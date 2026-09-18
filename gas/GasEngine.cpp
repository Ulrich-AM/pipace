#include "GasEngine.h"

#include "fluid/FluidEngine.h"
#include "fluid/DiagOutput.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/PhaseTransfer.h"
#include "thermal/ThermalTypes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace {

constexpr int kDx[4] = {-1, 1, 0, 0};
constexpr int kDy[4] = {0, 0, -1, 1};

int chunkIndex(int x, int y) {
    return (y / CHUNK) * CHUNK_W + (x / CHUNK);
}

} // namespace

GasEngine::GasEngine()
    : amount(GW * GH, 0.0f)
    , heat(GW * GH, 0.0f)
    , heatNext(GW * GH, 0.0f)
    , volume(GW * GH, 0.0f)
    , pressure(GW * GH, 0.0f)
    , u((GW + 1) * GH, 0.0f)
    , v(GW * (GH + 1), 0.0f)
    , fluxU((GW + 1) * GH, 0.0f)
    , fluxV(GW * (GH + 1), 0.0f)
    , outgoing(GW * GH, 0.0f)
    , chunkActivity(CHUNK_W * CHUNK_H, 0)
    , chunkQuietTicks(CHUNK_W * CHUNK_H, 0)
    , chunkSolveMask(CHUNK_W * CHUNK_H, 0)
    , thermalChunkWake(CHUNK_W * CHUNK_H, 0)
    , gasCompId(static_cast<size_t>(GW * GH) * kMaxGasComponents, SUBSTANCE_NONE)
    , gasCompAmt(static_cast<size_t>(GW * GH) * kMaxGasComponents, 0.0f)
    , gasCompCount(GW * GH, 0)
    , nextGasCompId(static_cast<size_t>(GW * GH) * kMaxGasComponents, SUBSTANCE_NONE)
    , nextGasCompAmt(static_cast<size_t>(GW * GH) * kMaxGasComponents, 0.0f)
    , nextGasCompCount(GW * GH, 0)
    , relocateStamp(GW * GH, 0)
    , prevVolume(GW * GH, 0.0f)
{
    relocateQueue.reserve(512);
}

float GasEngine::pressurePa(int index) const {
    return pressureAtm(index) * config.referencePressurePa;
}

float GasEngine::cellU(int x, int y) const {
    if (x < 0 || x >= GW || y < 0 || y >= GH) return 0.0f;
    return 0.5f * (u[static_cast<size_t>(ui(x, y))] + u[static_cast<size_t>(ui(x + 1, y))]);
}

float GasEngine::cellV(int x, int y) const {
    if (x < 0 || x >= GW || y < 0 || y >= GH) return 0.0f;
    return 0.5f * (v[static_cast<size_t>(vi(x, y))] + v[static_cast<size_t>(vi(x, y + 1))]);
}

bool GasEngine::isAccessible(FluidEngine const &fluid, int x, int y) const {
    return FluidEngine::inside(x, y) && availableVolume(fluid, x, y) >= GAS_MIN_VOLUME;
}

float GasEngine::availableVolume(FluidEngine const &fluid, int x, int y) const {
    if (!FluidEngine::inside(x, y)) return 0.0f;
    int i = ci(x, y);
    if (fluid.solid[static_cast<size_t>(i)] || fluid.dynamicSolid[static_cast<size_t>(i)]) return 0.0f;
    float fill = std::clamp(fluid.fill[static_cast<size_t>(i)], 0.0f, 1.0f);
    // Liquid is a gas barrier until bubbles exist. Partial-fill cells must not
    // compress leftover Air into a high-pressure cushion under floating bodies.
    if (fill >= MIN_PRESSURE_FILL) return 0.0f;
    return 1.0f;
}

double GasEngine::sumAmount() const {
    double s = 0.0;
    for (float a : amount) s += a;
    return s;
}

int GasEngine::compositionSlot(int cell, int slot) {
    return cell * kMaxGasComponents + slot;
}

void GasEngine::clearGasComposition(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        gasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_NONE;
        gasCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    gasCompCount[static_cast<size_t>(index)] = 0;
}

void GasEngine::compactGasComposition(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    int w = 0;
    for (int s = 0; s < n; ++s) {
        SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
        float amt = gasCompAmt[static_cast<size_t>(base + s)];
        if (!std::isfinite(amt) || amt < 0.0f) amt = 0.0f;
        if (amt > kMinGasComponent && validGasComponentId(id)) {
            if (w != s) {
                gasCompId[static_cast<size_t>(base + w)] = id;
                gasCompAmt[static_cast<size_t>(base + w)] = amt;
            }
            ++w;
        }
    }
    for (int s = w; s < n; ++s) {
        gasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_NONE;
        gasCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    gasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(w);
}

void GasEngine::copyGasCompToNext(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    nextGasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n);
    for (int s = 0; s < n; ++s) {
        nextGasCompId[static_cast<size_t>(base + s)] = gasCompId[static_cast<size_t>(base + s)];
        nextGasCompAmt[static_cast<size_t>(base + s)] = gasCompAmt[static_cast<size_t>(base + s)];
    }
    for (int s = n; s < kMaxGasComponents; ++s) {
        nextGasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_NONE;
        nextGasCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
}

void GasEngine::commitGasCompFromNext(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = nextGasCompCount[static_cast<size_t>(index)];
    if (n < 0) n = 0;
    if (n > kMaxGasComponents) n = kMaxGasComponents;
    gasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n);
    for (int s = 0; s < n; ++s) {
        gasCompId[static_cast<size_t>(base + s)] = nextGasCompId[static_cast<size_t>(base + s)];
        float amt = nextGasCompAmt[static_cast<size_t>(base + s)];
        if (!std::isfinite(amt) || amt < 0.0f) amt = 0.0f;
        gasCompAmt[static_cast<size_t>(base + s)] = amt;
    }
    for (int s = n; s < kMaxGasComponents; ++s) {
        gasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_NONE;
        gasCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    compactGasComposition(index);
}

float GasEngine::addGasComponentUntracked(int index, SubstanceId id, float componentAmount) {
    if (index < 0 || index >= GW * GH || !validGasComponentId(id)) return std::max(0.0f, componentAmount);
    if (!(componentAmount > kMinGasComponent) || !std::isfinite(componentAmount)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (gasCompId[static_cast<size_t>(base + s)] == id) {
            gasCompAmt[static_cast<size_t>(base + s)] += componentAmount;
            return 0.0f;
        }
    }
    if (n >= kMaxGasComponents) return componentAmount;
    gasCompId[static_cast<size_t>(base + n)] = id;
    gasCompAmt[static_cast<size_t>(base + n)] = componentAmount;
    gasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n + 1);
    return 0.0f;
}

float GasEngine::addNextGasComponentUntracked(int index, SubstanceId id, float componentAmount) {
    if (index < 0 || index >= GW * GH || !validGasComponentId(id)) return std::max(0.0f, componentAmount);
    if (!(componentAmount > kMinGasComponent) || !std::isfinite(componentAmount)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = nextGasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (nextGasCompId[static_cast<size_t>(base + s)] == id) {
            nextGasCompAmt[static_cast<size_t>(base + s)] += componentAmount;
            return 0.0f;
        }
    }
    if (n >= kMaxGasComponents) return componentAmount;
    nextGasCompId[static_cast<size_t>(base + n)] = id;
    nextGasCompAmt[static_cast<size_t>(base + n)] = componentAmount;
    nextGasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n + 1);
    return 0.0f;
}

void GasEngine::writePureGas(int index, SubstanceId id, float componentAmount) {
    clearGasComposition(index);
    if (!(componentAmount > GAS_MIN_AMOUNT) || !std::isfinite(componentAmount) || !validGasComponentId(id)) {
        amount[static_cast<size_t>(index)] = 0.0f;
        return;
    }
    addGasComponentUntracked(index, id, componentAmount);
    amount[static_cast<size_t>(index)] = componentAmount;
}

void GasEngine::scaleGasComposition(int index, float frac) {
    if (index < 0 || index >= GW * GH) return;
    if (!(frac > 0.0f) || !std::isfinite(frac)) {
        clearGasComposition(index);
        return;
    }
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s)
        gasCompAmt[static_cast<size_t>(base + s)] *= frac;
    compactGasComposition(index);
}

void GasEngine::syncAmountFromComposition(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    float sum = 0.0f;
    for (int s = 0; s < n; ++s) sum += gasCompAmt[static_cast<size_t>(base + s)];
    if (!std::isfinite(sum) || sum < 0.0f) sum = 0.0f;
    amount[static_cast<size_t>(index)] = sum;
}

float GasEngine::gasComponentAmount(int index, SubstanceId id) const {
    if (index < 0 || index >= GW * GH || !validGasComponentId(id)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (gasCompId[static_cast<size_t>(base + s)] == id) {
            float amt = gasCompAmt[static_cast<size_t>(base + s)];
            return (std::isfinite(amt) && amt > 0.0f) ? amt : 0.0f;
        }
    }
    return 0.0f;
}

float GasEngine::gasComponentFraction(int index, SubstanceId id) const {
    if (index < 0 || index >= GW * GH) return 0.0f;
    float a = amount[static_cast<size_t>(index)];
    if (!(a > GAS_MIN_AMOUNT)) return 0.0f;
    return std::clamp(gasComponentAmount(index, id) / a, 0.0f, 1.0f);
}

SubstanceId GasEngine::dominantGasSubstance(int index) const {
    if (index < 0 || index >= GW * GH) return SUBSTANCE_NONE;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    SubstanceId best = SUBSTANCE_NONE;
    float bestAmt = 0.0f;
    for (int s = 0; s < n; ++s) {
        SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
        float amt = gasCompAmt[static_cast<size_t>(base + s)];
        if (validGasComponentId(id) && amt > bestAmt) {
            bestAmt = amt;
            best = id;
        }
    }
    return best;
}

bool GasEngine::gasCompositionValid(int index) const {
    if (index < 0 || index >= GW * GH) return false;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    if (n < 0 || n > kMaxGasComponents) return false;
    float sum = 0.0f;
    for (int s = 0; s < n; ++s) {
        SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
        float amt = gasCompAmt[static_cast<size_t>(base + s)];
        if (!(amt > kMinGasComponent)) continue;
        if (!validGasComponentId(id)) return false;
        for (int b = 0; b < s; ++b) {
            if (gasCompId[static_cast<size_t>(base + b)] == id) return false;
        }
        sum += amt;
    }
    float a = amount[static_cast<size_t>(index)];
    if (!(a > GAS_MIN_AMOUNT))
        return !(sum > GAS_MIN_AMOUNT);
    if (!(sum > GAS_MIN_AMOUNT)) return false; // positive amount must have composition
    return std::abs(sum - a) <= 1.0e-4f * std::max(1.0f, a) + 1.0e-6f;
}

GasComponentView GasEngine::gasComponents(int index) const {
    GasComponentView view;
    if (index < 0 || index >= GW * GH) return view;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n && view.count < kMaxGasComponents; ++s) {
        SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
        float amt = gasCompAmt[static_cast<size_t>(base + s)];
        if (amt > kMinGasComponent && validGasComponentId(id))
            view.items[view.count++] = {id, amt};
    }
    return view;
}

void GasEngine::setGasComponentAmount(int index, SubstanceId id, float componentAmount) {
    if (index < 0 || index >= GW * GH || !validGasComponentId(id)) return;
    if (!std::isfinite(componentAmount) || componentAmount < 0.0f) componentAmount = 0.0f;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[static_cast<size_t>(index)];
    int slot = -1;
    for (int s = 0; s < n; ++s) {
        if (gasCompId[static_cast<size_t>(base + s)] == id) { slot = s; break; }
    }
    if (componentAmount <= kMinGasComponent) {
        if (slot >= 0) {
            gasCompAmt[static_cast<size_t>(base + slot)] = 0.0f;
            compactGasComposition(index);
        }
        syncAmountFromComposition(index);
        return;
    }
    if (slot < 0) {
        if (n >= kMaxGasComponents) return;
        slot = n;
        gasCompId[static_cast<size_t>(base + slot)] = id;
        gasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n + 1);
    }
    gasCompAmt[static_cast<size_t>(base + slot)] = componentAmount;
    compactGasComposition(index);
    syncAmountFromComposition(index);
}

void GasEngine::addGasComponentAmount(int index, SubstanceId id, float delta) {
    if (!(delta > 0.0f) || !std::isfinite(delta)) {
        if (delta < 0.0f) (void)takeGasComponentAmount(index, id, -delta);
        return;
    }
    if (addGasComponentUntracked(index, id, delta) > 0.0f) return;
    amount[static_cast<size_t>(index)] += delta;
}

float GasEngine::takeGasComponentAmount(int index, SubstanceId id, float da) {
    if (index < 0 || index >= GW * GH || !validGasComponentId(id) || !(da > 0.0f)) return 0.0f;
    float have = gasComponentAmount(index, id);
    float take = std::min(have, da);
    if (take <= kMinGasComponent) return 0.0f;
    size_t i = static_cast<size_t>(index);
    float a = amount[i];
    if (a > GAS_MIN_AMOUNT) heat[i] *= std::max(0.0f, (a - take) / a);
    setGasComponentAmount(index, id, have - take);
    if (amount[i] <= GAS_MIN_AMOUNT) heat[i] = 0.0f;
    return take;
}

bool GasEngine::tryCommitGasOccupancy(int index, GasComponentView const &view) {
    if (index < 0 || index >= GW * GH) return false;
    if (view.count < 0 || view.count > kMaxGasComponents) return false;
    GasComponent packed[kMaxGasComponents]{};
    int packedCount = 0;
    float sum = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId id = view.items[n].id;
        float amt = view.items[n].amount;
        if (!(amt > kMinGasComponent)) continue;
        if (!validGasComponentId(id)) return false;
        if (addGasPayload(packed, packedCount, id, amt) > 0.0f) return false;
        sum += amt;
    }
    if (!std::isfinite(sum) || sum < 0.0f) return false;
    int base = compositionSlot(index, 0);
    for (int s = 0; s < kMaxGasComponents; ++s) {
        if (s < packedCount) {
            gasCompId[static_cast<size_t>(base + s)] = packed[s].id;
            gasCompAmt[static_cast<size_t>(base + s)] = packed[s].amount;
        } else {
            gasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_NONE;
            gasCompAmt[static_cast<size_t>(base + s)] = 0.0f;
        }
    }
    gasCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(packedCount);
    amount[static_cast<size_t>(index)] = sum;
    if (!(sum > GAS_MIN_AMOUNT)) heat[static_cast<size_t>(index)] = 0.0f;
    return true;
}

float GasEngine::gasPartialPressurePa(int index, SubstanceId id) const {
    if (index < 0 || index >= GW * GH) return 0.0f;
    float vol = volume[static_cast<size_t>(index)];
    if (!(vol >= GAS_MIN_VOLUME)) return 0.0f;
    float amt = gasComponentAmount(index, id);
    if (!(amt > GAS_MIN_AMOUNT)) return 0.0f;
    return (amt / vol) * config.referencePressurePa;
}

double GasEngine::sumWaterVapor() const {
    double s = 0.0;
    for (int i = 0; i < GW * GH; ++i)
        s += static_cast<double>(gasComponentAmount(i, SUBSTANCE_WATER));
    return s;
}

float GasEngine::vaporAmount(int index) const {
    return gasComponentAmount(index, SUBSTANCE_WATER);
}

float GasEngine::airAmount(int index) const {
    return gasComponentAmount(index, SUBSTANCE_AIR);
}

float GasEngine::vaporFraction(int index) const {
    return gasComponentFraction(index, SUBSTANCE_WATER);
}

void GasEngine::clampSpecies(int index) {
    if (index < 0 || index >= GW * GH) return;
    size_t i = static_cast<size_t>(index);
    if (!std::isfinite(amount[i]) || amount[i] < 0.0f) amount[i] = 0.0f;
    compactGasComposition(index);
    float sum = 0.0f;
    int base = compositionSlot(index, 0);
    int n = gasCompCount[i];
    for (int s = 0; s < n; ++s) sum += gasCompAmt[static_cast<size_t>(base + s)];
    if (!std::isfinite(sum) || sum < 0.0f) sum = 0.0f;
    // Species are conserved mass. If bookkeeping drifted, keep components —
    // do not invent missing Air and do not delete leftover vapor.
    if (sum > amount[i]) amount[i] = sum;
    if (!(amount[i] > GAS_MIN_AMOUNT) && !(sum > GAS_MIN_AMOUNT)) {
        amount[i] = 0.0f;
        clearGasComposition(index);
        return;
    }
    if (!(sum > GAS_MIN_AMOUNT)) {
        // Positive amount with empty composition is invalid. Do not invent Air.
        amount[i] = 0.0f;
        clearGasComposition(index);
    }
}

void GasEngine::addWaterVapor(int index, float da) {
    addGasComponentAmount(index, SUBSTANCE_WATER, da);
}

float GasEngine::takeWaterVapor(int index, float da) {
    return takeGasComponentAmount(index, SUBSTANCE_WATER, da);
}

void GasEngine::wakeAt(int x, int y) {
    wakeChunkAtCell(x, y);
    wakeThermalAt(x, y);
}

void GasEngine::transferSpecies(int donor, int receiver, float q) {
    if (q <= GAS_MIN_AMOUNT || donor == receiver) return;
    float a0 = amount[static_cast<size_t>(donor)];
    if (a0 <= GAS_MIN_AMOUNT) return;
    GasComponentView src = gasComponents(donor);
    GasComponentView dst = gasComponents(receiver);
    if (!gasPayloadCanMerge(dst.items, dst.count, src.items, src.count)) return;
    float frac = q / a0;
    int base = compositionSlot(donor, 0);
    int n = gasCompCount[static_cast<size_t>(donor)];
    for (int s = 0; s < n; ++s) {
        SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
        float d = gasCompAmt[static_cast<size_t>(base + s)] * frac;
        if (!(d > kMinGasComponent) || !validGasComponentId(id)) continue;
        gasCompAmt[static_cast<size_t>(base + s)] -= d;
        float unplaced = addGasComponentUntracked(receiver, id, d);
        if (unplaced > kMinGasComponent)
            gasCompAmt[static_cast<size_t>(base + s)] += unplaced;
    }
    compactGasComposition(donor);
    compactGasComposition(receiver);
}

void GasEngine::commitExpected() {
    currentAmount = sumAmount();
    expectedAmount = currentAmount;
    amountError = 0.0;
    currentWaterVapor = sumWaterVapor();
    expectedWaterVapor = currentWaterVapor;
}

void GasEngine::wakeChunkAtCell(int x, int y, bool resetQuiet) {
    if (!FluidEngine::inside(x, y)) return;
    int c = chunkIndex(x, y);
    chunkActivity[static_cast<size_t>(c)] = 1;
    if (resetQuiet) chunkQuietTicks[static_cast<size_t>(c)] = 0;
}

void GasEngine::wakeThermalAt(int x, int y) {
    if (!FluidEngine::inside(x, y)) return;
    int c = chunkIndex(x, y);
    if (c >= 0 && c < static_cast<int>(thermalChunkWake.size())) thermalChunkWake[static_cast<size_t>(c)] = 1;
}

void GasEngine::wakeAll() {
    std::fill(chunkActivity.begin(), chunkActivity.end(), uint8_t{1});
    std::fill(chunkQuietTicks.begin(), chunkQuietTicks.end(), uint8_t{0});
}

void GasEngine::sleepAll() {
    std::fill(chunkActivity.begin(), chunkActivity.end(), uint8_t{0});
    std::fill(chunkQuietTicks.begin(), chunkQuietTicks.end(), uint8_t{255});
}

void GasEngine::recomputePressure() {
    for (int i = 0; i < GW * GH; ++i) {
        float vol = volume[static_cast<size_t>(i)];
        if (vol < GAS_MIN_VOLUME) pressure[static_cast<size_t>(i)] = 0.0f;
        else pressure[static_cast<size_t>(i)] = amount[static_cast<size_t>(i)] / vol;
    }
}

void GasEngine::rebuildVolumes(FluidEngine const &fluid, bool &volumeChanged, bool fullGrid) {
    volumeChanged = false;
    auto writeVol = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return;
        int i = ci(x, y);
        float vol = availableVolume(fluid, x, y);
        if (std::abs(vol - volume[static_cast<size_t>(i)]) > 1.0e-6f) {
            volumeChanged = true;
            wakeChunkAtCell(x, y);
        }
        prevVolume[static_cast<size_t>(i)] = volume[static_cast<size_t>(i)];
        volume[static_cast<size_t>(i)] = vol;
    };

    ++volumeScanTick;
    bool scanAll = fullGrid || (volumeScanTick % 16) == 0;
    if (scanAll) {
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) writeVol(x, y);
        lastOccupancyScan = fluid.dynamicOccupiedCells;
        return;
    }

    if (fluid.hasActiveSolveRegion) {
        int x0 = std::max(0, fluid.solveX0 - 2);
        int y0 = std::max(0, fluid.solveY0 - 2);
        int x1 = std::min(GW - 1, fluid.solveX1 + 2);
        int y1 = std::min(GH - 1, fluid.solveY1 + 2);
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) writeVol(x, y);
    }
    for (int i : lastOccupancyScan) writeVol(i % GW, i / GW);
    for (int i : fluid.dynamicOccupiedCells) writeVol(i % GW, i / GW);
    lastOccupancyScan = fluid.dynamicOccupiedCells;
}

float GasEngine::relocateAmount(FluidEngine const &fluid, int x, int y, GasComponentView parcel,
    float leftoverHeatPerAmount, float maxAtm)
{
    compactGasPayload(parcel.items, parcel.count);
    float leftover = gasPayloadSum(parcel.items, parcel.count);
    if (leftover <= GAS_MIN_AMOUNT) return 0.0f;
    if (++relocateEpoch == 0) {
        std::fill(relocateStamp.begin(), relocateStamp.end(), 0);
        relocateEpoch = 1;
    }
    relocateQueue.clear();
    auto consider = [&](int nx, int ny) {
        if (!isAccessible(fluid, nx, ny)) return;
        int ni = ci(nx, ny);
        if (relocateStamp[static_cast<size_t>(ni)] == relocateEpoch) return;
        relocateStamp[static_cast<size_t>(ni)] = relocateEpoch;
        relocateQueue.push_back(ni);
    };
    for (int n = 0; n < 4; ++n) consider(x + kDx[n], y + kDy[n]);
    size_t limit = static_cast<size_t>(GW * GH);
    for (size_t head = 0; head < relocateQueue.size() && leftover > GAS_MIN_AMOUNT && head < limit; ++head) {
        int i = relocateQueue[head];
        int cx = i % GW, cy = i / GW;
        float vol = volume[static_cast<size_t>(i)];
        if (vol >= GAS_MIN_VOLUME) {
            GasComponentView dest = gasComponents(i);
            if (!gasPayloadCanMerge(dest.items, dest.count, parcel.items, parcel.count)) {
                for (int n = 0; n < 4; ++n) consider(cx + kDx[n], cy + kDy[n]);
                continue;
            }
            float room = leftover;
            if (maxAtm > 0.0f) {
                float cap = vol * maxAtm;
                room = std::max(0.0f, cap - amount[static_cast<size_t>(i)]);
            }
            float placed = std::min(leftover, room);
            if (placed > GAS_MIN_AMOUNT) {
                float frac = placed / leftover;
                float accepted = 0.0f;
                for (int s = 0; s < parcel.count; ++s) {
                    float d = parcel.items[s].amount * frac;
                    float unplaced = addGasComponentUntracked(i, parcel.items[s].id, d);
                    float got = d - unplaced;
                    parcel.items[s].amount -= got;
                    accepted += got;
                }
                compactGasPayload(parcel.items, parcel.count);
                amount[static_cast<size_t>(i)] += accepted;
                clampSpecies(i);
                heat[static_cast<size_t>(i)] += accepted * leftoverHeatPerAmount;
                leftover = gasPayloadSum(parcel.items, parcel.count);
                wakeChunkAtCell(cx, cy);
                if (std::abs(leftoverHeatPerAmount) > 1.0e-6f) wakeThermalAt(cx, cy);
            }
        }
        for (int n = 0; n < 4; ++n) consider(cx + kDx[n], cy + kDy[n]);
    }
    return leftover;
}

void GasEngine::displaceBlocked(FluidEngine const &fluid) {
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = ci(x, y);
        if (volume[static_cast<size_t>(i)] >= GAS_MIN_VOLUME) continue;
        clampSpecies(i);
        float a = amount[static_cast<size_t>(i)];
        if (a <= GAS_MIN_AMOUNT) {
            amount[static_cast<size_t>(i)] = 0.0f;
            clearGasComposition(i);
            heat[static_cast<size_t>(i)] = 0.0f;
            continue;
        }
        float h = heat[static_cast<size_t>(i)];
        GasComponentView parcel = gasComponents(i);
        amount[static_cast<size_t>(i)] = 0.0f;
        clearGasComposition(i);
        heat[static_cast<size_t>(i)] = 0.0f;
        float maxAtm = 0.0f;
        if (!fluid.solid[static_cast<size_t>(i)] && !fluid.dynamicSolid[static_cast<size_t>(i)]
            && fluid.fill[static_cast<size_t>(i)] >= MIN_PRESSURE_FILL)
            maxAtm = config.ambientPressureAtm * 1.12f;
        float rem = relocateAmount(fluid, x, y, parcel, a > GAS_MIN_AMOUNT ? h / a : 0.0f, maxAtm);
        if (rem > 0.0f && std::isfinite(rem)) {
            for (int s = 0; s < parcel.count; ++s)
                (void)addGasComponentUntracked(i, parcel.items[s].id, parcel.items[s].amount);
            amount[static_cast<size_t>(i)] += rem;
            clampSpecies(i);
            heat[static_cast<size_t>(i)] += rem * (a > GAS_MIN_AMOUNT ? h / a : 0.0f);
        }
    }
}

void GasEngine::displaceLiquidOverflow(FluidEngine const &fluid) {
    // Partial liquid: push gas out of lost volume so we do not form trapped bubbles.
    // Solid occupancy uses compression instead (piston / sealed-cavity tests).
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = ci(x, y);
        float oldV = prevVolume[static_cast<size_t>(i)];
        float newV = volume[static_cast<size_t>(i)];
        if (newV >= oldV - 1.0e-6f) continue;
        if (fluid.solid[static_cast<size_t>(i)] || fluid.dynamicSolid[static_cast<size_t>(i)]) continue;
        float a = amount[static_cast<size_t>(i)];
        if (a <= GAS_MIN_AMOUNT || oldV < GAS_MIN_VOLUME) continue;
        float keep = a * (newV / oldV);
        float pushed = a - keep;
        if (pushed <= GAS_MIN_AMOUNT) continue;
        float h = heat[static_cast<size_t>(i)];
        float hKeep = a > GAS_MIN_AMOUNT ? h * (keep / a) : 0.0f;
        float hPush = h - hKeep;
        GasComponentView view = gasComponents(i);
        GasComponentView parcel;
        copyGasPayload(parcel.items, parcel.count, view.items, view.count);
        scaleGasPayload(parcel.items, parcel.count, pushed / a);
        compactGasPayload(parcel.items, parcel.count);
        scaleGasComposition(i, keep / a);
        amount[static_cast<size_t>(i)] = keep;
        clampSpecies(i);
        heat[static_cast<size_t>(i)] = hKeep;
        float rem = relocateAmount(fluid, x, y, parcel,
            pushed > GAS_MIN_AMOUNT ? hPush / pushed : 0.0f, config.ambientPressureAtm * 1.12f);
        if (rem > GAS_MIN_AMOUNT) {
            for (int s = 0; s < parcel.count; ++s)
                (void)addGasComponentUntracked(i, parcel.items[s].id, parcel.items[s].amount);
            amount[static_cast<size_t>(i)] += rem;
            clampSpecies(i);
            heat[static_cast<size_t>(i)] += rem * (hPush / pushed);
        }
        wakeChunkAtCell(x, y);
    }
}

void GasEngine::applyMovingBoundary(FluidEngine const &fluid) {
    for (int i : fluid.dynamicOccupiedCells) {
        int x = i % GW, y = i / GW;
        float vx = fluid.dynamicVelX[static_cast<size_t>(i)];
        float vy = fluid.dynamicVelY[static_cast<size_t>(i)];
        u[static_cast<size_t>(ui(x, y))] = vx;
        u[static_cast<size_t>(ui(x + 1, y))] = vx;
        v[static_cast<size_t>(vi(x, y))] = vy;
        v[static_cast<size_t>(vi(x, y + 1))] = vy;
        wakeChunkAtCell(x, y);
        for (int n = 0; n < 4; ++n) wakeChunkAtCell(x + kDx[n], y + kDy[n]);
    }
}

void GasEngine::accumulateTransfer(FluidEngine const &fluid, float dt) {
    (void)fluid;
    std::fill(fluxU.begin(), fluxU.end(), 0.0f);
    std::fill(fluxV.begin(), fluxV.end(), 0.0f);
    std::fill(outgoing.begin(), outgoing.end(), 0.0f);

    auto faceFlux = [&](int i0, int i1, float adv) {
        float vol0 = volume[static_cast<size_t>(i0)];
        float vol1 = volume[static_cast<size_t>(i1)];
        if (vol0 < GAS_MIN_VOLUME || vol1 < GAS_MIN_VOLUME) return 0.0f;
        float p0 = pressure[static_cast<size_t>(i0)];
        float p1 = pressure[static_cast<size_t>(i1)];
        float volFace = 2.0f * vol0 * vol1 / (vol0 + vol1);
        float flux = config.conductivity * (p0 - p1) * volFace * dt;
        float dens0 = amount[static_cast<size_t>(i0)] / vol0;
        float dens1 = amount[static_cast<size_t>(i1)] / vol1;
        float up = adv > 0.0f ? dens0 : dens1;
        flux += adv * dt * up * volFace;
        flux = std::clamp(flux, -config.maxFaceTransfer, config.maxFaceTransfer);
        return flux;
    };

    int x0 = 0, y0 = 0, x1 = GW - 1, y1 = GH - 1;
    bool sparse = false;
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) if (chunkActivity[static_cast<size_t>(c)] || chunkSolveMask[static_cast<size_t>(c)]) {
        sparse = true;
        break;
    }
    if (sparse) {
        x0 = GW; y0 = GH; x1 = -1; y1 = -1;
        for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
            if (!chunkActivity[static_cast<size_t>(c)] && !chunkSolveMask[static_cast<size_t>(c)]) continue;
            int cx = c % CHUNK_W, cy = c / CHUNK_W;
            x0 = std::min(x0, std::max(0, cx * CHUNK - 1));
            y0 = std::min(y0, std::max(0, cy * CHUNK - 1));
            x1 = std::max(x1, std::min(GW - 1, (cx + 1) * CHUNK));
            y1 = std::max(y1, std::min(GH - 1, (cy + 1) * CHUNK));
        }
        if (x1 < x0) return;
    }

    for (int y = y0; y <= y1; ++y) for (int x = x0 + 1; x <= x1; ++x) {
        int iL = ci(x - 1, y), iR = ci(x, y);
        float adv = u[static_cast<size_t>(ui(x, y))];
        float flux = faceFlux(iL, iR, adv);
        fluxU[static_cast<size_t>(ui(x, y))] = flux;
        if (flux > 0.0f) outgoing[static_cast<size_t>(iL)] += flux;
        else outgoing[static_cast<size_t>(iR)] += -flux;
    }
    for (int y = y0 + 1; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int iT = ci(x, y - 1), iB = ci(x, y);
        float adv = v[static_cast<size_t>(vi(x, y))];
        float flux = faceFlux(iT, iB, adv);
        fluxV[static_cast<size_t>(vi(x, y))] = flux;
        if (flux > 0.0f) outgoing[static_cast<size_t>(iT)] += flux;
        else outgoing[static_cast<size_t>(iB)] += -flux;
    }

    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        float cap = amount[static_cast<size_t>(i)] * config.maxExtractFraction;
        float out = outgoing[static_cast<size_t>(i)];
        if (out <= cap || out <= GAS_MIN_AMOUNT) continue;
        float s = cap / out;
        if (x > 0) {
            float &f = fluxU[static_cast<size_t>(ui(x, y))];
            if (f > 0.0f) f *= s;
        }
        if (x + 1 < GW) {
            float &f = fluxU[static_cast<size_t>(ui(x + 1, y))];
            if (f < 0.0f) f *= s;
        }
        if (y > 0) {
            float &f = fluxV[static_cast<size_t>(vi(x, y))];
            if (f > 0.0f) f *= s;
        }
        if (y + 1 < GH) {
            float &f = fluxV[static_cast<size_t>(vi(x, y + 1))];
            if (f < 0.0f) f *= s;
        }
    }
}

void GasEngine::applyBoundaryFlux(float dt) {
    if (config.boundary != GasBoundary::OpenAmbient) return;
    float pAmb = config.ambientPressureAtm;
    auto edge = [&](int i) {
        float vol = volume[static_cast<size_t>(i)];
        if (vol < GAS_MIN_VOLUME) return;
        float p = pressure[static_cast<size_t>(i)];
        float flux = config.conductivity * (p - pAmb) * vol * dt;
        flux = std::clamp(flux, -config.maxFaceTransfer, config.maxFaceTransfer);
        if (flux > 0.0f) flux = std::min(flux, amount[static_cast<size_t>(i)] * config.maxExtractFraction);
        if (flux > 0.0f) {
            float a = amount[static_cast<size_t>(i)];
            float dq = (a > GAS_MIN_AMOUNT) ? heat[static_cast<size_t>(i)] * (flux / a) : 0.0f;
            float dv = (a > GAS_MIN_AMOUNT) ? gasComponentAmount(i, SUBSTANCE_WATER) * (flux / a) : 0.0f;
            heat[static_cast<size_t>(i)] -= dq;
            scaleGasComposition(i, a > GAS_MIN_AMOUNT ? (a - flux) / a : 0.0f);
            escapedWaterVapor += dv;
            expectedWaterVapor -= dv;
            escapedHeat += dq;
        } else if (flux < 0.0f) {
            float cap = thermalCapacity(gasMassKg(-flux), thermalForSubstance(SUBSTANCE_AIR).specificHeat);
            float added = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
            heat[static_cast<size_t>(i)] += added;
            escapedHeat -= added;
            (void)addGasComponentUntracked(i, SUBSTANCE_AIR, -flux);
        }
        amount[static_cast<size_t>(i)] -= flux;
        if (amount[static_cast<size_t>(i)] < 0.0f) {
            flux += amount[static_cast<size_t>(i)];
            amount[static_cast<size_t>(i)] = 0.0f;
        }
        clampSpecies(i);
        escapedAmount += flux;
        expectedAmount -= flux;
        if (std::abs(flux) > 1.0e-6f) {
            int x = i % GW, y = i / GW;
            wakeChunkAtCell(x, y);
            wakeThermalAt(x, y);
        }
    };
    for (int y = 0; y < GH; ++y) {
        edge(ci(0, y));
        edge(ci(GW - 1, y));
    }
    for (int x = 1; x < GW - 1; ++x) {
        edge(ci(x, 0));
        edge(ci(x, GH - 1));
    }
}

void GasEngine::applyFluxes() {
    heatNext = heat;
    for (int i = 0; i < GW * GH; ++i) copyGasCompToNext(i);
    auto nextView = [&](int index) {
        GasComponentView view;
        int base = compositionSlot(index, 0);
        int n = nextGasCompCount[static_cast<size_t>(index)];
        for (int s = 0; s < n && view.count < kMaxGasComponents; ++s) {
            SubstanceId id = nextGasCompId[static_cast<size_t>(base + s)];
            float amt = nextGasCompAmt[static_cast<size_t>(base + s)];
            if (amt > kMinGasComponent)
                view.items[view.count++] = {id, amt};
        }
        return view;
    };
    auto mixFlux = [&](int donor, int receiver, float q) -> float {
        if (q <= 0.0f) return 0.0f;
        float a0 = amount[static_cast<size_t>(donor)];
        if (a0 <= GAS_MIN_AMOUNT) return 0.0f;
        GasComponentView src = gasComponents(donor);
        GasComponentView dst = nextView(receiver);
        if (!gasPayloadCanMerge(dst.items, dst.count, src.items, src.count))
            return 0.0f;
        float frac = q / a0;
        float dq = heat[static_cast<size_t>(donor)] * frac;
        heatNext[static_cast<size_t>(donor)] -= dq;
        heatNext[static_cast<size_t>(receiver)] += dq;
        int base = compositionSlot(donor, 0);
        int n = gasCompCount[static_cast<size_t>(donor)];
        for (int s = 0; s < n; ++s) {
            SubstanceId id = gasCompId[static_cast<size_t>(base + s)];
            float d = gasCompAmt[static_cast<size_t>(base + s)] * frac;
            if (!(d > 0.0f) || !validGasComponentId(id)) continue;
            float unplaced = addNextGasComponentUntracked(receiver, id, d);
            float got = d - unplaced;
            int nbase = compositionSlot(donor, 0);
            for (int t = 0; t < nextGasCompCount[static_cast<size_t>(donor)]; ++t) {
                if (nextGasCompId[static_cast<size_t>(nbase + t)] == id) {
                    nextGasCompAmt[static_cast<size_t>(nbase + t)] -= got;
                    break;
                }
            }
        }
        return q;
    };
    for (int y = 0; y < GH; ++y) for (int x = 1; x < GW; ++x) {
        float &flux = fluxU[static_cast<size_t>(ui(x, y))];
        if (flux == 0.0f) continue;
        int iL = ci(x - 1, y), iR = ci(x, y);
        int donor = flux > 0.0f ? iL : iR;
        int receiver = flux > 0.0f ? iR : iL;
        float accepted = mixFlux(donor, receiver, std::abs(flux));
        if (!(accepted > 0.0f)) {
            flux = 0.0f;
            continue;
        }
        if (std::abs(flux) > 1.0e-3f) { wakeThermalAt(x - 1, y); wakeThermalAt(x, y); }
        if (accepted + 1.0e-8f < std::abs(flux))
            flux = (flux > 0.0f ? accepted : -accepted);
    }
    for (int y = 1; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        float &flux = fluxV[static_cast<size_t>(vi(x, y))];
        if (flux == 0.0f) continue;
        int iT = ci(x, y - 1), iB = ci(x, y);
        int donor = flux > 0.0f ? iT : iB;
        int receiver = flux > 0.0f ? iB : iT;
        float accepted = mixFlux(donor, receiver, std::abs(flux));
        if (!(accepted > 0.0f)) {
            flux = 0.0f;
            continue;
        }
        if (std::abs(flux) > 1.0e-3f) { wakeThermalAt(x, y - 1); wakeThermalAt(x, y); }
        if (accepted + 1.0e-8f < std::abs(flux))
            flux = (flux > 0.0f ? accepted : -accepted);
    }
    for (int y = 0; y < GH; ++y) for (int x = 1; x < GW; ++x) {
        float flux = fluxU[static_cast<size_t>(ui(x, y))];
        if (flux == 0.0f) continue;
        int iL = ci(x - 1, y), iR = ci(x, y);
        amount[static_cast<size_t>(iL)] -= flux;
        amount[static_cast<size_t>(iR)] += flux;
        if (std::abs(flux) > 1.0e-5f) {
            wakeChunkAtCell(x - 1, y, false);
            wakeChunkAtCell(x, y, false);
        }
    }
    for (int y = 1; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        float flux = fluxV[static_cast<size_t>(vi(x, y))];
        if (flux == 0.0f) continue;
        int iT = ci(x, y - 1), iB = ci(x, y);
        amount[static_cast<size_t>(iT)] -= flux;
        amount[static_cast<size_t>(iB)] += flux;
        if (std::abs(flux) > 1.0e-5f) {
            wakeChunkAtCell(x, y - 1, false);
            wakeChunkAtCell(x, y, false);
        }
    }
    for (int i = 0; i < GW * GH; ++i) {
        if (amount[static_cast<size_t>(i)] < 0.0f) amount[static_cast<size_t>(i)] = 0.0f;
        commitGasCompFromNext(i);
        clampSpecies(i);
        if (amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT) heat[static_cast<size_t>(i)] = 0.0f;
        else {
            heat[static_cast<size_t>(i)] = heatNext[static_cast<size_t>(i)];
            if (!std::isfinite(heat[static_cast<size_t>(i)]) || heat[static_cast<size_t>(i)] < 0.0f)
                heat[static_cast<size_t>(i)] = 0.0f;
        }
    }
}

void GasEngine::integrateVelocity(FluidEngine const &fluid, float dt) {
    float damp = std::clamp(1.0f - config.velocityDamping, 0.0f, 1.0f);
    float gravity = fluid.gridGravity();
    float buoy = config.buoyancyScale;
    auto gasT = [&](int index) -> float {
        float a = amount[static_cast<size_t>(index)];
        if (a <= GAS_MIN_AMOUNT) return AMBIENT_TEMPERATURE_K;
        float h = heat[static_cast<size_t>(index)];
        if (!(h > 0.0f)) return AMBIENT_TEMPERATURE_K;
        float cap = gasMixtureThermalCapacity(gasComponents(index), fluid.config.cellsPerMeter);
        if (cap < MIN_THERMAL_CAPACITY) return AMBIENT_TEMPERATURE_K;
        return tempFromEnergy(h, cap);
    };
    for (int y = 0; y < GH; ++y) for (int x = 1; x < GW; ++x) {
        int iL = ci(x - 1, y), iR = ci(x, y);
        float &vel = u[static_cast<size_t>(ui(x, y))];
        if (volume[static_cast<size_t>(iL)] < GAS_MIN_VOLUME || volume[static_cast<size_t>(iR)] < GAS_MIN_VOLUME) {
            vel = 0.0f;
            continue;
        }
        float grad = pressure[static_cast<size_t>(iR)] - pressure[static_cast<size_t>(iL)];
        vel += -config.pressureAccel * grad * dt;
        vel *= damp;
        vel = std::clamp(vel, -config.maxVelocity, config.maxVelocity);
    }
    for (int y = 1; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int iT = ci(x, y - 1), iB = ci(x, y);
        float &vel = v[static_cast<size_t>(vi(x, y))];
        if (volume[static_cast<size_t>(iT)] < GAS_MIN_VOLUME || volume[static_cast<size_t>(iB)] < GAS_MIN_VOLUME) {
            vel = 0.0f;
            continue;
        }
        float grad = pressure[static_cast<size_t>(iB)] - pressure[static_cast<size_t>(iT)];
        vel += -config.pressureAccel * grad * dt;
        if (gravity > 0.0f) {
            int cT = chunkIndex(x, y - 1);
            int cB = chunkIndex(x, y);
            bool live = (cT >= 0 && cT < CHUNK_W * CHUNK_H
                    && (chunkActivity[static_cast<size_t>(cT)] || chunkSolveMask[static_cast<size_t>(cT)]))
                || (cB >= 0 && cB < CHUNK_W * CHUNK_H
                    && (chunkActivity[static_cast<size_t>(cB)] || chunkSolveMask[static_cast<size_t>(cB)]));
            if (live) {
                if (buoy > 1.0e-8f) {
                    float Tface = 0.5f * (gasT(iT) + gasT(iB));
                    float dT = std::clamp(Tface - AMBIENT_TEMPERATURE_K, -400.0f, 400.0f);
                    float accel = -gravity * buoy * (dT / AMBIENT_TEMPERATURE_K);
                    accel = std::clamp(accel, -3.0f * gravity, 3.0f * gravity);
                    vel += accel * dt;
                }
                float compScale = config.compositionBuoyancyScale;
                if (compScale > 1.0e-8f) {
                    auto ambientAir = [&](int idx) {
                        return gasCompCount[static_cast<size_t>(idx)] == 1
                            && gasCompId[static_cast<size_t>(compositionSlot(idx, 0))] == SUBSTANCE_AIR;
                    };
                    if (!(ambientAir(iT) && ambientAir(iB))) {
                        float rhoAir = AIR_DENSITY_KG_M3;
                        float rhoT = gasMixtureReferenceDensityKgM3(gasComponents(iT));
                        float rhoB = gasMixtureReferenceDensityKgM3(gasComponents(iB));
                        float rhoFace = 0.5f * (rhoT + rhoB);
                        if (rhoAir > 1.0e-8f && rhoFace > 0.0f && std::isfinite(rhoFace)) {
                            float ratio = rhoFace / rhoAir;
                            if (!std::isfinite(ratio)) ratio = 1.0f;
                            ratio = std::clamp(ratio, 0.05f, 8.0f);
                            float accel = gravity * compScale * (ratio - 1.0f);
                            accel = std::clamp(accel, -2.0f * gravity, 2.0f * gravity);
                            vel += accel * dt;
                        }
                    }
                }
            }
        }
        vel *= damp;
        vel = std::clamp(vel, -config.maxVelocity, config.maxVelocity);
    }
    applyMovingBoundary(fluid);
}

void GasEngine::rebuildActivity(bool advanceSleep) {
    std::array<float, CHUNK_W * CHUNK_H> maxDp{};
    std::array<float, CHUNK_W * CHUNK_H> maxSpeed{};
    activeGasCells = 0;
    currentAmount = sumAmount();
    amountError = currentAmount - expectedAmount;

    std::fill(chunkSolveMask.begin(), chunkSolveMask.end(), 0);
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) if (chunkActivity[static_cast<size_t>(c)]) {
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
            int nx = cx + ox, ny = cy + oy;
            if (nx >= 0 && nx < CHUNK_W && ny >= 0 && ny < CHUNK_H) chunkSolveMask[static_cast<size_t>(ny * CHUNK_W + nx)] = 1;
        }
    }

    for (int chunk = 0; chunk < CHUNK_W * CHUNK_H; ++chunk) {
        if (!chunkActivity[static_cast<size_t>(chunk)] && !chunkSolveMask[static_cast<size_t>(chunk)]) continue;
        int cx = chunk % CHUNK_W, cy = chunk / CHUNK_W;
        int beginX = cx * CHUNK, endX = std::min(GW, beginX + CHUNK);
        int beginY = cy * CHUNK, endY = std::min(GH, beginY + CHUNK);
        for (int y = beginY; y < endY; ++y) for (int x = beginX; x < endX; ++x) {
            int i = ci(x, y);
            if (volume[static_cast<size_t>(i)] < GAS_MIN_VOLUME) continue;
            if (chunkActivity[static_cast<size_t>(chunk)]) ++activeGasCells;
            float speed = std::abs(cellU(x, y)) + std::abs(cellV(x, y));
            maxSpeed[static_cast<size_t>(chunk)] = std::max(maxSpeed[static_cast<size_t>(chunk)], speed);
            for (int n = 0; n < 4; ++n) {
                int nx = x + kDx[n], ny = y + kDy[n];
                if (!FluidEngine::inside(nx, ny) || volume[static_cast<size_t>(ci(nx, ny))] < GAS_MIN_VOLUME) continue;
                maxDp[static_cast<size_t>(chunk)] = std::max(maxDp[static_cast<size_t>(chunk)],
                    std::abs(pressure[static_cast<size_t>(i)] - pressure[static_cast<size_t>(ci(nx, ny))]));
            }
        }
    }

    if (advanceSleep) {
        for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
            bool quiet = maxSpeed[static_cast<size_t>(c)] < config.sleepSpeed
                && maxDp[static_cast<size_t>(c)] < config.sleepPressureDelta;
            if (chunkActivity[static_cast<size_t>(c)]) {
                if (quiet) chunkQuietTicks[static_cast<size_t>(c)] = static_cast<uint8_t>(std::min(255, int(chunkQuietTicks[static_cast<size_t>(c)]) + 1));
                else chunkQuietTicks[static_cast<size_t>(c)] = 0;
                if (chunkQuietTicks[static_cast<size_t>(c)] >= config.sleepQuietTicks) chunkActivity[static_cast<size_t>(c)] = 0;
            }
        }
    }
    activeChunks = static_cast<int>(std::count(chunkActivity.begin(), chunkActivity.end(), uint8_t{1}));
    lastWasSleeping = activeChunks == 0;
}

void GasEngine::resetAmbient(FluidEngine &fluid) {
    bool changed = false;
    rebuildVolumes(fluid, changed);
    for (int i = 0; i < GW * GH; ++i) {
        float vol = volume[static_cast<size_t>(i)];
        amount[static_cast<size_t>(i)] = vol * config.ambientPressureAtm;
        writePureGas(i, SUBSTANCE_AIR, amount[static_cast<size_t>(i)]);
        float cap = thermalCapacity(gasMassKg(amount[static_cast<size_t>(i)]), thermalForSubstance(SUBSTANCE_AIR).specificHeat);
        heat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
    }
    std::fill(u.begin(), u.end(), 0.0f);
    std::fill(v.begin(), v.end(), 0.0f);
    escapedAmount = 0.0;
    escapedHeat = 0.0;
    escapedWaterVapor = 0.0;
    recomputePressure();
    commitExpected();
    sleepAll();
    rebuildActivity(false);
}

void GasEngine::handleWorldEdit(FluidEngine &fluid) {
    bool changed = false;
    rebuildVolumes(fluid, changed);
    displaceLiquidOverflow(fluid);
    displaceBlocked(fluid);
    recomputePressure();
    currentAmount = sumAmount();
    amountError = currentAmount - expectedAmount;
    currentWaterVapor = sumWaterVapor();
}

void GasEngine::applyPressureBrush(FluidEngine &fluid, int cx, int cy, int brushRadius,
    float signedAtmPerSec, float dt, BrushShape shape)
{
    float dAtm = signedAtmPerSec * std::max(dt, 1.0f / 30.0f);
    if (std::abs(dAtm) < 1.0e-8f) return;
    double net = 0.0;
    float maxAtm = std::max(config.ambientPressureAtm, config.brushMaxAtm);
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if (!brushContains(shape, cx, cy, x, y, brushRadius)) continue;
            int i = ci(x, y);
            float vol = availableVolume(fluid, x, y);
            volume[static_cast<size_t>(i)] = vol;
            if (vol < GAS_MIN_VOLUME) continue;
            float a = amount[static_cast<size_t>(i)];
            float p = a / vol;
            float np = std::clamp(p + dAtm, 0.0f, maxAtm);
            float na = np * vol;
            float da = na - a;
            if (std::abs(da) <= GAS_MIN_AMOUNT) continue;
            float capOld = (a > GAS_MIN_AMOUNT)
                ? gasMixtureThermalCapacity(gasComponents(i), fluid.config.cellsPerMeter)
                : 0.0f;
            float t = (a > GAS_MIN_AMOUNT && capOld >= MIN_THERMAL_CAPACITY)
                ? tempFromEnergy(heat[static_cast<size_t>(i)], capOld)
                : AMBIENT_TEMPERATURE_K;
            if (da > 0.0f) {
                float capAdd = thermalCapacity(gasMassKg(da), thermalForSubstance(substanceForGasSpecies()).specificHeat);
                heat[static_cast<size_t>(i)] += energyFromTemp(capAdd, t);
            } else if (a > GAS_MIN_AMOUNT) {
                heat[static_cast<size_t>(i)] *= (na / a);
            }
            amount[static_cast<size_t>(i)] = na;
            if (da > 0.0f)
                (void)addGasComponentUntracked(i, SUBSTANCE_AIR, da);
            else if (da < 0.0f && a > GAS_MIN_AMOUNT)
                scaleGasComposition(i, na / a);
            else if (na <= GAS_MIN_AMOUNT)
                clearGasComposition(i);
            clampSpecies(i);
            if (na <= GAS_MIN_AMOUNT) {
                amount[static_cast<size_t>(i)] = 0.0f;
                clearGasComposition(i);
                heat[static_cast<size_t>(i)] = 0.0f;
                na = 0.0f;
            }
            pressure[static_cast<size_t>(i)] = na / vol;
            net += static_cast<double>(na - a);
            wakeChunkAtCell(x, y);
            wakeThermalAt(x, y);
        }
    expectedAmount += net;
    currentAmount += net;
    amountError = currentAmount - expectedAmount;
}

void GasEngine::applyGasBrush(FluidEngine &fluid, int cx, int cy, int brushRadius,
    SubstanceId gasId, float amountPerSec, float dt, BrushShape shape)
{
    if (!validGasComponentId(gasId)) return;
    float dAtm = amountPerSec * std::max(dt, 1.0f / 30.0f);
    if (!(dAtm > 1.0e-8f) || !std::isfinite(dAtm)) return;
    double net = 0.0;
    float maxAtm = std::max(config.ambientPressureAtm, config.brushMaxAtm);
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if (!brushContains(shape, cx, cy, x, y, brushRadius)) continue;
            int i = ci(x, y);
            float vol = availableVolume(fluid, x, y);
            volume[static_cast<size_t>(i)] = vol;
            if (vol < GAS_MIN_VOLUME) continue;
            float a = amount[static_cast<size_t>(i)];
            float p = a / vol;
            float np = std::clamp(p + dAtm, 0.0f, maxAtm);
            float na = np * vol;
            float da = na - a;
            if (!(da > GAS_MIN_AMOUNT) || !std::isfinite(da)) continue;
            GasComponentView cur = gasComponents(i);
            GasComponent add[1] = {{gasId, da}};
            if (!gasPayloadCanMerge(cur.items, cur.count, add, 1)) continue;
            if (addGasComponentUntracked(i, gasId, da) > 0.0f) continue;
            GasComponentView parcel;
            parcel.count = 1;
            parcel.items[0] = {gasId, da};
            float capAdd = gasMixtureThermalCapacity(parcel, fluid.config.cellsPerMeter);
            heat[static_cast<size_t>(i)] += energyFromTemp(capAdd, AMBIENT_TEMPERATURE_K);
            if (!(heat[static_cast<size_t>(i)] >= 0.0f) || !std::isfinite(heat[static_cast<size_t>(i)]))
                heat[static_cast<size_t>(i)] = 0.0f;
            amount[static_cast<size_t>(i)] = na;
            clampSpecies(i);
            if (amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT) {
                amount[static_cast<size_t>(i)] = 0.0f;
                clearGasComposition(i);
                heat[static_cast<size_t>(i)] = 0.0f;
                na = 0.0f;
            }
            pressure[static_cast<size_t>(i)] = (vol > GAS_MIN_VOLUME) ? (na / vol) : 0.0f;
            net += static_cast<double>(na - a);
            wakeAt(x, y);
        }
    expectedAmount += net;
    currentAmount += net;
    amountError = currentAmount - expectedAmount;
}

void GasEngine::eraseAmountBrush(FluidEngine &fluid, int cx, int cy, int brushRadius, BrushShape shape)
{
    (void)fluid;
    double net = 0.0;
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if (!brushContains(shape, cx, cy, x, y, brushRadius)) continue;
            int i = ci(x, y);
            float a = amount[static_cast<size_t>(i)];
            if (a <= GAS_MIN_AMOUNT) continue;
            net -= static_cast<double>(a);
            expectedWaterVapor -= gasComponentAmount(i, SUBSTANCE_WATER);
            amount[static_cast<size_t>(i)] = 0.0f;
            clearGasComposition(i);
            heat[static_cast<size_t>(i)] = 0.0f;
            pressure[static_cast<size_t>(i)] = 0.0f;
            wakeChunkAtCell(x, y);
            wakeThermalAt(x, y);
        }
    expectedAmount += net;
    currentAmount += net;
    amountError = currentAmount - expectedAmount;
}

void GasEngine::simulationTick(FluidEngine &fluid) {
    if (config.simMode == GasSimMode::Off) {
        lastWasSleeping = true;
        lastSubsteps = 0;
        lastStepMs = 0.0;
        timing.physics = 0.0;
        return;
    }
    if (config.simMode == GasSimMode::Half && (fluid.tickNo & 1u) != 0u) {
        lastStepMs = 0.0;
        return;
    }
    auto start = Clock::now();
    bool volumeChanged = false;
    rebuildVolumes(fluid, volumeChanged, false);
    timing.occupancy = std::chrono::duration<double, std::milli>(Clock::now() - start).count();

    auto d0 = Clock::now();
    displaceLiquidOverflow(fluid);
    displaceBlocked(fluid);
    recomputePressure();
    timing.displace = std::chrono::duration<double, std::milli>(Clock::now() - d0).count();

    bool anyActive = false;
    for (uint8_t a : chunkActivity) if (a) { anyActive = true; break; }
    if (!anyActive && !volumeChanged) {
        currentAmount = sumAmount();
        amountError = currentAmount - expectedAmount;
        activeChunks = 0;
        activeGasCells = 0;
        lastSubsteps = 0;
        lastWasSleeping = true;
        lastStepMs = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        timing.physics = lastStepMs;
        return;
    }

    float maxSpeed = 0.0f;
    for (float s : u) maxSpeed = std::max(maxSpeed, std::abs(s));
    for (float s : v) maxSpeed = std::max(maxSpeed, std::abs(s));
    int sub = std::clamp(static_cast<int>(std::ceil(maxSpeed * PHYSICS_DT / config.maxTravelPerSubstep)), 1, config.maxSubsteps);
    lastSubsteps = sub;
    float dt = PHYSICS_DT / static_cast<float>(sub);

    auto t0 = Clock::now();
    for (int step = 0; step < sub; ++step) {
        integrateVelocity(fluid, dt);
        recomputePressure();
        accumulateTransfer(fluid, dt);
        applyFluxes();
        applyBoundaryFlux(dt);
        recomputePressure();
    }
    timing.transfer = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

    auto c0 = Clock::now();
    rebuildActivity(true);
    timing.chunks = std::chrono::duration<double, std::milli>(Clock::now() - c0).count();
    lastStepMs = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    timing.physics = lastStepMs;
}

void GasEngine::applyPressureForces(RigidBodyEngine &rigid, FluidEngine const &fluid) const {
    if (config.simMode == GasSimMode::Off) return;
    constexpr int ndx[4] = {-1, 1, 0, 0};
    constexpr int ndy[4] = {0, 0, -1, 1};
    float gravity = fluid.gridGravity();
    float ambient = config.ambientPressureAtm;
    for (size_t bi = 0; bi < rigid.bodies.size(); ++bi) {
        RigidBody &b = rigid.bodies[bi];
        if (b.anchored) continue;
        int x0 = std::max(1, static_cast<int>(std::floor(b.aabbX0)) - 1);
        int y0 = std::max(1, static_cast<int>(std::floor(b.aabbY0)) - 1);
        int x1 = std::min(GW - 2, static_cast<int>(std::ceil(b.aabbX1)) + 1);
        int y1 = std::min(GH - 2, static_cast<int>(std::ceil(b.aabbY1)) + 1);
        float gfx = 0.0f, gfy = 0.0f, gtq = 0.0f;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            if (rigid.occupant[static_cast<size_t>(ci(x, y))] != static_cast<int>(bi)) continue;
            for (int n = 0; n < 4; ++n) {
                int nx = x + ndx[n], ny = y + ndy[n];
                if (!FluidEngine::inside(nx, ny)) continue;
                int ni = ci(nx, ny);
                if (volume[static_cast<size_t>(ni)] < 0.25f) continue;
                if (fluid.fill[static_cast<size_t>(ni)] >= MIN_ACTIVE_FILL) continue;
                float gauge = pressure[static_cast<size_t>(ni)] - ambient;
                gauge = std::clamp(gauge, -2.0f, 2.0f);
                if (std::abs(gauge) < 0.02f) continue;
                float p = gauge * config.rigidForceScale;
                float fx = -p * static_cast<float>(ndx[n]);
                float fy = -p * static_cast<float>(ndy[n]);
                float rx = (x + 0.5f) - b.x, ry = (y + 0.5f) - b.y;
                gfx += fx; gfy += fy;
                gtq += rx * fy - ry * fx;
            }
        }
        b.fx += gfx; b.fy += gfy; b.torque += gtq;
        float accel = std::sqrt(b.fx * b.fx + b.fy * b.fy) * b.invMass;
        float aMax = 3.0f * gravity;
        if (accel > aMax && accel > 1e-5f) {
            float s = aMax / accel;
            b.fx *= s; b.fy *= s; b.torque *= s;
            accel = aMax;
        }
        if (b.dormant) {
            if (accel > 2.5f * gravity) {
                b.dormant = false;
                b.sleeping = false;
                b.quietTicks = 0;
            } else {
                b.fx = b.fy = b.torque = 0.0f;
            }
            continue;
        }
        if (b.sleeping && accel > 2.0f * gravity) { b.sleeping = false; b.quietTicks = 0; }
    }
}

void GasEngine::wallRect(FluidEngine &fluid, int x0, int y0, int x1, int y1) {
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
        if (FluidEngine::inside(x, y)) fluid.solid[static_cast<size_t>(ci(x, y))] = 1;
}

void GasEngine::vacuumAll() {
    std::fill(amount.begin(), amount.end(), 0.0f);
    std::fill(gasCompId.begin(), gasCompId.end(), SUBSTANCE_NONE);
    std::fill(gasCompAmt.begin(), gasCompAmt.end(), 0.0f);
    std::fill(gasCompCount.begin(), gasCompCount.end(), uint8_t{0});
    std::fill(heat.begin(), heat.end(), 0.0f);
    std::fill(u.begin(), u.end(), 0.0f);
    std::fill(v.begin(), v.end(), 0.0f);
    escapedAmount = 0.0;
    escapedHeat = 0.0;
    escapedWaterVapor = 0.0;
}

void GasEngine::fillRectAmount(FluidEngine const &fluid, int x0, int y0, int x1, int y1, float atm) {
    (void)fluid;
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        if (!FluidEngine::inside(x, y)) continue;
        int i = ci(x, y);
        float vol = volume[static_cast<size_t>(i)];
        if (vol < GAS_MIN_VOLUME) continue;
        amount[static_cast<size_t>(i)] = vol * atm;
        writePureGas(i, SUBSTANCE_AIR, amount[static_cast<size_t>(i)]);
        float cap = thermalCapacity(gasMassKg(amount[static_cast<size_t>(i)]),
            thermalForSubstance(SUBSTANCE_AIR).specificHeat);
        heat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
    }
}

void GasEngine::loadTestScene(FluidEngine &fluid, RigidBodyEngine &rigid, int scene) {
    rigid.clear();
    fluid.clearWorld();
    config.boundary = GasBoundary::Sealed;
    bool changed = false;
    rebuildVolumes(fluid, changed);
    vacuumAll();

    auto box = [&](int x0, int y0, int x1, int y1) {
        wallRect(fluid, x0, y0, x1, y0);
        wallRect(fluid, x0, y1, x1, y1);
        wallRect(fluid, x0, y0, x0, y1);
        wallRect(fluid, x1, y0, x1, y1);
    };

    if (scene == 1) {
        // A: uniform ambient air, open edges at 1 atm
        config.boundary = GasBoundary::OpenAmbient;
        rebuildVolumes(fluid, changed);
        resetAmbient(fluid);
        return;
    }
    if (scene == 2) {
        // B: sealed chamber of 1 atm in vacuum
        box(60, 30, 140, 90);
        rebuildVolumes(fluid, changed);
        fillRectAmount(fluid, 61, 31, 139, 89, 1.0f);
    } else if (scene == 3) {
        // C: high / low chambers with a closed wall (open one cell after load in diagnostics)
        box(40, 25, 160, 95);
        wallRect(fluid, 100, 26, 100, 94);
        rebuildVolumes(fluid, changed);
        fillRectAmount(fluid, 41, 26, 99, 94, 2.0f);
        fillRectAmount(fluid, 101, 26, 159, 94, 0.25f);
    } else if (scene == 4) {
        // D: solid block; cavity is created by diagnostics / player erase
        wallRect(fluid, 70, 35, 130, 85);
        rebuildVolumes(fluid, changed);
    } else if (scene == 5) {
        // E: sealed chamber ready to compress
        box(70, 40, 130, 90);
        rebuildVolumes(fluid, changed);
        fillRectAmount(fluid, 71, 41, 129, 89, 1.0f);
    } else if (scene == 6) {
        // F: sealed chamber with extra interior solid that can be erased to expand
        box(70, 40, 130, 90);
        wallRect(fluid, 71, 41, 90, 89);
        rebuildVolumes(fluid, changed);
        fillRectAmount(fluid, 91, 41, 129, 89, 1.0f);
    } else if (scene == 7) {
        // G: sealed air plus a falling body (rigid scene clears the world first)
        config.boundary = GasBoundary::Sealed;
        rigid.loadTestScene(fluid, 1);
        resetAmbient(fluid);
        config.boundary = GasBoundary::Sealed;
        handleWorldEdit(fluid);
        wakeAll();
        commitExpected();
        return;
    } else {
        // H: sealed uniform air for long conservation
        config.boundary = GasBoundary::Sealed;
        rebuildVolumes(fluid, changed);
        resetAmbient(fluid);
        config.boundary = GasBoundary::Sealed;
        wakeAll();
        commitExpected();
        return;
    }

    recomputePressure();
    std::fill(u.begin(), u.end(), 0.0f);
    std::fill(v.begin(), v.end(), 0.0f);
    escapedAmount = 0.0;
    escapedHeat = 0.0;
    commitExpected();
    wakeAll();
}

void GasEngine::runDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid) {
    std::ofstream out(miscFile("gas_diag.tsv"), std::ios::trunc);
    out << "case\ttick\tamount\texpected\terror\tescaped\tactive_chunks\tactive_cells\tmean_p\tmin_p\tmax_p\tms\tpass\n";

    auto stats = [&](float &meanP, float &minP, float &maxP, int &n) {
        double sum = 0.0;
        minP = 1.0e9f; maxP = -1.0e9f; n = 0;
        for (int i = 0; i < GW * GH; ++i) {
            if (volume[static_cast<size_t>(i)] < GAS_MIN_VOLUME) continue;
            float p = pressure[static_cast<size_t>(i)];
            sum += p; minP = std::min(minP, p); maxP = std::max(maxP, p); ++n;
        }
        meanP = n ? static_cast<float>(sum / n) : 0.0f;
        if (n == 0) { minP = 0.0f; maxP = 0.0f; }
    };

    auto emit = [&](char const *name, int tick, bool pass) {
        float meanP, minP, maxP; int n;
        stats(meanP, minP, maxP, n);
        out << name << '\t' << tick << '\t' << currentAmount << '\t' << expectedAmount << '\t'
            << amountError << '\t' << escapedAmount << '\t' << activeChunks << '\t' << activeGasCells << '\t'
            << meanP << '\t' << minP << '\t' << maxP << '\t' << lastStepMs << '\t' << (pass ? 1 : 0) << '\n';
    };

    auto regionAmount = [&](int x0, int y0, int x1, int y1) {
        double s = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            if (FluidEngine::inside(x, y)) s += amount[static_cast<size_t>(ci(x, y))];
        return s;
    };

    auto tickN = [&](int n) {
        for (int i = 0; i < n; ++i) {
            rigid.step(fluid, PHYSICS_DT);
            fluid.simulationTick();
            simulationTick(fluid);
            rigid.gatherFluidForces(fluid);
            applyPressureForces(rigid, fluid);
        }
    };

    // A uniform
    loadTestScene(fluid, rigid, 1);
    tickN(90);
    {
        float meanP, minP, maxP; int n; stats(meanP, minP, maxP, n);
        bool pass = std::abs(amountError) < 1.0e-3 && (maxP - minP) < 0.02f && activeChunks == 0;
        emit("A_uniform", 90, pass);
    }

    // B sealed chamber
    loadTestScene(fluid, rigid, 2);
    double b0 = currentAmount;
    tickN(90);
    {
        float pIn = pressure[static_cast<size_t>(ci(100, 60))];
        float pOut = pressure[static_cast<size_t>(ci(30, 20))];
        bool pass = std::abs(currentAmount - b0) < 1.0e-3 && pIn > 0.9f && pIn < 1.1f && pOut < 0.05f;
        emit("B_sealed", 90, pass);
    }

    // C open a hole
    loadTestScene(fluid, rigid, 3);
    double c0 = currentAmount;
    double left0a = regionAmount(41, 26, 99, 94);
    double right0a = regionAmount(101, 26, 159, 94);
    fluid.solid[static_cast<size_t>(ci(100, 60))] = 0;
    handleWorldEdit(fluid);
    tickN(150);
    {
        double left1a = regionAmount(41, 26, 99, 94);
        double right1a = regionAmount(101, 26, 159, 94);
        bool pass = std::abs(currentAmount - c0) < 1.0e-2 && left1a < left0a - 8.0 && right1a > right0a + 8.0;
        emit("C_open_hole", 150, pass);
    }

    // D erase cavity in solid
    loadTestScene(fluid, rigid, 4);
    for (int y = 40; y <= 80; ++y) for (int x = 80; x <= 120; ++x)
        fluid.solid[static_cast<size_t>(ci(x, y))] = 0;
    handleWorldEdit(fluid);
    double cavity = 0.0;
    for (int y = 40; y <= 80; ++y) for (int x = 80; x <= 120; ++x)
        cavity += amount[static_cast<size_t>(ci(x, y))];
    tickN(30);
    double cavityLater = 0.0;
    for (int y = 40; y <= 80; ++y) for (int x = 80; x <= 120; ++x)
        cavityLater += amount[static_cast<size_t>(ci(x, y))];
    emit("D_erased_cavity", 30, cavity < 0.05 && cavityLater < 0.05);

    // E compress sealed chamber by adding a solid slab
    loadTestScene(fluid, rigid, 5);
    double eAmt = currentAmount;
    wallRect(fluid, 71, 41, 129, 50);
    handleWorldEdit(fluid);
    tickN(24);
    {
        float maxInside = 0.0f;
        for (int y = 51; y <= 89; ++y) for (int x = 71; x <= 129; ++x)
            maxInside = std::max(maxInside, pressure[static_cast<size_t>(ci(x, y))]);
        bool pass = std::abs(currentAmount - eAmt) < 1.0e-2 && maxInside > 1.15f;
        emit("E_compress", 24, pass);
    }

    // F expand by erasing interior solid
    loadTestScene(fluid, rigid, 6);
    double fAmt = currentAmount;
    for (int y = 41; y <= 89; ++y) for (int x = 71; x <= 90; ++x)
        fluid.solid[static_cast<size_t>(ci(x, y))] = 0;
    handleWorldEdit(fluid);
    tickN(60);
    {
        double entered = regionAmount(71, 41, 90, 89);
        bool pass = std::abs(currentAmount - fAmt) < 1.0e-2 && entered > 80.0;
        emit("F_expand", 60, pass);
    }

    // G rigid displacement conservation
    loadTestScene(fluid, rigid, 7);
    double g0 = currentAmount;
    tickN(90);
    emit("G_rigid_displace", 90, std::abs(currentAmount - g0) < 0.05);

    // H long sealed conservation
    loadTestScene(fluid, rigid, 8);
    double h0 = currentAmount;
    tickN(240);
    emit("H_conservation", 240, std::abs(currentAmount - h0) < 1.0e-2 && std::abs(amountError) < 1.0e-2);
}

void GasEngine::runCompositionSanityCheck(FluidEngine &fluid, RigidBodyEngine &rigid) {
    std::ofstream out(miscFile("gas_composition_sanity.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };
    auto near = [](float a, float b) {
        return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 1.0e-4f;
    };

    GasComponent full[kMaxGasComponents] = {
        {SUBSTANCE_AIR, 0.25f}, {SUBSTANCE_AIR, 0.25f},
        {SUBSTANCE_AIR, 0.25f}, {SUBSTANCE_AIR, 0.25f}
    };
    GasComponent waterSrc[1] = {{SUBSTANCE_WATER, 0.5f}};
    emit("merge_reject_new_species_into_full_slots",
        !gasPayloadCanMerge(full, kMaxGasComponents, waterSrc, 1), "");
    emit("merge_accept_existing_species_into_full_slots",
        gasPayloadCanMerge(full, kMaxGasComponents, full, 1), "");
    GasComponent airOnly[1] = {{SUBSTANCE_AIR, 1.0f}};
    emit("merge_accept_water_into_air",
        gasPayloadCanMerge(airOnly, 1, waterSrc, 1), "");

    rigid.clear();
    fluid.clearWorld();
    config.boundary = GasBoundary::Sealed;
    bool changed = false;
    rebuildVolumes(fluid, changed);
    vacuumAll();
    int ia = ci(10, 10);
    int ib = ci(11, 10);
    volume[static_cast<size_t>(ia)] = 1.0f;
    volume[static_cast<size_t>(ib)] = 1.0f;
    writePureGas(ia, SUBSTANCE_WATER, 1.0f);
    writePureGas(ib, SUBSTANCE_AIR, 1.0f);
    int base = compositionSlot(ib, 0);
    for (int s = 0; s < kMaxGasComponents; ++s) {
        gasCompId[static_cast<size_t>(base + s)] = SUBSTANCE_AIR;
        gasCompAmt[static_cast<size_t>(base + s)] = 0.25f;
    }
    gasCompCount[static_cast<size_t>(ib)] = static_cast<uint8_t>(kMaxGasComponents);
    amount[static_cast<size_t>(ib)] = 1.0f;
    heat[static_cast<size_t>(ia)] = 1.0f;
    heat[static_cast<size_t>(ib)] = 1.0f;
    std::fill(fluxU.begin(), fluxU.end(), 0.0f);
    std::fill(fluxV.begin(), fluxV.end(), 0.0f);
    fluxU[static_cast<size_t>(ui(11, 10))] = 0.25f;
    applyFluxes();
    emit("overflow_face_rejects_amount_and_species",
        near(amount[static_cast<size_t>(ia)], 1.0f)
            && near(amount[static_cast<size_t>(ib)], 1.0f)
            && near(gasComponentAmount(ia, SUBSTANCE_WATER), 1.0f)
            && gasComponentAmount(ib, SUBSTANCE_WATER) <= GAS_MIN_AMOUNT
            && gasCompositionValid(ia),
        "a=" + std::to_string(amount[static_cast<size_t>(ia)])
            + " b=" + std::to_string(amount[static_cast<size_t>(ib)])
            + " bWater=" + std::to_string(gasComponentAmount(ib, SUBSTANCE_WATER)));

    writePureGas(ia, SUBSTANCE_WATER, 1.0f);
    writePureGas(ib, SUBSTANCE_AIR, 1.0f);
    heat[static_cast<size_t>(ia)] = 1.0f;
    heat[static_cast<size_t>(ib)] = 1.0f;
    std::fill(fluxU.begin(), fluxU.end(), 0.0f);
    std::fill(fluxV.begin(), fluxV.end(), 0.0f);
    fluxU[static_cast<size_t>(ui(11, 10))] = 0.25f;
    applyFluxes();
    emit("compatible_face_transfers_composition",
        near(amount[static_cast<size_t>(ia)] + amount[static_cast<size_t>(ib)], 2.0f)
            && gasCompositionValid(ia) && gasCompositionValid(ib)
            && gasComponentAmount(ib, SUBSTANCE_WATER) > GAS_MIN_AMOUNT,
        "a=" + std::to_string(amount[static_cast<size_t>(ia)])
            + " b=" + std::to_string(amount[static_cast<size_t>(ib)])
            + " bWater=" + std::to_string(gasComponentAmount(ib, SUBSTANCE_WATER)));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
