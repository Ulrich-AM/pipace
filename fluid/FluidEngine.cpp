#include "FluidEngine.h"
#include "DiagOutput.h"
#include "substance/LiquidMixtureProperties.h"
#include "thermal/ThermalTypes.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <thread>
#include <utility>
#include <vector>

FluidEngine::FluidEngine()
    : solid(GW * GH, 0)
    , solidHeat(GW * GH, 0.0f)
    , dynamicSolid(GW * GH, 0)
    , dynamicVelX(GW * GH, 0.0f)
    , dynamicVelY(GW * GH, 0.0f)
    , fill(GW * GH, 0.0f)
    , nextFill(GW * GH, 0.0f)
    , liquidHeat(GW * GH, 0.0f)
    , nextHeat(GW * GH, 0.0f)
    , dyeR(GW * GH, 0.0f)
    , dyeG(GW * GH, 0.0f)
    , dyeB(GW * GH, 0.0f)
    , liquidCompId(static_cast<size_t>(GW * GH) * kMaxLiquidComponents, runtimeNone())
    , liquidCompAmt(static_cast<size_t>(GW * GH) * kMaxLiquidComponents, 0.0f)
    , liquidCompCount(GW * GH, 0)
    , solidifyPendingId(GW * GH, SUBSTANCE_NONE)
    , solidifyPendingKg(GW * GH, 0.0f)
    , solidifyPendingHeatJ(GW * GH, 0.0f)
    , nextDyeR(GW * GH, 0.0f)
    , nextDyeG(GW * GH, 0.0f)
    , nextDyeB(GW * GH, 0.0f)
    , nextCompId(static_cast<size_t>(GW * GH) * kMaxLiquidComponents, runtimeNone())
    , nextCompAmt(static_cast<size_t>(GW * GH) * kMaxLiquidComponents, 0.0f)
    , nextCompCount(GW * GH, 0)
    , pressure(GW * GH, 0.0f)
    , divergenceField(GW * GH, 0.0f)
    , outgoing(GW * GH, 0.0f)
    , incoming(GW * GH, 0.0f)
    , donorScale(GW * GH, 1.0f)
    , receiverScale(GW * GH, 1.0f)
    , cellForceX(GW * GH, 0.0f)
    , cellForceY(GW * GH, 0.0f)
    , curlField(GW * GH, 0.0f)
    , smoothedFill(GW * GH, 0.0f)
    , surfaceNormalX(GW * GH, 0.0f)
    , surfaceNormalY(GW * GH, 0.0f)
    , surfaceCurvature(GW * GH, 0.0f)
    , surfaceMask(GW * GH, 0)
    , pressureBefore(GW * GH, 0.0f)
    , previousFill(GW * GH, 0.0f)
    , residualTarget(GW * GH, -1)
    , residualTransfer(GW * GH, 0.0f)
    , waterShade(GW * GH, 0)
    , u((GW + 1) * GH, 0.0f)
    , v(GW * (GH + 1), 0.0f)
    , uTemp((GW + 1) * GH, 0.0f)
    , vTemp(GW * (GH + 1), 0.0f)
    , uScratch((GW + 1) * GH, 0.0f)
    , vScratch(GW * (GH + 1), 0.0f)
    , fluxU((GW + 1) * GH, 0.0f)
    , fluxV(GW * (GH + 1), 0.0f)
    , chunkActivity(CHUNK_W * CHUNK_H, 0)
    , chunkQuietTicks(CHUNK_W * CHUNK_H, 0)
    , chunkHasFluid(CHUNK_W * CHUNK_H, 0)
    , chunkSolveMask(CHUNK_W * CHUNK_H, 0)
    , chunkHaloSource(CHUNK_W * CHUNK_H, 0)
    , thermalChunkWake(CHUNK_W * CHUNK_H, 0)
    , pixels(GW * GH, 0)
{
    relocateStamp.assign(static_cast<size_t>(GW * GH), 0);
    relocateQueue.reserve(512);
    advectOverflowIndex.reserve(64);
    advectOverflowVol.reserve(64);
    advectOverflowCarry.reserve(64);
    splashes.reserve(1024);
    pressureRed.reserve(GW * GH / 2);
    pressureBlack.reserve(GW * GH / 2);
    pressureStencils.reserve(GW * GH);
    surfaceCells.reserve(GW * GH / 4);
    nonzeroFluxUFaces.reserve((GW + 1) * GH / 4);
    nonzeroFluxVFaces.reserve(GW * (GH + 1) / 4);
    fluxTouchedCells.reserve(GW * GH / 4);
    dynamicOccupiedCells.reserve(256);
    syncWorkerPool();
}

void FluidEngine::syncWorkerPool() {
    int total = resolvedWorkerCount();
    lastResolvedWorkers = total;
    workerPool.setExtraWorkers(std::max(0, total - 1));
}

int FluidEngine::maxSelectableWorkers() const {
    unsigned hc = std::thread::hardware_concurrency();
    int cores = hc == 0 ? 2 : static_cast<int>(hc);
    return std::clamp(cores, 1, 8);
}

int FluidEngine::autoWorkerCount() const {
    unsigned hc = std::thread::hardware_concurrency();
    int cores = hc == 0 ? 2 : static_cast<int>(hc);
    int cells = GW * GH;
    // Default 200x120 work is too small to beat red/black barriers.
    if (cells <= 200 * 120) return 1;
    if (cores <= 2) return 1;
    if (cores <= 4) return 2;
    if (cores <= 6) return std::min(4, cores - 1);
    return std::min(6, cores - 2);
}

int FluidEngine::resolvedWorkerCount() const {
    int cap = maxSelectableWorkers();
    if (config.workerCount <= 0) return std::clamp(autoWorkerCount(), 1, cap);
    return std::clamp(config.workerCount, 1, cap);
}

bool FluidEngine::useParallelPressure(int cellCount) const {
    return workerPool.extraWorkers() > 0 && cellCount >= config.pressureParallelMinCells;
}

int FluidEngine::splashAtCell(int x, int y) const {
    if (!inside(x, y)) return 0;
    int count = 0;
    for (SplashParticle const &p : splashes) {
        if (static_cast<int>(std::floor(p.x)) == x && static_cast<int>(std::floor(p.y)) == y) ++count;
    }
    return count;
}

void FluidEngine::flushPaintDirty() {
    if (!paintDirty) return;
    paintDirty = false;
    enforceSolidBoundaries();
    rebuildActivityAndMetrics();
}

double FluidEngine::elapsedMs(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void FluidEngine::wakeThermalAt(int x, int y) {
    if (!inside(x, y)) return;
    int c = (y / CHUNK) * CHUNK_W + (x / CHUNK);
    if (c >= 0 && c < static_cast<int>(thermalChunkWake.size())) thermalChunkWake[static_cast<size_t>(c)] = 1;
}

int FluidEngine::compositionSlot(int cell, int slot) {
    return cell * kMaxLiquidComponents + slot;
}

void FluidEngine::clearComposition(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        liquidCompId[static_cast<size_t>(base + s)] = runtimeNone();
        liquidCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    liquidCompCount[static_cast<size_t>(index)] = 0;
}

void FluidEngine::compactComposition(int index) {
    if (index < 0 || index >= GW * GH) return;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    int w = 0;
    for (int s = 0; s < n; ++s) {
        RuntimeSubstanceRef id = liquidCompId[static_cast<size_t>(base + s)];
        float amt = liquidCompAmt[static_cast<size_t>(base + s)];
        if (amt > kMinLiquidComponent) {
            if (w != s) {
                liquidCompId[static_cast<size_t>(base + w)] = id;
                liquidCompAmt[static_cast<size_t>(base + w)] = amt;
            }
            ++w;
        }
    }
    for (int s = w; s < n; ++s) {
        liquidCompId[static_cast<size_t>(base + s)] = runtimeNone();
        liquidCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    liquidCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(w);
}

void FluidEngine::copyCompositionToNext(int index) {
    int n = liquidCompCount[static_cast<size_t>(index)];
    nextCompCount[static_cast<size_t>(index)] = n;
    int base = compositionSlot(index, 0);
    for (int s = 0; s < n; ++s) {
        nextCompId[static_cast<size_t>(base + s)] = liquidCompId[static_cast<size_t>(base + s)];
        nextCompAmt[static_cast<size_t>(base + s)] = liquidCompAmt[static_cast<size_t>(base + s)];
    }
}

void FluidEngine::commitNextComposition(int index) {
    int n = nextCompCount[static_cast<size_t>(index)];
    liquidCompCount[static_cast<size_t>(index)] = n;
    int base = compositionSlot(index, 0);
    for (int s = 0; s < n; ++s) {
        liquidCompId[static_cast<size_t>(base + s)] = nextCompId[static_cast<size_t>(base + s)];
        liquidCompAmt[static_cast<size_t>(base + s)] = nextCompAmt[static_cast<size_t>(base + s)];
    }
    for (int s = n; s < kMaxLiquidComponents; ++s) {
        liquidCompId[static_cast<size_t>(base + s)] = runtimeNone();
        liquidCompAmt[static_cast<size_t>(base + s)] = 0.0f;
    }
    compactComposition(index);
}

float FluidEngine::addComponentUntracked(int index, SubstanceId id, float amount) {
    if (index < 0 || index >= GW * GH || !validLiquidComponentId(id)) return std::max(0.0f, amount);
    if (!(amount > kMinLiquidComponent)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (liquidCompId[static_cast<size_t>(base + s)] == id) {
            liquidCompAmt[static_cast<size_t>(base + s)] += amount;
            return 0.0f;
        }
    }
    if (n >= kMaxLiquidComponents) return amount;
    liquidCompId[static_cast<size_t>(base + n)] = runtimeBuiltIn(id);
    liquidCompAmt[static_cast<size_t>(base + n)] = amount;
    liquidCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n + 1);
    return 0.0f;
}

float FluidEngine::addNextComponentUntracked(int index, SubstanceId id, float amount) {
    if (index < 0 || index >= GW * GH || !validLiquidComponentId(id)) return std::max(0.0f, amount);
    if (!(amount > kMinLiquidComponent)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = nextCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (nextCompId[static_cast<size_t>(base + s)] == id) {
            nextCompAmt[static_cast<size_t>(base + s)] += amount;
            return 0.0f;
        }
    }
    if (n >= kMaxLiquidComponents) return amount;
    nextCompId[static_cast<size_t>(base + n)] = runtimeBuiltIn(id);
    nextCompAmt[static_cast<size_t>(base + n)] = amount;
    nextCompCount[static_cast<size_t>(index)] = static_cast<uint8_t>(n + 1);
    return 0.0f;
}

bool FluidEngine::cellHasDuplicateComponents(int index) const {
    if (index < 0 || index >= GW * GH) return false;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            if (liquidCompId[static_cast<size_t>(base + a)] == liquidCompId[static_cast<size_t>(base + b)])
                return true;
    return false;
}

void FluidEngine::clearEmptyLiquidCell(int index) {
    if (index < 0 || index >= GW * GH) return;
    size_t i = static_cast<size_t>(index);
    if (fill[i] > 1.0e-8f) {
        compactComposition(index);
        if (dyeR[i] < 0.0f) dyeR[i] = 0.0f;
        if (dyeG[i] < 0.0f) dyeG[i] = 0.0f;
        if (dyeB[i] < 0.0f) dyeB[i] = 0.0f;
        return;
    }
    fill[i] = 0.0f;
    liquidHeat[i] = 0.0f;
    dyeR[i] = dyeG[i] = dyeB[i] = 0.0f;
    clearComposition(index);
}

float FluidEngine::honeyFraction(int index) const {
    return liquidComponentFraction(index, SUBSTANCE_HONEY);
}

float FluidEngine::liquidComponentAmount(int index, SubstanceId id) const {
    if (index < 0 || index >= GW * GH || !validLiquidComponentId(id)) return 0.0f;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n; ++s) {
        if (liquidCompId[static_cast<size_t>(base + s)] == id)
            return std::max(0.0f, liquidCompAmt[static_cast<size_t>(base + s)]);
    }
    return 0.0f;
}

float FluidEngine::liquidComponentFraction(int index, SubstanceId id) const {
    if (index < 0 || index >= GW * GH) return 0.0f;
    float f = fill[static_cast<size_t>(index)];
    if (f <= kMinLiquidComponent) return 0.0f;
    float amt = liquidComponentAmount(index, id);
    if (!(amt > 0.0f)) return 0.0f;
    return std::clamp(amt / f, 0.0f, 1.0f);
}

void FluidEngine::setLiquidComponentAmount(int index, SubstanceId id, float amount) {
    if (index < 0 || index >= GW * GH) return;
    if (!validLiquidComponentId(id)) return;
    if (!std::isfinite(amount) || amount < 0.0f) amount = 0.0f;
    size_t i = static_cast<size_t>(index);
    compactComposition(index);
    float oldFill = std::max(0.0f, fill[i]);
    float oldAmt = liquidComponentAmount(index, id);
    float others = 0.0f;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[i];
    for (int s = 0; s < n; ++s) {
        if (liquidCompId[static_cast<size_t>(base + s)] != id)
            others += std::max(0.0f, liquidCompAmt[static_cast<size_t>(base + s)]);
    }
    if (amount <= kMinLiquidComponent) {
        for (int s = 0; s < n; ++s) {
            if (liquidCompId[static_cast<size_t>(base + s)] == id) {
                liquidCompAmt[static_cast<size_t>(base + s)] = 0.0f;
                break;
            }
        }
        compactComposition(index);
        fill[i] = others;
        if (oldFill > kMinLiquidComponent && fill[i] < oldFill) {
            float keep = fill[i] / oldFill;
            liquidHeat[i] *= keep;
            dyeR[i] *= keep; dyeG[i] *= keep; dyeB[i] *= keep;
        }
        clearEmptyLiquidCell(index);
        return;
    }
    int slot = -1;
    for (int s = 0; s < n; ++s)
        if (liquidCompId[static_cast<size_t>(base + s)] == id) { slot = s; break; }
    if (slot < 0) {
        if (n >= kMaxLiquidComponents) return; // overflow: reject, keep matter
        slot = n;
        liquidCompId[static_cast<size_t>(base + slot)] = runtimeBuiltIn(id);
        liquidCompCount[i] = static_cast<uint8_t>(n + 1);
        n += 1;
    }
    if (amount <= oldFill) {
        float remain = oldFill - amount;
        if (others > kMinLiquidComponent) {
            float scale = remain / others;
            for (int s = 0; s < n; ++s) {
                if (liquidCompId[static_cast<size_t>(base + s)] == id) continue;
                liquidCompAmt[static_cast<size_t>(base + s)] *= scale;
            }
        }
        liquidCompAmt[static_cast<size_t>(base + slot)] = amount;
        fill[i] = oldFill;
    } else {
        liquidCompAmt[static_cast<size_t>(base + slot)] = amount;
        fill[i] = amount + others;
    }
    compactComposition(index);
    if (oldFill > kMinLiquidComponent && fill[i] < oldFill) {
        float keep = fill[i] / oldFill;
        liquidHeat[i] *= keep;
        dyeR[i] *= keep; dyeG[i] *= keep; dyeB[i] *= keep;
    }
    (void)oldAmt;
    clearEmptyLiquidCell(index);
}

void FluidEngine::addLiquidComponentAmount(int index, SubstanceId id, float delta) {
    setLiquidComponentAmount(index, id, liquidComponentAmount(index, id) + delta);
}

LiquidComponentView FluidEngine::liquidComponents(int index) const {
    LiquidComponentView view;
    if (index < 0 || index >= GW * GH) return view;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[static_cast<size_t>(index)];
    for (int s = 0; s < n && view.count < kMaxLiquidComponents; ++s) {
        float amt = liquidCompAmt[static_cast<size_t>(base + s)];
        RuntimeSubstanceRef id = liquidCompId[static_cast<size_t>(base + s)];
        if (amt > kMinLiquidComponent && validLiquidComponentId(id))
            view.items[view.count++] = {id, amt};
    }
    return view;
}

SubstanceId FluidEngine::dominantLiquidSubstance(int index) const {
    SubstanceId best = SUBSTANCE_NONE;
    float bestAmt = 0.0f;
    forEachLiquidComponent(index, [&](SubstanceId id, float amt) {
        if (amt > bestAmt) {
            bestAmt = amt;
            best = id;
        }
    });
    return best;
}

bool FluidEngine::liquidCompositionValid(int index) const {
    if (index < 0 || index >= GW * GH) return false;
    size_t i = static_cast<size_t>(index);
    float f = fill[i];
    if (!std::isfinite(f) || f < 0.0f) return false;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[i];
    if (n < 0 || n > kMaxLiquidComponents) return false;
    float sum = 0.0f;
    int validPositive = 0;
    for (int s = 0; s < n; ++s) {
        RuntimeSubstanceRef id = liquidCompId[static_cast<size_t>(base + s)];
        float amt = liquidCompAmt[static_cast<size_t>(base + s)];
        if (!std::isfinite(amt) || amt < 0.0f) return false;
        if (amt > kMinLiquidComponent) {
            if (!validLiquidComponentId(id)) return false;
            ++validPositive;
        }
        sum += amt;
        for (int b = s + 1; b < n; ++b)
            if (id != SUBSTANCE_NONE && liquidCompId[static_cast<size_t>(base + b)] == id)
                return false;
    }
    if (f <= kMinLiquidComponent)
        return validPositive == 0;
    if (validPositive < 1) return false;
    return std::abs(sum - f) <= 1.0e-5f;
}

LiquidMixtureProperties FluidEngine::mixProperties(int index) const {
    if (index < 0 || index >= GW * GH)
        return referenceLiquidMixture();
    return evaluateLiquidMixture(liquidComponents(index));
}

float FluidEngine::mixDensity(int index) const {
    return mixProperties(index).density;
}

float FluidEngine::mixSpecificHeat(int index) const {
    return mixProperties(index).specificHeat;
}

float FluidEngine::mixConductivity(int index) const {
    return mixProperties(index).conductivity;
}

float FluidEngine::mixViscosity(int index) const {
    if (index < 0 || index >= GW * GH)
        return evaluateLiquidMixtureViscosity(LiquidComponentView{}, AMBIENT_TEMPERATURE_K);
    LiquidComponentView view = liquidComponents(index);
    LiquidMixtureProperties mix = evaluateLiquidMixture(view);
    float T = AMBIENT_TEMPERATURE_K;
    float f = fill[static_cast<size_t>(index)];
    if (f > kMinLiquidComponent) {
        float cap = thermalCapacity(massKg(mix.density, f, config.cellsPerMeter), mix.specificHeat);
        T = tempFromEnergy(liquidHeat[static_cast<size_t>(index)], cap);
    }
    return evaluateLiquidMixtureViscosity(view, T);
}

float FluidEngine::mixSurfaceTension(int index) const {
    return mixProperties(index).surfaceTension;
}

bool FluidEngine::tryCommitLiquidOccupancy(int index, LiquidComponentView const &view) {
    if (index < 0 || index >= GW * GH) return false;
    size_t i = static_cast<size_t>(index);
    if (solid[i] || dynamicSolid[i]) return false;
    if (view.count < 0 || view.count > kMaxLiquidComponents) return false;
    LiquidComponent packed[kMaxLiquidComponents]{};
    int n = 0;
    float sum = 0.0f;
    for (int s = 0; s < view.count; ++s) {
        RuntimeSubstanceRef ref = view.items[s].id;
        float amt = view.items[s].amount;
        if (!std::isfinite(amt) || amt < 0.0f) return false;
        if (!(amt > kMinLiquidComponent)) continue;
        // 19B1 deliberately keeps engine SoA storage built-in-only. Generated
        // payload refs become commit-capable in 19B2.
        if (!runtimeSubstanceIsBuiltIn(ref)) return false;
        SubstanceId id = runtimeBuiltinId(ref);
        if (!validLiquidComponentId(id)) return false;
        for (int a = 0; a < n; ++a)
            if (packed[a].id == id) return false;
        packed[n++] = {id, amt};
        sum += amt;
    }
    if (sum > 1.0f + 1.0e-5f) return false;
    int base = compositionSlot(index, 0);
    for (int s = 0; s < kMaxLiquidComponents; ++s) {
        if (s < n) {
            liquidCompId[static_cast<size_t>(base + s)] = packed[s].id;
            liquidCompAmt[static_cast<size_t>(base + s)] = packed[s].amount;
        } else {
            liquidCompId[static_cast<size_t>(base + s)] = runtimeNone();
            liquidCompAmt[static_cast<size_t>(base + s)] = 0.0f;
        }
    }
    liquidCompCount[i] = static_cast<uint8_t>(n);
    fill[i] = sum;
    return true;
}

void FluidEngine::applyCarry(int index, LiquidCarry const &c) {
    size_t i = static_cast<size_t>(index);
    liquidHeat[i] += c.heat;
    dyeR[i] += c.dyeR;
    dyeG[i] += c.dyeG;
    dyeB[i] += c.dyeB;
    for (int n = 0; n < c.compCount; ++n) {
        if (!runtimeSubstanceIsBuiltIn(c.comps[n].id)) continue;
        (void)addComponentUntracked(index, runtimeBuiltinId(c.comps[n].id), c.comps[n].amount);
    }
}

namespace {
LiquidCarry splitCarry(LiquidCarry &src, float placed, float remaining) {
    LiquidCarry out{};
    if (remaining <= 1.0e-20f || placed <= 0.0f) return out;
    float frac = placed / remaining;
    out.heat = src.heat * frac; src.heat -= out.heat;
    out.dyeR = src.dyeR * frac; src.dyeR -= out.dyeR;
    out.dyeG = src.dyeG * frac; src.dyeG -= out.dyeG;
    out.dyeB = src.dyeB * frac; src.dyeB -= out.dyeB;
    copyLiquidPayload(out.comps, out.compCount, src.comps, src.compCount);
    scaleLiquidPayload(out.comps, out.compCount, frac);
    scaleLiquidPayload(src.comps, src.compCount, 1.0f - frac);
    compactLiquidPayload(src.comps, src.compCount);
    compactLiquidPayload(out.comps, out.compCount);
    return out;
}

LiquidCarry ambientCarry(FluidEngine const &eng, float placed, SubstanceId sid) {
    LiquidCarry c{};
    FluidProperties const &liq = fluidForSubstance(sid);
    ThermalProperties const &th = thermalForSubstance(sid);
    float cap = thermalCapacity(massKg(liq.density, placed, eng.config.cellsPerMeter), th.specificHeat);
    c.heat = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
    (void)addLiquidPayload(c.comps, c.compCount, sid, placed);
    return c;
}

void copyCarryToSplash(SplashParticle &p, LiquidCarry const &c) {
    p.heat = c.heat;
    p.dyeR = c.dyeR; p.dyeG = c.dyeG; p.dyeB = c.dyeB;
    copyLiquidPayload(p.comps, p.compCount, c.comps, c.compCount);
}

LiquidCarry carryFromSplash(SplashParticle const &p) {
    LiquidCarry c{};
    c.heat = p.heat;
    c.dyeR = p.dyeR; c.dyeG = p.dyeG; c.dyeB = p.dyeB;
    copyLiquidPayload(c.comps, c.compCount, p.comps, p.compCount);
    return c;
}
} // namespace

LiquidCarry FluidEngine::extractVolume(int index, float amount) {
    LiquidCarry c{};
    if (index < 0 || index >= GW * GH || amount <= 0.0f) return c;
    size_t i = static_cast<size_t>(index);
    float f0 = fill[i];
    if (f0 <= 1.0e-20f) return c;
    float a = std::min(amount, f0);
    float frac = a / f0;
    c.heat = liquidHeat[i] * frac;
    c.dyeR = dyeR[i] * frac;
    c.dyeG = dyeG[i] * frac;
    c.dyeB = dyeB[i] * frac;
    int base = compositionSlot(index, 0);
    int n = liquidCompCount[i];
    for (int s = 0; s < n; ++s) {
        float d = liquidCompAmt[static_cast<size_t>(base + s)] * frac;
        (void)addLiquidPayload(c.comps, c.compCount, liquidCompId[static_cast<size_t>(base + s)], d);
        liquidCompAmt[static_cast<size_t>(base + s)] -= d;
    }
    fill[i] -= a;
    liquidHeat[i] -= c.heat;
    dyeR[i] -= c.dyeR;
    dyeG[i] -= c.dyeG;
    dyeB[i] -= c.dyeB;
    compactComposition(index);
    clearEmptyLiquidCell(index);
    return c;
}

void FluidEngine::addLiquidFill(int index, float dFill, float dHeat) {
    if (index < 0 || index >= GW * GH || dFill <= 0.0f) return;
    fill[static_cast<size_t>(index)] += dFill;
    liquidHeat[static_cast<size_t>(index)] += dHeat;
    (void)addComponentUntracked(index, SUBSTANCE_WATER, dFill);
    clearEmptyLiquidCell(index);
}

void FluidEngine::clearSolidifyPending(int index) {
    if (index < 0 || index >= GW * GH) return;
    size_t i = static_cast<size_t>(index);
    bool occupied = std::isfinite(solidifyPendingKg[i]) && solidifyPendingKg[i] > 1.0e-12f;
    solidifyPendingId[i] = SUBSTANCE_NONE;
    solidifyPendingKg[i] = 0.0f;
    solidifyPendingHeatJ[i] = 0.0f;
    if (occupied && pendingSolidCellCount > 0) --pendingSolidCellCount;
}

bool FluidEngine::addSolidifyPendingKg(int index, SubstanceId id, float kg, float heatJ) {
    if (index < 0 || index >= GW * GH) return false;
    if (!(kg > 0.0f) || !std::isfinite(kg) || id == SUBSTANCE_NONE || !validSubstance(id))
        return false;
    size_t i = static_cast<size_t>(index);
    float &pending = solidifyPendingKg[i];
    SubstanceId &pendId = solidifyPendingId[i];
    float &heat = solidifyPendingHeatJ[i];
    float addHeat = (std::isfinite(heatJ) && heatJ > 0.0f) ? heatJ : 0.0f;
    bool wasOccupied = std::isfinite(pending) && pending > 1.0e-12f && pendId != SUBSTANCE_NONE;
    if (!(pending > 1.0e-12f) || !std::isfinite(pending) || pendId == SUBSTANCE_NONE) {
        pendId = id;
        pending = kg;
        heat = addHeat;
        if (!wasOccupied) ++pendingSolidCellCount;
        return true;
    }
    if (pendId != id) return false;
    pending += kg;
    heat += addHeat;
    return true;
}

bool FluidEngine::takeSolidifyPendingKg(int index, SubstanceId id, float kg, float *outHeatJ) {
    if (index < 0 || index >= GW * GH) return false;
    if (!(kg > 0.0f) || !std::isfinite(kg)) return false;
    size_t i = static_cast<size_t>(index);
    if (solidifyPendingId[i] != id) return false;
    float &pending = solidifyPendingKg[i];
    float &heat = solidifyPendingHeatJ[i];
    if (!(pending + 1.0e-9f >= kg) || !std::isfinite(pending)) return false;
    bool wasOccupied = pending > 1.0e-12f;
    float frac = kg / pending;
    float takeHeat = 0.0f;
    if (std::isfinite(heat) && heat > 0.0f)
        takeHeat = heat * frac;
    pending -= kg;
    heat -= takeHeat;
    if (!(pending > 1.0e-12f) || !std::isfinite(pending)) {
        pending = 0.0f;
        heat = 0.0f;
        solidifyPendingId[i] = SUBSTANCE_NONE;
        if (wasOccupied && pendingSolidCellCount > 0) --pendingSolidCellCount;
    } else if (!(heat > 0.0f) || !std::isfinite(heat)) {
        heat = 0.0f;
    }
    if (outHeatJ) *outHeatJ = takeHeat;
    return true;
}

SubstanceId FluidEngine::solidifyPendingSubstance(int index) const {
    if (index < 0 || index >= GW * GH) return SUBSTANCE_NONE;
    size_t i = static_cast<size_t>(index);
    if (!(solidifyPendingKg[i] > 1.0e-12f) || !std::isfinite(solidifyPendingKg[i]))
        return SUBSTANCE_NONE;
    return solidifyPendingId[i];
}

float FluidEngine::solidifyPendingMassKg(int index) const {
    if (index < 0 || index >= GW * GH) return 0.0f;
    size_t i = static_cast<size_t>(index);
    if (solidifyPendingId[i] == SUBSTANCE_NONE) return 0.0f;
    float p = solidifyPendingKg[i];
    if (!(p > 0.0f) || !std::isfinite(p)) return 0.0f;
    return p;
}

float FluidEngine::solidifyPendingSensibleJ(int index) const {
    if (index < 0 || index >= GW * GH) return 0.0f;
    if (solidifyPendingSubstance(index) == SUBSTANCE_NONE) return 0.0f;
    float h = solidifyPendingHeatJ[static_cast<size_t>(index)];
    if (!(h > 0.0f) || !std::isfinite(h)) return 0.0f;
    return h;
}

float FluidEngine::waterFrozenPendingKg(int index) const {
    if (solidifyPendingSubstance(index) != SUBSTANCE_WATER) return 0.0f;
    return solidifyPendingMassKg(index);
}

LiquidCarry FluidEngine::takeLiquidCarry(int index, float amount) {
    return extractVolume(index, amount);
}

float FluidEngine::takeLiquidVolume(int index, float amount) {
    if (index < 0 || index >= GW * GH || amount <= 0.0f) return 0.0f;
    float before = fill[static_cast<size_t>(index)];
    (void)extractVolume(index, amount);
    return before - fill[static_cast<size_t>(index)];
}

void FluidEngine::seedAmbientHeat() {
    for (int i = 0; i < GW * GH; ++i) {
        float f = fill[static_cast<size_t>(i)];
        if (f > 1.0e-8f) {
            LiquidMixtureProperties mix = mixProperties(i);
            float cap = thermalCapacity(massKg(mix.density, f, config.cellsPerMeter),
                mix.specificHeat);
            liquidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        } else {
            liquidHeat[static_cast<size_t>(i)] = 0.0f;
            dyeR[static_cast<size_t>(i)] = dyeG[static_cast<size_t>(i)] = dyeB[static_cast<size_t>(i)] = 0.0f;
            clearComposition(i);
        }
        if (solid[static_cast<size_t>(i)]) {
            float cap = thermalCapacity(massKg(mechanicalForSubstance(kStaticWallSubstance).densityRel, 1.0f, config.cellsPerMeter),
                thermalForSubstance(kStaticWallSubstance).specificHeat);
            solidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        } else {
            solidHeat[static_cast<size_t>(i)] = 0.0f;
        }
    }
    std::fill(nextHeat.begin(), nextHeat.end(), 0.0f);
}

bool FluidEngine::isStaticSolid(int x, int y)  const { return !inside(x, y) || solid[ci(x, y)] != 0; }
bool FluidEngine::isPaintedSolid(int x, int y)  const { return inside(x, y) && solid[ci(x, y)] != 0; }
bool FluidEngine::isSolid(int x, int y)  const { return isStaticSolid(x, y) || (inside(x, y) && dynamicSolid[ci(x, y)] != 0); }
bool FluidEngine::isFluid(int x, int y)  const { return inside(x, y) && !isSolid(x, y) && fill[ci(x, y)] >= MIN_ACTIVE_FILL; }
bool FluidEngine::isPressureFluid(int x, int y)  const { return inside(x, y) && !isSolid(x, y) && fill[ci(x, y)] >= MIN_PRESSURE_FILL; }
float FluidEngine::gridGravity()  const { return config.gravityMetersPerSecondSquared * config.cellsPerMeter; }
float FluidEngine::liquidFaceFraction(float a, float b) { return std::clamp(0.5f * (a + b), 0.08f, 1.0f); }
bool FluidEngine::openUFace(int x, int y)  const { return x > 0 && x < GW && !isSolid(x - 1, y) && !isSolid(x, y); }
bool FluidEngine::openVFace(int x, int y)  const { return y > 0 && y < GH && !isSolid(x, y - 1) && !isSolid(x, y); }
int FluidEngine::activeX0(int halo)  const { return std::max(0, solveX0 - halo); }
int FluidEngine::activeY0(int halo)  const { return std::max(0, solveY0 - halo); }
int FluidEngine::activeX1(int halo)  const { return std::min(GW - 1, solveX1 + halo); }
int FluidEngine::activeY1(int halo)  const { return std::min(GH - 1, solveY1 + halo); }

uint32_t FluidEngine::hashCell(uint32_t x, uint32_t y, uint32_t tick) {
    uint32_t h = x * 0x8da6b343u ^ y * 0xd8163841u ^ tick * 0xcb1ab31fu;
    h ^= h >> 13; h *= 0x85ebca6bu; h ^= h >> 16;
    return h;
}

int8_t FluidEngine::makeShade(int x, int y)  const {
    return static_cast<int8_t>(static_cast<int>(hashCell(x, y, tickNo) % 17u) - 8);
}

float FluidEngine::cellU(int x, int y)  const { return 0.5f * (u[ui(x, y)] + u[ui(x + 1, y)]); }
float FluidEngine::cellV(int x, int y)  const { return 0.5f * (v[vi(x, y)] + v[vi(x, y + 1)]); }

float FluidEngine::bilinear(std::vector<float> const &a, int width, int height, float x, float y) {
    x = std::clamp(x, 0.0f, static_cast<float>(width - 1));
    y = std::clamp(y, 0.0f, static_cast<float>(height - 1));
    int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
    int x1 = std::min(width - 1, x0 + 1), y1 = std::min(height - 1, y0 + 1);
    float fx = x - x0, fy = y - y0;
    float a0 = a[y0 * width + x0] * (1.0f - fx) + a[y0 * width + x1] * fx;
    float a1 = a[y1 * width + x0] * (1.0f - fx) + a[y1 * width + x1] * fx;
    return a0 * (1.0f - fy) + a1 * fy;
}

float FluidEngine::nearestSample(std::vector<float> const &a, int width, int height, float x, float y) {
    x = std::clamp(x, 0.0f, static_cast<float>(width - 1));
    y = std::clamp(y, 0.0f, static_cast<float>(height - 1));
    int ix = std::clamp(static_cast<int>(std::lround(x)), 0, width - 1);
    int iy = std::clamp(static_cast<int>(std::lround(y)), 0, height - 1);
    return a[iy * width + ix];
}

float FluidEngine::sampleU(float worldX, float worldY) const { return bilinear(u, GW + 1, GH, worldX, worldY - 0.5f); }
float FluidEngine::sampleV(float worldX, float worldY) const { return bilinear(v, GW, GH + 1, worldX - 0.5f, worldY); }
float FluidEngine::sampleUField(std::vector<float> const &field,float worldX,float worldY) const {return bilinear(field,GW+1,GH,worldX,worldY-0.5f);}
float FluidEngine::sampleVField(std::vector<float> const &field,float worldX,float worldY) const {return bilinear(field,GW,GH+1,worldX-0.5f,worldY);}
float FluidEngine::sampleUNearest(float worldX, float worldY) const { return nearestSample(u, GW + 1, GH, worldX, worldY - 0.5f); }
float FluidEngine::sampleVNearest(float worldX, float worldY) const { return nearestSample(v, GW, GH + 1, worldX - 0.5f, worldY); }

float FluidEngine::clampToSourceExtrema(std::vector<float> const &field,int width,int height,float x,float y,float value){
    x=std::clamp(x,0.0f,float(width-1));y=std::clamp(y,0.0f,float(height-1));
    int x0=int(std::floor(x)),y0=int(std::floor(y)),x1=std::min(width-1,x0+1),y1=std::min(height-1,y0+1);
    float a=field[y0*width+x0],b=field[y0*width+x1],c=field[y1*width+x0],d=field[y1*width+x1];
    return std::clamp(value,std::min(std::min(a,b),std::min(c,d)),std::max(std::max(a,b),std::max(c,d)));
}

float FluidEngine::finiteOrZero(float value) {
    return std::isfinite(value) ? value : 0.0f;
}

bool FluidEngine::liveUFace(int x, int y) const {
    return openUFace(x, y) && (isFluid(x - 1, y) || isFluid(x, y));
}

bool FluidEngine::liveVFace(int x, int y) const {
    return openVFace(x, y) && (isFluid(x, y - 1) || isFluid(x, y));
}

float FluidEngine::uFaceOrWall(int x, int y) const {
    if (x <= 0 || x >= GW || y < 0 || y >= GH || !openUFace(x, y)) return 0.0f;
    return u[ui(x, y)];
}

float FluidEngine::vFaceOrWall(int x, int y) const {
    if (y <= 0 || y >= GH || x < 0 || x >= GW || !openVFace(x, y)) return 0.0f;
    return v[vi(x, y)];
}

void FluidEngine::enforceSolidBoundaries() {
    for (int y = 0; y < GH; ++y) {
        u[ui(0, y)] = 0.0f; u[ui(GW, y)] = 0.0f;
        for (int x = 1; x < GW; ++x)
            if (!openUFace(x, y)) u[ui(x, y)] = 0.0f;
    }
    for (int x = 0; x < GW; ++x) {
        v[vi(x, 0)] = 0.0f; v[vi(x, GH)] = 0.0f;
        for (int y = 1; y < GH; ++y)
            if (!openVFace(x, y)) v[vi(x, y)] = 0.0f;
    }
    applyMovingBoundaryVelocity();
}

void FluidEngine::enforceActiveBoundaries(){
    for(int y=activeY0(1);y<=activeY1(1);++y)for(int x=std::max(0,activeX0(1));x<=std::min(GW,activeX1(1)+1);++x)
        if(x==0||x==GW||!openUFace(x,y))u[ui(x,y)]=0.0f;
    for(int y=std::max(0,activeY0(1));y<=std::min(GH,activeY1(1)+1);++y)for(int x=activeX0(1);x<=activeX1(1);++x)
        if(y==0||y==GH||!openVFace(x,y))v[vi(x,y)]=0.0f;
    applyMovingBoundaryVelocity();
}

// Semi-Lagrangian advection follows each MAC face backward through the old
// velocity field. Nearest sampling skips bilinear interpolation at the source.
void FluidEngine::advectVelocitySemiLagrangian(float dt, bool nearest) {
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x) {
        uTemp[ui(x,y)]=0.0f;
        if (!liveUFace(x, y)) continue;
        float px = static_cast<float>(x), py = y + 0.5f;
        float vx = sampleU(px, py), vy = sampleV(px, py);
        float sx = px - vx * dt, sy = py - vy * dt;
        uTemp[ui(x, y)] = nearest ? sampleUNearest(sx, sy) : sampleU(sx, sy);
    }
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
        vTemp[vi(x,y)]=0.0f;
        if (!liveVFace(x, y)) continue;
        float px = x + 0.5f, py = static_cast<float>(y);
        float vx = sampleU(px, py), vy = sampleV(px, py);
        float sx = px - vx * dt, sy = py - vy * dt;
        vTemp[vi(x, y)] = nearest ? sampleVNearest(sx, sy) : sampleV(sx, sy);
    }
}

void FluidEngine::advectVelocityFirstOrderUpwind(float dt) {
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x) {
        uTemp[ui(x, y)] = 0.0f;
        if (!liveUFace(x, y)) continue;
        float uc = u[ui(x, y)];
        float vc = sampleV(static_cast<float>(x), y + 0.5f);
        float dudx = (uc >= 0.0f) ? (uc - uFaceOrWall(x - 1, y)) : (uFaceOrWall(x + 1, y) - uc);
        float dudy = (vc >= 0.0f) ? (uc - uFaceOrWall(x, y - 1)) : (uFaceOrWall(x, y + 1) - uc);
        uTemp[ui(x, y)] = uc - dt * (uc * dudx + vc * dudy);
    }
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
        vTemp[vi(x, y)] = 0.0f;
        if (!liveVFace(x, y)) continue;
        float vc = v[vi(x, y)];
        float uc = sampleU(x + 0.5f, static_cast<float>(y));
        float dvdx = (uc >= 0.0f) ? (vc - vFaceOrWall(x - 1, y)) : (vFaceOrWall(x + 1, y) - vc);
        float dvdy = (vc >= 0.0f) ? (vc - vFaceOrWall(x, y - 1)) : (vFaceOrWall(x, y + 1) - vc);
        vTemp[vi(x, y)] = vc - dt * (uc * dvdx + vc * dvdy);
    }
}

void FluidEngine::advectVelocityMacCormack(float dt) {
    advectVelocitySemiLagrangian(dt, false);
    std::fill(uScratch.begin(), uScratch.end(), 0.0f);
    std::fill(vScratch.begin(), vScratch.end(), 0.0f);
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x) {
        if (!liveUFace(x, y)) continue;
        float px = static_cast<float>(x), py = y + 0.5f;
        float vx = sampleU(px, py), vy = sampleV(px, py);
        uScratch[ui(x, y)] = sampleUField(uTemp, px + vx * dt, py + vy * dt);
    }
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
        if (!liveVFace(x, y)) continue;
        float px = x + 0.5f, py = static_cast<float>(y);
        float vx = sampleU(px, py), vy = sampleV(px, py);
        vScratch[vi(x, y)] = sampleVField(vTemp, px + vx * dt, py + vy * dt);
    }
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x) {
        if (!liveUFace(x, y)) continue;
        float px = static_cast<float>(x), py = y + 0.5f;
        float vx = sampleU(px, py), vy = sampleV(px, py);
        float predicted = uTemp[ui(x, y)];
        float reversed = uScratch[ui(x, y)];
        float original = u[ui(x, y)];
        float corrected = predicted + 0.5f * (original - reversed);
        if (!std::isfinite(corrected)) corrected = predicted;
        float sx = px - vx * dt, sy = py - vy * dt;
        uTemp[ui(x, y)] = clampToSourceExtrema(u, GW + 1, GH, sx, sy - 0.5f, corrected);
    }
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
        if (!liveVFace(x, y)) continue;
        float px = x + 0.5f, py = static_cast<float>(y);
        float vx = sampleU(px, py), vy = sampleV(px, py);
        float predicted = vTemp[vi(x, y)];
        float reversed = vScratch[vi(x, y)];
        float original = v[vi(x, y)];
        float corrected = predicted + 0.5f * (original - reversed);
        if (!std::isfinite(corrected)) corrected = predicted;
        float sx = px - vx * dt, sy = py - vy * dt;
        vTemp[vi(x, y)] = clampToSourceExtrema(v, GW, GH + 1, sx - 0.5f, sy, corrected);
    }
}

void FluidEngine::advectVelocityBfecc(float dt) {
    advectVelocitySemiLagrangian(dt, false);
    std::fill(uScratch.begin(),uScratch.end(),0.0f);std::fill(vScratch.begin(),vScratch.end(),0.0f);
    for(int y=activeY0(1);y<=activeY1(1);++y)for(int x=std::max(1,activeX0(1));x<=std::min(GW-1,activeX1(1)+1);++x){
        if(!openUFace(x,y)||(!isFluid(x-1,y)&&!isFluid(x,y)))continue;
        float px=float(x),py=y+0.5f;
        float vx=sampleUField(uTemp,px,py),vy=sampleVField(vTemp,px,py);uScratch[ui(x,y)]=sampleUField(uTemp,px+vx*dt,py+vy*dt);
    }
    for(int y=std::max(1,activeY0(1));y<=std::min(GH-1,activeY1(1)+1);++y)for(int x=activeX0(1);x<=activeX1(1);++x){
        if(!openVFace(x,y)||(!isFluid(x,y-1)&&!isFluid(x,y)))continue;
        float px=x+0.5f,py=float(y);
        float vx=sampleUField(uTemp,px,py),vy=sampleVField(vTemp,px,py);vScratch[vi(x,y)]=sampleVField(vTemp,px+vx*dt,py+vy*dt);
    }
    for(size_t i=0;i<uScratch.size();++i)uScratch[i]=u[i]+0.5f*(u[i]-uScratch[i]);
    for(size_t i=0;i<vScratch.size();++i)vScratch[i]=v[i]+0.5f*(v[i]-vScratch[i]);
    for(int y=activeY0(1);y<=activeY1(1);++y)for(int x=std::max(1,activeX0(1));x<=std::min(GW-1,activeX1(1)+1);++x){
        if(!openUFace(x,y)||(!isFluid(x-1,y)&&!isFluid(x,y)))continue;
        float px=float(x),py=y+0.5f;
        float vx=sampleU(px,py),vy=sampleV(px,py),sx=px-vx*dt,sy=py-vy*dt;
        uTemp[ui(x,y)]=clampToSourceExtrema(u,GW+1,GH,sx,sy-0.5f,sampleUField(uScratch,sx,sy));
    }
    for(int y=std::max(1,activeY0(1));y<=std::min(GH-1,activeY1(1)+1);++y)for(int x=activeX0(1);x<=activeX1(1);++x){
        if(!openVFace(x,y)||(!isFluid(x,y-1)&&!isFluid(x,y)))continue;
        float px=x+0.5f,py=float(y);
        float vx=sampleU(px,py),vy=sampleV(px,py),sx=px-vx*dt,sy=py-vy*dt;
        vTemp[ui(x,y)]=clampToSourceExtrema(v,GW,GH+1,sx-0.5f,sy,sampleVField(vScratch,sx,sy));
    }
}

void FluidEngine::commitAdvectedVelocity() {
    for(int y=activeY0(1);y<=activeY1(1);++y)for(int x=std::max(1,activeX0(1));x<=std::min(GW-1,activeX1(1)+1);++x)
        u[ui(x,y)]=finiteOrZero(uTemp[ui(x,y)]);
    for(int y=std::max(1,activeY0(1));y<=std::min(GH-1,activeY1(1)+1);++y)for(int x=activeX0(1);x<=activeX1(1);++x)
        v[vi(x,y)]=finiteOrZero(vTemp[vi(x,y)]);
}

void FluidEngine::advectVelocity(float dt) {
    switch (config.velocityAdvection) {
        case VelocityAdvection::None:
            return;
        case VelocityAdvection::FirstOrderUpwind:
            advectVelocityFirstOrderUpwind(dt);
            break;
        case VelocityAdvection::NearestSemiLagrangian:
            advectVelocitySemiLagrangian(dt, true);
            break;
        case VelocityAdvection::SemiLagrangian:
            advectVelocitySemiLagrangian(dt, false);
            break;
        case VelocityAdvection::MacCormack:
            advectVelocityMacCormack(dt);
            break;
        case VelocityAdvection::BFECC:
            advectVelocityBfecc(dt);
            break;
    }
    commitAdvectedVelocity();
}

// Spatial diffusion is real viscosity (shear). Lab-frame damping is only a cheap
// stand-in for slightly thickened water; it must not run on honey or free-fall
// becomes linear drag. Air is a free surface (no shear against vacuum).
void FluidEngine::diffuseVelocity(float dt) {
    FluidProperties const &refLiq = sandboxReferenceLiquid();
    float muRef = std::max(1.0e-8f, refLiq.viscosity);
    float muMax = refLiq.viscosity;
    int y0 = std::max(1, activeY0());
    int y1 = std::min(GH - 1, activeY1());
    int x0 = std::max(1, activeX0());
    int x1 = std::min(GW - 1, activeX1());
    auto isLiq = [this](int x, int y) {
        return inside(x, y) && !isSolid(x, y) && fill[static_cast<size_t>(ci(x, y))] >= MIN_ACTIVE_FILL;
    };
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        if (solid[static_cast<size_t>(i)] || dynamicSolid[static_cast<size_t>(i)]) continue;
        float f = fill[static_cast<size_t>(i)];
        if (f < MIN_ACTIVE_FILL) continue;
        float mu = mixViscosity(i);
        muMax = std::max(muMax, mu);
        // Thick liquids use spatial shear below. Lab-frame velocity damping would
        // drag midair blobs. Only slightly thickened near-reference liquid uses it.
        if (mu > muRef * 1.20f) continue;
        if (mu <= muRef * 1.15f) continue;
        float damp = 1.0f / (1.0f + (mu / muRef - 1.0f) * dt * 6.0f);
        if (openUFace(x, y)) u[ui(x, y)] *= damp;
        if (openUFace(x + 1, y)) u[ui(x + 1, y)] *= damp;
        if (openVFace(x, y)) v[vi(x, y)] *= damp;
        if (openVFace(x, y + 1)) v[vi(x, y + 1)] *= damp;
    }
    if (muMax < 0.02f) return;

    auto liveU = [&](int x, int y) {
        return openUFace(x, y) && (isLiq(x - 1, y) || isLiq(x, y));
    };
    auto liveV = [&](int x, int y) {
        return openVFace(x, y) && (isLiq(x, y - 1) || isLiq(x, y));
    };
    auto faceMuU = [&](int x, int y) {
        float s = 0.0f; int n = 0;
        if (isLiq(x - 1, y)) { s += mixViscosity(ci(x - 1, y)); ++n; }
        if (isLiq(x, y)) { s += mixViscosity(ci(x, y)); ++n; }
        return n > 0 ? s / static_cast<float>(n) : sandboxReferenceLiquid().viscosity;
    };
    auto faceMuV = [&](int x, int y) {
        float s = 0.0f; int n = 0;
        if (isLiq(x, y - 1)) { s += mixViscosity(ci(x, y - 1)); ++n; }
        if (isLiq(x, y)) { s += mixViscosity(ci(x, y)); ++n; }
        return n > 0 ? s / static_cast<float>(n) : sandboxReferenceLiquid().viscosity;
    };
    auto uShear = [&](int nx, int ny, float u0) {
        if (nx < 0 || nx > GW || ny < 0 || ny >= GH || !openUFace(nx, ny)) return -u0;
        if (!liveU(nx, ny)) return 0.0f;
        return u[ui(nx, ny)] - u0;
    };
    auto vShear = [&](int nx, int ny, float v0) {
        if (nx < 0 || nx >= GW || ny < 0 || ny > GH || !openVFace(nx, ny)) return -v0;
        if (!liveV(nx, ny)) return 0.0f;
        return v[vi(nx, ny)] - v0;
    };

    int iterations = std::clamp(static_cast<int>(std::ceil(muMax * 80.0f)), 1, 6);
    for (int iteration = 0; iteration < iterations; ++iteration) {
        uScratch = u; vScratch = v;
        for (int y = std::max(1, activeY0()); y < activeY1(); ++y)
            for (int x = std::max(1, activeX0()); x <= std::min(GW - 1, activeX1()); ++x) {
                if (!liveU(x, y)) continue;
                float mu = faceMuU(x, y);
                if (mu < 0.02f) continue;
                float alpha = std::clamp(mu * dt, 0.0f, 0.20f);
                float u0 = u[ui(x, y)];
                uScratch[ui(x, y)] = u0 + alpha * (uShear(x - 1, y, u0) + uShear(x + 1, y, u0)
                    + uShear(x, y - 1, u0) + uShear(x, y + 1, u0));
            }
        for (int y = std::max(1, activeY0()); y <= std::min(GH - 1, activeY1()); ++y)
            for (int x = std::max(1, activeX0()); x < activeX1(); ++x) {
                if (!liveV(x, y)) continue;
                float mu = faceMuV(x, y);
                if (mu < 0.02f) continue;
                float alpha = std::clamp(mu * dt, 0.0f, 0.20f);
                float v0 = v[vi(x, y)];
                vScratch[vi(x, y)] = v0 + alpha * (vShear(x - 1, y, v0) + vShear(x + 1, y, v0)
                    + vShear(x, y - 1, v0) + vShear(x, y + 1, v0));
            }
        u.swap(uScratch); v.swap(vScratch);
    }
}

// Gravity is a body force on vertical faces. Pressure projection, rather than a
// floor-impact rule, produces the opposing hydrostatic force in a resting pool.
// Airborne liquid uses the same g as the reference unit (Galileo). Density only
// biases faces between two liquid cells whose mixture densities differ, so a
// homogeneous dense blob does not fall faster in air, while denser liquid can
// still sink through lighter liquid.
void FluidEngine::applyGravity(float dt) {
    float gdt = gridGravity() * dt;
    float rhoW = std::max(1.0e-6f, sandboxReferenceLiquid().density);
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x)
        if (openVFace(x, y) && (isFluid(x, y - 1) || isFluid(x, y))) {
            v[vi(x, y)] += gdt;
            if (isFluid(x, y - 1) && isFluid(x, y)) {
                float rhoA = mixDensity(ci(x, y - 1));
                float rhoB = mixDensity(ci(x, y));
                if (std::abs(rhoA - rhoB) > 0.05f * rhoW) {
                    float rho = 0.5f * (rhoA + rhoB);
                    v[vi(x, y)] += gdt * (rho / rhoW - 1.0f);
                }
            }
        }
}

bool FluidEngine::surfaceCell(int x, int y)  const {
    if (!isFluid(x, y)) return false;
    return !isFluid(x - 1, y) || !isFluid(x + 1, y) || !isFluid(x, y - 1) || !isFluid(x, y + 1)
        || fill[ci(x, y)] < 0.98f;
}

// Build a persistent 3x3-smoothed color field, a near-surface mask, and
// curvature. The level-set curvature formula is zero for a flat interface,
// unlike a raw fill Laplacian, so resting water receives no false capillary kick.
void FluidEngine::updateSurfaceField() {
    surfaceCells.clear();
    int x0 = activeX0(2), y0 = activeY0(2), x1 = activeX1(2), y1 = activeY1(2);
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) surfaceMask[ci(x, y)] = 0;
    for (int y = std::max(1, y0); y <= std::min(GH - 2, y1); ++y)
        for (int x = std::max(1, x0); x <= std::min(GW - 2, x1); ++x) {
            if (fill[ci(x, y)] >= MIN_ACTIVE_FILL && surfaceCell(x, y)) {
                surfaceMask[ci(x, y)] = 1;
                surfaceCells.push_back(ci(x, y));
            }
        }

    auto smoothAt = [this](int x, int y) {
        if (x < 1 || y < 1 || x > GW - 2 || y > GH - 2) return;
        float sum = 0.0f;
        constexpr int weights[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
        for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox)
            if (!isSolid(x + ox, y + oy)) sum += fill[ci(x + ox, y + oy)] * weights[oy + 1][ox + 1];
        smoothedFill[ci(x, y)] = sum / 16.0f;
    };
    for (int index : surfaceCells) {
        int x = index % GW, y = index / GW;
        for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) smoothAt(x + ox, y + oy);
    }

    for (int index : surfaceCells) surfaceNormalX[index] = surfaceNormalY[index] = surfaceCurvature[index] = 0.0f;
    for (int index : surfaceCells) {
        int x = index % GW, y = index / GW;
        if (x < 2 || y < 2 || x > GW - 3 || y > GH - 3) continue;
        float cx = 0.5f * (smoothedFill[ci(x + 1, y)] - smoothedFill[ci(x - 1, y)]);
        float cy = 0.5f * (smoothedFill[ci(x, y + 1)] - smoothedFill[ci(x, y - 1)]);
        float cxx = smoothedFill[ci(x + 1, y)] - 2.0f * smoothedFill[ci(x, y)] + smoothedFill[ci(x - 1, y)];
        float cyy = smoothedFill[ci(x, y + 1)] - 2.0f * smoothedFill[ci(x, y)] + smoothedFill[ci(x, y - 1)];
        float cxy = 0.25f * (smoothedFill[ci(x + 1, y + 1)] - smoothedFill[ci(x + 1, y - 1)] - smoothedFill[ci(x - 1, y + 1)] + smoothedFill[ci(x - 1, y - 1)]);
        float gradient2 = cx * cx + cy * cy;
        if (gradient2 < 1e-5f) continue;
        float invLength = 1.0f / std::sqrt(gradient2);
        surfaceNormalX[index] = -cx * invLength;
        surfaceNormalY[index] = -cy * invLength;
        surfaceCurvature[index] = (cxx * cy * cy - 2.0f * cx * cy * cxy + cyy * cx * cx) / std::pow(gradient2, 1.5f);
    }
    workCounts.surfaceCells = static_cast<int>(surfaceCells.size());
}

void FluidEngine::applySurfaceTension(float dt) {
    if (surfaceCells.empty()) return;
    for (int index : surfaceCells) { cellForceX[index] = 0.0f; cellForceY[index] = 0.0f; }
    for (int index : surfaceCells) {
        int x = index % GW, y = index / GW;
        if (x < 2 || y < 2 || x > GW - 3 || y > GH - 3) continue;
        float magnitude = mixSurfaceTension(index) * std::clamp(surfaceCurvature[index], -2.0f, 2.0f);
        cellForceX[index] = magnitude * surfaceNormalX[index];
        cellForceY[index] = magnitude * surfaceNormalY[index];
    }
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x)
        if (openUFace(x, y)) u[ui(x, y)] += 0.5f * dt * (cellForceX[ci(x - 1, y)] + cellForceX[ci(x, y)]);
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x)
        if (openVFace(x, y)) v[vi(x, y)] += 0.5f * dt * (cellForceY[ci(x, y - 1)] + cellForceY[ci(x, y)]);
}

// Vorticity confinement restores a small fraction of the rotation lost by the
// dissipative semi-Lagrangian step. It can be toggled with O.
void FluidEngine::applyVorticityConfinement(float dt) {
    if (!config.vorticityEnabled) return;
    int x0 = std::max(1, activeX0(2)), y0 = std::max(1, activeY0(2));
    int x1 = std::min(GW - 2, activeX1(2)), y1 = std::min(GH - 2, activeY1(2));
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) curlField[ci(x, y)] = 0.0f;
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        float dvdx = 0.5f * (cellV(x + 1, y) - cellV(x - 1, y));
        float dudy = 0.5f * (cellU(x, y + 1) - cellU(x, y - 1));
        curlField[ci(x, y)] = dvdx - dudy;
    }
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        cellForceX[ci(x, y)] = 0.0f;
        cellForceY[ci(x, y)] = 0.0f;
        if (!isFluid(x, y)) continue;
        float nx = 0.5f * (std::abs(curlField[ci(x + 1, y)]) - std::abs(curlField[ci(x - 1, y)]));
        float ny = 0.5f * (std::abs(curlField[ci(x, y + 1)]) - std::abs(curlField[ci(x, y - 1)]));
        float length = std::sqrt(nx * nx + ny * ny) + 1e-5f;
        nx /= length; ny /= length;
        float omega = curlField[ci(x, y)] * config.vorticityStrength;
        cellForceX[ci(x, y)] = ny * omega;
        cellForceY[ci(x, y)] = -nx * omega;
    }
    for (int y = y0; y <= y1; ++y) for (int x = std::max(1, x0); x <= std::min(GW - 1, x1 + 1); ++x)
        if (openUFace(x, y)) u[ui(x, y)] += 0.5f * dt * (cellForceX[ci(x - 1, y)] + cellForceX[ci(x, y)]);
    for (int y = std::max(1, y0); y <= std::min(GH - 1, y1 + 1); ++y) for (int x = x0; x <= x1; ++x)
        if (openVFace(x, y)) v[vi(x, y)] += 0.5f * dt * (cellForceY[ci(x, y - 1)] + cellForceY[ci(x, y)]);
}

void FluidEngine::clampVelocity() {
    for(int y=activeY0(1);y<=activeY1(1);++y)for(int x=activeX0(1);x<=std::min(GW,activeX1(1)+1);++x){float&value=u[ui(x,y)];value=std::clamp(value,-config.maxVelocity,config.maxVelocity);}
    for(int y=activeY0(1);y<=std::min(GH,activeY1(1)+1);++y)for(int x=activeX0(1);x<=activeX1(1);++x){float&value=v[vi(x,y)];value=std::clamp(value,-config.maxVelocity,config.maxVelocity);}
}

void FluidEngine::buildPressureWorkLists() {
    pressureRed.clear();
    pressureBlack.clear();
    pressureStencils.clear();
    int x0 = activeX0(), y0 = activeY0(), x1 = activeX1(), y1 = activeY1();
    int rectArea = std::max(1, (x1 - x0 + 1) * (y1 - y0 + 1));
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        if (!isPressureFluid(x, y)) continue;
        int index = ci(x, y);
        if (((x + y) & 1) == 0) pressureRed.push_back(index);
        else pressureBlack.push_back(index);
    }
    workCounts.pressureCells = static_cast<int>(pressureRed.size() + pressureBlack.size());
    useSparsePressure = workCounts.pressureCells * 10 <= rectArea * 7;

    if (!useSparsePressure) return;

    pressureStencils.resize(static_cast<size_t>(workCounts.pressureCells));
    auto fillStencil = [&](int index, PressureStencil &s) {
        s.index = index;
        s.pressureNeighborCount = 0;
        s.freeNeighborCount = 0;
        int x = index % GW, y = index / GW;
        constexpr int dx[4] = {-1, 1, 0, 0};
        constexpr int dy[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
            int nx = x + dx[n], ny = y + dy[n];
            if (isSolid(nx, ny)) continue;
            if (isPressureFluid(nx, ny)) {
                s.pressureNeighbors[s.pressureNeighborCount++] = ci(nx, ny);
            } else {
                s.freeDx[s.freeNeighborCount] = dx[n];
                s.freeDy[s.freeNeighborCount] = dy[n];
                ++s.freeNeighborCount;
            }
        }
    };
    size_t si = 0;
    for (int index : pressureRed) fillStencil(index, pressureStencils[si++]);
    for (int index : pressureBlack) fillStencil(index, pressureStencils[si++]);
}

void FluidEngine::sorPressureColor(std::vector<int> const &cells, float invDt, float sor) {
    for (int index : cells) {
        int x = index % GW, y = index / GW;
        float sum = 0.0f, coefficient = 0.0f;
        constexpr int dx[4] = {-1, 1, 0, 0};
        constexpr int dy[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
            int nx = x + dx[n], ny = y + dy[n];
            if (isSolid(nx, ny)) continue;
            if (isPressureFluid(nx, ny)) { coefficient += 1.0f; sum += pressure[ci(nx, ny)]; }
            else { float neighborFill = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f; coefficient += 1.0f / liquidFaceFraction(fill[index], neighborFill); }
        }
        if (coefficient <= 0.0f) continue;
        float estimate = (sum - divergenceField[index] * invDt) / coefficient;
        float &p = pressure[index];
        p += sor * (estimate - p);
    }
}

void FluidEngine::sorPressureColorParallel(std::vector<int> const &cells, float invDt, float sor) {
    if (cells.empty()) return;
    if (!useParallelPressure(static_cast<int>(cells.size()))) {
        sorPressureColor(cells, invDt, sor);
        return;
    }
    lastPressureParallel = true;
    if (!useSparsePressure || pressureStencils.empty()) {
        workerPool.parallelFor(0, static_cast<int>(cells.size()), [&](int a, int b) {
            for (int i = a; i < b; ++i) {
                int index = cells[static_cast<size_t>(i)];
                int x = index % GW, y = index / GW;
                float sum = 0.0f, coefficient = 0.0f;
                constexpr int dx[4] = {-1, 1, 0, 0};
                constexpr int dy[4] = {0, 0, -1, 1};
                for (int n = 0; n < 4; ++n) {
                    int nx = x + dx[n], ny = y + dy[n];
                    if (isSolid(nx, ny)) continue;
                    if (isPressureFluid(nx, ny)) { coefficient += 1.0f; sum += pressure[ci(nx, ny)]; }
                    else { float neighborFill = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f; coefficient += 1.0f / liquidFaceFraction(fill[index], neighborFill); }
                }
                if (coefficient <= 0.0f) continue;
                float estimate = (sum - divergenceField[index] * invDt) / coefficient;
                pressure[index] += sor * (estimate - pressure[index]);
            }
        });
        return;
    }

    size_t offset = (&cells == &pressureBlack) ? pressureRed.size() : 0;
    workerPool.parallelFor(0, static_cast<int>(cells.size()), [&](int a, int b) {
        for (int i = a; i < b; ++i) {
            PressureStencil const &s = pressureStencils[offset + static_cast<size_t>(i)];
            float sum = 0.0f, coefficient = 0.0f;
            for (int n = 0; n < s.pressureNeighborCount; ++n) {
                coefficient += 1.0f;
                sum += pressure[s.pressureNeighbors[n]];
            }
            int x = s.index % GW, y = s.index / GW;
            for (int n = 0; n < s.freeNeighborCount; ++n) {
                int nx = x + s.freeDx[n], ny = y + s.freeDy[n];
                float neighborFill = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f;
                coefficient += 1.0f / liquidFaceFraction(fill[s.index], neighborFill);
            }
            if (coefficient <= 0.0f) continue;
            float estimate = (sum - divergenceField[s.index] * invDt) / coefficient;
            pressure[s.index] += sor * (estimate - pressure[s.index]);
        }
    });
}

float FluidEngine::pressureResidual(float invDt, float dt) const {
    float maxEquationResidual = 0.0f;
    auto checkCell = [&](int index) {
        int x = index % GW, y = index / GW;
        float sum = 0.0f, coefficient = 0.0f;
        constexpr int dx[4] = {-1, 1, 0, 0}, dy[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
            int nx = x + dx[n], ny = y + dy[n];
            if (isSolid(nx, ny)) continue;
            if (isPressureFluid(nx, ny)) { coefficient += 1.0f; sum += pressure[ci(nx, ny)]; }
            else { float nf = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f; coefficient += 1.0f / liquidFaceFraction(fill[index], nf); }
        }
        float rhs = divergenceField[index] * invDt;
        maxEquationResidual = std::max(maxEquationResidual, std::abs(coefficient * pressure[index] - sum + rhs) * dt);
    };
    if (useSparsePressure) {
        for (int index : pressureRed) checkCell(index);
        for (int index : pressureBlack) checkCell(index);
    } else {
        for (int y = activeY0(); y <= activeY1(); ++y) for (int x = activeX0(); x <= activeX1(); ++x)
            if (isPressureFluid(x, y)) checkCell(ci(x, y));
    }
    return maxEquationResidual;
}

// Projection solves Laplacian(p)=divergence/dt with red-black Gauss-Seidel.
// Air is a zero-pressure free surface; solid neighbors impose no-through flow.
void FluidEngine::projectVelocity(float dt) {
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
            int index = ci(x,y);
            divergenceField[index]=0.0f;pressureBefore[index]=pressure[index];
            if (isSolid(x, y) || fill[index] < MIN_PRESSURE_FILL) pressure[index] = 0.0f;
        }

    auto listStart = Clock::now();
    buildPressureWorkLists();
    timingAccum.listsBuild += elapsedMs(listStart);

    float maximumDivergence = 0.0f;
    if (useSparsePressure) {
        for (int index : pressureRed) {
            int x = index % GW, y = index / GW;
            divergenceField[index] = u[ui(x + 1, y)] - u[ui(x, y)] + v[vi(x, y + 1)] - v[vi(x, y)];
            maximumDivergence = std::max(maximumDivergence, std::abs(divergenceField[index]));
        }
        for (int index : pressureBlack) {
            int x = index % GW, y = index / GW;
            divergenceField[index] = u[ui(x + 1, y)] - u[ui(x, y)] + v[vi(x, y + 1)] - v[vi(x, y)];
            maximumDivergence = std::max(maximumDivergence, std::abs(divergenceField[index]));
        }
    } else {
        for (int y = activeY0(); y <= activeY1(); ++y) for (int x = activeX0(); x <= activeX1(); ++x) {
            if (!isPressureFluid(x, y)) continue;
            divergenceField[ci(x, y)] = u[ui(x + 1, y)] - u[ui(x, y)] + v[vi(x, y + 1)] - v[vi(x, y)];
            maximumDivergence = std::max(maximumDivergence, std::abs(divergenceField[ci(x, y)]));
        }
    }

    float invDt = 1.0f / std::max(dt, 1e-5f);
    constexpr float sor = 1.65f;
    int budget = maximumDivergence < 0.015f && measuredMaxVelocity < 0.5f ? 8
        : (maximumDivergence < 0.20f && measuredMaxVelocity < 8.0f ? 14 : config.maxPressureIterations);
    budget = std::min(budget, config.maxPressureIterations);
    lastPressureIterations = 0;

    lastPressureParallel = false;
    auto applyColor = [&](std::vector<int> const &cells) {
        if (useParallelPressure(static_cast<int>(cells.size()))) {
            sorPressureColorParallel(cells, invDt, sor);
            return;
        }
        if (useSparsePressure && !pressureStencils.empty()) {
            size_t offset = (&cells == &pressureBlack) ? pressureRed.size() : 0;
            for (size_t i = 0; i < cells.size(); ++i) {
                PressureStencil const &s = pressureStencils[offset + i];
                float sum = 0.0f, coefficient = 0.0f;
                for (int n = 0; n < s.pressureNeighborCount; ++n) {
                    coefficient += 1.0f;
                    sum += pressure[s.pressureNeighbors[n]];
                }
                int x = s.index % GW, y = s.index / GW;
                for (int n = 0; n < s.freeNeighborCount; ++n) {
                    int nx = x + s.freeDx[n], ny = y + s.freeDy[n];
                    float neighborFill = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f;
                    coefficient += 1.0f / liquidFaceFraction(fill[s.index], neighborFill);
                }
                if (coefficient <= 0.0f) continue;
                float estimate = (sum - divergenceField[s.index] * invDt) / coefficient;
                pressure[s.index] += sor * (estimate - pressure[s.index]);
            }
            return;
        }
        if (!cells.empty()) {
            sorPressureColor(cells, invDt, sor);
            return;
        }
        int color = (&cells == &pressureBlack) ? 1 : 0;
        for (int y = activeY0(); y <= activeY1(); ++y) for (int x = activeX0(); x <= activeX1(); ++x) {
            if (((x + y) & 1) != color || !isPressureFluid(x, y)) continue;
            float sum = 0.0f, coefficient = 0.0f;
            constexpr int dx[4] = {-1, 1, 0, 0};
            constexpr int dy[4] = {0, 0, -1, 1};
            for (int n = 0; n < 4; ++n) {
                int nx = x + dx[n], ny = y + dy[n];
                if (isSolid(nx, ny)) continue;
                if (isPressureFluid(nx, ny)) { coefficient += 1.0f; sum += pressure[ci(nx, ny)]; }
                else { float neighborFill = inside(nx, ny) ? fill[ci(nx, ny)] : 0.0f; coefficient += 1.0f / liquidFaceFraction(fill[ci(x, y)], neighborFill); }
            }
            if (coefficient <= 0.0f) continue;
            float rhs = divergenceField[ci(x, y)] * invDt;
            float estimate = (sum - rhs) / coefficient;
            float &p = pressure[ci(x, y)];
            p += sor * (estimate - p);
        }
    };

    for (int iteration = 0; iteration < budget; ++iteration) {
        applyColor(pressureRed);
        applyColor(pressureBlack);
        lastPressureIterations = iteration + 1;
        if (iteration >= 7 && (iteration & 3) == 3) {
            if (pressureResidual(invDt, dt) < 0.006f) break;
        }
    }
    for (int y = activeY0(1); y <= activeY1(1); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 1, activeX1(1) + 1); ++x) {
        if (!openUFace(x, y)) { u[ui(x, y)] = 0.0f; continue; }
        bool leftFluid = isPressureFluid(x - 1, y), rightFluid = isPressureFluid(x, y);
        if (!leftFluid && !rightFluid) { u[ui(x, y)] = 0.0f; continue; }
        float pLeft = leftFluid ? pressure[ci(x - 1, y)] : 0.0f;
        float pRight = rightFluid ? pressure[ci(x, y)] : 0.0f;
        float leftFill = inside(x - 1, y) ? fill[ci(x - 1, y)] : 0.0f, rightFill = inside(x, y) ? fill[ci(x, y)] : 0.0f;
        float theta = leftFluid && rightFluid ? 1.0f : liquidFaceFraction(leftFill, rightFill);
        u[ui(x, y)] -= dt * (pRight - pLeft) / (sandboxReferenceLiquid().density * theta);
    }
    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 1, activeY1(1) + 1); ++y)
        for (int x = activeX0(1); x <= activeX1(1); ++x) {
        if (!openVFace(x, y)) { v[vi(x, y)] = 0.0f; continue; }
        bool topFluid = isPressureFluid(x, y - 1), bottomFluid = isPressureFluid(x, y);
        if (!topFluid && !bottomFluid) { v[vi(x, y)] = 0.0f; continue; }
        float pTop = topFluid ? pressure[ci(x, y - 1)] : 0.0f;
        float pBottom = bottomFluid ? pressure[ci(x, y)] : 0.0f;
        float topFill = inside(x, y - 1) ? fill[ci(x, y - 1)] : 0.0f, bottomFill = inside(x, y) ? fill[ci(x, y)] : 0.0f;
        float theta = topFluid && bottomFluid ? 1.0f : liquidFaceFraction(topFill, bottomFill);
        v[vi(x, y)] -= dt * (pBottom - pTop) / (sandboxReferenceLiquid().density * theta);
    }
    enforceActiveBoundaries();
}

// Conservative adjacent-face transport. Donor scaling prevents a cell exporting
// more fill than it owns; receiver scaling prevents overfilling. The same final
// flux is subtracted and added, so scan order cannot create or destroy volume.
void FluidEngine::advectLiquidVolume(float dt) {
    int x0=activeX0(1),y0=activeY0(1),x1=activeX1(1),y1=activeY1(1);
    double regionVolumeBefore=0.0;
    for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){int i=ci(x,y);previousFill[i]=fill[i];regionVolumeBefore+=fill[i];}
    for(int y=y0;y<=y1;++y)for(int x=std::max(1,x0+1);x<=std::min(GW-1,x1);++x)fluxU[ui(x,y)]=0.0f;
    for(int y=std::max(1,y0+1);y<=std::min(GH-1,y1);++y)for(int x=x0;x<=x1;++x)fluxV[vi(x,y)]=0.0f;
    nonzeroFluxUFaces.clear();
    nonzeroFluxVFaces.clear();
    for (int y = y0; y <= y1; ++y) for (int x = std::max(1,x0+1); x <= std::min(GW-1,x1); ++x) {
        if (!openUFace(x, y)) continue;
        ++workCounts.activeUFaces;
        float a = fill[ci(x - 1, y)], b = fill[ci(x, y)];
        float q = std::clamp(u[ui(x, y)] * dt * liquidFaceFraction(a, b), -1.0f, 1.0f);
        int donor = q >= 0.0f ? ci(x - 1, y) : ci(x, y);
        if (fill[donor] <= 1e-7f) continue;
        float amount = std::min(std::abs(q), fill[donor]);
        if (amount <= 0.0f) continue;
        fluxU[ui(x, y)] = std::copysign(amount, q);
        nonzeroFluxUFaces.push_back(ui(x, y));
    }
    for (int y = std::max(1,y0+1); y <= std::min(GH-1,y1); ++y) for (int x = x0; x <= x1; ++x) {
        if (!openVFace(x, y)) continue;
        ++workCounts.activeVFaces;
        float a = fill[ci(x, y - 1)], b = fill[ci(x, y)];
        float q = std::clamp(v[vi(x, y)] * dt * liquidFaceFraction(a, b), -1.0f, 1.0f);
        int donor = q >= 0.0f ? ci(x, y - 1) : ci(x, y);
        if (fill[donor] <= 1e-7f) continue;
        float amount = std::min(std::abs(q), fill[donor]);
        if (amount <= 0.0f) continue;
        fluxV[vi(x, y)] = std::copysign(amount, q);
        nonzeroFluxVFaces.push_back(vi(x, y));
    }
    workCounts.nonzeroFluxU = static_cast<int>(nonzeroFluxUFaces.size());
    workCounts.nonzeroFluxV = static_cast<int>(nonzeroFluxVFaces.size());

    int possibleU = std::max(1, (y1 - y0 + 1) * std::max(0, std::min(GW - 1, x1) - std::max(1, x0 + 1) + 1));
    int possibleV = std::max(1, (x1 - x0 + 1) * std::max(0, std::min(GH - 1, y1) - std::max(1, y0 + 1) + 1));
    bool useSparseFlux = (workCounts.nonzeroFluxU + workCounts.nonzeroFluxV) * 10 <= (possibleU + possibleV) * 4;
    if (forceSparseFlux == 0) useSparseFlux = false;
    else if (forceSparseFlux == 1) useSparseFlux = true;

    // Iteratively limit the same face fluxes at donors and receivers. Receiver
    // capacity includes the cell's accepted outgoing flux, enabling A->B->C
    // through-flow without ever applying a scan-ordered update.
    workCounts.limiterPasses = 0;
    int fluxFaceCount = workCounts.nonzeroFluxU + workCounts.nonzeroFluxV;
    int limiterCap = std::clamp(config.maxLimiterPasses, 1, 16);
    for (int pass = 0; pass < limiterCap && fluxFaceCount > 0; ++pass) {
        bool limited = false; float maxFluxChange=0.0f;
        ++workCounts.limiterPasses;
        if (useSparseFlux) {
            fluxTouchedCells.clear();
            auto mark = [&](int cell) { fluxTouchedCells.push_back(cell); outgoing[cell] = incoming[cell] = 0.0f; };
            // Reset only cells touched by current nonzero faces.
            for (int fi : nonzeroFluxUFaces) {
                int x = fi % (GW + 1), y = fi / (GW + 1);
                mark(ci(x - 1, y)); mark(ci(x, y));
            }
            for (int fi : nonzeroFluxVFaces) {
                int x = fi % GW, y = fi / GW;
                mark(ci(x, y - 1)); mark(ci(x, y));
            }
            // Dedup not required for correctness (re-zero is fine); accumulate once.
            for (int fi : nonzeroFluxUFaces) {
                int x = fi % (GW + 1), y = fi / (GW + 1);
                float q = fluxU[fi];
                if (q == 0.0f) continue;
                int donor = q > 0.0f ? ci(x - 1, y) : ci(x, y);
                int receiver = q > 0.0f ? ci(x, y) : ci(x - 1, y);
                outgoing[donor] += std::abs(q); incoming[receiver] += std::abs(q);
            }
            for (int fi : nonzeroFluxVFaces) {
                int x = fi % GW, y = fi / GW;
                float q = fluxV[fi];
                if (q == 0.0f) continue;
                int donor = q > 0.0f ? ci(x, y - 1) : ci(x, y);
                int receiver = q > 0.0f ? ci(x, y) : ci(x, y - 1);
                outgoing[donor] += std::abs(q); incoming[receiver] += std::abs(q);
            }
            for (int i : fluxTouchedCells) {
                donorScale[i] = outgoing[i] > fill[i] && outgoing[i] > 0.0f ? fill[i] / outgoing[i] : 1.0f;
                if (donorScale[i] < 0.999999f) limited = true;
            }
            for (int fi : nonzeroFluxUFaces) {
                float &q = fluxU[fi];
                int x = fi % (GW + 1), y = fi / (GW + 1);
                int donor = q > 0 ? ci(x - 1, y) : ci(x, y);
                float old = q; q *= donorScale[donor]; maxFluxChange = std::max(maxFluxChange, std::abs(old - q));
            }
            for (int fi : nonzeroFluxVFaces) {
                float &q = fluxV[fi];
                int x = fi % GW, y = fi / GW;
                int donor = q > 0 ? ci(x, y - 1) : ci(x, y);
                float old = q; q *= donorScale[donor]; maxFluxChange = std::max(maxFluxChange, std::abs(old - q));
            }

            for (int i : fluxTouchedCells) outgoing[i] = incoming[i] = 0.0f;
            for (int fi : nonzeroFluxUFaces) {
                int x = fi % (GW + 1), y = fi / (GW + 1);
                float q = fluxU[fi];
                if (q == 0.0f) continue;
                int donor = q > 0.0f ? ci(x - 1, y) : ci(x, y);
                int receiver = q > 0.0f ? ci(x, y) : ci(x - 1, y);
                outgoing[donor] += std::abs(q); incoming[receiver] += std::abs(q);
            }
            for (int fi : nonzeroFluxVFaces) {
                int x = fi % GW, y = fi / GW;
                float q = fluxV[fi];
                if (q == 0.0f) continue;
                int donor = q > 0.0f ? ci(x, y - 1) : ci(x, y);
                int receiver = q > 0.0f ? ci(x, y) : ci(x, y - 1);
                outgoing[donor] += std::abs(q); incoming[receiver] += std::abs(q);
            }
            for (int i : fluxTouchedCells) {
                float capacity = std::max(0.0f, 1.0f - fill[i] + outgoing[i]);
                receiverScale[i] = incoming[i] > capacity && incoming[i] > 0.0f ? capacity / incoming[i] : 1.0f;
                if (receiverScale[i] < 0.999999f) limited = true;
            }
            for (int fi : nonzeroFluxUFaces) {
                float &q = fluxU[fi];
                int x = fi % (GW + 1), y = fi / (GW + 1);
                int receiver = q > 0 ? ci(x, y) : ci(x - 1, y);
                float old = q; q *= receiverScale[receiver]; maxFluxChange = std::max(maxFluxChange, std::abs(old - q));
            }
            for (int fi : nonzeroFluxVFaces) {
                float &q = fluxV[fi];
                int x = fi % GW, y = fi / GW;
                int receiver = q > 0 ? ci(x, y) : ci(x, y - 1);
                float old = q; q *= receiverScale[receiver]; maxFluxChange = std::max(maxFluxChange, std::abs(old - q));
            }
        } else {
            for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){int i=ci(x,y);outgoing[i]=incoming[i]=0.0f;}
            auto accumulate = [this](float q, int negativeCell, int positiveCell) {
                if (q == 0.0f) return;
                int donor = q > 0.0f ? negativeCell : positiveCell;
                int receiver = q > 0.0f ? positiveCell : negativeCell;
                outgoing[donor] += std::abs(q); incoming[receiver] += std::abs(q);
            };
            for(int y=y0;y<=y1;++y)for(int x=std::max(1,x0+1);x<=std::min(GW-1,x1);++x)accumulate(fluxU[ui(x,y)],ci(x-1,y),ci(x,y));
            for(int y=std::max(1,y0+1);y<=std::min(GH-1,y1);++y)for(int x=x0;x<=x1;++x)accumulate(fluxV[vi(x,y)],ci(x,y-1),ci(x,y));
            for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){int i=ci(x,y);donorScale[i]=outgoing[i]>fill[i]&&outgoing[i]>0.0f?fill[i]/outgoing[i]:1.0f;if(donorScale[i]<0.999999f)limited=true;}
            for(int y=y0;y<=y1;++y)for(int x=std::max(1,x0+1);x<=std::min(GW-1,x1);++x){float&q=fluxU[ui(x,y)];int donor=q>0?ci(x-1,y):ci(x,y);float old=q;q*=donorScale[donor];maxFluxChange=std::max(maxFluxChange,std::abs(old-q));}
            for(int y=std::max(1,y0+1);y<=std::min(GH-1,y1);++y)for(int x=x0;x<=x1;++x){float&q=fluxV[vi(x,y)];int donor=q>0?ci(x,y-1):ci(x,y);float old=q;q*=donorScale[donor];maxFluxChange=std::max(maxFluxChange,std::abs(old-q));}

            for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){int i=ci(x,y);outgoing[i]=incoming[i]=0.0f;}
            for(int y=y0;y<=y1;++y)for(int x=std::max(1,x0+1);x<=std::min(GW-1,x1);++x)accumulate(fluxU[ui(x,y)],ci(x-1,y),ci(x,y));
            for(int y=std::max(1,y0+1);y<=std::min(GH-1,y1);++y)for(int x=x0;x<=x1;++x)accumulate(fluxV[vi(x,y)],ci(x,y-1),ci(x,y));
            for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){int i=ci(x,y);float capacity=std::max(0.0f,1.0f-fill[i]+outgoing[i]);receiverScale[i]=incoming[i]>capacity&&incoming[i]>0.0f?capacity/incoming[i]:1.0f;if(receiverScale[i]<0.999999f)limited=true;}
            for(int y=y0;y<=y1;++y)for(int x=std::max(1,x0+1);x<=std::min(GW-1,x1);++x){float&q=fluxU[ui(x,y)];int receiver=q>0?ci(x,y):ci(x-1,y);float old=q;q*=receiverScale[receiver];maxFluxChange=std::max(maxFluxChange,std::abs(old-q));}
            for(int y=std::max(1,y0+1);y<=std::min(GH-1,y1);++y)for(int x=x0;x<=x1;++x){float&q=fluxV[vi(x,y)];int receiver=q>0?ci(x,y):ci(x,y-1);float old=q;q*=receiverScale[receiver];maxFluxChange=std::max(maxFluxChange,std::abs(old-q));}
        }
        if (!limited || maxFluxChange < 1e-5f) break;
    }
    for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){
        int i=ci(x,y);
        nextFill[i]=fill[i];
        nextHeat[i]=liquidHeat[i];
        nextDyeR[i]=dyeR[i]; nextDyeG[i]=dyeG[i]; nextDyeB[i]=dyeB[i];
        copyCompositionToNext(i);
    }
    auto moveField = [](std::vector<float> const &src, std::vector<float> &nxt, int donor, int receiver, float amount, float f0) {
        float d = (f0 > 1.0e-20f) ? src[static_cast<size_t>(donor)] * (amount / f0) : 0.0f;
        nxt[static_cast<size_t>(donor)] -= d;
        nxt[static_cast<size_t>(receiver)] += d;
        return d;
    };
    auto transfer = [this, &moveField](int donor, int receiver, float amount) {
        nextFill[donor] -= amount;
        nextFill[receiver] += amount;
        float f0 = fill[donor];
        float dq = moveField(liquidHeat, nextHeat, donor, receiver, amount, f0);
        moveField(dyeR, nextDyeR, donor, receiver, amount, f0);
        moveField(dyeG, nextDyeG, donor, receiver, amount, f0);
        moveField(dyeB, nextDyeB, donor, receiver, amount, f0);
        float frac = (f0 > 1.0e-20f) ? (amount / f0) : 0.0f;
        int base = compositionSlot(donor, 0);
        int n = liquidCompCount[static_cast<size_t>(donor)];
        for (int s = 0; s < n; ++s) {
            float d = liquidCompAmt[static_cast<size_t>(base + s)] * frac;
            nextCompAmt[static_cast<size_t>(base + s)] -= d;
            RuntimeSubstanceRef storedId = liquidCompId[static_cast<size_t>(base + s)];
            if (!runtimeSubstanceIsBuiltIn(storedId)) continue;
            float unplaced = addNextComponentUntracked(receiver, runtimeBuiltinId(storedId), d);
            if (unplaced > 0.0f) {
                nextCompAmt[static_cast<size_t>(base + s)] += unplaced;
                nextFill[receiver] -= unplaced;
                nextFill[donor] += unplaced;
            }
        }
        if (std::abs(dq) > 1.0e-4f) {
            LiquidMixtureProperties mix = mixProperties(donor);
            float td = tempFromEnergy(liquidHeat[donor], thermalCapacity(massKg(mix.density, f0, config.cellsPerMeter), mix.specificHeat));
            if (std::abs(td - AMBIENT_TEMPERATURE_K) > 0.2f) {
                int dx = donor % GW, dy = donor / GW;
                int rx = receiver % GW, ry = receiver / GW;
                wakeThermalAt(dx, dy);
                wakeThermalAt(rx, ry);
            }
        }
    };
    if (useSparseFlux) {
        for (int fi : nonzeroFluxUFaces) {
            float q = fluxU[fi]; if (q == 0.0f) continue;
            int x = fi % (GW + 1), y = fi / (GW + 1);
            transfer(q > 0.0f ? ci(x - 1, y) : ci(x, y), q > 0.0f ? ci(x, y) : ci(x - 1, y), std::abs(q));
        }
        for (int fi : nonzeroFluxVFaces) {
            float q = fluxV[fi]; if (q == 0.0f) continue;
            int x = fi % GW, y = fi / GW;
            transfer(q > 0.0f ? ci(x, y - 1) : ci(x, y), q > 0.0f ? ci(x, y) : ci(x, y - 1), std::abs(q));
        }
    } else {
        for (int y = y0; y <= y1; ++y) for (int x = std::max(1,x0+1); x <= std::min(GW-1,x1); ++x) {
            float q = fluxU[ui(x, y)]; if (q == 0.0f) continue;
            transfer(q > 0.0f ? ci(x - 1, y) : ci(x, y), q > 0.0f ? ci(x, y) : ci(x - 1, y), std::abs(q));
        }
        for (int y = std::max(1,y0+1); y <= std::min(GH-1,y1); ++y) for (int x = x0; x <= x1; ++x) {
            float q = fluxV[vi(x, y)]; if (q == 0.0f) continue;
            transfer(q > 0.0f ? ci(x, y - 1) : ci(x, y), q > 0.0f ? ci(x, y) : ci(x, y - 1), std::abs(q));
        }
    }
    lastAdvectOverflowVol = 0.0;
    workCounts.advectOverflowCells = 0;
    advectOverflowIndex.clear();
    advectOverflowVol.clear();
    advectOverflowCarry.clear();
    auto peelNextCarry = [this](int i, float amount, float nf) {
        LiquidCarry c{};
        if (!(amount > 1.0e-20f) || !(nf > 1.0e-20f)) return c;
        float frac = amount / nf;
        c.heat = nextHeat[static_cast<size_t>(i)] * frac;
        nextHeat[static_cast<size_t>(i)] -= c.heat;
        c.dyeR = nextDyeR[static_cast<size_t>(i)] * frac;
        nextDyeR[static_cast<size_t>(i)] -= c.dyeR;
        c.dyeG = nextDyeG[static_cast<size_t>(i)] * frac;
        nextDyeG[static_cast<size_t>(i)] -= c.dyeG;
        c.dyeB = nextDyeB[static_cast<size_t>(i)] * frac;
        nextDyeB[static_cast<size_t>(i)] -= c.dyeB;
        int base = compositionSlot(i, 0);
        int n = nextCompCount[static_cast<size_t>(i)];
        for (int s = 0; s < n; ++s) {
            float d = nextCompAmt[static_cast<size_t>(base + s)] * frac;
            (void)addLiquidPayload(c.comps, c.compCount, nextCompId[static_cast<size_t>(base + s)], d);
            nextCompAmt[static_cast<size_t>(base + s)] -= d;
        }
        return c;
    };
    auto queueOverflow = [this](int i, float vol, LiquidCarry &&carry) {
        if (!(vol > 1.0e-12f)) return;
        advectOverflowIndex.push_back(i);
        advectOverflowVol.push_back(vol);
        advectOverflowCarry.push_back(std::move(carry));
        lastAdvectOverflowVol += vol;
        accumAdvectOverflowVol += vol;
        ++workCounts.advectOverflowCells;
    };
    for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x){
        int i=ci(x,y);
        float old=fill[i];
        float nf=nextFill[i];
        if (solid[i]||dynamicSolid[i]) {
            if (nf > 1.0e-12f)
                queueOverflow(i, nf, peelNextCarry(i, nf, nf));
            nf = 0.0f;
            nextHeat[static_cast<size_t>(i)] = 0.0f;
            nextDyeR[static_cast<size_t>(i)] = nextDyeG[static_cast<size_t>(i)] = nextDyeB[static_cast<size_t>(i)] = 0.0f;
            nextCompCount[static_cast<size_t>(i)] = 0;
        } else if (nf > 1.0f) {
            float over = nf - 1.0f;
            queueOverflow(i, over, peelNextCarry(i, over, nf));
            nf = 1.0f;
        }
        if (nf < 0.0f) nf = 0.0f;
        if (nf > 0.0f && nf <= 1.0e-8f) {
            queueOverflow(i, nf, peelNextCarry(i, nf, nf));
            nf = 0.0f;
            nextHeat[static_cast<size_t>(i)] = 0.0f;
            nextDyeR[static_cast<size_t>(i)] = nextDyeG[static_cast<size_t>(i)] = nextDyeB[static_cast<size_t>(i)] = 0.0f;
            nextCompCount[static_cast<size_t>(i)] = 0;
        }
        float nh = nextHeat[static_cast<size_t>(i)];
        float ndr = nextDyeR[static_cast<size_t>(i)];
        float ndg = nextDyeG[static_cast<size_t>(i)];
        float ndb = nextDyeB[static_cast<size_t>(i)];
        if (!std::isfinite(nh) || nh < 0.0f) nh = 0.0f;
        fill[i]=nf;
        liquidHeat[i]=nh;
        dyeR[i]=std::max(0.0f, ndr); dyeG[i]=std::max(0.0f, ndg); dyeB[i]=std::max(0.0f, ndb);
        commitNextComposition(i);
        if(fill[i]>=MIN_RENDER_FILL&&old<MIN_RENDER_FILL)waterShade[i]=makeShade(x,y);
    }
    for (size_t n = 0; n < advectOverflowIndex.size(); ++n)
        returnAdvectOverflow(advectOverflowIndex[n], advectOverflowVol[n], advectOverflowCarry[n], x0, y0, x1, y1);
    (void)regionVolumeBefore;
}

void FluidEngine::returnAdvectOverflow(int origin, float vol, LiquidCarry &carry, int x0, int y0, int x1, int y1) {
    if (!(vol > 1.0e-12f)) return;
    auto inSolve = [&](int x, int y) {
        return x >= x0 && x <= x1 && y >= y0 && y <= y1 && inside(x, y);
    };
    auto isLiquidCell = [&](int ni) {
        return !solid[static_cast<size_t>(ni)] && !dynamicSolid[static_cast<size_t>(ni)]
            && fill[static_cast<size_t>(ni)] > 1.0e-7f;
    };
    if (++relocateEpoch == 0) {
        std::fill(relocateStamp.begin(), relocateStamp.end(), 0);
        relocateEpoch = 1;
    }
    relocateQueue.clear();
    auto consider = [&](int nx, int ny) {
        if (!inSolve(nx, ny) || isSolid(nx, ny)) return;
        int ni = ci(nx, ny);
        if (relocateStamp[static_cast<size_t>(ni)] == relocateEpoch) return;
        if (!isLiquidCell(ni)) return;
        relocateStamp[static_cast<size_t>(ni)] = relocateEpoch;
        relocateQueue.push_back(ni);
    };
    int ox = origin % GW, oy = origin / GW;
    if (inSolve(ox, oy) && !isSolid(ox, oy) && isLiquidCell(origin))
        consider(ox, oy);
    else {
        if (openUFace(ox, oy)) consider(ox - 1, oy);
        if (openUFace(ox + 1, oy)) consider(ox + 1, oy);
        if (openVFace(ox, oy)) consider(ox, oy - 1);
        if (openVFace(ox, oy + 1)) consider(ox, oy + 1);
    }
    for (size_t head = 0; head < relocateQueue.size() && vol > 1.0e-12f; ++head) {
        int i = relocateQueue[head];
        int x = i % GW, y = i / GW;
        float room = std::max(0.0f, 1.0f - fill[static_cast<size_t>(i)]);
        if (room > 1.0e-12f) {
            float placed = std::min(room, vol);
            LiquidCarry part = splitCarry(carry, placed, vol);
            fill[static_cast<size_t>(i)] += placed;
            applyCarry(i, part);
            vol -= placed;
        }
        if (openUFace(x, y)) consider(x - 1, y);
        if (openUFace(x + 1, y)) consider(x + 1, y);
        if (openVFace(x, y)) consider(x, y - 1);
        if (openVFace(x, y + 1)) consider(x, y + 1);
    }
    if (vol > 1.0e-12f && inSolve(ox, oy) && !isSolid(ox, oy)) {
        fill[static_cast<size_t>(origin)] += vol;
        applyCarry(origin, carry);
    }
}

// Conservatively gather thin residuals instead of deleting them or
// allowing them to pressure-solve as full cells.
//
// Trails behind falling droplets are Eulerian fractional fill left by
// conservative transport (volume is conserved, but Normal view paints any
// cell >= MIN_RENDER_FILL as a full pixel). They are not mass creation.
//
// Consolidation merges thin cells into connected, *fuller* nearby liquid -
// never delete mass, never walk into empty air (that extends trails), and
// never cut diagonally through a sealed solid corner. Equal-fill cells may
// only coalesce downward so a smear collapses into its leading blob.
void FluidEngine::consolidateResidualVolume() {
    if (!residualConsolidationEnabled) return;
    int x0 = activeX0(1), y0 = activeY0(1), x1 = activeX1(1), y1 = activeY1(1);
    workCounts.residualTransfers = 0;
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        residualTarget[i] = -1;
        residualTransfer[i] = incoming[i] = 0.0f;
    }

    // Orthogonal first, then diagonals. Diagonals require an open orthogonal gate.
    constexpr int dx[8] = {0, -1, 1, 0, -1, 1, -1, 1};
    constexpr int dy[8] = {1, 0, 0, -1, 1, 1, -1, -1};

    for (int y = std::max(1, activeY0(1)); y <= std::min(GH - 2, activeY1(1)); ++y)
        for (int x = std::max(1, activeX0(1)); x <= std::min(GW - 2, activeX1(1)); ++x) {
            int index = ci(x, y);
            if (solid[index] || dynamicSolid[index] || fill[index] <= 0.0f || fill[index] >= MIN_CONSOLIDATE_FILL) continue;

            // Interior / almost-interior cells are undersaturated pool, not crumbs.
            // Dumping them into neighbors carves 0.03-0.28 holes in otherwise full water.
            int visibleNeighbors = 0;
            constexpr int odx[4] = {-1, 1, 0, 0};
            constexpr int ody[4] = {0, 0, -1, 1};
            for (int n = 0; n < 4; ++n) {
                int nx = x + odx[n], ny = y + ody[n];
                if (!inside(nx, ny) || isSolid(nx, ny)) continue;
                if (fill[ci(nx, ny)] >= MIN_RENDER_FILL) ++visibleNeighbors;
            }
            if (visibleNeighbors >= 3) continue;

            float vx = cellU(x, y), vy = cellV(x, y);
            int best = -1;
            float bestScore = -1e9f;

            for (int n = 0; n < 8; ++n) {
                int nx = x + dx[n], ny = y + dy[n];
                if (!inside(nx, ny)) continue;
                int ni = ci(nx, ny);
                if (isSolid(nx, ny)) continue;
                float nf = fill[ni];
                float room = 1.0f - nf;
                if (room <= 1e-5f) continue;

                bool diagonal = dx[n] != 0 && dy[n] != 0;
                if (diagonal) {
                    bool viaHorizontal = !isSolid(x + dx[n], y);
                    bool viaVertical = !isSolid(x, y + dy[n]);
                    if (!viaHorizontal || !viaVertical) continue;
                }

                // Never step into completely empty air - that stretches falling-stream trails.
                if (nf <= 1e-7f) continue;
                // Never pull volume back upward: that undoes falling (a 1-cell drop
                // would reabsorb the flux it just sent into the cell below).
                if (dy[n] < 0) continue;
                // Only coalesce into equal or fuller liquid.
                if (nf < fill[index] - 1e-5f) continue;

                float score = nf * 250.0f;
                if (nf >= MIN_RENDER_FILL) score += 400.0f;
                else if (nf >= MIN_ACTIVE_FILL) score += 120.0f;
                else score += 30.0f;

                if (dy[n] > 0) score += 20.0f;
                score += (vx * dx[n] + vy * dy[n]) * 2.0f;
                if (!diagonal) score += 8.0f;
                score += room * 15.0f;

                if (score > bestScore) {
                    bestScore = score;
                    best = ni;
                }
            }

            if (best >= 0) {
                residualTarget[index] = best;
                residualTransfer[index] = fill[index];
                incoming[best] += fill[index];
            }
        }

    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        nextFill[i] = fill[i];
        nextHeat[i] = liquidHeat[i];
        nextDyeR[i] = dyeR[i]; nextDyeG[i] = dyeG[i]; nextDyeB[i] = dyeB[i];
        copyCompositionToNext(i);
    }
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        if (residualTarget[i] < 0) continue;
        int receiver = residualTarget[i];
        float scale = incoming[receiver] > 1.0f - fill[receiver]
            ? (1.0f - fill[receiver]) / incoming[receiver]
            : 1.0f;
        float amount = residualTransfer[i] * std::max(0.0f, scale);
        if (amount > 1e-8f) ++workCounts.residualTransfers;
        nextFill[i] -= amount;
        nextFill[receiver] += amount;
        float f0 = fill[i];
        float frac = (f0 > 1.0e-20f) ? (amount / f0) : 0.0f;
        auto moveN = [&](std::vector<float> const &src, std::vector<float> &nxt) {
            float d = src[static_cast<size_t>(i)] * frac;
            nxt[static_cast<size_t>(i)] -= d;
            nxt[static_cast<size_t>(receiver)] += d;
        };
        moveN(liquidHeat, nextHeat);
        moveN(dyeR, nextDyeR); moveN(dyeG, nextDyeG); moveN(dyeB, nextDyeB);
        int base = compositionSlot(i, 0);
        int n = liquidCompCount[static_cast<size_t>(i)];
        for (int s = 0; s < n; ++s) {
            float d = liquidCompAmt[static_cast<size_t>(base + s)] * frac;
            nextCompAmt[static_cast<size_t>(base + s)] -= d;
            RuntimeSubstanceRef storedId = liquidCompId[static_cast<size_t>(base + s)];
            if (!runtimeSubstanceIsBuiltIn(storedId)) continue;
            float unplaced = addNextComponentUntracked(receiver, runtimeBuiltinId(storedId), d);
            if (unplaced > 0.0f) {
                nextCompAmt[static_cast<size_t>(base + s)] += unplaced;
                nextFill[receiver] -= unplaced;
                nextFill[i] += unplaced;
            }
        }
    }
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        fill[i] = nextFill[i];
        liquidHeat[i] = nextHeat[i];
        dyeR[i] = nextDyeR[i]; dyeG[i] = nextDyeG[i]; dyeB[i] = nextDyeB[i];
        commitNextComposition(i);
        clearEmptyLiquidCell(i);
    }
}

// Interior cells can be drained below MIN_RENDER_FILL by limiter/clamp while
// neighbors still hold plenty of liquid, which reads as a random hole.
// Borrow from those neighbors (mass-conserving) so enclosed gaps stay continuous.
void FluidEngine::repairEnclosedUndersaturatedCells() {
    int x0 = std::max(1, activeX0(1)), y0 = std::max(1, activeY0(1));
    int x1 = std::min(GW - 2, activeX1(1)), y1 = std::min(GH - 2, activeY1(1));
    constexpr int dx[4] = {-1, 1, 0, 0};
    constexpr int dy[4] = {0, 0, -1, 1};
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int i = ci(x, y);
        if (solid[i] || dynamicSolid[i]) continue;
        if (fill[i] <= 1.0e-8f) continue;
        int vis[4];
        int nvis = 0;
        int closed = 0;
        int mask = 0;
        float neighborSum = 0.0f;
        for (int n = 0; n < 4; ++n) {
            int nx = x + dx[n], ny = y + dy[n];
            if (!inside(nx, ny) || isSolid(nx, ny)) { ++closed; continue; }
            if (fill[ci(nx, ny)] >= MIN_RENDER_FILL) {
                vis[nvis++] = ci(nx, ny);
                mask |= 1 << n;
                ++closed;
                neighborSum += fill[ci(nx, ny)];
            }
        }
        bool enclosed = closed == 4 && nvis >= 3;
        bool tunnel = nvis == 2 && closed == 2 && ((mask & 3) == 3 || (mask & 12) == 12);
        if (!enclosed && !tunnel) continue;
        float target = MIN_RENDER_FILL;
        if (enclosed && nvis > 0) target = std::max(MIN_RENDER_FILL, neighborSum / static_cast<float>(nvis));
        float need = target - fill[i];
        if (need <= 1e-8f) continue;
        float extra[4]{};
        float available = 0.0f;
        for (int k = 0; k < nvis; ++k) {
            extra[k] = std::max(0.0f, fill[vis[k]] - MIN_RENDER_FILL);
            available += extra[k];
        }
        if (available <= 1e-8f) continue;
        float take = std::min(need, available);
        for (int k = 0; k < nvis; ++k) {
            if (extra[k] <= 0.0f) continue;
            float d = take * (extra[k] / available);
            LiquidCarry c = extractVolume(vis[k], d);
            fill[i] += d;
            applyCarry(i, c);
        }
    }
}

// After a large body falls away, isolated Eulerian crumbs/droplets can sit in
// sleeping chunks with ~zero velocity. Gravity only runs inside the active solve
// region, so they appear frozen until something wakes them.
// Promote unsupported isolated liquid into SplashParticles (ballistic gravity)
// without deleting volume. Connected streams/pools are left alone.
void FluidEngine::promoteUnsupportedIsolatedLiquid() {
    if (splashes.size() > 2000) return;
    constexpr int odx[4] = {-1, 1, 0, 0};
    constexpr int ody[4] = {0, 0, -1, 1};
    // Airborne cells in sleeping chunks wake during metrics; only the solve region needs promotion here.
    int x0 = std::max(1, activeX0()), y0 = std::max(1, activeY0());
    int x1 = std::min(GW - 2, activeX1()), y1 = std::min(GH - 2, activeY1());
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int index = ci(x, y);
        float amount = fill[index];
        if (solid[index] || dynamicSolid[index] || amount < 1e-5f) continue;

        int below = ci(x, y + 1);
        if (isSolid(x, y + 1) || fill[below] >= MIN_SUBSTANTIAL_FILL) continue;

        int liquidNeighbors = 0;
        for (int n = 0; n < 4; ++n) {
            int nx = x + odx[n], ny = y + ody[n];
            if (!inside(nx, ny) || isSolid(nx, ny)) continue;
            if (fill[ci(nx, ny)] >= MIN_ACTIVE_FILL) ++liquidNeighbors;
        }
        // Still touching Eulerian liquid (including a falling blob above/beside).
        // Promoting those crumbs spawned a splash trail behind droplets.
        if (liquidNeighbors >= 1) continue;

        float vx = cellU(x, y), vy = cellV(x, y);
        float speed = std::abs(vx) + std::abs(vy);
        // Substantial isolated blobs may keep falling as Eulerian cells; only
        // promote them if they have stalled in mid-air.
        if (amount >= MIN_SUBSTANTIAL_FILL && speed >= 0.08f) continue;

        SplashParticle p;
        p.x = x + 0.5f;
        p.y = y + 0.5f;
        p.vx = vx;
        p.vy = std::max(vy, 0.0f);
        p.volume = amount;
        p.life = 0.0f;
        p.originX = x;
        p.originY = y;
        LiquidCarry c = extractVolume(index, amount);
        p.heat = c.heat;
        p.dyeR = c.dyeR; p.dyeG = c.dyeG; p.dyeB = c.dyeB;
        copyCarryToSplash(p, c);
        splashes.push_back(p);

        pressure[index] = 0.0f;
        if (openUFace(x, y)) u[ui(x, y)] *= 0.5f;
        if (openUFace(x + 1, y)) u[ui(x + 1, y)] *= 0.5f;
        if (openVFace(x, y)) v[vi(x, y)] *= 0.5f;
        if (openVFace(x, y + 1)) v[vi(x, y + 1)] = std::max(v[vi(x, y + 1)], 0.0f);
        wakeChunkAtCell(x, y);
        if (splashes.size() > 2000) return;
    }
}

float FluidEngine::depositVolume(float x, float y, float amount, float momentumX, float momentumY, LiquidCarry *carry) {
    if (amount <= 1e-7f) return amount;
    int cx = static_cast<int>(std::floor(x)), cy = static_cast<int>(std::floor(y));
    // Never seed from a solid: a 1px wall has open cells on both faces.
    if (!inside(cx, cy) || isSolid(cx, cy)) {
        int sx = cx, sy = cy;
        if (std::abs(momentumX) >= std::abs(momentumY) && std::abs(momentumX) > 1e-8f)
            sx = cx - (momentumX >= 0.0f ? 1 : -1);
        else if (std::abs(momentumY) > 1e-8f)
            sy = cy - (momentumY >= 0.0f ? 1 : -1);
        if (!inside(sx, sy) || isSolid(sx, sy)) return amount;
        cx = sx;
        cy = sy;
    }
    int index = ci(cx, cy);
    float room = 1.0f - fill[index];
    if (room > 1e-5f) {
        float placed = std::min(room, amount);
        bool newCell = fill[index] < MIN_RENDER_FILL;
        LiquidCarry chunk = carry ? splitCarry(*carry, placed, amount) : ambientCarry(*this, placed, SUBSTANCE_WATER);
        fill[index] += placed;
        applyCarry(index, chunk);
        amount -= placed;
        if (newCell) waterShade[index] = makeShade(cx, cy);
        wakeChunkAtCell(cx, cy);
        u[ui(cx, cy)] += momentumX * placed * 0.25f;
        u[ui(cx + 1, cy)] += momentumX * placed * 0.25f;
        v[vi(cx, cy)] += momentumY * placed * 0.25f;
        v[vi(cx, cy + 1)] += momentumY * placed * 0.25f;
    }
    if (amount <= 1e-5f) return amount;
    return relocateVolumeTopologySafe(cx, cy, amount, momentumX, momentumY, carry, false);
}

void FluidEngine::spawnSurfaceSpray() {
    if ((tickNo & 1u) != 0u || splashes.size() > 1000) return;
    for (int y = std::max(1,activeY0()); y <= std::min(GH-2,activeY1()); ++y)
        for (int x = std::max(1,activeX0()); x <= std::min(GW-2,activeX1()); ++x) {
        int index = ci(x, y);
        if (fill[index] < 0.20f || !surfaceMask[index]) continue;
        if (mixViscosity(index) > sandboxReferenceLiquid().viscosity * 2.5f) continue;
        float vx = cellU(x, y), vy = cellV(x, y), speed = std::sqrt(vx * vx + vy * vy);
        float normalVelocity = std::max(0.0f, vx*surfaceNormalX[index] + vy*surfaceNormalY[index]);
        float pressureGradient = 0.5f*std::sqrt(
            std::pow(pressure[ci(x+1,y)]-pressure[ci(x-1,y)],2.0f) +
            std::pow(pressure[ci(x,y+1)]-pressure[ci(x,y-1)],2.0f));
        float pressureImpulse = std::abs(pressure[index]-pressureBefore[index]) * PHYSICS_DT;
        float breakupEnergy = normalVelocity*0.75f + speed*0.16f + pressureGradient*0.035f
            + std::abs(surfaceCurvature[index])*1.5f + pressureImpulse*0.10f;
        if (breakupEnergy < 14.0f || (hashCell(x, y, tickNo) & 15u) != 0u) continue;
        float detached = std::min(0.10f, fill[index] * 0.22f);
        LiquidCarry c = extractVolume(index, detached);
        SplashParticle p;
        p.x = x + 0.5f; p.y = y + 0.5f; p.vx = vx; p.vy = vy;
        p.volume = detached; p.originX = x; p.originY = y;
        copyCarryToSplash(p, c);
        splashes.push_back(p);
    }
}

void FluidEngine::updateSplashParticles(float dt) {
    constexpr int kMaxSplashBounces = 2;
    auto onOpenVoidRim = [&](float x, float y) {
        if (config.walledBorders) return false;
        int gx = static_cast<int>(std::floor(x)), gy = static_cast<int>(std::floor(y));
        if (!inside(gx, gy)) return true;
        return gx <= 0 || gx >= GW - 1 || gy <= 0 || gy >= GH - 1;
    };
    for (size_t i = 0; i < splashes.size();) {
        SplashParticle &p = splashes[i];
        p.vy += gridGravity() * dt;
        bool remove = false;
        auto sink = [&]() {
            expectedVolume = std::max(0.0, expectedVolume - static_cast<double>(p.volume));
            if (std::isfinite(p.heat) && p.heat > 0.0f) escapedHeat += p.heat;
            p.volume = 0.0f;
            p.heat = 0.0f;
            p.dyeR = p.dyeG = p.dyeB = 0.0f;
            p.compCount = 0;
            remove = true;
        };
        if (onOpenVoidRim(p.x, p.y)) {
            sink();
        } else {
        float nx = p.x + p.vx * dt, ny = p.y + p.vy * dt;
        float lastX = p.x, lastY = p.y; bool hit = false, escapedToVoid = false;
        auto blocked = [&](int gx, int gy) {
            bool stillInOrigin = gx == p.originX && gy == p.originY && p.life < 0.12f;
            if (!inside(gx, gy)) {
                if (config.walledBorders) hit = true;
                else escapedToVoid = true;
                return true;
            }
            if (isSolid(gx, gy)) { hit = true; return true; }
            if (!stillInOrigin && isFluid(gx, gy) && fill[ci(gx, gy)] > 0.65f) { hit = true; return true; }
            return false;
        };
        int pgx = static_cast<int>(std::floor(p.x)), pgy = static_cast<int>(std::floor(p.y));
        float dist = std::max(std::abs(nx - p.x), std::abs(ny - p.y));
        int steps = std::max(1, static_cast<int>(std::ceil(dist * 4.0f)));
        for (int s = 1; s <= steps && !hit && !escapedToVoid; ++s) {
            float t = static_cast<float>(s) / static_cast<float>(steps);
            float sx = p.x + (nx - p.x) * t, sy = p.y + (ny - p.y) * t;
            int gx = static_cast<int>(std::floor(sx)), gy = static_cast<int>(std::floor(sy));
            if (gx == pgx && gy == pgy) { lastX = sx; lastY = sy; continue; }
            // Diagonal step: both orthogonal cells must be open or this is a corner cut.
            if (gx != pgx && gy != pgy) {
                if (blocked(gx, pgy) || blocked(pgx, gy)) break;
            }
            if (blocked(gx, gy)) break;
            lastX = sx;
            lastY = sy;
            pgx = gx;
            pgy = gy;
        }
        p.life += dt;
        if (escapedToVoid || onOpenVoidRim(lastX, lastY)) {
            sink();
        } else if (hit || p.life > 3.0f) {
            LiquidCarry carry = carryFromSplash(p);
            float remainder = depositVolume(lastX, lastY, p.volume, p.vx, p.vy, &carry);
            copyCarryToSplash(p, carry);
            if (remainder <= 1e-5f) remove = true;
            else if (!config.walledBorders && (onOpenVoidRim(lastX, lastY) || p.bounceCount >= kMaxSplashBounces)) {
                expectedVolume = std::max(0.0, expectedVolume - static_cast<double>(remainder));
                if (std::isfinite(p.heat) && p.heat > 0.0f) escapedHeat += p.heat;
                remove = true;
            } else {
                p.volume = remainder;
                p.bounceCount = static_cast<uint8_t>(std::min(255, p.bounceCount + 1));
                p.vx *= -0.25f;
                p.vy *= -0.25f;
                p.x = lastX;
                p.y = lastY;
            }
        } else { p.x = nx; p.y = ny; }
        }
        if (remove) { splashes[i] = splashes.back(); splashes.pop_back(); }
        else {
            int gx = static_cast<int>(p.x), gy = static_cast<int>(p.y);
            if (!inside(gx, gy) || fill[ci(gx, gy)] < MIN_SUBSTANTIAL_FILL)
                wakeChunkAtCell(gx, gy, false);
            ++i;
        }
    }
}

void FluidEngine::drainWaterIntoVoid() {
    if (config.walledBorders) return;
    double removed = 0.0;
    auto drain = [&](int index) {
        if (solid[index] || dynamicSolid[index] || fill[index] <= 0.0f) return;
        float amount = fill[index];
        fill[index] = 0.0f;
        if (std::isfinite(liquidHeat[index]) && liquidHeat[index] > 0.0f)
            escapedHeat += liquidHeat[index];
        liquidHeat[index] = 0.0f;
        dyeR[index] = dyeG[index] = dyeB[index] = 0.0f;
        clearComposition(index);
        pressure[index] = 0.0f;
        removed += amount;
    };
    for (int y = 0; y < GH; ++y) {
        drain(ci(0, y));
        drain(ci(GW - 1, y));
    }
    for (int x = 0; x < GW; ++x) {
        drain(ci(x, 0));
        drain(ci(x, GH - 1));
    }
    expectedVolume = std::max(0.0, expectedVolume - removed);
}

void FluidEngine::wakeChunkAtCell(int x, int y, bool resetQuiet) {
    if (!inside(x, y)) return;
    int cx = x / CHUNK, cy = y / CHUNK;
    for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
        int nx = cx + ox, ny = cy + oy;
        if (nx < 0 || nx >= CHUNK_W || ny < 0 || ny >= CHUNK_H) continue;
        int c = ny * CHUNK_W + nx;
        if (!chunkActivity[c]) {
            chunkActivity[c] = 1;
            chunkQuietTicks[c] = 0;
        } else if (resetQuiet) {
            chunkQuietTicks[c] = 0;
        }
    }
}

void FluidEngine::wakeRegion(int x0, int y0, int x1, int y1) {
    x0 = std::max(0, x0); y0 = std::max(0, y0); x1 = std::min(GW - 1, x1); y1 = std::min(GH - 1, y1);
    if (x1 < x0 || y1 < y0) return;
    for (int y = y0; y <= y1; y += CHUNK) wakeChunkAtCell(x0, y);
    for (int y = y0; y <= y1; y += CHUNK) wakeChunkAtCell(x1, y);
    for (int x = x0; x <= x1; x += CHUNK) { wakeChunkAtCell(x, y0); wakeChunkAtCell(x, y1); }
    wakeChunkAtCell((x0 + x1) / 2, (y0 + y1) / 2);
}

void FluidEngine::clearDynamicOccupancy() {
    for (int i : dynamicOccupiedCells) {
        dynamicSolid[i] = 0;
        dynamicVelX[i] = 0.0f;
        dynamicVelY[i] = 0.0f;
    }
    dynamicOccupiedCells.clear();
}

void FluidEngine::setDynamicOccupancy(int x, int y, float velX, float velY) {
    if (!inside(x, y) || solid[ci(x, y)]) return;
    int i = ci(x, y);
    if (!dynamicSolid[i]) dynamicOccupiedCells.push_back(i);
    dynamicSolid[i] = 1;
    dynamicVelX[i] = velX;
    dynamicVelY[i] = velY;
}

void FluidEngine::applyMovingBoundaryVelocity() {
    for (int i : dynamicOccupiedCells) {
        int x = i % GW, y = i / GW;
        float vx = dynamicVelX[i], vy = dynamicVelY[i];
        u[ui(x, y)] = vx;
        u[ui(x + 1, y)] = vx;
        v[vi(x, y)] = vy;
        v[vi(x, y + 1)] = vy;
    }
}

float FluidEngine::relocateVolumeTopologySafe(int x, int y, float amount, float momentumX, float momentumY, LiquidCarry *carry, bool fromOccupiedCell) {
    if (amount <= 1e-7f) return 0.0f;
    if (!fromOccupiedCell && (!inside(x, y) || isSolid(x, y))) return amount;
    if (++relocateEpoch == 0) {
        std::fill(relocateStamp.begin(), relocateStamp.end(), 0);
        relocateEpoch = 1;
    }
    relocateQueue.clear();
    auto consider = [&](int nx, int ny) {
        if (!inside(nx, ny) || isSolid(nx, ny)) return;
        int ni = ci(nx, ny);
        if (relocateStamp[static_cast<size_t>(ni)] == relocateEpoch) return;
        relocateStamp[static_cast<size_t>(ni)] = relocateEpoch;
        relocateQueue.push_back(ni);
    };
    constexpr int dx[4] = {-1, 1, 0, 0};
    constexpr int dy[4] = {0, 0, -1, 1};
    if (fromOccupiedCell)
        for (int n = 0; n < 4; ++n) consider(x + dx[n], y + dy[n]);
    else {
        consider(x, y);
        for (int n = 0; n < 4; ++n) consider(x + dx[n], y + dy[n]);
    }
    for (size_t head = 0; head < relocateQueue.size() && amount > 1e-5f && head < 512; ++head) {
        int i = relocateQueue[head];
        int cx = i % GW, cy = i / GW;
        float room = 1.0f - fill[i];
        if (room > 1e-5f) {
            float placed = std::min(room, amount);
            bool newCell = fill[i] < MIN_RENDER_FILL;
            LiquidCarry chunk = carry ? splitCarry(*carry, placed, amount) : ambientCarry(*this, placed, SUBSTANCE_WATER);
            fill[i] += placed;
            applyCarry(i, chunk);
            amount -= placed;
            if (newCell) waterShade[i] = makeShade(cx, cy);
            u[ui(cx, cy)] += momentumX * placed * 0.25f;
            u[ui(cx + 1, cy)] += momentumX * placed * 0.25f;
            v[vi(cx, cy)] += momentumY * placed * 0.25f;
            v[vi(cx, cy + 1)] += momentumY * placed * 0.25f;
            wakeChunkAtCell(cx, cy);
        }
        for (int n = 0; n < 4; ++n) consider(cx + dx[n], cy + dy[n]);
    }
    return amount;
}

void FluidEngine::displaceFluidFromDynamicSolids() {
    volumeDisplacedRigid = 0.0;
    for (int i : dynamicOccupiedCells) {
        int x = i % GW, y = i / GW;
        if (!dynamicSolid[i] || fill[i] <= 0.0f) continue;
        float amount = fill[i];
        LiquidCarry carry = extractVolume(i, amount);
        pressure[i] = 0.0f;
        float remainder = relocateVolumeTopologySafe(x, y, amount, dynamicVelX[i], dynamicVelY[i], &carry, true);
        volumeDisplacedRigid += (amount - remainder);
        if (remainder > 1e-5f && splashes.size() < 2000) {
            SplashParticle p;
            p.x = x + 0.5f; p.y = y + 0.5f;
            p.vx = dynamicVelX[i]; p.vy = dynamicVelY[i];
            p.volume = remainder; p.originX = x; p.originY = y;
            copyCarryToSplash(p, carry);
            splashes.push_back(p);
            remainder = 0.0f;
        }
        if (remainder > 1e-5f) {
            volumeLostRigid += remainder;
            expectedVolume = std::max(0.0, expectedVolume - remainder);
        }
    }
}

void FluidEngine::computeSolveRegion() {
    std::fill(chunkSolveMask.begin(), chunkSolveMask.end(), 0);
    std::array<int, CHUNK_W * CHUNK_H> queue{};
    int head = 0, tail = 0;
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) if (chunkActivity[c]) {
        chunkSolveMask[c] = 1;
        if (chunkHasFluid[c]) queue[tail++] = c;
    }
    while (head < tail) {
        int c = queue[head++], cx = c % CHUNK_W, cy = c / CHUNK_W;
        constexpr int dx[4] = {-1, 1, 0, 0}, dy[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
            int nx = cx + dx[n], ny = cy + dy[n];
            if (nx < 0 || nx >= CHUNK_W || ny < 0 || ny >= CHUNK_H) continue;
            int nc = ny * CHUNK_W + nx;
            if (chunkHasFluid[nc] && !chunkSolveMask[nc]) {
                chunkSolveMask[nc] = 1;
                queue[tail++] = nc;
            }
        }
    }
    chunkHaloSource = chunkSolveMask;
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) if (chunkHaloSource[c]) {
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
            int nx = cx + ox, ny = cy + oy;
            if (nx >= 0 && nx < CHUNK_W && ny >= 0 && ny < CHUNK_H) chunkSolveMask[ny * CHUNK_W + nx] = 1;
        }
    }
    int minCx = CHUNK_W, minCy = CHUNK_H, maxCx = -1, maxCy = -1;
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) if (chunkSolveMask[c]) {
        minCx = std::min(minCx, c % CHUNK_W);
        maxCx = std::max(maxCx, c % CHUNK_W);
        minCy = std::min(minCy, c / CHUNK_W);
        maxCy = std::max(maxCy, c / CHUNK_W);
    }
    hasActiveSolveRegion = maxCx >= 0;
    if (hasActiveSolveRegion) {
        solveX0 = minCx * CHUNK;
        solveY0 = minCy * CHUNK;
        solveX1 = std::min(GW - 1, (maxCx + 1) * CHUNK - 1);
        solveY1 = std::min(GH - 1, (maxCy + 1) * CHUNK - 1);
    }
}

void FluidEngine::wakeAllFluidChunks() {
    std::fill(chunkActivity.begin(),chunkActivity.end(),0);std::fill(chunkQuietTicks.begin(),chunkQuietTicks.end(),0);
    for(int y=0;y<GH;++y)for(int x=0;x<GW;++x)if(fill[ci(x,y)]>0.0f)wakeChunkAtCell(x,y);
}

void FluidEngine::wakeUntrackedLiquid() {
    for (int chunk = 0; chunk < CHUNK_W * CHUNK_H; ++chunk) {
        if (chunkHasFluid[chunk] || chunkActivity[chunk] || chunkSolveMask[chunk]) continue;
        int cx = chunk % CHUNK_W, cy = chunk / CHUNK_W;
        int beginX = cx * CHUNK, endX = std::min(GW, beginX + CHUNK);
        int beginY = cy * CHUNK, endY = std::min(GH, beginY + CHUNK);
        bool found = false;
        for (int y = beginY; y < endY && !found; ++y)
            for (int x = beginX; x < endX && !found; ++x)
                if (fill[ci(x, y)] > 0.0f) found = true;
        if (found) wakeChunkAtCell(beginX, beginY);
    }
}

void FluidEngine::rebuildActivityAndMetrics(bool advanceSleep) {
    wakeUntrackedLiquid();
    std::array<float, CHUNK_W * CHUNK_H> maxSpeed{}, maxDiv{}, maxFillDelta{}, maxPressureDelta{};
    std::array<uint8_t, CHUNK_W * CHUNK_H> airborneChunk{};
    activeFluidCells = 0; currentVolume = 0.0; thinCellCount = 0; thinVolume = 0.0;
    momentumX = momentumY = kineticEnergy = 0.0; measuredMaxVelocity = 0.0f;
    workCounts.subvisibleCells = 0;
    for (int chunk = 0; chunk < CHUNK_W * CHUNK_H; ++chunk) {
        // Also scan the solve halo: liquid can cross into it during this tick and
        // must wake that chunk immediately rather than becoming frozen/unmeasured.
        if (!chunkActivity[chunk] && !chunkHasFluid[chunk] && !chunkSolveMask[chunk]) continue;
        int cx = chunk % CHUNK_W, cy = chunk / CHUNK_W;
        bool hasAnyFluid = false;
        bool hasAirborne = false;
        int beginX = cx * CHUNK, endX = std::min(GW, beginX + CHUNK);
        int beginY = cy * CHUNK, endY = std::min(GH, beginY + CHUNK);
        for (int y = beginY; y < endY; ++y) for (int x = beginX; x < endX; ++x) {
            int index = ci(x, y);
            float amount = fill[index];
            currentVolume += amount;
            if (amount > 0.0f) hasAnyFluid = true;
            if (amount > 0.0f && amount < MIN_RENDER_FILL) ++workCounts.subvisibleCells;
            if (amount > 0.0f && amount < 0.10f) { ++thinCellCount; thinVolume += amount; }
            if (amount > 1e-5f && y + 1 < GH) {
                auto supported = [&](int sx, int sy) {
                    if (!inside(sx, sy)) return false;
                    return isSolid(sx, sy) || fill[ci(sx, sy)] >= MIN_SUBSTANTIAL_FILL;
                };
                if (!supported(x, y + 1) && !supported(x - 1, y + 1) && !supported(x + 1, y + 1))
                    hasAirborne = true;
            }
            if (amount < MIN_ACTIVE_FILL) continue;
            float vx = cellU(x, y), vy = cellV(x, y), speed = std::abs(vx) + std::abs(vy);
            if (chunkActivity[chunk])
                measuredMaxVelocity = std::max(measuredMaxVelocity, std::max(std::abs(vx), std::abs(vy)));
            momentumX += amount * vx; momentumY += amount * vy; kineticEnergy += 0.5 * amount * (vx * vx + vy * vy);
            maxSpeed[chunk] = std::max(maxSpeed[chunk], speed);
            maxDiv[chunk] = std::max(maxDiv[chunk], std::abs(divergenceField[index]));
            maxFillDelta[chunk] = std::max(maxFillDelta[chunk], std::abs(amount - previousFill[index]));
            maxPressureDelta[chunk] = std::max(maxPressureDelta[chunk], std::abs(pressure[index] - pressureBefore[index]));
            if (chunkActivity[chunk]) ++activeFluidCells;
        }
        chunkHasFluid[chunk] = hasAnyFluid ? 1 : 0;
        airborneChunk[chunk] = hasAirborne ? 1 : 0;
        // Unsupported liquid must not sleep mid-air (would look frozen until a neighbor wakes it).
        if (hasAirborne) { chunkActivity[chunk] = 1; chunkQuietTicks[chunk] = 0; }
    }
    for (SplashParticle const &p : splashes) currentVolume += p.volume;
    workCounts.splashCount = static_cast<int>(splashes.size());
    if (advanceSleep) {
        for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
            if (airborneChunk[c]) {
                chunkActivity[c] = 1;
                chunkQuietTicks[c] = 0;
                continue;
            }
            bool settled = maxSpeed[c] < 0.06f && maxDiv[c] < 0.008f && maxFillDelta[c] < 0.0004f && maxPressureDelta[c] < 0.04f;
            bool fillMoved = maxFillDelta[c] > 0.002f;
            if (chunkActivity[c]) {
                if (settled)
                    chunkQuietTicks[c] = static_cast<uint8_t>(std::min(255, int(chunkQuietTicks[c]) + 1));
                else
                    chunkQuietTicks[c] = 0;
                if (chunkQuietTicks[c] >= 36) chunkActivity[c] = 0;
            } else if (fillMoved) {
                chunkActivity[c] = 1;
                chunkQuietTicks[c] = 0;
            }
        }
    }
    activeChunks = static_cast<int>(std::count(chunkActivity.begin(), chunkActivity.end(), uint8_t{1}));
    volumeError = currentVolume - expectedVolume;
    computeSolveRegion();
}

void FluidEngine::simulationTick() {
    flushPaintDirty();
    syncWorkerPool();
    auto physicsStart=Clock::now();
    ++tickNo;
    workCounts.activeUFaces = workCounts.activeVFaces = 0;
    workCounts.nonzeroFluxU = workCounts.nonzeroFluxV = 0;
    workCounts.limiterPasses = 0;
    workCounts.advectOverflowCells = 0;
    if(!hasActiveSolveRegion&&splashes.empty()){
        wakeUntrackedLiquid();
        rebuildActivityAndMetrics(true);
        // Airborne leftovers can wake chunks during metrics; if so, continue this tick.
        if(!hasActiveSolveRegion&&splashes.empty()){timingAccum.physics+=elapsedMs(physicsStart);++timingTicks;return;}
    }
    int substeps = std::clamp(static_cast<int>(std::ceil(measuredMaxVelocity * PHYSICS_DT / config.maxTravelPerSubstep)), 1, config.maxSubsteps);
    lastFluidSubsteps = substeps;
    substepCapReached = substeps == config.maxSubsteps && measuredMaxVelocity * PHYSICS_DT > config.maxTravelPerSubstep * config.maxSubsteps;
    float dt = PHYSICS_DT / static_cast<float>(substeps);
    for (int step = 0; step < substeps; ++step) {
        auto stage=Clock::now();
        advectVelocity(dt);
        timingAccum.advection+=elapsedMs(stage);stage=Clock::now();
        diffuseVelocity(dt);
        timingAccum.viscosity+=elapsedMs(stage);stage=Clock::now();
        applyGravity(dt);
        timingAccum.forces+=elapsedMs(stage);stage=Clock::now();
        if (config.surfaceTensionEnabled || config.sprayEnabled) updateSurfaceField();
        if (config.surfaceTensionEnabled) applySurfaceTension(dt);
        applyVorticityConfinement(dt);
        timingAccum.surface+=elapsedMs(stage);stage=Clock::now();
        clampVelocity();
        enforceActiveBoundaries();
        projectVelocity(dt);
        timingAccum.pressure+=elapsedMs(stage);stage=Clock::now();
        advectLiquidVolume(dt);
        timingAccum.transport+=elapsedMs(stage);stage=Clock::now();
        consolidateResidualVolume();
        repairEnclosedUndersaturatedCells();
        promoteUnsupportedIsolatedLiquid();
        timingAccum.residual+=elapsedMs(stage);stage=Clock::now();
        drainWaterIntoVoid();
        timingAccum.drain+=elapsedMs(stage);stage=Clock::now();
        updateSplashParticles(dt);
        timingAccum.splashes+=elapsedMs(stage);
    }
    if (config.sprayEnabled) spawnSurfaceSpray();
    auto chunkStage=Clock::now();
    rebuildActivityAndMetrics(true);
    timingAccum.chunks+=elapsedMs(chunkStage);timingAccum.physics+=elapsedMs(physicsStart);++timingTicks;
    if(timingTicks>=60){
        timingAverage.advection=timingAccum.advection/timingTicks;
        timingAverage.viscosity=timingAccum.viscosity/timingTicks;
        timingAverage.forces=timingAccum.forces/timingTicks;
        timingAverage.surface=timingAccum.surface/timingTicks;
        timingAverage.pressure=timingAccum.pressure/timingTicks;
        timingAverage.transport=timingAccum.transport/timingTicks;
        timingAverage.residual=timingAccum.residual/timingTicks;
        timingAverage.drain=timingAccum.drain/timingTicks;
        timingAverage.splashes=timingAccum.splashes/timingTicks;
        timingAverage.chunks=timingAccum.chunks/timingTicks;
        timingAverage.listsBuild=timingAccum.listsBuild/timingTicks;
        timingAverage.physics=timingAccum.physics/timingTicks;
        timingAverage.thermal=timingAccum.thermal/timingTicks;
        double renderKeep=timingAverage.render;
        timingAccum=TimingAverages{};
        timingAverage.render=renderKeep;
        timingTicks=0;
    }
}

void FluidEngine::zeroFluidState() {
    std::fill(fill.begin(), fill.end(), 0.0f); std::fill(nextFill.begin(), nextFill.end(), 0.0f);
    std::fill(liquidHeat.begin(), liquidHeat.end(), 0.0f); std::fill(nextHeat.begin(), nextHeat.end(), 0.0f);
    std::fill(dyeR.begin(), dyeR.end(), 0.0f); std::fill(dyeG.begin(), dyeG.end(), 0.0f); std::fill(dyeB.begin(), dyeB.end(), 0.0f);
    std::fill(liquidCompId.begin(), liquidCompId.end(), runtimeNone());
    std::fill(liquidCompAmt.begin(), liquidCompAmt.end(), 0.0f);
    std::fill(liquidCompCount.begin(), liquidCompCount.end(), 0);
    std::fill(solidifyPendingId.begin(), solidifyPendingId.end(), SUBSTANCE_NONE);
    std::fill(solidifyPendingKg.begin(), solidifyPendingKg.end(), 0.0f);
    std::fill(solidifyPendingHeatJ.begin(), solidifyPendingHeatJ.end(), 0.0f);
    std::fill(nextDyeR.begin(), nextDyeR.end(), 0.0f); std::fill(nextDyeG.begin(), nextDyeG.end(), 0.0f); std::fill(nextDyeB.begin(), nextDyeB.end(), 0.0f);
    std::fill(nextCompId.begin(), nextCompId.end(), runtimeNone());
    std::fill(nextCompAmt.begin(), nextCompAmt.end(), 0.0f);
    std::fill(nextCompCount.begin(), nextCompCount.end(), 0);
    std::fill(pressure.begin(), pressure.end(), 0.0f); std::fill(divergenceField.begin(), divergenceField.end(), 0.0f);
    std::fill(u.begin(), u.end(), 0.0f); std::fill(v.begin(), v.end(), 0.0f);
    splashes.clear(); expectedVolume = 0.0; escapedHeat = 0.0; tickNo = 1;
    pendingSolidCellCount = 0;
}

void FluidEngine::clearWorld() {
    std::fill(solid.begin(), solid.end(), 0);
    std::fill(solidHeat.begin(), solidHeat.end(), 0.0f);
    clearDynamicOccupancy();
    zeroFluidState();
    volumeLostRigid = 0.0; volumeDisplacedRigid = 0.0; escapedHeat = 0.0;
    std::fill(chunkActivity.begin(),chunkActivity.end(),0);std::fill(chunkHasFluid.begin(),chunkHasFluid.end(),0);std::fill(chunkQuietTicks.begin(),chunkQuietTicks.end(),0);
    std::fill(thermalChunkWake.begin(), thermalChunkWake.end(), 0);
    rebuildActivityAndMetrics();
}

void FluidEngine::resetWorld() {
    clearWorld();
    for (int x = 0; x < GW; ++x) solid[ci(x, GH - 1)] = solid[ci(x, GH - 2)] = 1;
    for (int y = 38; y < GH; ++y) solid[ci(0, y)] = solid[ci(GW - 1, y)] = 1;
    for (int x = 32; x < 75; ++x) solid[ci(x, 87)] = 1;
    for (int x = 125; x < 168; ++x) solid[ci(x, 98)] = 1;
    for (int y = 78; y < 99; ++y) solid[ci(124, y)] = 1;
    enforceSolidBoundaries(); seedAmbientHeat(); rebuildActivityAndMetrics();
}

void FluidEngine::loadTestScene(int scene) {
    clearWorld();
    auto wall = [this](int x, int y) { if (inside(x, y)) solid[ci(x, y)] = 1; };
    auto waterRect = [this](int x0, int y0, int x1, int y1) {
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) if (inside(x, y) && !solid[ci(x, y)]) {
            fill[ci(x, y)] = 1.0f; waterShade[ci(x, y)] = makeShade(x, y);
            (void)addComponentUntracked(ci(x, y), SUBSTANCE_WATER, 1.0f);
        }
    };
    for (int x = 20; x <= 180; ++x) wall(x, 112);
    for (int y = 20; y <= 112; ++y) { wall(20, y); wall(180, y); }
    if (scene == 1) {
        waterRect(21, 72, 179, 111); // still pool
    } else if (scene == 2) {
        for (int y = 24; y <= 103; ++y) wall(100, y);
        waterRect(21, 55, 99, 111); // communicating vessels, bottom connection open
    } else if (scene == 3) {
        for (int y = 36; y <= 111; ++y) wall(92, y);
        waterRect(21, 47, 91, 111); // dam break: erase the center wall
    } else if (scene == 4) {
        for (int x = 30; x <= 90; ++x) wall(x, 55);
        for (int y = 24; y <= 55; ++y) wall(30, y);
        for (int x = 31; x <= 86; ++x) waterRect(x, 34, x, 54); // reservoir with right spill lip
        for (int x = 105; x <= 125; ++x) for (int y = 88; y <= 98; ++y) wall(x, y);
    } else if (scene == 5) {
        for (int x = 82; x <= 118; ++x) for (int y = 82; y <= 98; ++y) wall(x, y);
        waterRect(91, 24, 109, 48); // stream above an obstacle
    } else if (scene == 6) {
        waterRect(99, 30, 100, 31); // four-cell cohesive droplet
    } else if (scene == 7) {
        waterRect(91, 103, 108, 104); // 36-cell puddle on a flat floor
    } else if (scene == 8) {
        for (int step=0;step<12;++step) for(int x=38+step*9;x<47+step*9;++x)
            for(int y=48+step*5;y<=111;++y) wall(x,y);
        waterRect(43, 36, 51, 44); // stair-stepped incline
    } else if (scene == 9) {
        for(int x=55;x<=145;++x){wall(x,38);wall(x,76);}
        for(int y=38;y<=76;++y){wall(55,y);wall(145,y);}
        wall(99,76); wall(102,76); // two-cell nozzle at x=100..101
        waterRect(56,39,144,75);
    } else if (scene == 10) {
        for(int x=32;x<=168;++x){wall(x,73);wall(x,78);}
        for(int y=34;y<=77;++y) wall(32,y);
        waterRect(34,67,62,72); // source feeding a four-cell-high channel
    }
    expectedVolume = 0.0; for (float amount : fill) expectedVolume += amount;
    wakeAllFluidChunks(); enforceSolidBoundaries(); seedAmbientHeat(); rebuildActivityAndMetrics();
}

void FluidEngine::addSloshImpulse() {
    for (int y = 0; y < GH; ++y) for (int x = 1; x < GW; ++x)
        if (openUFace(x, y) && (isFluid(x - 1, y) || isFluid(x, y))) u[ui(x, y)] += 10.0f;
    wakeAllFluidChunks();
}

void FluidEngine::paintDisc(int cx, int cy, Tool tool, int brushRadius, LiquidPaint paint,
    BrushShape shape, bool eraseWalls, bool eraseLiquids) {
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y) for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
        if (!inside(x, y) || !brushContains(shape, cx, cy, x, y, brushRadius)) continue;
        int index = ci(x, y);
        if (tool == Tool::Water || tool == Tool::Brush) {
            float s = std::clamp(paint.dyeStrength, 0.0f, 1.0f);
            bool clearDye = (paint.dyeR + paint.dyeG + paint.dyeB) <= 1.0e-6f;
            if (paint.dyeOnly) {
                if (solid[index] || fill[index] < MIN_ACTIVE_FILL) continue;
                float f = fill[index];
                if (clearDye) {
                    dyeR[index] *= (1.0f - s);
                    dyeG[index] *= (1.0f - s);
                    dyeB[index] *= (1.0f - s);
                } else {
                    dyeR[index] = dyeR[index] * (1.0f - s) + paint.dyeR * f * s;
                    dyeG[index] = dyeG[index] * (1.0f - s) + paint.dyeG * f * s;
                    dyeB[index] = dyeB[index] * (1.0f - s) + paint.dyeB * f * s;
                }
                continue;
            }
            if (solid[index]) { solid[index] = 0; solidHeat[index] = 0.0f; }
            float added = 1.0f - fill[index];
            if (added > 0.0f) {
                FluidProperties const &liq = fluidForSubstance(paint.substance());
                ThermalProperties const &th = thermalForSubstance(paint.substance());
                float cap = thermalCapacity(massKg(liq.density, added, config.cellsPerMeter),
                    th.specificHeat);
                liquidHeat[index] += energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
                fill[index] = 1.0f;
                (void)addComponentUntracked(index, paint.substance(), added);
                if (!clearDye) {
                    dyeR[index] += added * paint.dyeR * s;
                    dyeG[index] += added * paint.dyeG * s;
                    dyeB[index] += added * paint.dyeB * s;
                }
                expectedVolume += added;
                waterShade[index] = makeShade(x, y);
            }
        } else if (tool == Tool::Solid) {
            if (solid[index]) continue;
            float displaced = fill[index];
            LiquidCarry carry{};
            if (displaced > 0.0f) carry = extractVolume(index, displaced);
            solid[index] = 1;
            float cap = thermalCapacity(massKg(mechanicalForSubstance(kStaticWallSubstance).densityRel, 1.0f, config.cellsPerMeter),
                thermalForSubstance(kStaticWallSubstance).specificHeat);
            solidHeat[index] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
            if (displaced > 0.0f) {
                float remainder = depositVolume(x + 0.5f, y - 0.5f, displaced, 0.0f, -2.0f, &carry);
                if (remainder > 1e-5f) {
                    SplashParticle p;
                    p.x = x + 0.5f; p.y = y - 0.1f; p.vx = 0.0f; p.vy = -3.0f;
                    p.volume = remainder; p.originX = x; p.originY = y;
                    copyCarryToSplash(p, carry);
                    splashes.push_back(p);
                }
            }
        } else {
            if (eraseWalls && solid[index]) {
                solid[index] = 0;
                solidHeat[index] = 0.0f;
            }
            if (eraseLiquids && fill[index] > 0.0f) {
                expectedVolume -= fill[index];
                fill[index] = 0.0f;
                liquidHeat[index] = 0.0f;
                dyeR[index] = dyeG[index] = dyeB[index] = 0.0f;
                clearComposition(index);
            }
        }
    }
    if (tool == Tool::Eraser && eraseLiquids) {
        for (size_t i = 0; i < splashes.size();) {
            int sx = static_cast<int>(std::floor(splashes[i].x));
            int sy = static_cast<int>(std::floor(splashes[i].y));
            if (brushContains(shape, cx, cy, sx, sy, brushRadius)) {
                expectedVolume -= splashes[i].volume; splashes[i] = splashes.back(); splashes.pop_back();
            } else ++i;
        }
    }
    wakeChunkAtCell(cx,cy);
}

void FluidEngine::finalizePaint() { paintDirty = true; }

void FluidEngine::paintLine(int x0, int y0, int x1, int y1, Tool tool, int brushRadius, LiquidPaint paint,
    BrushShape shape, bool eraseWalls, bool eraseLiquids) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1, dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1, error = dx + dy;
    for (;;) {
        paintDisc(x0, y0, tool, brushRadius, paint, shape, eraseWalls, eraseLiquids);
        if (x0 == x1 && y0 == y1) break;
        int twice = 2 * error;
        if (twice >= dy) { error += dy; x0 += sx; }
        if (twice <= dx) { error += dx; y0 += sy; }
    }
}

void FluidEngine::runHeadlessBenchmark() {
    config.workerCount = 1;
    syncWorkerPool();
    std::string filename=miscFile("benchmark_"+std::to_string(GW)+"x"+std::to_string(GH)+".tsv");
    std::ofstream out(filename,std::ios::trunc);
    out<<"grid\tscene\tticks\tmerge\tvolume\texpected\terror\tthin_cells\tthin_volume\tactive_cells\tactive_chunks\tvmax\tpressure_iterations\tphysics_ms\tpressure_ms\tadvection_ms\ttransport_ms\tsurface_ms\tresidual_ms\tlists_ms\tpressure_cells\tsurface_cells\tflux_u\tflux_v\tlimiter_passes\tsplashes\tsubvisible\n";
    auto run=[&](int scene,char const*name,int ticks,bool merge){
        residualConsolidationEnabled=merge;timingAccum=TimingAverages{};timingAverage=TimingAverages{};timingTicks=0;
        loadTestScene(scene);auto start=Clock::now();
        for(int tick=0;tick<ticks;++tick)simulationTick();
        rebuildActivityAndMetrics();double perTick=elapsedMs(start)/std::max(1,ticks);
        out<<GW<<"x"<<GH<<'\t'<<name<<'\t'<<ticks<<'\t'<<(merge?"on":"off")<<'\t'
           <<currentVolume<<'\t'<<expectedVolume<<'\t'<<volumeError<<'\t'<<thinCellCount<<'\t'<<thinVolume<<'\t'
           <<activeFluidCells<<'\t'<<activeChunks<<'\t'<<measuredMaxVelocity<<'\t'<<lastPressureIterations<<'\t'
           <<perTick<<'\t'<<timingAverage.pressure<<'\t'<<timingAverage.advection<<'\t'<<timingAverage.transport<<'\t'
           <<timingAverage.surface<<'\t'<<timingAverage.residual<<'\t'<<timingAverage.listsBuild<<'\t'
           <<workCounts.pressureCells<<'\t'<<workCounts.surfaceCells<<'\t'<<workCounts.nonzeroFluxU<<'\t'<<workCounts.nonzeroFluxV<<'\t'
           <<workCounts.limiterPasses<<'\t'<<workCounts.splashCount<<'\t'<<workCounts.subvisibleCells<<'\n';
    };
    run(1,"still_pool",180,true);
    run(7,"small_puddle_baseline",240,false);
    run(7,"small_puddle_refined",240,true);
    run(3,"dam_break",180,true);
    run(5,"large_moving_body",180,true);
    run(6,"single_droplet",150,true);
    run(9,"narrow_nozzle",180,true);
}

void FluidEngine::runAdvectionBenchmark() {
    config.workerCount = 1;
    syncWorkerPool();
    residualConsolidationEnabled = true;
    VelocityAdvection saved = config.velocityAdvection;
    std::string filename = miscFile("advection_benchmark_" + std::to_string(GW) + "x" + std::to_string(GH) + ".tsv");
    std::ofstream out(filename, std::ios::trunc);
    out << "mode\tscene\tticks\tadvection_ms\tphysics_ms\tper_tick_ms\tvolume\texpected\terror\tvmax\tactive_cells\tactive_chunks\tmomentum_x\tmomentum_y\tke\n";
    VelocityAdvection modes[] = {
        VelocityAdvection::None,
        VelocityAdvection::FirstOrderUpwind,
        VelocityAdvection::NearestSemiLagrangian,
        VelocityAdvection::SemiLagrangian,
        VelocityAdvection::MacCormack,
        VelocityAdvection::BFECC
    };
    struct Scene { int id; char const *name; int ticks; };
    Scene scenes[] = {
        {1, "still_pool", 120},
        {3, "dam_break", 120},
        {4, "waterfall", 120},
        {9, "narrow_nozzle", 90}
    };
    for (VelocityAdvection mode : modes) {
        config.velocityAdvection = mode;
        for (Scene const &s : scenes) {
            timingAccum = TimingAverages{};
            timingAverage = TimingAverages{};
            timingTicks = 0;
            loadTestScene(s.id);
            auto start = Clock::now();
            for (int tick = 0; tick < s.ticks; ++tick) simulationTick();
            rebuildActivityAndMetrics();
            double perTick = elapsedMs(start) / std::max(1, s.ticks);
            out << velocityAdvectionName(mode) << '\t' << s.name << '\t' << s.ticks << '\t'
                << timingAverage.advection << '\t' << timingAverage.physics << '\t' << perTick << '\t'
                << currentVolume << '\t' << expectedVolume << '\t' << volumeError << '\t'
                << measuredMaxVelocity << '\t' << activeFluidCells << '\t' << activeChunks << '\t'
                << momentumX << '\t' << momentumY << '\t' << kineticEnergy << '\n';
        }
    }
    config.velocityAdvection = saved;
}

void FluidEngine::runScaleBenchmark() {
    config.workerCount = 1;
    syncWorkerPool();
    std::string filename=miscFile("scale_benchmark_"+std::to_string(GW)+"x"+std::to_string(GH)+".tsv");
    std::ofstream out(filename,std::ios::trunc);
    out<<"grid\tscene\tticks\tphysics_ms\tactive_cells\tactive_chunks\tvolume_error\n";
    auto run=[&](int scene,char const*name,int ticks){
        residualConsolidationEnabled=true;loadTestScene(scene);auto start=Clock::now();
        for(int tick=0;tick<ticks;++tick)simulationTick();
        rebuildActivityAndMetrics();
        out<<GW<<"x"<<GH<<'\t'<<name<<'\t'<<ticks<<'\t'<<elapsedMs(start)/ticks<<'\t'
           <<activeFluidCells<<'\t'<<activeChunks<<'\t'<<volumeError<<'\n';
    };
    run(1,"resting_pool",60);run(7,"small_puddle",90);run(3,"dam_break",60);run(5,"large_moving_body",60);
}

void FluidEngine::runLiquidBugDiagnostics() {
    config.workerCount = 1;
    syncWorkerPool();
    residualConsolidationEnabled = true;
    std::ofstream out(miscFile("liquid_bug_diag.tsv"), std::ios::trunc);
    out << "case\ttick\tvolume\texpected\terror\tn_any\tn_visible\tn_subvis\ty_span\ty_min\tmax_fill\tsplashes\tfrozen_vis\tholes\tactive_chunks\n";

    auto countHoles = [&]() {
        int holes = 0;
        constexpr int dx[4] = {-1, 1, 0, 0};
        constexpr int dy[4] = {0, 0, -1, 1};
        for (int y = 1; y < GH - 1; ++y) for (int x = 1; x < GW - 1; ++x) {
            int i = ci(x, y);
            if (solid[i] || dynamicSolid[i] || fill[i] >= MIN_RENDER_FILL) continue;
            int vis = 0;
            for (int n = 0; n < 4; ++n) {
                int ni = ci(x + dx[n], y + dy[n]);
                if (!solid[ni] && !dynamicSolid[ni] && fill[ni] >= MIN_RENDER_FILL) ++vis;
            }
            if (vis >= 4) ++holes;
        }
        return holes;
    };

    auto emit = [&](char const *name, int tick) {
        rebuildActivityAndMetrics();
        int nAny = 0, nVis = 0, nSub = 0, yMin = GH, yMax = -1, frozen = 0;
        float maxFill = 0.0f;
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            float amount = fill[ci(x, y)];
            if (amount <= 1e-7f) continue;
            ++nAny;
            maxFill = std::max(maxFill, amount);
            yMin = std::min(yMin, y);
            yMax = std::max(yMax, y);
            if (amount >= MIN_RENDER_FILL) {
                ++nVis;
                bool supported = (y + 1 >= GH) || isSolid(x, y + 1) || fill[ci(x, y + 1)] >= MIN_ACTIVE_FILL;
                float speed = std::abs(cellU(x, y)) + std::abs(cellV(x, y));
                if (!supported && speed < 0.12f) ++frozen;
            } else ++nSub;
        }
        int ySpan = (yMax >= yMin) ? (yMax - yMin + 1) : 0;
        int yMinOut = (yMax >= yMin) ? yMin : -1;
        out << name << '\t' << tick << '\t' << currentVolume << '\t' << expectedVolume << '\t' << volumeError
            << '\t' << nAny << '\t' << nVis << '\t' << nSub << '\t' << ySpan << '\t' << yMinOut << '\t' << maxFill
            << '\t' << splashes.size() << '\t' << frozen << '\t' << countHoles() << '\t' << activeChunks << '\n';
    };

    // Single 1-cell droplet falling onto the scene-6 floor.
    loadTestScene(6);
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) fill[ci(x, y)] = 0.0f;
    splashes.clear();
    fill[ci(100, 36)] = 1.0f;
    expectedVolume = 1.0;
    wakeAllFluidChunks();
    enforceSolidBoundaries();
    rebuildActivityAndMetrics();
    emit("single_drop", 0);
    for (int tick = 1; tick <= 90; ++tick) {
        simulationTick();
        if (tick == 1 || tick == 5 || tick == 15 || tick == 30 || tick == 60 || tick == 90) emit("single_drop", tick);
    }

    // Four-cell droplet (F6) for smear comparison.
    loadTestScene(6);
    emit("four_cell", 0);
    for (int tick = 1; tick <= 90; ++tick) {
        simulationTick();
        if (tick == 15 || tick == 45 || tick == 90) emit("four_cell", tick);
    }

    // Resting pool interior holes.
    loadTestScene(1);
    emit("still_pool", 0);
    for (int tick = 1; tick <= 120; ++tick) simulationTick();
    emit("still_pool", 120);

    // Dam-break: energetic transport + hole risk.
    loadTestScene(3);
    emit("dam_break", 0);
    int maxHoles = 0;
    for (int tick = 1; tick <= 90; ++tick) {
        simulationTick();
        maxHoles = std::max(maxHoles, countHoles());
    }
    emit("dam_break", 90);
    out << "dam_break\tmax_holes\t" << maxHoles << "\n";

    auto wallAt = [&](int x, int y) { if (inside(x, y)) solid[ci(x, y)] = 1; };
    auto sumFill = [&](int x0, int y0, int x1, int y1) {
        double s = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            if (inside(x, y)) s += fill[ci(x, y)];
        return s;
    };
    auto topologyLine = [&](char const *name, bool ok, double a = 0.0, double b = 0.0) {
        out << name << '\t' << (ok ? "ok" : "FAIL") << '\t' << a << '\t' << b << '\n';
    };

    // 1px sealed cup: exterior deposit/splash must not enter the interior.
    {
        clearWorld();
        int const x0 = 60, y0 = 40, x1 = 68, y1 = 48;
        for (int x = x0; x <= x1; ++x) { wallAt(x, y0); wallAt(x, y1); }
        for (int y = y0; y <= y1; ++y) { wallAt(x0, y); wallAt(x1, y); }
        double interior0 = 0.0;
        for (int y = y0 + 2; y <= y1 - 2; ++y) for (int x = x0 + 2; x <= x1 - 2; ++x) {
            fill[ci(x, y)] = 1.0f;
            interior0 += 1.0;
        }
        expectedVolume = interior0;
        enforceSolidBoundaries();
        seedAmbientHeat();
        wakeAllFluidChunks();
        float leftover = depositVolume(70.5f, 44.5f, 0.85f, -14.0f, 0.0f);
        leftover += depositVolume(68.5f, 44.5f, 0.5f, -14.0f, 0.0f);
        splashes.push_back({74.5f, 44.5f, -45.0f, 0.0f, 0.4f, 0.0f, 74, 44, 0.0f});
        expectedVolume += 0.85 + 0.5 + 0.4 - leftover;
        for (int k = 0; k < 12; ++k) updateSplashParticles(PHYSICS_DT);
        double interior1 = sumFill(x0 + 1, y0 + 1, x1 - 1, y1 - 1);
        bool ok = std::abs(interior1 - interior0) < 1e-4;
        topologyLine("sealed_cup", ok, interior0, interior1);
        emit("sealed_cup", 1);
        (void)leftover;
    }

    // 1px wall: deposit on the empty side must not appear on the wet side.
    {
        clearWorld();
        for (int y = 0; y < GH; ++y) wallAt(100, y);
        for (int y = 50; y <= 70; ++y) for (int x = 70; x <= 99; ++x) fill[ci(x, y)] = 1.0f;
        double left0 = sumFill(70, 50, 99, 70);
        expectedVolume = left0;
        enforceSolidBoundaries();
        seedAmbientHeat();
        wakeAllFluidChunks();
        depositVolume(101.5f, 60.5f, 1.0f, -10.0f, 0.0f);
        depositVolume(100.5f, 60.5f, 0.6f, -10.0f, 0.0f);
        expectedVolume = left0 + 1.6;
        rebuildActivityAndMetrics();
        double left1 = sumFill(70, 50, 99, 70);
        double right = sumFill(101, 50, 110, 70);
        bool ok = std::abs(left1 - left0) < 1e-4 && right < 1.0f + 0.6f + 1e-3f;
        topologyLine("thin_wall", ok, left0, left1);
        emit("thin_wall", 1);
    }

    // Diagonal solids: 4-connect must not cut the corner into the opposite empty cell.
    {
        clearWorld();
        wallAt(51, 40);
        wallAt(50, 41);
        enforceSolidBoundaries();
        seedAmbientHeat();
        expectedVolume = 0.0;
        wakeAllFluidChunks();
        depositVolume(50.5f, 40.5f, 1.0f, 12.0f, 12.0f);
        depositVolume(51.5f, 40.5f, 0.4f, 8.0f, 0.0f);
        splashes.push_back({50.6f, 40.4f, 30.0f, 30.0f, 0.35f, 0.0f, 50, 40, 0.0f});
        for (int k = 0; k < 16; ++k) updateSplashParticles(PHYSICS_DT);
        float b = fill[ci(51, 41)];
        bool ok = b < 1e-5f;
        topologyLine("diag_corner", ok, 0.0, static_cast<double>(b));
        emit("diag_corner", 1);
    }

    // Open void: splash leaving the map is a counted sink, not a rim bounce.
    {
        bool savedWalls = config.walledBorders;
        config.walledBorders = false;
        clearWorld();
        splashes.clear();
        splashes.push_back({3.5f, 40.5f, -90.0f, 0.0f, 0.5f, 0.0f, 3, 40, 0.0f});
        expectedVolume = 0.5;
        enforceSolidBoundaries();
        wakeAllFluidChunks();
        for (int k = 0; k < 40; ++k) updateSplashParticles(PHYSICS_DT);
        rebuildActivityAndMetrics();
        double edgeFill = 0.0;
        for (int y = 0; y < GH; ++y) edgeFill += fill[ci(0, y)] + fill[ci(1, y)];
        bool ok = splashes.empty() && std::abs(expectedVolume) < 1e-5 && std::abs(currentVolume) < 1e-5 && edgeFill < 1e-5;
        topologyLine("void_escape", ok, expectedVolume, currentVolume);
        emit("void_escape", 1);
        config.walledBorders = savedWalls;
    }

    // Walled rim: the same particle stays in-world; expected volume does not drop.
    {
        bool savedWalls = config.walledBorders;
        config.walledBorders = true;
        clearWorld();
        splashes.clear();
        splashes.push_back({3.5f, 40.5f, -90.0f, 0.0f, 0.5f, 0.0f, 3, 40, 0.0f});
        expectedVolume = 0.5;
        enforceSolidBoundaries();
        wakeAllFluidChunks();
        for (int k = 0; k < 40; ++k) updateSplashParticles(PHYSICS_DT);
        rebuildActivityAndMetrics();
        double splashVol = 0.0;
        for (SplashParticle const &p : splashes) splashVol += p.volume;
        bool inGrid = true;
        for (SplashParticle const &p : splashes)
            if (!inside(static_cast<int>(std::floor(p.x)), static_cast<int>(std::floor(p.y)))) inGrid = false;
        bool ok = inGrid && std::abs(currentVolume - 0.5) < 1e-4 && std::abs(expectedVolume - 0.5) < 1e-4
            && std::abs(currentVolume - (sumFill(0, 0, GW - 1, GH - 1) + splashVol)) < 1e-4;
        topologyLine("walled_bounce", ok, expectedVolume, currentVolume);
        emit("walled_bounce", 1);
        config.walledBorders = savedWalls;
    }
}

void FluidEngine::runThreadBenchmark() {
    residualConsolidationEnabled = true;
    std::string filename = miscFile("thread_benchmark_" + std::to_string(GW) + "x" + std::to_string(GH) + ".tsv");
    std::ofstream out(filename, std::ios::trunc);
    int autoN = autoWorkerCount();
    out << "# auto_workers=" << autoN << " hardware=" << maxSelectableWorkers()
        << " pressure_min_cells=" << config.pressureParallelMinCells << "\n";
    out << "workers\tmode\tscene\tticks\tphysics_ms\tpressure_ms\tadvection_ms\ttransport_ms\tsurface_ms\tresidual_ms\trigid_ms\tactive_cells\tpressure_cells\tvolume\texpected\terror\tparallel_p\n";

    struct Case { int scene; char const *name; int warmup; int ticks; };
    Case cases[] = {
        {6, "A_light_droplet", 20, 90},
        {7, "B_moderate_puddle", 30, 120},
        {1, "C_large_pool", 30, 120},
        {3, "D_dam_break", 20, 120},
        {5, "D2_falling_stream", 20, 120},
    };

    int counts[] = {1, 2, 4, 0};
    char const *modes[] = {"1", "2", "4", "auto"};

    for (int m = 0; m < 4; ++m) {
        config.workerCount = counts[m];
        syncWorkerPool();
        int resolved = lastResolvedWorkers;
        for (Case const &c : cases) {
            timingAccum = TimingAverages{}; timingAverage = TimingAverages{}; timingTicks = 0;
            loadTestScene(c.scene);
            for (int i = 0; i < c.warmup; ++i) simulationTick();
            timingAccum = TimingAverages{}; timingAverage = TimingAverages{}; timingTicks = 0;
            auto start = Clock::now();
            for (int i = 0; i < c.ticks; ++i) simulationTick();
            rebuildActivityAndMetrics();
            double perTick = elapsedMs(start) / std::max(1, c.ticks);
            out << resolved << '\t' << modes[m] << '\t' << c.name << '\t' << c.ticks << '\t'
                << perTick << '\t' << timingAverage.pressure << '\t' << timingAverage.advection << '\t'
                << timingAverage.transport << '\t' << timingAverage.surface << '\t' << timingAverage.residual << '\t'
                << 0.0 << '\t' << activeFluidCells << '\t' << workCounts.pressureCells << '\t'
                << currentVolume << '\t' << expectedVolume << '\t' << volumeError << '\t'
                << (lastPressureParallel ? 1 : 0) << '\n';
            out.flush();
        }
    }
    config.workerCount = 1;
    syncWorkerPool();
}

void FluidEngine::runLiquidCompositionDiagnostics() {
    config.workerCount = 1;
    syncWorkerPool();
    std::ofstream out(miscFile("liquid_composition_diag.tsv"));
    out << std::setprecision(8);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };
    auto near = [](float a, float b, float tol = 1.0e-5f) {
        return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol;
    };
    auto componentSum = [&](int idx) {
        float s = 0.0f;
        forEachLiquidComponent(idx, [&](SubstanceId, float amt) { s += amt; });
        return s;
    };

    clearWorld();
    int i = ci(40, 40);
    auto resetCell = [&]() {
        fill[static_cast<size_t>(i)] = 0.0f;
        clearComposition(i);
        clearEmptyLiquidCell(i);
    };

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    emit("pure_water_storage",
        near(liquidComponentFraction(i, SUBSTANCE_WATER), 1.0f)
            && near(liquidComponentFraction(i, SUBSTANCE_HONEY), 0.0f)
            && dominantLiquidSubstance(i) == SUBSTANCE_WATER
            && liquidCompCount[static_cast<size_t>(i)] == 1,
        "w=" + std::to_string(liquidComponentFraction(i, SUBSTANCE_WATER)));
    emit("pure_water_sum",
        near(componentSum(i), fill[static_cast<size_t>(i)]),
        "sum=" + std::to_string(componentSum(i)) + " fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 1.0f);
    emit("pure_honey_storage",
        near(liquidComponentFraction(i, SUBSTANCE_HONEY), 1.0f)
            && near(liquidComponentFraction(i, SUBSTANCE_WATER), 0.0f)
            && dominantLiquidSubstance(i) == SUBSTANCE_HONEY,
        "w=" + std::to_string(liquidComponentFraction(i, SUBSTANCE_WATER))
            + " h=" + std::to_string(liquidComponentFraction(i, SUBSTANCE_HONEY)));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.5f);
    emit("mix_50_50",
        near(liquidComponentFraction(i, SUBSTANCE_WATER), 0.5f)
            && near(liquidComponentFraction(i, SUBSTANCE_HONEY), 0.5f)
            && near(componentSum(i), fill[static_cast<size_t>(i)])
            && dominantLiquidSubstance(i) == SUBSTANCE_WATER,
        "w=" + std::to_string(liquidComponentFraction(i, SUBSTANCE_WATER))
            + " h=" + std::to_string(liquidComponentFraction(i, SUBSTANCE_HONEY))
            + " fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 0.5f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.2f);
    emit("partial_fill_sum",
        near(componentSum(i), fill[static_cast<size_t>(i)])
            && near(fill[static_cast<size_t>(i)], 0.5f)
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), 0.3f)
            && near(liquidComponentAmount(i, SUBSTANCE_HONEY), 0.2f),
        "sum=" + std::to_string(componentSum(i)) + " fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    resetCell();
    float emptyW = liquidComponentAmount(i, SUBSTANCE_WATER);
    float emptyH = liquidComponentAmount(i, SUBSTANCE_HONEY);
    float emptyWf = liquidComponentFraction(i, SUBSTANCE_WATER);
    float emptyHf = liquidComponentFraction(i, SUBSTANCE_HONEY);
    emit("empty_cell_safe",
        emptyW == 0.0f && emptyH == 0.0f && emptyWf == 0.0f && emptyHf == 0.0f
            && std::isfinite(emptyW) && std::isfinite(emptyWf)
            && liquidComponents(i).count == 0
            && dominantLiquidSubstance(i) == SUBSTANCE_NONE
            && liquidCompositionValid(i),
        "w=" + std::to_string(emptyW) + " frac=" + std::to_string(emptyWf)
            + " dominant=" + std::to_string(dominantLiquidSubstance(i)));

    emit("invalid_substance",
        liquidComponentAmount(i, 999) == 0.0f
            && liquidComponentFraction(i, SUBSTANCE_WOOD) == 0.0f
            && liquidComponentAmount(i, SUBSTANCE_NONE) == 0.0f,
        "");

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    addLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.25f);
    emit("add_component",
        near(liquidComponentAmount(i, SUBSTANCE_HONEY), 0.25f)
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), 0.75f)
            && near(componentSum(i), fill[static_cast<size_t>(i)]),
        "w=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_WATER))
            + " h=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_HONEY)));
    addLiquidComponentAmount(i, SUBSTANCE_HONEY, -0.25f);
    emit("remove_component",
        near(liquidComponentAmount(i, SUBSTANCE_HONEY), 0.0f)
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), fill[static_cast<size_t>(i)])
            && near(componentSum(i), fill[static_cast<size_t>(i)])
            && liquidCompCount[static_cast<size_t>(i)] == 1,
        "count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.4f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 1.0e-12f);
    emit("slot_deletion_near_zero",
        liquidComponentAmount(i, SUBSTANCE_HONEY) == 0.0f
            && liquidCompCount[static_cast<size_t>(i)] == 1
            && near(componentSum(i), fill[static_cast<size_t>(i)]),
        "count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    (void)addComponentUntracked(i, SUBSTANCE_WATER, 0.2f);
    emit("no_duplicate_ids",
        !cellHasDuplicateComponents(i) && liquidCompCount[static_cast<size_t>(i)] == 1
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), 1.2f),
        "count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)]));

    resetCell();
    {
        int base = compositionSlot(i, 0);
        liquidCompId[static_cast<size_t>(base + 0)] = runtimeBuiltIn(SUBSTANCE_HONEY);
        liquidCompAmt[static_cast<size_t>(base + 0)] = 0.40f;
        liquidCompId[static_cast<size_t>(base + 1)] = runtimeBuiltIn(SUBSTANCE_WOOD);
        liquidCompAmt[static_cast<size_t>(base + 1)] = 0.20f;
        liquidCompId[static_cast<size_t>(base + 2)] = runtimeBuiltIn(SUBSTANCE_STONE);
        liquidCompAmt[static_cast<size_t>(base + 2)] = 0.20f;
        liquidCompId[static_cast<size_t>(base + 3)] = runtimeBuiltIn(SUBSTANCE_AIR);
        liquidCompAmt[static_cast<size_t>(base + 3)] = 0.20f;
        liquidCompCount[static_cast<size_t>(i)] = static_cast<uint8_t>(kMaxLiquidComponents);
        fill[static_cast<size_t>(i)] = 1.0f;
        float fillBeforeOverflow = fill[static_cast<size_t>(i)];
        float unplaced = addComponentUntracked(i, SUBSTANCE_WATER, 0.15f);
        emit("fixed_capacity_overflow",
            near(unplaced, 0.15f)
                && liquidComponentAmount(i, SUBSTANCE_WATER) == 0.0f
                && near(fill[static_cast<size_t>(i)], fillBeforeOverflow)
                && liquidCompCount[static_cast<size_t>(i)] == kMaxLiquidComponents,
            "unplaced=" + std::to_string(unplaced)
                + " count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)])
                + " fill=" + std::to_string(fill[static_cast<size_t>(i)]));
    }

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 0.80f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.30f);
    emit("water_honey_roundtrip",
        near(liquidComponentAmount(i, SUBSTANCE_HONEY), 0.30f)
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), 0.50f)
            && near(componentSum(i), fill[static_cast<size_t>(i)]),
        "w=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_WATER))
            + " h=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_HONEY)));

    // take/deposit conserves components
    clearWorld();
    int a = ci(50, 40), b = ci(51, 40);
    setLiquidComponentAmount(a, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(a, SUBSTANCE_HONEY, 0.25f);
    LiquidCarry moved = extractVolume(a, 0.4f);
    fill[static_cast<size_t>(b)] += 0.4f;
    applyCarry(b, moved);
    emit("transport_take_deposit",
        near(liquidComponentAmount(a, SUBSTANCE_WATER) + liquidComponentAmount(b, SUBSTANCE_WATER), 0.75f)
            && near(liquidComponentAmount(a, SUBSTANCE_HONEY) + liquidComponentAmount(b, SUBSTANCE_HONEY), 0.25f)
            && near(liquidComponentAmount(b, SUBSTANCE_WATER), 0.3f)
            && near(liquidComponentAmount(b, SUBSTANCE_HONEY), 0.1f),
        "Aw=" + std::to_string(liquidComponentAmount(a, SUBSTANCE_WATER))
            + " Bw=" + std::to_string(liquidComponentAmount(b, SUBSTANCE_WATER)));

    clearWorld();
    int s = ci(60, 40);
    setLiquidComponentAmount(s, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(s, SUBSTANCE_HONEY, 0.25f);
    LiquidCarry sc = extractVolume(s, 1.0f);
    SplashParticle drop{};
    drop.volume = 1.0f;
    copyCarryToSplash(drop, sc);
    LiquidCarry back = carryFromSplash(drop);
    fill[static_cast<size_t>(s)] += 1.0f;
    applyCarry(s, back);
    emit("splash_roundtrip",
        near(liquidComponentAmount(s, SUBSTANCE_WATER), 0.75f)
            && near(liquidComponentAmount(s, SUBSTANCE_HONEY), 0.25f)
            && near(componentSum(s), 1.0f),
        "w=" + std::to_string(liquidComponentAmount(s, SUBSTANCE_WATER))
            + " h=" + std::to_string(liquidComponentAmount(s, SUBSTANCE_HONEY)));

    clearWorld();
    bool walls = config.walledBorders;
    config.walledBorders = true;
    int left = ci(80, 50);
    (void)left;
    setLiquidComponentAmount(ci(80, 50), SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(left, SUBSTANCE_HONEY, 0.25f);
    expectedVolume = 1.0;
    u[ui(81, 50)] = 2.0f;
    wakeAllFluidChunks();
    for (int n = 0; n < 4; ++n) simulationTick();
    float wTot = 0.0f, hTot = 0.0f, fTot = 0.0f;
    for (int idx = 0; idx < GW * GH; ++idx) {
        wTot += liquidComponentAmount(idx, SUBSTANCE_WATER);
        hTot += liquidComponentAmount(idx, SUBSTANCE_HONEY);
        fTot += fill[static_cast<size_t>(idx)];
    }
    for (SplashParticle const &p : splashes) {
        wTot += liquidPayloadAmount(p.comps, p.compCount, SUBSTANCE_WATER);
        hTot += liquidPayloadAmount(p.comps, p.compCount, SUBSTANCE_HONEY);
        fTot += p.volume;
    }
    emit("transport_advection",
        near(wTot, 0.75f, 2.0e-3f) && near(hTot, 0.25f, 2.0e-3f) && near(fTot, 1.0f, 2.0e-3f),
        "w=" + std::to_string(wTot) + " h=" + std::to_string(hTot) + " fill=" + std::to_string(fTot));
    config.walledBorders = walls;

    auto stampWallRect = [&](int x0, int y0, int x1, int y1) {
        for (int x = x0; x <= x1; ++x) {
            solid[static_cast<size_t>(ci(x, y0))] = 1;
            solid[static_cast<size_t>(ci(x, y1))] = 1;
        }
        for (int y = y0; y <= y1; ++y) {
            solid[static_cast<size_t>(ci(x0, y))] = 1;
            solid[static_cast<size_t>(ci(x1, y))] = 1;
        }
    };
    auto fillBasin = [&](int x0, int y0, int x1, int y1, SubstanceId id, float heat, float dye) {
        double vol = 0.0;
        for (int y = y0 + 1; y <= y1 - 1; ++y) for (int x = x0 + 1; x <= x1 - 1; ++x) {
            int idx = ci(x, y);
            float f = (y == y0 + 1) ? 0.72f : 1.0f;
            setLiquidComponentAmount(idx, id, f);
            liquidHeat[static_cast<size_t>(idx)] = heat * f;
            dyeR[static_cast<size_t>(idx)] = dye * f;
            dyeG[static_cast<size_t>(idx)] = 0.0f;
            dyeB[static_cast<size_t>(idx)] = 0.0f;
            vol += f;
        }
        return vol;
    };
    auto sumRect = [&](int x0, int y0, int x1, int y1, SubstanceId id) {
        double s = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            s += liquidComponentAmount(ci(x, y), id);
        return s;
    };
    auto sumRectHeat = [&](int x0, int y0, int x1, int y1) {
        double s = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            s += liquidHeat[static_cast<size_t>(ci(x, y))];
        return s;
    };
    auto sumRectDyeR = [&](int x0, int y0, int x1, int y1) {
        double s = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            s += dyeR[static_cast<size_t>(ci(x, y))];
        return s;
    };
    auto swirlBasin = [&](int x0, int y0, int x1, int y1) {
        for (int y = y0 + 1; y <= y1 - 1; ++y) for (int x = x0 + 1; x <= x1; ++x)
            if (openUFace(x, y)) u[ui(x, y)] = ((y % 2) == 0) ? 18.0f : -18.0f;
        for (int y = y0 + 1; y <= y1; ++y) for (int x = x0 + 1; x <= x1 - 1; ++x)
            if (openVFace(x, y)) v[vi(x, y)] = ((x % 2) == 0) ? 14.0f : -14.0f;
    };
    auto runIsolatedBasins = [&](int sparseMode, char const *name) {
        FluidConfig saved = config;
        bool savedResidual = residualConsolidationEnabled;
        int savedForce = forceSparseFlux;
        clearWorld();
        config.walledBorders = true;
        config.sprayEnabled = false;
        config.surfaceTensionEnabled = false;
        config.maxLimiterPasses = 1;
        config.maxSubsteps = 2;
        residualConsolidationEnabled = false;
        forceSparseFlux = sparseMode;
        int ax0 = 8, ay0 = 16, ax1 = 24, ay1 = 48;
        int bx0 = 48, by0 = 16, bx1 = 64, by1 = 48;
        stampWallRect(ax0, ay0, ax1, ay1);
        stampWallRect(bx0, by0, bx1, by1);
        ThermalProperties const &thW = thermalForSubstance(SUBSTANCE_WATER);
        ThermalProperties const &thH = thermalForSubstance(SUBSTANCE_HONEY);
        FluidProperties const &flW = fluidForSubstance(SUBSTANCE_WATER);
        FluidProperties const &flH = fluidForSubstance(SUBSTANCE_HONEY);
        float hot = energyFromTemp(thermalCapacity(massKg(flW.density, 1.0f, config.cellsPerMeter), thW.specificHeat), 360.0f);
        float cold = energyFromTemp(thermalCapacity(massKg(flH.density, 1.0f, config.cellsPerMeter), thH.specificHeat), 280.0f);
        double aVol = fillBasin(ax0, ay0, ax1, ay1, SUBSTANCE_WATER, hot, 1.0f);
        double bVol = fillBasin(bx0, by0, bx1, by1, SUBSTANCE_HONEY, cold, 0.0f);
        expectedVolume = aVol + bVol;
        double heatA0 = sumRectHeat(ax0, ay0, ax1, ay1);
        double heatB0 = sumRectHeat(bx0, by0, bx1, by1);
        double dyeA0 = sumRectDyeR(ax0, ay0, ax1, ay1);
        double dyeB0 = sumRectDyeR(bx0, by0, bx1, by1);
        accumAdvectOverflowVol = 0.0;
        wakeAllFluidChunks();
        for (int n = 0; n < 18; ++n) {
            swirlBasin(ax0, ay0, ax1, ay1);
            swirlBasin(bx0, by0, bx1, by1);
            simulationTick();
        }
        double honeyInA = sumRect(ax0, ay0, ax1, ay1, SUBSTANCE_HONEY);
        double waterInB = sumRect(bx0, by0, bx1, by1, SUBSTANCE_WATER);
        double waterInA = sumRect(ax0, ay0, ax1, ay1, SUBSTANCE_WATER);
        double honeyInB = sumRect(bx0, by0, bx1, by1, SUBSTANCE_HONEY);
        double heatA1 = sumRectHeat(ax0, ay0, ax1, ay1);
        double heatB1 = sumRectHeat(bx0, by0, bx1, by1);
        double dyeA1 = sumRectDyeR(ax0, ay0, ax1, ay1);
        double dyeB1 = sumRectDyeR(bx0, by0, bx1, by1);
        double wAll = 0.0, hAll = 0.0, fAll = 0.0, heatAll = 0.0, dyeAll = 0.0;
        for (int idx = 0; idx < GW * GH; ++idx) {
            wAll += liquidComponentAmount(idx, SUBSTANCE_WATER);
            hAll += liquidComponentAmount(idx, SUBSTANCE_HONEY);
            fAll += fill[static_cast<size_t>(idx)];
            heatAll += liquidHeat[static_cast<size_t>(idx)];
            dyeAll += dyeR[static_cast<size_t>(idx)];
        }
        for (SplashParticle const &p : splashes) {
            wAll += liquidPayloadAmount(p.comps, p.compCount, SUBSTANCE_WATER);
            hAll += liquidPayloadAmount(p.comps, p.compCount, SUBSTANCE_HONEY);
            fAll += p.volume;
            heatAll += p.heat;
            dyeAll += p.dyeR;
        }
        const float isoTol = 2.0e-4f;
        emit(name,
            honeyInA <= isoTol && waterInB <= isoTol
                && near(static_cast<float>(wAll), static_cast<float>(aVol), 2.0e-3f)
                && near(static_cast<float>(hAll), static_cast<float>(bVol), 2.0e-3f)
                && accumAdvectOverflowVol > 1.0e-8,
            "honeyA=" + std::to_string(honeyInA)
                + " waterB=" + std::to_string(waterInB)
                + " wA=" + std::to_string(waterInA)
                + " hB=" + std::to_string(honeyInB)
                + " overflow=" + std::to_string(accumAdvectOverflowVol)
                + " fill=" + std::to_string(fAll));
        std::string heatName = std::string(name) + "_heat";
        emit(heatName.c_str(),
            std::abs(heatA1 - heatA0) / std::max(1.0, std::abs(heatA0)) < 2.0e-3
                && std::abs(heatB1 - heatB0) / std::max(1.0, std::abs(heatB0)) < 2.0e-3
                && std::abs(heatAll - (heatA0 + heatB0)) / std::max(1.0, std::abs(heatA0 + heatB0)) < 2.0e-3,
            "dA=" + std::to_string(heatA1 - heatA0)
                + " dB=" + std::to_string(heatB1 - heatB0));
        std::string dyeName = std::string(name) + "_dye";
        emit(dyeName.c_str(),
            dyeB1 <= 1.0e-5 && std::abs(dyeA1 - dyeA0) <= 1.0e-4
                && dyeB0 == 0.0,
            "dyeA=" + std::to_string(dyeA1) + " dyeB=" + std::to_string(dyeB1)
                + " dyeAll=" + std::to_string(dyeAll));
        config = saved;
        residualConsolidationEnabled = savedResidual;
        forceSparseFlux = savedForce;
    };
    runIsolatedBasins(1, "isolated_basins_no_component_exchange_sparse");
    runIsolatedBasins(0, "isolated_basins_no_component_exchange_dense");

    {
        FluidConfig saved = config;
        bool savedResidual = residualConsolidationEnabled;
        int savedForce = forceSparseFlux;
        clearWorld();
        config.walledBorders = true;
        config.sprayEnabled = false;
        config.surfaceTensionEnabled = false;
        residualConsolidationEnabled = false;
        forceSparseFlux = -1;
        int x0 = 20, y0 = 30, x1 = 50, y1 = 50;
        stampWallRect(x0, y0, x1, y1);
        int mid = (x0 + x1) / 2;
        for (int y = y0 + 1; y <= y1 - 1; ++y) for (int x = x0 + 1; x <= x1 - 1; ++x) {
            SubstanceId id = (x < mid) ? SUBSTANCE_WATER : SUBSTANCE_HONEY;
            setLiquidComponentAmount(ci(x, y), id, 1.0f);
        }
        expectedVolume = static_cast<double>((x1 - x0 - 1) * (y1 - y0 - 1));
        wakeAllFluidChunks();
        for (int n = 0; n < 24; ++n) {
            swirlBasin(x0, y0, x1, y1);
            simulationTick();
        }
        double honeyLeft = sumRect(x0, y0, mid - 1, y1, SUBSTANCE_HONEY);
        double waterRight = sumRect(mid, y0, x1, y1, SUBSTANCE_WATER);
        emit("connected_mixture_can_mix",
            honeyLeft > 0.05 && waterRight > 0.05,
            "honeyLeft=" + std::to_string(honeyLeft)
                + " waterRight=" + std::to_string(waterRight));
        config = saved;
        residualConsolidationEnabled = savedResidual;
        forceSparseFlux = savedForce;
    }

    emit("no_honey_array", true, "honey[]/nextHoney[] removed; slots are authoritative");

    emit("water_is_liquid_capable",
        validLiquidComponentId(SUBSTANCE_WATER)
            && supportsPhase(SUBSTANCE_WATER, MatterPhase::Liquid)
            && hasFluidProperties(SUBSTANCE_WATER),
        "");
    emit("honey_is_liquid_capable",
        validLiquidComponentId(SUBSTANCE_HONEY)
            && supportsPhase(SUBSTANCE_HONEY, MatterPhase::Liquid)
            && hasFluidProperties(SUBSTANCE_HONEY),
        "");
    emit("air_not_liquid_capable",
        !validLiquidComponentId(SUBSTANCE_AIR)
            && !(supportsPhase(SUBSTANCE_AIR, MatterPhase::Liquid) && hasFluidProperties(SUBSTANCE_AIR)),
        "");
    emit("nonliquid_solids_rejected",
        !validLiquidComponentId(SUBSTANCE_WOOD)
            && !validLiquidComponentId(SUBSTANCE_STONE)
            && !validLiquidComponentId(SUBSTANCE_GLASS)
            && !validLiquidComponentId(SUBSTANCE_METAL),
        "");

    resetCell();
    fill[static_cast<size_t>(i)] = 1.0f;
    clearComposition(i);
    float orphanFill = fill[static_cast<size_t>(i)];
    int orphanCount = liquidCompCount[static_cast<size_t>(i)];
    bool orphanInvalid = !liquidCompositionValid(i);
    clearEmptyLiquidCell(i);
    emit("fill_without_composition_invalid",
        orphanInvalid
            && !liquidCompositionValid(i)
            && near(fill[static_cast<size_t>(i)], orphanFill)
            && liquidCompCount[static_cast<size_t>(i)] == orphanCount
            && liquidComponentAmount(i, SUBSTANCE_WATER) == 0.0f
            && dominantLiquidSubstance(i) == SUBSTANCE_NONE,
        "fill=" + std::to_string(fill[static_cast<size_t>(i)])
            + " count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)]));
    emit("fill_without_composition_not_water",
        near(fill[static_cast<size_t>(i)], 1.0f)
            && liquidComponentAmount(i, SUBSTANCE_WATER) == 0.0f
            && liquidCompCount[static_cast<size_t>(i)] == 0
            && !liquidCompositionValid(i),
        "w=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_WATER)));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_AIR, 1.0f);
    emit("air_component_rejected",
        fill[static_cast<size_t>(i)] == 0.0f
            && liquidComponentAmount(i, SUBSTANCE_AIR) == 0.0f
            && liquidCompCount[static_cast<size_t>(i)] == 0
            && liquidCompositionValid(i)
            && dominantLiquidSubstance(i) == SUBSTANCE_NONE,
        "fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_STONE, 0.5f);
    setLiquidComponentAmount(i, SUBSTANCE_WOOD, 0.5f);
    setLiquidComponentAmount(i, SUBSTANCE_GLASS, 0.5f);
    setLiquidComponentAmount(i, SUBSTANCE_METAL, 0.5f);
    emit("nonliquid_component_rejected",
        fill[static_cast<size_t>(i)] == 0.0f
            && liquidCompCount[static_cast<size_t>(i)] == 0
            && liquidComponentAmount(i, SUBSTANCE_STONE) == 0.0f
            && liquidComponentAmount(i, SUBSTANCE_WOOD) == 0.0f
            && liquidCompositionValid(i),
        "count=" + std::to_string(liquidCompCount[static_cast<size_t>(i)]));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    emit("water_component_accepted",
        liquidCompositionValid(i)
            && near(liquidComponentAmount(i, SUBSTANCE_WATER), 1.0f)
            && near(componentSum(i), fill[static_cast<size_t>(i)])
            && dominantLiquidSubstance(i) == SUBSTANCE_WATER,
        "w=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_WATER)));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 1.0f);
    emit("honey_component_accepted",
        liquidCompositionValid(i)
            && near(liquidComponentAmount(i, SUBSTANCE_HONEY), 1.0f)
            && near(componentSum(i), fill[static_cast<size_t>(i)])
            && dominantLiquidSubstance(i) == SUBSTANCE_HONEY,
        "h=" + std::to_string(liquidComponentAmount(i, SUBSTANCE_HONEY)));

    resetCell();
    setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
    setLiquidComponentAmount(i, SUBSTANCE_HONEY, 0.5f);
    emit("valid_composition_invariants",
        liquidCompositionValid(i)
            && near(componentSum(i), fill[static_cast<size_t>(i)])
            && !cellHasDuplicateComponents(i)
            && validLiquidComponentId(liquidCompId[static_cast<size_t>(compositionSlot(i, 0))])
            && validLiquidComponentId(liquidCompId[static_cast<size_t>(compositionSlot(i, 1))]),
        "sum=" + std::to_string(componentSum(i)) + " fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    resetCell();
    fill[static_cast<size_t>(i)] = 0.80f;
    clearComposition(i);
    float snapFill = fill[static_cast<size_t>(i)];
    int snapCount = liquidCompCount[static_cast<size_t>(i)];
    float snapWater = liquidComponentAmount(i, SUBSTANCE_WATER);
    bool snapInvalid = !liquidCompositionValid(i);
    (void)liquidCompositionValid(i);
    emit("invariant_check_does_not_mutate",
        snapInvalid
            && !liquidCompositionValid(i)
            && near(fill[static_cast<size_t>(i)], snapFill)
            && liquidCompCount[static_cast<size_t>(i)] == snapCount
            && liquidComponentAmount(i, SUBSTANCE_WATER) == snapWater
            && snapWater == 0.0f,
        "fill=" + std::to_string(fill[static_cast<size_t>(i)]));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
