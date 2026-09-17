#include "ThermalEngine.h"

#include "fluid/FluidEngine.h"
#include "fluid/DiagOutput.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/SubstanceRegistry.h"
#include "substance/PhaseTransfer.h"
#include "substance/LiquidMixtureProperties.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace {

bool chunkInside(int cx, int cy) {
    return cx >= 0 && cx < CHUNK_W && cy >= 0 && cy < CHUNK_H;
}

void addEnergy(float &energy, float capacity, float dQ) {
    energy += dQ;
    if (!std::isfinite(energy) || energy < 0.0f) energy = 0.0f;
    if (capacity > MIN_THERMAL_CAPACITY) {
        float t = tempFromEnergy(energy, capacity);
        energy = energyFromTemp(capacity, t);
    }
}

} // namespace

ThermalEngine::ThermalEngine()
    : chunkActivity(CHUNK_W * CHUNK_H, 0)
    , chunkQuietTicks(CHUNK_W * CHUNK_H, 0)
{
}

void ThermalEngine::clear() {
    std::fill(chunkActivity.begin(), chunkActivity.end(), 0);
    std::fill(chunkQuietTicks.begin(), chunkQuietTicks.end(), 0);
    work = ThermalWorkCounts{};
    thermalEnergyEscaped = 0.0;
    thermalEnergyEntered = 0.0;
    tickNo = 0;
}

void ThermalEngine::wakeCell(int x, int y) {
    if (x < 0 || x >= GW || y < 0 || y >= GH) return;
    int cx = x / CHUNK, cy = y / CHUNK;
    for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
        int nx = cx + ox, ny = cy + oy;
        if (!chunkInside(nx, ny)) continue;
        int c = ny * CHUNK_W + nx;
        chunkActivity[static_cast<size_t>(c)] = 1;
        chunkQuietTicks[static_cast<size_t>(c)] = 0;
    }
}

void ThermalEngine::wakeRect(int x0, int y0, int x1, int y1) {
    x0 = std::max(0, x0); y0 = std::max(0, y0);
    x1 = std::min(GW - 1, x1); y1 = std::min(GH - 1, y1);
    if (x1 < x0 || y1 < y0) return;
    for (int y = y0; y <= y1; y += CHUNK) wakeCell(x0, y);
    for (int y = y0; y <= y1; y += CHUNK) wakeCell(x1, y);
    for (int x = x0; x <= x1; x += CHUNK) { wakeCell(x, y0); wakeCell(x, y1); }
    wakeCell((x0 + x1) / 2, (y0 + y1) / 2);
}

void ThermalEngine::ingestEngineWakes(FluidEngine &fluid, GasEngine &gas) {
    auto ingest = [&](std::vector<uint8_t> &src) {
        if (src.size() != chunkActivity.size()) return;
        for (size_t c = 0; c < src.size(); ++c) {
            if (!src[c]) continue;
            int cx = static_cast<int>(c) % CHUNK_W, cy = static_cast<int>(c) / CHUNK_W;
            for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
                int nx = cx + ox, ny = cy + oy;
                if (!chunkInside(nx, ny)) continue;
                int nc = ny * CHUNK_W + nx;
                chunkActivity[static_cast<size_t>(nc)] = 1;
                chunkQuietTicks[static_cast<size_t>(nc)] = 0;
            }
            src[c] = 0;
        }
    };
    ingest(fluid.thermalChunkWake);
    ingest(gas.thermalChunkWake);
}

float ThermalEngine::liquidCapacity(FluidEngine const &fluid, int index) {
    if (index < 0 || index >= GW * GH) return 0.0f;
    float fill = fluid.fill[static_cast<size_t>(index)];
    if (fill <= 1.0e-8f) return 0.0f;
    // Mixture density + mass-weighted Cp; fill scales mass only. Invalid
    // composition uses reference solver numbers without writing Water identity.
    LiquidMixtureProperties mix = evaluateLiquidMixture(fluid.liquidComponents(index));
    float mass = massKg(mix.density, fill, fluid.config.cellsPerMeter);
    return thermalCapacity(mass, mix.specificHeat);
}

float ThermalEngine::wallCapacity(FluidEngine const &fluid, int index) {
    if (index < 0 || index >= GW * GH) return 0.0f;
    if (!fluid.solid[static_cast<size_t>(index)]) return 0.0f;
    float mass = massKg(mechanicalForSubstance(kStaticWallSubstance).densityRel, 1.0f, fluid.config.cellsPerMeter);
    return thermalCapacity(mass, thermalForSubstance(kStaticWallSubstance).specificHeat);
}

float ThermalEngine::gasCapacity(GasEngine const &gas, int index) {
    if (index < 0 || index >= GW * GH) return 0.0f;
    if (!(gas.amount[static_cast<size_t>(index)] > GAS_MIN_AMOUNT)) return 0.0f;
    GasComponentView view = gas.gasComponents(index);
    if (view.count <= 0) return 0.0f;
    if (view.count == 1 && view.items[0].id == SUBSTANCE_AIR)
        return thermalCapacity(gasMassKg(view.items[0].amount),
            thermalForSubstance(SUBSTANCE_AIR).specificHeat);
    float cap = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId id = view.items[n].id;
        float amt = view.items[n].amount;
        if (!(amt > GAS_MIN_AMOUNT) || !validGasComponentId(id)) continue;
        float mass = 0.0f;
        if (id == SUBSTANCE_AIR)
            mass = gasMassKg(amt);
        else {
            mass = static_cast<float>(gasAmountToMassKg(id, amt, 4.0));
            // Numerical safety only: missing molar/density data uses Air's
            // reference mass scale. Does not relabel the component as Air.
            if (!(mass > 0.0f) || !std::isfinite(mass))
                mass = gasMassKg(amt);
        }
        cap += thermalCapacity(mass, gasPhaseSpecificHeat(id));
    }
    return cap;
}

float ThermalEngine::rigidPixelCapacity(RigidBody const &b, int localIndex) {
    if (localIndex < 0 || localIndex >= static_cast<int>(b.mask.size())) return 0.0f;
    MaterialId id = b.mask[static_cast<size_t>(localIndex)];
    if (id == MATERIAL_EMPTY) return 0.0f;
    MaterialDefinition const &mat = materialDef(id);
    float remain = 1.0f;
    if (localIndex < static_cast<int>(b.solidRemain.size()))
        remain = std::max(0.0f, b.solidRemain[static_cast<size_t>(localIndex)]);
    float mass = massKg(mat.density, remain, 4.0f);
    return thermalCapacity(mass, solidPhaseSpecificHeat(substanceForMaterialId(id)));
}

bool ThermalEngine::isChunkActive(int x, int y) const {
    if (x < 0 || x >= GW || y < 0 || y >= GH) return false;
    int c = (y / CHUNK) * CHUNK_W + (x / CHUNK);
    if (c < 0 || c >= static_cast<int>(chunkActivity.size())) return false;
    return chunkActivity[static_cast<size_t>(c)] != 0;
}

float ThermalEngine::liquidTempK(FluidEngine const &fluid, int index) {
    return tempFromEnergy(fluid.liquidHeat[static_cast<size_t>(index)], liquidCapacity(fluid, index));
}
float ThermalEngine::wallTempK(FluidEngine const &fluid, int index) {
    return tempFromEnergy(fluid.solidHeat[static_cast<size_t>(index)], wallCapacity(fluid, index));
}
float ThermalEngine::gasTempK(GasEngine const &gas, int index) {
    return tempFromEnergy(gas.heat[static_cast<size_t>(index)], gasCapacity(gas, index));
}
float ThermalEngine::rigidPixelTempK(RigidBody const &b, int localIndex) {
    if (localIndex < 0 || localIndex >= static_cast<int>(b.heat.size())) return AMBIENT_TEMPERATURE_K;
    return tempFromEnergy(b.heat[static_cast<size_t>(localIndex)], rigidPixelCapacity(b, localIndex));
}

void ThermalEngine::seedAmbient(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas) {
    float cpm = fluid.config.cellsPerMeter;
    for (int i = 0; i < GW * GH; ++i) {
        float fill = fluid.fill[static_cast<size_t>(i)];
        if (fill > 1.0e-8f) {
            LiquidMixtureProperties mix = evaluateLiquidMixture(fluid.liquidComponents(i));
            float cap = thermalCapacity(massKg(mix.density, fill, cpm), mix.specificHeat);
            fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        } else {
            fluid.clearEmptyLiquidCell(i);
        }
        if (fluid.solid[static_cast<size_t>(i)]) {
            float cap = thermalCapacity(massKg(mechanicalForSubstance(kStaticWallSubstance).densityRel, 1.0f, cpm),
                thermalForSubstance(kStaticWallSubstance).specificHeat);
            fluid.solidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        } else {
            fluid.solidHeat[static_cast<size_t>(i)] = 0.0f;
        }
        float a = gas.amount[static_cast<size_t>(i)];
        if (a > GAS_MIN_AMOUNT) {
            float cap = gasCapacity(gas, i);
            gas.heat[static_cast<size_t>(i)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        } else {
            gas.heat[static_cast<size_t>(i)] = 0.0f;
        }
    }
    for (RigidBody &b : rigid.bodies) {
        if (b.heat.size() != b.mask.size()) b.heat.assign(b.mask.size(), 0.0f);
        for (int li : b.occupiedLocal) {
            float cap = rigidPixelCapacity(b, li);
            b.heat[static_cast<size_t>(li)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        }
    }
    clear();
}

ThermalCellSample ThermalEngine::sampleCell(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas, int x, int y)
{
    ThermalCellSample s;
    if (!FluidEngine::inside(x, y)) return s;
    int i = FluidEngine::ci(x, y);
    int body = rigid.occupant[static_cast<size_t>(i)];
    if (body >= 0 && body < static_cast<int>(rigid.bodies.size())) {
        RigidBody const &b = rigid.bodies[static_cast<size_t>(body)];
        float lx, ly;
        RigidBodyEngine::worldToLocal(b, x + 0.5f, y + 0.5f, lx, ly);
        int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
        if (RigidBodyEngine::maskOccupied(b, ix, iy)) {
            int li = iy * b.maskW + ix;
            s.kind = ThermalSampleKind::Rigid;
            s.hasMatter = true;
            s.materialId = b.mask[static_cast<size_t>(li)];
            s.capacityJK = rigidPixelCapacity(b, li);
            s.energyJ = (li >= 0 && li < static_cast<int>(b.heat.size())) ? b.heat[static_cast<size_t>(li)] : 0.0f;
            s.temperatureK = rigidPixelTempK(b, li);
            return s;
        }
    }
    if (fluid.solid[static_cast<size_t>(i)]) {
        s.kind = ThermalSampleKind::Wall;
        s.hasMatter = true;
        s.materialId = rigidMaterialForSubstance(kStaticWallSubstance);
        s.capacityJK = wallCapacity(fluid, i);
        s.energyJ = fluid.solidHeat[static_cast<size_t>(i)];
        s.temperatureK = wallTempK(fluid, i);
        return s;
    }
    float fill = fluid.fill[static_cast<size_t>(i)];
    float gasAmt = gas.amount[static_cast<size_t>(i)];
    float liqCap = (fill > 1.0e-8f) ? liquidCapacity(fluid, i) : 0.0f;
    float gCap = (gasAmt > GAS_MIN_AMOUNT) ? gasCapacity(gas, i) : 0.0f;
    if (fill >= MIN_ACTIVE_FILL || (liqCap > gCap && liqCap > MIN_THERMAL_CAPACITY)) {
        s.kind = ThermalSampleKind::Liquid;
        s.hasMatter = true;
        s.capacityJK = liqCap;
        s.energyJ = fluid.liquidHeat[static_cast<size_t>(i)];
        s.temperatureK = liquidTempK(fluid, i);
        return s;
    }
    if (gCap > MIN_THERMAL_CAPACITY) {
        s.kind = ThermalSampleKind::Gas;
        s.hasMatter = true;
        s.capacityJK = gCap;
        s.energyJ = gas.heat[static_cast<size_t>(i)];
        s.temperatureK = gasTempK(gas, i);
        return s;
    }
    return s;
}

float ThermalEngine::sampleTemperatureK(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas, int x, int y)
{
    return sampleCell(fluid, rigid, gas, x, y).temperatureK;
}

ThermalWorldStats ThermalEngine::collectStats(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas) const
{
    ThermalWorldStats st;
    auto consider = [](float t, float &tMin, float &tMax, int &n) {
        if (n == 0) { tMin = tMax = t; }
        else { tMin = std::min(tMin, t); tMax = std::max(tMax, t); }
        ++n;
    };
    double tSum = 0.0;
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        ThermalCellSample s = sampleCell(fluid, rigid, gas, x, y);
        float t = s.temperatureK;
        if (!std::isfinite(t)) {
            if (std::isnan(t)) ++st.nNan;
            else ++st.nInf;
            continue;
        }
        if (t < 0.0f) ++st.nNegK;
        if (t <= 0.05f) ++st.nZeroK;
        consider(t, st.tMin, st.tMax, st.nSampled);
        tSum += t;
        st.maxAbsDeltaAmbient = std::max(st.maxAbsDeltaAmbient, std::abs(t - AMBIENT_TEMPERATURE_K));
        if (std::abs(t - AMBIENT_TEMPERATURE_K) > TEMP_VIZ_DEADBAND_K) ++st.nAwayFromAmbient;
        int i = FluidEngine::ci(x, y);
        switch (s.kind) {
            case ThermalSampleKind::Gas: consider(t, st.tMinGas, st.tMaxGas, st.nGas); break;
            case ThermalSampleKind::Liquid:
                consider(t, st.tMinLiquid, st.tMaxLiquid, st.nLiquid);
                if (st.nLiquid == 1 || t <= st.tMinLiquid) st.tMinLiquidFill = fluid.fill[static_cast<size_t>(i)];
                if (st.nLiquid == 1 || t >= st.tMaxLiquid) st.tMaxLiquidFill = fluid.fill[static_cast<size_t>(i)];
                break;
            case ThermalSampleKind::Wall: consider(t, st.tMinSolid, st.tMaxSolid, st.nWall); break;
            case ThermalSampleKind::Rigid: consider(t, st.tMinRigid, st.tMaxRigid, st.nRigid); break;
            default: ++st.nEmpty; break;
        }
        if (fluid.fill[static_cast<size_t>(i)] <= 1.0e-8f && fluid.liquidHeat[static_cast<size_t>(i)] > 0.0f) {
            ++st.leftoverLiquidCells;
            st.leftoverLiquidEnergy += fluid.liquidHeat[static_cast<size_t>(i)];
        }
        if (!fluid.solid[static_cast<size_t>(i)] && fluid.solidHeat[static_cast<size_t>(i)] > 0.0f) {
            ++st.leftoverSolidCells;
            st.leftoverSolidEnergy += fluid.solidHeat[static_cast<size_t>(i)];
        }
    }
    st.energy = totalThermalEnergy(fluid, rigid, gas);
    st.tMean = (st.nSampled > 0) ? static_cast<float>(tSum / st.nSampled) : AMBIENT_TEMPERATURE_K;
    return st;
}

double ThermalEngine::totalThermalEnergy(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas) const
{
    double e = 0.0;
    for (int i = 0; i < GW * GH; ++i) {
        e += fluid.liquidHeat[static_cast<size_t>(i)];
        e += fluid.solidHeat[static_cast<size_t>(i)];
        e += gas.heat[static_cast<size_t>(i)];
    }
    for (RigidBody const &b : rigid.bodies)
        for (float h : b.heat) e += h;
    for (SplashParticle const &p : fluid.splashes) e += p.heat;
    return e;
}

double ThermalEngine::netExternalHeat(FluidEngine const &fluid, GasEngine const &gas) const
{
    return thermalEnergyEscaped - thermalEnergyEntered + gas.escapedHeat + fluid.escapedHeat;
}

void ThermalEngine::applyBrush(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    int cx, int cy, int brushRadius, float signedStrength, float dt, BrushShape shape)
{
    float watts = config.heatToolWatts * signedStrength;
    float dQcell = watts * std::max(dt, 1.0f / 30.0f);
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if (!brushContains(shape, cx, cy, x, y, brushRadius)) continue;
            int i = FluidEngine::ci(x, y);
            bool any = false;
            int body = rigid.occupant[static_cast<size_t>(i)];
            if (body >= 0 && body < static_cast<int>(rigid.bodies.size())) {
                RigidBody &b = rigid.bodies[static_cast<size_t>(body)];
                float lx, ly;
                RigidBodyEngine::worldToLocal(b, x + 0.5f, y + 0.5f, lx, ly);
                int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
                if (RigidBodyEngine::maskOccupied(b, ix, iy)) {
                    int li = iy * b.maskW + ix;
                    if (li >= 0 && li < static_cast<int>(b.heat.size())) {
                        addEnergy(b.heat[static_cast<size_t>(li)], rigidPixelCapacity(b, li), dQcell);
                        any = true;
                    }
                }
            }
            if (!any && fluid.solid[static_cast<size_t>(i)]) {
                addEnergy(fluid.solidHeat[static_cast<size_t>(i)], wallCapacity(fluid, i), dQcell);
                any = true;
            }
            if (!any && fluid.fill[static_cast<size_t>(i)] > 1.0e-6f) {
                addEnergy(fluid.liquidHeat[static_cast<size_t>(i)], liquidCapacity(fluid, i), dQcell);
                any = true;
            }
            if (!any && gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT
                && !fluid.solid[static_cast<size_t>(i)] && !fluid.dynamicSolid[static_cast<size_t>(i)]) {
                addEnergy(gas.heat[static_cast<size_t>(i)], gasCapacity(gas, i), dQcell);
                any = true;
            }
            if (any) wakeCell(x, y);
            if (any) gas.wakeAt(x, y);
        }
}

void ThermalEngine::wakeMovingBodies(RigidBodyEngine const &rigid) {
    for (RigidBody const &b : rigid.bodies) {
        bool moved = std::abs(b.aabbX0 - b.prevAabbX0) > 0.05f
            || std::abs(b.aabbY0 - b.prevAabbY0) > 0.05f
            || std::abs(b.aabbX1 - b.prevAabbX1) > 0.05f
            || std::abs(b.aabbY1 - b.prevAabbY1) > 0.05f
            || std::abs(b.vx) + std::abs(b.vy) + std::abs(b.omega) > 0.05f;
        if (!moved) continue;
        bool interesting = false;
        for (int li : b.occupiedLocal) {
            float t = rigidPixelTempK(b, li);
            if (std::abs(t - AMBIENT_TEMPERATURE_K) > config.sleepAmbientEps) { interesting = true; break; }
        }
        if (!interesting) continue;
        int x0 = std::max(0, static_cast<int>(std::floor(b.aabbX0)));
        int y0 = std::max(0, static_cast<int>(std::floor(b.aabbY0)));
        int x1 = std::min(GW - 1, static_cast<int>(std::ceil(b.aabbX1)));
        int y1 = std::min(GH - 1, static_cast<int>(std::ceil(b.aabbY1)));
        wakeRect(x0, y0, x1, y1);
    }
}

void ThermalEngine::conductRigidBodies(RigidBodyEngine &rigid, float dt, float areaOverDx) {
    work.activeBodies = 0;
    for (RigidBody &b : rigid.bodies) {
        if (b.heat.size() != b.mask.size()) b.heat.assign(b.mask.size(), 0.0f);
        if (b.occupiedLocal.size() < 2) continue;
        bool any = false;
        float tMin = 1.0e9f, tMax = -1.0e9f;
        for (int li : b.occupiedLocal) {
            float t = rigidPixelTempK(b, li);
            tMin = std::min(tMin, t);
            tMax = std::max(tMax, t);
        }
        if (tMax - tMin < config.sleepTempEps &&
            std::abs(tMax - AMBIENT_TEMPERATURE_K) < config.sleepAmbientEps &&
            std::abs(tMin - AMBIENT_TEMPERATURE_K) < config.sleepAmbientEps)
            continue;
        ++work.activeBodies;
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW, ly = li / b.maskW;
            auto neigh = [&](int nx, int ny) {
                if (nx < 0 || ny < 0 || nx >= b.maskW || ny >= b.maskH) return;
                if (!RigidBodyEngine::maskOccupied(b, nx, ny)) return;
                int ni = ny * b.maskW + nx;
                if (ni <= li) return;
                MaterialId a = b.mask[static_cast<size_t>(li)];
                MaterialId c = b.mask[static_cast<size_t>(ni)];
                float Ca = rigidPixelCapacity(b, li);
                float Cb = rigidPixelCapacity(b, ni);
                float dQ = exchangeThermalEnergy(
                    b.heat[static_cast<size_t>(li)], Ca, solidPhaseConductivity(substanceForMaterialId(a)),
                    b.heat[static_cast<size_t>(ni)], Cb, solidPhaseConductivity(substanceForMaterialId(c)),
                    dt, areaOverDx, config.conductivityScale);
                if (std::abs(dQ) > 1.0e-8f) {
                    ++work.conductionPairs;
                    any = true;
                }
            };
            neigh(lx + 1, ly);
            neigh(lx, ly + 1);
        }
        if (any) {
            int x0 = std::max(0, static_cast<int>(std::floor(b.aabbX0)));
            int y0 = std::max(0, static_cast<int>(std::floor(b.aabbY0)));
            int x1 = std::min(GW - 1, static_cast<int>(std::ceil(b.aabbX1)));
            int y1 = std::min(GH - 1, static_cast<int>(std::ceil(b.aabbY1)));
            wakeRect(x0, y0, x1, y1);
        }
    }
}

void ThermalEngine::conductActive(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas, float dt) {
    float dx = cellLengthM(fluid.config.cellsPerMeter);
    float areaOverDx = dx;
    work.conductionPairs = 0;
    work.activeCells = 0;

    auto exchangePair = [&](int x0, int y0, int x1, int y1) {
        float *eA = nullptr, *eB = nullptr, cA = 0, cB = 0, kA = 0, kB = 0;
        int bA = -1, bB = -1;
        bool gasA = false, gasB = false;
        if (!resolveNode(fluid, rigid, gas, x0, y0, eA, cA, kA, bA, gasA)) return;
        if (!resolveNode(fluid, rigid, gas, x1, y1, eB, cB, kB, bB, gasB)) return;
        if (bA >= 0 && bA == bB) return;
        if (eA == eB) return;
        float scale = (gasA || gasB) ? config.gasConductivityScale : config.conductivityScale;
        float dQ = exchangeThermalEnergy(*eA, cA, kA, *eB, cB, kB, dt, areaOverDx, scale);
        if (std::abs(dQ) > 1.0e-8f) {
            ++work.conductionPairs;
            if (std::abs(dQ) > 1.0e-3f) {
                wakeCell(x0, y0);
                wakeCell(x1, y1);
                if (gasA || gasB) {
                    gas.wakeAt(x0, y0);
                    gas.wakeAt(x1, y1);
                }
            }
        }
    };

    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
        if (!chunkActivity[static_cast<size_t>(c)]) continue;
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        int x0 = cx * CHUNK, y0 = cy * CHUNK;
        int x1 = std::min(GW - 1, x0 + CHUNK - 1);
        int y1 = std::min(GH - 1, y0 + CHUNK - 1);
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            ++work.activeCells;
            if (x + 1 < GW) exchangePair(x, y, x + 1, y);
            if (y + 1 < GH) exchangePair(x, y, x, y + 1);
        }
    }
    conductRigidBodies(rigid, dt, areaOverDx);
}

bool ThermalEngine::resolveNode(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    int x, int y, float *&energy, float &cap, float &k, int &bodyId, bool &isGas)
{
    energy = nullptr; cap = 0.0f; k = 0.0f; bodyId = -1; isGas = false;
    if (!FluidEngine::inside(x, y)) return false;
    int i = FluidEngine::ci(x, y);
    int body = rigid.occupant[static_cast<size_t>(i)];
    if (body >= 0 && body < static_cast<int>(rigid.bodies.size())) {
        RigidBody &b = rigid.bodies[static_cast<size_t>(body)];
        float lx, ly;
        RigidBodyEngine::worldToLocal(b, x + 0.5f, y + 0.5f, lx, ly);
        int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
        if (RigidBodyEngine::maskOccupied(b, ix, iy)) {
            int li = iy * b.maskW + ix;
            if (li >= 0 && li < static_cast<int>(b.heat.size())) {
                energy = &b.heat[static_cast<size_t>(li)];
                cap = rigidPixelCapacity(b, li);
                k = solidPhaseConductivity(substanceForMaterialId(b.mask[static_cast<size_t>(li)]));
                bodyId = body;
                return cap > MIN_THERMAL_CAPACITY;
            }
        }
    }
    if (fluid.solid[static_cast<size_t>(i)]) {
        energy = &fluid.solidHeat[static_cast<size_t>(i)];
        cap = wallCapacity(fluid, i);
        k = thermalForSubstance(kStaticWallSubstance).conductivity;
        return cap > MIN_THERMAL_CAPACITY;
    }
    if (fluid.fill[static_cast<size_t>(i)] > 1.0e-6f) {
        energy = &fluid.liquidHeat[static_cast<size_t>(i)];
        cap = liquidCapacity(fluid, i);
        k = fluid.mixConductivity(i);
        return cap > MIN_THERMAL_CAPACITY;
    }
    if (gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT) {
        energy = &gas.heat[static_cast<size_t>(i)];
        cap = gasCapacity(gas, i);
        k = thermalForSubstance(substanceForGasSpecies()).conductivity;
        isGas = true;
        return cap > MIN_THERMAL_CAPACITY;
    }
    return false;
}

void ThermalEngine::conductOpenBoundary(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas, float dt)
{
    if (fluid.config.walledBorders) return;
    float dx = cellLengthM(fluid.config.cellsPerMeter);
    float areaOverDx = dx;
    float kAir = thermalForSubstance(SUBSTANCE_AIR).conductivity;
    auto expose = [&](int x, int y, int faces) {
        if (faces <= 0) return;
        float *e = nullptr, cap = 0.0f, k = 0.0f;
        int bodyId = -1;
        bool isGas = false;
        if (!resolveNode(fluid, rigid, gas, x, y, e, cap, k, bodyId, isGas)) return;
        float dQ = exchangeThermalEnergyWithAmbient(*e, cap, k, kAir, dt,
            areaOverDx * static_cast<float>(faces), config.gasConductivityScale);
        if (!std::isfinite(dQ) || dQ == 0.0f) return;
        if (dQ > 0.0f) thermalEnergyEscaped += dQ;
        else thermalEnergyEntered += static_cast<double>(-dQ);
        ++work.conductionPairs;
        if (std::abs(dQ) > 1.0e-3f) {
            wakeCell(x, y);
            if (isGas) gas.wakeAt(x, y);
        }
    };
    for (int y = 0; y < GH; ++y) {
        int facesL = 1 + (y == 0 || y == GH - 1 ? 1 : 0);
        int facesR = 1 + (y == 0 || y == GH - 1 ? 1 : 0);
        expose(0, y, facesL);
        expose(GW - 1, y, facesR);
    }
    for (int x = 1; x < GW - 1; ++x) {
        expose(x, 0, 1);
        expose(x, GH - 1, 1);
    }
}

void ThermalEngine::sleepChunks(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas) {
    work.activeChunks = 0;
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
        if (!chunkActivity[static_cast<size_t>(c)]) continue;
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        int x0 = cx * CHUNK, y0 = cy * CHUNK;
        int x1 = std::min(GW - 1, x0 + CHUNK - 1);
        int y1 = std::min(GH - 1, y0 + CHUNK - 1);
        float maxD = 0.0f;
        float maxAmb = 0.0f;
        auto consider = [&](float t) {
            maxAmb = std::max(maxAmb, std::abs(t - AMBIENT_TEMPERATURE_K));
        };
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            float t = sampleTemperatureK(fluid, rigid, gas, x, y);
            consider(t);
            if (x + 1 <= x1) maxD = std::max(maxD, std::abs(t - sampleTemperatureK(fluid, rigid, gas, x + 1, y)));
            if (y + 1 <= y1) maxD = std::max(maxD, std::abs(t - sampleTemperatureK(fluid, rigid, gas, x, y + 1)));
            auto halo = [&](int nx, int ny) {
                if (!FluidEngine::inside(nx, ny)) return;
                float tn = sampleTemperatureK(fluid, rigid, gas, nx, ny);
                float d = std::abs(t - tn);
                maxD = std::max(maxD, d);
                if (d >= config.sleepTempEps) wakeCell(nx, ny);
            };
            if (x == x0) halo(x - 1, y);
            if (x == x1) halo(x + 1, y);
            if (y == y0) halo(x, y - 1);
            if (y == y1) halo(x, y + 1);
        }
        if (maxD < config.sleepTempEps && maxAmb < config.sleepAmbientEps) {
            uint8_t q = chunkQuietTicks[static_cast<size_t>(c)];
            if (q < 250) ++q;
            chunkQuietTicks[static_cast<size_t>(c)] = q;
            if (q >= config.sleepQuietTicks) chunkActivity[static_cast<size_t>(c)] = 0;
        } else {
            chunkQuietTicks[static_cast<size_t>(c)] = 0;
        }
        if (chunkActivity[static_cast<size_t>(c)]) ++work.activeChunks;
    }
}

void ThermalEngine::publishTiming(FluidEngine &fluid) {
    fluid.timingAverage.thermal = timingAverage;
    fluid.workCounts.thermalCells = work.activeCells;
    fluid.workCounts.thermalChunks = work.activeChunks;
    fluid.workCounts.thermalPairs = work.conductionPairs;
}

void ThermalEngine::scrubMasslessHeat(FluidEngine &fluid) {
    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
        if (!chunkActivity[static_cast<size_t>(c)]) continue;
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        int x0 = cx * CHUNK, y0 = cy * CHUNK;
        int x1 = std::min(GW - 1, x0 + CHUNK - 1);
        int y1 = std::min(GH - 1, y0 + CHUNK - 1);
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            int i = FluidEngine::ci(x, y);
            if (fluid.fill[static_cast<size_t>(i)] <= 1.0e-8f)
                fluid.clearEmptyLiquidCell(i);
            if (!fluid.solid[static_cast<size_t>(i)])
                fluid.solidHeat[static_cast<size_t>(i)] = 0.0f;
        }
    }
}

void ThermalEngine::simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas, float dt) {
    auto start = std::chrono::steady_clock::now();
    ++tickNo;
    ingestEngineWakes(fluid, gas);
    wakeMovingBodies(rigid);

    if (!config.enabled) {
        lastStepMs = 0.0;
        work = ThermalWorkCounts{};
        publishTiming(fluid);
        return;
    }

    for (int c = 0; c < CHUNK_W * CHUNK_H; ++c) {
        if (!chunkActivity[static_cast<size_t>(c)]) continue;
        int cx = c % CHUNK_W, cy = c / CHUNK_W;
        int x0 = cx * CHUNK, y0 = cy * CHUNK;
        int x1 = std::min(GW - 1, x0 + CHUNK - 1);
        int y1 = std::min(GH - 1, y0 + CHUNK - 1);
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            int i = FluidEngine::ci(x, y);
            if (fluid.solid[static_cast<size_t>(i)] || fluid.fill[static_cast<size_t>(i)] > 1.0e-6f) continue;
            if (gas.amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT) continue;
            float t = gasTempK(gas, i);
            if (std::abs(t - AMBIENT_TEMPERATURE_K) > config.sleepAmbientEps)
                gas.wakeAt(x, y);
        }
    }

    int interval = std::max(1, config.intervalTicks);
    bool runConduction = (tickNo % static_cast<uint32_t>(interval)) == 0u;
    if (runConduction) {
        bool any = false;
        for (uint8_t a : chunkActivity) if (a) { any = true; break; }
        float dtCond = dt * static_cast<float>(interval);
        if (any) {
            conductActive(fluid, rigid, gas, dtCond);
            scrubMasslessHeat(fluid);
        } else {
            work.activeCells = 0;
            work.conductionPairs = 0;
            work.activeBodies = 0;
        }
        if (!fluid.config.walledBorders)
            conductOpenBoundary(fluid, rigid, gas, dtCond);
        sleepChunks(fluid, rigid, gas);
    }

    lastStepMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    timingAccum += lastStepMs;
    ++timingTicks;
    if (timingTicks >= 60) {
        timingAverage = timingAccum / static_cast<double>(timingTicks);
        timingAccum = 0.0;
        timingTicks = 0;
    }
    fluid.timingAccum.thermal += lastStepMs;
    publishTiming(fluid);
}

void ThermalEngine::runDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas) {
    std::ofstream out(miscFile("thermal_diag.tsv"), std::ios::trunc);
    out << "case\ttick\tenergy\tactive_chunks\tactive_cells\tpairs\tthermal_ms\textra\tpass\n";

    auto worldTick = [&]() {
        rigid.step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gas.simulationTick(fluid);
        rigid.gatherFluidForces(fluid);
        gas.applyPressureForces(rigid, fluid);
        simulationTick(fluid, rigid, gas, PHYSICS_DT);
    };
    auto tickN = [&](int n) { for (int i = 0; i < n; ++i) worldTick(); };

    auto emit = [&](char const *name, int tick, double energy, double extra, bool pass) {
        out << name << '\t' << tick << '\t' << energy << '\t' << work.activeChunks << '\t'
            << work.activeCells << '\t' << work.conductionPairs << '\t' << lastStepMs << '\t'
            << extra << '\t' << (pass ? 1 : 0) << '\n';
    };

    auto resetBare = [&]() {
        rigid.clear();
        fluid.clearWorld();
        fluid.config.walledBorders = true;
        gas.config.boundary = GasBoundary::Sealed;
        gas.config.simMode = GasSimMode::Off;
        gas.resetAmbient(fluid);
        std::fill(gas.amount.begin(), gas.amount.end(), 0.0f);
        std::fill(gas.heat.begin(), gas.heat.end(), 0.0f);
        config.enabled = true;
        config.intervalTicks = 1;
        config.sleepTempEps = 0.05f;
        config.sleepAmbientEps = 0.15f;
        seedAmbient(fluid, rigid, gas);
        std::fill(gas.amount.begin(), gas.amount.end(), 0.0f);
        std::fill(gas.heat.begin(), gas.heat.end(), 0.0f);
    };

    auto paintBodyRect = [&](MaterialId mat, int x0, int y0, int x1, int y1, bool anchored) {
        rigid.drawMaterial = mat;
        rigid.placeAnchored = anchored;
        rigid.placeSleeping = false;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            if (FluidEngine::inside(x, y))
                rigid.pending[static_cast<size_t>(FluidEngine::ci(x, y))] = mat;
        }
        rigid.commitPending(fluid);
        rigid.syncOccupancy(fluid);
    };

    auto setLiquid = [&](int x, int y, float fill, float tempK) {
        if (!FluidEngine::inside(x, y) || fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y))]) return;
        int i = FluidEngine::ci(x, y);
        fluid.fill[static_cast<size_t>(i)] = fill;
        float cap = liquidCapacity(fluid, i);
        fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, tempK);
        fluid.expectedVolume += fill;
        wakeCell(x, y);
    };

    // 1. Hot/cold liquid mixing
    resetBare();
    for (int x = 40; x <= 160; ++x) fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 90))] = 1;
    for (int y = 40; y <= 90; ++y) {
        fluid.solid[static_cast<size_t>(FluidEngine::ci(40, y))] = 1;
        fluid.solid[static_cast<size_t>(FluidEngine::ci(160, y))] = 1;
    }
    for (int x = 40; x <= 160; ++x) fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 40))] = 1;
    seedAmbient(fluid, rigid, gas);
    int nHot = 0, nCold = 0;
    for (int y = 50; y <= 80; ++y) for (int x = 50; x <= 95; ++x) { setLiquid(x, y, 1.0f, 350.0f); ++nHot; }
    for (int y = 50; y <= 80; ++y) for (int x = 105; x <= 150; ++x) { setLiquid(x, y, 1.0f, 280.0f); ++nCold; }
    fluid.rebuildActivityAndMetrics();
    double eMix0 = 0.0, cMix0 = 0.0;
    for (int i = 0; i < GW * GH; ++i) {
        eMix0 += fluid.liquidHeat[static_cast<size_t>(i)];
        cMix0 += liquidCapacity(fluid, i);
    }
    tickN(90);
    double eMix1 = 0.0, cMix1 = 0.0;
    for (int i = 0; i < GW * GH; ++i) {
        eMix1 += fluid.liquidHeat[static_cast<size_t>(i)];
        cMix1 += liquidCapacity(fluid, i);
    }
    float tMean = (cMix1 > MIN_THERMAL_CAPACITY) ? static_cast<float>(eMix1 / cMix1) : 0.0f;
    float tExpect = (350.0f * static_cast<float>(nHot) + 280.0f * static_cast<float>(nCold))
        / static_cast<float>(std::max(1, nHot + nCold));
    bool passMix = std::abs(tMean - tExpect) < 12.0f
        && std::abs(eMix1 - eMix0) / std::max(1.0, std::abs(eMix0)) < 0.04;
    emit("mix_hot_cold", 90, eMix1, tMean, passMix);

    // 2. Metal vs wood conduction
    resetBare();
    paintBodyRect(MATERIAL_METAL, 30, 50, 41, 53, true);
    paintBodyRect(MATERIAL_WOOD, 90, 50, 101, 53, true);
    seedAmbient(fluid, rigid, gas);
    auto heatEnd = [&](int bodyIndex, int localX0) {
        if (bodyIndex < 0 || bodyIndex >= static_cast<int>(rigid.bodies.size())) return;
        RigidBody &b = rigid.bodies[static_cast<size_t>(bodyIndex)];
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW;
            if (lx <= localX0) {
                float cap = rigidPixelCapacity(b, li);
                b.heat[static_cast<size_t>(li)] = energyFromTemp(cap, 450.0f);
            }
        }
        wakeRect(static_cast<int>(b.aabbX0), static_cast<int>(b.aabbY0),
            static_cast<int>(b.aabbX1), static_cast<int>(b.aabbY1));
    };
    int metalI = -1, woodI = -1;
    for (size_t i = 0; i < rigid.bodies.size(); ++i) {
        MaterialId m = MATERIAL_EMPTY;
        for (int li : rigid.bodies[i].occupiedLocal) {
            m = rigid.bodies[i].mask[static_cast<size_t>(li)];
            break;
        }
        if (m == MATERIAL_METAL) metalI = static_cast<int>(i);
        if (m == MATERIAL_WOOD) woodI = static_cast<int>(i);
    }
    heatEnd(metalI, 1);
    heatEnd(woodI, 1);
    tickN(60);
    auto probeTemp = [&](int bi, int xTarget) {
        if (bi < 0) return AMBIENT_TEMPERATURE_K;
        RigidBody const &b = rigid.bodies[static_cast<size_t>(bi)];
        float sum = 0.0f;
        int n = 0;
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW;
            if (lx == xTarget || lx == xTarget + 1) {
                sum += rigidPixelTempK(b, li);
                ++n;
            }
        }
        return n ? sum / static_cast<float>(n) : AMBIENT_TEMPERATURE_K;
    };
    float tMetalFar = probeTemp(metalI, 6);
    float tWoodFar = probeTemp(woodI, 6);
    bool passCond = tMetalFar > tWoodFar + 12.0f;
    emit("metal_vs_wood", 60, tMetalFar, tWoodFar, passCond);

    // 3. Moving hot liquid — heat follows water
    resetBare();
    for (int x = 20; x <= 180; ++x) fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 100))] = 1;
    for (int y = 20; y <= 100; ++y) {
        fluid.solid[static_cast<size_t>(FluidEngine::ci(20, y))] = 1;
        fluid.solid[static_cast<size_t>(FluidEngine::ci(180, y))] = 1;
    }
    seedAmbient(fluid, rigid, gas);
    for (int y = 40; y <= 90; ++y) for (int x = 30; x <= 70; ++x) setLiquid(x, y, 1.0f, 360.0f);
    fluid.rebuildActivityAndMetrics();
    tickN(90);
    double trailHeat = 0.0, waterHeat = 0.0;
    for (int y = 40; y <= 90; ++y) for (int x = 30; x <= 70; ++x) {
        int i = FluidEngine::ci(x, y);
        if (fluid.fill[static_cast<size_t>(i)] < MIN_ACTIVE_FILL)
            trailHeat += fluid.liquidHeat[static_cast<size_t>(i)];
    }
    for (int i = 0; i < GW * GH; ++i) {
        if (fluid.fill[static_cast<size_t>(i)] >= MIN_ACTIVE_FILL) waterHeat += fluid.liquidHeat[static_cast<size_t>(i)];
    }
    bool passMoveLiq = trailHeat < 0.02 * std::max(1.0, waterHeat);
    emit("moving_hot_liquid", 90, waterHeat, trailHeat, passMoveLiq);

    // 4. Moving rigid body — heat stays on local pixels
    resetBare();
    paintBodyRect(MATERIAL_METAL, 80, 40, 95, 55, false);
    seedAmbient(fluid, rigid, gas);
    int hotLocal = -1;
    if (!rigid.bodies.empty()) {
        RigidBody &b = rigid.bodies[0];
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW, ly = li / b.maskW;
            if (lx <= 1 && ly <= 1) {
                b.heat[static_cast<size_t>(li)] = energyFromTemp(rigidPixelCapacity(b, li), 500.0f);
                hotLocal = li;
            }
        }
        b.omega = 4.0f;
        b.vx = 6.0f;
        wakeRect(static_cast<int>(b.aabbX0) - 2, static_cast<int>(b.aabbY0) - 2,
            static_cast<int>(b.aabbX1) + 2, static_cast<int>(b.aabbY1) + 2);
    }
    tickN(8);
    float tHot = AMBIENT_TEMPERATURE_K, tOther = 0.0f;
    int nOther = 0;
    if (!rigid.bodies.empty() && hotLocal >= 0) {
        RigidBody const &b = rigid.bodies[0];
        int idx = std::min(hotLocal, static_cast<int>(b.heat.size()) - 1);
        tHot = rigidPixelTempK(b, idx);
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW, ly = li / b.maskW;
            if (lx > 3 && ly > 3) { tOther += rigidPixelTempK(b, li); ++nOther; }
        }
    }
    float tO = nOther ? tOther / static_cast<float>(nOther) : AMBIENT_TEMPERATURE_K;
    bool passMoveRigid = tHot > tO + 30.0f;
    emit("moving_rigid", 8, tHot, tO, passMoveRigid);

    // 5. Fracture inherits pixel temperatures
    resetBare();
    paintBodyRect(MATERIAL_GLASS, 70, 50, 100, 58, true);
    seedAmbient(fluid, rigid, gas);
    if (!rigid.bodies.empty()) {
        RigidBody &b = rigid.bodies[0];
        for (int li : b.occupiedLocal) {
            int lx = li % b.maskW;
            float t = lx < b.maskW / 2 ? 420.0f : 293.15f;
            b.heat[static_cast<size_t>(li)] = energyFromTemp(rigidPixelCapacity(b, li), t);
        }
        int midX = static_cast<int>(b.x);
        int y0 = static_cast<int>(std::floor(b.aabbY0));
        int y1 = static_cast<int>(std::ceil(b.aabbY1));
        rigid.eraseLine(midX, y0, midX, y1, 0, fluid);
        rigid.syncOccupancy(fluid);
    }
    tickN(5);
    float tFragHot = 0, tFragCold = 0;
    int nH = 0, nC = 0;
    for (RigidBody const &b : rigid.bodies) {
        for (int li : b.occupiedLocal) {
            float t = rigidPixelTempK(b, li);
            if (t > 350.0f) { tFragHot += t; ++nH; }
            else { tFragCold += t; ++nC; }
        }
    }
    bool passFrac = nH > 0 && nC > 0 && (tFragHot / nH) > 370.0f && (tFragCold / nC) < 330.0f;
    emit("rigid_fracture", 5, nH ? tFragHot / nH : 0, nC ? tFragCold / nC : 0, passFrac);

    // 6. Thermal sleeping at uniform ambient
    resetBare();
    for (int x = 20; x <= 180; ++x) fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 110))] = 1;
    seedAmbient(fluid, rigid, gas);
    tickN(40);
    bool passSleep = work.activeChunks == 0 && lastStepMs < 0.4;
    emit("thermal_sleep", 40, lastStepMs, work.activeChunks, passSleep);

    // 7. Closed-system energy conservation
    resetBare();
    for (int x = 50; x <= 120; ++x) {
        fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 70))] = 1;
        fluid.solid[static_cast<size_t>(FluidEngine::ci(x, 40))] = 1;
    }
    for (int y = 40; y <= 70; ++y) {
        fluid.solid[static_cast<size_t>(FluidEngine::ci(50, y))] = 1;
        fluid.solid[static_cast<size_t>(FluidEngine::ci(120, y))] = 1;
    }
    seedAmbient(fluid, rigid, gas);
    for (int y = 45; y <= 65; ++y) for (int x = 55; x <= 85; ++x) setLiquid(x, y, 1.0f, 330.0f);
    for (int y = 45; y <= 65; ++y) for (int x = 90; x <= 115; ++x) setLiquid(x, y, 1.0f, 280.0f);
    fluid.rebuildActivityAndMetrics();
    double eClosed0 = totalThermalEnergy(fluid, rigid, gas);
    tickN(120);
    double eClosed1 = totalThermalEnergy(fluid, rigid, gas);
    double rel = std::abs(eClosed1 - eClosed0) / std::max(1.0, std::abs(eClosed0));
    ThermalWorldStats closedStats = collectStats(fluid, rigid, gas);
    bool passCons = rel < 0.02 && closedStats.nNan == 0 && closedStats.nInf == 0 && closedStats.nNegK == 0;
    emit("energy_conservation", 120, eClosed1, rel, passCons);
    emit("energy_conservation_T", 120, closedStats.tMin, closedStats.tMax, passCons);

    // 8. One-step conduction must not leap past equilibrium.
    {
        float cap = 1000.0f;
        float eA = energyFromTemp(cap, 400.0f);
        float eB = energyFromTemp(cap, 280.0f);
        float sum0 = eA + eB;
        exchangeThermalEnergy(eA, cap, 50.0f, eB, cap, 50.0f, 50.0f, 1.0f, 25000.0f);
        float Ta = tempFromEnergy(eA, cap);
        float Tb = tempFromEnergy(eB, cap);
        float sum1 = eA + eB;
        bool noCross = Ta + 1.0e-3f >= Tb;
        bool cons = std::abs(sum1 - sum0) / std::max(1.0f, std::abs(sum0)) < 1.0e-5f;
        bool nearEq = std::abs(Ta - Tb) < 0.05f;
        bool passOver = noCross && cons && nearEq;
        emit("conduction_overshoot", 1, Ta, Tb, passOver);
    }

    // 9. Ambient world with gas: measure real T spread (visualization diagnosis).
    {
        rigid.clear();
        fluid.resetWorld();
        gas.config.simMode = GasSimMode::Full;
        gas.resetAmbient(fluid);
        config.enabled = true;
        config.intervalTicks = 1;
        seedAmbient(fluid, rigid, gas);
        ThermalWorldStats seeded = collectStats(fluid, rigid, gas);
        tickN(60);
        ThermalWorldStats idle = collectStats(fluid, rigid, gas);
        bool passSeed = seeded.nNan == 0 && seeded.nInf == 0 && seeded.nNegK == 0
            && seeded.leftoverLiquidCells == 0 && seeded.leftoverSolidCells == 0
            && seeded.maxAbsDeltaAmbient < 0.05f;
        bool passIdle = idle.nNan == 0 && idle.nInf == 0 && idle.nNegK == 0
            && idle.leftoverLiquidCells == 0 && idle.maxAbsDeltaAmbient < 0.25f;
        emit("ambient_seed_spread", 0, seeded.maxAbsDeltaAmbient, seeded.tMean, passSeed);
        emit("ambient_idle_spread", 60, idle.maxAbsDeltaAmbient, idle.tMean, passIdle);

        std::ofstream stats(miscFile("thermal_stats.txt"), std::ios::trunc);
        stats << "PIPACE thermal numerical dump\n";
        stats << "ambient_K " << AMBIENT_TEMPERATURE_K << "\n";
        stats << "viz_deadband_K " << TEMP_VIZ_DEADBAND_K << "\n\n";
        auto dump = [&](char const *label, ThermalWorldStats const &s) {
            stats << "[" << label << "]\n";
            stats << "energy_J " << s.energy << "\n";
            stats << "t_min " << s.tMin << " t_max " << s.tMax << " t_mean " << s.tMean << "\n";
            stats << "max_|dT_amb| " << s.maxAbsDeltaAmbient << " n_away_deadband " << s.nAwayFromAmbient << "\n";
            stats << "n_gas " << s.nGas << " n_liquid " << s.nLiquid << " n_wall " << s.nWall
                  << " n_rigid " << s.nRigid << " n_empty " << s.nEmpty << "\n";
            if (s.nGas) stats << "gas_T " << s.tMinGas << " .. " << s.tMaxGas << "\n";
            if (s.nLiquid) stats << "liquid_T " << s.tMinLiquid << " .. " << s.tMaxLiquid
                << " fill_at_min " << s.tMinLiquidFill << " fill_at_max " << s.tMaxLiquidFill << "\n";
            if (s.nWall) stats << "wall_T " << s.tMinSolid << " .. " << s.tMaxSolid << "\n";
            if (s.nRigid) stats << "rigid_T " << s.tMinRigid << " .. " << s.tMaxRigid << "\n";
            stats << "nan " << s.nNan << " inf " << s.nInf << " negK " << s.nNegK << " near0K " << s.nZeroK << "\n";
            stats << "leftover_liquid_cells " << s.leftoverLiquidCells
                  << " leftover_liquid_J " << s.leftoverLiquidEnergy << "\n";
            stats << "leftover_solid_cells " << s.leftoverSolidCells
                  << " leftover_solid_J " << s.leftoverSolidEnergy << "\n\n";
        };
        dump("seeded_ambient_with_gas", seeded);
        dump("idle_60_ticks", idle);
        dump("closed_conservation_after_120", closedStats);
        stats << "closed_rel_energy_error " << rel << "\n";
        stats << "closed_energy_0 " << eClosed0 << " closed_energy_1 " << eClosed1 << "\n";
        stats.flush();
    }

    out.flush();
}

void ThermalEngine::runSpreadDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas) {
    std::ofstream out(miscFile("thermal_spread_diag.tsv"), std::ios::trunc);
    out << std::setprecision(8);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto f8 = [](double v) {
        std::ostringstream o;
        o << std::setprecision(8) << v;
        return o.str();
    };

    auto worldTick = [&]() {
        rigid.step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gas.simulationTick(fluid);
        rigid.gatherFluidForces(fluid);
        gas.applyPressureForces(rigid, fluid);
        simulationTick(fluid, rigid, gas, PHYSICS_DT);
    };

    auto box = [&](int x0, int y0, int x1, int y1) {
        for (int x = x0; x <= x1; ++x) {
            fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y0))] = 1;
            fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y1))] = 1;
        }
        for (int y = y0; y <= y1; ++y) {
            fluid.solid[static_cast<size_t>(FluidEngine::ci(x0, y))] = 1;
            fluid.solid[static_cast<size_t>(FluidEngine::ci(x1, y))] = 1;
        }
    };

    auto resetSealed = [&]() {
        rigid.clear();
        fluid.clearWorld();
        fluid.config.walledBorders = true;
        gas.config.boundary = GasBoundary::Sealed;
        gas.config.simMode = GasSimMode::Full;
        config.enabled = true;
        config.intervalTicks = 1;
        box(40, 18, 90, 95);
        gas.resetAmbient(fluid);
        seedAmbient(fluid, rigid, gas);
    };

    auto setGasT = [&](int x, int y, float tK) {
        if (!FluidEngine::inside(x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (fluid.solid[static_cast<size_t>(i)]) return;
        float cap = gasCapacity(gas, i);
        if (cap > MIN_THERMAL_CAPACITY)
            gas.heat[static_cast<size_t>(i)] = energyFromTemp(cap, tK);
        wakeCell(x, y);
        gas.wakeAt(x, y);
    };
    auto setWallT = [&](int x, int y, float tK) {
        if (!FluidEngine::inside(x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (!fluid.solid[static_cast<size_t>(i)]) return;
        float cap = wallCapacity(fluid, i);
        fluid.solidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, tK);
        wakeCell(x, y);
    };
    auto setLiqT = [&](int x, int y, float tK, float fillAmt) {
        if (!FluidEngine::inside(x, y) || fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y))]) return;
        int i = FluidEngine::ci(x, y);
        fluid.fill[static_cast<size_t>(i)] = fillAmt;
        float cap = liquidCapacity(fluid, i);
        fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, tK);
        fluid.expectedVolume += fillAmt;
        wakeCell(x, y);
        fluid.wakeChunkAtCell(x, y);
    };

    struct Probe {
        float center = AMBIENT_TEMPERATURE_K;
        float ring1 = AMBIENT_TEMPERATURE_K;
        float ring2 = AMBIENT_TEMPERATURE_K;
        double energy = 0.0;
        float tMax = AMBIENT_TEMPERATURE_K;
        float tMean = AMBIENT_TEMPERATURE_K;
        int nAway = 0;
        float comY = 0.0f;
        int activeChunks = 0;
        int nRing1 = 0, nRing2 = 0;
    };
    auto probeAt = [&](int cx, int cy, int x0, int y0, int x1, int y1) -> Probe {
        Probe p;
        p.center = sampleTemperatureK(fluid, rigid, gas, cx, cy);
        p.energy = totalThermalEnergy(fluid, rigid, gas);
        p.activeChunks = work.activeChunks;
        ThermalWorldStats st = collectStats(fluid, rigid, gas);
        p.tMax = st.tMax;
        p.tMean = st.tMean;
        p.nAway = st.nAwayFromAmbient;
        double r1 = 0.0, r2 = 0.0;
        int n1 = 0, n2 = 0;
        double wsum = 0.0, wy = 0.0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if (fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y))]) continue;
            float t = sampleTemperatureK(fluid, rigid, gas, x, y);
            int cheb = std::max(std::abs(x - cx), std::abs(y - cy));
            if (cheb == 1) { r1 += t; ++n1; }
            else if (cheb >= 2 && cheb <= 4) { r2 += t; ++n2; }
            float dT = t - AMBIENT_TEMPERATURE_K;
            if (dT > TEMP_VIZ_DEADBAND_K) {
                wsum += dT;
                wy += dT * static_cast<double>(y);
            }
        }
        p.nRing1 = n1;
        p.nRing2 = n2;
        p.ring1 = n1 ? static_cast<float>(r1 / n1) : AMBIENT_TEMPERATURE_K;
        p.ring2 = n2 ? static_cast<float>(r2 / n2) : AMBIENT_TEMPERATURE_K;
        p.comY = (wsum > 1.0e-6) ? static_cast<float>(wy / wsum) : static_cast<float>(cy);
        return p;
    };
    auto meanGasBand = [&](int x0, int y0, int x1, int y1) -> float {
        double s = 0.0;
        int n = 0;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            int i = FluidEngine::ci(x, y);
            if (fluid.solid[static_cast<size_t>(i)] || fluid.fill[static_cast<size_t>(i)] > 1.0e-6f) continue;
            if (gas.amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT) continue;
            s += gasTempK(gas, i);
            ++n;
        }
        return n ? static_cast<float>(s / n) : AMBIENT_TEMPERATURE_K;
    };
    auto row = [&](char const *name, int tick, Probe const &p) {
        out << "case\t" << name << '\t' << tick << '\t' << p.center << '\t' << p.ring1 << '\t'
            << p.ring2 << '\t' << p.energy << '\t' << p.tMax << '\t' << p.tMean << '\t'
            << p.nAway << '\t' << p.comY << '\t' << p.activeChunks << '\n';
    };

    out << "case\tname\ttick\tcenter_T\tring1_T\tring2_T\tenergy\tmax_T\tmean_T\tn_away\tcom_y\tactive_chunks\n";

    // 1. Localized hot gas in sealed ambient air
    resetSealed();
    int gx = 65, gy = 78;
    for (int y = gy - 3; y <= gy + 3; ++y) for (int x = gx - 3; x <= gx + 3; ++x) setGasT(x, y, 450.0f);
    Probe g0 = probeAt(gx, gy, 41, 19, 89, 94);
    float above0 = meanGasBand(gx - 6, gy - 16, gx + 6, gy - 6);
    row("hot_gas_patch", 0, g0);
    Probe g15{}, g30{}, g90{};
    for (int n = 1; n <= 90; ++n) {
        worldTick();
        if (n == 15) { g15 = probeAt(gx, gy, 41, 19, 89, 94); row("hot_gas_patch", 15, g15); }
        if (n == 30) { g30 = probeAt(gx, gy, 41, 19, 89, 94); row("hot_gas_patch", 30, g30); }
        if (n == 90) { g90 = probeAt(gx, gy, 41, 19, 89, 94); row("hot_gas_patch", 90, g90); }
    }
    float above90 = meanGasBand(gx - 6, gy - 16, gx + 6, gy - 6);
    double gRel = std::abs(g90.energy - g0.energy) / std::max(1.0, std::abs(g0.energy));
    emit("hot_gas_center_cools", g90.center + 15.0f < g0.center,
        "T0=" + f8(g0.center) + " T90=" + f8(g90.center));
    emit("hot_gas_neighbors_warm", above90 > above0 + 4.0f,
        "above0=" + f8(above0) + " above90=" + f8(above90));
    emit("hot_gas_region_expands", g90.nAway > g0.nAway + 8,
        "n0=" + std::to_string(g0.nAway) + " n90=" + std::to_string(g90.nAway));
    emit("hot_gas_sealed_energy", gRel < 0.02, "rel=" + f8(gRel));
    emit("hot_gas_rises", g90.comY + 1.5f < g0.comY,
        "com0=" + f8(g0.comY) + " com90=" + f8(g90.comY));
    emit("hot_gas_sleep_does_not_freeze",
        g15.activeChunks >= 1 && g30.activeChunks >= 2 && g90.nAway > g0.nAway,
        "c15=" + std::to_string(g15.activeChunks) + " c30=" + std::to_string(g30.activeChunks)
            + " nAway90=" + std::to_string(g90.nAway));

    // 2. Hot wall next to air
    resetSealed();
    int wx = 40, wy = 55;
    for (int y = wy - 2; y <= wy + 2; ++y) setWallT(wx, y, 480.0f);
    int ax = wx + 1, ay = wy;
    float wall0 = wallTempK(fluid, FluidEngine::ci(wx, wy));
    float air0 = gasTempK(gas, GasEngine::ci(ax, ay));
    double eWall0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float wall1 = wallTempK(fluid, FluidEngine::ci(wx, wy));
    float air1 = gasTempK(gas, GasEngine::ci(ax, ay));
    double eWall1 = totalThermalEnergy(fluid, rigid, gas);
    double wallRel = std::abs(eWall1 - eWall0) / std::max(1.0, std::abs(eWall0));
    emit("hot_wall_warms_air", air1 > air0 + 8.0f,
        "air0=" + f8(air0) + " air90=" + f8(air1));
    emit("hot_wall_cools", wall1 < wall0 - 0.002f,
        "wall0=" + f8(wall0) + " wall90=" + f8(wall1));
    emit("hot_wall_sealed_energy", wallRel < 0.02, "rel=" + f8(wallRel));

    // 3. Hot liquid next to air (below boiling). Rest on the sealed floor so it stays put.
    resetSealed();
    int lx = 65, ly = 94;
    for (int x = lx - 6; x <= lx + 6; ++x) setLiqT(x, ly, 340.0f, 1.0f);
    fluid.rebuildActivityAndMetrics();
    gas.handleWorldEdit(fluid);
    int gxL = lx, gyL = ly - 1;
    float liq0 = liquidTempK(fluid, FluidEngine::ci(lx, ly));
    float lair0 = gasTempK(gas, GasEngine::ci(gxL, gyL));
    double eLiq0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float liq1 = liquidTempK(fluid, FluidEngine::ci(lx, ly));
    float lair1 = gasTempK(gas, GasEngine::ci(gxL, gyL));
    if (!(lair1 > lair0 + 4.0f))
        lair1 = meanGasBand(lx - 6, ly - 4, lx + 6, ly - 1);
    double eLiq1 = totalThermalEnergy(fluid, rigid, gas);
    double liqRel = std::abs(eLiq1 - eLiq0) / std::max(1.0, std::abs(eLiq0));
    emit("hot_liquid_warms_air", lair1 > lair0 + 4.0f,
        "air0=" + f8(lair0) + " air90=" + f8(lair1));
    emit("hot_liquid_cools", liq1 < liq0 - 0.01f,
        "liq0=" + f8(liq0) + " liq90=" + f8(liq1));
    emit("hot_liquid_sealed_energy", liqRel < 0.02, "rel=" + f8(liqRel));

    // 4. Cold gas patch
    resetSealed();
    int cxC = 65, cyC = 50;
    for (int y = cyC - 3; y <= cyC + 3; ++y) for (int x = cxC - 3; x <= cxC + 3; ++x)
        setGasT(x, y, 250.0f);
    Probe c0 = probeAt(cxC, cyC, 41, 19, 89, 94);
    float below0 = meanGasBand(cxC - 6, cyC + 6, cxC + 6, cyC + 16);
    row("cold_gas_patch", 0, c0);
    for (int n = 0; n < 90; ++n) worldTick();
    Probe c90 = probeAt(cxC, cyC, 41, 19, 89, 94);
    float below90 = meanGasBand(cxC - 6, cyC + 6, cxC + 6, cyC + 16);
    row("cold_gas_patch", 90, c90);
    double cRel = std::abs(c90.energy - c0.energy) / std::max(1.0, std::abs(c0.energy));
    emit("cold_gas_center_warms", c90.center > c0.center + 8.0f,
        "T0=" + f8(c0.center) + " T90=" + f8(c90.center));
    emit("cold_gas_neighbors_cool", below90 < below0 - 2.0f,
        "below0=" + f8(below0) + " below90=" + f8(below90));
    emit("cold_gas_sealed_energy", cRel < 0.02, "rel=" + f8(cRel));

    // 5. Open-boundary gas mass heat + conduction accounting
    rigid.clear();
    fluid.clearWorld();
    fluid.config.walledBorders = false;
    gas.config.boundary = GasBoundary::OpenAmbient;
    gas.config.simMode = GasSimMode::Full;
    config.enabled = true;
    config.intervalTicks = 1;
    gas.resetAmbient(fluid);
    seedAmbient(fluid, rigid, gas);
    int ox = GW / 2, oy = 2;
    for (int y = oy; y <= oy + 5; ++y) for (int x = ox - 4; x <= ox + 4; ++x) setGasT(x, y, 420.0f);
    double eOpen0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 60; ++n) worldTick();
    double eOpen1 = totalThermalEnergy(fluid, rigid, gas);
    double netOpen = netExternalHeat(fluid, gas);
    double openAccounted = eOpen1 + netOpen;
    double openRel = std::abs(openAccounted - eOpen0) / std::max(1.0, std::abs(eOpen0));
    emit("open_boundary_heat_accounted",
        openRel < 0.03,
        "E0=" + f8(eOpen0) + " E1=" + f8(eOpen1)
            + " condEsc=" + f8(thermalEnergyEscaped) + " condEnt=" + f8(thermalEnergyEntered)
            + " gasEsc=" + f8(gas.escapedHeat) + " fluidEsc=" + f8(fluid.escapedHeat)
            + " rel=" + f8(openRel));

    auto meanEdgeWall = [&](int x, int y0, int y1) -> float {
        double s = 0.0;
        int n = 0;
        for (int y = y0; y <= y1; ++y) {
            int i = FluidEngine::ci(x, y);
            if (!fluid.solid[static_cast<size_t>(i)]) continue;
            s += wallTempK(fluid, i);
            ++n;
        }
        return n ? static_cast<float>(s / n) : AMBIENT_TEMPERATURE_K;
    };
    auto meanEdgeGas = [&](int x, int y0, int y1) -> float {
        double s = 0.0;
        int n = 0;
        for (int y = y0; y <= y1; ++y) {
            int i = FluidEngine::ci(x, y);
            if (fluid.solid[static_cast<size_t>(i)] || fluid.fill[static_cast<size_t>(i)] > 1.0e-6f) continue;
            if (gas.amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT) continue;
            s += gasTempK(gas, i);
            ++n;
        }
        return n ? static_cast<float>(s / n) : AMBIENT_TEMPERATURE_K;
    };
    auto placeEdgeWall = [&](int x, int y0, int y1, float tK) {
        for (int y = y0; y <= y1; ++y) {
            fluid.solid[static_cast<size_t>(FluidEngine::ci(x, y))] = 1;
        }
        seedAmbient(fluid, rigid, gas);
        for (int y = y0; y <= y1; ++y) setWallT(x, y, tK);
    };
    auto resetOpenBare = [&]() {
        rigid.clear();
        fluid.clearWorld();
        fluid.config.walledBorders = false;
        gas.config.boundary = GasBoundary::OpenAmbient;
        gas.config.simMode = GasSimMode::Full;
        config.enabled = true;
        config.intervalTicks = 1;
        gas.resetAmbient(fluid);
        seedAmbient(fluid, rigid, gas);
    };
    auto resetWalledBare = [&]() {
        rigid.clear();
        fluid.clearWorld();
        fluid.config.walledBorders = true;
        gas.config.boundary = GasBoundary::Sealed;
        gas.config.simMode = GasSimMode::Full;
        config.enabled = true;
        config.intervalTicks = 1;
        gas.resetAmbient(fluid);
        seedAmbient(fluid, rigid, gas);
    };
    auto accountOk = [&](double e0, double e1) {
        double rel = std::abs((e1 + netExternalHeat(fluid, gas)) - e0) / std::max(1.0, std::abs(e0));
        return rel;
    };

    // 6. OPEN BORDER HOT PATCH — wall on the map rim
    resetOpenBare();
    int ehx = 0, ehy0 = 50, ehy1 = 70;
    placeEdgeWall(ehx, ehy0, ehy1, 480.0f);
    float edgeHot0 = meanEdgeWall(ehx, ehy0, ehy1);
    double eEh0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float edgeHot1 = meanEdgeWall(ehx, ehy0, ehy1);
    double eEh1 = totalThermalEnergy(fluid, rigid, gas);
    double ehRel = accountOk(eEh0, eEh1);
    double ehCondEsc = thermalEnergyEscaped;
    emit("open_border_hot_cools", edgeHot1 < edgeHot0 - 0.002f && edgeHot1 > MIN_SAFE_TEMPERATURE_K,
        "T0=" + f8(edgeHot0) + " T90=" + f8(edgeHot1));
    emit("open_border_hot_energy_leaves",
        thermalEnergyEscaped > 1.0 && eEh1 < eEh0 - 1.0,
        "condEsc=" + f8(thermalEnergyEscaped) + " dE=" + f8(eEh1 - eEh0));
    emit("open_border_hot_accounted", ehRel < 0.03,
        "E0=" + f8(eEh0) + " E1=" + f8(eEh1) + " netExt=" + f8(netExternalHeat(fluid, gas))
            + " rel=" + f8(ehRel));

    // 7. OPEN BORDER COLD PATCH — cold gas on the rim
    resetOpenBare();
    int ecx = 0, ecy0 = 50, ecy1 = 70;
    for (int y = ecy0; y <= ecy1; ++y) setGasT(ecx, y, 250.0f);
    float edgeCold0 = meanEdgeGas(ecx, ecy0, ecy1);
    double eEc0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float edgeCold1 = meanEdgeGas(ecx, ecy0, ecy1);
    double eEc1 = totalThermalEnergy(fluid, rigid, gas);
    double ecRel = accountOk(eEc0, eEc1);
    emit("open_border_cold_warms", edgeCold1 > edgeCold0 + 2.0f,
        "T0=" + f8(edgeCold0) + " T90=" + f8(edgeCold1));
    emit("open_border_cold_energy_enters",
        thermalEnergyEntered > 0.1 && eEc1 > eEc0,
        "condEnt=" + f8(thermalEnergyEntered) + " dE=" + f8(eEc1 - eEc0)
            + " gasEsc=" + f8(gas.escapedHeat));
    emit("open_border_cold_accounted", ecRel < 0.03,
        "E0=" + f8(eEc0) + " E1=" + f8(eEc1) + " netExt=" + f8(netExternalHeat(fluid, gas))
            + " rel=" + f8(ecRel));

    // 8. WALLED BORDER HOT PATCH — same wall, no direct external loss
    resetWalledBare();
    placeEdgeWall(ehx, ehy0, ehy1, 480.0f);
    float walledHot0 = meanEdgeWall(ehx, ehy0, ehy1);
    double eWh0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float walledHot1 = meanEdgeWall(ehx, ehy0, ehy1);
    double eWh1 = totalThermalEnergy(fluid, rigid, gas);
    double walledRel = std::abs(eWh1 - eWh0) / std::max(1.0, std::abs(eWh0));
    emit("walled_border_hot_no_external",
        thermalEnergyEscaped < 1.0e-6 && thermalEnergyEntered < 1.0e-6
            && std::abs(gas.escapedHeat) < 1.0e-6 && std::abs(fluid.escapedHeat) < 1.0e-6,
        "condEsc=" + f8(thermalEnergyEscaped) + " condEnt=" + f8(thermalEnergyEntered)
            + " gasEsc=" + f8(gas.escapedHeat));
    emit("walled_border_hot_conserved", walledRel < 0.02,
        "rel=" + f8(walledRel) + " T0=" + f8(walledHot0) + " T90=" + f8(walledHot1));

    // 9. INTERIOR CONTROL — identical hot wall away from the rim
    resetOpenBare();
    int ihx = 100, ihy0 = 50, ihy1 = 70;
    placeEdgeWall(ihx, ihy0, ihy1, 480.0f);
    float intHot0 = meanEdgeWall(ihx, ihy0, ihy1);
    double eIh0 = totalThermalEnergy(fluid, rigid, gas);
    for (int n = 0; n < 90; ++n) worldTick();
    float intHot1 = meanEdgeWall(ihx, ihy0, ihy1);
    double eIh1 = totalThermalEnergy(fluid, rigid, gas);
    double ihRel = accountOk(eIh0, eIh1);
    emit("interior_hot_cools_in_map", intHot1 < intHot0 - 0.002f,
        "T0=" + f8(intHot0) + " T90=" + f8(intHot1));
    emit("interior_no_direct_edge_dump",
        thermalEnergyEscaped < ehCondEsc * 0.05,
        "condEsc=" + f8(thermalEnergyEscaped) + " edgeCondEsc=" + f8(ehCondEsc));
    emit("interior_hot_accounted", ihRel < 0.03,
        "rel=" + f8(ihRel) + " netExt=" + f8(netExternalHeat(fluid, gas)));

    ThermalWorldStats endStats = collectStats(fluid, rigid, gas);
    emit("no_nan_inf_negK", endStats.nNan == 0 && endStats.nInf == 0 && endStats.nNegK == 0, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
