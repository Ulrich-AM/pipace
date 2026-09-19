#include "world/PhaseChangeEngine.h"

#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/LiquidMixtureProperties.h"
#include "substance/PhaseTransfer.h"
#include "thermal/ThermalEngine.h"
#include "world/WaterPhaseChange.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

constexpr int kDx[4] = {0, -1, 1, 0};
constexpr int kDy[4] = {-1, 0, 0, 1};
// Time-based liquid/gas conversion cap. 10.5 fill/s matches 0.35 fill/tick at 30 Hz.
constexpr float kMaxLiquidGasFillPerSec = 10.5f;
constexpr float kSolverSafetyAtm = 1.0e4f;
constexpr float kMinFillMove = 1.0e-6f;
constexpr float kBoilEpsK = 0.05f;
constexpr float kSatHystRel = 0.002f;
constexpr int kSearchLimit = 96;

struct PhaseScratch {
    std::vector<int> q;
    std::vector<uint32_t> seen;
    uint32_t gen = 1;

    void prepareSeen() {
        size_t n = static_cast<size_t>(GW * GH);
        if (seen.size() != n) {
            seen.assign(n, 0);
            gen = 1;
        } else {
            ++gen;
            if (gen == 0) {
                std::fill(seen.begin(), seen.end(), 0);
                gen = 1;
            }
        }
        q.clear();
    }

    bool mark(int i) {
        uint32_t &slot = seen[static_cast<size_t>(i)];
        if (slot == gen) return false;
        slot = gen;
        return true;
    }
};

PhaseScratch &liquidGasScratch() {
    static PhaseScratch s;
    return s;
}

bool isBlockedSolid(FluidEngine const &fluid, RigidBodyEngine const &rigid, int x, int y) {
    if (!FluidEngine::inside(x, y)) return true;
    int i = FluidEngine::ci(x, y);
    if (fluid.solid[static_cast<size_t>(i)]) return true;
    int body = rigid.occupant[static_cast<size_t>(i)];
    return body >= 0;
}

void wakeAllEngines(FluidEngine &fluid, GasEngine &gas, ThermalEngine &thermal, int x, int y) {
    fluid.wakeChunkAtCell(x, y);
    fluid.wakeThermalAt(x, y);
    gas.wakeAt(x, y);
    thermal.wakeCell(x, y);
}

uint32_t destTieKey(int x, int y, uint32_t salt) {
    return static_cast<uint32_t>(x) * 73856093u
        ^ static_cast<uint32_t>(y) * 19349663u
        ^ salt * 83492791u;
}

bool destBetter(float score, uint32_t key, float bestScore, uint32_t bestKey) {
    if (score > bestScore + 1.0e-5f) return true;
    if (score < bestScore - 1.0e-5f) return false;
    return key < bestKey;
}

float cellTempK(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y)
{
    if (!FluidEngine::inside(x, y)) return AMBIENT_TEMPERATURE_K;
    return ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
}

bool cellIsSolidSurface(FluidEngine const &fluid, RigidBodyEngine const &rigid, int x, int y) {
    if (!FluidEngine::inside(x, y)) return true;
    int i = FluidEngine::ci(x, y);
    return fluid.solid[static_cast<size_t>(i)] || rigid.occupant[static_cast<size_t>(i)] >= 0;
}

float coolSurfaceBonus(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y, float Tref)
{
    if (!(Tref > 1.0f)) return 0.0f;
    float bonus = 0.0f;
    for (int n = 0; n < 4; ++n) {
        int nx = x + kDx[n], ny = y + kDy[n];
        if (!cellIsSolidSurface(fluid, rigid, nx, ny)) continue;
        float Ts = cellTempK(fluid, rigid, gas, nx, ny);
        if (!(Ts + 0.5f < Tref)) continue;
        float span = std::clamp((Tref - Ts) / 40.0f, 0.0f, 1.0f);
        bonus += 3.2f * span;
    }
    return bonus;
}

float gasCellPressureAtm(GasEngine const &gas, int i, float vol) {
    if (!(vol >= GAS_MIN_VOLUME)) return 0.0f;
    return gasPressureAtmFromState(gas.amount[static_cast<size_t>(i)], vol, gas.cellTemperatureK(i));
}

float gasAmountRoomForPressureAtm(GasEngine const &gas, int i, float vol, float maxAtm) {
    if (!(vol >= GAS_MIN_VOLUME)) return 0.0f;
    float capAmt = gasAmountFromPressureAtm(maxAtm, vol, gas.cellTemperatureK(i));
    return capAmt - gas.amount[static_cast<size_t>(i)];
}

bool gasCellCanAcceptComponent(GasEngine const &gas, int index, SubstanceId id) {
    GasComponentView view = gas.gasComponents(index);
    GasComponent src{id, 1.0f};
    return gasPayloadCanMerge(view.items, view.count, &src, 1);
}

bool liquidCellAcceptsCondensate(FluidEngine const &fluid, int i, SubstanceId id) {
    float fill = fluid.fill[static_cast<size_t>(i)];
    if (fill < kMinFillMove) return true;
    SubstanceId pureId = SUBSTANCE_NONE;
    float mainAmt = 0.0f;
    float otherAmt = 0.0f;
    LiquidComponentView view = fluid.liquidComponents(i);
    for (int n = 0; n < view.count; ++n) {
        SubstanceId sid = view.items[n].id;
        float amt = view.items[n].amount;
        if (!validLiquidComponentId(sid) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        if (pureId == SUBSTANCE_NONE || sid == pureId) {
            pureId = sid;
            mainAmt += amt;
        } else {
            otherAmt += amt;
        }
    }
    if (pureId == SUBSTANCE_NONE) return true;
    if (otherAmt > kMeaningfulLiquidComponent) return false;
    return pureId == id;
}

bool resolvePureLiquidForVaporization(FluidEngine const &fluid, int i,
    SubstanceId &outId, float &outFill)
{
    outId = SUBSTANCE_NONE;
    outFill = 0.0f;
    float fill = fluid.fill[static_cast<size_t>(i)];
    if (fill < kMinFillMove) return false;
    LiquidComponentView view = fluid.liquidComponents(i);
    SubstanceId found = SUBSTANCE_NONE;
    float mainAmt = 0.0f;
    float otherAmt = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId sid = view.items[n].id;
        float amt = view.items[n].amount;
        if (!validLiquidComponentId(sid) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        if (found == SUBSTANCE_NONE || sid == found) {
            found = sid;
            mainAmt += amt;
        } else {
            otherAmt += amt;
        }
    }
    if (found == SUBSTANCE_NONE) return false;
    if (otherAmt > kMeaningfulLiquidComponent) return false;
    if (mainAmt + kMeaningfulLiquidComponent < fill) return false;
    if (!supportsLiveLiquidGasTransition(found)) return false;
    outId = found;
    outFill = fluid.liquidComponentAmount(i, found);
    return outFill >= kMinFillMove;
}

float localAmbientPressurePa(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y)
{
    (void)rigid;
    float pref = gas.config.referencePressurePa;
    float bestAtm = -1.0f;
    auto consider = [&](int cx, int cy) {
        if (!gas.isAccessible(fluid, cx, cy)) return;
        int i = GasEngine::ci(cx, cy);
        float vol = gas.availableVolume(fluid, cx, cy);
        if (!(vol >= GAS_MIN_VOLUME)) return;
        float pAtm = gasCellPressureAtm(gas, i, vol);
        if (bestAtm < 0.0f || pAtm < bestAtm) bestAtm = pAtm;
    };
    consider(x, y);
    for (int n = 0; n < 4; ++n) consider(x + kDx[n], y + kDy[n]);
    if (bestAtm < 0.0f) return 1.0f;
    return std::max(1.0f, bestAtm * pref);
}

int findVaporDest(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    SubstanceId id, int sx, int sy, float needAmt, uint32_t salt, float maxAtm,
    bool *hitSafety, bool *hitEq)
{
    if (needAmt <= GAS_MIN_AMOUNT) return -1;
    if (hitSafety) *hitSafety = false;
    if (hitEq) *hitEq = false;
    if (!(maxAtm > 0.0f) || !std::isfinite(maxAtm)) maxAtm = kSolverSafetyAtm;
    auto room = [&](int x, int y) -> float {
        if (isBlockedSolid(fluid, rigid, x, y)) return 0.0f;
        if (!gas.isAccessible(fluid, x, y)) return 0.0f;
        int i = GasEngine::ci(x, y);
        if (!gasCellCanAcceptComponent(gas, i, id)) return 0.0f;
        float vol = gas.availableVolume(fluid, x, y);
        if (vol < GAS_MIN_VOLUME) return 0.0f;
        float thermo = gasAmountRoomForPressureAtm(gas, i, vol, maxAtm);
        float safety = gasAmountRoomForPressureAtm(gas, i, vol, kSolverSafetyAtm);
        return std::min(thermo, safety);
    };
    int best = -1;
    float bestRoom = 0.0f;
    uint32_t bestKey = ~0u;
    auto consider = [&](int x, int y) {
        float r = room(x, y);
        if (!(r > GAS_MIN_AMOUNT)) return;
        uint32_t key = destTieKey(x, y, salt);
        if (destBetter(r, key, bestRoom, bestKey)) {
            bestRoom = r;
            bestKey = key;
            best = GasEngine::ci(x, y);
        }
    };
    for (int n = 0; n < 4; ++n) consider(sx + kDx[n], sy + kDy[n]);
    consider(sx, sy);
    if (best >= 0 && bestRoom >= std::min(needAmt, 1.0e-4f)) return best;

    PhaseScratch &sc = liquidGasScratch();
    sc.prepareSeen();
    sc.q.reserve(kSearchLimit);
    auto push = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return;
        if (isBlockedSolid(fluid, rigid, x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (!sc.mark(i)) return;
        sc.q.push_back(i);
    };
    if (FluidEngine::inside(sx, sy)) {
        int si = FluidEngine::ci(sx, sy);
        if (sc.mark(si)) sc.q.push_back(si);
    }
    for (size_t head = 0; head < sc.q.size() && static_cast<int>(head) < kSearchLimit; ++head) {
        int i = sc.q[head];
        int x = i % GW, y = i / GW;
        consider(x, y);
        if (best >= 0 && bestRoom >= needAmt) return best;
        for (int n = 0; n < 4; ++n) push(x + kDx[n], y + kDy[n]);
    }
    if (best >= 0 && bestRoom > GAS_MIN_AMOUNT) return best;
    bool anyAccess = gas.isAccessible(fluid, sx, sy);
    for (int n = 0; n < 4 && !anyAccess; ++n)
        anyAccess = gas.isAccessible(fluid, sx + kDx[n], sy + kDy[n]);
    if (anyAccess) {
        bool anySafety = false;
        auto probe = [&](int x, int y) {
            if (!gas.isAccessible(fluid, x, y)) return;
            float vol = gas.availableVolume(fluid, x, y);
            if (vol < GAS_MIN_VOLUME) return;
            float amt = gas.amount[static_cast<size_t>(GasEngine::ci(x, y))];
            if (amt >= gasAmountFromPressureAtm(kSolverSafetyAtm, vol, gas.cellTemperatureK(GasEngine::ci(x, y)))
                - GAS_MIN_AMOUNT) anySafety = true;
        };
        probe(sx, sy);
        for (int n = 0; n < 4; ++n) probe(sx + kDx[n], sy + kDy[n]);
        if (anySafety && hitSafety) *hitSafety = true;
        else if (hitEq) *hitEq = true;
    }
    return -1;
}

int findLiquidDest(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    SubstanceId id, int sx, int sy, uint32_t salt, bool preferCoolSurfaces, float Tref)
{
    auto score = [&](int x, int y) -> float {
        if (isBlockedSolid(fluid, rigid, x, y)) return -1.0f;
        int i = FluidEngine::ci(x, y);
        float fill = fluid.fill[static_cast<size_t>(i)];
        float room = 1.0f - fill;
        if (room < kMinFillMove) return -1.0f;
        if (!liquidCellAcceptsCondensate(fluid, i, id)) return -1.0f;
        float s = room;
        if (fill > kMinFillMove) s += 4.0f + fill;
        if (preferCoolSurfaces) s += coolSurfaceBonus(fluid, rigid, gas, x, y, Tref);
        if (y > sy) s += 1.5f;
        return s;
    };
    int best = -1;
    float bestS = -1.0f;
    uint32_t bestKey = ~0u;
    auto consider = [&](int x, int y) {
        float s = score(x, y);
        if (s < 0.0f) return;
        uint32_t key = destTieKey(x, y, salt);
        if (destBetter(s, key, bestS, bestKey)) {
            bestS = s;
            bestKey = key;
            best = FluidEngine::ci(x, y);
        }
    };
    consider(sx, sy);
    for (int n = 0; n < 4; ++n) consider(sx + kDx[n], sy + kDy[n]);
    if (best >= 0 && bestS > 0.0f) return best;

    PhaseScratch &sc = liquidGasScratch();
    sc.prepareSeen();
    sc.q.reserve(kSearchLimit);
    auto push = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return;
        if (isBlockedSolid(fluid, rigid, x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (!sc.mark(i)) return;
        sc.q.push_back(i);
    };
    if (FluidEngine::inside(sx, sy)) {
        int si = FluidEngine::ci(sx, sy);
        if (sc.mark(si)) sc.q.push_back(si);
    }
    for (size_t head = 0; head < sc.q.size() && static_cast<int>(head) < kSearchLimit; ++head) {
        int i = sc.q[head];
        int x = i % GW, y = i / GW;
        consider(x, y);
        for (int n = 0; n < 4; ++n) push(x + kDx[n], y + kDy[n]);
    }
    return (best >= 0 && bestS > 0.0f) ? best : -1;
}

float heatRoomToK(float energy, float cap, float tMax) {
    if (cap < MIN_THERMAL_CAPACITY) return 0.0f;
    float T = tempFromEnergy(energy, cap);
    if (!(T < tMax)) return 0.0f;
    return cap * (tMax - T);
}

void addHeatToK(float &energy, float cap, float tMax, float &budget) {
    if (!(budget > 0.0f)) return;
    float room = heatRoomToK(energy, cap, tMax);
    float add = std::min(budget, room);
    if (!(add > 0.0f)) return;
    energy += add;
    budget -= add;
}

float heatAvailableAboveK(float energy, float cap, float tMin) {
    if (cap < MIN_THERMAL_CAPACITY) return 0.0f;
    float T = tempFromEnergy(energy, cap);
    if (!(T > tMin)) return 0.0f;
    return cap * (T - tMin);
}

void takeHeatAboveK(float &energy, float cap, float tMin, float &budget) {
    if (!(budget > 0.0f)) return;
    float avail = heatAvailableAboveK(energy, cap, tMin);
    float take = std::min(budget, avail);
    if (!(take > 0.0f)) return;
    energy -= take;
    budget -= take;
    if (!std::isfinite(energy) || energy < 0.0f) energy = 0.0f;
}

void dumpHeatBudget(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int x, int y, float tMax, float &budget)
{
    if (!(budget > 0.0f) || !FluidEngine::inside(x, y)) return;
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
                addHeatToK(b.heat[static_cast<size_t>(li)],
                    ThermalEngine::rigidPixelCapacity(b, li), tMax, budget);
                thermal.wakeCell(x, y);
            }
        }
        return;
    }
    if (fluid.solid[static_cast<size_t>(i)]) {
        addHeatToK(fluid.solidHeat[static_cast<size_t>(i)],
            ThermalEngine::wallCapacity(fluid, i), tMax, budget);
        thermal.wakeCell(x, y);
        return;
    }
    if (fluid.fill[static_cast<size_t>(i)] > 1.0e-8f) {
        addHeatToK(fluid.liquidHeat[static_cast<size_t>(i)],
            ThermalEngine::liquidCapacity(fluid, i), tMax, budget);
        return;
    }
    if (gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT) {
        addHeatToK(gas.heat[static_cast<size_t>(i)],
            ThermalEngine::gasCapacity(gas, i), tMax, budget);
        gas.wakeAt(x, y);
    }
}

void takeHeatAtCell(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int x, int y, float tMin, float &budget)
{
    if (!(budget > 0.0f) || !FluidEngine::inside(x, y)) return;
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
                takeHeatAboveK(b.heat[static_cast<size_t>(li)],
                    ThermalEngine::rigidPixelCapacity(b, li), tMin, budget);
                thermal.wakeCell(x, y);
            }
        }
        return;
    }
    if (fluid.solid[static_cast<size_t>(i)]) {
        takeHeatAboveK(fluid.solidHeat[static_cast<size_t>(i)],
            ThermalEngine::wallCapacity(fluid, i), tMin, budget);
        thermal.wakeCell(x, y);
        return;
    }
    if (fluid.fill[static_cast<size_t>(i)] > 1.0e-8f) {
        takeHeatAboveK(fluid.liquidHeat[static_cast<size_t>(i)],
            ThermalEngine::liquidCapacity(fluid, i), tMin, budget);
        return;
    }
    if (gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT) {
        takeHeatAboveK(gas.heat[static_cast<size_t>(i)],
            ThermalEngine::gasCapacity(gas, i), tMin, budget);
        gas.wakeAt(x, y);
    }
}

float heatAvailableAtCell(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y, float tMin)
{
    if (!FluidEngine::inside(x, y)) return 0.0f;
    int i = FluidEngine::ci(x, y);
    int body = rigid.occupant[static_cast<size_t>(i)];
    if (body >= 0 && body < static_cast<int>(rigid.bodies.size())) {
        RigidBody const &b = rigid.bodies[static_cast<size_t>(body)];
        float lx, ly;
        RigidBodyEngine::worldToLocal(b, x + 0.5f, y + 0.5f, lx, ly);
        int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
        if (RigidBodyEngine::maskOccupied(b, ix, iy)) {
            int li = iy * b.maskW + ix;
            if (li >= 0 && li < static_cast<int>(b.heat.size()))
                return heatAvailableAboveK(b.heat[static_cast<size_t>(li)],
                    ThermalEngine::rigidPixelCapacity(b, li), tMin);
        }
        return 0.0f;
    }
    if (fluid.solid[static_cast<size_t>(i)])
        return heatAvailableAboveK(fluid.solidHeat[static_cast<size_t>(i)],
            ThermalEngine::wallCapacity(fluid, i), tMin);
    if (fluid.fill[static_cast<size_t>(i)] > 1.0e-8f)
        return heatAvailableAboveK(fluid.liquidHeat[static_cast<size_t>(i)],
            ThermalEngine::liquidCapacity(fluid, i), tMin);
    if (gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT)
        return heatAvailableAboveK(gas.heat[static_cast<size_t>(i)],
            ThermalEngine::gasCapacity(gas, i), tMin);
    return 0.0f;
}

float heatAvailableAround(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y, float tMin, bool includeCenter)
{
    float q = 0.0f;
    if (includeCenter) q += heatAvailableAtCell(fluid, rigid, gas, x, y, tMin);
    for (int n = 0; n < 4; ++n)
        q += heatAvailableAtCell(fluid, rigid, gas, x + kDx[n], y + kDy[n], tMin);
    return q;
}

void takeHeatAround(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int x, int y, float tMin, float &budget, bool includeCenter)
{
    if (includeCenter) takeHeatAtCell(fluid, rigid, gas, thermal, x, y, tMin, budget);
    for (int n = 0; n < 4 && budget > 0.0f; ++n)
        takeHeatAtCell(fluid, rigid, gas, thermal, x + kDx[n], y + kDy[n], tMin, budget);
}

} // namespace

bool supportsLiveLiquidGasTransition(SubstanceId id) {
    if (id == SUBSTANCE_NONE || !validSubstance(id)) return false;
    if (!supportsPhase(id, MatterPhase::Liquid) || !supportsPhase(id, MatterPhase::Gas))
        return false;
    if (!hasFluidProperties(id) || !hasGasProperties(id)) return false;
    if (!canTransition(id, MatterPhase::Liquid, MatterPhase::Gas)) return false;
    float M = chemicalForSubstance(id).molarMass;
    if (!(M > 0.0f) || !std::isfinite(M)) return false;
    PhaseProperties const &p = phaseForSubstance(id);
    if (!p.valid) return false;
    if (!(p.boilingPointK > 1.0f) || !std::isfinite(p.boilingPointK)) return false;
    if (!(p.latentHeatVaporization > 0.0f) || !std::isfinite(p.latentHeatVaporization)) return false;
    if (!(p.referencePressurePa > 0.0f) || !std::isfinite(p.referencePressurePa)) return false;
    if (!(saturationVaporPressurePa(id, p.boilingPointK) > 0.0)) return false;
    if (!(saturationTemperatureK(id, p.referencePressurePa) > 1.0)) return false;
    if (!(referenceGasDensityKgM3(id) > 0.0)) return false;
    if (!(liquidFillToMassKg(id, 1.0, 4.0) > 0.0)) return false;
    if (!(massKgToGasAmount(id, liquidFillToMassKg(id, 1.0, 4.0), 4.0) > 0.0)) return false;
    return true;
}

float liquidGasTransferRateFillPerSec() {
    return kMaxLiquidGasFillPerSec;
}

void applyLiquidGasStatsToWaterCompat(LiquidGasPhaseTickStats const &lg, WaterPhaseTickStats &st) {
    st.massBoiledKg += lg.massVaporizedKg;
    st.massCondensedKg += lg.massCondensedKg;
    st.latentAbsorbedJ += lg.latentVaporizationAbsorbedJ;
    st.latentReleasedJ += lg.latentCondensationReleasedJ;
    st.vaporPlacedAmount += lg.vaporPlacedAmount;
    st.liquidFillRemoved += lg.liquidFillRemoved;
    st.liquidFillAdded += lg.liquidFillAdded;
    st.boilCells += lg.vaporizedCells;
    st.condenseCells += lg.condensedCells;
    st.blockedBoilNoGasSpace += lg.blockedNoGasSpace;
    st.blockedCondenseNoLiquidSpace += lg.blockedNoLiquidSpace;
    st.blockedBoilEquilibrium += lg.blockedVaporizationEquilibrium;
    st.blockedCondenseEquilibrium += lg.blockedCondensationEquilibrium;
    st.blockedBoilEnergy += lg.blockedVaporizationEnergy;
    st.blockedCondenseEnergy += lg.blockedCondensationEnergy;
    st.blockedBoilSafetyLimit += lg.blockedSafetyLimit;
    st.blockedBoil += lg.blockedNoGasSpace + lg.blockedVaporizationEquilibrium
        + lg.blockedVaporizationEnergy + lg.blockedSafetyLimit;
    st.blockedCondense += lg.blockedNoLiquidSpace + lg.blockedCondensationEquilibrium
        + lg.blockedCondensationEnergy;
}

LiquidGasPhaseTickStats stepLiquidGasPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt, uint32_t salt)
{
    LiquidGasPhaseTickStats st;
    float dtSafe = (dt > 1.0e-6f && std::isfinite(dt)) ? dt : PHYSICS_DT;
    double cpm = fluid.config.cellsPerMeter;
    float pref = gas.config.referencePressurePa;
    if (!(pref > 0.0f) || !std::isfinite(pref)) pref = GAS_REFERENCE_PRESSURE_PA;
    double mRateFillRef = liquidFillToMassKg(SUBSTANCE_WATER,
        static_cast<double>(kMaxLiquidGasFillPerSec) * static_cast<double>(dtSafe), cpm);

    bool anyLiquid = fluid.expectedVolume > 1.0e-10;
    // Compatibility/accounting only: Water vapor counters are a cheap skip while
    // Water is the only live Liquid/Gas substance. Identity still comes from
    // stored components, not these counters.
    bool anyVapor = gas.currentWaterVapor > GAS_MIN_AMOUNT
        || gas.expectedWaterVapor > GAS_MIN_AMOUNT;
    if (!anyLiquid && !anyVapor) return st;

    if (anyLiquid) {
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int i = FluidEngine::ci(x, y);
            if (isBlockedSolid(fluid, rigid, x, y)) continue;
            float fill = fluid.fill[static_cast<size_t>(i)];
            if (fill < kMinFillMove) continue;
            SubstanceId id = SUBSTANCE_NONE;
            float liquidFill = 0.0f;
            if (!resolvePureLiquidForVaporization(fluid, i, id, liquidFill)) continue;
            PhaseProperties const &phase = phaseForSubstance(id);
            float Lv = phase.latentHeatVaporization;
            float cpLiquid = thermalForSubstance(id).specificHeat;
            float cpVapor = gasPhaseSpecificHeat(id);
            float C = ThermalEngine::liquidCapacity(fluid, i);
            float E = fluid.liquidHeat[static_cast<size_t>(i)];
            float T = (C > MIN_THERMAL_CAPACITY) ? tempFromEnergy(E, C)
                : ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
            float Tplace = T;
            // Tpay is at least Tm+1 when Solid is supported, so colder liquid cannot
            // pay latent heat. Skip dest search (no world mutation) in that case.
            if (supportsPhase(id, MatterPhase::Solid)) {
                float Tm = phase.meltingPointK;
                if (Tm > 1.0f && std::isfinite(Tm) && !(T > Tm + 1.0f + kBoilEpsK))
                    continue;
            }
            double PsatT = saturationVaporPressurePa(id, T);
            float maxAtm = static_cast<float>(PsatT / static_cast<double>(pref));
            if (!(maxAtm > 0.0f) || !std::isfinite(maxAtm)) maxAtm = kSolverSafetyAtm;
            maxAtm = std::min(maxAtm, kSolverSafetyAtm);
            double mRateFill = liquidFillToMassKg(id,
                static_cast<double>(kMaxLiquidGasFillPerSec) * static_cast<double>(dtSafe), cpm);
            if (!(mRateFill > 0.0)) mRateFill = mRateFillRef;
            double mFill = liquidFillToMassKg(id, liquidFill, cpm);
            double vaporWantProbe = massKgToGasAmount(id, std::min(mFill, mRateFill), cpm);
            bool hitSafety = false, hitEq = false;
            int dest = findVaporDest(fluid, rigid, gas, id, x, y,
                static_cast<float>(std::max(vaporWantProbe, 1.0e-4)),
                salt, maxAtm, &hitSafety, &hitEq);
            if (dest < 0) {
                if (hitSafety) ++st.blockedSafetyLimit;
                else if (hitEq) ++st.blockedVaporizationEquilibrium;
                else ++st.blockedNoGasSpace;
                continue;
            }
            if (!gasCellCanAcceptComponent(gas, dest, id)) {
                ++st.blockedNoGasSpace;
                continue;
            }
            int dx = dest % GW, dy = dest / GW;
            float vol = gas.availableVolume(fluid, dx, dy);
            float Pdest = std::max(1.0f, gasCellPressureAtm(gas, dest, vol) * pref);
            float Tsat = static_cast<float>(saturationTemperatureK(id, Pdest));
            float Tpay = Tsat;
            if (supportsPhase(id, MatterPhase::Solid)) {
                float Tm = phase.meltingPointK;
                if (Tm > 1.0f && std::isfinite(Tm))
                    Tpay = std::max(Tsat, Tm + 1.0f);
            }
            if (!(T > Tpay + kBoilEpsK)) {
                ++st.blockedVaporizationEquilibrium;
                continue;
            }
            wakeAllEngines(fluid, gas, thermal, x, y);
            float Eplat = (C > MIN_THERMAL_CAPACITY) ? energyFromTemp(C, Tpay) : 0.0f;
            float excess = (C > MIN_THERMAL_CAPACITY) ? std::max(0.0f, E - Eplat) : 0.0f;
            float neighborQ = heatAvailableAround(fluid, rigid, gas, x, y, Tpay, false);
            float totalQ = excess + neighborQ;
            double denom = static_cast<double>(Lv)
                + static_cast<double>(cpVapor) * static_cast<double>(Tplace)
                - static_cast<double>(cpLiquid) * static_cast<double>(Tpay);
            if (denom < 1.0) denom = static_cast<double>(Lv);
            double mEnergy = (denom > 1.0) ? static_cast<double>(totalQ) / denom : mFill;
            double mWant = std::min({mEnergy, mFill, mRateFill});
            if (!(mWant > 1.0e-9)) {
                ++st.blockedVaporizationEnergy;
                continue;
            }
            double vaporWant = massKgToGasAmount(id, mWant, cpm);
            float roomThermo = gasAmountRoomForPressureAtm(gas, dest, vol, maxAtm);
            float roomSafety = gasAmountRoomForPressureAtm(gas, dest, vol, kSolverSafetyAtm);
            float room = std::min(roomThermo, roomSafety);
            if (room <= GAS_MIN_AMOUNT) {
                if (roomSafety <= GAS_MIN_AMOUNT) ++st.blockedSafetyLimit;
                else ++st.blockedVaporizationEquilibrium;
                continue;
            }
            double vaporGot = std::min(vaporWant, static_cast<double>(room));
            double m = gasAmountToMassKg(id, vaporGot, cpm);
            if (!(m > 1.0e-9)) {
                ++st.blockedVaporizationEquilibrium;
                continue;
            }
            float dFill = static_cast<float>(massKgToLiquidFill(id, m, cpm));
            dFill = std::min(dFill, liquidFill);
            if (dFill < kMinFillMove) continue;
            m = liquidFillToMassKg(id, dFill, cpm);
            vaporGot = massKgToGasAmount(id, m, cpm);
            if (!gasCellCanAcceptComponent(gas, dest, id)) {
                ++st.blockedNoGasSpace;
                continue;
            }

            float Eremoved = 0.0f;
            if (fill > 1.0e-8f) Eremoved = E * (dFill / fill);
            float Edest = static_cast<float>(m) * cpVapor * Tplace;
            float Qneed = static_cast<float>(m * static_cast<double>(Lv)) + Edest - Eremoved;

            (void)fluid.takeLiquidVolume(i, dFill);
            float remainFill = fluid.fill[static_cast<size_t>(i)];
            float remainE = std::max(0.0f, E - Eremoved);
            if (remainFill <= 1.0e-8f) {
                fluid.liquidHeat[static_cast<size_t>(i)] = 0.0f;
                if (remainE > 0.0f) Qneed -= remainE;
            } else {
                fluid.liquidHeat[static_cast<size_t>(i)] = remainE;
            }
            if (Qneed > 0.0f) {
                float pay = Qneed;
                float newC = ThermalEngine::liquidCapacity(fluid, i);
                takeHeatAboveK(fluid.liquidHeat[static_cast<size_t>(i)], newC, Tpay, pay);
                if (pay > 0.0f)
                    takeHeatAround(fluid, rigid, gas, thermal, x, y, Tpay, pay, false);
                if (pay > 0.0f)
                    takeHeatAround(fluid, rigid, gas, thermal, dx, dy, Tpay, pay, true);
                if (pay > 0.0f)
                    takeHeatAboveK(gas.heat[static_cast<size_t>(dest)],
                        ThermalEngine::gasCapacity(gas, dest), Tpay, pay);
            } else if (Qneed < 0.0f) {
                gas.heat[static_cast<size_t>(dest)] += -Qneed;
            }

            gas.addGasComponentAmount(dest, id, static_cast<float>(vaporGot));
            gas.heat[static_cast<size_t>(dest)] += Edest;
            if (!std::isfinite(gas.heat[static_cast<size_t>(dest)]) || gas.heat[static_cast<size_t>(dest)] < 0.0f)
                gas.heat[static_cast<size_t>(dest)] = 0.0f;

            fluid.expectedVolume -= static_cast<double>(dFill);
            gas.expectedAmount += vaporGot;
            // Compatibility/accounting only: GasEngine Water vapor ledger.
            if (id == SUBSTANCE_WATER)
                gas.expectedWaterVapor += vaporGot;
            st.massVaporizedKg += m;
            st.latentVaporizationAbsorbedJ += m * static_cast<double>(Lv);
            st.vaporPlacedAmount += vaporGot;
            st.liquidFillRemoved += dFill;
            ++st.vaporizedCells;
            wakeAllEngines(fluid, gas, thermal, x, y);
            wakeAllEngines(fluid, gas, thermal, dx, dy);
        }
    }

    if (anyVapor) {
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int gi = GasEngine::ci(x, y);
            if (isBlockedSolid(fluid, rigid, x, y)) continue;
            float vol = gas.availableVolume(fluid, x, y);
            if (vol < GAS_MIN_VOLUME) continue;
            GasComponentView view = gas.gasComponents(gi);
            float Tgas = ThermalEngine::gasTempK(gas, gi);
            for (int c = 0; c < view.count && c < kMaxGasComponents; ++c) {
                SubstanceId id = view.items[c].id;
                float vap = view.items[c].amount;
                if (!(vap > GAS_MIN_AMOUNT) || !supportsLiveLiquidGasTransition(id)) continue;
                PhaseProperties const &phase = phaseForSubstance(id);
                float Lv = phase.latentHeatVaporization;
                float cpLiquid = thermalForSubstance(id).specificHeat;
                double Pv = static_cast<double>(gas.gasPartialPressurePaFromCurrentState(gi, id));
                double Psat = saturationVaporPressurePa(id, Tgas);
                if (!(Pv > Psat * (1.0 + static_cast<double>(kSatHystRel)))) {
                    if (Pv > Psat * 0.9 && vap > 1.0e-3f)
                        ++st.blockedCondensationEquilibrium;
                    continue;
                }
                float psatAtm = static_cast<float>(Psat / static_cast<double>(pref));
                float satAmt = gasAmountFromPressureAtm(psatAtm, vol, Tgas);
                double excessAmt = static_cast<double>(vap) - static_cast<double>(satAmt);
                if (!(excessAmt > GAS_MIN_AMOUNT) || !std::isfinite(excessAmt)) continue;
                excessAmt = std::min(excessAmt, static_cast<double>(vap));
                double mRateFill = liquidFillToMassKg(id,
                    static_cast<double>(kMaxLiquidGasFillPerSec) * static_cast<double>(dtSafe), cpm);
                if (!(mRateFill > 0.0)) mRateFill = mRateFillRef;
                double mAvail = gasAmountToMassKg(id, vap, cpm);
                double mExcess = gasAmountToMassKg(id, excessAmt, cpm);
                double mWant = std::min({mAvail, mExcess, mRateFill});
                if (!(mWant > 1.0e-9)) continue;
                int dest = findLiquidDest(fluid, rigid, gas, id, x, y, salt, true, Tgas);
                if (dest < 0) {
                    ++st.blockedNoLiquidSpace;
                    continue;
                }
                int dx = dest % GW, dy = dest / GW;
                if (fluid.solid[static_cast<size_t>(dest)]) {
                    ++st.blockedNoLiquidSpace;
                    continue;
                }
                if (!liquidCellAcceptsCondensate(fluid, dest, id)) {
                    ++st.blockedNoLiquidSpace;
                    continue;
                }
                float roomFill = 1.0f - fluid.fill[static_cast<size_t>(dest)];
                if (roomFill < kMinFillMove) {
                    ++st.blockedNoLiquidSpace;
                    continue;
                }
                float Pdest = localAmbientPressurePa(fluid, rigid, gas, dx, dy);
                float TsatDest = static_cast<float>(saturationTemperatureK(id, Pdest));
                float Cdest = ThermalEngine::liquidCapacity(fluid, dest);
                float Edest0 = fluid.liquidHeat[static_cast<size_t>(dest)];
                float Tdest = (Cdest > MIN_THERMAL_CAPACITY) ? tempFromEnergy(Edest0, Cdest) : Tgas;
                if (Tdest > TsatDest - kBoilEpsK) {
                    ++st.blockedCondensationEquilibrium;
                    continue;
                }
                double m = mWant;
                float dFill = static_cast<float>(massKgToLiquidFill(id, m, cpm));
                dFill = std::min(dFill, roomFill);
                if (dFill < kMinFillMove) continue;
                m = liquidFillToMassKg(id, dFill, cpm);
                float vaporTake = static_cast<float>(massKgToGasAmount(id, m, cpm));
                float heatBefore = gas.heat[static_cast<size_t>(gi)];
                float taken = gas.takeGasComponentAmount(gi, id, vaporTake);
                if (taken <= GAS_MIN_AMOUNT) {
                    ++st.blockedCondensationEnergy;
                    continue;
                }
                float Eremoved = std::max(0.0f, heatBefore - gas.heat[static_cast<size_t>(gi)]);
                m = gasAmountToMassKg(id, taken, cpm);
                dFill = static_cast<float>(massKgToLiquidFill(id, m, cpm));
                float TmPlace = 0.0f;
                if (supportsPhase(id, MatterPhase::Solid) && phase.meltingPointK > 1.0f)
                    TmPlace = phase.meltingPointK;
                float Tplace = std::min(Tgas, TsatDest);
                if (TmPlace > 1.0f) Tplace = std::max(Tplace, TmPlace);
                Tplace = std::clamp(Tplace, MIN_SAFE_TEMPERATURE_K,
                    (TsatDest > TmPlace) ? TsatDest : std::max(TsatDest, TmPlace));
                float Eliquid = static_cast<float>(m) * cpLiquid * Tplace;
                float leftover = static_cast<float>(m * static_cast<double>(Lv)) + Eremoved - Eliquid;
                if (id == SUBSTANCE_WATER)
                    fluid.addLiquidFill(dest, dFill, Eliquid);
                else {
                    fluid.addLiquidComponentAmount(dest, id, dFill);
                    fluid.liquidHeat[static_cast<size_t>(dest)] += Eliquid;
                }
                if (leftover > 0.0f) {
                    float TliqCap = std::max(TmPlace, TsatDest - 1.0f);
                    addHeatToK(fluid.liquidHeat[static_cast<size_t>(dest)],
                        ThermalEngine::liquidCapacity(fluid, dest), TliqCap, leftover);
                    addHeatToK(gas.heat[static_cast<size_t>(gi)], ThermalEngine::gasCapacity(gas, gi),
                        TsatDest + 40.0f, leftover);
                    for (int n = 0; n < 4; ++n) {
                        if (!(leftover > 0.0f)) break;
                        dumpHeatBudget(fluid, rigid, gas, thermal, dx + kDx[n], dy + kDy[n],
                            TsatDest + 80.0f, leftover);
                    }
                    for (int n = 0; n < 4 && leftover > 0.0f; ++n) {
                        int nx = dx + kDx[n], ny = dy + kDy[n];
                        if (!FluidEngine::inside(nx, ny)) continue;
                        int ni = FluidEngine::ci(nx, ny);
                        if (!fluid.solid[static_cast<size_t>(ni)]) continue;
                        fluid.solidHeat[static_cast<size_t>(ni)] += leftover;
                        leftover = 0.0f;
                        thermal.wakeCell(nx, ny);
                    }
                    if (leftover > 0.0f)
                        gas.heat[static_cast<size_t>(gi)] += leftover;
                } else if (leftover < 0.0f) {
                    float need = -leftover;
                    takeHeatAround(fluid, rigid, gas, thermal, dx, dy, Tplace, need, true);
                }
                int xw = dest % GW, yw = dest / GW;
                if (id == SUBSTANCE_WATER)
                    fluid.waterShade[static_cast<size_t>(dest)] = fluid.makeShade(xw, yw);
                fluid.expectedVolume += static_cast<double>(dFill);
                gas.expectedAmount -= taken;
                if (id == SUBSTANCE_WATER)
                    gas.expectedWaterVapor -= taken;
                st.massCondensedKg += m;
                st.latentCondensationReleasedJ += m * static_cast<double>(Lv);
                st.liquidFillAdded += dFill;
                ++st.condensedCells;
                wakeAllEngines(fluid, gas, thermal, x, y);
                wakeAllEngines(fluid, gas, thermal, dx, dy);
                view = gas.gasComponents(gi);
            }
        }
    }
    return st;
}

PhaseChangeTickStats stepPhaseChanges(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt)
{
    PhaseChangeTickStats st;
    WaterPhaseTickStats water = stepWaterPhaseChange(fluid, rigid, gas, thermal, dt);
    st.liquidGas.massVaporizedKg = water.massBoiledKg;
    st.liquidGas.massCondensedKg = water.massCondensedKg;
    st.liquidGas.latentVaporizationAbsorbedJ = water.latentAbsorbedJ;
    st.liquidGas.latentCondensationReleasedJ = water.latentReleasedJ;
    st.liquidGas.vaporPlacedAmount = water.vaporPlacedAmount;
    st.liquidGas.liquidFillRemoved = water.liquidFillRemoved;
    st.liquidGas.liquidFillAdded = water.liquidFillAdded;
    st.liquidGas.vaporizedCells = water.boilCells;
    st.liquidGas.condensedCells = water.condenseCells;
    st.liquidGas.blockedNoGasSpace = water.blockedBoilNoGasSpace;
    st.liquidGas.blockedNoLiquidSpace = water.blockedCondenseNoLiquidSpace;
    st.liquidGas.blockedVaporizationEquilibrium = water.blockedBoilEquilibrium;
    st.liquidGas.blockedCondensationEquilibrium = water.blockedCondenseEquilibrium;
    st.liquidGas.blockedVaporizationEnergy = water.blockedBoilEnergy;
    st.liquidGas.blockedCondensationEnergy = water.blockedCondenseEnergy;
    st.liquidGas.blockedSafetyLimit = water.blockedBoilSafetyLimit;
    return st;
}
