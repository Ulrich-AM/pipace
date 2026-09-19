#include "world/PhaseChangeEngine.h"

#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/LiquidMixtureProperties.h"
#include "substance/PhaseTransfer.h"
#include "substance/SubstanceRegistry.h"
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
constexpr int kSolidSearchLimit = 32;
constexpr float kFreezeEpsK = 0.05f;
constexpr float kMaxSolidPixelsPerSec = 2.5f;
constexpr float kMinSolidRemain = 0.03f;

struct PhaseScratch {
    std::vector<int> q;
    std::vector<uint32_t> seen;
    uint32_t gen = 1;
    std::vector<uint8_t> reserved;
    std::vector<MaterialId> reservedMat;
    std::vector<int> spawnCells;
    std::vector<int> worldCells;
    std::vector<int> spawnByMat[MATERIAL_COUNT];
    struct SpawnJob {
        int src = -1;
        int dest = -1;
        SubstanceId id = SUBSTANCE_NONE;
        MaterialId mat = MATERIAL_EMPTY;
        float massKg = 0.0f;
    };
    std::vector<SpawnJob> spawnJobs;
    std::vector<RigidBodyEngine::SpawnedSourcePixel> spawnUnique;

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

PhaseScratch &solidLiquidScratch() {
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

bool resolvePureLiquidForSolidification(FluidEngine const &fluid, int i,
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
    if (!supportsLiveSolidLiquidTransition(found)) return false;
    outId = found;
    outFill = fluid.liquidComponentAmount(i, found);
    return outFill >= kMinFillMove;
}

void addLiquidComponentWithHeat(FluidEngine &fluid, int dest, SubstanceId id, float dFill, float dHeat) {
    if (dest < 0 || dest >= GW * GH || !(dFill > 0.0f)) return;
    fluid.addLiquidComponentAmount(dest, id, dFill);
    fluid.liquidHeat[static_cast<size_t>(dest)] += dHeat;
    if (id == SUBSTANCE_WATER) {
        int x = dest % GW, y = dest / GW;
        fluid.waterShade[static_cast<size_t>(dest)] = fluid.makeShade(x, y);
    }
    fluid.clearEmptyLiquidCell(dest);
}

float pendingSolidCapacityJK(SubstanceId id, float kg) {
    if (!(kg > 0.0f) || !validSubstance(id)) return 0.0f;
    float cp = solidPhaseSpecificHeat(id);
    if (!(cp > 0.0f) || !std::isfinite(cp)) return 0.0f;
    return kg * cp;
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

float heatRoomAround(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int x, int y, float tMax)
{
    float room = 0.0f;
    auto add = [&](int cx, int cy) {
        if (!FluidEngine::inside(cx, cy)) return;
        int i = FluidEngine::ci(cx, cy);
        int body = rigid.occupant[static_cast<size_t>(i)];
        if (body >= 0 && body < static_cast<int>(rigid.bodies.size())) {
            RigidBody const &b = rigid.bodies[static_cast<size_t>(body)];
            float lx, ly;
            RigidBodyEngine::worldToLocal(b, cx + 0.5f, cy + 0.5f, lx, ly);
            int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
            if (RigidBodyEngine::maskOccupied(b, ix, iy)) {
                int li = iy * b.maskW + ix;
                if (li >= 0 && li < static_cast<int>(b.heat.size()))
                    room += heatRoomToK(b.heat[static_cast<size_t>(li)],
                        ThermalEngine::rigidPixelCapacity(b, li), tMax);
            }
            return;
        }
        if (fluid.solid[static_cast<size_t>(i)]) {
            room += heatRoomToK(fluid.solidHeat[static_cast<size_t>(i)],
                ThermalEngine::wallCapacity(fluid, i), tMax);
            return;
        }
        if (fluid.fill[static_cast<size_t>(i)] > 1.0e-8f) {
            room += heatRoomToK(fluid.liquidHeat[static_cast<size_t>(i)],
                ThermalEngine::liquidCapacity(fluid, i), tMax);
            return;
        }
        if (gas.amount[static_cast<size_t>(i)] > GAS_MIN_AMOUNT)
            room += heatRoomToK(gas.heat[static_cast<size_t>(i)],
                ThermalEngine::gasCapacity(gas, i), tMax);
    };
    add(x, y);
    for (int n = 0; n < 4; ++n) add(x + kDx[n], y + kDy[n]);
    return room;
}

bool solidDestCellFree(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    std::vector<uint8_t> const &reserved, int x, int y)
{
    if (!FluidEngine::inside(x, y)) return false;
    int i = FluidEngine::ci(x, y);
    if (fluid.solid[static_cast<size_t>(i)]) return false;
    if (rigid.occupant[static_cast<size_t>(i)] >= 0) return false;
    if (!reserved.empty() && reserved[static_cast<size_t>(i)]) return false;
    return true;
}

bool adjacentSameMaterial(RigidBodyEngine const &rigid, int x, int y, MaterialId mat) {
    if (mat == MATERIAL_EMPTY) return false;
    for (int n = 0; n < 4; ++n) {
        int nx = x + kDx[n], ny = y + kDy[n];
        if (!FluidEngine::inside(nx, ny)) continue;
        if (rigid.worldCellMaterial(nx, ny) == mat) return true;
    }
    return false;
}

bool adjacentReservedMaterial(std::vector<uint8_t> const &reserved,
    std::vector<MaterialId> const &reservedMat, int x, int y, MaterialId mat)
{
    if (mat == MATERIAL_EMPTY || reserved.empty()) return false;
    for (int n = 0; n < 4; ++n) {
        int nx = x + kDx[n], ny = y + kDy[n];
        if (!FluidEngine::inside(nx, ny)) continue;
        int ni = FluidEngine::ci(nx, ny);
        if (!reserved[static_cast<size_t>(ni)]) continue;
        if (!reservedMat.empty() && reservedMat[static_cast<size_t>(ni)] == mat)
            return true;
    }
    return false;
}

int findSolidDest(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    std::vector<uint8_t> const &reserved, std::vector<MaterialId> const &reservedMat,
    MaterialId mat, int sx, int sy)
{
    auto consider = [&](int x, int y) -> bool {
        return solidDestCellFree(fluid, rigid, reserved, x, y);
    };
    int best = -1;
    float bestScore = -1.0f;
    uint32_t bestKey = ~0u;
    auto take = [&](int x, int y, float bonus) {
        if (!consider(x, y)) return;
        float s = bonus;
        if (adjacentSameMaterial(rigid, x, y, mat)) s += 1000.0f;
        if (adjacentReservedMaterial(reserved, reservedMat, x, y, mat)) s += 500.0f;
        if (x == sx && y == sy) s += 100.0f;
        else {
            for (int n = 0; n < 4; ++n) {
                if (x == sx + kDx[n] && y == sy + kDy[n]) {
                    s += 50.0f;
                    break;
                }
            }
        }
        s -= static_cast<float>(std::abs(x - sx) + std::abs(y - sy)) * 2.0f;
        uint32_t key = destTieKey(x, y, 1u);
        if (destBetter(s, key, bestScore, bestKey)) {
            bestScore = s;
            bestKey = key;
            best = FluidEngine::ci(x, y);
        }
    };

    // Priority 1: free cells adjacent to existing same-material rigid solids.
    if (mat != MATERIAL_EMPTY) {
        for (RigidBody const &b : rigid.bodies) {
            for (int li : b.occupiedLocal) {
                if (li < 0 || li >= static_cast<int>(b.mask.size())) continue;
                if (b.mask[static_cast<size_t>(li)] != mat) continue;
                int lx = li % b.maskW, ly = li / b.maskW;
                float wx, wy;
                RigidBodyEngine::localToWorld(b, static_cast<float>(lx) + 0.5f,
                    static_cast<float>(ly) + 0.5f, wx, wy);
                int gx = static_cast<int>(std::floor(wx));
                int gy = static_cast<int>(std::floor(wy));
                for (int n = 0; n < 4; ++n) {
                    int nx = gx + kDx[n], ny = gy + kDy[n];
                    if (std::abs(nx - sx) + std::abs(ny - sy) > kSolidSearchLimit) continue;
                    take(nx, ny, 2000.0f);
                }
            }
        }
        if (best >= 0) return best;
    }

    take(sx, sy, 0.0f);
    for (int n = 0; n < 4; ++n)
        take(sx + kDx[n], sy + kDy[n], 0.0f);
    PhaseScratch &sc = solidLiquidScratch();
    sc.prepareSeen();
    sc.q.reserve(kSolidSearchLimit);
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
    for (size_t head = 0; head < sc.q.size() && static_cast<int>(head) < kSolidSearchLimit; ++head) {
        int i = sc.q[head];
        int x = i % GW, y = i / GW;
        take(x, y, 0.0f);
        for (int n = 0; n < 4; ++n) push(x + kDx[n], y + kDy[n]);
    }
    if (best < 0 && consider(sx, sy))
        best = FluidEngine::ci(sx, sy);
    return best;
}

bool anyEligibleGasVapor(GasEngine const &gas) {
    for (int i = 0; i < GW * GH; ++i) {
        if (gas.amount[static_cast<size_t>(i)] <= GAS_MIN_AMOUNT
            && gas.gasCompCount[static_cast<size_t>(i)] == 0)
            continue;
        if (gasViewHasLiveLiquidGasVapor(gas.gasComponents(i)))
            return true;
    }
    return false;
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

float solidLiquidTransferRatePixelsPerSec() {
    return kMaxSolidPixelsPerSec;
}

bool gasViewHasLiveLiquidGasVapor(GasComponentView const &view) {
    int n = view.count;
    if (n > kMaxGasComponents) n = kMaxGasComponents;
    for (int i = 0; i < n; ++i) {
        SubstanceId id = view.items[i].id;
        float amt = view.items[i].amount;
        if (!(amt > GAS_MIN_AMOUNT) || !std::isfinite(amt)) continue;
        if (supportsLiveLiquidGasTransition(id)) return true;
    }
    return false;
}

bool supportsLiveSolidLiquidTransition(SubstanceId id) {
    if (id == SUBSTANCE_NONE || !validSubstance(id)) return false;
    if (!supportsPhase(id, MatterPhase::Solid) || !supportsPhase(id, MatterPhase::Liquid))
        return false;
    if (!hasMechanicalProperties(id) || !hasFluidProperties(id)) return false;
    if (!canTransition(id, MatterPhase::Solid, MatterPhase::Liquid)) return false;
    if (rigidMaterialForSubstance(id) == MATERIAL_EMPTY) return false;
    PhaseProperties const &p = phaseForSubstance(id);
    if (!p.valid) return false;
    if (!(p.meltingPointK > 1.0f) || !std::isfinite(p.meltingPointK)) return false;
    if (!(p.latentHeatFusion > 0.0f) || !std::isfinite(p.latentHeatFusion)) return false;
    if (!(solidFractionToMassKg(id, 1.0, 4.0) > 0.0)) return false;
    if (!(liquidFillToMassKg(id, 1.0, 4.0) > 0.0)) return false;
    return true;
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

void applySolidLiquidStatsToWaterCompat(SolidLiquidPhaseTickStats const &sl, WaterPhaseTickStats &st) {
    st.massFrozenKg += sl.massSolidifiedKg;
    st.massMeltedKg += sl.massMeltedKg;
    st.latentFusionReleasedJ += sl.latentFusionReleasedJ;
    st.latentFusionAbsorbedJ += sl.latentFusionAbsorbedJ;
    st.liquidFillRemoved += sl.liquidFillRemoved;
    st.liquidFillAdded += sl.liquidFillAdded;
    st.freezeCells += sl.solidifyCells;
    st.meltPixels += sl.meltPixels;
    st.icePixelsSpawned += sl.solidPixelsSpawned;
    st.blockedFreeze += sl.blockedSolidify;
    st.blockedMelt += sl.blockedMelt;
    st.blockedFreezeNoDest += sl.blockedSolidifyNoDest;
    st.blockedMeltNoDest += sl.blockedMeltNoDest;
    st.blockedMeltEnergy += sl.blockedMeltEnergy;
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
    bool anyVapor = anyEligibleGasVapor(gas);
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

SolidLiquidPhaseTickStats stepSolidLiquidPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt, uint32_t salt)
{
    SolidLiquidPhaseTickStats st;
    float dtSafe = (dt > 1.0e-6f && std::isfinite(dt)) ? dt : PHYSICS_DT;
    double cpm = fluid.config.cellsPerMeter;

    PhaseScratch &sc = solidLiquidScratch();
    size_t nCells = static_cast<size_t>(GW * GH);
    if (sc.reserved.size() != nCells) {
        sc.reserved.assign(nCells, 0);
        sc.reservedMat.assign(nCells, MATERIAL_EMPTY);
    } else {
        std::fill(sc.reserved.begin(), sc.reserved.end(), 0);
        std::fill(sc.reservedMat.begin(), sc.reservedMat.end(), MATERIAL_EMPTY);
    }
    for (int m = 0; m < MATERIAL_COUNT; ++m)
        sc.spawnByMat[m].clear();
    sc.spawnCells.clear();
    sc.spawnJobs.clear();
    sc.spawnUnique.clear();

    auto meltPendingAt = [&](int x, int y) {
        int i = FluidEngine::ci(x, y);
        SubstanceId id = fluid.solidifyPendingSubstance(i);
        if (id == SUBSTANCE_NONE || !supportsLiveSolidLiquidTransition(id)) return;
        PhaseProperties const &phase = phaseForSubstance(id);
        float Tm = phase.meltingPointK;
        float Lf = phase.latentHeatFusion;
        double pixelKg = solidFractionToMassKg(id, 1.0, cpm);
        if (!(pixelKg > 1.0e-9) || !(Lf > 1.0f) || !(Tm > 1.0f)) return;
        float pendingKg = fluid.solidifyPendingMassKg(i);
        if (!(pendingKg > 1.0e-9f)) return;
        float cap = pendingSolidCapacityJK(id, pendingKg);
        float Epend = fluid.solidifyPendingSensibleJ(i);
        float Tpend = (cap > MIN_THERMAL_CAPACITY) ? tempFromEnergy(Epend, cap)
            : ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
        if (!(Tpend + kFreezeEpsK >= Tm)) return;
        float Eplat = energyFromTemp(std::max(cap, MIN_THERMAL_CAPACITY), Tm);
        float excess = std::max(0.0f, Epend - Eplat);
        float neighborQ = heatAvailableAround(fluid, rigid, gas, x, y, Tm, false);
        float totalQ = excess + neighborQ;
        if (!(totalQ > 0.0f)) {
            ++st.blockedMelt;
            ++st.blockedMeltEnergy;
            return;
        }
        double mRate = pixelKg * static_cast<double>(kMaxSolidPixelsPerSec) * static_cast<double>(dtSafe);
        double mEnergy = static_cast<double>(totalQ) / static_cast<double>(Lf);
        double mWant = std::min({mEnergy, static_cast<double>(pendingKg), mRate});
        if (!(mWant > 1.0e-9)) {
            ++st.blockedMelt;
            ++st.blockedMeltEnergy;
            return;
        }
        int dest = -1;
        if (!isBlockedSolid(fluid, rigid, x, y) && liquidCellAcceptsCondensate(fluid, i, id)
            && (1.0f - fluid.fill[static_cast<size_t>(i)]) >= kMinFillMove)
            dest = i;
        if (dest < 0)
            dest = findLiquidDest(fluid, rigid, gas, id, x, y, salt, false, Tm);
        if (dest < 0) {
            ++st.blockedMelt;
            ++st.blockedMeltNoDest;
            return;
        }
        int dx = dest % GW, dy = dest / GW;
        if (isBlockedSolid(fluid, rigid, dx, dy) || !liquidCellAcceptsCondensate(fluid, dest, id)) {
            ++st.blockedMelt;
            ++st.blockedMeltNoDest;
            return;
        }
        float roomFill = 1.0f - fluid.fill[static_cast<size_t>(dest)];
        if (roomFill < kMinFillMove) {
            ++st.blockedMelt;
            ++st.blockedMeltNoDest;
            return;
        }
        float dFill = static_cast<float>(massKgToLiquidFill(id, mWant, cpm));
        dFill = std::min(dFill, roomFill);
        if (dFill < kMinFillMove) {
            ++st.blockedMelt;
            ++st.blockedMeltNoDest;
            return;
        }
        double m = liquidFillToMassKg(id, dFill, cpm);
        float takeHeat = 0.0f;
        if (!fluid.takeSolidifyPendingKg(i, id, static_cast<float>(m), &takeHeat)) {
            ++st.blockedMelt;
            return;
        }
        float pay = static_cast<float>(m * static_cast<double>(Lf));
        float iceTake = std::min(pay, excess);
        pay -= iceTake;
        if (pay > 0.0f)
            takeHeatAround(fluid, rigid, gas, thermal, x, y, Tm, pay, false);
        float leftoverEx = std::max(0.0f, excess - iceTake);
        float remainPend = fluid.solidifyPendingMassKg(i);
        if (remainPend > 1.0e-12f) {
            float remainCap = pendingSolidCapacityJK(id, remainPend);
            fluid.solidifyPendingHeatJ[static_cast<size_t>(i)] =
                energyFromTemp(std::max(remainCap, MIN_THERMAL_CAPACITY), Tm) + leftoverEx;
        }
        float sensible = static_cast<float>(m) * thermalForSubstance(id).specificHeat * Tm;
        addLiquidComponentWithHeat(fluid, dest, id, dFill, sensible);
        fluid.expectedVolume += static_cast<double>(dFill);
        st.massMeltedKg += m;
        st.latentFusionAbsorbedJ += m * static_cast<double>(Lf);
        st.liquidFillAdded += dFill;
        ++st.meltPixels;
        wakeAllEngines(fluid, gas, thermal, x, y);
        wakeAllEngines(fluid, gas, thermal, dx, dy);
    };

    auto trySpawnPending = [&](int x, int y) {
        int i = FluidEngine::ci(x, y);
        SubstanceId id = fluid.solidifyPendingSubstance(i);
        if (id == SUBSTANCE_NONE || !supportsLiveSolidLiquidTransition(id)) return;
        MaterialId mat = rigidMaterialForSubstance(id);
        if (mat == MATERIAL_EMPTY || mat >= MATERIAL_COUNT) {
            ++st.blockedSolidify;
            ++st.blockedSolidifyNoDest;
            return;
        }
        double pixelKg = solidFractionToMassKg(id, 1.0, cpm);
        if (!(pixelKg > 1.0e-9)) return;
        double avail = static_cast<double>(fluid.solidifyPendingMassKg(i));
        while (avail + 1.0e-9 >= pixelKg) {
            int dest = findSolidDest(fluid, rigid, sc.reserved, sc.reservedMat, mat, x, y);
            if (dest < 0) {
                ++st.blockedSolidify;
                ++st.blockedSolidifyNoDest;
                break;
            }
            int dx = dest % GW, dy = dest / GW;
            float leftoverFill = fluid.fill[static_cast<size_t>(dest)];
            if (leftoverFill > kMinFillMove && dest != i) {
                LiquidCarry carry = fluid.takeLiquidCarry(dest, leftoverFill);
                float rem = fluid.relocateVolumeTopologySafe(dx, dy, leftoverFill, 0.0f, 0.0f, &carry, true);
                if (rem > 1.0e-5f) {
                    for (int n = 0; n < carry.compCount; ++n) {
                        if (carry.comps[n].amount > 0.0f && validLiquidComponentId(carry.comps[n].id))
                            fluid.addLiquidComponentAmount(dest, carry.comps[n].id, carry.comps[n].amount);
                    }
                    fluid.liquidHeat[static_cast<size_t>(dest)] += carry.heat;
                    fluid.dyeR[static_cast<size_t>(dest)] += carry.dyeR;
                    fluid.dyeG[static_cast<size_t>(dest)] += carry.dyeG;
                    fluid.dyeB[static_cast<size_t>(dest)] += carry.dyeB;
                    fluid.clearEmptyLiquidCell(dest);
                    sc.reserved[static_cast<size_t>(dest)] = 1;
                    sc.reservedMat[static_cast<size_t>(dest)] = mat;
                    ++st.blockedSolidify;
                    continue;
                }
            }
            sc.reserved[static_cast<size_t>(dest)] = 1;
            sc.reservedMat[static_cast<size_t>(dest)] = mat;
            sc.spawnByMat[mat].push_back(dest);
            PhaseScratch::SpawnJob job;
            job.src = i;
            job.dest = dest;
            job.id = id;
            job.mat = mat;
            job.massKg = static_cast<float>(pixelKg);
            sc.spawnJobs.push_back(job);
            avail -= pixelKg;
        }
    };

    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = FluidEngine::ci(x, y);
        meltPendingAt(x, y);
        if (isBlockedSolid(fluid, rigid, x, y)) {
            trySpawnPending(x, y);
            continue;
        }
        float fill = fluid.fill[static_cast<size_t>(i)];
        SubstanceId id = SUBSTANCE_NONE;
        float liquidFill = 0.0f;
        if (fill >= kMinFillMove && resolvePureLiquidForSolidification(fluid, i, id, liquidFill)) {
            PhaseProperties const &phase = phaseForSubstance(id);
            float Tm = phase.meltingPointK;
            float Lf = phase.latentHeatFusion;
            double pixelKg = solidFractionToMassKg(id, 1.0, cpm);
            double mRate = pixelKg * static_cast<double>(kMaxSolidPixelsPerSec) * static_cast<double>(dtSafe);
            float C = ThermalEngine::liquidCapacity(fluid, i);
            float E = fluid.liquidHeat[static_cast<size_t>(i)];
            float T = (C > MIN_THERMAL_CAPACITY) ? tempFromEnergy(E, C)
                : ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
            if (liquidFill >= kMinFillMove && T <= Tm + kFreezeEpsK) {
                SubstanceId pendId = fluid.solidifyPendingSubstance(i);
                if (pendId != SUBSTANCE_NONE && pendId != id) {
                    ++st.blockedSolidify;
                } else {
                    double pendingNow = static_cast<double>(fluid.solidifyPendingMassKg(i));
                    double roomMass = pixelKg - pendingNow;
                    if (!(roomMass > 1.0e-9)) {
                        // Placement-limited: keep existing pending, do not grow a hidden reservoir.
                    } else {
                    float roomJ = heatRoomAround(fluid, rigid, gas, x, y, Tm);
                    double mEnergy = static_cast<double>(roomJ) / static_cast<double>(Lf);
                    double mFill = liquidFillToMassKg(id, liquidFill, cpm);
                    double mWant = std::min({mEnergy, mFill, mRate, roomMass});
                    if (mWant > 1.0e-9) {
                        float dFill = static_cast<float>(massKgToLiquidFill(id, mWant, cpm));
                        dFill = std::min(dFill, liquidFill);
                        if (dFill >= kMinFillMove) {
                            double m = liquidFillToMassKg(id, dFill, cpm);
                            float fill0 = fluid.fill[static_cast<size_t>(i)];
                            float E0 = fluid.liquidHeat[static_cast<size_t>(i)];
                            (void)fluid.takeLiquidVolume(i, dFill);
                            float fill1 = fluid.fill[static_cast<size_t>(i)];
                            float Ekeep = (fill0 > 1.0e-8f) ? E0 * (fill1 / fill0) : 0.0f;
                            if (!(Ekeep > 0.0f) || !std::isfinite(Ekeep)) Ekeep = 0.0f;
                            fluid.liquidHeat[static_cast<size_t>(i)] = Ekeep;
                            float Eremoved = std::max(0.0f, E0 - Ekeep);
                            float Tplace = std::min(T, Tm);
                            if (!(Tplace > MIN_SAFE_TEMPERATURE_K)) Tplace = MIN_SAFE_TEMPERATURE_K;
                            float pendingHeat = static_cast<float>(m) * solidPhaseSpecificHeat(id) * Tplace;
                            if (!(pendingHeat > 0.0f) || !std::isfinite(pendingHeat)) pendingHeat = 0.0f;
                            float leftoverSensible = Eremoved - pendingHeat;
                            fluid.expectedVolume -= static_cast<double>(dFill);
                            if (!fluid.addSolidifyPendingKg(i, id, static_cast<float>(m), pendingHeat)) {
                                addLiquidComponentWithHeat(fluid, i, id, dFill, Eremoved);
                                fluid.expectedVolume += static_cast<double>(dFill);
                                ++st.blockedSolidify;
                            } else {
                                if (leftoverSensible > 0.0f) {
                                    float dumpT = std::min(Tm, T);
                                    for (int n = 0; n < 4 && leftoverSensible > 0.0f; ++n)
                                        dumpHeatBudget(fluid, rigid, gas, thermal,
                                            x + kDx[n], y + kDy[n], dumpT, leftoverSensible);
                                    if (leftoverSensible > 0.0f) {
                                        for (int n = 0; n < 4 && leftoverSensible > 0.0f; ++n) {
                                            int nx = x + kDx[n], ny = y + kDy[n];
                                            if (!FluidEngine::inside(nx, ny)) continue;
                                            int ni = FluidEngine::ci(nx, ny);
                                            if (!fluid.solid[static_cast<size_t>(ni)]) continue;
                                            dumpHeatBudget(fluid, rigid, gas, thermal, nx, ny, Tm, leftoverSensible);
                                        }
                                    }
                                } else if (leftoverSensible < 0.0f) {
                                    float need = -leftoverSensible;
                                    takeHeatAround(fluid, rigid, gas, thermal, x, y, Tplace, need, false);
                                }
                                float latent = static_cast<float>(m * static_cast<double>(Lf));
                                dumpHeatBudget(fluid, rigid, gas, thermal, x, y, Tm, latent);
                                for (int n = 0; n < 4 && latent > 0.0f; ++n)
                                    dumpHeatBudget(fluid, rigid, gas, thermal, x + kDx[n], y + kDy[n], Tm, latent);
                                if (latent > 0.0f)
                                    dumpHeatBudget(fluid, rigid, gas, thermal, x, y, Tm, latent);
                                st.massSolidifiedKg += m;
                                st.latentFusionReleasedJ += m * static_cast<double>(Lf);
                                st.liquidFillRemoved += dFill;
                                ++st.solidifyCells;
                                wakeAllEngines(fluid, gas, thermal, x, y);
                            }
                        }
                    } else if (mFill > 1.0e-9 && mEnergy <= 1.0e-9) {
                        ++st.blockedSolidify;
                        ++st.blockedSolidifyEnergy;
                    }
                    }
                }
            }
        }
        trySpawnPending(x, y);
    }

    for (int m = 1; m < MATERIAL_COUNT; ++m) {
        if (sc.spawnByMat[m].empty()) continue;
        SubstanceId sid = substanceForMaterialId(static_cast<MaterialId>(m));
        float Tm = phaseForSubstance(sid).meltingPointK;
        if (!(Tm > 1.0f)) Tm = AMBIENT_TEMPERATURE_K;
        sc.spawnUnique.clear();
        rigid.addSameMaterialWorldCells(fluid, sc.spawnByMat[m], static_cast<MaterialId>(m), Tm,
            &sc.spawnUnique);
        size_t k = 0;
        for (PhaseScratch::SpawnJob const &job : sc.spawnJobs) {
            if (static_cast<int>(job.mat) != m) continue;
            bool unique = k < sc.spawnUnique.size() && sc.spawnUnique[k].uniqueCreated;
            if (unique) {
                float heatTake = 0.0f;
                if (fluid.takeSolidifyPendingKg(job.src, job.id, job.massKg, &heatTake)) {
                    if (heatTake > 1.0e-6f)
                        rigid.commitSourcePixelState(sc.spawnUnique[k].bodyId,
                            sc.spawnUnique[k].localIndex, 1.0f, heatTake);
                    ++st.solidPixelsSpawned;
                }
            }
            int dx = job.dest % GW, dy = job.dest / GW;
            wakeAllEngines(fluid, gas, thermal, dx, dy);
            ++k;
        }
        for (RigidBodyEngine::SpawnedSourcePixel const &u : sc.spawnUnique) {
            if (!u.uniqueCreated || u.bodyId == 0) continue;
            int bi = rigid.indexOfId(u.bodyId);
            if (bi < 0 || bi >= static_cast<int>(rigid.bodies.size())) continue;
            RigidBody &b = rigid.bodies[static_cast<size_t>(bi)];
            b.vx = b.vy = b.omega = 0.0f;
            b.sleeping = true;
            b.quietTicks = rigid.sleepQuietTicks;
        }
    }

    bool meltedDirty = false;
    sc.worldCells.clear();
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        MaterialId mat = rigid.worldCellMaterial(x, y);
        if (mat == MATERIAL_EMPTY) continue;
        SubstanceId sid = substanceForMaterialId(mat);
        if (!supportsLiveSolidLiquidTransition(sid)) continue;
        sc.worldCells.push_back(FluidEngine::ci(x, y));
    }
    for (int index : sc.worldCells) {
        int x = index % GW, y = index / GW;
        RigidBodyEngine::SourcePixel site = rigid.resolveSourcePixel(fluid, x, y);
        if (!site.valid) continue;
        SubstanceId id = site.substance;
        if (!supportsLiveSolidLiquidTransition(id)) continue;
        if (!rigid.sourcePixelStillValid(site)) continue;
        PhaseProperties const &phase = phaseForSubstance(id);
        float Tm = phase.meltingPointK;
        float Lf = phase.latentHeatFusion;
        double pixelKg = solidFractionToMassKg(id, 1.0, cpm);
        if (!(pixelKg > 1.0e-9) || !(Lf > 1.0f) || !(Tm > 1.0f)) continue;
        double mRate = pixelKg * static_cast<double>(kMaxSolidPixelsPerSec) * static_cast<double>(dtSafe);
        float remain = site.fraction;
        if (remain <= 1.0e-6f) continue;
        int bi = rigid.indexOfId(site.bodyId);
        if (bi < 0 || bi >= static_cast<int>(rigid.bodies.size())) continue;
        RigidBody &b = rigid.bodies[static_cast<size_t>(bi)];
        int li = site.localIndex;
        float cap = ThermalEngine::rigidPixelCapacity(b, li);
        float Eice = site.heatJ;
        float Tice = tempFromEnergy(Eice, cap);
        wakeAllEngines(fluid, gas, thermal, x, y);
        if (!(Tice + kFreezeEpsK >= Tm)) continue;
        float Eplat = energyFromTemp(std::max(cap, MIN_THERMAL_CAPACITY), Tm);
        float excess = std::max(0.0f, Eice - Eplat);
        float neighborQ = heatAvailableAround(fluid, rigid, gas, x, y, Tm, false);
        float totalQ = excess + neighborQ;
        if (!(totalQ > 0.0f)) {
            ++st.blockedMelt;
            ++st.blockedMeltEnergy;
            continue;
        }
        double mEnergy = static_cast<double>(totalQ) / static_cast<double>(Lf);
        double mRemain = pixelKg * static_cast<double>(remain);
        double mWant = std::min({mEnergy, mRemain, mRate});
        if (!(mWant > 1.0e-9)) {
            ++st.blockedMelt;
            ++st.blockedMeltEnergy;
            continue;
        }
        int dest = findLiquidDest(fluid, rigid, gas, id, x, y, salt, false, Tm);
        if (dest < 0) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
        int dx = dest % GW, dy = dest / GW;
        if (dx == x && dy == y) {
            dest = -1;
            uint32_t bestKey = ~0u;
            float bestS = -1.0f;
            for (int n = 0; n < 4; ++n) {
                int nx = x + kDx[n], ny = y + kDy[n];
                if (!FluidEngine::inside(nx, ny) || isBlockedSolid(fluid, rigid, nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                float room = 1.0f - fluid.fill[static_cast<size_t>(ni)];
                if (room < kMinFillMove) continue;
                if (!liquidCellAcceptsCondensate(fluid, ni, id)) continue;
                uint32_t key = destTieKey(nx, ny, salt);
                if (destBetter(room, key, bestS, bestKey)) {
                    bestS = room;
                    bestKey = key;
                    dest = ni;
                }
            }
            if (dest < 0) dest = findLiquidDest(fluid, rigid, gas, id, x, y, salt, false, Tm);
            if (dest < 0) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
            dx = dest % GW; dy = dest / GW;
            if (dx == x && dy == y) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
        }
        if (isBlockedSolid(fluid, rigid, dx, dy)) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
        if (!liquidCellAcceptsCondensate(fluid, dest, id)) {
            ++st.blockedMelt;
            ++st.blockedMeltNoDest;
            continue;
        }
        float roomFill = 1.0f - fluid.fill[static_cast<size_t>(dest)];
        if (roomFill < kMinFillMove) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
        float dFill = static_cast<float>(massKgToLiquidFill(id, mWant, cpm));
        dFill = std::min(dFill, roomFill);
        if (dFill < kMinFillMove) { ++st.blockedMelt; ++st.blockedMeltNoDest; continue; }
        double m = liquidFillToMassKg(id, dFill, cpm);
        float pay = static_cast<float>(m * static_cast<double>(Lf));
        float iceTake = std::min(pay, excess);
        pay -= iceTake;
        if (pay > 0.0f)
            takeHeatAround(fluid, rigid, gas, thermal, x, y, Tm, pay, false);
        float leftoverEx = std::max(0.0f, excess - iceTake);
        float sensible = static_cast<float>(m) * thermalForSubstance(id).specificHeat * Tm;
        addLiquidComponentWithHeat(fluid, dest, id, dFill, sensible);
        remain -= static_cast<float>(m / pixelKg);
        bool removePixel = remain < kMinSolidRemain;
        if (removePixel) {
            double leftover = pixelKg * static_cast<double>(std::max(0.0f, remain));
            if (leftover > 1.0e-9) {
                float extraFill = static_cast<float>(massKgToLiquidFill(id, leftover, cpm));
                float extraRoom = 1.0f - fluid.fill[static_cast<size_t>(dest)];
                extraFill = std::min(extraFill, extraRoom);
                if (extraFill > kMinFillMove) {
                    double m2 = liquidFillToMassKg(id, extraFill, cpm);
                    addLiquidComponentWithHeat(fluid, dest, id, extraFill,
                        static_cast<float>(m2) * thermalForSubstance(id).specificHeat * Tm);
                    m += m2;
                    dFill += extraFill;
                    leftover -= m2;
                }
                if (leftover > 1.0e-9) {
                    float leftoverHeat = static_cast<float>(leftover) * solidPhaseSpecificHeat(id) * Tm
                        + leftoverEx;
                    if (fluid.addSolidifyPendingKg(dest, id, static_cast<float>(leftover), leftoverHeat)
                        || fluid.addSolidifyPendingKg(index, id, static_cast<float>(leftover), leftoverHeat)) {
                        leftoverEx = 0.0f;
                    } else {
                        remain = static_cast<float>(leftover / pixelKg);
                        removePixel = false;
                    }
                }
            }
        }
        if (removePixel) {
            remain = 0.0f;
            rigid.commitSourcePixelState(site.bodyId, li, 0.0f, 0.0f);
            meltedDirty = true;
        } else {
            float massRemain = static_cast<float>(solidFractionToMassKg(id, remain, cpm));
            float newCap = thermalCapacity(massRemain, solidPhaseSpecificHeat(id));
            float newHeat = energyFromTemp(newCap, Tm) + leftoverEx;
            rigid.commitSourcePixelState(site.bodyId, li, remain, newHeat);
            meltedDirty = true;
        }
        fluid.expectedVolume += static_cast<double>(dFill);
        st.massMeltedKg += m;
        st.latentFusionAbsorbedJ += m * static_cast<double>(Lf);
        st.liquidFillAdded += dFill;
        ++st.meltPixels;
        wakeAllEngines(fluid, gas, thermal, x, y);
        wakeAllEngines(fluid, gas, thermal, dx, dy);
    }

    if (meltedDirty) rigid.finalizeChemistryEdits(fluid);
    return st;
}

PhaseChangeTickStats stepPhaseChanges(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt)
{
    PhaseChangeTickStats st;
    float dtSafe = (dt > 1.0e-6f && std::isfinite(dt)) ? dt : PHYSICS_DT;
    static uint32_t phaseSalt = 0;
    ++phaseSalt;
    uint32_t salt = (thermal.tickNo + phaseSalt) * 0x9E3779B9u;

    bool anyLiquid = fluid.expectedVolume > 1.0e-10;
    bool anyVapor = anyEligibleGasVapor(gas);
    bool anyPending = false;
    if (!anyLiquid) {
        for (int i = 0; i < GW * GH; ++i) {
            if (fluid.solidifyPendingMassKg(i) > 1.0e-12f) { anyPending = true; break; }
        }
    }
    bool anyEligibleSolid = false;
    for (RigidBody const &b : rigid.bodies) {
        for (int li : b.occupiedLocal) {
            if (li < 0 || li >= static_cast<int>(b.mask.size())) continue;
            MaterialId mat = b.mask[static_cast<size_t>(li)];
            if (mat == MATERIAL_EMPTY) continue;
            if (supportsLiveSolidLiquidTransition(substanceForMaterialId(mat))) {
                anyEligibleSolid = true;
                break;
            }
        }
        if (anyEligibleSolid) break;
    }
    if (!anyLiquid && !anyVapor && !anyEligibleSolid && !anyPending) return st;

    if (anyLiquid || anyVapor)
        st.liquidGas = stepLiquidGasPhaseChange(fluid, rigid, gas, thermal, dtSafe, salt);
    if (anyLiquid || anyEligibleSolid || anyPending)
        st.solidLiquid = stepSolidLiquidPhaseChange(fluid, rigid, gas, thermal, dtSafe, salt);

    double vaporBeforeEdit = gas.sumWaterVapor();
    gas.handleWorldEdit(fluid);
    gas.currentWaterVapor = gas.sumWaterVapor();
    st.worldEditVaporDelta = gas.currentWaterVapor - vaporBeforeEdit;
    return st;
}
