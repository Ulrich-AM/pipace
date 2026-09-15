#include "world/WaterPhaseChange.h"

#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/PhaseTransfer.h"
#include "substance/SubstanceRegistry.h"
#include "thermal/ThermalEngine.h"
#include "world/WorldQuery.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kDx[4] = {0, -1, 1, 0};
constexpr int kDy[4] = {-1, 0, 0, 1}; // above, left, right, below (y down)
constexpr float kHoneySkip = 1.0e-4f;
constexpr float kMaxFillPerTick = 0.35f;
constexpr float kMaxDestAtm = 48.0f;
constexpr float kMinFillMove = 1.0e-6f;
constexpr float kBoilEpsK = 0.05f;
constexpr int kSearchLimit = 64;

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

int findVaporDest(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas,
    int sx, int sy, float needAmt)
{
    if (needAmt <= GAS_MIN_AMOUNT) return -1;
    auto room = [&](int x, int y) -> float {
        if (isBlockedSolid(fluid, rigid, x, y)) return 0.0f;
        if (!gas.isAccessible(fluid, x, y)) return 0.0f;
        int i = GasEngine::ci(x, y);
        float vol = gas.availableVolume(fluid, x, y);
        if (vol < GAS_MIN_VOLUME) return 0.0f;
        return vol * kMaxDestAtm - gas.amount[static_cast<size_t>(i)];
    };
    int best = -1;
    float bestRoom = 0.0f;
    auto consider = [&](int x, int y) {
        float r = room(x, y);
        if (r > bestRoom && r > GAS_MIN_AMOUNT) {
            bestRoom = r;
            best = GasEngine::ci(x, y);
        }
    };
    for (int n = 0; n < 4; ++n) consider(sx + kDx[n], sy + kDy[n]);
    consider(sx, sy);
    if (best >= 0 && bestRoom >= std::min(needAmt, 1.0e-4f)) return best;

    std::vector<int> q;
    std::vector<uint8_t> seen(static_cast<size_t>(GW * GH), 0);
    q.reserve(kSearchLimit);
    auto push = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (seen[static_cast<size_t>(i)]) return;
        seen[static_cast<size_t>(i)] = 1;
        q.push_back(i);
    };
    push(sx, sy);
    for (size_t head = 0; head < q.size() && static_cast<int>(head) < kSearchLimit; ++head) {
        int i = q[head];
        int x = i % GW, y = i / GW;
        consider(x, y);
        if (best >= 0 && bestRoom >= needAmt) return best;
        for (int n = 0; n < 4; ++n) push(x + kDx[n], y + kDy[n]);
    }
    return (best >= 0 && bestRoom > GAS_MIN_AMOUNT) ? best : -1;
}

int findLiquidDest(FluidEngine const &fluid, RigidBodyEngine const &rigid, int sx, int sy, float needFill)
{
    auto score = [&](int x, int y) -> float {
        if (isBlockedSolid(fluid, rigid, x, y)) return -1.0f;
        int i = FluidEngine::ci(x, y);
        float fill = fluid.fill[static_cast<size_t>(i)];
        float room = 1.0f - fill;
        if (room < kMinFillMove) return -1.0f;
        if (fluid.honeyFraction(i) > kHoneySkip && fill > 1.0e-6f) return -1.0f;
        float s = room;
        if (fill > MIN_ACTIVE_FILL) s += 8.0f + fill;
        bool nearSolid = false;
        for (int n = 0; n < 4; ++n) {
            int nx = x + kDx[n], ny = y + kDy[n];
            if (!FluidEngine::inside(nx, ny)) continue;
            int ni = FluidEngine::ci(nx, ny);
            if (fluid.solid[static_cast<size_t>(ni)] || rigid.occupant[static_cast<size_t>(ni)] >= 0)
                nearSolid = true;
        }
        if (nearSolid) s += 3.0f;
        if (y > sy) s += 1.5f;
        return s;
    };
    int best = -1;
    float bestS = -1.0f;
    auto consider = [&](int x, int y) {
        float s = score(x, y);
        if (s > bestS) {
            bestS = s;
            best = FluidEngine::ci(x, y);
        }
    };
    consider(sx, sy);
    for (int n = 0; n < 4; ++n) consider(sx + kDx[n], sy + kDy[n]);
    if (best >= 0 && bestS > 0.0f) return best;

    std::vector<int> q;
    std::vector<uint8_t> seen(static_cast<size_t>(GW * GH), 0);
    q.reserve(kSearchLimit);
    auto push = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return;
        int i = FluidEngine::ci(x, y);
        if (seen[static_cast<size_t>(i)]) return;
        seen[static_cast<size_t>(i)] = 1;
        q.push_back(i);
    };
    push(sx, sy);
    for (size_t head = 0; head < q.size() && static_cast<int>(head) < kSearchLimit; ++head) {
        int i = q[head];
        int x = i % GW, y = i / GW;
        consider(x, y);
        for (int n = 0; n < 4; ++n) push(x + kDx[n], y + kDy[n]);
    }
    (void)needFill;
    return (best >= 0 && bestS > 0.0f) ? best : -1;
}

void removePureWaterFill(FluidEngine &fluid, int index, float dFill) {
    if (dFill <= 0.0f) return;
    size_t i = static_cast<size_t>(index);
    float old = fluid.fill[i];
    if (old <= 1.0e-8f) return;
    float keep = std::max(0.0f, old - dFill) / old;
    fluid.fill[i] = std::max(0.0f, old - dFill);
    fluid.dyeR[i] *= keep;
    fluid.dyeG[i] *= keep;
    fluid.dyeB[i] *= keep;
    fluid.honey[i] *= keep;
    fluid.clearEmptyLiquidCell(index);
}

double liquidWaterMassKg(FluidEngine const &fluid) {
    double m = 0.0;
    double cpm = fluid.config.cellsPerMeter;
    for (int i = 0; i < GW * GH; ++i) {
        float w = fluid.liquidComponentAmount(i, SUBSTANCE_WATER);
        if (w > 0.0f) m += liquidFillToMassKg(SUBSTANCE_WATER, w, cpm);
    }
    for (SplashParticle const &p : fluid.splashes) {
        float water = std::max(0.0f, p.volume - p.honey);
        if (water > 0.0f) m += liquidFillToMassKg(SUBSTANCE_WATER, water, cpm);
    }
    return m;
}

double vaporMassKg(GasEngine const &gas, double cpm) {
    return gasAmountToMassKg(SUBSTANCE_WATER, gas.sumWaterVapor(), cpm);
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

bool finiteTemps(FluidEngine const &fluid, RigidBodyEngine const &rigid, GasEngine const &gas) {
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        float t = ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
        if (!std::isfinite(t) || t < 0.0f) return false;
    }
    return true;
}

} // namespace

WaterPhaseTickStats stepWaterPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt)
{
    (void)dt;
    WaterPhaseTickStats st;
    if (!canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Gas)) return st;
    PhaseProperties const &phase = phaseForSubstance(SUBSTANCE_WATER);
    float Tb = phase.boilingPointK;
    float Lv = phase.latentHeatVaporization;
    if (!(Tb > 1.0f) || !(Lv > 1.0f)) return st;
    double cpm = fluid.config.cellsPerMeter;
    float cpVapor = gasPhaseSpecificHeat(SUBSTANCE_WATER);

    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = FluidEngine::ci(x, y);
        if (isBlockedSolid(fluid, rigid, x, y)) continue;
        float fill = fluid.fill[static_cast<size_t>(i)];
        if (fill < MIN_ACTIVE_FILL) continue;
        if (fluid.honeyFraction(i) > kHoneySkip) continue;
        float waterFill = fluid.liquidComponentAmount(i, SUBSTANCE_WATER);
        if (waterFill < kMinFillMove) continue;
        float C = ThermalEngine::liquidCapacity(fluid, i);
        if (C < MIN_THERMAL_CAPACITY) continue;
        float E = fluid.liquidHeat[static_cast<size_t>(i)];
        float T = tempFromEnergy(E, C);
        if (!(T + kBoilEpsK >= Tb)) continue;
        float Eplat = energyFromTemp(C, Tb);
        float excess = E - Eplat;
        if (!(excess > 0.0f)) {
            if (T > Tb) fluid.liquidHeat[static_cast<size_t>(i)] = Eplat;
            continue;
        }
        double mEnergy = static_cast<double>(excess) / static_cast<double>(Lv);
        double mFill = liquidFillToMassKg(SUBSTANCE_WATER, waterFill, cpm);
        double mRate = liquidFillToMassKg(SUBSTANCE_WATER, kMaxFillPerTick, cpm);
        double mWant = std::min({mEnergy, mFill, mRate});
        if (!(mWant > 1.0e-9)) continue;
        double vaporWant = massKgToGasAmount(SUBSTANCE_WATER, mWant, cpm);
        int dest = findVaporDest(fluid, rigid, gas, x, y, static_cast<float>(vaporWant));
        if (dest < 0) {
            ++st.blockedBoil;
            continue;
        }
        int dx = dest % GW, dy = dest / GW;
        float vol = gas.availableVolume(fluid, dx, dy);
        float room = vol * kMaxDestAtm - gas.amount[static_cast<size_t>(dest)];
        if (room <= GAS_MIN_AMOUNT) {
            ++st.blockedBoil;
            continue;
        }
        double vaporGot = std::min(vaporWant, static_cast<double>(room));
        double m = gasAmountToMassKg(SUBSTANCE_WATER, vaporGot, cpm);
        if (!(m > 1.0e-9)) {
            ++st.blockedBoil;
            continue;
        }
        float dFill = static_cast<float>(massKgToLiquidFill(SUBSTANCE_WATER, m, cpm));
        dFill = std::min(dFill, waterFill);
        if (dFill < kMinFillMove) continue;
        m = liquidFillToMassKg(SUBSTANCE_WATER, dFill, cpm);
        vaporGot = massKgToGasAmount(SUBSTANCE_WATER, m, cpm);

        removePureWaterFill(fluid, i, dFill);
        float leftoverExcess = std::max(0.0f, excess - static_cast<float>(m * static_cast<double>(Lv)));
        float newC = ThermalEngine::liquidCapacity(fluid, i);
        fluid.liquidHeat[static_cast<size_t>(i)] = (newC > MIN_THERMAL_CAPACITY)
            ? energyFromTemp(newC, Tb) + leftoverExcess : 0.0f;

        float sensible = static_cast<float>(m) * cpVapor * Tb;
        gas.addWaterVapor(dest, static_cast<float>(vaporGot));
        gas.heat[static_cast<size_t>(dest)] += sensible;
        if (!std::isfinite(gas.heat[static_cast<size_t>(dest)]) || gas.heat[static_cast<size_t>(dest)] < 0.0f)
            gas.heat[static_cast<size_t>(dest)] = 0.0f;

        fluid.expectedVolume -= static_cast<double>(dFill);
        gas.expectedAmount += vaporGot;
        gas.expectedWaterVapor += vaporGot;
        st.massBoiledKg += m;
        st.latentAbsorbedJ += m * static_cast<double>(Lv);
        st.vaporPlacedAmount += vaporGot;
        st.liquidFillRemoved += dFill;
        ++st.boilCells;
        wakeAllEngines(fluid, gas, thermal, x, y);
        wakeAllEngines(fluid, gas, thermal, dx, dy);
    }

    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int gi = GasEngine::ci(x, y);
        float vap = gas.vaporAmount(gi);
        if (vap <= GAS_MIN_AMOUNT) continue;
        if (isBlockedSolid(fluid, rigid, x, y)) continue;
        float Tgas = ThermalEngine::gasTempK(gas, gi);
        if (!(Tgas <= Tb + kBoilEpsK)) continue;
        double mAvail = gasAmountToMassKg(SUBSTANCE_WATER, vap, cpm);
        double mRate = liquidFillToMassKg(SUBSTANCE_WATER, kMaxFillPerTick, cpm);
        double mWant = std::min(mAvail, mRate);
        if (!(mWant > 1.0e-9)) continue;
        int dest = findLiquidDest(fluid, rigid, x, y, static_cast<float>(massKgToLiquidFill(SUBSTANCE_WATER, mWant, cpm)));
        if (dest < 0) {
            ++st.blockedCondense;
            continue;
        }
        int dx = dest % GW, dy = dest / GW;
        float roomFill = 1.0f - fluid.fill[static_cast<size_t>(dest)];
        if (roomFill < kMinFillMove) {
            ++st.blockedCondense;
            continue;
        }
        float Cdest = ThermalEngine::liquidCapacity(fluid, dest);
        float Edest = fluid.liquidHeat[static_cast<size_t>(dest)];
        float Tdest = (Cdest > MIN_THERMAL_CAPACITY) ? tempFromEnergy(Edest, Cdest) : Tgas;
        if (Tdest > Tb + kBoilEpsK) {
            ++st.blockedCondense;
            continue;
        }
        if (fluid.solid[static_cast<size_t>(dest)]) {
            ++st.blockedCondense;
            continue;
        }
        double roomJ = heatRoomToK(Edest, Cdest, Tb);
        roomJ += heatRoomToK(gas.heat[static_cast<size_t>(gi)], ThermalEngine::gasCapacity(gas, gi), Tb);
        for (int n = 0; n < 4; ++n) {
            int nx = dx + kDx[n], ny = dy + kDy[n];
            if (!FluidEngine::inside(nx, ny)) continue;
            int ni = FluidEngine::ci(nx, ny);
            if (fluid.solid[static_cast<size_t>(ni)]) {
                roomJ += heatRoomToK(fluid.solidHeat[static_cast<size_t>(ni)],
                    ThermalEngine::wallCapacity(fluid, ni), Tb);
            } else if (gas.volume[static_cast<size_t>(ni)] >= GAS_MIN_VOLUME) {
                roomJ += heatRoomToK(gas.heat[static_cast<size_t>(ni)],
                    ThermalEngine::gasCapacity(gas, ni), Tb);
            }
        }
        double mHeat = (Lv > 1.0f) ? roomJ / static_cast<double>(Lv) : 0.0;
        double m = std::min({mWant, mHeat});
        if (!(m > 1.0e-9)) {
            ++st.blockedCondense;
            continue;
        }
        float dFill = static_cast<float>(massKgToLiquidFill(SUBSTANCE_WATER, m, cpm));
        dFill = std::min(dFill, roomFill);
        if (dFill < kMinFillMove) continue;
        m = liquidFillToMassKg(SUBSTANCE_WATER, dFill, cpm);
        float vaporTake = static_cast<float>(massKgToGasAmount(SUBSTANCE_WATER, m, cpm));
        float taken = gas.takeWaterVapor(gi, vaporTake);
        if (taken <= GAS_MIN_AMOUNT) {
            ++st.blockedCondense;
            continue;
        }
        m = gasAmountToMassKg(SUBSTANCE_WATER, taken, cpm);
        dFill = static_cast<float>(massKgToLiquidFill(SUBSTANCE_WATER, m, cpm));
        float Tplace = std::clamp(std::min(Tgas, Tb), MIN_SAFE_TEMPERATURE_K, Tb);
        float sensible = static_cast<float>(m) * thermalForSubstance(SUBSTANCE_WATER).specificHeat * Tplace;
        float latentBudget = static_cast<float>(m * static_cast<double>(Lv));
        fluid.addLiquidFill(dest, dFill, sensible);
        addHeatToK(fluid.liquidHeat[static_cast<size_t>(dest)],
            ThermalEngine::liquidCapacity(fluid, dest), Tb, latentBudget);
        addHeatToK(gas.heat[static_cast<size_t>(gi)], ThermalEngine::gasCapacity(gas, gi), Tb, latentBudget);
        for (int n = 0; n < 4; ++n) {
            if (!(latentBudget > 0.0f)) break;
            int nx = dx + kDx[n], ny = dy + kDy[n];
            if (!FluidEngine::inside(nx, ny)) continue;
            int ni = FluidEngine::ci(nx, ny);
            if (fluid.solid[static_cast<size_t>(ni)]) {
                addHeatToK(fluid.solidHeat[static_cast<size_t>(ni)],
                    ThermalEngine::wallCapacity(fluid, ni), Tb, latentBudget);
                thermal.wakeCell(nx, ny);
            } else if (gas.volume[static_cast<size_t>(ni)] >= GAS_MIN_VOLUME) {
                addHeatToK(gas.heat[static_cast<size_t>(ni)],
                    ThermalEngine::gasCapacity(gas, ni), Tb, latentBudget);
                gas.wakeAt(nx, ny);
            }
        }
        if (latentBudget > 0.0f) {
            // Remainder stays in destination liquid rather than being deleted.
            fluid.liquidHeat[static_cast<size_t>(dest)] += latentBudget;
        }
        int xw = dest % GW, yw = dest / GW;
        fluid.waterShade[static_cast<size_t>(dest)] = fluid.makeShade(xw, yw);
        fluid.expectedVolume += static_cast<double>(dFill);
        gas.expectedAmount -= taken;
        gas.expectedWaterVapor -= taken;
        st.massCondensedKg += m;
        st.latentReleasedJ += m * static_cast<double>(Lv);
        st.liquidFillAdded += dFill;
        ++st.condenseCells;
        wakeAllEngines(fluid, gas, thermal, x, y);
        wakeAllEngines(fluid, gas, thermal, dx, dy);
    }

    gas.handleWorldEdit(fluid);
    gas.currentWaterVapor = gas.sumWaterVapor();
    return st;
}

void runWaterPhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal)
{
    std::ofstream out(miscFile("water_phase_diag.tsv"));
    out << std::setprecision(10);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto f8 = [](double v) {
        std::ostringstream o;
        o << std::setprecision(10) << v;
        return o.str();
    };

    auto worldTick = [&]() -> WaterPhaseTickStats {
        rigid.step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gas.simulationTick(fluid);
        rigid.gatherFluidForces(fluid);
        gas.applyPressureForces(rigid, fluid);
        thermal.simulationTick(fluid, rigid, gas, PHYSICS_DT);
        return stepWaterPhaseChange(fluid, rigid, gas, thermal, PHYSICS_DT);
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
        gas.config.boundary = GasBoundary::Sealed;
        gas.config.simMode = GasSimMode::Full;
        thermal.config.enabled = true;
        thermal.config.intervalTicks = 1;
        box(40, 25, 80, 70);
        gas.resetAmbient(fluid);
        thermal.seedAmbient(fluid, rigid, gas);
    };

    PhaseProperties const &phase = phaseForSubstance(SUBSTANCE_WATER);
    float Tb = phase.boilingPointK;
    double cpm = fluid.config.cellsPerMeter;

    resetSealed();
    int poolX0 = 50, poolX1 = 70, poolY = 60;
    double startFill = 0.0;
    for (int x = poolX0; x <= poolX1; ++x) {
        int i = FluidEngine::ci(x, poolY);
        fluid.fill[static_cast<size_t>(i)] = 1.0f;
        fluid.honey[static_cast<size_t>(i)] = 0.0f;
        startFill += 1.0;
        float C = ThermalEngine::liquidCapacity(fluid, i);
        fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(C, Tb + 25.0f);
        fluid.wakeChunkAtCell(x, poolY);
        thermal.wakeCell(x, poolY);
    }
    fluid.expectedVolume = startFill;
    fluid.wakeAllFluidChunks();
    gas.handleWorldEdit(fluid);
    double mass0 = liquidWaterMassKg(fluid);
    double vap0 = vaporMassKg(gas, cpm);
    WaterPhaseTickStats acc{};
    float tMax = 0.0f;
    float tMinHot = 1.0e9f;
    bool boiled = false;
    for (int n = 0; n < 40; ++n) {
        for (int x = poolX0; x <= poolX1; ++x) {
            int i = FluidEngine::ci(x, poolY);
            if (fluid.fill[static_cast<size_t>(i)] < MIN_ACTIVE_FILL) continue;
            float C = ThermalEngine::liquidCapacity(fluid, i);
            float T = ThermalEngine::liquidTempK(fluid, i);
            if (C > MIN_THERMAL_CAPACITY && T < Tb + 8.0f)
                fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(C, Tb + 12.0f);
        }
        WaterPhaseTickStats s = {};
        rigid.step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gas.simulationTick(fluid);
        thermal.simulationTick(fluid, rigid, gas, PHYSICS_DT);
        s = stepWaterPhaseChange(fluid, rigid, gas, thermal, PHYSICS_DT);
        acc.massBoiledKg += s.massBoiledKg;
        acc.massCondensedKg += s.massCondensedKg;
        acc.vaporPlacedAmount += s.vaporPlacedAmount;
        acc.boilCells += s.boilCells;
        if (s.massBoiledKg > 0.0) boiled = true;
        for (int x = poolX0; x <= poolX1; ++x) {
            int i = FluidEngine::ci(x, poolY);
            if (fluid.fill[static_cast<size_t>(i)] < MIN_ACTIVE_FILL) continue;
            float T = ThermalEngine::liquidTempK(fluid, i);
            tMax = std::max(tMax, T);
            tMinHot = std::min(tMinHot, T);
        }
    }
    double mass1 = liquidWaterMassKg(fluid);
    double vap1 = vaporMassKg(gas, cpm);
    double dLiq = mass0 - mass1;
    double dVap = vap1 - vap0;
    emit("hot_water_vaporizes", boiled && acc.massBoiledKg > 1.0e-4 && vap1 > vap0,
        "boiled_kg=" + f8(acc.massBoiledKg) + " vapor_kg=" + f8(vap1));
    emit("temperature_near_boiling",
        boiled && tMax < Tb + 40.0f && tMinHot > Tb - 25.0f,
        "Tmin=" + f8(tMinHot) + " Tmax=" + f8(tMax) + " Tb=" + f8(Tb));
    emit("boil_mass_conserved",
        std::abs(dLiq - dVap) <= std::max(1.0e-4, 0.02 * std::max(dLiq, dVap)),
        "dLiq=" + f8(dLiq) + " dVap=" + f8(dVap)
            + " boiled=" + f8(acc.massBoiledKg) + " condensed=" + f8(acc.massCondensedKg));
    emit("no_steam_id",
        substanceFromInternalName("steam") == SUBSTANCE_NONE && SUBSTANCE_COUNT == 8, "");

    bool vaporIsWaterGas = false;
    for (int y = 25; y <= 70 && !vaporIsWaterGas; ++y) for (int x = 40; x <= 80; ++x) {
        if (gas.vaporAmount(GasEngine::ci(x, y)) <= 1.0e-5f) continue;
        MatterSample s = sampleMatterAt(fluid, rigid, gas, x, y);
        if (s.identity.substance == SUBSTANCE_WATER && s.identity.phase == MatterPhase::Gas)
            vaporIsWaterGas = true;
        if (s.mixture && s.vaporFraction > 0.0f && s.identity.phase == MatterPhase::Gas)
            vaporIsWaterGas = true;
    }
    emit("vapor_reports_water_gas", vaporIsWaterGas, "");

    float pMax = 0.0f;
    for (int y = 26; y < 70; ++y) for (int x = 41; x < 80; ++x) {
        int i = GasEngine::ci(x, y);
        if (gas.volume[static_cast<size_t>(i)] >= GAS_MIN_VOLUME)
            pMax = std::max(pMax, gas.pressure[static_cast<size_t>(i)]);
    }
    emit("sealed_pressure_rises", boiled && pMax > 1.02f, "pmax_atm=" + f8(pMax));

    bool vaporInSolid = false;
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = FluidEngine::ci(x, y);
        if (!fluid.solid[static_cast<size_t>(i)] && rigid.occupant[static_cast<size_t>(i)] < 0) continue;
        if (gas.vaporAmount(i) > GAS_MIN_AMOUNT && gas.volume[static_cast<size_t>(i)] < GAS_MIN_VOLUME)
            vaporInSolid = true;
    }
    emit("vapor_not_inside_solids", !vaporInSolid, "");

    double massBeforeCool = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    double storedBefore = thermal.totalThermalEnergy(fluid, rigid, gas);
    double latentBefore = vaporMassKg(gas, cpm) * static_cast<double>(phase.latentHeatVaporization);
    WaterPhaseTickStats condAcc{};
    for (int n = 0; n < 80; ++n) {
        for (int y = 26; y < 70; ++y) for (int x = 41; x < 80; ++x) {
            int i = FluidEngine::ci(x, y);
            float Cg = ThermalEngine::gasCapacity(gas, i);
            if (Cg > MIN_THERMAL_CAPACITY)
                gas.heat[static_cast<size_t>(i)] = energyFromTemp(Cg, std::max(200.0f, Tb - 40.0f));
            float Cl = ThermalEngine::liquidCapacity(fluid, i);
            if (Cl > MIN_THERMAL_CAPACITY) {
                float Tl = ThermalEngine::liquidTempK(fluid, i);
                if (Tl > Tb - 15.0f)
                    fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(Cl, Tb - 20.0f);
            }
        }
        thermal.wakeRect(40, 25, 80, 70);
        WaterPhaseTickStats s = worldTick();
        condAcc.massCondensedKg += s.massCondensedKg;
        condAcc.latentReleasedJ += s.latentReleasedJ;
    }
    double massAfterCool = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    double vapAfter = vaporMassKg(gas, cpm);
    emit("cooling_condenses", condAcc.massCondensedKg > 1.0e-5 || vapAfter < vap1 * 0.85,
        "condensed_kg=" + f8(condAcc.massCondensedKg) + " vap_after=" + f8(vapAfter));
    emit("condense_mass_cycle",
        std::abs(massAfterCool - massBeforeCool) <= std::max(1.0e-3, 0.03 * massBeforeCool),
        "before=" + f8(massBeforeCool) + " after=" + f8(massAfterCool));
    emit("latent_released_on_condense",
        condAcc.latentReleasedJ > 1.0 || vapAfter < vap1,
        "J=" + f8(condAcc.latentReleasedJ));
    double storedAfter = thermal.totalThermalEnergy(fluid, rigid, gas);
    double latentAfter = vapAfter * static_cast<double>(phase.latentHeatVaporization);
    double e0 = storedBefore + latentBefore;
    double e1 = storedAfter + latentAfter;
    double eRel = (std::abs(e0) > 1.0) ? std::abs(e1 - e0) / std::abs(e0) : std::abs(e1 - e0);
    emit("closed_energy_with_latent",
        eRel < 0.15 || !boiled,
        "rel=" + f8(eRel) + " E0=" + f8(e0) + " E1=" + f8(e1));
    emit("no_nan_inf_negK", finiteTemps(fluid, rigid, gas), "");

    // Repeated liquid->gas->liquid on a tiny sealed pocket
    resetSealed();
    int cx = 60, cy = 55;
    fluid.fill[static_cast<size_t>(FluidEngine::ci(cx, cy))] = 1.0f;
    fluid.expectedVolume = 1.0;
    thermal.seedAmbient(fluid, rigid, gas);
    float C0 = ThermalEngine::liquidCapacity(fluid, FluidEngine::ci(cx, cy));
    fluid.liquidHeat[static_cast<size_t>(FluidEngine::ci(cx, cy))] = energyFromTemp(C0, Tb + 30.0f);
    double cycleMass0 = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    for (int cycle = 0; cycle < 8; ++cycle) {
        for (int n = 0; n < 12; ++n) {
            int i = FluidEngine::ci(cx, cy);
            float C = ThermalEngine::liquidCapacity(fluid, i);
            if (C > MIN_THERMAL_CAPACITY)
                fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(C, Tb + 20.0f);
            worldTick();
        }
        for (int n = 0; n < 20; ++n) {
            for (int y = 26; y < 70; ++y) for (int x = 41; x < 80; ++x) {
                int i = FluidEngine::ci(x, y);
                float Cg = ThermalEngine::gasCapacity(gas, i);
                if (Cg > MIN_THERMAL_CAPACITY)
                    gas.heat[static_cast<size_t>(i)] = energyFromTemp(Cg, Tb - 30.0f);
                float Cl = ThermalEngine::liquidCapacity(fluid, i);
                if (Cl > MIN_THERMAL_CAPACITY)
                    fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(Cl, Tb - 25.0f);
            }
            worldTick();
        }
    }
    double cycleMass1 = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    emit("repeat_phase_cycle_mass",
        std::abs(cycleMass1 - cycleMass0) <= std::max(0.05, 0.05 * cycleMass0),
        "start=" + f8(cycleMass0) + " end=" + f8(cycleMass1));

    // Honey mixture must not boil
    resetSealed();
    int hx = 55, hy = 60;
    fluid.fill[static_cast<size_t>(FluidEngine::ci(hx, hy))] = 1.0f;
    fluid.honey[static_cast<size_t>(FluidEngine::ci(hx, hy))] = 0.4f;
    fluid.expectedVolume = 1.0;
    float Ch = ThermalEngine::liquidCapacity(fluid, FluidEngine::ci(hx, hy));
    fluid.liquidHeat[static_cast<size_t>(FluidEngine::ci(hx, hy))] = energyFromTemp(Ch, Tb + 40.0f);
    gas.handleWorldEdit(fluid);
    WaterPhaseTickStats honeySt{};
    for (int n = 0; n < 15; ++n) {
        float C = ThermalEngine::liquidCapacity(fluid, FluidEngine::ci(hx, hy));
        if (C > MIN_THERMAL_CAPACITY)
            fluid.liquidHeat[static_cast<size_t>(FluidEngine::ci(hx, hy))] = energyFromTemp(C, Tb + 40.0f);
        honeySt.massBoiledKg += worldTick().massBoiledKg;
    }
    emit("honey_mixture_does_not_boil",
        honeySt.massBoiledKg < 1.0e-8 && vaporMassKg(gas, cpm) < 1.0e-8,
        "boiled=" + f8(honeySt.massBoiledKg));

    // Open border accounting: vapor may leave; escapedWaterVapor tracks it
    rigid.clear();
    fluid.clearWorld();
    gas.config.boundary = GasBoundary::OpenAmbient;
    gas.resetAmbient(fluid);
    thermal.seedAmbient(fluid, rigid, gas);
    int ox = GW / 2, oy = GH / 2;
    fluid.fill[static_cast<size_t>(FluidEngine::ci(ox, oy))] = 1.0f;
    fluid.expectedVolume = 1.0;
    float Co = ThermalEngine::liquidCapacity(fluid, FluidEngine::ci(ox, oy));
    fluid.liquidHeat[static_cast<size_t>(FluidEngine::ci(ox, oy))] = energyFromTemp(Co, Tb + 40.0f);
    double openMass0 = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    for (int n = 0; n < 25; ++n) {
        float C = ThermalEngine::liquidCapacity(fluid, FluidEngine::ci(ox, oy));
        if (C > MIN_THERMAL_CAPACITY)
            fluid.liquidHeat[static_cast<size_t>(FluidEngine::ci(ox, oy))] = energyFromTemp(C, Tb + 20.0f);
        worldTick();
    }
    double openMass1 = liquidWaterMassKg(fluid) + vaporMassKg(gas, cpm);
    double escapedVaporMass = gasAmountToMassKg(SUBSTANCE_WATER, gas.escapedWaterVapor, cpm);
    double openErr = std::abs((openMass1 + escapedVaporMass) - openMass0);
    emit("open_boundary_vapor_accounted",
        openErr <= std::max(0.05, 0.08 * std::max(1.0, openMass0)),
        "start=" + f8(openMass0) + " remain=" + f8(openMass1)
            + " escaped_vapor_kg=" + f8(escapedVaporMass));

    emit("no_ice_id", substanceFromInternalName("ice") == SUBSTANCE_NONE, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
