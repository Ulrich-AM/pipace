#include "sim/SimulationQuality.h"

#include "chemistry/ReactionEngine.h"
#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/SubstanceTypes.h"
#include "thermal/ThermalEngine.h"
#include "world/WaterPhaseChange.h"

#include <cmath>
#include <fstream>
#include <sstream>

SimulationQualityProfile profileForQualityLevel(int level) {
    SimulationQualityProfile p;
    if (level <= 0) {
        p.level = 0;
        p.fluid.maxPressureIterations = 8;
        p.fluid.maxSubsteps = 2;
        p.fluid.maxLimiterPasses = 4;
        p.fluid.vorticityEnabled = false;
        p.fluid.surfaceTensionEnabled = false;
        p.fluid.sprayEnabled = false;
        p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
        p.fluid.physicsHz = 20;
        p.fluid.catchUpTicks = 1;
        p.gas.simMode = GasSimMode::Half;
        p.gas.maxSubsteps = 2;
        p.gas.sleepQuietTicks = 14;
        p.gas.sleepPressureDelta = 0.008f;
        p.rigid.velocityContactIters = 4;
        p.rigid.positionalIters = 1;
        p.rigid.sleepQuietTicks = 12;
        p.rigid.sleepLin = 0.22f;
        p.thermal.intervalTicks = 2;
        p.thermal.sleepTempEps = 0.20f;
        p.thermal.sleepAmbientEps = 0.50f;
        p.thermal.sleepQuietTicks = 8;
        p.chemistry.intervalTicks = 2;
        p.phase.intervalTicks = 2;
        return p;
    }
    if (level == 1) {
        p.level = 1;
        p.fluid.maxPressureIterations = 24;
        p.fluid.maxSubsteps = 6;
        p.fluid.maxLimiterPasses = 16;
        p.fluid.vorticityEnabled = false;
        p.fluid.surfaceTensionEnabled = true;
        p.fluid.sprayEnabled = true;
        p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
        p.fluid.physicsHz = 30;
        p.fluid.catchUpTicks = 2;
        p.gas.simMode = GasSimMode::Full;
        p.gas.maxSubsteps = 4;
        p.gas.sleepQuietTicks = 22;
        p.gas.sleepPressureDelta = 0.004f;
        p.rigid.velocityContactIters = 8;
        p.rigid.positionalIters = 1;
        p.rigid.sleepQuietTicks = 18;
        p.rigid.sleepLin = 0.16f;
        p.thermal.intervalTicks = 1;
        p.thermal.sleepTempEps = 0.08f;
        p.thermal.sleepAmbientEps = 0.25f;
        p.thermal.sleepQuietTicks = 12;
        p.chemistry.intervalTicks = 1;
        p.phase.intervalTicks = 1;
        return p;
    }
    p.level = 2;
    p.fluid.maxPressureIterations = 30;
    p.fluid.maxSubsteps = 6;
    p.fluid.maxLimiterPasses = 16;
    p.fluid.vorticityEnabled = false;
    p.fluid.surfaceTensionEnabled = true;
    p.fluid.sprayEnabled = true;
    p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
    p.fluid.physicsHz = 30;
    p.fluid.catchUpTicks = 2;
    p.gas.simMode = GasSimMode::Full;
    p.gas.maxSubsteps = 5;
    p.gas.sleepQuietTicks = 28;
    p.gas.sleepPressureDelta = 0.003f;
    p.rigid.velocityContactIters = 10;
    p.rigid.positionalIters = 2;
    p.rigid.sleepQuietTicks = 24;
    p.rigid.sleepLin = 0.12f;
    p.thermal.intervalTicks = 1;
    p.thermal.sleepTempEps = 0.05f;
    p.thermal.sleepAmbientEps = 0.15f;
    p.thermal.sleepQuietTicks = 16;
    p.chemistry.intervalTicks = 1;
    p.phase.intervalTicks = 1;
    return p;
}

void applySimulationQuality(SimulationQualityProfile const &profile,
    FluidEngine &fluid, GasEngine &gas, ThermalEngine &thermal,
    RigidBodyEngine &rigid, SimulationScheduleState &schedule)
{
    FluidQualitySettings const &f = profile.fluid;
    fluid.config.maxPressureIterations = f.maxPressureIterations;
    fluid.config.maxSubsteps = f.maxSubsteps;
    fluid.config.maxLimiterPasses = f.maxLimiterPasses;
    fluid.config.vorticityEnabled = f.vorticityEnabled;
    fluid.config.surfaceTensionEnabled = f.surfaceTensionEnabled;
    fluid.config.sprayEnabled = f.sprayEnabled;
    fluid.config.velocityAdvection = f.velocityAdvection;
    fluid.config.physicsHz = f.physicsHz;
    fluid.config.catchUpTicks = f.catchUpTicks;
    fluid.residualConsolidationEnabled = true;

    GasQualitySettings const &g = profile.gas;
    if (gas.config.simMode != GasSimMode::Off)
        gas.config.simMode = g.simMode;
    gas.config.maxSubsteps = g.maxSubsteps;
    gas.config.sleepQuietTicks = g.sleepQuietTicks;
    gas.config.sleepPressureDelta = g.sleepPressureDelta;

    ThermalQualitySettings const &th = profile.thermal;
    // enabled is a user toggle. intervalTicks is displayed; conduction skip is
    // owned by SimulationScheduleState so ThermalEngine does not double-skip.
    thermal.config.sleepTempEps = th.sleepTempEps;
    thermal.config.sleepAmbientEps = th.sleepAmbientEps;
    thermal.config.sleepQuietTicks = th.sleepQuietTicks;
    thermal.config.intervalTicks = th.intervalTicks;

    rigid.velocityContactIters = profile.rigid.velocityContactIters;
    rigid.positionalIters = profile.rigid.positionalIters;
    rigid.sleepQuietTicks = profile.rigid.sleepQuietTicks;
    rigid.sleepLin = profile.rigid.sleepLin;

    schedule.setIntervals(profile.thermal.intervalTicks,
        profile.chemistry.intervalTicks, profile.phase.intervalTicks);
}

std::string describeSimulationQuality(SimulationQualityProfile const &profile,
    SimulationScheduleState const &schedule, GasSimMode gasMode, bool thermalEnabled)
{
    char const *gasName = "Off";
    if (gasMode == GasSimMode::Half) gasName = "Half";
    else if (gasMode == GasSimMode::Full) gasName = "Full";
    char const *levelName = profile.level <= 0 ? "Performance"
        : (profile.level == 1 ? "Balanced" : "Accurate");
    std::ostringstream out;
    out << "sim quality:\n"
        << "  level = " << levelName << " (" << profile.level << ")\n"
        << "  fluid Piter=" << profile.fluid.maxPressureIterations
        << " sub=" << profile.fluid.maxSubsteps
        << " lim=" << profile.fluid.maxLimiterPasses << "\n"
        << "  gas mode = " << gasName
        << " sub=" << profile.gas.maxSubsteps << "\n"
        << "  thermal interval = " << schedule.thermalInterval
        << (thermalEnabled ? "\n" : " (Thermal Off)\n")
        << "  chemistry interval = " << schedule.chemistryInterval << "\n"
        << "  phase interval = " << schedule.phaseInterval << "\n"
        << "  rigid iters = " << profile.rigid.velocityContactIters;
    return out.str();
}

namespace {

bool nearRel(double a, double b, double rel = 0.12) {
    double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= rel * scale + 1.0e-4;
}

} // namespace

void runSimulationQualitySanity() {
    std::ofstream out(miscFile("quality_profile_sanity.tsv"), std::ios::trunc);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };

    SimulationQualityProfile perf = profileForQualityLevel(0);
    SimulationQualityProfile bal = profileForQualityLevel(1);
    SimulationQualityProfile acc = profileForQualityLevel(2);
    emit("performance_maps_knobs",
        perf.level == 0 && perf.gas.simMode == GasSimMode::Half
            && perf.thermal.intervalTicks == 2 && perf.chemistry.intervalTicks == 2
            && perf.phase.intervalTicks == 2 && perf.fluid.maxPressureIterations == 8
            && perf.rigid.velocityContactIters == 4
            && perf.gas.maxSubsteps == 2, "");
    emit("balanced_maps_current_like",
        bal.level == 1 && bal.gas.simMode == GasSimMode::Full
            && bal.thermal.intervalTicks == 1 && bal.chemistry.intervalTicks == 1
            && bal.phase.intervalTicks == 1 && bal.fluid.maxPressureIterations == 24
            && bal.rigid.velocityContactIters == 8
            && bal.gas.maxSubsteps == 4, "");
    emit("accurate_ge_balanced_effort",
        acc.level == 2 && acc.gas.simMode == GasSimMode::Full
            && acc.fluid.maxPressureIterations >= bal.fluid.maxPressureIterations
            && acc.gas.maxSubsteps >= bal.gas.maxSubsteps
            && acc.rigid.velocityContactIters >= bal.rigid.velocityContactIters
            && acc.thermal.intervalTicks == 1 && acc.chemistry.intervalTicks == 1, "");

    FluidEngine fOff;
    GasEngine gOff;
    ThermalEngine tOff;
    RigidBodyEngine rOff;
    SimulationScheduleState sched;
    gOff.config.simMode = GasSimMode::Off;
    tOff.config.enabled = false;
    applySimulationQuality(perf, fOff, gOff, tOff, rOff, sched);
    applySimulationQuality(bal, fOff, gOff, tOff, rOff, sched);
    applySimulationQuality(acc, fOff, gOff, tOff, rOff, sched);
    emit("gas_off_survives_quality", gOff.config.simMode == GasSimMode::Off, "");
    emit("thermal_off_survives_quality", tOff.config.enabled == false, "");

    FluidEngine fOn;
    GasEngine gOn;
    gOn.config.simMode = GasSimMode::Full;
    applySimulationQuality(perf, fOn, gOn, tOff, rOff, sched);
    emit("performance_sets_half_when_gas_on", gOn.config.simMode == GasSimMode::Half, "");

    SimulationScheduleState sChange;
    sChange.thermalInterval = 2;
    sChange.thermalCounter = 2;
    sChange.thermalAccumDt = 0.20f;
    sChange.setIntervals(1, 1, 1);
    emit("profile_change_clamps_accum",
        sChange.thermalAccumDt <= PHYSICS_DT + 1.0e-6f
            && sChange.thermalInterval == 1 && sChange.thermalCounter == 0,
        "accum=" + std::to_string(sChange.thermalAccumDt));

    auto runChemSecond = [](int interval) {
        FluidEngine fluid;
        RigidBodyEngine rigid;
        GasEngine gas;
        ThermalEngine thermal;
        ReactionEngine rx;
        fluid.config.walledBorders = true;
        gas.config.simMode = GasSimMode::Off;
        gas.config.boundary = GasBoundary::Sealed;
        thermal.config.enabled = false;
        gas.resetAmbient(fluid);
        int i = FluidEngine::ci(40, 40);
        gas.volume[static_cast<size_t>(i)] = 1.0f;
        GasComponentView mix{};
        mix.count = 2;
        mix.items[0] = {SUBSTANCE_HYDROGEN, 0.5f};
        mix.items[1] = {SUBSTANCE_OXYGEN, 0.25f};
        (void)gas.tryCommitGasOccupancy(i, mix);
        float cap = ThermalEngine::gasCapacity(gas, i);
        gas.heat[static_cast<size_t>(i)] = energyFromTemp(cap, 950.0f);
        SimulationScheduleState local;
        local.setIntervals(1, interval, 1);
        float extent = 0.0f;
        int ticks = 30;
        for (int n = 0; n < ticks; ++n) {
            local.beginPhysicsTick(PHYSICS_DT);
            float dtChem = 0.0f;
            if (local.takeChemistry(dtChem)) {
                rx.simulationTick(fluid, rigid, gas, thermal, dtChem);
                extent += rx.extentLastTick;
            }
        }
        return extent;
    };
    float e1 = runChemSecond(1);
    float e2 = runChemSecond(2);
    emit("chemistry_interval_2_same_extent_per_second",
        e1 > 0.05f && e2 > 0.05f && nearRel(static_cast<double>(e1), static_cast<double>(e2), 0.20),
        "e1=" + std::to_string(e1) + " e2=" + std::to_string(e2));

    auto runPhaseSecond = [](int interval) {
        FluidEngine fluid;
        RigidBodyEngine rigid;
        GasEngine gas;
        ThermalEngine thermal;
        fluid.config.walledBorders = true;
        gas.config.simMode = GasSimMode::Off;
        gas.config.boundary = GasBoundary::Sealed;
        thermal.config.enabled = false;
        gas.resetAmbient(fluid);
        int i = FluidEngine::ci(50, 50);
        fluid.fill[static_cast<size_t>(i)] = 1.0f;
        fluid.setLiquidComponentAmount(i, SUBSTANCE_WATER, 1.0f);
        float cap = ThermalEngine::liquidCapacity(fluid, i);
        fluid.liquidHeat[static_cast<size_t>(i)] = energyFromTemp(cap, 400.0f);
        SimulationScheduleState local;
        local.setIntervals(1, 1, interval);
        double boiled = 0.0;
        for (int n = 0; n < 30; ++n) {
            local.beginPhysicsTick(PHYSICS_DT);
            float dtPhase = 0.0f;
            if (local.takePhase(dtPhase))
                boiled += stepWaterPhaseChange(fluid, rigid, gas, thermal, dtPhase).massBoiledKg;
        }
        return boiled;
    };
    double b1 = runPhaseSecond(1);
    double b2 = runPhaseSecond(2);
    bool phaseOk = (b1 <= 1.0e-8 && b2 <= 1.0e-8) || nearRel(b1, b2, 0.25);
    emit("phase_interval_2_same_water_transfer",
        phaseOk, "b1=" + std::to_string(b1) + " b2=" + std::to_string(b2));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
