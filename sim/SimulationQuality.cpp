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

namespace {

char const *levelLetter(int level) {
    if (level <= 0) return "P";
    if (level == 1) return "B";
    return "A";
}

char const *levelName(int level) {
    if (level <= 0) return "Performance";
    if (level == 1) return "Balanced";
    return "Accurate";
}

bool sameFluid(FluidQualitySettings const &a, FluidQualitySettings const &b) {
    return a.maxPressureIterations == b.maxPressureIterations
        && a.maxSubsteps == b.maxSubsteps
        && a.maxLimiterPasses == b.maxLimiterPasses
        && a.vorticityEnabled == b.vorticityEnabled
        && a.surfaceTensionEnabled == b.surfaceTensionEnabled
        && a.sprayEnabled == b.sprayEnabled
        && a.velocityAdvection == b.velocityAdvection;
}

bool sameGas(GasQualitySettings const &a, GasQualitySettings const &b) {
    return a.simMode == b.simMode
        && a.maxSubsteps == b.maxSubsteps
        && a.sleepQuietTicks == b.sleepQuietTicks
        && a.sleepPressureDelta == b.sleepPressureDelta;
}

bool sameRigid(RigidQualitySettings const &a, RigidQualitySettings const &b) {
    return a.velocityContactIters == b.velocityContactIters
        && a.positionalIters == b.positionalIters
        && a.sleepQuietTicks == b.sleepQuietTicks
        && a.sleepLin == b.sleepLin;
}

bool sameThermal(ThermalQualitySettings const &a, ThermalQualitySettings const &b) {
    return a.intervalTicks == b.intervalTicks
        && a.sleepTempEps == b.sleepTempEps
        && a.sleepAmbientEps == b.sleepAmbientEps
        && a.sleepQuietTicks == b.sleepQuietTicks;
}

bool nearRel(double a, double b, double rel = 0.12) {
    double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= rel * scale + 1.0e-4;
}

} // namespace

SimulationQualityProfile profileForQualityLevel(int level) {
    SimulationQualityProfile p;
    int lv = std::max(0, std::min(2, level));
    p.level = lv;
    p.custom = false;
    p.levels = SimulationQualityLevels{lv, lv, lv, lv, lv, lv};
    if (lv == 0) {
        p.fluid.maxPressureIterations = 8;
        p.fluid.maxSubsteps = 2;
        p.fluid.maxLimiterPasses = 4;
        p.fluid.vorticityEnabled = false;
        p.fluid.surfaceTensionEnabled = false;
        p.fluid.sprayEnabled = false;
        p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
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
    if (lv == 1) {
        p.fluid.maxPressureIterations = 24;
        p.fluid.maxSubsteps = 6;
        p.fluid.maxLimiterPasses = 16;
        p.fluid.vorticityEnabled = false;
        p.fluid.surfaceTensionEnabled = true;
        p.fluid.sprayEnabled = true;
        p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
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
    p.fluid.maxPressureIterations = 30;
    p.fluid.maxSubsteps = 6;
    p.fluid.maxLimiterPasses = 16;
    p.fluid.vorticityEnabled = false;
    p.fluid.surfaceTensionEnabled = true;
    p.fluid.sprayEnabled = true;
    p.fluid.velocityAdvection = VelocityAdvection::SemiLagrangian;
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

SimulationQualityProfile profileForCustomLevels(SimulationQualityLevels levels) {
    levels.clampAll();
    SimulationQualityProfile fluidP = profileForQualityLevel(levels.fluid);
    SimulationQualityProfile gasP = profileForQualityLevel(levels.gas);
    SimulationQualityProfile rigidP = profileForQualityLevel(levels.rigid);
    SimulationQualityProfile thermalP = profileForQualityLevel(levels.thermal);
    SimulationQualityProfile chemP = profileForQualityLevel(levels.chemistry);
    SimulationQualityProfile phaseP = profileForQualityLevel(levels.phase);
    SimulationQualityProfile p;
    p.level = -1;
    p.custom = true;
    p.levels = levels;
    p.fluid = fluidP.fluid;
    p.gas = gasP.gas;
    p.rigid = rigidP.rigid;
    p.thermal = thermalP.thermal;
    p.chemistry = chemP.chemistry;
    p.phase = phaseP.phase;
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
    std::ostringstream out;
    if (profile.custom) {
        out << "sim quality: Custom\n"
            << "fluid=" << levelLetter(profile.levels.fluid)
            << " gas=" << levelLetter(profile.levels.gas)
            << " rigid=" << levelLetter(profile.levels.rigid)
            << " thermal=" << levelLetter(profile.levels.thermal)
            << " chemistry=" << levelLetter(profile.levels.chemistry)
            << " phase=" << levelLetter(profile.levels.phase);
    } else {
        out << "sim quality: " << levelName(profile.level);
    }
    if (gasMode == GasSimMode::Off) out << "\nGas Off";
    else if (gasMode == GasSimMode::Half) out << "\nGas Half";
    if (!thermalEnabled) out << "\nThermal Off";
    out << "\n  fluid Piter=" << profile.fluid.maxPressureIterations
        << " sub=" << profile.fluid.maxSubsteps
        << " lim=" << profile.fluid.maxLimiterPasses
        << "\n  thermal interval = " << schedule.thermalInterval
        << "  chemistry = " << schedule.chemistryInterval
        << "  phase = " << schedule.phaseInterval;
    return out.str();
}

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
        perf.level == 0 && !perf.custom
            && perf.gas.simMode == GasSimMode::Half
            && perf.thermal.intervalTicks == 2 && perf.chemistry.intervalTicks == 2
            && perf.phase.intervalTicks == 2 && perf.fluid.maxPressureIterations == 8
            && perf.rigid.velocityContactIters == 4
            && perf.gas.maxSubsteps == 2, "");
    emit("balanced_maps_current_like",
        bal.level == 1 && !bal.custom
            && bal.gas.simMode == GasSimMode::Full
            && bal.thermal.intervalTicks == 1 && bal.chemistry.intervalTicks == 1
            && bal.phase.intervalTicks == 1 && bal.fluid.maxPressureIterations == 24
            && bal.rigid.velocityContactIters == 8
            && bal.gas.maxSubsteps == 4, "");
    emit("accurate_ge_balanced_effort",
        acc.level == 2 && !acc.custom
            && acc.gas.simMode == GasSimMode::Full
            && acc.fluid.maxPressureIterations >= bal.fluid.maxPressureIterations
            && acc.gas.maxSubsteps >= bal.gas.maxSubsteps
            && acc.rigid.velocityContactIters >= bal.rigid.velocityContactIters
            && acc.thermal.intervalTicks == 1 && acc.chemistry.intervalTicks == 1, "");

    emit("global_performance_composes_all",
        sameFluid(perf.fluid, profileForQualityLevel(0).fluid)
            && sameGas(perf.gas, profileForQualityLevel(0).gas)
            && sameRigid(perf.rigid, profileForQualityLevel(0).rigid)
            && sameThermal(perf.thermal, profileForQualityLevel(0).thermal)
            && perf.chemistry.intervalTicks == 2 && perf.phase.intervalTicks == 2, "");
    emit("global_balanced_composes_all",
        sameFluid(bal.fluid, profileForQualityLevel(1).fluid)
            && sameGas(bal.gas, profileForQualityLevel(1).gas)
            && sameRigid(bal.rigid, profileForQualityLevel(1).rigid)
            && sameThermal(bal.thermal, profileForQualityLevel(1).thermal)
            && bal.chemistry.intervalTicks == 1 && bal.phase.intervalTicks == 1, "");
    emit("global_accurate_composes_all",
        sameFluid(acc.fluid, profileForQualityLevel(2).fluid)
            && sameGas(acc.gas, profileForQualityLevel(2).gas)
            && sameRigid(acc.rigid, profileForQualityLevel(2).rigid)
            && sameThermal(acc.thermal, profileForQualityLevel(2).thermal)
            && acc.chemistry.intervalTicks == 1 && acc.phase.intervalTicks == 1, "");

    SimulationQualityLevels mix;
    mix.fluid = 2;
    mix.gas = 0;
    mix.rigid = 1;
    mix.thermal = 0;
    mix.chemistry = 0;
    mix.phase = 1;
    SimulationQualityProfile custom = profileForCustomLevels(mix);
    emit("custom_mixed_matches_tables",
        custom.custom
            && sameFluid(custom.fluid, acc.fluid)
            && sameGas(custom.gas, perf.gas)
            && sameRigid(custom.rigid, bal.rigid)
            && sameThermal(custom.thermal, perf.thermal)
            && custom.chemistry.intervalTicks == perf.chemistry.intervalTicks
            && custom.phase.intervalTicks == bal.phase.intervalTicks, "");

    SimulationQualityLevels saved = mix;
    (void)profileForQualityLevel(1);
    SimulationQualityProfile customAgain = profileForCustomLevels(saved);
    emit("custom_levels_survive_leave",
        saved.fluid == 2 && saved.gas == 0 && saved.rigid == 1
            && saved.thermal == 0 && saved.chemistry == 0 && saved.phase == 1
            && customAgain.custom && sameFluid(customAgain.fluid, acc.fluid)
            && sameGas(customAgain.gas, perf.gas), "");

    FluidEngine fOff;
    GasEngine gOff;
    ThermalEngine tOff;
    RigidBodyEngine rOff;
    SimulationScheduleState sched;
    fOff.config.physicsHz = 30;
    fOff.config.catchUpTicks = 2;
    gOff.config.simMode = GasSimMode::Off;
    tOff.config.enabled = false;
    applySimulationQuality(perf, fOff, gOff, tOff, rOff, sched);
    applySimulationQuality(bal, fOff, gOff, tOff, rOff, sched);
    applySimulationQuality(acc, fOff, gOff, tOff, rOff, sched);
    applySimulationQuality(custom, fOff, gOff, tOff, rOff, sched);
    emit("gas_off_survives_quality", gOff.config.simMode == GasSimMode::Off, "");
    emit("thermal_off_survives_quality", tOff.config.enabled == false, "");
    emit("quality_leaves_physicsHz", fOff.config.physicsHz == 30,
        "hz=" + std::to_string(fOff.config.physicsHz));
    emit("quality_leaves_catchUpTicks", fOff.config.catchUpTicks == 2,
        "catch=" + std::to_string(fOff.config.catchUpTicks));

    fOff.config.physicsHz = 20;
    fOff.config.catchUpTicks = 1;
    applySimulationQuality(perf, fOff, gOff, tOff, rOff, sched);
    emit("quality_leaves_manual_20hz",
        fOff.config.physicsHz == 20 && fOff.config.catchUpTicks == 1, "");

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

    SimulationScheduleState sceneSched;
    sceneSched.thermalAccumDt = 0.5f;
    sceneSched.chemistryAccumDt = 0.4f;
    sceneSched.phaseAccumDt = 0.3f;
    sceneSched.thermalCounter = 3;
    sceneSched.chemistryCounter = 2;
    sceneSched.phaseCounter = 4;
    sceneSched.thermalInterval = 2;
    sceneSched.reset();
    emit("scene_reset_clears_accum",
        sceneSched.thermalAccumDt == 0.0f && sceneSched.chemistryAccumDt == 0.0f
            && sceneSched.phaseAccumDt == 0.0f
            && sceneSched.thermalCounter == 0 && sceneSched.chemistryCounter == 0
            && sceneSched.phaseCounter == 0
            && sceneSched.thermalInterval == 2, "");

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

    SimulationScheduleState phaseSched;
    phaseSched.setIntervals(1, 1, 2);
    float phaseSum = 0.0f;
    int phaseFires = 0;
    for (int n = 0; n < 4; ++n) {
        phaseSched.beginPhysicsTick(PHYSICS_DT);
        float dtPhase = 0.0f;
        if (phaseSched.takePhase(dtPhase)) {
            phaseSum += dtPhase;
            ++phaseFires;
        }
    }
    float expect = 4.0f * PHYSICS_DT;
    emit("phase_scheduler_interval_2_accum",
        phaseFires == 2 && std::abs(phaseSum - expect) < 1.0e-6f,
        "fires=" + std::to_string(phaseFires) + " sum=" + std::to_string(phaseSum));

    SimulationScheduleState phase1;
    phase1.setIntervals(1, 1, 1);
    float p1sum = 0.0f;
    int p1fires = 0;
    for (int n = 0; n < 4; ++n) {
        phase1.beginPhysicsTick(PHYSICS_DT);
        float dtPhase = 0.0f;
        if (phase1.takePhase(dtPhase)) {
            p1sum += dtPhase;
            ++p1fires;
        }
    }
    emit("phase_scheduler_interval_1_accum",
        p1fires == 4 && std::abs(p1sum - expect) < 1.0e-6f
            && std::abs(p1sum - phaseSum) < 1.0e-6f,
        "fires=" + std::to_string(p1fires) + " sum=" + std::to_string(p1sum));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
