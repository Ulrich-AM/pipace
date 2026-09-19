#pragma once

#include "fluid/FluidConfig.h"
#include "gas/GasConfig.h"
#include "thermal/ThermalConfig.h"

#include <algorithm>
#include <string>

struct FluidEngine;
struct GasEngine;
struct ThermalEngine;
struct RigidBodyEngine;

// Quality knobs only. Material constants, reaction ΔH, stoichiometry, gravity,
// and ignition gates are not part of a profile.
//
// Simulation fidelity is independent of simulation rate. physicsHz / catchUpTicks
// stay on FluidConfig as user Settings. PHYSICS_DT remains 1/30 even when the
// menu is 20 Hz (known limitation; not repaired here).

struct FluidQualitySettings {
    int maxPressureIterations = 24;
    int maxSubsteps = 6;
    int maxLimiterPasses = 16;
    bool vorticityEnabled = false;
    bool surfaceTensionEnabled = true;
    bool sprayEnabled = true;
    VelocityAdvection velocityAdvection = VelocityAdvection::SemiLagrangian;
};

struct GasQualitySettings {
    GasSimMode simMode = GasSimMode::Full; // never Off; Off is a user toggle
    int maxSubsteps = 4;
    int sleepQuietTicks = 22;
    float sleepPressureDelta = 0.004f;
};

struct RigidQualitySettings {
    int velocityContactIters = 8;
    int positionalIters = 1;
    int sleepQuietTicks = 18;
    float sleepLin = 0.16f;
};

struct ThermalQualitySettings {
    int intervalTicks = 1;
    float sleepTempEps = 0.08f;
    float sleepAmbientEps = 0.25f;
    int sleepQuietTicks = 12;
};

struct ChemistryQualitySettings {
    int intervalTicks = 1;
};

struct PhaseQualitySettings {
    int intervalTicks = 1;
};

// 0 = Performance, 1 = Balanced, 2 = Accurate. Default all Balanced.
struct SimulationQualityLevels {
    int fluid = 1;
    int gas = 1;
    int rigid = 1;
    int thermal = 1;
    int chemistry = 1;
    int phase = 1;

    void clampAll() {
        auto clampLevel = [](int v) { return std::max(0, std::min(2, v)); };
        fluid = clampLevel(fluid);
        gas = clampLevel(gas);
        rigid = clampLevel(rigid);
        thermal = clampLevel(thermal);
        chemistry = clampLevel(chemistry);
        phase = clampLevel(phase);
    }
};

struct SimulationQualityProfile {
    int level = 1; // 0 Performance, 1 Balanced, 2 Accurate; -1 Custom mixed
    bool custom = false;
    SimulationQualityLevels levels;
    FluidQualitySettings fluid;
    GasQualitySettings gas;
    RigidQualitySettings rigid;
    ThermalQualitySettings thermal;
    ChemistryQualitySettings chemistry;
    PhaseQualitySettings phase;
};

struct SimulationScheduleState {
    int thermalInterval = 1;
    int chemistryInterval = 1;
    int phaseInterval = 1;
    int thermalCounter = 0;
    int chemistryCounter = 0;
    int phaseCounter = 0;
    float thermalAccumDt = 0.0f;
    float chemistryAccumDt = 0.0f;
    float phaseAccumDt = 0.0f;

    void reset() {
        thermalCounter = chemistryCounter = phaseCounter = 0;
        thermalAccumDt = chemistryAccumDt = phaseAccumDt = 0.0f;
    }

    void setIntervals(int thermal, int chemistry, int phase) {
        thermalInterval = std::max(1, thermal);
        chemistryInterval = std::max(1, chemistry);
        phaseInterval = std::max(1, phase);
        float tCap = PHYSICS_DT * static_cast<float>(thermalInterval);
        float cCap = PHYSICS_DT * static_cast<float>(chemistryInterval);
        float pCap = PHYSICS_DT * static_cast<float>(phaseInterval);
        if (thermalAccumDt > tCap) thermalAccumDt = tCap;
        if (chemistryAccumDt > cCap) chemistryAccumDt = cCap;
        if (phaseAccumDt > pCap) phaseAccumDt = pCap;
        if (thermalCounter >= thermalInterval) thermalCounter = 0;
        if (chemistryCounter >= chemistryInterval) chemistryCounter = 0;
        if (phaseCounter >= phaseInterval) phaseCounter = 0;
    }

    void beginPhysicsTick(float dt) {
        thermalAccumDt += dt;
        chemistryAccumDt += dt;
        phaseAccumDt += dt;
        ++thermalCounter;
        ++chemistryCounter;
        ++phaseCounter;
    }

    bool takeDue(int &counter, int interval, float &accum, float &outDt) {
        if (counter < interval) return false;
        counter = 0;
        outDt = accum;
        accum = 0.0f;
        return outDt > 0.0f;
    }

    bool takeThermal(float &outDt) {
        return takeDue(thermalCounter, thermalInterval, thermalAccumDt, outDt);
    }
    bool takeChemistry(float &outDt) {
        return takeDue(chemistryCounter, chemistryInterval, chemistryAccumDt, outDt);
    }
    bool takePhase(float &outDt) {
        return takeDue(phaseCounter, phaseInterval, phaseAccumDt, outDt);
    }
};

// 0 = Performance, 1 = Balanced, 2 = Accurate. Auto selects 0 or 1 only.
SimulationQualityProfile profileForQualityLevel(int level);
SimulationQualityProfile profileForCustomLevels(SimulationQualityLevels levels);

void applySimulationQuality(SimulationQualityProfile const &profile,
    FluidEngine &fluid, GasEngine &gas, ThermalEngine &thermal,
    RigidBodyEngine &rigid, SimulationScheduleState &schedule);

std::string describeSimulationQuality(SimulationQualityProfile const &profile,
    SimulationScheduleState const &schedule, GasSimMode gasMode, bool thermalEnabled);

void runSimulationQualitySanity();
