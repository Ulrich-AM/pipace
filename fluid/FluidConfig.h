#pragma once

#include "FluidTypes.h"

inline LiquidProperties kHoneyLiquid() {
    LiquidProperties p{1.42f, 0.35f, 0.080f};
    p.thermal = kHoneyThermal();
    p.viscArrheniusK = 2500.0f;
    return p;
}

struct FluidConfig {
    // Compatibility (Prompt 1): these are the current liquid-phase tables for
    // SUBSTANCE_WATER / SUBSTANCE_HONEY. Identity is SubstanceId.
    LiquidProperties water{1.0f, 0.006f, 0.055f};
    LiquidProperties honey = kHoneyLiquid();
    float cellsPerMeter = 4.0f;
    float gravityMetersPerSecondSquared = 9.81f;
    float maxVelocity = 75.0f;             // prevents pathological CFL spikes
    float maxTravelPerSubstep = 1.0f;      // CFL-like transport limit
    float vorticityStrength = 0.035f;
    int maxPressureIterations = 24;
    int maxSubsteps = 6;
    int maxLimiterPasses = 16;
    bool vorticityEnabled = false;
    bool surfaceTensionEnabled = true;
    bool sprayEnabled = true;
    VelocityAdvection velocityAdvection = VelocityAdvection::SemiLagrangian;
    QualityPreset quality = QualityPreset::Medium;
    int autoQualityLevel = 1; // 0 = Low knobs, 1 = Medium (Auto never climbs to High)
    int physicsHz = 30;
    int catchUpTicks = 2;
    // 0 = Auto (conservative hardware pick). 1 = deterministic sequential path.
    // 2/4/6/8 = requested total workers (capped to the machine).
    int workerCount = 0;
    // Skip the pool when a red or black pressure list is smaller than this.
    int pressureParallelMinCells = 20000;
    // When true, grid edges are closed walls. When false they are a void sink.
    bool walledBorders = false;
};

inline void applyFluidQualityKnobs(FluidConfig &c, int level) {
    if (level <= 0) {
        c.maxPressureIterations = 8;
        c.maxSubsteps = 2;
        c.maxLimiterPasses = 4;
        c.vorticityEnabled = false;
        c.surfaceTensionEnabled = false;
        c.sprayEnabled = false;
        c.velocityAdvection = VelocityAdvection::SemiLagrangian;
        c.physicsHz = 20;
        c.catchUpTicks = 1;
    } else if (level == 1) {
        c.maxPressureIterations = 24;
        c.maxSubsteps = 6;
        c.maxLimiterPasses = 16;
        c.surfaceTensionEnabled = true;
        c.sprayEnabled = true;
        c.physicsHz = 30;
        c.catchUpTicks = 2;
    } else {
        c.maxPressureIterations = 30;
        c.maxSubsteps = 6;
        c.maxLimiterPasses = 16;
        c.surfaceTensionEnabled = true;
        c.sprayEnabled = true;
        c.physicsHz = 30;
        c.catchUpTicks = 2;
    }
}
