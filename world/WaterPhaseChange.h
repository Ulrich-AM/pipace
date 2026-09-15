#pragma once

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

struct WaterPhaseTickStats {
    double massBoiledKg = 0.0;
    double massCondensedKg = 0.0;
    double latentAbsorbedJ = 0.0;
    double latentReleasedJ = 0.0;
    double vaporPlacedAmount = 0.0;
    double liquidFillRemoved = 0.0;
    double liquidFillAdded = 0.0;
    int boilCells = 0;
    int condenseCells = 0;
    int blockedBoil = 0;
    int blockedCondense = 0;
};

// Live WATER liquid ⇄ gas transfer. Honey mixtures are skipped.
// Occupancy stays one primary medium: vapor is placed in neighboring gas cells.
WaterPhaseTickStats stepWaterPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt);

void runWaterPhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal);
