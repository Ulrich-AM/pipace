#pragma once

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

struct WaterPhaseTickStats {
    double massBoiledKg = 0.0;
    double massCondensedKg = 0.0;
    double massFrozenKg = 0.0;
    double massMeltedKg = 0.0;
    double latentAbsorbedJ = 0.0;
    double latentReleasedJ = 0.0;
    double latentFusionAbsorbedJ = 0.0;
    double latentFusionReleasedJ = 0.0;
    double vaporPlacedAmount = 0.0;
    double liquidFillRemoved = 0.0;
    double liquidFillAdded = 0.0;
    int boilCells = 0;
    int condenseCells = 0;
    int freezeCells = 0;
    int meltPixels = 0;
    int icePixelsSpawned = 0;
    int blockedBoil = 0;
    int blockedCondense = 0;
    int blockedFreeze = 0;
    int blockedMelt = 0;
    int blockedBoilNoGasSpace = 0;
    int blockedBoilEquilibrium = 0;
    int blockedBoilEnergy = 0;
    int blockedBoilSafetyLimit = 0;
    int blockedCondenseNoLiquidSpace = 0;
    int blockedCondenseEquilibrium = 0;
    int blockedCondenseEnergy = 0;
    int blockedMeltNoDest = 0;
    int blockedMeltEnergy = 0;
    int blockedFreezeNoDest = 0;
    double worldEditVaporDelta = 0.0;
};

void runWaterPhaseStabilityDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal);

// Live WATER diagnostics / compatibility wrappers. Authoritative Liquid ⇄ Gas
// and Solid ⇄ Liquid physics live in PhaseChangeEngine. Honey mixtures are
// skipped. Solid water uses rigid MATERIAL_WATER_SOLID (still SUBSTANCE_WATER).
WaterPhaseTickStats stepWaterPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt);

void runWaterPhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal);

void runWaterSolidPhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal);

// Extended empirical validation: tick-rate sensitivity, incremental phase cost,
// stress cases, and storage-integrity audits. Writes misc/water_phase_validation.tsv.
void runWaterPhaseValidation(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal);
