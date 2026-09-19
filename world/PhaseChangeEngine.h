#pragma once

#include "substance/SubstanceTypes.h"

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;
struct WaterPhaseTickStats;

// Live Liquid <-> Gas from SubstanceId + PhaseProperties.
// Solid <-> Liquid remains Water-specific in WaterPhaseChange.cpp.

struct LiquidGasPhaseTickStats {
    double massVaporizedKg = 0.0;
    double massCondensedKg = 0.0;
    double latentVaporizationAbsorbedJ = 0.0;
    double latentCondensationReleasedJ = 0.0;
    double vaporPlacedAmount = 0.0;
    double liquidFillRemoved = 0.0;
    double liquidFillAdded = 0.0;
    int vaporizedCells = 0;
    int condensedCells = 0;
    int blockedNoGasSpace = 0;
    int blockedNoLiquidSpace = 0;
    int blockedVaporizationEquilibrium = 0;
    int blockedCondensationEquilibrium = 0;
    int blockedVaporizationEnergy = 0;
    int blockedCondensationEnergy = 0;
    int blockedSafetyLimit = 0;
};

struct PhaseChangeTickStats {
    LiquidGasPhaseTickStats liquidGas;
};

// Eligible only when both liquid and gas endpoints are supported, properties
// are valid, and saturation helpers can run. Does not infer from names.
bool supportsLiveLiquidGasTransition(SubstanceId id);

// One engine-level sandbox transfer-rate cap (fill/s). Not a per-substance
// kinetic coefficient. Equilibrium and energy still limit actual transfer.
float liquidGasTransferRateFillPerSec();

LiquidGasPhaseTickStats stepLiquidGasPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt, uint32_t salt);

// Generic Liquid <-> Gas plus the existing Water Solid <-> Liquid path.
PhaseChangeTickStats stepPhaseChanges(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt);

void applyLiquidGasStatsToWaterCompat(LiquidGasPhaseTickStats const &lg, WaterPhaseTickStats &st);
