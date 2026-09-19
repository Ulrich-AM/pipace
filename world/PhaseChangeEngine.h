#pragma once

#include "gas/GasTypes.h"
#include "substance/SubstanceTypes.h"

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;
struct WaterPhaseTickStats;

// Live Liquid <-> Gas and Solid <-> Liquid from SubstanceId + PhaseProperties.
// WaterPhaseChange is a diagnostics/compatibility wrapper only.

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

struct SolidLiquidPhaseTickStats {
    double massSolidifiedKg = 0.0;
    double massMeltedKg = 0.0;
    double latentFusionReleasedJ = 0.0;
    double latentFusionAbsorbedJ = 0.0;
    double liquidFillRemoved = 0.0;
    double liquidFillAdded = 0.0;
    int solidifyCells = 0;
    int meltPixels = 0;
    int solidPixelsSpawned = 0;
    int blockedSolidify = 0;
    int blockedMelt = 0;
    int blockedSolidifyNoDest = 0;
    int blockedSolidifyEnergy = 0;
    int blockedMeltNoDest = 0;
    int blockedMeltEnergy = 0;
};

struct PhaseChangeTickStats {
    LiquidGasPhaseTickStats liquidGas;
    SolidLiquidPhaseTickStats solidLiquid;
    double worldEditVaporDelta = 0.0;
};

// Eligible only when both liquid and gas endpoints are supported, properties
// are valid, and saturation helpers can run. Does not infer from names.
bool supportsLiveLiquidGasTransition(SubstanceId id);

// Eligible only when both solid and liquid endpoints are supported, mechanical
// and fluid properties are valid, a rigid MaterialId exists, melting/fusion
// data are finite, and mass helpers can run. Does not infer from names.
bool supportsLiveSolidLiquidTransition(SubstanceId id);

// One engine-level sandbox transfer-rate cap (fill/s). Not a per-substance
// kinetic coefficient. Equilibrium and energy still limit actual transfer.
float liquidGasTransferRateFillPerSec();

// One engine-level solid-pixel formation cap (pixels/s). Converted through
// each substance's full solid-pixel mass. Not per-substance kinetics.
float solidLiquidTransferRatePixelsPerSec();

// True if a stored gas-component view contains a live Liquid/Gas volatile.
// Inspects GasComponentView identity only — never Water vapor counters.
bool gasViewHasLiveLiquidGasVapor(GasComponentView const &view);

LiquidGasPhaseTickStats stepLiquidGasPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt, uint32_t salt);

SolidLiquidPhaseTickStats stepSolidLiquidPhaseChange(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt, uint32_t salt);

// Authoritative live phase-change step: generic Liquid <-> Gas plus
// generic Solid <-> Liquid. Does not call stepWaterPhaseChange.
PhaseChangeTickStats stepPhaseChanges(FluidEngine &fluid, RigidBodyEngine &rigid,
    GasEngine &gas, ThermalEngine &thermal, float dt);

void applyLiquidGasStatsToWaterCompat(LiquidGasPhaseTickStats const &lg, WaterPhaseTickStats &st);
void applySolidLiquidStatsToWaterCompat(SolidLiquidPhaseTickStats const &sl, WaterPhaseTickStats &st);
