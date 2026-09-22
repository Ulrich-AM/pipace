#pragma once

#include "chemistry/SaceCatalog.h"

#include <cstdint>

// SACE Phase 13: generated simulation-readiness / compilation preflight.
// Answers whether a generated SACE record has enough coherent data to
// eventually compile into live liquid/gas simulation properties.
// Not identity. Not world simulation. Not SubstanceId allocation.
// Computed from trusted graph + coherent structural models; not stored
// on SaceGeneratedRecord. Do not call from physics ticks.

enum class SaceSimulationRequirement : uint64_t {
    MolecularGraph = 1ull << 0,
    MolarMass = 1ull << 1,
    NormalBoilingPoint = 1ull << 2,
    CriticalTemperature = 1ull << 3,
    CriticalPressure = 1ull << 4,
    AcentricFactor = 1ull << 5,
    VaporPressureModel = 1ull << 6,
    VaporizationEnthalpyModel = 1ull << 7,
    GasHeatCapacityModel = 1ull << 8,
    LiquidHeatCapacityModel = 1ull << 9,
    LiquidDensityModel = 1ull << 10,
    LiquidViscosity = 1ull << 11,
    LiquidSurfaceTension = 1ull << 12,
    LiquidThermalConductivity = 1ull << 13,
    GasThermalConductivity = 1ull << 14
};

struct SaceSimulationReadiness {
    bool gasThermoReady = false;
    bool liquidThermoReady = false;
    bool liquidGasEquilibriumReady = false;
    bool liveGasReady = false;
    bool liveLiquidReady = false;
    bool liveLiquidGasPhaseChangeReady = false;
    uint64_t missingMask = 0;
};

SaceSimulationReadiness saceAssessSimulationReadiness(SaceGeneratedRecord const &record);

inline bool saceReadinessMissing(
    SaceSimulationReadiness const &r,
    SaceSimulationRequirement requirement)
{
    return (r.missingMask & static_cast<uint64_t>(requirement)) != 0;
}

char const *saceSimulationRequirementKey(SaceSimulationRequirement requirement);

// Headless: --sace-readiness-diag -> misc/sace_readiness_diag.tsv
void runSaceReadinessDiagnostics();
