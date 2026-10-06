#pragma once

#include "substance/RuntimeSubstance.h"

// Phase 19A: read-only engine property bridge for RuntimeSubstanceRef.
//
// Built-ins delegate to the existing SubstanceDefinition tables.
// Generated refs resolve through GeneratedMaterialRegistry and evaluate the
// compiled SACE models at the requested temperature.
//
// This is a property-query layer only. It does not change liquid/gas component
// storage, does not spawn generated matter, and does not allocate SubstanceId.

struct RuntimeLiquidPropertySample {
    bool valid = false;
    float temperatureK = 0.0f;
    float densityRelativeToWater = 0.0f;
    float specificHeatJPerKgK = 0.0f;
    float thermalConductivityWPerMK = 0.0f;
    float solverViscosity = 0.0f;
    float solverSurfaceTension = 0.0f;
};

struct RuntimeGasPropertySample {
    bool valid = false;
    float temperatureK = 0.0f;
    float molarMassGPerMol = 0.0f;
    float specificHeatJPerKgK = 0.0f;
    float thermalConductivityWPerMK = 0.0f;
};

// Samples phase-appropriate properties at temperatureK.
// Invalid/stale refs, unsupported phases, invalid temperatures, and generated
// model-domain failures return false with an all-zero output.
// Generated queries never borrow Water/Air material properties as fallbacks.
bool sampleRuntimeLiquidProperties(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeLiquidPropertySample &out);

bool sampleRuntimeGasProperties(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeGasPropertySample &out);

// Headless: --runtime-property-diag -> misc/runtime_property_diag.tsv
void runRuntimeSubstancePropertyDiagnostics();
