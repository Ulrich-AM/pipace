#pragma once

#include "chemistry/SaceCatalog.h"
#include "chemistry/SaceEstimation.h"
#include "chemistry/SaceSimulationReadiness.h"

// SACE Phase 17: generated simulation-property compiler.
// Translates coherent SACE physical models into an engine-facing profile.
// Not identity. Not SubstanceId. Not world storage. Not spawnable matter.
// Do not call from physics ticks. Compilation is record-prep work, not a
// per-cell query. Direction: identity -> physical models -> readiness ->
// compilation. Never copy solver scalars back into SaceGeneratedProperties.
//
// profile.valid means at least one coherent subprofile compiled. Gas, liquid,
// and liquid-gas equilibrium compile independently when their models exist.
// Missing data fails that subprofile; Water/Air properties are never borrowed
// as material fallbacks. Solver viscosity/surface-tension use Water-relative
// unit-calibration anchors only.

constexpr float kSaceCompilerReferenceTemperatureK = 298.15f;

// Physical Water anchors for sandbox unit bridges. Not SACE estimates and
// not replacements for built-in Water FluidProperties.
constexpr double kPhysicalWaterViscosityAt298KPaS = 8.90e-4;
constexpr double kPhysicalWaterSurfaceTensionAt298KNPerM = 0.07197;

struct SaceCompiledLiquidProfile {
    bool valid = false;
    float referenceTemperatureK = kSaceCompilerReferenceTemperatureK;
    float densityRelativeToWater = 0.0f;
    float specificHeatJPerKgK = 0.0f;
    float thermalConductivityWPerMK = 0.0f;
    float solverViscosity = 0.0f;
    float solverSurfaceTension = 0.0f;
    float viscArrheniusK = 0.0f;
    SaceCostaldLiquidDensityModel densityModel{};
    SaceRowlinsonPolingLiquidCpModel heatCapacityModel{};
    SaceJobackLiquidViscosityModel viscosityModel{};
    SaceSastriRaoSurfaceTensionModel surfaceTensionModel{};
    SaceSatoRiedelLiquidConductivityModel conductivityModel{};
};

struct SaceCompiledGasProfile {
    bool valid = false;
    float referenceTemperatureK = kSaceCompilerReferenceTemperatureK;
    float molarMassGPerMol = 0.0f;
    float specificHeatJPerKgK = 0.0f;
    float thermalConductivityWPerMK = 0.0f;
    SaceJobackIdealGasCpModel heatCapacityModel{};
    SaceGharagheiziGasConductivityModel conductivityModel{};
};

struct SaceCompiledPhaseProfile {
    bool liquidGasValid = false;
    float normalBoilingPointK = 0.0f;
    float criticalTemperatureK = 0.0f;
    float criticalPressurePa = 0.0f;
    float referencePressurePa = 101325.0f;
    double acentricFactor = 0.0;
    double molarMassGPerMol = 0.0;
    SaceLeeKeslerVaporModel vaporPressureModel{};
    SaceWatsonVaporizationModel vaporizationModel{};
};

struct SaceCompiledSimulationProfile {
    SaceRecordId sourceRecord = kSaceRecordNone;
    bool gasReady = false;
    bool liquidReady = false;
    bool liquidGasPhaseChangeReady = false;
    SaceCompiledLiquidProfile liquid{};
    SaceCompiledGasProfile gas{};
    SaceCompiledPhaseProfile phase{};
    bool valid = false;
};

void saceResetCompiledSimulationProfile(SaceCompiledSimulationProfile &out);

bool sacePhysicalDensityToWaterRelative(double densityKgPerM3, float &outRelative);

bool saceCalibrateSolverViscosityFromPhysicalPaS(double physicalPaS, float &outSolver);
bool saceCalibrateSolverSurfaceTensionFromPhysicalNPerM(double physicalNPerM, float &outSolver);

bool saceCompiledSolverViscosityAtTemperature(
    SaceCompiledLiquidProfile const &liquid,
    float temperatureK,
    float &outSolver);

// Rebuilds one coherent Joback path from the record graph. Does not mutate
// the record or catalog. Does not allocate a SubstanceId.
bool saceCompileSimulationProfile(
    SaceGeneratedRecord const &record,
    SaceCompiledSimulationProfile &out);

// Headless: --sace-compiler-diag -> misc/sace_compiler_diag.tsv
void runSaceCompilerDiagnostics();
