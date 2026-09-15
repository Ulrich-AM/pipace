#pragma once

#include "gas/GasTypes.h"
#include "substance/SubstanceTypes.h"
#include "thermal/ThermalTypes.h"

#include <cmath>

// Phase-transfer accounting. Conversion is mathematical only — no engine writes,
// no temperature-triggered boiling/melting/freezing/condensation.
//
// Current sandbox units (do not invent a second system):
//
//   cellsPerMeter          default 4 → cell edge dx = 0.25 m
//   implied 2D depth       one cell, so V = dx³ = 0.015625 m³
//   FluidEngine fill       fractional occupancy of that cell volume (typically 0..1;
//                          overfill is allowed by the solver but mass helpers treat
//                          fill as a linear scale). fill = 1 is one full liquid cell.
//   liquid density         FluidProperties.density, relative to water = 1.0
//   liquid mass            fill * densityRel * 1000 kg/m³ * V
//                          full water cell = 1 * 1 * 1000 * 0.015625 = 15.625 kg
//   GasEngine amount       conserved "cell-atmospheres". 1.0 amount in 1.0 available
//                          volume = 1 atm (isothermal P = amount / volume).
//                          1.0 ≈ one cell of 1 atm of the stored gas species.
//                          Today the engine stores only Air (GasSpecies::Air).
//   air mass               amount * 1.204 kg/m³ * V  (gasMassKg; ambient-T air)
//   water vapor mass       amount * ρ_vapor(id) * V
//                          ρ_vapor from ideal gas at phase.referencePressurePa and
//                          AMBIENT_TEMPERATURE_K using chemical.molarMass.
//                          Not a Steam SubstanceId. Not yet a GasEngine species.
//   thermal energy         Joules on liquidHeat / gas.heat / rigid pixels
//   latent heats           PhaseProperties (copied from thermal authoring), J/kg
//
// Solid water: SubstanceId stays WATER. There is no mechanical.densityRel for ice.
// Solid mass helpers return 0 until a water-specific solid table exists.
// Do not borrow stone.

struct PhaseTransferResult {
    SubstanceId substance = SUBSTANCE_NONE;
    MatterPhase from = MatterPhase::None;
    MatterPhase to = MatterPhase::None;
    double massTransferred = 0.0;           // kg
    double energyTransferred = 0.0;         // J; positive = latent heat absorbed
    double sourceAmountRemoved = 0.0;       // fill, gas amount, or solid fraction
    double destinationAmountAdded = 0.0;
    bool success = false;
};

inline double sandboxCellLengthM(double cellsPerMeter) {
    double cpm = cellsPerMeter > 0.1 ? cellsPerMeter : 4.0;
    return 1.0 / cpm;
}

inline double sandboxCellVolumeM3(double cellsPerMeter) {
    double dx = sandboxCellLengthM(cellsPerMeter);
    return dx * dx * dx;
}

// Metadata only: both endpoints are capability-supported. Not "transition now".
inline bool canTransition(SubstanceId id, MatterPhase from, MatterPhase to) {
    if (id == SUBSTANCE_NONE || !validSubstance(id)) return false;
    if (from == to) return false;
    if (from == MatterPhase::None || to == MatterPhase::None) return false;
    if (from == MatterPhase::Plasma || to == MatterPhase::Plasma) return false;
    return supportsPhase(id, from) && supportsPhase(id, to);
}

inline double fusionEnergyJ(SubstanceId id, double massKg) {
    if (!(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    return massKg * static_cast<double>(phaseForSubstance(id).latentHeatFusion);
}

inline double vaporizationEnergyJ(SubstanceId id, double massKg) {
    if (!(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    return massKg * static_cast<double>(phaseForSubstance(id).latentHeatVaporization);
}

// J/(kg·K). 0 if molar mass is unknown.
inline double specificGasConstantJPerKgK(SubstanceId id) {
    float M = chemicalForSubstance(id).molarMass;
    if (!(M > 0.0f) || !std::isfinite(M)) return 0.0;
    return static_cast<double>(UNIVERSAL_GAS_R_J_MOL_K) * 1000.0 / static_cast<double>(M);
}

// Reference gas density at 1 atm-equivalent and ambient T.
// Air uses the tabulated sandbox density so thermal/gasMassKg stay identical.
// Other gas-capable substances with molarMass use the ideal-gas law.
inline double referenceGasDensityKgM3(SubstanceId id) {
    if (!supportsPhase(id, MatterPhase::Gas)) return 0.0;
    if (id == SUBSTANCE_AIR) return static_cast<double>(AIR_DENSITY_KG_M3);
    float M = chemicalForSubstance(id).molarMass;
    if (!(M > 0.0f) || !std::isfinite(M)) return 0.0;
    float P = phaseForSubstance(id).referencePressurePa;
    if (!(P > 0.0f) || !std::isfinite(P)) P = GAS_REFERENCE_PRESSURE_PA;
    double T = static_cast<double>(AMBIENT_TEMPERATURE_K);
    if (!(T > 1.0)) T = 293.15;
    double Mkg = static_cast<double>(M) * 0.001;
    double rho = static_cast<double>(P) * Mkg / (static_cast<double>(UNIVERSAL_GAS_R_J_MOL_K) * T);
    if (!std::isfinite(rho) || rho <= 0.0) return 0.0;
    return rho;
}

inline double liquidFillToMassKg(SubstanceId id, double fill, double cellsPerMeter) {
    if (!hasFluidProperties(id) || !(fill > 0.0) || !std::isfinite(fill)) return 0.0;
    double rhoRel = static_cast<double>(fluidForSubstance(id).density);
    if (!(rhoRel > 0.0) || !std::isfinite(rhoRel)) return 0.0;
    return rhoRel * static_cast<double>(WATER_DENSITY_KG_M3) * sandboxCellVolumeM3(cellsPerMeter) * fill;
}

inline double massKgToLiquidFill(SubstanceId id, double massKg, double cellsPerMeter) {
    double one = liquidFillToMassKg(id, 1.0, cellsPerMeter);
    if (!(one > 0.0) || !(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    double fill = massKg / one;
    return std::isfinite(fill) ? fill : 0.0;
}

inline double gasAmountToMassKg(SubstanceId id, double amount, double cellsPerMeter) {
    if (!(amount > 0.0) || !std::isfinite(amount)) return 0.0;
    double rho = referenceGasDensityKgM3(id);
    if (!(rho > 0.0)) return 0.0;
    return amount * rho * sandboxCellVolumeM3(cellsPerMeter);
}

inline double massKgToGasAmount(SubstanceId id, double massKg, double cellsPerMeter) {
    double one = gasAmountToMassKg(id, 1.0, cellsPerMeter);
    if (!(one > 0.0) || !(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    double amount = massKg / one;
    return std::isfinite(amount) ? amount : 0.0;
}

inline double solidFractionToMassKg(SubstanceId id, double fraction, double cellsPerMeter) {
    if (!hasMechanicalProperties(id) || !(fraction > 0.0) || !std::isfinite(fraction)) return 0.0;
    double rhoRel = static_cast<double>(mechanicalForSubstance(id).densityRel);
    if (!(rhoRel > 0.0) || !std::isfinite(rhoRel)) return 0.0;
    return rhoRel * static_cast<double>(WATER_DENSITY_KG_M3) * sandboxCellVolumeM3(cellsPerMeter) * fraction;
}

inline double massKgToSolidFraction(SubstanceId id, double massKg, double cellsPerMeter) {
    double one = solidFractionToMassKg(id, 1.0, cellsPerMeter);
    if (!(one > 0.0) || !(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    double f = massKg / one;
    return std::isfinite(f) ? f : 0.0;
}

inline double latentEnergyForTransitionJ(SubstanceId id, MatterPhase from, MatterPhase to, double massKg) {
    if (!(massKg > 0.0)) return 0.0;
    double Lf = fusionEnergyJ(id, massKg);
    double Lv = vaporizationEnergyJ(id, massKg);
    if (from == MatterPhase::Solid && to == MatterPhase::Liquid) return Lf;
    if (from == MatterPhase::Liquid && to == MatterPhase::Solid) return -Lf;
    if (from == MatterPhase::Liquid && to == MatterPhase::Gas) return Lv;
    if (from == MatterPhase::Gas && to == MatterPhase::Liquid) return -Lv;
    if (from == MatterPhase::Solid && to == MatterPhase::Gas) return Lf + Lv;
    if (from == MatterPhase::Gas && to == MatterPhase::Solid) return -(Lf + Lv);
    return 0.0;
}

inline double amountToMassKg(SubstanceId id, MatterPhase phase, double amount, double cellsPerMeter) {
    switch (phase) {
        case MatterPhase::Liquid: return liquidFillToMassKg(id, amount, cellsPerMeter);
        case MatterPhase::Gas: return gasAmountToMassKg(id, amount, cellsPerMeter);
        case MatterPhase::Solid: return solidFractionToMassKg(id, amount, cellsPerMeter);
        default: return 0.0;
    }
}

inline double massKgToAmount(SubstanceId id, MatterPhase phase, double massKg, double cellsPerMeter) {
    switch (phase) {
        case MatterPhase::Liquid: return massKgToLiquidFill(id, massKg, cellsPerMeter);
        case MatterPhase::Gas: return massKgToGasAmount(id, massKg, cellsPerMeter);
        case MatterPhase::Solid: return massKgToSolidFraction(id, massKg, cellsPerMeter);
        default: return 0.0;
    }
}

// Closed accounting only. Does not mutate FluidEngine / GasEngine / rigid storage.
inline PhaseTransferResult convertPhaseAmount(SubstanceId id, MatterPhase from, MatterPhase to,
    double sourceAmount, double cellsPerMeter)
{
    PhaseTransferResult r;
    r.substance = id;
    r.from = from;
    r.to = to;
    if (!canTransition(id, from, to)) return r;
    if (!(sourceAmount > 0.0) || !std::isfinite(sourceAmount)) return r;
    double mass = amountToMassKg(id, from, sourceAmount, cellsPerMeter);
    if (!(mass > 0.0) || !std::isfinite(mass)) return r;
    double dest = massKgToAmount(id, to, mass, cellsPerMeter);
    if (!(dest > 0.0) || !std::isfinite(dest)) return r;
    r.massTransferred = mass;
    r.energyTransferred = latentEnergyForTransitionJ(id, from, to, mass);
    r.sourceAmountRemoved = sourceAmount;
    r.destinationAmountAdded = dest;
    r.success = true;
    return r;
}

void runPhaseTransferDiagnostics();
