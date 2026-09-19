#pragma once

#include "gas/GasTypes.h"
#include "substance/SubstanceTypes.h"
#include "thermal/ThermalTypes.h"

#include <algorithm>
#include <cmath>

// Phase-transfer accounting helpers (fill ↔ mass ↔ gas amount, latent energy,
// Clausius–Clapeyron saturation). Live Liquid ⇄ Gas is in
// world/PhaseChangeEngine.cpp. Live water solid ⇄ liquid is in
// world/WaterPhaseChange.cpp.
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
//   GasEngine amount       conserved "cell-atmospheres" at Tref = AMBIENT_TEMPERATURE_K.
//                          Pressure is P_atm = (amount / volume) * (T / Tref).
//                          At ambient T this matches 1.0 amount in 1.0 volume = 1 atm.
//                          There is no PdV work; pressure is derived and does not
//                          modify thermal energy.
//                          1.0 ≈ one cell of 1 atm of the stored gas species.
//                          Generic composition: amount[] is total cell-atmospheres;
//                          per-cell SubstanceId slots (Air, Water vapor, …) sum to amount.
//                          Water vapor is SUBSTANCE_WATER + MatterPhase::Gas.
//   air mass               amount * 1.204 kg/m³ * V  (gasMassKg; ambient-T air)
//   water vapor mass       amount * ρ_vapor(id) * V
//                          ρ_vapor from ideal gas at phase.referencePressurePa and
//                          AMBIENT_TEMPERATURE_K using chemical.molarMass.
//                          Not a Steam SubstanceId.
//   thermal energy         Joules on liquidHeat / gas.heat / rigid pixels
//   latent heats           PhaseProperties (copied from thermal authoring), J/kg
//
// Solid water: SubstanceId stays WATER. mechanical.densityRel is ice (~0.917).
// Do not borrow stone. There is no SUBSTANCE_ICE.

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

// Numerical window for the integrated Clausius–Clapeyron curve.
// This is not a steam table and not substance-specific critical-point data.
// Low-T clamp keeps 1/T and log finite; high-T clamp keeps exp(Lv/R Δ(1/T)) finite.
inline double phaseEquilibriumTminK(SubstanceId id) {
    float Tm = phaseForSubstance(id).meltingPointK;
    double lo = 180.0;
    if (Tm > 1.0f) lo = std::min(lo, std::max(50.0, static_cast<double>(Tm) * 0.5));
    return lo;
}

inline double phaseEquilibriumTmaxK(SubstanceId id) {
    float Tb = phaseForSubstance(id).boilingPointK;
    double hi = (Tb > 1.0f) ? static_cast<double>(Tb) * 2.5 : 2000.0;
    return std::min(std::max(hi, 400.0), 2000.0);
}

// P_sat(T) = P0 * exp( -(Lv/R_spec) * (1/T - 1/Tb) )
// Uses PhaseProperties.referencePressurePa / boilingPointK / latentHeatVaporization
// and chemical.molarMass. Valid near the reference boiling point.
inline double saturationVaporPressurePa(SubstanceId id, double temperatureK) {
    PhaseProperties const &p = phaseForSubstance(id);
    double Tb = static_cast<double>(p.boilingPointK);
    double P0 = static_cast<double>(p.referencePressurePa);
    double Lv = static_cast<double>(p.latentHeatVaporization);
    double R = specificGasConstantJPerKgK(id);
    if (!(Tb > 1.0) || !(P0 > 0.0) || !(Lv > 0.0) || !(R > 0.0)) return 0.0;
    if (!std::isfinite(temperatureK) || temperatureK <= 1.0) return 0.0;
    double Tmin = phaseEquilibriumTminK(id);
    double Tmax = phaseEquilibriumTmaxK(id);
    double T = std::clamp(temperatureK, Tmin, Tmax);
    double expo = -(Lv / R) * (1.0 / T - 1.0 / Tb);
    expo = std::clamp(expo, -50.0, 50.0);
    double P = P0 * std::exp(expo);
    if (!std::isfinite(P) || P < 0.0) return 0.0;
    return P;
}

// Inverse: T_sat(P) from the same integrated Clausius–Clapeyron relation.
inline double saturationTemperatureK(SubstanceId id, double pressurePa) {
    PhaseProperties const &p = phaseForSubstance(id);
    double Tb = static_cast<double>(p.boilingPointK);
    double P0 = static_cast<double>(p.referencePressurePa);
    double Lv = static_cast<double>(p.latentHeatVaporization);
    double R = specificGasConstantJPerKgK(id);
    if (!(Tb > 1.0) || !(P0 > 0.0) || !(Lv > 0.0) || !(R > 0.0))
        return Tb > 1.0 ? Tb : 0.0;
    double Tmin = phaseEquilibriumTminK(id);
    double Tmax = phaseEquilibriumTmaxK(id);
    double P = pressurePa;
    if (!std::isfinite(P) || P < 1.0) P = 1.0;
    double ln = std::log(P / P0);
    ln = std::clamp(ln, -50.0, 50.0);
    double inv = 1.0 / Tb - (R / Lv) * ln;
    if (!(inv > 1.0 / Tmax)) return Tmax;
    if (!(inv < 1.0 / Tmin)) return Tmin;
    double T = 1.0 / inv;
    if (!std::isfinite(T)) return Tb;
    return std::clamp(T, Tmin, Tmax);
}

// Reference gas density at 1 atm-equivalent and ambient T.
// Air uses the tabulated sandbox density so thermal/gasMassKg stay identical.
// Other gas-capable substances with molarMass use the ideal-gas law.
inline double referenceGasDensityKgM3(SubstanceId id) {
    if (!supportsPhase(id, MatterPhase::Gas)) return 0.0;
    if (id == SUBSTANCE_AIR) return static_cast<double>(AIR_DENSITY_KG_M3);
    float M = chemicalForSubstance(id).molarMass;
    if (M > 0.0f && std::isfinite(M)) {
        float P = phaseForSubstance(id).referencePressurePa;
        if (!(P > 0.0f) || !std::isfinite(P)) P = GAS_REFERENCE_PRESSURE_PA;
        double T = static_cast<double>(AMBIENT_TEMPERATURE_K);
        if (!(T > 1.0)) T = 293.15;
        double Mkg = static_cast<double>(M) * 0.001;
        double rho = static_cast<double>(P) * Mkg / (static_cast<double>(UNIVERSAL_GAS_R_J_MOL_K) * T);
        if (std::isfinite(rho) && rho > 0.0) return rho;
    }
    // Numerical safety only. Does not relabel the gas as Air.
    return static_cast<double>(AIR_DENSITY_KG_M3);
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

inline bool substanceHasMolarMass(SubstanceId id) {
    ChemicalProperties const &c = chemicalForSubstance(id);
    return c.valid && c.molarMass > 0.0f && std::isfinite(c.molarMass);
}

inline double massKgToMoles(SubstanceId id, double massKg) {
    if (!substanceHasMolarMass(id) || !(massKg > 0.0) || !std::isfinite(massKg)) return 0.0;
    double kgPerMol = static_cast<double>(chemicalForSubstance(id).molarMass) * 1.0e-3;
    if (!(kgPerMol > 0.0)) return 0.0;
    double n = massKg / kgPerMol;
    return std::isfinite(n) ? n : 0.0;
}

inline double molesToMassKg(SubstanceId id, double moles) {
    if (!substanceHasMolarMass(id) || !(moles > 0.0) || !std::isfinite(moles)) return 0.0;
    double kgPerMol = static_cast<double>(chemicalForSubstance(id).molarMass) * 1.0e-3;
    if (!(kgPerMol > 0.0)) return 0.0;
    double m = moles * kgPerMol;
    return std::isfinite(m) ? m : 0.0;
}

// Storage amount (liquid fill, gas cell-atmospheres, or solid fraction) -> moles.
// Returns 0 if molar mass or the phase mass path is unavailable. Does not invent M.
inline double storageAmountToMoles(SubstanceId id, MatterPhase phase, double amount, double cellsPerMeter) {
    return massKgToMoles(id, amountToMassKg(id, phase, amount, cellsPerMeter));
}

inline double molesToStorageAmount(SubstanceId id, MatterPhase phase, double moles, double cellsPerMeter) {
    return massKgToAmount(id, phase, molesToMassKg(id, moles), cellsPerMeter);
}

// Shared gas heat capacity: Σ mass_i · Cp_gas,i. Pure Air uses gasMassKg × air Cp.
// Missing molar/density data uses Air's mass scale as a numerical fallback only.
inline float gasMixtureThermalCapacity(GasComponentView const &view, float cellsPerMeter = 4.0f) {
    if (view.count <= 0) return 0.0f;
    if (view.count == 1 && view.items[0].id == SUBSTANCE_AIR)
        return thermalCapacity(gasMassKg(view.items[0].amount, cellsPerMeter),
            thermalForSubstance(SUBSTANCE_AIR).specificHeat);
    float cap = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId id = view.items[n].id;
        float amt = view.items[n].amount;
        if (!(amt > GAS_MIN_AMOUNT) || !validGasComponentId(id)) continue;
        float mass = 0.0f;
        if (id == SUBSTANCE_AIR)
            mass = gasMassKg(amt, cellsPerMeter);
        else {
            mass = static_cast<float>(gasAmountToMassKg(id, amt, static_cast<double>(cellsPerMeter)));
            if (!(mass > 0.0f) || !std::isfinite(mass))
                mass = gasMassKg(amt, cellsPerMeter);
        }
        cap += thermalCapacity(mass, gasPhaseSpecificHeat(id));
    }
    return cap;
}

// Amount-fraction mixture of referenceGasDensityKgM3. Pure-component fast path.
// Missing/invalid density uses Air's tabulated rho as a numerical fallback only;
// it does not relabel the substance as Air.
inline float gasMixtureReferenceDensityKgM3(GasComponentView const &view) {
    float airRho = AIR_DENSITY_KG_M3;
    if (!(airRho > 0.0f) || !std::isfinite(airRho)) airRho = 1.204f;
    auto rhoOf = [&](SubstanceId id) {
        double rho = referenceGasDensityKgM3(id);
        if (rho > 0.0 && std::isfinite(rho)) return static_cast<float>(rho);
        return airRho;
    };
    if (view.count <= 0) return airRho;
    if (view.count == 1) return rhoOf(view.items[0].id);
    float tot = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        if (view.items[n].amount > GAS_MIN_AMOUNT && validGasComponentId(view.items[n].id))
            tot += view.items[n].amount;
    }
    if (!(tot > GAS_MIN_AMOUNT)) return airRho;
    float rho = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId id = view.items[n].id;
        float amt = view.items[n].amount;
        if (!(amt > GAS_MIN_AMOUNT) || !validGasComponentId(id)) continue;
        rho += (amt / tot) * rhoOf(id);
    }
    if (!(rho > 0.0f) || !std::isfinite(rho)) return airRho;
    return rho;
}

// Amount-fraction mixture of thermal.conductivity. Pure-component fast path.
// Missing/invalid k uses Air conductivity as a numerical fallback only.
inline float gasMixtureConductivity(GasComponentView const &view) {
    float airK = thermalForSubstance(SUBSTANCE_AIR).conductivity;
    if (!(airK > 0.0f) || !std::isfinite(airK)) airK = 0.026f;
    auto kOf = [&](SubstanceId id) {
        float k = thermalForSubstance(id).conductivity;
        if (k > 0.0f && std::isfinite(k)) return k;
        return airK;
    };
    if (view.count <= 0) return airK;
    if (view.count == 1) return kOf(view.items[0].id);
    float tot = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        if (view.items[n].amount > GAS_MIN_AMOUNT && validGasComponentId(view.items[n].id))
            tot += view.items[n].amount;
    }
    if (!(tot > GAS_MIN_AMOUNT)) return airK;
    float k = 0.0f;
    for (int n = 0; n < view.count; ++n) {
        SubstanceId id = view.items[n].id;
        float amt = view.items[n].amount;
        if (!(amt > GAS_MIN_AMOUNT) || !validGasComponentId(id)) continue;
        k += (amt / tot) * kOf(id);
    }
    if (!(k > 0.0f) || !std::isfinite(k)) return airK;
    return k;
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
