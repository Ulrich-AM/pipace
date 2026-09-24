#pragma once

#include "chemistry/SaceMolecule.h"
#include "chemistry/SaceProperties.h"
#include "chemistry/SaceDescriptors.h"

#include <cstdint>

// SACE Phase 7–16: Joback-Reid 1987 subset + Lee-Kesler Psat + Watson Hvap + Cp + COSTALD + viscosity + Sastri-Rao sigma + thermal conductivity.
// Phase 7: Tb[K] = 198.2 + SUM(Tb groups).
// Phase 8: Tc, Pc, Vc from the same six groups and the Joback Tb (not stored Tb).
// Phase 9: Lee-Kesler omega and Psat(T) from the same Joback Tb/Tc/Pc tuple.
// Phase 10: Joback ΔHvap(Tb) [J/mol] + Watson Hvap(T). Not live phase physics.
// Phase 11: Joback ideal-gas Cp(T) + Rowlinson-Poling liquid Cp(T). 298–1000 K.
// Phase 12: COSTALD saturated-liquid Vs(T)/rho(T). Joback Vc as V* fallback.
// Phase 14: Joback liquid dynamic viscosity mu = MW[g/mol] * exp(A/T + B) in Pa*s.
// Phase 15: Sastri-Rao liquid surface tension in N/m from coherent Joback Tb/Tc/Pc
// plus functional hydroxyl classification (Alcohol vs GeneralOrganic).
// Phase 16: Sato-Riedel liquid kL(T) and Gharagheizi gas kG(T) in W/(m*K).
// Not live ThermalEngine. No phase-specific conductivity mapping yet.
// No published universal Joback viscosity temperature interval is claimed.
// StructuralEstimate / Low. Unsupported chemistry stays Unknown.
// Not canonical identity. Not live PhaseProperties. Do not call from physics ticks.
//
// Applicability: acyclic, net-neutral, fully single-bonded H/C/O small molecules
// coverable by -CH3, -CH2-, >CH-, >C<, alcohol -OH, non-ring -O-.
// Critical molar volume is not ordinary liquid molar volume.

constexpr double kJobackTbInterceptK = 198.2;
constexpr double kJobackTbCH3K = 23.58;
constexpr double kJobackTbCH2K = 22.88;
constexpr double kJobackTbCHK = 21.74;
constexpr double kJobackTbCK = 18.25;
constexpr double kJobackTbAlcoholOHK = 92.88;
constexpr double kJobackTbEtherOK = 22.42;

constexpr double kJobackTcCH3 = 0.0141;
constexpr double kJobackTcCH2 = 0.0189;
constexpr double kJobackTcCH = 0.0164;
constexpr double kJobackTcC = 0.0067;
constexpr double kJobackTcAlcoholOH = 0.0741;
constexpr double kJobackTcEtherO = 0.0168;

constexpr double kJobackPcCH3 = -0.0012;
constexpr double kJobackPcCH2 = 0.0000;
constexpr double kJobackPcCH = 0.0020;
constexpr double kJobackPcC = 0.0043;
constexpr double kJobackPcAlcoholOH = 0.0112;
constexpr double kJobackPcEtherO = 0.0015;

constexpr double kJobackVcCH3Cm3 = 65.0;
constexpr double kJobackVcCH2Cm3 = 56.0;
constexpr double kJobackVcCHCm3 = 41.0;
constexpr double kJobackVcCCm3 = 27.0;
constexpr double kJobackVcAlcoholOHCm3 = 28.0;
constexpr double kJobackVcEtherOCm3 = 18.0;

constexpr double kJobackVcInterceptCm3 = 17.5;
constexpr double kJobackBarToPa = 100000.0;
constexpr double kJobackCm3ToM3 = 1.0e-6;

constexpr double kJobackHvapInterceptKJPerMol = 15.30;
constexpr double kJobackHvapCH3KJPerMol = 2.373;
constexpr double kJobackHvapCH2KJPerMol = 2.226;
constexpr double kJobackHvapCHKJPerMol = 1.691;
constexpr double kJobackHvapCKJPerMol = 0.636;
constexpr double kJobackHvapAlcoholOHKJPerMol = 16.826;
constexpr double kJobackHvapEtherOKJPerMol = 2.410;
constexpr double kJobackKJToJ = 1000.0;
constexpr double kWatsonHvapExponent = 0.38;
constexpr double kSaceHeatCapacityTMinK = 298.0;
constexpr double kSaceHeatCapacityTMaxK = 1000.0;
constexpr double kSaceHeatCapacityReferenceK = 298.15;
constexpr double kSaceLiquidCpTrReject = 0.98;
constexpr double kSaceCostaldTrMinExclusive = 0.25;
constexpr double kSaceCostaldTrMaxExclusive = 0.95;

constexpr double kJobackViscAOffset = 597.82;
constexpr double kJobackViscBOffset = 11.202;
constexpr double kJobackViscMuACH3 = 548.29;
constexpr double kJobackViscMuBCH3 = -1.719;
constexpr double kJobackViscMuACH2 = 94.16;
constexpr double kJobackViscMuBCH2 = -0.199;
constexpr double kJobackViscMuACH = -322.15;
constexpr double kJobackViscMuBCH = 1.187;
constexpr double kJobackViscMuAC = -573.56;
constexpr double kJobackViscMuBC = 2.307;
constexpr double kJobackViscMuAAlcoholOH = 2173.72;
constexpr double kJobackViscMuBAlcoholOH = -5.057;
constexpr double kJobackViscMuAEtherO = 122.09;
constexpr double kJobackViscMuBEtherO = -0.386;

constexpr double kSastriRaoPcPaToBar = 1.0e-5;
constexpr double kSastriRaoMNPerMToNPerM = 1.0e-3;
constexpr double kSastriRaoGeneralK = 0.158;
constexpr double kSastriRaoGeneralX = 0.50;
constexpr double kSastriRaoGeneralY = -1.5;
constexpr double kSastriRaoGeneralZ = 1.85;
constexpr double kSastriRaoGeneralM = 11.0 / 9.0;
constexpr double kSastriRaoAlcoholK = 2.28;
constexpr double kSastriRaoAlcoholX = 0.25;
constexpr double kSastriRaoAlcoholY = 0.175;
constexpr double kSastriRaoAlcoholZ = 0.0;
constexpr double kSastriRaoAlcoholM = 0.8;

constexpr double kSatoRiedelPrefactor = 1.1053;
constexpr double kGharagheiziPcScale = 1.0e-4;
constexpr double kGharagheiziConst32825 = 3.2825;
constexpr double kGharagheiziX0Omega = 3.9752;
constexpr double kGharagheiziX0P = 0.1;
constexpr double kGharagheiziX0B = 1.9876;
constexpr double kGharagheiziX0Const = 6.5243;
constexpr double kGharagheiziK0 = 7.9505e-4;
constexpr double kGharagheiziKT = 3.989e-5;
constexpr double kGharagheiziKMW = 5.419e-5;
constexpr double kGharagheiziKA = 3.989e-5;

constexpr double kJobackCpAIntercept = -37.93;
constexpr double kJobackCpBIntercept = 0.210;
constexpr double kJobackCpCIntercept = -3.91e-4;
constexpr double kJobackCpDIntercept = 2.06e-7;

constexpr double kJobackCpACH3 = 19.5;
constexpr double kJobackCpBCH3 = -8.08e-3;
constexpr double kJobackCpCCH3 = 1.53e-4;
constexpr double kJobackCpDCH3 = -9.67e-8;

constexpr double kJobackCpACH2 = -0.909;
constexpr double kJobackCpBCH2 = 9.50e-2;
constexpr double kJobackCpCCH2 = -5.44e-5;
constexpr double kJobackCpDCH2 = 1.19e-8;

constexpr double kJobackCpACH = -23.0;
constexpr double kJobackCpBCH = 0.204;
constexpr double kJobackCpCCH = -2.65e-4;
constexpr double kJobackCpDCH = 1.20e-7;

constexpr double kJobackCpAC = -66.2;
constexpr double kJobackCpBC = 0.427;
constexpr double kJobackCpCC = -6.41e-4;
constexpr double kJobackCpDC = 3.01e-7;

constexpr double kJobackCpAAlcoholOH = 25.7;
constexpr double kJobackCpBAlcoholOH = -6.91e-2;
constexpr double kJobackCpCAlcoholOH = 1.77e-4;
constexpr double kJobackCpDAlcoholOH = -9.88e-8;

constexpr double kJobackCpAEtherO = 25.5;
constexpr double kJobackCpBEtherO = -6.32e-2;
constexpr double kJobackCpCEtherO = 1.11e-4;
constexpr double kJobackCpDEtherO = -5.48e-8;

enum class SaceJobackGroup : uint8_t {
    CarbonCH3 = 0,
    CarbonCH2,
    CarbonCH,
    CarbonC,
    HydroxylAlcohol,
    EtherNonRing
};

struct SaceJobackGroupCounts {
    int carbonCH3 = 0;
    int carbonCH2 = 0;
    int carbonCH = 0;
    int carbonC = 0;
    int hydroxylAlcohol = 0;
    int etherNonRing = 0;
};

struct SaceJobackEstimateBundle {
    SaceJobackGroupCounts groups{};
    SaceScalarProperty normalBoilingPointK{};
    SaceScalarProperty criticalTemperatureK{};
    SaceScalarProperty criticalPressurePa{};
    SaceScalarProperty criticalMolarVolumeM3PerMol{};
    SaceScalarProperty acentricFactor{};
    SaceScalarProperty enthalpyVaporizationAtNormalBoilingJPerMol{};
    SaceScalarProperty idealGasHeatCapacityAt298KJPerMolK{};
    SaceScalarProperty saturatedLiquidHeatCapacityAt298KJPerMolK{};
    SaceScalarProperty saturatedLiquidDensityAt298KKgPerM3{};
    SaceScalarProperty liquidDynamicViscosityAt298KPaS{};
    SaceScalarProperty liquidSurfaceTensionAt298KNPerM{};
    SaceScalarProperty liquidThermalConductivityAt298KWPerMK{};
    SaceScalarProperty gasThermalConductivityAt298KWPerMK{};
};

// Lee-Kesler 1975 corresponding-states vapor-pressure model.
// Built from a self-consistent Joback Tb/Tc/Pc/omega tuple, not mixed
// independently replaced SaceGeneratedProperties. Not canonical identity.
// Not live phase behavior. Independent Reference/Structural fields must not
// be assumed to form a coherent Psat(T) set.
struct SaceLeeKeslerVaporModel {
    double normalBoilingPointK = 0.0;
    double criticalTemperatureK = 0.0;
    double criticalPressurePa = 0.0;
    double acentricFactor = 0.0;
    bool valid = false;
};

constexpr double kLeeKeslerAtmPa = 101325.0;

// Watson temperature scaling of Joback ΔHvap(Tb). Coherent structural set only:
// graph -> one Joback bundle -> Watson model. Not mixed stored properties.
// Not canonical identity. Not live PhaseChangeEngine. Exponent is 0.38.
// Values below Tb are mathematical extrapolations, not proof of liquid stability.
struct SaceWatsonVaporizationModel {
    double referenceTemperatureK = 0.0;
    double criticalTemperatureK = 0.0;
    double referenceEnthalpyJPerMol = 0.0;
    double molarMassKgPerMol = 0.0;
    bool valid = false;
};

struct SaceJobackIdealGasCpModel {
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    bool valid = false;
};

// Rowlinson-Poling liquid Cp from coherent Joback Cp + Joback Tc + Lee-Kesler omega.
// Alcohols (Graph A) are a difficult associating case; no empirical correction.
// Fail closed for Tr >= 0.98. Not live ThermalEngine. Not mixed stored properties.
struct SaceRowlinsonPolingLiquidCpModel {
    SaceJobackIdealGasCpModel idealGasCp{};
    double criticalTemperatureK = 0.0;
    double acentricFactor = 0.0;
    double molarMassKgPerMol = 0.0;
    bool valid = false;
};

// COSTALD saturated-liquid volume (Hankinson–Thomson 1979).
// Phase 12 has no fitted COSTALD V* / omega_SRK; Joback Vc is the V* fallback
// and SACE Lee-Kesler omega is the available acentric estimate. Low confidence.
// Saturated density at Psat(T), not compressed-liquid rho(T,P).
// Strict domain: 0.25 < Tr < 0.95. Not live FluidEngine. Not kg/m3 -> relative density.
struct SaceCostaldLiquidDensityModel {
    double criticalTemperatureK = 0.0;
    double characteristicVolumeM3PerMol = 0.0;
    double acentricFactor = 0.0;
    double molarMassKgPerMol = 0.0;
    bool valid = false;
};

// Joback liquid dynamic viscosity. MW is the published numerical g/mol
// convention in mu = MW * exp(A/T + B). Output is physical Pa*s, not
// FluidProperties::viscosity sandbox units. No invented T interval.
// Liquid-domain validity is assessed separately (COSTALD gate for the
// stored 298.15 K scalar). Not live FluidEngine.
struct SaceJobackLiquidViscosityModel {
    double aKelvin = 0.0;
    double bDimensionless = 0.0;
    double molarMassGPerMol = 0.0;
    bool valid = false;
};

enum class SaceSurfaceTensionClass : uint8_t {
    GeneralOrganic = 0,
    Alcohol
};

// Sastri-Rao 1995 liquid surface tension from one coherent Joback Tb/Tc/Pc
// tuple plus SACE functional hydroxyl classification. Output is physical N/m,
// not FluidProperties::surfaceTension sandbox units. No acid branch.
// A value below an unknown melting/triple point is not proof of liquid stability.
// sigma(Tc) = 0. Not live FluidEngine.
struct SaceSastriRaoSurfaceTensionModel {
    double normalBoilingPointK = 0.0;
    double criticalTemperatureK = 0.0;
    double criticalPressurePa = 0.0;
    SaceSurfaceTensionClass chemicalClass = SaceSurfaceTensionClass::GeneralOrganic;
    bool valid = false;
};

// Sato-Riedel corresponding-states liquid thermal conductivity.
// kL = (1.1053 / sqrt(MW[g/mol])) * [3+20(1-Tr)^(2/3)] / [3+20(1-Tbr)^(2/3)]
// Output W/(m*K). Approximate corresponding-states estimate, not reference data.
// Not high-pressure liquid k. Not proof of liquid stability. Not live ThermalEngine.
struct SaceSatoRiedelLiquidConductivityModel {
    double normalBoilingPointK = 0.0;
    double criticalTemperatureK = 0.0;
    double molarMassGPerMol = 0.0;
    bool valid = false;
};

// Gharagheizi pure-gas thermal conductivity. Pc stored in Pa; internal scale
// P = PcPa * 1e-4 is the corrected implementation convention, not ordinary bar.
// Output W/(m*K). Not a pressure-dependent transport model. Not live ThermalEngine.
struct SaceGharagheiziGasConductivityModel {
    double normalBoilingPointK = 0.0;
    double criticalPressurePa = 0.0;
    double acentricFactor = 0.0;
    double molarMassGPerMol = 0.0;
    bool valid = false;
};

enum class SaceJobackFragmentationResult : uint8_t {
    Ok = 0,
    InvalidGraph,
    UnsupportedElement,
    Charged,
    Cyclic,
    Unsaturated,
    UnsupportedCarbonEnvironment,
    UnsupportedOxygenEnvironment,
    UnsupportedHydrogenEnvironment,
    IncompleteCoverage
};

inline char const *saceJobackFragmentationResultKey(SaceJobackFragmentationResult r) {
    switch (r) {
        case SaceJobackFragmentationResult::Ok: return "ok";
        case SaceJobackFragmentationResult::InvalidGraph: return "invalid-graph";
        case SaceJobackFragmentationResult::UnsupportedElement: return "unsupported-element";
        case SaceJobackFragmentationResult::Charged: return "charged";
        case SaceJobackFragmentationResult::Cyclic: return "cyclic";
        case SaceJobackFragmentationResult::Unsaturated: return "unsaturated";
        case SaceJobackFragmentationResult::UnsupportedCarbonEnvironment: return "unsupported-carbon";
        case SaceJobackFragmentationResult::UnsupportedOxygenEnvironment: return "unsupported-oxygen";
        case SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment: return "unsupported-hydrogen";
        case SaceJobackFragmentationResult::IncompleteCoverage: return "incomplete-coverage";
    }
    return "unknown";
}

inline double saceJobackTcContributionSum(SaceJobackGroupCounts const &g) {
    return kJobackTcCH3 * g.carbonCH3
        + kJobackTcCH2 * g.carbonCH2
        + kJobackTcCH * g.carbonCH
        + kJobackTcC * g.carbonC
        + kJobackTcAlcoholOH * g.hydroxylAlcohol
        + kJobackTcEtherO * g.etherNonRing;
}

inline double saceJobackPcContributionSum(SaceJobackGroupCounts const &g) {
    return kJobackPcCH3 * g.carbonCH3
        + kJobackPcCH2 * g.carbonCH2
        + kJobackPcCH * g.carbonCH
        + kJobackPcC * g.carbonC
        + kJobackPcAlcoholOH * g.hydroxylAlcohol
        + kJobackPcEtherO * g.etherNonRing;
}

inline double saceJobackVcContributionSumCm3PerMol(SaceJobackGroupCounts const &g) {
    return kJobackVcCH3Cm3 * g.carbonCH3
        + kJobackVcCH2Cm3 * g.carbonCH2
        + kJobackVcCHCm3 * g.carbonCH
        + kJobackVcCCm3 * g.carbonC
        + kJobackVcAlcoholOHCm3 * g.hydroxylAlcohol
        + kJobackVcEtherOCm3 * g.etherNonRing;
}

inline double saceJobackViscosityMuASum(SaceJobackGroupCounts const &g) {
    return kJobackViscMuACH3 * g.carbonCH3
        + kJobackViscMuACH2 * g.carbonCH2
        + kJobackViscMuACH * g.carbonCH
        + kJobackViscMuAC * g.carbonC
        + kJobackViscMuAAlcoholOH * g.hydroxylAlcohol
        + kJobackViscMuAEtherO * g.etherNonRing;
}

inline double saceJobackViscosityMuBSum(SaceJobackGroupCounts const &g) {
    return kJobackViscMuBCH3 * g.carbonCH3
        + kJobackViscMuBCH2 * g.carbonCH2
        + kJobackViscMuBCH * g.carbonCH
        + kJobackViscMuBC * g.carbonC
        + kJobackViscMuBAlcoholOH * g.hydroxylAlcohol
        + kJobackViscMuBEtherO * g.etherNonRing;
}

inline double saceJobackHvapContributionSumKJPerMol(SaceJobackGroupCounts const &g) {
    return kJobackHvapCH3KJPerMol * g.carbonCH3
        + kJobackHvapCH2KJPerMol * g.carbonCH2
        + kJobackHvapCHKJPerMol * g.carbonCH
        + kJobackHvapCKJPerMol * g.carbonC
        + kJobackHvapAlcoholOHKJPerMol * g.hydroxylAlcohol
        + kJobackHvapEtherOKJPerMol * g.etherNonRing;
}

SaceJobackFragmentationResult saceFragmentJobackGroups(
    SaceMolecularGraph const &graph, SaceJobackGroupCounts &counts);

bool saceJobackGroupCountsEqual(SaceJobackGroupCounts const &a, SaceJobackGroupCounts const &b);

void saceResetJobackEstimateBundle(SaceJobackEstimateBundle &out);
bool saceEstimateJobackBundle(SaceMolecularGraph const &graph, SaceJobackEstimateBundle &out);

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out);
bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out, SaceJobackGroupCounts &counts);

bool saceLeeKeslerSaturationPressurePa(
    double temperatureK,
    double criticalTemperatureK,
    double criticalPressurePa,
    double acentricFactor,
    double &outPressurePa);

bool saceBuildLeeKeslerVaporModelFromJoback(
    SaceJobackEstimateBundle const &bundle,
    SaceLeeKeslerVaporModel &out);

bool saceBuildWatsonVaporizationModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceWatsonVaporizationModel &out);

bool saceWatsonEnthalpyVaporizationJPerMol(
    SaceWatsonVaporizationModel const &model,
    double temperatureK,
    double &outJPerMol);

bool saceWatsonLatentHeatVaporizationJPerKg(
    SaceWatsonVaporizationModel const &model,
    double temperatureK,
    double &outJPerKg);

bool saceBuildJobackIdealGasCpModel(
    SaceJobackGroupCounts const &groups,
    SaceJobackIdealGasCpModel &out);
bool saceJobackIdealGasHeatCapacityJPerMolK(
    SaceJobackIdealGasCpModel const &model,
    double temperatureK,
    double &outJPerMolK);

bool saceBuildRowlinsonPolingLiquidCpModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceRowlinsonPolingLiquidCpModel &out);
bool saceRowlinsonPolingLiquidHeatCapacityJPerMolK(
    SaceRowlinsonPolingLiquidCpModel const &model,
    double temperatureK,
    double &outJPerMolK);

bool saceMolarHeatCapacityToSpecificJPerKgK(
    double cpJPerMolK,
    double molarMassKgPerMol,
    double &outJPerKgK);

bool saceBuildCostaldLiquidDensityModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceCostaldLiquidDensityModel &out);
bool saceCostaldSaturatedLiquidMolarVolumeM3PerMol(
    SaceCostaldLiquidDensityModel const &model,
    double temperatureK,
    double &outM3PerMol);
bool saceCostaldSaturatedLiquidDensityKgPerM3(
    SaceCostaldLiquidDensityModel const &model,
    double temperatureK,
    double &outKgPerM3);

bool saceBuildJobackLiquidViscosityModel(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceJobackLiquidViscosityModel &out);
bool saceJobackLiquidViscosityPaS(
    SaceJobackLiquidViscosityModel const &model,
    double temperatureK,
    double &outPaS);

bool saceBuildSastriRaoSurfaceTensionModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceSastriRaoSurfaceTensionModel &out);
bool saceSastriRaoSurfaceTensionNPerM(
    SaceSastriRaoSurfaceTensionModel const &model,
    double temperatureK,
    double &outNPerM);

bool saceBuildSatoRiedelLiquidConductivityModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceSatoRiedelLiquidConductivityModel &out);
bool saceSatoRiedelLiquidThermalConductivityWPerMK(
    SaceSatoRiedelLiquidConductivityModel const &model,
    double temperatureK,
    double &outWPerMK);

bool saceBuildGharagheiziGasConductivityModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceGharagheiziGasConductivityModel &out);
bool saceGharagheiziGasThermalConductivityWPerMK(
    SaceGharagheiziGasConductivityModel const &model,
    double temperatureK,
    double &outWPerMK);

void runSaceEstimationDiagnostics();
