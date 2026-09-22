#pragma once

#include "chemistry/SaceMolecule.h"
#include "chemistry/SaceProperties.h"
#include "chemistry/SaceDescriptors.h"

#include <cstdint>

// SACE Phase 7–10: Joback-Reid 1987 subset + Lee-Kesler Psat + Watson Hvap.
// Phase 7: Tb[K] = 198.2 + SUM(Tb groups).
// Phase 8: Tc, Pc, Vc from the same six groups and the Joback Tb (not stored Tb).
// Phase 9: Lee-Kesler omega and Psat(T) from the same Joback Tb/Tc/Pc tuple.
// Phase 10: Joback ΔHvap(Tb) [J/mol] + Watson Hvap(T). Not live phase physics.
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

void runSaceEstimationDiagnostics();
