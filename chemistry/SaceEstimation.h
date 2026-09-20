#pragma once

#include "chemistry/SaceMolecule.h"
#include "chemistry/SaceProperties.h"
#include "chemistry/SaceDescriptors.h"

#include <cstdint>

// SACE Phase 7: Joback-Reid 1987 normal boiling point, published subset only.
// Tb[K] = 198.2 + SUM(group contributions). Coefficients are not fitted.
// StructuralEstimate / Low. Unsupported chemistry stays Unknown.
// Not canonical identity. Not live PhaseProperties. Do not call from physics ticks.
//
// Applicability: acyclic, net-neutral, fully single-bonded H/C/O small molecules
// coverable by -CH3, -CH2-, >CH-, >C<, alcohol -OH, non-ring -O-.

constexpr double kJobackTbInterceptK = 198.2;
constexpr double kJobackTbCH3K = 23.58;
constexpr double kJobackTbCH2K = 22.88;
constexpr double kJobackTbCHK = 21.74;
constexpr double kJobackTbCK = 18.25;
constexpr double kJobackTbAlcoholOHK = 92.88;
constexpr double kJobackTbEtherOK = 22.42;

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

SaceJobackFragmentationResult saceFragmentJobackGroups(
    SaceMolecularGraph const &graph, SaceJobackGroupCounts &counts);

bool saceJobackGroupCountsEqual(SaceJobackGroupCounts const &a, SaceJobackGroupCounts const &b);

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out);
bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out, SaceJobackGroupCounts &counts);

void runSaceEstimationDiagnostics();
