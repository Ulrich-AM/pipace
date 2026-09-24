#pragma once

#include "chemistry/SaceTypes.h"

#include <cstdint>

// SACE Phase 3: property provenance + identity-derived molar mass.
// Identity may be exact while properties remain unknown.
// Properties never participate in canonical identity.
// Do not call from physics ticks.

enum class SacePropertySource : uint8_t {
    Unknown = 0,
    Reference,
    IdentityDerived,
    EmpiricalEstimate,
    StructuralEstimate,
    MixtureEstimate,
    Fallback
};

enum class SaceConfidence : uint8_t {
    Unknown = 0,
    Low,
    Medium,
    High
};

// Future precedence (higher wins; never silently overwrite with worse data):
//   Reference > IdentityDerived > StructuralEstimate / EmpiricalEstimate
//   > MixtureEstimate > Fallback
// Equal source rank: higher confidence wins. Exact source-rank + confidence
// tie preserves the existing value.
inline int sacePropertySourceRank(SacePropertySource source) {
    switch (source) {
        case SacePropertySource::Reference: return 5;
        case SacePropertySource::IdentityDerived: return 4;
        case SacePropertySource::StructuralEstimate: return 3;
        case SacePropertySource::EmpiricalEstimate: return 3;
        case SacePropertySource::MixtureEstimate: return 2;
        case SacePropertySource::Fallback: return 1;
        case SacePropertySource::Unknown: return 0;
    }
    return 0;
}

inline int saceConfidenceRank(SaceConfidence c) {
    switch (c) {
        case SaceConfidence::High: return 3;
        case SaceConfidence::Medium: return 2;
        case SaceConfidence::Low: return 1;
        case SaceConfidence::Unknown: return 0;
    }
    return 0;
}

struct SaceScalarProperty {
    float value = 0.0f;
    bool known = false;
    SacePropertySource source = SacePropertySource::Unknown;
    SaceConfidence confidence = SaceConfidence::Unknown;
};

inline SaceScalarProperty saceUnknownScalarProperty() {
    return {};
}

struct SaceGeneratedProperties {
    SaceScalarProperty molarMassGPerMol{};
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
};

inline char const *sacePropertySourceKey(SacePropertySource source) {
    switch (source) {
        case SacePropertySource::Reference: return "reference";
        case SacePropertySource::IdentityDerived: return "identity-derived";
        case SacePropertySource::EmpiricalEstimate: return "empirical-estimate";
        case SacePropertySource::StructuralEstimate: return "structural-estimate";
        case SacePropertySource::MixtureEstimate: return "mixture-estimate";
        case SacePropertySource::Fallback: return "fallback";
        case SacePropertySource::Unknown: return "unknown";
    }
    return "unknown";
}

inline char const *saceConfidenceKey(SaceConfidence c) {
    switch (c) {
        case SaceConfidence::Low: return "low";
        case SaceConfidence::Medium: return "medium";
        case SaceConfidence::High: return "high";
        case SaceConfidence::Unknown: return "unknown";
    }
    return "unknown";
}

// Tiny H/C/O reference table. Unknown Z fails; no invented averages.
bool saceAtomicMassGPerMol(AtomicNumber z, double &out);

// Composition arithmetic only. Any missing atomic mass => unknown, not a partial sum.
bool saceDeriveMolarMass(ElementalComposition const &elemental, SaceScalarProperty &out);

// Assign only if src is known and not worse than an already-known destination.
bool saceAssignScalarProperty(SaceScalarProperty &dst, SaceScalarProperty const &src);

void initializeGeneratedProperties(SaceGeneratedProperties &props, ElementalComposition const &elemental);

// Headless: --sace-property-diag -> misc/sace_property_diag.tsv
void runSacePropertyDiagnostics();
