#pragma once

#include "substance/SubstanceProperties.h"

#include <cstdint>

// Engine identity for a substance. Compact, stable, compile-time for built-ins.
// Not a player-facing Element #HEX id (that is a future SACE identifier).
// Not a phase: the same SubstanceId may later exist as solid, liquid, or gas.
using SubstanceId = uint16_t;

constexpr SubstanceId SUBSTANCE_NONE = 0;
constexpr SubstanceId SUBSTANCE_WATER = 1;
constexpr SubstanceId SUBSTANCE_HONEY = 2;
constexpr SubstanceId SUBSTANCE_WOOD = 3;
constexpr SubstanceId SUBSTANCE_STONE = 4;
constexpr SubstanceId SUBSTANCE_GLASS = 5;
constexpr SubstanceId SUBSTANCE_METAL = 6;
constexpr SubstanceId SUBSTANCE_AIR = 7;
constexpr SubstanceId SUBSTANCE_COUNT = 8;

enum class SubstanceClass : uint8_t {
    Unknown = 0,
    PureSubstance,
    Mixture,
    Composite,
    Biological
};

// Placeholder for later SACE/composition work. Do not invent formulas.
enum class CompositionKind : uint8_t {
    Unknown = 0,
    PureChemical,
    Mixture,
    Composite
};

struct SubstanceDefinition {
    SubstanceId id = SUBSTANCE_NONE;
    char const *internalName = "none";
    char const *displayName = "(none)";
    char const *displayNameKey = "ins_mat_none";
    SubstanceClass classification = SubstanceClass::Unknown;
    CompositionKind compositionKind = CompositionKind::Unknown;
    char const *formulaHint = nullptr; // only when genuinely known (e.g. water H2O)
    bool canExistAsSolid = false;
    bool canExistAsLiquid = false;
    bool canExistAsGas = false;

    MechanicalProperties mechanical{};
    FluidProperties fluid{};
    ThermalProperties thermal{};
    PhaseProperties phase{};
    PorousProperties porous{};
    ChemicalProperties chemical{};       // placeholder; not simulated
    ElectricalProperties electrical{};   // placeholder; not simulated
    SubstanceVisualMetadata visual{};    // RGB seed; renderer still owns MaterialVisual
};

// Function-local static table in SubstanceRegistry.cpp (safe first-call init).
SubstanceDefinition const *builtinSubstanceTable();

inline bool validSubstance(SubstanceId id) {
    return id < SUBSTANCE_COUNT;
}

inline SubstanceDefinition const &substanceDef(SubstanceId id) {
    SubstanceDefinition const *table = builtinSubstanceTable();
    if (id >= SUBSTANCE_COUNT) return table[SUBSTANCE_NONE];
    return table[id];
}

inline MechanicalProperties const &mechanicalForSubstance(SubstanceId id) {
    return substanceDef(id).mechanical;
}
inline FluidProperties const &fluidForSubstance(SubstanceId id) {
    return substanceDef(id).fluid;
}
inline ThermalProperties const &thermalForSubstance(SubstanceId id) {
    return substanceDef(id).thermal;
}
inline PorousProperties const &porousForSubstance(SubstanceId id) {
    return substanceDef(id).porous;
}
inline PhaseProperties const &phaseForSubstance(SubstanceId id) {
    return substanceDef(id).phase;
}
inline ChemicalProperties const &chemicalForSubstance(SubstanceId id) {
    return substanceDef(id).chemical;
}
inline ElectricalProperties const &electricalForSubstance(SubstanceId id) {
    return substanceDef(id).electrical;
}

// MaterialId lives in rigid/RigidBodyTypes.h (0 empty, 1 wood, 2 stone, 3 glass, 4 metal).
// Numeric bridge so thermal/fluid can map without including rigid headers.
inline SubstanceId substanceForMaterialId(uint16_t materialId) {
    switch (materialId) {
        case 1: return SUBSTANCE_WOOD;
        case 2: return SUBSTANCE_STONE;
        case 3: return SUBSTANCE_GLASS;
        case 4: return SUBSTANCE_METAL;
        default: return SUBSTANCE_NONE;
    }
}

inline ThermalProperties const &thermalForMaterial(uint16_t materialId) {
    return thermalForSubstance(substanceForMaterialId(materialId));
}

inline SubstanceId substanceForLiquidPaint(bool asHoney) {
    return asHoney ? SUBSTANCE_HONEY : SUBSTANCE_WATER;
}

inline SubstanceId substanceForHoneyFraction(float honeyFrac) {
    return honeyFrac > 0.5f ? SUBSTANCE_HONEY : SUBSTANCE_WATER;
}

inline SubstanceId substanceForAmbientGas() {
    return SUBSTANCE_AIR;
}
