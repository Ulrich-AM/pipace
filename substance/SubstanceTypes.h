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
constexpr SubstanceId SUBSTANCE_HYDROGEN = 8;
constexpr SubstanceId SUBSTANCE_OXYGEN = 9;
constexpr SubstanceId SUBSTANCE_COUNT = 10;

// Static world walls currently use the stone solid table. This is an explicit
// world-boundary identity, not a rule that every future solid is stone.
constexpr SubstanceId kStaticWallSubstance = SUBSTANCE_STONE;

// Current world representation. Not intrinsic identity. Same SubstanceId may
// later be Solid, Liquid, or Gas. Plasma is a placeholder; do not simulate it.
// Do not store this on SubstanceDefinition — phase is world state.
enum class MatterPhase : uint8_t {
    None = 0,
    Solid,
    Liquid,
    Gas,
    Plasma
};

struct MatterIdentity {
    SubstanceId substance = SUBSTANCE_NONE;
    MatterPhase phase = MatterPhase::None;
};

// Read-only world sample. Not for hot sim loops. Mixtures report the dominant
// component plus water/honey fractions; no generated mixture SubstanceId.
struct MatterSample {
    MatterIdentity identity{};
    bool hasMatter = false;
    bool mixture = false;
    float amount = 0.0f;
    float temperatureK = AMBIENT_TEMPERATURE_K;
    float waterFraction = 0.0f;
    float honeyFraction = 0.0f;
    float vaporFraction = 0.0f;
};

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
    // Phase capability lives on phase.*Capable. Current world phase is NOT here.

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

inline float gasPhaseSpecificHeat(SubstanceId id) {
    ThermalProperties const &t = thermalForSubstance(id);
    if (t.gasSpecificHeat > 1.0f) return t.gasSpecificHeat;
    return t.specificHeat;
}

inline bool supportsPhase(SubstanceId id, MatterPhase phase) {
    PhaseProperties const &p = phaseForSubstance(id);
    switch (phase) {
        case MatterPhase::Solid: return p.solidCapable;
        case MatterPhase::Liquid: return p.liquidCapable;
        case MatterPhase::Gas: return p.gasCapable;
        case MatterPhase::Plasma: return false;
        case MatterPhase::None: return false;
    }
    return false;
}

inline bool hasMechanicalProperties(SubstanceId id) {
    return substanceDef(id).mechanical.valid;
}
inline bool hasFluidProperties(SubstanceId id) {
    return substanceDef(id).fluid.valid;
}
inline bool hasGasProperties(SubstanceId id) {
    return phaseForSubstance(id).gasCapable;
}

// Engine selection must not be inferred from SubstanceId alone.
// These queries describe which property groups exist today for a phase.
inline bool hasPropertiesForPhase(SubstanceId id, MatterPhase phase) {
    switch (phase) {
        case MatterPhase::Solid: return hasMechanicalProperties(id);
        case MatterPhase::Liquid: return hasFluidProperties(id);
        case MatterPhase::Gas: return hasGasProperties(id);
        default: return false;
    }
}

inline char const *matterPhaseKey(MatterPhase phase) {
    switch (phase) {
        case MatterPhase::Solid: return "ins_phase_solid";
        case MatterPhase::Liquid: return "ins_phase_liquid";
        case MatterPhase::Gas: return "ins_phase_gas";
        case MatterPhase::Plasma: return "ins_phase_plasma";
        case MatterPhase::None: return "ins_phase_none";
    }
    return "ins_phase_none";
}

inline MatterIdentity makeMatterIdentity(SubstanceId substance, MatterPhase phase) {
    MatterIdentity id;
    id.substance = substance;
    id.phase = phase;
    return id;
}

inline MatterIdentity staticWallIdentity() {
    return makeMatterIdentity(kStaticWallSubstance, MatterPhase::Solid);
}

inline MatterIdentity ambientGasIdentity() {
    return makeMatterIdentity(SUBSTANCE_AIR, MatterPhase::Gas);
}

// Solver reference liquid: relative density 1.0 and the thin-liquid viscosity
// origin. Instantiated as SUBSTANCE_WATER today, but this is a unit reference,
// not a claim that every liquid cell is water.
inline FluidProperties const &sandboxReferenceLiquid() {
    return fluidForSubstance(SUBSTANCE_WATER);
}

// MaterialId lives in rigid/RigidBodyTypes.h (0 empty, 1 wood, 2 stone, 3 glass, 4 metal,
// 5 water-solid). Numeric bridge so thermal/fluid can map without including rigid headers.
inline SubstanceId substanceForMaterialId(uint16_t materialId) {
    switch (materialId) {
        case 1: return SUBSTANCE_WOOD;
        case 2: return SUBSTANCE_STONE;
        case 3: return SUBSTANCE_GLASS;
        case 4: return SUBSTANCE_METAL;
        case 5: return SUBSTANCE_WATER; // MATERIAL_WATER_SOLID — still WATER
        default: return SUBSTANCE_NONE;
    }
}

inline float solidPhaseSpecificHeat(SubstanceId id) {
    float v = thermalForSubstance(id).solidSpecificHeat;
    if (v > 1.0f) return v;
    return thermalForSubstance(id).specificHeat;
}

inline float solidPhaseConductivity(SubstanceId id) {
    float v = thermalForSubstance(id).solidConductivity;
    if (v > 0.0f) return v;
    return thermalForSubstance(id).conductivity;
}

inline MatterIdentity rigidIdentityForMaterial(uint16_t materialId) {
    SubstanceId sid = substanceForMaterialId(materialId);
    if (sid == SUBSTANCE_NONE) return makeMatterIdentity(SUBSTANCE_NONE, MatterPhase::None);
    return makeMatterIdentity(sid, MatterPhase::Solid);
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
