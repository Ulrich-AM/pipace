#pragma once

#include "substance/SubstanceTypes.h"
#include "rigid/RigidBodyTypes.h"

// Canonical built-in table. Indexed by SubstanceId. Prompt 1: compile-time IDs only.
extern SubstanceDefinition const kBuiltinSubstances[SUBSTANCE_COUNT];

inline bool validSubstance(SubstanceId id) {
    return id < SUBSTANCE_COUNT;
}

inline SubstanceDefinition const &substanceDef(SubstanceId id) {
    if (id >= SUBSTANCE_COUNT) return kBuiltinSubstances[SUBSTANCE_NONE];
    return kBuiltinSubstances[id];
}

inline SubstanceId substanceForMaterial(MaterialId material) {
    return substanceForMaterialId(material);
}

// Compatibility: rigid masks still store MaterialId. Water/honey/air have no rigid mask.
inline MaterialId rigidMaterialForSubstance(SubstanceId substance) {
    switch (substance) {
        case SUBSTANCE_WOOD: return MATERIAL_WOOD;
        case SUBSTANCE_STONE: return MATERIAL_STONE;
        case SUBSTANCE_GLASS: return MATERIAL_GLASS;
        case SUBSTANCE_METAL: return MATERIAL_METAL;
        default: return MATERIAL_EMPTY;
    }
}

SubstanceId substanceFromInternalName(char const *name);
