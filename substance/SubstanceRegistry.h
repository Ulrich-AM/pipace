#pragma once

#include "substance/SubstanceTypes.h"
#include "rigid/RigidBodyTypes.h"

// Canonical built-in table lives in SubstanceRegistry.cpp (function-local static).
// Access: substanceDef(id) and the grouped accessors in SubstanceTypes.h.

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

// Headless: --substance-registry-diag  -> misc/substance_registry_diag.tsv
void runSubstanceRegistryDiagnostics();
