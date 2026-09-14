#include "substance/SubstanceRegistry.h"

#include <cstring>

SubstanceDefinition const kBuiltinSubstances[SUBSTANCE_COUNT] = {
    {SUBSTANCE_NONE,  "none",  "(none)", "ins_mat_none",
        SubstanceClass::Unknown, CompositionKind::Unknown, nullptr,
        false, false, false},
    {SUBSTANCE_WATER, "water", "Water", "ins_mat_water",
        SubstanceClass::PureSubstance, CompositionKind::PureChemical, "H2O",
        true, true, true},
    {SUBSTANCE_HONEY, "honey", "Honey", "ins_mat_honey",
        SubstanceClass::Mixture, CompositionKind::Mixture, nullptr,
        true, true, false},
    {SUBSTANCE_WOOD,  "wood",  "wood", "wood",
        SubstanceClass::Biological, CompositionKind::Composite, nullptr,
        true, false, false},
    {SUBSTANCE_STONE, "stone", "stone", "stone",
        SubstanceClass::Mixture, CompositionKind::Mixture, nullptr,
        true, false, false},
    {SUBSTANCE_GLASS, "glass", "glass", "glass",
        SubstanceClass::Mixture, CompositionKind::Mixture, nullptr,
        true, true, false},
    {SUBSTANCE_METAL, "metal", "metal", "metal",
        SubstanceClass::Mixture, CompositionKind::Mixture, nullptr,
        true, true, false},
    {SUBSTANCE_AIR,   "air",   "Air", "ins_mat_air",
        SubstanceClass::Mixture, CompositionKind::Mixture, nullptr,
        false, false, true},
};

SubstanceId substanceFromInternalName(char const *name) {
    if (!name || !name[0]) return SUBSTANCE_NONE;
    for (SubstanceId id = 0; id < SUBSTANCE_COUNT; ++id) {
        char const *internal = kBuiltinSubstances[id].internalName;
        if (internal && std::strcmp(internal, name) == 0) return id;
    }
    return SUBSTANCE_NONE;
}
