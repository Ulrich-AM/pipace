#include "chemistry/ReactionRegistry.h"

#include "substance/SubstanceTypes.h"

#include <cmath>
#include <cstring>

namespace {

ReactionDefinition makeNoneReaction() {
    ReactionDefinition r;
    r.id = REACTION_NONE;
    r.internalName = "none";
    r.reactantCount = 0;
    r.productCount = 0;
    r.energyChangeJPerExtent = 0.0f;
    return r;
}

// No scientifically defensible chemical reaction exists among current built-ins
// (Water, Honey, Wood, Stone, Glass, Metal, Air). Do not invent H2/O2, wood
// ash, or honey combustion here. Slot 0 is the empty sentinel only.

constexpr int kReactionTableSize = 1;

} // namespace

ReactionDefinition const *builtinReactionTable() {
    static ReactionDefinition const table[kReactionTableSize] = {
        makeNoneReaction(),
    };
    return table;
}

int reactionTableSize() {
    return kReactionTableSize;
}

int reactionCount() {
    int n = 0;
    ReactionDefinition const *table = builtinReactionTable();
    for (int i = 0; i < kReactionTableSize; ++i)
        if (table[i].id != REACTION_NONE) ++n;
    return n;
}

ReactionId reactionFromInternalName(char const *name) {
    if (!name || !name[0]) return REACTION_NONE;
    ReactionDefinition const *table = builtinReactionTable();
    for (int i = 0; i < kReactionTableSize; ++i) {
        char const *internal = table[i].internalName;
        if (internal && std::strcmp(internal, name) == 0) return table[i].id;
    }
    return REACTION_NONE;
}

namespace {

bool participantMolarMassKgMol(ReactionParticipant const &p, double &massKgPerMol) {
    if (!reactionParticipantUsed(p)) {
        massKgPerMol = 0.0;
        return true;
    }
    ChemicalProperties const &chem = chemicalForSubstance(p.substance);
    if (!chem.valid || !(chem.molarMass > 0.0f) || !std::isfinite(chem.molarMass))
        return false;
    massKgPerMol = static_cast<double>(p.coefficient) * (static_cast<double>(chem.molarMass) * 1.0e-3);
    return true;
}

} // namespace

ReactionMassConservation reactionMassConservation(ReactionId id, float relativeTolerance) {
    if (!validReaction(id)) return ReactionMassConservation::Unavailable;
    ReactionDefinition const &def = reactionDef(id);
    double react = 0.0;
    double prod = 0.0;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        double m = 0.0;
        if (!participantMolarMassKgMol(def.reactants[n], m))
            return ReactionMassConservation::Unavailable;
        react += m;
    }
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n) {
        double m = 0.0;
        if (!participantMolarMassKgMol(def.products[n], m))
            return ReactionMassConservation::Unavailable;
        prod += m;
    }
    if (!(react > 0.0) || !(prod > 0.0) || !std::isfinite(react) || !std::isfinite(prod))
        return ReactionMassConservation::Unavailable;
    double scale = std::max(react, prod);
    float tol = (std::isfinite(relativeTolerance) && relativeTolerance > 0.0f)
        ? relativeTolerance : 0.02f;
    if (std::abs(react - prod) <= static_cast<double>(tol) * scale)
        return ReactionMassConservation::Balanced;
    return ReactionMassConservation::Imbalanced;
}
