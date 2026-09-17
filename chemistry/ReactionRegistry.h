#pragma once

#include "chemistry/ReactionTypes.h"

// Built-in table lives in ReactionRegistry.cpp (function-local static), same
// pattern as SubstanceRegistry. No runtime registration yet.
// There is no ReactionEngine in this layer: definitions only.

ReactionDefinition const *builtinReactionTable();
int reactionTableSize(); // includes REACTION_NONE at slot 0
int reactionCount();     // defined reactions excluding NONE

inline bool validReaction(ReactionId id) {
    return id != REACTION_NONE && id < static_cast<ReactionId>(reactionTableSize());
}

inline ReactionDefinition const &reactionDef(ReactionId id) {
    ReactionDefinition const *table = builtinReactionTable();
    int n = reactionTableSize();
    if (id >= static_cast<ReactionId>(n)) return table[REACTION_NONE];
    return table[id];
}

ReactionId reactionFromInternalName(char const *name);
ReactionMassConservation reactionMassConservation(ReactionId id, float relativeTolerance = 0.02f);
