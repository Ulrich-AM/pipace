#pragma once

#include "substance/SubstanceTypes.h"

#include <cstdint>

// Reaction identity only. Not progress, temperature, or local world state.
using ReactionId = uint16_t;

constexpr ReactionId REACTION_NONE = 0;

// Fixed participant slots. Covers 2 H2 + O2 -> 2 H2O with room for extras.
// No per-cell or per-definition heap.
constexpr int kMaxReactionParticipants = 8;

// One participant: a substance in an optional required phase, with a
// stoichiometric coefficient in moles per reaction extent. Coefficients need
// not sum to 1. MatterPhase::None means any supported phase.
// Same SubstanceId may appear as liquid and as gas in different slots.
struct ReactionParticipant {
    SubstanceId substance = SUBSTANCE_NONE;
    MatterPhase requiredPhase = MatterPhase::None;
    float coefficient = 0.0f; // stoichiometric moles of this participant per extent
};

// Declarative bounds. Unset flags mean unbounded (no fake 0 K / 0 Pa limits).
// catalyst == SUBSTANCE_NONE means no catalyst is required.
struct ReactionConditions {
    bool minTemperatureValid = false;
    float minTemperatureK = 0.0f;
    bool maxTemperatureValid = false;
    float maxTemperatureK = 0.0f;
    bool minPressureValid = false;
    float minPressurePa = 0.0f;
    bool maxPressureValid = false;
    float maxPressurePa = 0.0f;
    SubstanceId catalyst = SUBSTANCE_NONE;
};

// energyChangeJPerExtent is ΔH-like heat of reaction per unit extent:
//   < 0  exothermic (releases heat to the thermal system)
//   > 0  endothermic (consumes heat from the thermal system)
//   = 0  thermally neutral / unspecified
//
// Extent is one mole of the written stoichiometric reaction. Example:
//   2 H2 + O2 -> 2 H2O
// extent = 1 means consume 2 mol H2, 1 mol O2, produce 2 mol H2O.
// energyChangeJPerExtent is then J per mole of that reaction as written.
// Diagnostic synthetic reactions may use a separate unit path (not moles).
// This layer does not move heat; ReactionEngine applies the sign later.
struct ReactionDefinition {
    ReactionId id = REACTION_NONE;
    char const *internalName = "none";
    ReactionParticipant reactants[kMaxReactionParticipants]{};
    uint8_t reactantCount = 0;
    ReactionParticipant products[kMaxReactionParticipants]{};
    uint8_t productCount = 0;
    ReactionConditions conditions{};
    float energyChangeJPerExtent = 0.0f;
    // Authorable first-order cap. Not Arrhenius kinetics. <= 0 means the
    // engine default (see ReactionEngine). Extent 1 consumes the written
    // coefficients in one second at this rate when matter allows.
    float maxExtentPerSecond = 0.0f;
};

inline bool reactionParticipantUsed(ReactionParticipant const &p) {
    return p.substance != SUBSTANCE_NONE && p.coefficient > 0.0f;
}

inline bool reactionParticipantPhaseOk(ReactionParticipant const &p, MatterPhase actual) {
    if (p.requiredPhase == MatterPhase::None) return true;
    return p.requiredPhase == actual;
}

inline bool reactionConditionsMatch(ReactionConditions const &c, float temperatureK, float pressurePa) {
    if (c.minTemperatureValid && !(temperatureK >= c.minTemperatureK)) return false;
    if (c.maxTemperatureValid && !(temperatureK <= c.maxTemperatureK)) return false;
    if (c.minPressureValid && !(pressurePa >= c.minPressurePa)) return false;
    if (c.maxPressureValid && !(pressurePa <= c.maxPressurePa)) return false;
    return true;
}

// Heat into the world for a given extent. Positive = warming.
inline float reactionHeatReleasedJ(ReactionDefinition const &def, float extent) {
    return -def.energyChangeJPerExtent * extent;
}

enum class ReactionMassConservation : uint8_t {
    Unavailable = 0, // at least one participant lacks chemical.valid + molarMass
    Balanced,
    Imbalanced
};
