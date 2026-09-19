#include "chemistry/ReactionRegistry.h"

#include "substance/SubstanceTypes.h"

#include <algorithm>
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
    r.topology = ReactionTopology::HomogeneousCell;
    return r;
}

// 2 H2(g) + O2(g) -> 2 H2O(g)
//
// energyChangeJPerExtent is a simplified sandbox heat of reaction for the
// gas-water product (~483.6 kJ per mole of the written reaction). It is not a
// temperature-dependent thermochemical model. Do not substitute the
// liquid-water formation enthalpy here.
//
// minTemperatureK = 850 K is a practical ignition gate until reaction kinetics
// exist. It is NOT an Arrhenius model and NOT a universal autoignition
// temperature.
ReactionDefinition makeHydrogenCombustion() {
    ReactionDefinition r;
    r.id = REACTION_HYDROGEN_COMBUSTION;
    r.internalName = "hydrogen_combustion";
    r.reactants[0] = {SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0f};
    r.reactants[1] = {SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0f};
    r.reactantCount = 2;
    r.products[0] = {SUBSTANCE_WATER, MatterPhase::Gas, 2.0f};
    r.productCount = 1;
    r.conditions.minTemperatureValid = true;
    r.conditions.minTemperatureK = 850.0f;
    r.energyChangeJPerExtent = -483600.0f;
    r.topology = ReactionTopology::HomogeneousCell;
    // 1.0 mol/s/cell was a first-pass correctness cap. A 1-atm stoichiometric
    // H2/O2 cell holds ~0.2 mol of written reaction, so 1 mol/s burned for
    // several ticks before neighbors reached 850 K by conduction. 4 mol/s/cell
    // still cannot eat a large cloud in one 30 Hz tick (~0.13 mol/tick) while
    // dumping heat fast enough for thermal transport to ignite adjacent mix.
    // Ignition threshold stays 850 K.
    r.maxExtentPerSecond = 4.0f;
    return r;
}

// C(s) + O2(g) -> CO2(g)
//
// energyChangeJPerExtent is simplified graphite-to-CO2 formation enthalpy
// (~393.5 kJ per mole of the written reaction). Not temperature-dependent.
//
// minTemperatureK = 900 K is a practical sandbox ignition gate, not a
// universal carbon autoignition temperature. Not Arrhenius kinetics.
//
// maxExtentPerSecond = 120 mol/s per exposed 4-neighbor interface face.
// A full Carbon pixel at 4 cells/m is ~2940 mol, so one face consumes it
// in tens of seconds: visible erosion, not an instant vanish. Multiple
// exposed faces burn faster. This is a sandbox rate, not a graphite burn
// model. Ambient SUBSTANCE_AIR is not oxygen and does not participate.
ReactionDefinition makeCarbonCombustion() {
    ReactionDefinition r;
    r.id = REACTION_CARBON_COMBUSTION;
    r.internalName = "carbon_combustion";
    r.reactants[0] = {SUBSTANCE_CARBON, MatterPhase::Solid, 1.0f};
    r.reactants[1] = {SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0f};
    r.reactantCount = 2;
    r.products[0] = {SUBSTANCE_CARBON_DIOXIDE, MatterPhase::Gas, 1.0f};
    r.productCount = 1;
    r.conditions.minTemperatureValid = true;
    r.conditions.minTemperatureK = 900.0f;
    r.energyChangeJPerExtent = -393500.0f;
    r.topology = ReactionTopology::SolidGasSurface;
    r.maxExtentPerSecond = 120.0f;
    return r;
}

constexpr int kReactionTableSize = 3;

} // namespace

ReactionDefinition const *builtinReactionTable() {
    static ReactionDefinition const table[kReactionTableSize] = {
        makeNoneReaction(),
        makeHydrogenCombustion(),
        makeCarbonCombustion(),
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
