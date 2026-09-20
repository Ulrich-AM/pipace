#pragma once

#include "chemistry/SaceTypes.h"
#include "chemistry/ReactionTypes.h"

// SACE Phase 1 identity helpers. Registry / diagnostic / registration-time.
// Do not call from per-cell ReactionEngine ticks.

CompositionNormalizeResult normalizeElementalComposition(ElementalComposition &c);

bool addScaledComposition(ElementalInventory &acc, ElementalComposition const &src, int64_t scale);
bool elementalInventoriesEqual(ElementalInventory const &a, ElementalInventory const &b);

// Deterministic signature. Does not include formula text or display names.
// Format: kind|Z:count,Z:count|structureKey
bool writeChemicalSignature(ChemicalIdentity const &id, char *out, int cap);

bool chemicalIdentitiesEquivalent(ChemicalIdentity const &a, ChemicalIdentity const &b);

// Resolves an exact identity to a built-in SubstanceId. NONE if no exact match.
// Does not allocate a new SubstanceId.
SubstanceId findBuiltInByChemicalIdentity(ChemicalIdentity const &query);

ConservationCheck reactionAtomConservation(ReactionDefinition const &def);

ChemicalIdentity saceUnknownIdentity(ChemicalRepresentationKind kind = ChemicalRepresentationKind::Unknown);
ChemicalIdentity saceBuiltinWaterIdentity();
ChemicalIdentity saceBuiltinHydrogenIdentity();
ChemicalIdentity saceBuiltinOxygenIdentity();
ChemicalIdentity saceBuiltinCarbonIdentity();
ChemicalIdentity saceBuiltinCarbonDioxideIdentity();

// Headless: --sace-identity-diag -> misc/sace_identity_diag.tsv
void runSaceIdentityDiagnostics();
