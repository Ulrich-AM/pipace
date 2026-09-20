#pragma once

#include "chemistry/SaceTypes.h"
#include "chemistry/SaceProperties.h"
#include "chemistry/SaceMolecule.h"
#include "chemistry/SaceDescriptors.h"
#include "chemistry/SaceFunctional.h"
#include "substance/SubstanceTypes.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <string>
#include <unordered_map>

// SACE Phase 2/3: in-memory generated identity catalog + identity-derived properties.
// SaceRecordId is not a SubstanceId. Generated records are not spawnable
// and are not live world matter. Session-local until persistence exists.
// Identity-defining fields are immutable after insert. Properties may be
// unknown even when identity is exact. Do not call from physics ticks.

using SaceRecordId = uint32_t;
constexpr SaceRecordId kSaceRecordNone = 0;

enum class SaceResolutionKind : uint8_t {
    Invalid = 0,
    BuiltIn,
    Generated
};

struct SaceSubstanceRef {
    SaceResolutionKind kind = SaceResolutionKind::Invalid;
    SubstanceId builtinId = SUBSTANCE_NONE;
    SaceRecordId generatedId = kSaceRecordNone;
};

struct SaceGeneratedRecord {
    SaceRecordId recordId = kSaceRecordNone;
    ChemicalRepresentationKind kind = ChemicalRepresentationKind::Unknown;
    ElementalComposition elemental{};
    std::string formula;
    std::string structureKey;
    std::string canonicalSignature;
    uint32_t displayOrdinal = 0;
    bool exactIdentity = false;
    bool spawnable = false;
    SaceGeneratedProperties properties{};
    bool hasMolecularGraph = false;
    SaceMolecularGraph molecularGraph{};
    std::string molecularGraphStructureKey;
    bool hasMolecularDescriptors = false;
    SaceMolecularDescriptors molecularDescriptors{};
    bool hasFunctionalProfile = false;
    SaceFunctionalProfile functionalProfile{};
};

// Player-facing display only. Uppercase hex, never part of canonical identity.
bool formatGeneratedElementName(uint32_t displayOrdinal, char *out, int cap);

class SaceCatalog {
public:
    void clear();
    std::size_t generatedCount() const { return records_.size(); }

    // Lookup-only when create == false. Built-in exact identities win.
    SaceSubstanceRef resolve(ChemicalIdentity const &query, bool create);

    // Pointer remains valid across later insertions. Invalid after clear().
    SaceGeneratedRecord const *record(SaceRecordId id) const;

    // Pointers in the returned view remain valid until that record is destroyed
    // (clear). Inserting other records does not invalidate them.
    ChemicalIdentity identityView(SaceGeneratedRecord const &rec) const;

    // Trusted provenance: claimedStructureKey must equal record.structureKey.
    // This is not graph-isomorphism proof. First exact attach wins.
    // Graph + binding key + descriptors + functional profile commit together.
    // Joback boiling-point estimate is attempted after local analysis; failure
    // does not fail attachment. Does not change canonical signature, record id,
    // or display ordinal. Molar mass stays identity-derived.
    bool attachMolecularGraph(SaceRecordId id, char const *claimedStructureKey,
        SaceMolecularGraph const &graph);

private:
    bool eligibleForGeneratedIdentity(ChemicalIdentity const &id, char *signature, int cap) const;
    SaceRecordId insertGenerated(ChemicalIdentity const &id, char const *signature);

    std::deque<SaceGeneratedRecord> records_;
    std::unordered_map<std::string, SaceRecordId> bySignature_;
    uint32_t nextDisplayOrdinal_ = 1;
};

SaceCatalog &saceGeneratedCatalog();

void runSaceCatalogDiagnostics();
