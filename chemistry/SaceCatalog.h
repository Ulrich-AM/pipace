#pragma once

#include "chemistry/SaceTypes.h"
#include "substance/SubstanceTypes.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// SACE Phase 2: in-memory generated identity catalog.
// SaceRecordId is not a SubstanceId. Generated records are not spawnable
// and are not live world matter. Session-local until persistence exists.
// Do not call from physics ticks.

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
};

// Player-facing display only. Uppercase hex, never part of canonical identity.
bool formatGeneratedElementName(uint32_t displayOrdinal, char *out, int cap);

class SaceCatalog {
public:
    void clear();
    std::size_t generatedCount() const { return records_.size(); }

    // Lookup-only when create == false. Built-in exact identities win.
    SaceSubstanceRef resolve(ChemicalIdentity const &query, bool create);

    SaceGeneratedRecord const *record(SaceRecordId id) const;

    // Pointers in the returned view are valid only until the next catalog mutation.
    ChemicalIdentity identityView(SaceGeneratedRecord const &rec) const;

private:
    bool eligibleForGeneratedIdentity(ChemicalIdentity const &id, char *signature, int cap) const;
    SaceRecordId insertGenerated(ChemicalIdentity const &id, char const *signature);

    std::vector<SaceGeneratedRecord> records_;
    std::unordered_map<std::string, SaceRecordId> bySignature_;
    uint32_t nextDisplayOrdinal_ = 1;
};

SaceCatalog &saceGeneratedCatalog();

void runSaceCatalogDiagnostics();
