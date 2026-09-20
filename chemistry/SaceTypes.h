#pragma once

#include <cstdint>
#include <cstring>

// SACE Phase 1: exact chemical identity and elemental composition.
// Properties may be approximate; identity/composition stay as exact as practical.
// Formula text is display/reference only. Canonical identity uses kind,
// elemental counts, and an optional structure discriminator.

using AtomicNumber = uint8_t;

constexpr AtomicNumber kAtomicNone = 0;
constexpr AtomicNumber kAtomicHydrogen = 1;
constexpr AtomicNumber kAtomicCarbon = 6;
constexpr AtomicNumber kAtomicOxygen = 8;
constexpr AtomicNumber kMinAtomicNumber = 1;
constexpr AtomicNumber kMaxAtomicNumber = 118;

// Small-molecule / near-term compound capacity. 8 covers C/H/O plus a few
// extra elements without heap. Overflow is an error, never silent truncation.
constexpr int kMaxElementalSpecies = 8;
constexpr int kSaceSignatureCap = 192;

inline bool validAtomicNumber(AtomicNumber z) {
    return z >= kMinAtomicNumber && z <= kMaxAtomicNumber;
}

enum class ChemicalRepresentationKind : uint8_t {
    Unknown = 0,
    AtomicSpecies,
    SmallMolecule,
    IonicMaterial,
    Mixture,
    Polymer,
    NetworkSolid,
    Composite,
    GeneratedUnknown
};

enum class CompositionNormalizeResult : uint8_t {
    Ok = 0,
    Invalid,
    Overflow
};

enum class ConservationCheck : uint8_t {
    Unknown = 0,
    Balanced,
    Unbalanced
};

struct ElementCount {
    AtomicNumber atomicNumber = kAtomicNone;
    uint16_t count = 0;
};

struct ElementalComposition {
    ElementCount entries[kMaxElementalSpecies]{};
    uint8_t count = 0;
};

// Wide inventory for stoichiometric atom totals. Diagnostic / registration-time.
constexpr int kMaxElementalInventory = 16;

struct ElementInventoryEntry {
    AtomicNumber atomicNumber = kAtomicNone;
    int64_t count = 0;
};

struct ElementalInventory {
    ElementInventoryEntry entries[kMaxElementalInventory]{};
    uint8_t count = 0;
};

// Exact identity record. Not ChemicalProperties (molar mass / flammability).
// Not MatterPhase. Not a player-facing display name.
struct ChemicalIdentity {
    ChemicalRepresentationKind kind = ChemicalRepresentationKind::Unknown;
    ElementalComposition elemental{};
    char const *formula = nullptr;
    char const *structureKey = nullptr;
    bool exact = false;
};

inline char const *chemicalRepresentationKindKey(ChemicalRepresentationKind kind) {
    switch (kind) {
        case ChemicalRepresentationKind::AtomicSpecies: return "atomic";
        case ChemicalRepresentationKind::SmallMolecule: return "molecule";
        case ChemicalRepresentationKind::IonicMaterial: return "ionic";
        case ChemicalRepresentationKind::Mixture: return "mixture";
        case ChemicalRepresentationKind::Polymer: return "polymer";
        case ChemicalRepresentationKind::NetworkSolid: return "network";
        case ChemicalRepresentationKind::Composite: return "composite";
        case ChemicalRepresentationKind::GeneratedUnknown: return "generated";
        case ChemicalRepresentationKind::Unknown: return "unknown";
    }
    return "unknown";
}

inline char const *conservationCheckKey(ConservationCheck c) {
    switch (c) {
        case ConservationCheck::Balanced: return "balanced";
        case ConservationCheck::Unbalanced: return "unbalanced";
        case ConservationCheck::Unknown: return "unknown";
    }
    return "unknown";
}

inline uint16_t elementalCountOf(ElementalComposition const &c, AtomicNumber z) {
    if (!validAtomicNumber(z)) return 0;
    for (int n = 0; n < c.count && n < kMaxElementalSpecies; ++n) {
        if (c.entries[n].atomicNumber == z)
            return c.entries[n].count;
    }
    return 0;
}

inline bool elementalCompositionsEqual(ElementalComposition const &a, ElementalComposition const &b) {
    if (a.count != b.count) return false;
    for (int n = 0; n < a.count && n < kMaxElementalSpecies; ++n) {
        if (a.entries[n].atomicNumber != b.entries[n].atomicNumber) return false;
        if (a.entries[n].count != b.entries[n].count) return false;
    }
    return true;
}

inline char const *saceStructureKeyOrEmpty(char const *key) {
    return (key && key[0]) ? key : "";
}

inline bool saceStructureKeysEqual(char const *a, char const *b) {
    return std::strcmp(saceStructureKeyOrEmpty(a), saceStructureKeyOrEmpty(b)) == 0;
}

inline bool chemicalIdentityExact(ChemicalIdentity const &id) {
    return id.exact && id.elemental.count > 0
        && (id.kind == ChemicalRepresentationKind::AtomicSpecies
            || id.kind == ChemicalRepresentationKind::SmallMolecule
            || id.kind == ChemicalRepresentationKind::IonicMaterial);
}
