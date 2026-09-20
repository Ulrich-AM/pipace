#pragma once

#include "chemistry/SaceTypes.h"
#include "substance/SubstanceTypes.h"

#include <cstdint>
#include <vector>

// SACE Phase 4: structured small-molecule graphs.
// Metadata only. Not world physics. Not canonical identity.
// structureKey remains the structural discriminator until graph
// isomorphism canonicalization exists. Debug serialization is NOT
// a catalog cache key. Do not call from physics ticks.

constexpr int kMaxSaceMoleculeAtoms = 64;
constexpr int kMaxSaceMoleculeBonds = 96;

enum class SaceBondOrder : uint8_t {
    Single = 0,
    Double,
    Triple,
    Aromatic
};

struct SaceAtom {
    AtomicNumber atomicNumber = kAtomicNone;
    int8_t formalCharge = 0;
};

struct SaceBond {
    uint16_t atomA = 0;
    uint16_t atomB = 0;
    SaceBondOrder order = SaceBondOrder::Single;
};

struct SaceMolecularGraph {
    std::vector<SaceAtom> atoms;
    std::vector<SaceBond> bonds;
};

enum class SaceGraphValidation : uint8_t {
    Ok = 0,
    Empty,
    TooManyAtoms,
    TooManyBonds,
    InvalidAtom,
    InvalidBond,
    SelfBond,
    DuplicateBond,
    Disconnected,
    CompositionOverflow
};

inline char const *saceGraphValidationKey(SaceGraphValidation v) {
    switch (v) {
        case SaceGraphValidation::Ok: return "ok";
        case SaceGraphValidation::Empty: return "empty";
        case SaceGraphValidation::TooManyAtoms: return "too-many-atoms";
        case SaceGraphValidation::TooManyBonds: return "too-many-bonds";
        case SaceGraphValidation::InvalidAtom: return "invalid-atom";
        case SaceGraphValidation::InvalidBond: return "invalid-bond";
        case SaceGraphValidation::SelfBond: return "self-bond";
        case SaceGraphValidation::DuplicateBond: return "duplicate-bond";
        case SaceGraphValidation::Disconnected: return "disconnected";
        case SaceGraphValidation::CompositionOverflow: return "composition-overflow";
    }
    return "invalid";
}

inline bool validSaceBondOrder(SaceBondOrder order) {
    switch (order) {
        case SaceBondOrder::Single:
        case SaceBondOrder::Double:
        case SaceBondOrder::Triple:
        case SaceBondOrder::Aromatic:
            return true;
    }
    return false;
}

SaceGraphValidation validateMolecularGraph(SaceMolecularGraph const &graph, bool requireConnected);

bool elementalCompositionFromGraph(SaceMolecularGraph const &graph, ElementalComposition &out);
bool graphMatchesChemicalIdentity(SaceMolecularGraph const &graph, ChemicalIdentity const &id);
bool molecularGraphsStoredEqual(SaceMolecularGraph const &a, SaceMolecularGraph const &b);

int32_t molecularGraphNetFormalCharge(SaceMolecularGraph const &graph);
int atomCount(SaceMolecularGraph const &graph);
int bondCount(SaceMolecularGraph const &graph);
int atomDegree(SaceMolecularGraph const &graph, int atomIndex);
int countElement(SaceMolecularGraph const &graph, AtomicNumber z);
bool hasBond(SaceMolecularGraph const &graph, int atomA, int atomB);

// Debug serialization of the CURRENT atom/bond array order.
// NOT canonical chemical identity. Do not use as a catalog cache key.
bool writeMolecularGraphDebug(SaceMolecularGraph const &graph, char *out, int cap);

// Advisory common-valence check for H/C/O only. Not universal chemistry.
bool saceCommonValenceLooksTypical(SaceMolecularGraph const &graph);

bool saceBuiltinMolecularGraph(SubstanceId id, SaceMolecularGraph &out);
// Associated built-in structureKey from ChemicalIdentity. nullptr if no graph.
char const *saceBuiltinMolecularStructureKey(SubstanceId id);
bool saceSyntheticC2H6OGraphA(SaceMolecularGraph &out); // C-C-O skeleton, explicit H
bool saceSyntheticC2H6OGraphB(SaceMolecularGraph &out); // C-O-C skeleton, explicit H

void runSaceMoleculeDiagnostics();
