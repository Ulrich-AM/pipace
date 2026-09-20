#pragma once

#include "chemistry/SaceMolecule.h"

#include <cstdint>

// SACE Phase 5: structure-derived descriptors from a valid molecular graph.
// NOT canonical identity. NOT physical-property estimates.
// Do not put these into canonicalSignature or SaceCatalog cache keys.
// Do not call from physics ticks.

constexpr int kMaxSaceBondClasses = 32;

struct SaceBondClass {
    AtomicNumber elementA = kAtomicNone; // min(Z,Z)
    AtomicNumber elementB = kAtomicNone; // max(Z,Z)
    SaceBondOrder order = SaceBondOrder::Single;
    uint16_t count = 0;
};

struct SaceMolecularDescriptors {
    int atomCount = 0;
    int heavyAtomCount = 0;     // non-hydrogen
    int hydrogenCount = 0;
    int heteroAtomCount = 0;    // neither H nor C (organic-style convention)

    int bondCount = 0;
    int singleBondCount = 0;
    int doubleBondCount = 0;
    int tripleBondCount = 0;
    int aromaticBondCount = 0;

    int branchAtomCount = 0;    // degree >= 3
    int terminalAtomCount = 0;  // degree == 1
    int cycleRank = 0;          // E - V + 1 for a connected graph

    int32_t netFormalCharge = 0;

    SaceBondClass bondClasses[kMaxSaceBondClasses]{};
    uint8_t bondClassCount = 0;
};

bool deriveMolecularDescriptors(SaceMolecularGraph const &graph, SaceMolecularDescriptors &out);
bool molecularDescriptorsEqual(SaceMolecularDescriptors const &a, SaceMolecularDescriptors const &b);
int countBondsBetweenElements(SaceMolecularDescriptors const &d,
    AtomicNumber a, AtomicNumber b, SaceBondOrder order);

void runSaceDescriptorDiagnostics();
