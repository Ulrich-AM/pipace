#pragma once

#include "chemistry/SaceMolecule.h"
#include "chemistry/SaceDescriptors.h"

#include <cstdint>

// SACE Phase 6: H/C/O functional motifs from a valid molecular graph.
// Graph-derived metadata only. NOT canonical identity. NOT physical properties.
// Motifs inspect actual connectivity, never formula or structureKey.
// Do not call from physics ticks.

enum class SaceHydrogenBondingClass : uint8_t {
    Unknown = 0,
    None,
    AcceptorOnly,
    DonorOnly,
    DonorAndAcceptor
};

inline char const *saceHydrogenBondingClassKey(SaceHydrogenBondingClass c) {
    switch (c) {
        case SaceHydrogenBondingClass::None: return "none";
        case SaceHydrogenBondingClass::AcceptorOnly: return "acceptor-only";
        case SaceHydrogenBondingClass::DonorOnly: return "donor-only";
        case SaceHydrogenBondingClass::DonorAndAcceptor: return "donor-and-acceptor";
        case SaceHydrogenBondingClass::Unknown: return "unknown";
    }
    return "unknown";
}

struct SaceFunctionalProfile {
    int hydroxylCount = 0;      // C-O-H oxygen sites; Water is not counted
    int etherOxygenCount = 0;   // C-O-C, no bonded H
    int carbonylCount = 0;      // number of C=O bonds, not a ketone classifier

    int hBondDonorCount = 0;    // oxygen donor sites (not each H)
    int hBondAcceptorCount = 0; // conservative H/C/O oxygen acceptors; O=O is not

    int carbonCarbonSingleCount = 0;
    int carbonCarbonDoubleCount = 0;
    int carbonCarbonTripleCount = 0;
    int carbonOxygenSingleCount = 0;
    int carbonOxygenDoubleCount = 0;
    int oxygenHydrogenSingleCount = 0;

    SaceHydrogenBondingClass hydrogenBondingClass = SaceHydrogenBondingClass::Unknown;
    bool supported = false;     // true only for H/C/O-only graphs
};

bool deriveFunctionalProfile(SaceMolecularGraph const &graph, SaceFunctionalProfile &out);
bool functionalProfilesEqual(SaceFunctionalProfile const &a, SaceFunctionalProfile const &b);

void runSaceFunctionalDiagnostics();
