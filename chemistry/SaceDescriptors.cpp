#include "chemistry/SaceDescriptors.h"

#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>

namespace {

bool isHydrogen(AtomicNumber z) { return z == kAtomicHydrogen; }
bool isCarbon(AtomicNumber z) { return z == kAtomicCarbon; }

int bondClassCmp(SaceBondClass const &a, SaceBondClass const &b) {
    if (a.elementA != b.elementA) return a.elementA < b.elementA ? -1 : 1;
    if (a.elementB != b.elementB) return a.elementB < b.elementB ? -1 : 1;
    uint8_t oa = static_cast<uint8_t>(a.order);
    uint8_t ob = static_cast<uint8_t>(b.order);
    if (oa != ob) return oa < ob ? -1 : 1;
    return 0;
}

bool addBondClass(SaceMolecularDescriptors &d, AtomicNumber za, AtomicNumber zb, SaceBondOrder order) {
    AtomicNumber a = za < zb ? za : zb;
    AtomicNumber b = za < zb ? zb : za;
    for (int i = 0; i < d.bondClassCount; ++i) {
        SaceBondClass &c = d.bondClasses[i];
        if (c.elementA == a && c.elementB == b && c.order == order) {
            uint32_t n = static_cast<uint32_t>(c.count) + 1u;
            if (n > 0xFFFFu) return false;
            c.count = static_cast<uint16_t>(n);
            return true;
        }
    }
    if (d.bondClassCount >= kMaxSaceBondClasses)
        return false;
    SaceBondClass &slot = d.bondClasses[d.bondClassCount];
    slot.elementA = a;
    slot.elementB = b;
    slot.order = order;
    slot.count = 1;
    ++d.bondClassCount;
    return true;
}

void sortBondClasses(SaceMolecularDescriptors &d) {
    for (int i = 0; i < d.bondClassCount; ++i) {
        int best = i;
        for (int j = i + 1; j < d.bondClassCount; ++j)
            if (bondClassCmp(d.bondClasses[j], d.bondClasses[best]) < 0) best = j;
        if (best != i) {
            SaceBondClass tmp = d.bondClasses[i];
            d.bondClasses[i] = d.bondClasses[best];
            d.bondClasses[best] = tmp;
        }
    }
}

} // namespace

int countBondsBetweenElements(SaceMolecularDescriptors const &d,
    AtomicNumber a, AtomicNumber b, SaceBondOrder order)
{
    AtomicNumber lo = a < b ? a : b;
    AtomicNumber hi = a < b ? b : a;
    for (int i = 0; i < d.bondClassCount && i < kMaxSaceBondClasses; ++i) {
        SaceBondClass const &c = d.bondClasses[i];
        if (c.elementA == lo && c.elementB == hi && c.order == order)
            return c.count;
    }
    return 0;
}

bool molecularDescriptorsEqual(SaceMolecularDescriptors const &a, SaceMolecularDescriptors const &b) {
    if (a.atomCount != b.atomCount) return false;
    if (a.heavyAtomCount != b.heavyAtomCount) return false;
    if (a.hydrogenCount != b.hydrogenCount) return false;
    if (a.heteroAtomCount != b.heteroAtomCount) return false;
    if (a.bondCount != b.bondCount) return false;
    if (a.singleBondCount != b.singleBondCount) return false;
    if (a.doubleBondCount != b.doubleBondCount) return false;
    if (a.tripleBondCount != b.tripleBondCount) return false;
    if (a.aromaticBondCount != b.aromaticBondCount) return false;
    if (a.branchAtomCount != b.branchAtomCount) return false;
    if (a.terminalAtomCount != b.terminalAtomCount) return false;
    if (a.cycleRank != b.cycleRank) return false;
    if (a.netFormalCharge != b.netFormalCharge) return false;
    if (a.bondClassCount != b.bondClassCount) return false;
    for (int i = 0; i < a.bondClassCount; ++i) {
        SaceBondClass const &ca = a.bondClasses[i];
        SaceBondClass const &cb = b.bondClasses[i];
        if (ca.elementA != cb.elementA || ca.elementB != cb.elementB) return false;
        if (ca.order != cb.order || ca.count != cb.count) return false;
    }
    return true;
}

bool deriveMolecularDescriptors(SaceMolecularGraph const &graph, SaceMolecularDescriptors &out) {
    out = {};
    if (validateMolecularGraph(graph, true) != SaceGraphValidation::Ok)
        return false;
    SaceMolecularDescriptors d{};
    int v = atomCount(graph);
    int e = bondCount(graph);
    d.atomCount = v;
    d.bondCount = e;
    d.netFormalCharge = molecularGraphNetFormalCharge(graph);
    int cycle = e - v + 1;
    if (cycle < 0)
        return false;
    d.cycleRank = cycle;

    for (int i = 0; i < v; ++i) {
        AtomicNumber z = graph.atoms[static_cast<size_t>(i)].atomicNumber;
        if (isHydrogen(z)) ++d.hydrogenCount;
        else ++d.heavyAtomCount;
        if (!isHydrogen(z) && !isCarbon(z)) ++d.heteroAtomCount;
        int deg = atomDegree(graph, i);
        if (deg == 1) ++d.terminalAtomCount;
        if (deg >= 3) ++d.branchAtomCount;
    }

    for (SaceBond const &bond : graph.bonds) {
        switch (bond.order) {
            case SaceBondOrder::Single: ++d.singleBondCount; break;
            case SaceBondOrder::Double: ++d.doubleBondCount; break;
            case SaceBondOrder::Triple: ++d.tripleBondCount; break;
            case SaceBondOrder::Aromatic: ++d.aromaticBondCount; break;
        }
        AtomicNumber za = graph.atoms[bond.atomA].atomicNumber;
        AtomicNumber zb = graph.atoms[bond.atomB].atomicNumber;
        if (!addBondClass(d, za, zb, bond.order))
            return false;
    }
    sortBondClasses(d);
    out = d;
    return true;
}

namespace {

SaceMolecularGraph permuteWaterGraph(SaceMolecularGraph const &src) {
    // Original: O=0 H=1 H=2. Permuted: H=0 H=1 O=2 with the same H-O-H topology.
    (void)src;
    SaceMolecularGraph g{};
    g.atoms.push_back({kAtomicHydrogen, 0});
    g.atoms.push_back({kAtomicHydrogen, 0});
    g.atoms.push_back({kAtomicOxygen, 0});
    g.bonds.push_back({2, 0, SaceBondOrder::Single});
    g.bonds.push_back({2, 1, SaceBondOrder::Single});
    return g;
}

SaceMolecularGraph reorderBonds(SaceMolecularGraph g) {
    if (g.bonds.size() >= 2)
        std::swap(g.bonds[0], g.bonds[g.bonds.size() - 1]);
    return g;
}

SaceMolecularGraph reverseBondEnds(SaceMolecularGraph g) {
    for (SaceBond &b : g.bonds) {
        uint16_t t = b.atomA;
        b.atomA = b.atomB;
        b.atomB = t;
    }
    return g;
}

bool nearC2H6O(ElementalComposition const &c) {
    return elementalCountOf(c, kAtomicCarbon) == 2
        && elementalCountOf(c, kAtomicHydrogen) == 6
        && elementalCountOf(c, kAtomicOxygen) == 1;
}

} // namespace

void runSaceDescriptorDiagnostics() {
    std::ofstream out(miscFile("sace_descriptor_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    SaceMolecularGraph empty{};
    SaceMolecularDescriptors junk{};
    junk.atomCount = 99;
    bool invalidFail = !deriveMolecularDescriptors(empty, junk)
        && junk.atomCount == 0 && junk.bondClassCount == 0;
    emit("invalid_graph_fails_descriptor_derivation", invalidFail, "");

    SaceMolecularGraph h2{}, o2{}, h2o{}, co2{};
    SaceMolecularDescriptors dH2{}, dO2{}, dH2O{}, dCO2{};
    emit("h2_descriptors_correct",
        saceBuiltinMolecularGraph(SUBSTANCE_HYDROGEN, h2)
            && deriveMolecularDescriptors(h2, dH2)
            && dH2.atomCount == 2 && dH2.heavyAtomCount == 0 && dH2.hydrogenCount == 2
            && dH2.bondCount == 1 && dH2.cycleRank == 0
            && dH2.singleBondCount == 1, "");
    emit("o2_descriptors_correct",
        saceBuiltinMolecularGraph(SUBSTANCE_OXYGEN, o2)
            && deriveMolecularDescriptors(o2, dO2)
            && dO2.atomCount == 2 && dO2.heavyAtomCount == 2 && dO2.hydrogenCount == 0
            && dO2.doubleBondCount == 1 && dO2.cycleRank == 0, "");
    emit("h2o_descriptors_correct",
        saceBuiltinMolecularGraph(SUBSTANCE_WATER, h2o)
            && deriveMolecularDescriptors(h2o, dH2O)
            && dH2O.atomCount == 3 && dH2O.heavyAtomCount == 1 && dH2O.hydrogenCount == 2
            && dH2O.heteroAtomCount == 1 && dH2O.bondCount == 2 && dH2O.cycleRank == 0, "");
    emit("co2_descriptors_correct",
        saceBuiltinMolecularGraph(SUBSTANCE_CARBON_DIOXIDE, co2)
            && deriveMolecularDescriptors(co2, dCO2)
            && dCO2.atomCount == 3 && dCO2.heavyAtomCount == 3 && dCO2.hydrogenCount == 0
            && dCO2.heteroAtomCount == 2 && dCO2.doubleBondCount == 2 && dCO2.cycleRank == 0, "");
    emit("water_oh_single_count_2",
        countBondsBetweenElements(dH2O, kAtomicOxygen, kAtomicHydrogen, SaceBondOrder::Single) == 2, "");
    emit("o2_oo_double_count_1",
        countBondsBetweenElements(dO2, kAtomicOxygen, kAtomicOxygen, SaceBondOrder::Double) == 1, "");
    emit("co2_co_double_count_2",
        countBondsBetweenElements(dCO2, kAtomicCarbon, kAtomicOxygen, SaceBondOrder::Double) == 2, "");

    SaceMolecularGraph gA{}, gB{};
    SaceMolecularDescriptors dA{}, dB{};
    emit("c2h6o_a_descriptors_valid",
        saceSyntheticC2H6OGraphA(gA) && deriveMolecularDescriptors(gA, dA), "");
    emit("c2h6o_b_descriptors_valid",
        saceSyntheticC2H6OGraphB(gB) && deriveMolecularDescriptors(gB, dB), "");
    ElementalComposition cA{}, cB{};
    emit("same_elemental_composition",
        elementalCompositionFromGraph(gA, cA) && elementalCompositionFromGraph(gB, cB)
            && nearC2H6O(cA) && elementalCompositionsEqual(cA, cB), "");
    emit("a_cc_single_1",
        countBondsBetweenElements(dA, kAtomicCarbon, kAtomicCarbon, SaceBondOrder::Single) == 1, "");
    emit("a_co_single_1",
        countBondsBetweenElements(dA, kAtomicCarbon, kAtomicOxygen, SaceBondOrder::Single) == 1, "");
    emit("a_oh_single_1",
        countBondsBetweenElements(dA, kAtomicOxygen, kAtomicHydrogen, SaceBondOrder::Single) == 1, "");
    emit("b_cc_single_0",
        countBondsBetweenElements(dB, kAtomicCarbon, kAtomicCarbon, SaceBondOrder::Single) == 0, "");
    emit("b_co_single_2",
        countBondsBetweenElements(dB, kAtomicCarbon, kAtomicOxygen, SaceBondOrder::Single) == 2, "");
    emit("b_oh_single_0",
        countBondsBetweenElements(dB, kAtomicOxygen, kAtomicHydrogen, SaceBondOrder::Single) == 0, "");
    emit("both_cycle_rank_0", dA.cycleRank == 0 && dB.cycleRank == 0, "");

    SaceMolecularGraph waterPerm = permuteWaterGraph(h2o);
    SaceMolecularDescriptors dPerm{};
    emit("atom_reorder_preserves_descriptors",
        deriveMolecularDescriptors(waterPerm, dPerm)
            && molecularDescriptorsEqual(dH2O, dPerm), "");

    SaceMolecularDescriptors dBondOrder{}, dEnds{};
    emit("bond_reorder_preserves_descriptors",
        deriveMolecularDescriptors(reorderBonds(h2o), dBondOrder)
            && molecularDescriptorsEqual(dH2O, dBondOrder), "");
    emit("reversing_bond_endpoints_preserves_descriptors",
        deriveMolecularDescriptors(reverseBondEnds(co2), dEnds)
            && molecularDescriptorsEqual(dCO2, dEnds), "");

    char sig[kSaceSignatureCap]{};
    bool wrote = writeChemicalSignature(saceBuiltinWaterIdentity(), sig, kSaceSignatureCap);
    emit("descriptors_not_in_canonical_signature",
        wrote && std::strstr(sig, "heavy") == nullptr
            && std::strstr(sig, "cycle") == nullptr
            && std::strstr(sig, "descriptor") == nullptr
            && std::strstr(sig, "debug|") == nullptr, sig);

    SaceMolecularDescriptors leftover{};
    leftover.atomCount = 7;
    leftover.bondClassCount = 3;
    leftover.bondClasses[0].count = 9;
    bool failClears = !deriveMolecularDescriptors(empty, leftover)
        && leftover.atomCount == 0 && leftover.bondClassCount == 0
        && leftover.bondClasses[0].count == 0;
    emit("descriptor_failure_resets_output", failClears, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
