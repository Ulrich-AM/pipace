#include "chemistry/SaceFunctional.h"

#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>

namespace {

bool isH(AtomicNumber z) { return z == kAtomicHydrogen; }
bool isC(AtomicNumber z) { return z == kAtomicCarbon; }
bool isO(AtomicNumber z) { return z == kAtomicOxygen; }
bool isHCO(AtomicNumber z) { return isH(z) || isC(z) || isO(z); }

SaceHydrogenBondingClass classifyHBond(int donors, int acceptors) {
    if (donors > 0 && acceptors > 0) return SaceHydrogenBondingClass::DonorAndAcceptor;
    if (donors > 0) return SaceHydrogenBondingClass::DonorOnly;
    if (acceptors > 0) return SaceHydrogenBondingClass::AcceptorOnly;
    return SaceHydrogenBondingClass::None;
}

void countOxygenNeighbors(SaceMolecularGraph const &g, int oIndex,
    int &hSingle, int &cSingle, int &cDouble)
{
    hSingle = 0;
    cSingle = 0;
    cDouble = 0;
    uint16_t idx = static_cast<uint16_t>(oIndex);
    for (SaceBond const &b : g.bonds) {
        int other = -1;
        if (b.atomA == idx) other = b.atomB;
        else if (b.atomB == idx) other = b.atomA;
        if (other < 0 || other >= atomCount(g)) continue;
        AtomicNumber z = g.atoms[static_cast<size_t>(other)].atomicNumber;
        if (b.order == SaceBondOrder::Single && isH(z)) ++hSingle;
        if (b.order == SaceBondOrder::Single && isC(z)) ++cSingle;
        if (b.order == SaceBondOrder::Double && isC(z)) ++cDouble;
    }
}

} // namespace

bool functionalProfilesEqual(SaceFunctionalProfile const &a, SaceFunctionalProfile const &b) {
    return a.hydroxylCount == b.hydroxylCount
        && a.etherOxygenCount == b.etherOxygenCount
        && a.carbonylCount == b.carbonylCount
        && a.hBondDonorCount == b.hBondDonorCount
        && a.hBondAcceptorCount == b.hBondAcceptorCount
        && a.carbonCarbonSingleCount == b.carbonCarbonSingleCount
        && a.carbonCarbonDoubleCount == b.carbonCarbonDoubleCount
        && a.carbonCarbonTripleCount == b.carbonCarbonTripleCount
        && a.carbonOxygenSingleCount == b.carbonOxygenSingleCount
        && a.carbonOxygenDoubleCount == b.carbonOxygenDoubleCount
        && a.oxygenHydrogenSingleCount == b.oxygenHydrogenSingleCount
        && a.hydrogenBondingClass == b.hydrogenBondingClass
        && a.supported == b.supported;
}

bool deriveFunctionalProfile(SaceMolecularGraph const &graph, SaceFunctionalProfile &out) {
    out = {};
    if (validateMolecularGraph(graph, true) != SaceGraphValidation::Ok)
        return false;

    bool onlyHCO = true;
    for (SaceAtom const &atom : graph.atoms) {
        if (!isHCO(atom.atomicNumber)) {
            onlyHCO = false;
            break;
        }
    }
    if (!onlyHCO) {
        out.supported = false;
        out.hydrogenBondingClass = SaceHydrogenBondingClass::Unknown;
        return true;
    }

    SaceMolecularDescriptors desc{};
    if (!deriveMolecularDescriptors(graph, desc))
        return false;

    SaceFunctionalProfile p{};
    p.supported = true;
    p.carbonCarbonSingleCount = countBondsBetweenElements(desc, kAtomicCarbon, kAtomicCarbon, SaceBondOrder::Single);
    p.carbonCarbonDoubleCount = countBondsBetweenElements(desc, kAtomicCarbon, kAtomicCarbon, SaceBondOrder::Double);
    p.carbonCarbonTripleCount = countBondsBetweenElements(desc, kAtomicCarbon, kAtomicCarbon, SaceBondOrder::Triple);
    p.carbonOxygenSingleCount = countBondsBetweenElements(desc, kAtomicCarbon, kAtomicOxygen, SaceBondOrder::Single);
    p.carbonOxygenDoubleCount = countBondsBetweenElements(desc, kAtomicCarbon, kAtomicOxygen, SaceBondOrder::Double);
    p.oxygenHydrogenSingleCount = countBondsBetweenElements(desc, kAtomicOxygen, kAtomicHydrogen, SaceBondOrder::Single);
    p.carbonylCount = p.carbonOxygenDoubleCount;

    int n = atomCount(graph);
    for (int i = 0; i < n; ++i) {
        if (!isO(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        int8_t charge = graph.atoms[static_cast<size_t>(i)].formalCharge;
        int hSingle = 0, cSingle = 0, cDouble = 0;
        countOxygenNeighbors(graph, i, hSingle, cSingle, cDouble);
        if (hSingle >= 1 && cSingle >= 1)
            ++p.hydroxylCount;
        if (charge == 0 && cSingle == 2 && hSingle == 0)
            ++p.etherOxygenCount;
        if (charge == 0 && hSingle >= 1)
            ++p.hBondDonorCount;
        // Conservative acceptor: neutral O bonded to C or H. O=O is excluded.
        if (charge == 0 && (cSingle + cDouble + hSingle) > 0)
            ++p.hBondAcceptorCount;
    }
    p.hydrogenBondingClass = classifyHBond(p.hBondDonorCount, p.hBondAcceptorCount);
    out = p;
    return true;
}

namespace {

SaceMolecularGraph permuteWater() {
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

SaceMolecularGraph reverseEnds(SaceMolecularGraph g) {
    for (SaceBond &b : g.bonds) {
        uint16_t t = b.atomA;
        b.atomA = b.atomB;
        b.atomB = t;
    }
    return g;
}

} // namespace

void runSaceFunctionalDiagnostics() {
    std::ofstream out(miscFile("sace_functional_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    SaceMolecularGraph empty{};
    SaceFunctionalProfile junk{};
    junk.hydroxylCount = 9;
    junk.supported = true;
    bool invalidFail = !deriveFunctionalProfile(empty, junk)
        && junk.hydroxylCount == 0 && !junk.supported;
    emit("invalid_graph_fails_and_resets", invalidFail, "");

    SaceMolecularGraph nitrogen{};
    nitrogen.atoms.push_back({7, 0});
    SaceFunctionalProfile nProf{};
    bool nOk = deriveFunctionalProfile(nitrogen, nProf);
    emit("unsupported_element_returns_unsupported",
        nOk && !nProf.supported
            && nProf.hydroxylCount == 0 && nProf.etherOxygenCount == 0
            && nProf.carbonylCount == 0
            && nProf.hydrogenBondingClass == SaceHydrogenBondingClass::Unknown, "");

    SaceMolecularGraph h2{}, o2{}, h2o{}, co2{}, gA{}, gB{};
    SaceFunctionalProfile pH2{}, pO2{}, pW{}, pCO2{}, pA{}, pB{};
    emit("h2_supported",
        saceBuiltinMolecularGraph(SUBSTANCE_HYDROGEN, h2)
            && deriveFunctionalProfile(h2, pH2) && pH2.supported, "");
    emit("h2_no_hbond",
        pH2.hBondDonorCount == 0 && pH2.hBondAcceptorCount == 0
            && pH2.hydrogenBondingClass == SaceHydrogenBondingClass::None, "");
    emit("o2_no_hbond",
        saceBuiltinMolecularGraph(SUBSTANCE_OXYGEN, o2)
            && deriveFunctionalProfile(o2, pO2) && pO2.supported
            && pO2.hBondDonorCount == 0 && pO2.hBondAcceptorCount == 0
            && pO2.hydrogenBondingClass == SaceHydrogenBondingClass::None, "");

    emit("water_hydroxyl_0",
        saceBuiltinMolecularGraph(SUBSTANCE_WATER, h2o)
            && deriveFunctionalProfile(h2o, pW) && pW.supported && pW.hydroxylCount == 0, "");
    emit("water_ether_0", pW.etherOxygenCount == 0, "");
    emit("water_donor_1", pW.hBondDonorCount == 1, "");
    emit("water_acceptor_1", pW.hBondAcceptorCount == 1, "");
    emit("water_donor_and_acceptor",
        pW.hydrogenBondingClass == SaceHydrogenBondingClass::DonorAndAcceptor, "");

    emit("co2_carbonyl_2",
        saceBuiltinMolecularGraph(SUBSTANCE_CARBON_DIOXIDE, co2)
            && deriveFunctionalProfile(co2, pCO2) && pCO2.supported && pCO2.carbonylCount == 2, "");
    emit("co2_donor_0", pCO2.hBondDonorCount == 0, "");
    emit("co2_acceptor_2",
        pCO2.hBondAcceptorCount == 2
            && pCO2.hydrogenBondingClass == SaceHydrogenBondingClass::AcceptorOnly, "");

    emit("c2h6o_a_hydroxyl_1",
        saceSyntheticC2H6OGraphA(gA) && deriveFunctionalProfile(gA, pA)
            && pA.supported && pA.hydroxylCount == 1, "");
    emit("c2h6o_a_ether_0", pA.etherOxygenCount == 0, "");
    emit("c2h6o_a_donor_1", pA.hBondDonorCount == 1, "");
    emit("c2h6o_a_acceptor_1", pA.hBondAcceptorCount == 1, "");
    emit("c2h6o_a_cc_single_1", pA.carbonCarbonSingleCount == 1, "");
    emit("c2h6o_a_co_single_1", pA.carbonOxygenSingleCount == 1, "");
    emit("c2h6o_a_oh_single_1", pA.oxygenHydrogenSingleCount == 1, "");

    emit("c2h6o_b_hydroxyl_0",
        saceSyntheticC2H6OGraphB(gB) && deriveFunctionalProfile(gB, pB)
            && pB.supported && pB.hydroxylCount == 0, "");
    emit("c2h6o_b_ether_1", pB.etherOxygenCount == 1, "");
    emit("c2h6o_b_donor_0", pB.hBondDonorCount == 0, "");
    emit("c2h6o_b_acceptor_1", pB.hBondAcceptorCount == 1, "");
    emit("c2h6o_b_cc_single_0", pB.carbonCarbonSingleCount == 0, "");
    emit("c2h6o_b_co_single_2", pB.carbonOxygenSingleCount == 2, "");
    emit("c2h6o_b_oh_single_0", pB.oxygenHydrogenSingleCount == 0, "");

    ElementalComposition cA{}, cB{};
    emit("same_elemental_composition",
        elementalCompositionFromGraph(gA, cA) && elementalCompositionFromGraph(gB, cB)
            && elementalCompositionsEqual(cA, cB), "");
    emit("a_and_b_profiles_differ", !functionalProfilesEqual(pA, pB), "");

    SaceFunctionalProfile pPerm{}, pBond{}, pEnds{};
    emit("atom_permutation_preserves_profile",
        deriveFunctionalProfile(permuteWater(), pPerm)
            && functionalProfilesEqual(pW, pPerm), "");
    emit("bond_reorder_preserves_profile",
        deriveFunctionalProfile(reorderBonds(h2o), pBond)
            && functionalProfilesEqual(pW, pBond), "");
    emit("reversing_bond_endpoints_preserves_profile",
        deriveFunctionalProfile(reverseEnds(co2), pEnds)
            && functionalProfilesEqual(pCO2, pEnds), "");

    char sig[kSaceSignatureCap]{};
    bool wrote = writeChemicalSignature(saceBuiltinWaterIdentity(), sig, kSaceSignatureCap);
    emit("functional_profile_absent_from_canonical_signature",
        wrote && std::strstr(sig, "hydroxyl") == nullptr
            && std::strstr(sig, "ether") == nullptr
            && std::strstr(sig, "acceptor") == nullptr
            && std::strstr(sig, "functional") == nullptr, sig);

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
