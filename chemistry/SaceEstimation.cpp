#include "chemistry/SaceEstimation.h"

#include "chemistry/SaceCatalog.h"
#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

bool isH(AtomicNumber z) { return z == kAtomicHydrogen; }
bool isC(AtomicNumber z) { return z == kAtomicCarbon; }
bool isO(AtomicNumber z) { return z == kAtomicOxygen; }
bool isHCO(AtomicNumber z) { return isH(z) || isC(z) || isO(z); }

void addAtom(SaceMolecularGraph &g, AtomicNumber z, int8_t charge = 0) {
    g.atoms.push_back({z, charge});
}
void addBond(SaceMolecularGraph &g, uint16_t a, uint16_t b, SaceBondOrder order = SaceBondOrder::Single) {
    g.bonds.push_back({a, b, order});
}
void addHs(SaceMolecularGraph &g, uint16_t center, int nH) {
    for (int i = 0; i < nH; ++i) {
        uint16_t h = static_cast<uint16_t>(g.atoms.size());
        addAtom(g, kAtomicHydrogen);
        addBond(g, center, h);
    }
}

SaceMolecularGraph methaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addHs(g, 0, 4);
    return g;
}

SaceMolecularGraph isobutaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon); // 0 CH
    addAtom(g, kAtomicCarbon); // 1 CH3
    addAtom(g, kAtomicCarbon); // 2 CH3
    addAtom(g, kAtomicCarbon); // 3 CH3
    addBond(g, 0, 1);
    addBond(g, 0, 2);
    addBond(g, 0, 3);
    addHs(g, 0, 1);
    addHs(g, 1, 3);
    addHs(g, 2, 3);
    addHs(g, 3, 3);
    return g;
}

SaceMolecularGraph neopentaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon); // 0 C
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1);
    addBond(g, 0, 2);
    addBond(g, 0, 3);
    addBond(g, 0, 4);
    addHs(g, 1, 3);
    addHs(g, 2, 3);
    addHs(g, 3, 3);
    addHs(g, 4, 3);
    return g;
}

SaceMolecularGraph cyclopropaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1);
    addBond(g, 1, 2);
    addBond(g, 2, 0);
    addHs(g, 0, 2);
    addHs(g, 1, 2);
    addHs(g, 2, 2);
    return g;
}

SaceMolecularGraph etheneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1, SaceBondOrder::Double);
    addHs(g, 0, 2);
    addHs(g, 1, 2);
    return g;
}

SaceMolecularGraph ethyneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1, SaceBondOrder::Triple);
    addHs(g, 0, 1);
    addHs(g, 1, 1);
    return g;
}

SaceMolecularGraph chargedWaterLike() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicOxygen, 1);
    addHs(g, 0, 2);
    return g;
}

SaceMolecularGraph nitrogenGraph() {
    SaceMolecularGraph g{};
    addAtom(g, 7);
    addHs(g, 0, 3);
    return g;
}

SaceMolecularGraph permuteGraph(SaceMolecularGraph const &src) {
    SaceMolecularGraph g{};
    int n = atomCount(src);
    if (n <= 0) return g;
    std::vector<int> map(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        map[static_cast<size_t>(i)] = n - 1 - i;
    g.atoms.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        g.atoms[static_cast<size_t>(map[static_cast<size_t>(i)])] = src.atoms[static_cast<size_t>(i)];
    for (SaceBond const &b : src.bonds) {
        SaceBond nb = b;
        nb.atomA = static_cast<uint16_t>(map[b.atomA]);
        nb.atomB = static_cast<uint16_t>(map[b.atomB]);
        g.bonds.push_back(nb);
    }
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

bool nearK(double a, double b, double tol = 0.02) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

} // namespace

bool saceJobackGroupCountsEqual(SaceJobackGroupCounts const &a, SaceJobackGroupCounts const &b) {
    return a.carbonCH3 == b.carbonCH3 && a.carbonCH2 == b.carbonCH2
        && a.carbonCH == b.carbonCH && a.carbonC == b.carbonC
        && a.hydroxylAlcohol == b.hydroxylAlcohol && a.etherNonRing == b.etherNonRing;
}

SaceJobackFragmentationResult saceFragmentJobackGroups(
    SaceMolecularGraph const &graph, SaceJobackGroupCounts &counts)
{
    counts = {};
    if (validateMolecularGraph(graph, true) != SaceGraphValidation::Ok)
        return SaceJobackFragmentationResult::InvalidGraph;

    for (SaceAtom const &atom : graph.atoms) {
        if (!isHCO(atom.atomicNumber))
            return SaceJobackFragmentationResult::UnsupportedElement;
    }
    if (molecularGraphNetFormalCharge(graph) != 0)
        return SaceJobackFragmentationResult::Charged;

    SaceMolecularDescriptors desc{};
    if (!deriveMolecularDescriptors(graph, desc))
        return SaceJobackFragmentationResult::InvalidGraph;
    if (desc.cycleRank != 0)
        return SaceJobackFragmentationResult::Cyclic;

    for (SaceBond const &b : graph.bonds) {
        if (b.order != SaceBondOrder::Single)
            return SaceJobackFragmentationResult::Unsaturated;
    }

    int n = atomCount(graph);
    for (int i = 0; i < n; ++i) {
        if (!isH(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 1)
            return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
        uint16_t idx = static_cast<uint16_t>(i);
        bool okH = false;
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            if (b.order != SaceBondOrder::Single)
                return SaceJobackFragmentationResult::Unsaturated;
            AtomicNumber z = graph.atoms[static_cast<size_t>(other)].atomicNumber;
            if (!isC(z) && !isO(z))
                return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
            okH = true;
        }
        if (!okH)
            return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
    }

    std::vector<char> covered(static_cast<size_t>(n), 0);

    auto cover = [&](int idx) {
        if (idx < 0 || idx >= n) return false;
        if (covered[static_cast<size_t>(idx)]) return false;
        covered[static_cast<size_t>(idx)] = 1;
        return true;
    };

    auto neighborZs = [&](int atom, int &hN, int &cN, int &oN, int &otherN) {
        hN = cN = oN = otherN = 0;
        uint16_t idx = static_cast<uint16_t>(atom);
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            AtomicNumber z = graph.atoms[static_cast<size_t>(other)].atomicNumber;
            if (isH(z)) ++hN;
            else if (isC(z)) ++cN;
            else if (isO(z)) ++oN;
            else ++otherN;
        }
    };

    for (int i = 0; i < n; ++i) {
        if (!isC(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        int hN = 0, cN = 0, oN = 0, otherN = 0;
        neighborZs(i, hN, cN, oN, otherN);
        if (otherN != 0 || hN + cN + oN != 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (hN == 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (hN == 3) ++counts.carbonCH3;
        else if (hN == 2) ++counts.carbonCH2;
        else if (hN == 1) ++counts.carbonCH;
        else if (hN == 0) ++counts.carbonC;
        else
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (!cover(i))
            return SaceJobackFragmentationResult::IncompleteCoverage;
        uint16_t idx = static_cast<uint16_t>(i);
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            if (!isH(graph.atoms[static_cast<size_t>(other)].atomicNumber))
                continue;
            if (!cover(other))
                return SaceJobackFragmentationResult::IncompleteCoverage;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!isO(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 2)
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        int hN = 0, cN = 0, oN = 0, otherN = 0;
        neighborZs(i, hN, cN, oN, otherN);
        if (otherN != 0 || oN != 0)
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        if (hN == 1 && cN == 1) {
            ++counts.hydroxylAlcohol;
            if (!cover(i))
                return SaceJobackFragmentationResult::IncompleteCoverage;
            uint16_t idx = static_cast<uint16_t>(i);
            for (SaceBond const &b : graph.bonds) {
                int other = -1;
                if (b.atomA == idx) other = b.atomB;
                else if (b.atomB == idx) other = b.atomA;
                if (other < 0) continue;
                if (!isH(graph.atoms[static_cast<size_t>(other)].atomicNumber))
                    continue;
                if (!cover(other))
                    return SaceJobackFragmentationResult::IncompleteCoverage;
            }
        } else if (hN == 0 && cN == 2) {
            ++counts.etherNonRing;
            if (!cover(i))
                return SaceJobackFragmentationResult::IncompleteCoverage;
        } else {
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!covered[static_cast<size_t>(i)])
            return SaceJobackFragmentationResult::IncompleteCoverage;
    }
    return SaceJobackFragmentationResult::Ok;
}

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out, SaceJobackGroupCounts &counts)
{
    out = saceUnknownScalarProperty();
    counts = {};
    SaceJobackFragmentationResult r = saceFragmentJobackGroups(graph, counts);
    if (r != SaceJobackFragmentationResult::Ok) {
        counts = {};
        return false;
    }
    double tb = kJobackTbInterceptK
        + kJobackTbCH3K * counts.carbonCH3
        + kJobackTbCH2K * counts.carbonCH2
        + kJobackTbCHK * counts.carbonCH
        + kJobackTbCK * counts.carbonC
        + kJobackTbAlcoholOHK * counts.hydroxylAlcohol
        + kJobackTbEtherOK * counts.etherNonRing;
    if (!(tb > 0.0) || !std::isfinite(tb)) {
        counts = {};
        return false;
    }
    out.value = static_cast<float>(tb);
    out.known = true;
    out.source = SacePropertySource::StructuralEstimate;
    out.confidence = SaceConfidence::Low;
    return true;
}

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out)
{
    SaceJobackGroupCounts counts{};
    return saceEstimateJobackNormalBoilingPoint(graph, out, counts);
}

void runSaceEstimationDiagnostics() {
    std::ofstream out(miscFile("sace_estimation_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    auto unsupported = [&](char const *name, SaceMolecularGraph const &g,
        SaceJobackFragmentationResult expect)
    {
        SaceJobackGroupCounts c{};
        SaceScalarProperty p = {};
        p.known = true;
        p.value = 99.0f;
        SaceJobackFragmentationResult r = saceFragmentJobackGroups(g, c);
        bool est = saceEstimateJobackNormalBoilingPoint(g, p);
        emit(name, r == expect && !est && !p.known
            && p.source == SacePropertySource::Unknown
            && p.confidence == SaceConfidence::Unknown, saceJobackFragmentationResultKey(r));
    };

    SaceMolecularGraph empty{};
    SaceScalarProperty junk{};
    junk.known = true;
    junk.value = 400.0f;
    junk.source = SacePropertySource::StructuralEstimate;
    emit("invalid_graph_unknown",
        !saceEstimateJobackNormalBoilingPoint(empty, junk) && !junk.known, "");

    SaceMolecularGraph water{}, h2{}, o2{}, co2{};
    saceBuiltinMolecularGraph(SUBSTANCE_WATER, water);
    saceBuiltinMolecularGraph(SUBSTANCE_HYDROGEN, h2);
    saceBuiltinMolecularGraph(SUBSTANCE_OXYGEN, o2);
    saceBuiltinMolecularGraph(SUBSTANCE_CARBON_DIOXIDE, co2);
    unsupported("water_unsupported", water, SaceJobackFragmentationResult::UnsupportedOxygenEnvironment);
    unsupported("h2_unsupported", h2, SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment);
    unsupported("o2_unsupported", o2, SaceJobackFragmentationResult::Unsaturated);
    unsupported("co2_unsupported", co2, SaceJobackFragmentationResult::Unsaturated);
    unsupported("methane_unsupported", methaneGraph(), SaceJobackFragmentationResult::UnsupportedCarbonEnvironment);
    unsupported("charged_unsupported", chargedWaterLike(), SaceJobackFragmentationResult::Charged);
    unsupported("ring_unsupported", cyclopropaneGraph(), SaceJobackFragmentationResult::Cyclic);
    unsupported("nitrogen_unsupported", nitrogenGraph(), SaceJobackFragmentationResult::UnsupportedElement);
    unsupported("cc_double_unsupported", etheneGraph(), SaceJobackFragmentationResult::Unsaturated);
    unsupported("cc_triple_unsupported", ethyneGraph(), SaceJobackFragmentationResult::Unsaturated);

    SaceMolecularGraph gA{}, gB{};
    saceSyntheticC2H6OGraphA(gA);
    saceSyntheticC2H6OGraphB(gB);
    SaceJobackGroupCounts cA{}, cB{};
    SaceScalarProperty tbA{}, tbB{};
    emit("c2h6o_a_fragmentation_ok",
        saceEstimateJobackNormalBoilingPoint(gA, tbA, cA)
            && saceFragmentJobackGroups(gA, cA) == SaceJobackFragmentationResult::Ok, "");
    emit("a_ch3_count_1", cA.carbonCH3 == 1, std::to_string(cA.carbonCH3));
    emit("a_ch2_count_1", cA.carbonCH2 == 1, std::to_string(cA.carbonCH2));
    emit("a_oh_count_1", cA.hydroxylAlcohol == 1, std::to_string(cA.hydroxylAlcohol));
    emit("a_ether_count_0", cA.etherNonRing == 0, std::to_string(cA.etherNonRing));
    emit("a_tb_337_54", tbA.known && nearK(tbA.value, 337.54), std::to_string(tbA.value));
    emit("a_source_structural_estimate",
        tbA.source == SacePropertySource::StructuralEstimate, sacePropertySourceKey(tbA.source));
    emit("a_confidence_low",
        tbA.confidence == SaceConfidence::Low, saceConfidenceKey(tbA.confidence));

    emit("c2h6o_b_fragmentation_ok",
        saceEstimateJobackNormalBoilingPoint(gB, tbB, cB), "");
    emit("b_ch3_count_2", cB.carbonCH3 == 2, std::to_string(cB.carbonCH3));
    emit("b_oh_count_0", cB.hydroxylAlcohol == 0, std::to_string(cB.hydroxylAlcohol));
    emit("b_ether_count_1", cB.etherNonRing == 1, std::to_string(cB.etherNonRing));
    emit("b_tb_267_78", tbB.known && nearK(tbB.value, 267.78), std::to_string(tbB.value));
    emit("b_source_structural_estimate",
        tbB.source == SacePropertySource::StructuralEstimate, "");
    emit("b_confidence_low", tbB.confidence == SaceConfidence::Low, "");

    SaceScalarProperty mmA{}, mmB{};
    ElementalComposition elA{}, elB{};
    elementalCompositionFromGraph(gA, elA);
    elementalCompositionFromGraph(gB, elB);
    saceDeriveMolarMass(elA, mmA);
    saceDeriveMolarMass(elB, mmB);
    emit("a_b_molar_mass_equal",
        mmA.known && mmB.known && nearK(mmA.value, mmB.value, 1e-4), "");
    emit("a_b_boiling_points_differ",
        tbA.known && tbB.known && !nearK(tbA.value, tbB.value, 0.5), "");
    emit("a_tb_greater_than_b", tbA.known && tbB.known && tbA.value > tbB.value, "");
    emit("tb_difference_69_76",
        tbA.known && tbB.known && nearK(tbA.value - tbB.value, 69.76),
        std::to_string(tbA.value - tbB.value));

    emit("ethanol_topology_vs_literature", true,
        "joback=" + std::to_string(tbA.value) + " ref=351.44 err="
            + std::to_string(tbA.value - 351.44f));
    emit("dme_topology_vs_literature", true,
        "joback=" + std::to_string(tbB.value) + " ref=248.3 err="
            + std::to_string(tbB.value - 248.3f));

    SaceJobackGroupCounts isoC{};
    SaceScalarProperty isoTb{};
    bool isoOk = saceEstimateJobackNormalBoilingPoint(isobutaneGraph(), isoTb, isoC)
        && isoC.carbonCH == 1 && isoC.carbonCH3 == 3;
    emit("isobutane_exercises_ch", isoOk, std::to_string(isoC.carbonCH));
    SaceJobackGroupCounts neoC{};
    SaceScalarProperty neoTb{};
    bool neoOk = saceEstimateJobackNormalBoilingPoint(neopentaneGraph(), neoTb, neoC)
        && neoC.carbonC == 1 && neoC.carbonCH3 == 4;
    emit("neopentane_exercises_quaternary_c", neoOk, std::to_string(neoC.carbonC));

    SaceJobackGroupCounts permC{}, bondC{}, endC{};
    SaceScalarProperty permTb{}, bondTb{}, endTb{};
    emit("atom_permutation_preserves_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(permuteGraph(gA), permTb, permC)
            && saceJobackGroupCountsEqual(cA, permC) && nearK(permTb.value, tbA.value), "");
    emit("bond_reorder_preserves_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(reorderBonds(gA), bondTb, bondC)
            && saceJobackGroupCountsEqual(cA, bondC) && nearK(bondTb.value, tbA.value), "");
    emit("reversed_endpoints_preserve_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(reverseEnds(gA), endTb, endC)
            && saceJobackGroupCountsEqual(cA, endC) && nearK(endTb.value, tbA.value), "");

    SaceScalarProperty resetP{};
    resetP.known = true;
    resetP.value = 12.0f;
    resetP.source = SacePropertySource::StructuralEstimate;
    resetP.confidence = SaceConfidence::Low;
    emit("failed_estimation_resets_unknown",
        !saceEstimateJobackNormalBoilingPoint(empty, resetP) && !resetP.known
            && resetP.source == SacePropertySource::Unknown, "");

    char sig[kSaceSignatureCap]{};
    writeChemicalSignature(saceBuiltinWaterIdentity(), sig, kSaceSignatureCap);
    emit("joback_not_in_canonical_signature",
        std::strstr(sig, "337") == nullptr && std::strstr(sig, "Joback") == nullptr
            && std::strstr(sig, "boiling") == nullptr, sig);

    SaceScalarProperty refHigh{};
    refHigh.known = true;
    refHigh.value = 351.44f;
    refHigh.source = SacePropertySource::Reference;
    refHigh.confidence = SaceConfidence::High;
    SaceScalarProperty jobackLow = tbA;
    emit("higher_quality_not_overwritten_by_joback",
        !saceAssignScalarProperty(refHigh, jobackLow)
            && refHigh.source == SacePropertySource::Reference
            && nearK(refHigh.value, 351.44), "");

    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();
    ElementCount c2h6o[] = {
        {kAtomicCarbon, 2}, {kAtomicHydrogen, 6}, {kAtomicOxygen, 1}
    };
    ChemicalIdentity idA = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-a", c2h6o, 3);
    ChemicalIdentity idB = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-b", c2h6o, 3);
    SaceSubstanceRef refA = cat.resolve(idA, true);
    SaceSubstanceRef refB = cat.resolve(idB, true);
    cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    cat.attachMolecularGraph(refB.generatedId, "synthetic-structure-b", gB);
    SaceGeneratedRecord const *recA = cat.record(refA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(refB.generatedId);
    emit("catalog_a_joback_tb",
        recA && recA->hasMolecularGraph && recA->hasMolecularDescriptors
            && recA->hasFunctionalProfile
            && recA->properties.normalBoilingPointK.known
            && recA->properties.normalBoilingPointK.source == SacePropertySource::StructuralEstimate
            && recA->properties.normalBoilingPointK.confidence == SaceConfidence::Low
            && nearK(recA->properties.normalBoilingPointK.value, 337.54),
        recA ? std::to_string(recA->properties.normalBoilingPointK.value) : "");
    emit("catalog_b_joback_tb",
        recB && recB->properties.normalBoilingPointK.known
            && nearK(recB->properties.normalBoilingPointK.value, 267.78),
        recB ? std::to_string(recB->properties.normalBoilingPointK.value) : "");

    ElementCount c2h4[] = {{kAtomicCarbon, 2}, {kAtomicHydrogen, 4}};
    ChemicalIdentity idE = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H4", "synthetic-ethene", c2h4, 2);
    SaceSubstanceRef refE = cat.resolve(idE, true);
    SaceMolecularGraph ethene = etheneGraph();
    bool attached = cat.attachMolecularGraph(refE.generatedId, "synthetic-ethene", ethene);
    SaceGeneratedRecord const *recE = cat.record(refE.generatedId);
    emit("unsupported_joback_still_attaches_graph",
        attached && recE && recE->hasMolecularGraph && recE->hasFunctionalProfile
            && !recE->properties.normalBoilingPointK.known, "");

    cat.clear();
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
