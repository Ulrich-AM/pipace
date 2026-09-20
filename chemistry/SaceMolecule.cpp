#include "chemistry/SaceMolecule.h"

#include "chemistry/SaceCatalog.h"
#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceProperties.h"
#include "fluid/DiagOutput.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <queue>
#include <string>

namespace {

void addBond(SaceMolecularGraph &g, uint16_t a, uint16_t b, SaceBondOrder order) {
    g.bonds.push_back({a, b, order});
}

void addAtom(SaceMolecularGraph &g, AtomicNumber z) {
    g.atoms.push_back({z, 0});
}

void addHydrogens(SaceMolecularGraph &g, uint16_t center, int nH) {
    for (int i = 0; i < nH; ++i) {
        uint16_t h = static_cast<uint16_t>(g.atoms.size());
        addAtom(g, kAtomicHydrogen);
        addBond(g, center, h, SaceBondOrder::Single);
    }
}

bool bondUndirectedEqual(SaceBond const &a, SaceBond const &b) {
    uint16_t a0 = std::min(a.atomA, a.atomB);
    uint16_t a1 = std::max(a.atomA, a.atomB);
    uint16_t b0 = std::min(b.atomA, b.atomB);
    uint16_t b1 = std::max(b.atomA, b.atomB);
    return a0 == b0 && a1 == b1 && a.order == b.order;
}

int appendDebug(char *out, int cap, int used, char const *s) {
    if (!out || cap <= 0 || !s) return used;
    while (*s) {
        if (used >= cap - 1) return -1;
        out[used++] = *s++;
    }
    out[used] = '\0';
    return used;
}

int appendDebugU32(char *out, int cap, int used, unsigned v) {
    char buf[16];
    int n = std::snprintf(buf, sizeof(buf), "%u", v);
    if (n <= 0) return -1;
    return appendDebug(out, cap, used, buf);
}

char bondOrderChar(SaceBondOrder order) {
    switch (order) {
        case SaceBondOrder::Single: return '1';
        case SaceBondOrder::Double: return '2';
        case SaceBondOrder::Triple: return '3';
        case SaceBondOrder::Aromatic: return 'a';
    }
    return '?';
}

int bondOrderWeight(SaceBondOrder order) {
    switch (order) {
        case SaceBondOrder::Single: return 1;
        case SaceBondOrder::Double: return 2;
        case SaceBondOrder::Triple: return 3;
        case SaceBondOrder::Aromatic: return 2;
    }
    return 0;
}

bool graphConnected(SaceMolecularGraph const &graph) {
    int n = static_cast<int>(graph.atoms.size());
    if (n <= 1) return n == 1;
    std::vector<char> seen(static_cast<size_t>(n), 0);
    std::queue<int> q;
    seen[0] = 1;
    q.push(0);
    int reached = 1;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (SaceBond const &b : graph.bonds) {
            int v = -1;
            if (b.atomA == u) v = b.atomB;
            else if (b.atomB == u) v = b.atomA;
            if (v < 0 || v >= n || seen[static_cast<size_t>(v)]) continue;
            seen[static_cast<size_t>(v)] = 1;
            q.push(v);
            ++reached;
        }
    }
    return reached == n;
}

} // namespace

SaceGraphValidation validateMolecularGraph(SaceMolecularGraph const &graph, bool requireConnected) {
    int nAtoms = static_cast<int>(graph.atoms.size());
    int nBonds = static_cast<int>(graph.bonds.size());
    if (nAtoms <= 0) return SaceGraphValidation::Empty;
    if (nAtoms > kMaxSaceMoleculeAtoms) return SaceGraphValidation::TooManyAtoms;
    if (nBonds > kMaxSaceMoleculeBonds) return SaceGraphValidation::TooManyBonds;
    for (SaceAtom const &atom : graph.atoms) {
        if (!validAtomicNumber(atom.atomicNumber))
            return SaceGraphValidation::InvalidAtom;
    }
    for (int i = 0; i < nBonds; ++i) {
        SaceBond const &b = graph.bonds[static_cast<size_t>(i)];
        if (!validSaceBondOrder(b.order))
            return SaceGraphValidation::InvalidBond;
        if (b.atomA >= static_cast<uint16_t>(nAtoms) || b.atomB >= static_cast<uint16_t>(nAtoms))
            return SaceGraphValidation::InvalidBond;
        if (b.atomA == b.atomB)
            return SaceGraphValidation::SelfBond;
        for (int j = 0; j < i; ++j) {
            SaceBond const &prev = graph.bonds[static_cast<size_t>(j)];
            uint16_t a0 = std::min(b.atomA, b.atomB);
            uint16_t a1 = std::max(b.atomA, b.atomB);
            uint16_t p0 = std::min(prev.atomA, prev.atomB);
            uint16_t p1 = std::max(prev.atomA, prev.atomB);
            if (a0 == p0 && a1 == p1)
                return SaceGraphValidation::DuplicateBond;
        }
    }
    if (requireConnected && !graphConnected(graph))
        return SaceGraphValidation::Disconnected;
    return SaceGraphValidation::Ok;
}

bool elementalCompositionFromGraph(SaceMolecularGraph const &graph, ElementalComposition &out) {
    out = {};
    if (validateMolecularGraph(graph, false) != SaceGraphValidation::Ok)
        return false;
    ElementalComposition raw{};
    for (SaceAtom const &atom : graph.atoms) {
        bool found = false;
        for (int i = 0; i < raw.count && i < kMaxElementalSpecies; ++i) {
            if (raw.entries[i].atomicNumber == atom.atomicNumber) {
                uint32_t sum = static_cast<uint32_t>(raw.entries[i].count) + 1u;
                if (sum > 0xFFFFu) return false;
                raw.entries[i].count = static_cast<uint16_t>(sum);
                found = true;
                break;
            }
        }
        if (!found) {
            if (raw.count >= kMaxElementalSpecies) return false;
            raw.entries[raw.count] = {atom.atomicNumber, 1};
            ++raw.count;
        }
    }
    if (normalizeElementalComposition(raw) != CompositionNormalizeResult::Ok)
        return false;
    out = raw;
    return true;
}

bool graphMatchesChemicalIdentity(SaceMolecularGraph const &graph, ChemicalIdentity const &id) {
    if (!id.exact) return false;
    ElementalComposition fromGraph{};
    if (!elementalCompositionFromGraph(graph, fromGraph))
        return false;
    ChemicalIdentity norm = id;
    if (norm.elemental.count > kMaxElementalSpecies)
        return false;
    if (normalizeElementalComposition(norm.elemental) != CompositionNormalizeResult::Ok)
        return false;
    return elementalCompositionsEqual(fromGraph, norm.elemental);
}

bool molecularGraphsStoredEqual(SaceMolecularGraph const &a, SaceMolecularGraph const &b) {
    if (a.atoms.size() != b.atoms.size() || a.bonds.size() != b.bonds.size())
        return false;
    for (size_t i = 0; i < a.atoms.size(); ++i) {
        if (a.atoms[i].atomicNumber != b.atoms[i].atomicNumber) return false;
        if (a.atoms[i].formalCharge != b.atoms[i].formalCharge) return false;
    }
    std::vector<char> used(b.bonds.size(), 0);
    for (SaceBond const &bond : a.bonds) {
        bool found = false;
        for (size_t i = 0; i < b.bonds.size(); ++i) {
            if (used[i]) continue;
            if (!bondUndirectedEqual(bond, b.bonds[i])) continue;
            used[i] = 1;
            found = true;
            break;
        }
        if (!found) return false;
    }
    return true;
}

int32_t molecularGraphNetFormalCharge(SaceMolecularGraph const &graph) {
    int32_t sum = 0;
    for (SaceAtom const &atom : graph.atoms)
        sum += atom.formalCharge;
    return sum;
}

int atomCount(SaceMolecularGraph const &graph) {
    return static_cast<int>(graph.atoms.size());
}

int bondCount(SaceMolecularGraph const &graph) {
    return static_cast<int>(graph.bonds.size());
}

int atomDegree(SaceMolecularGraph const &graph, int atomIndex) {
    if (atomIndex < 0 || atomIndex >= atomCount(graph)) return 0;
    int deg = 0;
    uint16_t idx = static_cast<uint16_t>(atomIndex);
    for (SaceBond const &b : graph.bonds)
        if (b.atomA == idx || b.atomB == idx) ++deg;
    return deg;
}

int countElement(SaceMolecularGraph const &graph, AtomicNumber z) {
    int n = 0;
    for (SaceAtom const &atom : graph.atoms)
        if (atom.atomicNumber == z) ++n;
    return n;
}

bool hasBond(SaceMolecularGraph const &graph, int atomA, int atomB) {
    if (atomA < 0 || atomB < 0) return false;
    uint16_t a = static_cast<uint16_t>(atomA);
    uint16_t b = static_cast<uint16_t>(atomB);
    for (SaceBond const &bond : graph.bonds) {
        if ((bond.atomA == a && bond.atomB == b) || (bond.atomA == b && bond.atomB == a))
            return true;
    }
    return false;
}

bool writeMolecularGraphDebug(SaceMolecularGraph const &graph, char *out, int cap) {
    if (!out || cap <= 0) return false;
    out[0] = '\0';
    int used = appendDebug(out, cap, 0, "debug|");
    if (used < 0) { out[0] = '\0'; return false; }
    for (size_t i = 0; i < graph.atoms.size(); ++i) {
        if (i > 0) { used = appendDebug(out, cap, used, ","); if (used < 0) { out[0] = '\0'; return false; } }
        used = appendDebugU32(out, cap, used, graph.atoms[i].atomicNumber);
        if (used < 0) { out[0] = '\0'; return false; }
        used = appendDebug(out, cap, used, ":");
        if (used < 0) { out[0] = '\0'; return false; }
        char chbuf[8];
        std::snprintf(chbuf, sizeof(chbuf), "%d", static_cast<int>(graph.atoms[i].formalCharge));
        used = appendDebug(out, cap, used, chbuf);
        if (used < 0) { out[0] = '\0'; return false; }
    }
    used = appendDebug(out, cap, used, "|");
    if (used < 0) { out[0] = '\0'; return false; }
    for (size_t i = 0; i < graph.bonds.size(); ++i) {
        if (i > 0) { used = appendDebug(out, cap, used, ","); if (used < 0) { out[0] = '\0'; return false; } }
        SaceBond const &b = graph.bonds[i];
        used = appendDebugU32(out, cap, used, b.atomA);
        if (used < 0) { out[0] = '\0'; return false; }
        used = appendDebug(out, cap, used, "-");
        if (used < 0) { out[0] = '\0'; return false; }
        used = appendDebugU32(out, cap, used, b.atomB);
        if (used < 0) { out[0] = '\0'; return false; }
        char oc[2] = {bondOrderChar(b.order), 0};
        used = appendDebug(out, cap, used, ":");
        if (used < 0) { out[0] = '\0'; return false; }
        used = appendDebug(out, cap, used, oc);
        if (used < 0) { out[0] = '\0'; return false; }
    }
    return true;
}

bool saceCommonValenceLooksTypical(SaceMolecularGraph const &graph) {
    if (validateMolecularGraph(graph, false) != SaceGraphValidation::Ok)
        return false;
    int n = atomCount(graph);
    for (int i = 0; i < n; ++i) {
        int sum = 0;
        uint16_t idx = static_cast<uint16_t>(i);
        for (SaceBond const &b : graph.bonds)
            if (b.atomA == idx || b.atomB == idx)
                sum += bondOrderWeight(b.order);
        AtomicNumber z = graph.atoms[static_cast<size_t>(i)].atomicNumber;
        if (z == kAtomicHydrogen && sum > 1) return false;
        if (z == kAtomicOxygen && sum > 2) return false;
        if (z == kAtomicCarbon && sum > 4) return false;
    }
    return true;
}

bool saceBuiltinMolecularGraph(SubstanceId id, SaceMolecularGraph &out) {
    out = {};
    if (id == SUBSTANCE_HYDROGEN) {
        addAtom(out, kAtomicHydrogen);
        addAtom(out, kAtomicHydrogen);
        addBond(out, 0, 1, SaceBondOrder::Single);
        return true;
    }
    if (id == SUBSTANCE_OXYGEN) {
        addAtom(out, kAtomicOxygen);
        addAtom(out, kAtomicOxygen);
        addBond(out, 0, 1, SaceBondOrder::Double);
        return true;
    }
    if (id == SUBSTANCE_WATER) {
        addAtom(out, kAtomicOxygen);
        addAtom(out, kAtomicHydrogen);
        addAtom(out, kAtomicHydrogen);
        addBond(out, 0, 1, SaceBondOrder::Single);
        addBond(out, 0, 2, SaceBondOrder::Single);
        return true;
    }
    if (id == SUBSTANCE_CARBON_DIOXIDE) {
        addAtom(out, kAtomicCarbon);
        addAtom(out, kAtomicOxygen);
        addAtom(out, kAtomicOxygen);
        addBond(out, 0, 1, SaceBondOrder::Double);
        addBond(out, 0, 2, SaceBondOrder::Double);
        return true;
    }
    return false;
}

char const *saceBuiltinMolecularStructureKey(SubstanceId id) {
    SaceMolecularGraph unused{};
    if (!saceBuiltinMolecularGraph(id, unused))
        return nullptr;
    char const *key = substanceDef(id).chemicalIdentity.structureKey;
    return saceStructureKeyOrEmpty(key)[0] ? key : nullptr;
}

bool saceSyntheticC2H6OGraphA(SaceMolecularGraph &out) {
    out = {};
    addAtom(out, kAtomicCarbon); // 0
    addAtom(out, kAtomicCarbon); // 1
    addAtom(out, kAtomicOxygen); // 2
    addBond(out, 0, 1, SaceBondOrder::Single);
    addBond(out, 1, 2, SaceBondOrder::Single);
    addHydrogens(out, 0, 3);
    addHydrogens(out, 1, 2);
    addHydrogens(out, 2, 1);
    return true;
}

bool saceSyntheticC2H6OGraphB(SaceMolecularGraph &out) {
    out = {};
    addAtom(out, kAtomicCarbon); // 0
    addAtom(out, kAtomicOxygen); // 1
    addAtom(out, kAtomicCarbon); // 2
    addBond(out, 0, 1, SaceBondOrder::Single);
    addBond(out, 1, 2, SaceBondOrder::Single);
    addHydrogens(out, 0, 3);
    addHydrogens(out, 2, 3);
    return true;
}

void runSaceMoleculeDiagnostics() {
    std::ofstream out(miscFile("sace_molecule_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    SaceMolecularGraph empty{};
    emit("empty_graph_invalid",
        validateMolecularGraph(empty, true) == SaceGraphValidation::Empty, "");

    SaceMolecularGraph badZ{};
    badZ.atoms.push_back({0, 0});
    emit("invalid_atomic_number_rejected",
        validateMolecularGraph(badZ, false) == SaceGraphValidation::InvalidAtom, "");

    SaceMolecularGraph selfBond{};
    addAtom(selfBond, kAtomicOxygen);
    addBond(selfBond, 0, 0, SaceBondOrder::Single);
    emit("self_bond_rejected",
        validateMolecularGraph(selfBond, false) == SaceGraphValidation::SelfBond, "");

    SaceMolecularGraph oor{};
    addAtom(oor, kAtomicHydrogen);
    addBond(oor, 0, 1, SaceBondOrder::Single);
    emit("out_of_range_bond_rejected",
        validateMolecularGraph(oor, false) == SaceGraphValidation::InvalidBond, "");

    SaceMolecularGraph dup{};
    addAtom(dup, kAtomicHydrogen);
    addAtom(dup, kAtomicHydrogen);
    addBond(dup, 0, 1, SaceBondOrder::Single);
    addBond(dup, 1, 0, SaceBondOrder::Single);
    emit("duplicate_reversed_bond_rejected",
        validateMolecularGraph(dup, false) == SaceGraphValidation::DuplicateBond, "");

    SaceMolecularGraph disconnected{};
    addAtom(disconnected, kAtomicHydrogen);
    addAtom(disconnected, kAtomicHydrogen);
    emit("disconnected_small_molecule_rejected",
        validateMolecularGraph(disconnected, true) == SaceGraphValidation::Disconnected, "");

    SaceMolecularGraph h2{}, o2{}, h2o{}, co2{};
    emit("h2_graph_valid",
        saceBuiltinMolecularGraph(SUBSTANCE_HYDROGEN, h2)
            && validateMolecularGraph(h2, true) == SaceGraphValidation::Ok, "");
    emit("o2_graph_valid",
        saceBuiltinMolecularGraph(SUBSTANCE_OXYGEN, o2)
            && validateMolecularGraph(o2, true) == SaceGraphValidation::Ok, "");
    emit("h2o_graph_valid",
        saceBuiltinMolecularGraph(SUBSTANCE_WATER, h2o)
            && validateMolecularGraph(h2o, true) == SaceGraphValidation::Ok, "");
    emit("co2_graph_valid",
        saceBuiltinMolecularGraph(SUBSTANCE_CARBON_DIOXIDE, co2)
            && validateMolecularGraph(co2, true) == SaceGraphValidation::Ok, "");

    emit("h2_graph_matches_identity",
        graphMatchesChemicalIdentity(h2, substanceDef(SUBSTANCE_HYDROGEN).chemicalIdentity), "");
    emit("o2_graph_matches_identity",
        graphMatchesChemicalIdentity(o2, substanceDef(SUBSTANCE_OXYGEN).chemicalIdentity), "");
    emit("h2o_graph_matches_identity",
        graphMatchesChemicalIdentity(h2o, substanceDef(SUBSTANCE_WATER).chemicalIdentity), "");
    emit("co2_graph_matches_identity",
        graphMatchesChemicalIdentity(co2, substanceDef(SUBSTANCE_CARBON_DIOXIDE).chemicalIdentity), "");
    emit("water_net_formal_charge_zero",
        molecularGraphNetFormalCharge(h2o) == 0, std::to_string(molecularGraphNetFormalCharge(h2o)));

    emit("builtin_structure_keys",
        saceBuiltinMolecularStructureKey(SUBSTANCE_WATER)
            && std::strcmp(saceBuiltinMolecularStructureKey(SUBSTANCE_WATER), "water") == 0
            && std::strcmp(saceBuiltinMolecularStructureKey(SUBSTANCE_HYDROGEN), "molecular-hydrogen") == 0
            && std::strcmp(saceBuiltinMolecularStructureKey(SUBSTANCE_OXYGEN), "molecular-oxygen") == 0
            && std::strcmp(saceBuiltinMolecularStructureKey(SUBSTANCE_CARBON_DIOXIDE), "carbon-dioxide") == 0
            && saceBuiltinMolecularStructureKey(SUBSTANCE_CARBON) == nullptr, "");

    SaceMolecularGraph carbonGraph{};
    emit("carbon_has_no_small_molecule_graph",
        !saceBuiltinMolecularGraph(SUBSTANCE_CARBON, carbonGraph), "");

    SaceMolecularGraph gA{}, gB{};
    emit("c2h6o_graph_a_valid",
        saceSyntheticC2H6OGraphA(gA)
            && validateMolecularGraph(gA, true) == SaceGraphValidation::Ok
            && saceCommonValenceLooksTypical(gA), "");
    emit("c2h6o_graph_b_valid",
        saceSyntheticC2H6OGraphB(gB)
            && validateMolecularGraph(gB, true) == SaceGraphValidation::Ok
            && saceCommonValenceLooksTypical(gB), "");

    ElementalComposition compA{}, compB{};
    bool derivedA = elementalCompositionFromGraph(gA, compA);
    bool derivedB = elementalCompositionFromGraph(gB, compB);
    emit("both_derive_c2h6o",
        derivedA && derivedB
            && elementalCountOf(compA, kAtomicCarbon) == 2
            && elementalCountOf(compA, kAtomicHydrogen) == 6
            && elementalCountOf(compA, kAtomicOxygen) == 1
            && elementalCompositionsEqual(compA, compB), "");
    emit("topologies_differ",
        hasBond(gA, 0, 1) && !hasBond(gB, 0, 2)
            && hasBond(gB, 0, 1) && hasBond(gB, 1, 2)
            && !molecularGraphsStoredEqual(gA, gB), "");

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
    SaceGeneratedRecord const *recA = cat.record(refA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(refB.generatedId);
    std::string sigA = recA ? recA->canonicalSignature : "";
    std::string sigB = recB ? recB->canonicalSignature : "";
    uint32_t ordA = recA ? recA->displayOrdinal : 0;
    uint32_t ordB = recB ? recB->displayOrdinal : 0;
    SaceRecordId ridA = recA ? recA->recordId : kSaceRecordNone;
    SaceRecordId ridB = recB ? recB->recordId : kSaceRecordNone;
    float massA = recA ? recA->properties.molarMassGPerMol.value : 0.0f;

    char const *keyA = "synthetic-structure-a";
    char const *keyB = "synthetic-structure-b";
    bool rejectBOnFreshA = cat.attachMolecularGraph(refA.generatedId, keyB, gB);
    recA = cat.record(refA.generatedId);
    emit("fresh_a_rejects_graph_b_key_b",
        !rejectBOnFreshA && recA && !recA->hasMolecularGraph && !recA->hasMolecularDescriptors, "");

    bool attachA = recA && cat.attachMolecularGraph(refA.generatedId, keyA, gA);
    recA = cat.record(refA.generatedId);
    emit("fresh_a_accepts_graph_a_key_a",
        attachA && recA && recA->hasMolecularGraph
            && recA->molecularGraphStructureKey == recA->structureKey
            && recA->structureKey == keyA
            && recA->hasMolecularDescriptors
            && molecularGraphsStoredEqual(recA->molecularGraph, gA), "");

    recB = cat.record(refB.generatedId);
    bool rejectAOnFreshB = recB && !recB->hasMolecularGraph
        && cat.attachMolecularGraph(refB.generatedId, keyA, gA);
    recB = cat.record(refB.generatedId);
    emit("fresh_b_rejects_graph_a_key_a",
        !rejectAOnFreshB && recB && !recB->hasMolecularGraph && !recB->hasMolecularDescriptors, "");

    bool attachB = recB && cat.attachMolecularGraph(refB.generatedId, keyB, gB);
    recB = cat.record(refB.generatedId);
    emit("fresh_b_accepts_graph_b_key_b",
        attachB && recB && recB->hasMolecularGraph
            && recB->molecularGraphStructureKey == recB->structureKey
            && recB->structureKey == keyB, "");

    bool attachBOnA = cat.attachMolecularGraph(refA.generatedId, keyB, gB);
    bool attachAOnB = cat.attachMolecularGraph(refB.generatedId, keyA, gA);
    recA = cat.record(refA.generatedId);
    recB = cat.record(refB.generatedId);
    emit("graph_a_attaches_only_to_identity_a",
        attachA && !attachBOnA && recA && recA->hasMolecularGraph
            && molecularGraphsStoredEqual(recA->molecularGraph, gA), "");
    emit("graph_b_attaches_only_to_identity_b",
        attachB && !attachAOnB && recB && recB->hasMolecularGraph
            && molecularGraphsStoredEqual(recB->molecularGraph, gB)
            && recA && recA->recordId != recB->recordId, "");

    bool mismatch = cat.attachMolecularGraph(refA.generatedId, keyA, h2o);
    emit("mismatched_graph_identity_rejected", !mismatch, "");

    recA = cat.record(refA.generatedId);
    recB = cat.record(refB.generatedId);
    emit("graph_attachment_does_not_change_signature",
        recA && recB && recA->canonicalSignature == sigA && recB->canonicalSignature == sigB, "");
    emit("graph_attachment_does_not_change_record_id",
        recA && recB && recA->recordId == ridA && recB->recordId == ridB, "");
    emit("graph_attachment_does_not_change_element_hex",
        recA && recB && recA->displayOrdinal == ordA && recB->displayOrdinal == ordB, "");

    SaceScalarProperty mmFromGraph{};
    bool massAgree = recA
        && recA->properties.molarMassGPerMol.known
        && recA->properties.molarMassGPerMol.source == SacePropertySource::IdentityDerived
        && recA->properties.molarMassGPerMol.confidence == SaceConfidence::High
        && saceDeriveMolarMass(compA, mmFromGraph)
        && std::fabs(mmFromGraph.value - recA->properties.molarMassGPerMol.value) < 1.0e-6f
        && std::fabs(recA->properties.molarMassGPerMol.value - massA) < 1.0e-6f;
    emit("graph_composition_molar_mass_agrees", massAgree, "");
    emit("attached_descriptors_from_stored_graph",
        recA && recA->hasMolecularDescriptors
            && recA->molecularDescriptors.atomCount == atomCount(recA->molecularGraph)
            && recA->properties.molarMassGPerMol.value == massA, "");

    emit("debug_serialization_is_not_catalog_key",
        recA && recA->canonicalSignature.find("debug|") == std::string::npos, "");

    cat.clear();
    emit("molecule_diag_clears_catalog", cat.generatedCount() == 0, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
