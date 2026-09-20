#include "chemistry/SaceIdentity.h"

#include "chemistry/ReactionRegistry.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceTypes.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

CompositionNormalizeResult normalizeElementalComposition(ElementalComposition &c) {
    if (c.count > kMaxElementalSpecies)
        return CompositionNormalizeResult::Overflow;
    ElementCount merged[kMaxElementalSpecies]{};
    int nMerged = 0;
    int nIn = static_cast<int>(c.count);
    for (int i = 0; i < nIn; ++i) {
        AtomicNumber z = c.entries[i].atomicNumber;
        uint16_t n = c.entries[i].count;
        if (n == 0) continue;
        if (!validAtomicNumber(z))
            return CompositionNormalizeResult::Invalid;
        int slot = -1;
        for (int m = 0; m < nMerged; ++m) {
            if (merged[m].atomicNumber == z) { slot = m; break; }
        }
        if (slot < 0) {
            if (nMerged >= kMaxElementalSpecies)
                return CompositionNormalizeResult::Overflow;
            merged[nMerged] = {z, n};
            ++nMerged;
        } else {
            uint32_t sum = static_cast<uint32_t>(merged[slot].count) + n;
            if (sum > 0xFFFFu)
                return CompositionNormalizeResult::Overflow;
            merged[slot].count = static_cast<uint16_t>(sum);
        }
    }
    for (int a = 0; a < nMerged; ++a) {
        int best = a;
        for (int b = a + 1; b < nMerged; ++b)
            if (merged[b].atomicNumber < merged[best].atomicNumber) best = b;
        if (best != a) {
            ElementCount tmp = merged[a];
            merged[a] = merged[best];
            merged[best] = tmp;
        }
    }
    c = {};
    c.count = static_cast<uint8_t>(nMerged);
    for (int i = 0; i < nMerged; ++i)
        c.entries[i] = merged[i];
    return CompositionNormalizeResult::Ok;
}

bool addScaledComposition(ElementalInventory &acc, ElementalComposition const &src, int64_t scale) {
    if (src.count > kMaxElementalSpecies) return false;
    if (acc.count > kMaxElementalInventory) return false;
    if (scale == 0) return true;
    for (int i = 0; i < src.count; ++i) {
        AtomicNumber z = src.entries[i].atomicNumber;
        uint16_t n = src.entries[i].count;
        if (n == 0) continue;
        if (!validAtomicNumber(z)) return false;
        int slot = -1;
        for (int m = 0; m < acc.count && m < kMaxElementalInventory; ++m) {
            if (acc.entries[m].atomicNumber == z) { slot = m; break; }
        }
        int64_t delta = scale * static_cast<int64_t>(n);
        if (slot < 0) {
            if (acc.count >= kMaxElementalInventory) return false;
            acc.entries[acc.count] = {z, delta};
            ++acc.count;
        } else {
            acc.entries[slot].count += delta;
        }
    }
    return true;
}

bool elementalInventoriesEqual(ElementalInventory const &a, ElementalInventory const &b) {
    if (a.count > kMaxElementalInventory || b.count > kMaxElementalInventory)
        return false;
    int64_t counts[kMaxAtomicNumber + 1]{};
    auto fold = [&](ElementalInventory const &inv, int sign) {
        for (int i = 0; i < inv.count; ++i) {
            AtomicNumber z = inv.entries[i].atomicNumber;
            if (!validAtomicNumber(z)) return false;
            counts[z] += static_cast<int64_t>(sign) * inv.entries[i].count;
        }
        return true;
    };
    if (!fold(a, 1) || !fold(b, -1)) return false;
    for (int z = kMinAtomicNumber; z <= kMaxAtomicNumber; ++z)
        if (counts[z] != 0) return false;
    return true;
}

namespace {

bool appendStr(char *out, int cap, int &used, char const *s) {
    if (!out || cap <= 0 || !s) return false;
    while (*s) {
        if (used >= cap - 1) return false;
        out[used++] = *s++;
    }
    out[used] = '\0';
    return true;
}

bool appendU32(char *out, int cap, int &used, unsigned v) {
    char buf[16];
    int n = std::snprintf(buf, sizeof(buf), "%u", v);
    if (n <= 0) return false;
    return appendStr(out, cap, used, buf);
}

} // namespace

ChemicalIdentity saceExactIdentity(ChemicalRepresentationKind kind, char const *formula,
    char const *structureKey, ElementCount const *items, int n)
{
    if (n < 0 || n > kMaxElementalSpecies)
        return saceUnknownIdentity();
    if (n > 0 && !items)
        return saceUnknownIdentity();
    ChemicalIdentity id{};
    id.kind = kind;
    id.formula = formula;
    id.structureKey = structureKey;
    id.exact = true;
    for (int i = 0; i < n; ++i)
        id.elemental.entries[i] = items[i];
    id.elemental.count = static_cast<uint8_t>(n);
    if (normalizeElementalComposition(id.elemental) != CompositionNormalizeResult::Ok)
        return saceUnknownIdentity();
    return id;
}

namespace {

bool stoichInteger(float coefficient, int64_t &out) {
    if (!std::isfinite(coefficient) || !(coefficient > 0.0f)) return false;
    double r = std::round(static_cast<double>(coefficient));
    if (std::abs(static_cast<double>(coefficient) - r) > 1.0e-6) return false;
    if (r < 1.0 || r > 1.0e9) return false;
    out = static_cast<int64_t>(r);
    return true;
}

bool accumulateReactionSide(ReactionParticipant const *parts, uint8_t nParts,
    ElementalInventory &inv)
{
    if (nParts > kMaxReactionParticipants) return false;
    for (int i = 0; i < nParts; ++i) {
        ReactionParticipant const &p = parts[i];
        if (!reactionParticipantUsed(p)) continue;
        if (!validSubstance(p.substance)) return false;
        ChemicalIdentity const &chem = substanceDef(p.substance).chemicalIdentity;
        if (!chemicalIdentityExact(chem)) return false;
        int64_t k = 0;
        if (!stoichInteger(p.coefficient, k)) return false;
        if (!addScaledComposition(inv, chem.elemental, k)) return false;
    }
    return true;
}

} // namespace

bool writeChemicalSignature(ChemicalIdentity const &id, char *out, int cap) {
    if (!out || cap <= 0) return false;
    out[0] = '\0';
    ChemicalIdentity norm = id;
    if (norm.elemental.count > kMaxElementalSpecies)
        return false;
    if (normalizeElementalComposition(norm.elemental) != CompositionNormalizeResult::Ok)
        return false;
    int used = 0;
    auto fail = [&]() {
        out[0] = '\0';
        return false;
    };
    if (!appendStr(out, cap, used, chemicalRepresentationKindKey(norm.kind))) return fail();
    if (!appendStr(out, cap, used, "|")) return fail();
    for (int i = 0; i < norm.elemental.count; ++i) {
        if (i > 0 && !appendStr(out, cap, used, ",")) return fail();
        if (!appendU32(out, cap, used, norm.elemental.entries[i].atomicNumber)) return fail();
        if (!appendStr(out, cap, used, ":")) return fail();
        if (!appendU32(out, cap, used, norm.elemental.entries[i].count)) return fail();
    }
    if (!appendStr(out, cap, used, "|")) return fail();
    if (!appendStr(out, cap, used, saceStructureKeyOrEmpty(norm.structureKey))) return fail();
    return true;
}

bool chemicalIdentitiesEquivalent(ChemicalIdentity const &a, ChemicalIdentity const &b) {
    if (a.elemental.count > kMaxElementalSpecies || b.elemental.count > kMaxElementalSpecies)
        return false;
    ChemicalIdentity na = a;
    ChemicalIdentity nb = b;
    if (normalizeElementalComposition(na.elemental) != CompositionNormalizeResult::Ok)
        return false;
    if (normalizeElementalComposition(nb.elemental) != CompositionNormalizeResult::Ok)
        return false;
    if (na.kind != nb.kind) return false;
    if (!elementalCompositionsEqual(na.elemental, nb.elemental)) return false;
    if (!saceStructureKeysEqual(na.structureKey, nb.structureKey)) return false;
    return true;
}

SubstanceId findBuiltInByChemicalIdentity(ChemicalIdentity const &query) {
    if (query.elemental.count > kMaxElementalSpecies)
        return SUBSTANCE_NONE;
    ChemicalIdentity q = query;
    if (normalizeElementalComposition(q.elemental) != CompositionNormalizeResult::Ok)
        return SUBSTANCE_NONE;
    if (!chemicalIdentityExact(q))
        return SUBSTANCE_NONE;
    for (SubstanceId id = 1; id < SUBSTANCE_COUNT; ++id) {
        ChemicalIdentity const &built = substanceDef(id).chemicalIdentity;
        if (!chemicalIdentityExact(built)) continue;
        if (chemicalIdentitiesEquivalent(q, built))
            return id;
    }
    return SUBSTANCE_NONE;
}

ConservationCheck reactionAtomConservation(ReactionDefinition const &def) {
    if (def.reactantCount > kMaxReactionParticipants
        || def.productCount > kMaxReactionParticipants)
        return ConservationCheck::Unknown;
    ElementalInventory react{};
    ElementalInventory prod{};
    if (!accumulateReactionSide(def.reactants, def.reactantCount, react))
        return ConservationCheck::Unknown;
    if (!accumulateReactionSide(def.products, def.productCount, prod))
        return ConservationCheck::Unknown;
    if (react.count == 0 && prod.count == 0)
        return ConservationCheck::Unknown;
    return elementalInventoriesEqual(react, prod)
        ? ConservationCheck::Balanced
        : ConservationCheck::Unbalanced;
}

ChemicalIdentity saceUnknownIdentity(ChemicalRepresentationKind kind) {
    ChemicalIdentity id{};
    id.kind = kind;
    id.exact = false;
    return id;
}

ChemicalIdentity saceBuiltinWaterIdentity() {
    ElementCount items[] = {{kAtomicHydrogen, 2}, {kAtomicOxygen, 1}};
    return saceExactIdentity(ChemicalRepresentationKind::SmallMolecule, "H2O", "water",
        items, 2);
}

ChemicalIdentity saceBuiltinHydrogenIdentity() {
    ElementCount items[] = {{kAtomicHydrogen, 2}};
    return saceExactIdentity(ChemicalRepresentationKind::SmallMolecule, "H2",
        "molecular-hydrogen", items, 1);
}

ChemicalIdentity saceBuiltinOxygenIdentity() {
    ElementCount items[] = {{kAtomicOxygen, 2}};
    return saceExactIdentity(ChemicalRepresentationKind::SmallMolecule, "O2",
        "molecular-oxygen", items, 1);
}

ChemicalIdentity saceBuiltinCarbonIdentity() {
    ElementCount items[] = {{kAtomicCarbon, 1}};
    return saceExactIdentity(ChemicalRepresentationKind::AtomicSpecies, "C",
        "elemental-carbon", items, 1);
}

ChemicalIdentity saceBuiltinCarbonDioxideIdentity() {
    ElementCount items[] = {{kAtomicCarbon, 1}, {kAtomicOxygen, 2}};
    return saceExactIdentity(ChemicalRepresentationKind::SmallMolecule, "CO2",
        "carbon-dioxide", items, 2);
}

void runSaceIdentityDiagnostics() {
    std::ofstream out(miscFile("sace_identity_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    ChemicalIdentity water = saceBuiltinWaterIdentity();
    ChemicalIdentity hydrogen = saceBuiltinHydrogenIdentity();
    ChemicalIdentity oxygen = saceBuiltinOxygenIdentity();
    ChemicalIdentity carbon = saceBuiltinCarbonIdentity();
    ChemicalIdentity co2 = saceBuiltinCarbonDioxideIdentity();

    emit("water_elemental_H2_O1",
        elementalCountOf(water.elemental, kAtomicHydrogen) == 2
            && elementalCountOf(water.elemental, kAtomicOxygen) == 1
            && water.elemental.count == 2
            && water.exact,
        "H=" + std::to_string(elementalCountOf(water.elemental, kAtomicHydrogen))
            + " O=" + std::to_string(elementalCountOf(water.elemental, kAtomicOxygen)));
    emit("hydrogen_elemental_H2",
        elementalCountOf(hydrogen.elemental, kAtomicHydrogen) == 2
            && hydrogen.elemental.count == 1 && hydrogen.exact, "");
    emit("oxygen_elemental_O2",
        elementalCountOf(oxygen.elemental, kAtomicOxygen) == 2
            && oxygen.elemental.count == 1 && oxygen.exact, "");
    emit("carbon_elemental_C1",
        elementalCountOf(carbon.elemental, kAtomicCarbon) == 1
            && carbon.elemental.count == 1 && carbon.exact, "");
    emit("co2_elemental_C1_O2",
        elementalCountOf(co2.elemental, kAtomicCarbon) == 1
            && elementalCountOf(co2.elemental, kAtomicOxygen) == 2
            && co2.elemental.count == 2 && co2.exact, "");

    char sigA[kSaceSignatureCap]{};
    char sigB[kSaceSignatureCap]{};
    bool wroteA = writeChemicalSignature(water, sigA, kSaceSignatureCap);
    bool wroteB = writeChemicalSignature(substanceDef(SUBSTANCE_WATER).chemicalIdentity, sigB, kSaceSignatureCap);
    emit("water_canonical_key_deterministic",
        wroteA && wroteB && std::strcmp(sigA, sigB) == 0 && std::strstr(sigA, "water") != nullptr
            && std::strstr(sigA, "H2O") == nullptr,
        sigA);

    ChemicalIdentity shuffled{};
    shuffled.kind = ChemicalRepresentationKind::SmallMolecule;
    shuffled.structureKey = "water";
    shuffled.exact = true;
    shuffled.elemental.entries[0] = {kAtomicOxygen, 1};
    shuffled.elemental.entries[1] = {kAtomicHydrogen, 2};
    shuffled.elemental.count = 2;
    emit("reordered_HO_normalizes_to_water",
        normalizeElementalComposition(shuffled.elemental) == CompositionNormalizeResult::Ok
            && findBuiltInByChemicalIdentity(shuffled) == SUBSTANCE_WATER
            && chemicalIdentitiesEquivalent(shuffled, water),
        "H=" + std::to_string(elementalCountOf(shuffled.elemental, kAtomicHydrogen))
            + " O=" + std::to_string(elementalCountOf(shuffled.elemental, kAtomicOxygen)));

    ChemicalIdentity isomer = water;
    isomer.structureKey = "not-water";
    isomer.formula = "H2O";
    emit("different_structure_key_is_not_water",
        findBuiltInByChemicalIdentity(isomer) == SUBSTANCE_NONE
            && !chemicalIdentitiesEquivalent(isomer, water),
        "");

    emit("lookup_water",
        findBuiltInByChemicalIdentity(water) == SUBSTANCE_WATER, "");
    emit("lookup_co2",
        findBuiltInByChemicalIdentity(co2) == SUBSTANCE_CARBON_DIOXIDE, "");
    emit("lookup_hydrogen",
        findBuiltInByChemicalIdentity(hydrogen) == SUBSTANCE_HYDROGEN, "");
    emit("lookup_oxygen",
        findBuiltInByChemicalIdentity(oxygen) == SUBSTANCE_OXYGEN, "");
    emit("lookup_carbon",
        findBuiltInByChemicalIdentity(carbon) == SUBSTANCE_CARBON, "");

    ChemicalIdentity const &honeyId = substanceDef(SUBSTANCE_HONEY).chemicalIdentity;
    emit("honey_not_exact_molecule",
        !honeyId.exact
            && honeyId.kind == ChemicalRepresentationKind::Mixture
            && honeyId.elemental.count == 0
            && findBuiltInByChemicalIdentity(honeyId) == SUBSTANCE_NONE,
        conservationCheckKey(ConservationCheck::Unknown));

    ChemicalIdentity const &airId = substanceDef(SUBSTANCE_AIR).chemicalIdentity;
    emit("air_is_mixture_not_o2_n2",
        !airId.exact
            && airId.kind == ChemicalRepresentationKind::Mixture
            && airId.elemental.count == 0
            && findBuiltInByChemicalIdentity(airId) == SUBSTANCE_NONE
            && findBuiltInByChemicalIdentity(oxygen) == SUBSTANCE_OXYGEN,
        "");

    ConservationCheck h2o2 = reactionAtomConservation(reactionDef(REACTION_HYDROGEN_COMBUSTION));
    emit("h2_o2_atom_balanced", h2o2 == ConservationCheck::Balanced, conservationCheckKey(h2o2));
    ConservationCheck co2rx = reactionAtomConservation(reactionDef(REACTION_CARBON_COMBUSTION));
    emit("c_o2_atom_balanced", co2rx == ConservationCheck::Balanced, conservationCheckKey(co2rx));

    ReactionDefinition unbalanced{};
    unbalanced.reactants[0] = {SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0f};
    unbalanced.reactants[1] = {SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0f};
    unbalanced.reactantCount = 2;
    unbalanced.products[0] = {SUBSTANCE_WATER, MatterPhase::Gas, 1.0f};
    unbalanced.productCount = 1;
    ConservationCheck unb = reactionAtomConservation(unbalanced);
    emit("synthetic_unbalanced", unb == ConservationCheck::Unbalanced, conservationCheckKey(unb));

    ReactionDefinition unknownRx{};
    unknownRx.reactants[0] = {SUBSTANCE_HONEY, MatterPhase::Liquid, 1.0f};
    unknownRx.reactantCount = 1;
    unknownRx.products[0] = {SUBSTANCE_WATER, MatterPhase::Liquid, 1.0f};
    unknownRx.productCount = 1;
    ConservationCheck unk = reactionAtomConservation(unknownRx);
    emit("unknown_composition_reaction", unk == ConservationCheck::Unknown, conservationCheckKey(unk));

    emit("formula_is_not_identity",
        water.formula && std::strcmp(water.formula, "H2O") == 0
            && isomer.formula && std::strcmp(isomer.formula, "H2O") == 0
            && findBuiltInByChemicalIdentity(isomer) == SUBSTANCE_NONE, "");
    emit("display_name_not_in_signature",
        std::strstr(sigA, "Water") == nullptr
            && std::strstr(sigA, "ins_mat_water") == nullptr, sigA);

    emit("no_runtime_generated_ids", SUBSTANCE_COUNT == 12, std::to_string(SUBSTANCE_COUNT));

    ElementalComposition overCap{};
    overCap.count = static_cast<uint8_t>(kMaxElementalSpecies + 1);
    overCap.entries[0] = {kAtomicHydrogen, 1};
    CompositionNormalizeResult overNorm = normalizeElementalComposition(overCap);
    emit("elemental_count_over_capacity",
        overNorm == CompositionNormalizeResult::Overflow
            && overCap.count == static_cast<uint8_t>(kMaxElementalSpecies + 1),
        "");

    ElementCount tooMany[kMaxElementalSpecies + 1]{};
    for (int i = 0; i < kMaxElementalSpecies + 1; ++i)
        tooMany[i] = {static_cast<AtomicNumber>(i + 1), 1};
    ChemicalIdentity overId = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "X", "synthetic-over-capacity", tooMany, kMaxElementalSpecies + 1);
    emit("exact_constructor_over_capacity",
        !overId.exact && overId.elemental.count == 0, "");

    char tiny[8]{};
    bool tinyOk = writeChemicalSignature(water, tiny, 8);
    emit("signature_small_buffer", !tinyOk && tiny[0] == '\0', "");

    ReactionDefinition malformedCount{};
    malformedCount.reactants[0] = {SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0f};
    malformedCount.reactants[1] = {SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0f};
    malformedCount.reactantCount = static_cast<uint8_t>(kMaxReactionParticipants + 1);
    malformedCount.products[0] = {SUBSTANCE_WATER, MatterPhase::Gas, 2.0f};
    malformedCount.productCount = 1;
    ConservationCheck malformedRx = reactionAtomConservation(malformedCount);
    emit("malformed_reaction_participant_count",
        malformedRx == ConservationCheck::Unknown, conservationCheckKey(malformedRx));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
