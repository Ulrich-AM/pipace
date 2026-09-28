#include "chemistry/SaceCatalog.h"

#include "chemistry/SaceEstimation.h"
#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"
#include "substance/GeneratedMaterialRegistry.h"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <fstream>

bool formatGeneratedElementName(uint32_t displayOrdinal, char *out, int cap) {
    if (!out || cap <= 0) return false;
    out[0] = '\0';
    if (displayOrdinal == 0) return false;
    int n = std::snprintf(out, static_cast<std::size_t>(cap), "Element #%X", displayOrdinal);
    return n > 0 && n < cap;
}

SaceCatalog &saceGeneratedCatalog() {
    static SaceCatalog catalog;
    return catalog;
}

void SaceCatalog::clear() {
    resetGeneratedRuntimeMaterials();
    records_.clear();
    bySignature_.clear();
    nextDisplayOrdinal_ = 1;
}

SaceGeneratedRecord const *SaceCatalog::record(SaceRecordId id) const {
    if (id == kSaceRecordNone || id > records_.size())
        return nullptr;
    SaceGeneratedRecord const &rec = records_[id - 1];
    if (rec.recordId != id) return nullptr;
    return &rec;
}

SaceGeneratedRecord *SaceCatalog::recordMutable(SaceRecordId id) {
    if (id == kSaceRecordNone || id > records_.size())
        return nullptr;
    SaceGeneratedRecord &rec = records_[id - 1];
    if (rec.recordId != id) return nullptr;
    return &rec;
}

void SaceCatalog::applyStructuralThermoEstimates(SaceGeneratedRecord &record,
    SaceJobackEstimateBundle const &bundle)
{
    saceAssignScalarProperty(record.properties.normalBoilingPointK, bundle.normalBoilingPointK);
    saceAssignScalarProperty(record.properties.criticalTemperatureK, bundle.criticalTemperatureK);
    saceAssignScalarProperty(record.properties.criticalPressurePa, bundle.criticalPressurePa);
    saceAssignScalarProperty(record.properties.criticalMolarVolumeM3PerMol,
        bundle.criticalMolarVolumeM3PerMol);
    saceAssignScalarProperty(record.properties.acentricFactor, bundle.acentricFactor);
    saceAssignScalarProperty(record.properties.enthalpyVaporizationAtNormalBoilingJPerMol,
        bundle.enthalpyVaporizationAtNormalBoilingJPerMol);
    saceAssignScalarProperty(record.properties.idealGasHeatCapacityAt298KJPerMolK,
        bundle.idealGasHeatCapacityAt298KJPerMolK);
    saceAssignScalarProperty(record.properties.saturatedLiquidHeatCapacityAt298KJPerMolK,
        bundle.saturatedLiquidHeatCapacityAt298KJPerMolK);
    saceAssignScalarProperty(record.properties.saturatedLiquidDensityAt298KKgPerM3,
        bundle.saturatedLiquidDensityAt298KKgPerM3);
    saceAssignScalarProperty(record.properties.liquidDynamicViscosityAt298KPaS,
        bundle.liquidDynamicViscosityAt298KPaS);
    saceAssignScalarProperty(record.properties.liquidSurfaceTensionAt298KNPerM,
        bundle.liquidSurfaceTensionAt298KNPerM);
    saceAssignScalarProperty(record.properties.liquidThermalConductivityAt298KWPerMK,
        bundle.liquidThermalConductivityAt298KWPerMK);
    saceAssignScalarProperty(record.properties.gasThermalConductivityAt298KWPerMK,
        bundle.gasThermalConductivityAt298KWPerMK);
}

ChemicalIdentity SaceCatalog::identityView(SaceGeneratedRecord const &rec) const {
    ChemicalIdentity id{};
    id.kind = rec.kind;
    id.elemental = rec.elemental;
    id.formula = rec.formula.empty() ? nullptr : rec.formula.c_str();
    id.structureKey = rec.structureKey.empty() ? nullptr : rec.structureKey.c_str();
    id.exact = rec.exactIdentity;
    return id;
}

bool SaceCatalog::eligibleForGeneratedIdentity(ChemicalIdentity const &id, char *signature, int cap) const {
    if (!signature || cap <= 0) return false;
    if (!id.exact) return false;
    if (id.elemental.count == 0 || id.elemental.count > kMaxElementalSpecies)
        return false;
    if (id.kind != ChemicalRepresentationKind::AtomicSpecies
        && id.kind != ChemicalRepresentationKind::SmallMolecule
        && id.kind != ChemicalRepresentationKind::IonicMaterial)
        return false;
    if ((id.kind == ChemicalRepresentationKind::SmallMolecule
            || id.kind == ChemicalRepresentationKind::IonicMaterial)
        && saceStructureKeyOrEmpty(id.structureKey)[0] == '\0')
        return false;
    ChemicalIdentity norm = id;
    if (normalizeElementalComposition(norm.elemental) != CompositionNormalizeResult::Ok)
        return false;
    if (!chemicalIdentityExact(norm))
        return false;
    if (!writeChemicalSignature(norm, signature, cap))
        return false;
    return signature[0] != '\0';
}

SaceRecordId SaceCatalog::insertGenerated(ChemicalIdentity const &id, char const *signature) {
    SaceGeneratedRecord rec{};
    rec.recordId = static_cast<SaceRecordId>(records_.size() + 1);
    rec.kind = id.kind;
    rec.elemental = id.elemental;
    if (normalizeElementalComposition(rec.elemental) != CompositionNormalizeResult::Ok)
        return kSaceRecordNone;
    rec.formula = id.formula ? id.formula : "";
    rec.structureKey = saceStructureKeyOrEmpty(id.structureKey);
    rec.canonicalSignature = signature ? signature : "";
    rec.displayOrdinal = nextDisplayOrdinal_++;
    rec.exactIdentity = true;
    rec.spawnable = false;
    initializeGeneratedProperties(rec.properties, rec.elemental);
    records_.push_back(std::move(rec));
    SaceGeneratedRecord &stored = records_.back();
    bySignature_.emplace(stored.canonicalSignature, stored.recordId);
    return stored.recordId;
}

SaceSubstanceRef SaceCatalog::resolve(ChemicalIdentity const &query, bool create) {
    SaceSubstanceRef result{};
    ChemicalIdentity norm = query;
    if (norm.elemental.count > kMaxElementalSpecies)
        return result;
    if (normalizeElementalComposition(norm.elemental) != CompositionNormalizeResult::Ok)
        return result;

    SubstanceId builtin = findBuiltInByChemicalIdentity(norm);
    if (builtin != SUBSTANCE_NONE) {
        result.kind = SaceResolutionKind::BuiltIn;
        result.builtinId = builtin;
        return result;
    }

    char signature[kSaceSignatureCap]{};
    if (!eligibleForGeneratedIdentity(norm, signature, kSaceSignatureCap))
        return result;

    auto found = bySignature_.find(signature);
    if (found != bySignature_.end()) {
        result.kind = SaceResolutionKind::Generated;
        result.generatedId = found->second;
        return result;
    }
    if (!create)
        return result;

    SaceRecordId id = insertGenerated(norm, signature);
    if (id == kSaceRecordNone)
        return result;
    result.kind = SaceResolutionKind::Generated;
    result.generatedId = id;
    return result;
}

bool SaceCatalog::attachMolecularGraph(SaceRecordId id, char const *claimedStructureKey,
    SaceMolecularGraph const &graph)
{
    if (id == kSaceRecordNone || id > records_.size())
        return false;
    SaceGeneratedRecord &rec = records_[id - 1];
    if (rec.recordId != id) return false;
    if (!rec.exactIdentity) return false;
    if (rec.kind != ChemicalRepresentationKind::SmallMolecule)
        return false;
    if (!claimedStructureKey || claimedStructureKey[0] == '\0')
        return false;
    if (rec.structureKey != claimedStructureKey)
        return false;
    if (validateMolecularGraph(graph, true) != SaceGraphValidation::Ok)
        return false;
    ChemicalIdentity view = identityView(rec);
    if (!graphMatchesChemicalIdentity(graph, view))
        return false;
    SaceMolecularDescriptors desc{};
    if (!deriveMolecularDescriptors(graph, desc))
        return false;
    SaceFunctionalProfile func{};
    if (!deriveFunctionalProfile(graph, func))
        return false;
    SaceJobackEstimateBundle joback{};
    (void)saceEstimateJobackBundle(graph, joback);
    if (rec.hasMolecularGraph) {
        if (rec.molecularGraphStructureKey != rec.structureKey)
            return false;
        if (!molecularGraphsStoredEqual(rec.molecularGraph, graph))
            return false;
        applyStructuralThermoEstimates(rec, joback);
        return true;
    }
    rec.molecularGraph = graph;
    rec.molecularGraphStructureKey = rec.structureKey;
    rec.molecularDescriptors = desc;
    rec.functionalProfile = func;
    rec.hasMolecularGraph = true;
    rec.hasMolecularDescriptors = true;
    rec.hasFunctionalProfile = true;
    applyStructuralThermoEstimates(rec, joback);
    return true;
}

void runSaceCatalogDiagnostics() {
    std::ofstream out(miscFile("sace_catalog_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();
    emit("fresh_catalog_empty", cat.generatedCount() == 0, "");

    ChemicalIdentity water = saceBuiltinWaterIdentity();
    std::size_t beforeWater = cat.generatedCount();
    SaceSubstanceRef waterRef = cat.resolve(water, true);
    emit("water_resolves_builtin",
        waterRef.kind == SaceResolutionKind::BuiltIn
            && waterRef.builtinId == SUBSTANCE_WATER
            && waterRef.generatedId == kSaceRecordNone, "");
    emit("water_creates_zero_generated",
        cat.generatedCount() == beforeWater, std::to_string(cat.generatedCount()));

    ChemicalIdentity co2 = saceBuiltinCarbonDioxideIdentity();
    SaceSubstanceRef co2Ref = cat.resolve(co2, true);
    emit("co2_resolves_builtin",
        co2Ref.kind == SaceResolutionKind::BuiltIn
            && co2Ref.builtinId == SUBSTANCE_CARBON_DIOXIDE, "");

    ElementCount c2h6o[] = {
        {kAtomicCarbon, 2}, {kAtomicHydrogen, 6}, {kAtomicOxygen, 1}
    };
    ChemicalIdentity synthA = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-a", c2h6o, 3);
    std::size_t beforeGen = cat.generatedCount();
    SaceSubstanceRef aRef = cat.resolve(synthA, true);
    SaceGeneratedRecord const *aRec = cat.record(aRef.generatedId);
    emit("synthetic_creates_one_generated",
        aRef.kind == SaceResolutionKind::Generated
            && aRef.generatedId != kSaceRecordNone
            && cat.generatedCount() == beforeGen + 1
            && aRec != nullptr, "");
    emit("generated_record_not_spawnable",
        aRec
            && aRec->recordId != kSaceRecordNone
            && aRec->displayOrdinal != 0
            && aRec->spawnable == false
            && aRec->exactIdentity, "");

    char display[32]{};
    bool named = aRec && formatGeneratedElementName(aRec->displayOrdinal, display, 32);
    emit("display_name_element_hex",
        named && std::strncmp(display, "Element #", 9) == 0
            && display[9] != '\0'
            && std::strchr(display, '#') != nullptr, display);

    std::size_t afterFirst = cat.generatedCount();
    SaceSubstanceRef aAgain = cat.resolve(synthA, true);
    emit("same_identity_same_handle",
        aAgain.kind == SaceResolutionKind::Generated
            && aAgain.generatedId == aRef.generatedId, "");
    emit("same_identity_does_not_increment",
        cat.generatedCount() == afterFirst, std::to_string(cat.generatedCount()));

    ChemicalIdentity shuffled = synthA;
    shuffled.elemental = {};
    shuffled.elemental.entries[0] = {kAtomicOxygen, 1};
    shuffled.elemental.entries[1] = {kAtomicCarbon, 2};
    shuffled.elemental.entries[2] = {kAtomicHydrogen, 6};
    shuffled.elemental.count = 3;
    SaceSubstanceRef shuffledRef = cat.resolve(shuffled, true);
    emit("reordered_entries_same_record",
        shuffledRef.kind == SaceResolutionKind::Generated
            && shuffledRef.generatedId == aRef.generatedId, "");

    ChemicalIdentity formulaVariant = synthA;
    formulaVariant.formula = "EtOH";
    SaceSubstanceRef formulaRef = cat.resolve(formulaVariant, true);
    emit("formula_text_ignored",
        formulaRef.kind == SaceResolutionKind::Generated
            && formulaRef.generatedId == aRef.generatedId, "");

    ChemicalIdentity synthB = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-b", c2h6o, 3);
    SaceSubstanceRef bRef = cat.resolve(synthB, true);
    SaceGeneratedRecord const *bRec = cat.record(bRef.generatedId);
    emit("different_structure_different_record",
        bRef.kind == SaceResolutionKind::Generated
            && bRef.generatedId != kSaceRecordNone
            && bRef.generatedId != aRef.generatedId
            && cat.generatedCount() == afterFirst + 1, "");

    char displayB[32]{};
    bool namedB = bRec && formatGeneratedElementName(bRec->displayOrdinal, displayB, 32);
    emit("different_element_hex",
        namedB && aRec && bRec
            && aRec->displayOrdinal != bRec->displayOrdinal
            && std::strcmp(display, displayB) != 0, displayB);

    ChemicalIdentity ambiguous = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "", c2h6o, 3);
    std::size_t beforeAmb = cat.generatedCount();
    SaceSubstanceRef ambRef = cat.resolve(ambiguous, true);
    emit("ambiguous_molecule_rejected",
        ambRef.kind == SaceResolutionKind::Invalid
            && ambRef.generatedId == kSaceRecordNone
            && cat.generatedCount() == beforeAmb, "");

    ChemicalIdentity inexact = saceUnknownIdentity(ChemicalRepresentationKind::Mixture);
    inexact.exact = false;
    inexact.elemental.entries[0] = {kAtomicCarbon, 1};
    inexact.elemental.count = 1;
    std::size_t beforeInexact = cat.generatedCount();
    SaceSubstanceRef inexactRef = cat.resolve(inexact, true);
    emit("inexact_identity_rejected",
        inexactRef.kind == SaceResolutionKind::Invalid
            && cat.generatedCount() == beforeInexact, "");

    ChemicalIdentity malformed = synthA;
    malformed.elemental.count = static_cast<uint8_t>(kMaxElementalSpecies + 1);
    std::size_t beforeBad = cat.generatedCount();
    SaceSubstanceRef badRef = cat.resolve(malformed, true);
    emit("malformed_elemental_count_rejected",
        badRef.kind == SaceResolutionKind::Invalid
            && cat.generatedCount() == beforeBad, "");

    char tiny[8]{};
    bool tinyOk = writeChemicalSignature(synthA, tiny, 8);
    emit("tiny_signature_buffer_fails", !tinyOk && tiny[0] == '\0', "");

    bool truncatedInCache = false;
    for (std::size_t i = 0; i < cat.generatedCount(); ++i) {
        SaceGeneratedRecord const *rec = cat.record(static_cast<SaceRecordId>(i + 1));
        if (!rec) { truncatedInCache = true; break; }
        char full[kSaceSignatureCap]{};
        ChemicalIdentity view = cat.identityView(*rec);
        if (!writeChemicalSignature(view, full, kSaceSignatureCap)
            || rec->canonicalSignature != full
            || rec->canonicalSignature.empty())
            truncatedInCache = true;
    }
    emit("no_truncated_signature_in_cache", !truncatedInCache, "");

    SaceGeneratedRecord const *stableA = cat.record(aRef.generatedId);
    SaceRecordId stableId = stableA ? stableA->recordId : kSaceRecordNone;
    uint32_t stableOrdinal = stableA ? stableA->displayOrdinal : 0;
    std::string stableSig = stableA ? stableA->canonicalSignature : "";
    for (int i = 0; i < 256; ++i) {
        char key[64];
        std::snprintf(key, sizeof(key), "synthetic-stability-%d", i);
        ChemicalIdentity extra = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
            "C2H6O", key, c2h6o, 3);
        cat.resolve(extra, true);
    }
    emit("record_pointer_stable_across_insertions",
        stableA
            && stableA->recordId == stableId
            && stableA->displayOrdinal == stableOrdinal
            && stableA->canonicalSignature == stableSig
            && cat.record(aRef.generatedId) == stableA
            && cat.generatedCount() >= 256, "");

    SaceGeneratedRecord const *massA = cat.record(aRef.generatedId);
    SaceGeneratedRecord const *massB = cat.record(bRef.generatedId);
    emit("c2h6o_a_molar_mass_known",
        massA && massA->properties.molarMassGPerMol.known
            && std::fabs(massA->properties.molarMassGPerMol.value - 46.069f) < 0.002f, "");
    emit("c2h6o_b_molar_mass_known",
        massB && massB->properties.molarMassGPerMol.known
            && std::fabs(massB->properties.molarMassGPerMol.value - 46.069f) < 0.002f, "");
    emit("c2h6o_masses_equal",
        massA && massB
            && massA->properties.molarMassGPerMol.known
            && massB->properties.molarMassGPerMol.known
            && std::fabs(massA->properties.molarMassGPerMol.value
                - massB->properties.molarMassGPerMol.value) < 1.0e-6f, "");
    emit("c2h6o_chemically_distinct",
        massA && massB
            && massA->recordId != massB->recordId
            && massA->canonicalSignature != massB->canonicalSignature, "");
    emit("molar_mass_source_identity_derived",
        massA
            && massA->properties.molarMassGPerMol.source == SacePropertySource::IdentityDerived, "");
    emit("molar_mass_confidence_high",
        massA
            && massA->properties.molarMassGPerMol.confidence == SaceConfidence::High, "");
    emit("generated_still_not_spawnable",
        massA && massB && !massA->spawnable && !massB->spawnable, "");

    ElementCount nitrogen[] = {{7, 1}};
    ChemicalIdentity nitrogenId = saceExactIdentity(ChemicalRepresentationKind::AtomicSpecies,
        "N", "elemental-nitrogen-synthetic", nitrogen, 1);
    char nitrogenSig[kSaceSignatureCap]{};
    bool nitrogenSigOk = writeChemicalSignature(nitrogenId, nitrogenSig, kSaceSignatureCap);
    SaceSubstanceRef nRef = cat.resolve(nitrogenId, true);
    SaceGeneratedRecord const *nRec = cat.record(nRef.generatedId);
    emit("unsupported_element_still_cataloged",
        nRef.kind == SaceResolutionKind::Generated
            && nRec && nRec->exactIdentity
            && nRec->kind == ChemicalRepresentationKind::AtomicSpecies, "");
    emit("unsupported_element_molar_mass_unknown",
        nRec && !nRec->properties.molarMassGPerMol.known
            && nRec->properties.molarMassGPerMol.source == SacePropertySource::Unknown
            && nRec->properties.molarMassGPerMol.confidence == SaceConfidence::Unknown, "");
    emit("unknown_molar_mass_does_not_change_signature",
        nRec && nitrogenSigOk
            && nRec->canonicalSignature == nitrogenSig
            && nRec->canonicalSignature.find("molar") == std::string::npos, "");

    cat.clear();
    emit("clear_returns_empty",
        cat.generatedCount() == 0
            && cat.record(aRef.generatedId) == nullptr, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
