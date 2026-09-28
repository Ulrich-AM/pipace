#include "substance/GeneratedMaterialRegistry.h"

#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceMolecule.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceTypes.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_map>

namespace {

struct RegistryState {
    std::unordered_map<uint32_t, GeneratedRuntimeMaterial> byHandle;
    std::unordered_map<SaceRecordId, uint32_t> byRecord;
    uint32_t nextHandle = 1;
};

RegistryState &state() {
    static RegistryState s;
    return s;
}

SaceMolecularGraph etheneGraph() {
    SaceMolecularGraph g{};
    g.atoms.push_back({kAtomicCarbon, 0});
    g.atoms.push_back({kAtomicCarbon, 0});
    g.bonds.push_back({0, 1, SaceBondOrder::Double});
    for (uint16_t c = 0; c < 2; ++c) {
        for (int i = 0; i < 2; ++i) {
            uint16_t h = static_cast<uint16_t>(g.atoms.size());
            g.atoms.push_back({kAtomicHydrogen, 0});
            g.bonds.push_back({c, h, SaceBondOrder::Single});
        }
    }
    return g;
}

bool near(double a, double b, double tol) {
    if (!std::isfinite(a) || !std::isfinite(b))
        return false;
    return std::fabs(a - b) <= tol;
}

} // namespace

void resetGeneratedRuntimeMaterials() {
    RegistryState &s = state();
    s.byHandle.clear();
    s.byRecord.clear();
}

std::size_t generatedRuntimeMaterialCount() {
    return state().byHandle.size();
}

bool registerGeneratedRuntimeMaterial(SaceRecordId recordId, RuntimeSubstanceRef &out) {
    out = runtimeNone();
    if (recordId == kSaceRecordNone)
        return false;

    RegistryState &s = state();
    auto existing = s.byRecord.find(recordId);
    if (existing != s.byRecord.end()) {
        auto it = s.byHandle.find(existing->second);
        if (it != s.byHandle.end() && it->second.valid) {
            out = it->second.runtimeRef;
            return true;
        }
    }

    SaceGeneratedRecord const *rec = saceGeneratedCatalog().record(recordId);
    if (!rec)
        return false;

    SaceCompiledSimulationProfile profile{};
    if (!saceCompileSimulationProfile(*rec, profile) || !profile.valid)
        return false;

    if (s.nextHandle == 0 || s.nextHandle > kRuntimeGeneratedHandleMax)
        return false;

    uint32_t handle = s.nextHandle++;
    RuntimeSubstanceRef ref = runtimeGeneratedHandle(handle);
    if (!runtimeSubstanceIsGenerated(ref))
        return false;

    GeneratedRuntimeMaterial entry{};
    entry.runtimeRef = ref;
    entry.sourceRecord = recordId;
    entry.profile = profile;
    entry.displayOrdinal = rec->displayOrdinal;
    entry.valid = true;
    s.byHandle[handle] = entry;
    s.byRecord[recordId] = handle;
    out = ref;
    return true;
}

RuntimeSubstanceRef runtimeGenerated(SaceRecordId recordId) {
    if (recordId == kSaceRecordNone)
        return runtimeNone();
    auto it = state().byRecord.find(recordId);
    if (it == state().byRecord.end())
        return runtimeNone();
    auto hit = state().byHandle.find(it->second);
    if (hit == state().byHandle.end() || !hit->second.valid)
        return runtimeNone();
    return hit->second.runtimeRef;
}

SaceRecordId runtimeGeneratedId(RuntimeSubstanceRef ref) {
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat)
        return kSaceRecordNone;
    return mat->sourceRecord;
}

RuntimeSubstanceRef runtimeRefFromSaceSubstanceRef(SaceSubstanceRef const &ref) {
    if (ref.kind == SaceResolutionKind::BuiltIn)
        return runtimeBuiltIn(ref.builtinId);
    if (ref.kind == SaceResolutionKind::Generated)
        return runtimeGenerated(ref.generatedId);
    return runtimeNone();
}

GeneratedRuntimeMaterial const *generatedRuntimeMaterial(RuntimeSubstanceRef ref) {
    if (!runtimeSubstanceIsGenerated(ref))
        return nullptr;
    auto it = state().byHandle.find(runtimeGeneratedHandleValue(ref));
    if (it == state().byHandle.end() || !it->second.valid)
        return nullptr;
    return &it->second;
}

SaceCompiledLiquidProfile const *runtimeLiquidProfile(RuntimeSubstanceRef ref) {
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat || !mat->profile.liquid.valid)
        return nullptr;
    return &mat->profile.liquid;
}

SaceCompiledGasProfile const *runtimeGasProfile(RuntimeSubstanceRef ref) {
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat || !mat->profile.gas.valid)
        return nullptr;
    return &mat->profile.gas;
}

SaceCompiledPhaseProfile const *runtimeLiquidGasPhaseProfile(RuntimeSubstanceRef ref) {
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat || !mat->profile.phase.liquidGasValid)
        return nullptr;
    return &mat->profile.phase;
}

bool runtimeSupportsPhase(RuntimeSubstanceRef ref, MatterPhase phase) {
    if (runtimeSubstanceIsBuiltIn(ref))
        return supportsPhase(runtimeBuiltinId(ref), phase);
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat)
        return false;
    switch (phase) {
        case MatterPhase::Liquid: return mat->profile.liquid.valid;
        case MatterPhase::Gas: return mat->profile.gas.valid;
        case MatterPhase::Solid: return false;
        case MatterPhase::Plasma: return false;
        case MatterPhase::None: return false;
    }
    return false;
}

bool runtimeFormatDisplayName(RuntimeSubstanceRef ref, char *out, int cap) {
    if (!out || cap <= 0)
        return false;
    out[0] = '\0';
    if (runtimeSubstanceIsBuiltIn(ref)) {
        char const *name = substanceDef(runtimeBuiltinId(ref)).displayName;
        if (!name)
            return false;
        std::snprintf(out, static_cast<std::size_t>(cap), "%s", name);
        return out[0] != '\0';
    }
    GeneratedRuntimeMaterial const *mat = generatedRuntimeMaterial(ref);
    if (!mat)
        return false;
    return formatGeneratedElementName(mat->displayOrdinal, out, cap);
}

void runSaceRuntimeRegistryDiagnostics() {
    std::ofstream out(miscFile("sace_runtime_registry_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    resetGeneratedRuntimeMaterials();
    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();

    RuntimeSubstanceRef water = runtimeBuiltIn(SUBSTANCE_WATER);
    RuntimeSubstanceRef air = runtimeBuiltIn(SUBSTANCE_AIR);
    RuntimeSubstanceRef oxygen = runtimeBuiltIn(SUBSTANCE_OXYGEN);
    RuntimeSubstanceRef co2 = runtimeBuiltIn(SUBSTANCE_CARBON_DIOXIDE);
    emit("built_in_water_roundtrip",
        runtimeSubstanceValid(water)
            && runtimeSubstanceIsBuiltIn(water)
            && runtimeBuiltinId(water) == SUBSTANCE_WATER
            && !runtimeSubstanceIsGenerated(water)
            && runtimeGeneratedId(water) == kSaceRecordNone, "");
    emit("built_in_air_roundtrip",
        runtimeSubstanceValid(air)
            && runtimeSubstanceIsBuiltIn(air)
            && runtimeBuiltinId(air) == SUBSTANCE_AIR
            && runtimeGeneratedId(air) == kSaceRecordNone, "");
    emit("built_in_oxygen_roundtrip",
        runtimeSubstanceValid(oxygen)
            && runtimeBuiltinId(oxygen) == SUBSTANCE_OXYGEN, "");
    emit("built_in_co2_roundtrip",
        runtimeSubstanceValid(co2)
            && runtimeBuiltinId(co2) == SUBSTANCE_CARBON_DIOXIDE, "");
    emit("built_in_wrong_generated_extract_rejected",
        runtimeGeneratedId(water) == kSaceRecordNone
            && generatedRuntimeMaterial(water) == nullptr
            && runtimeLiquidProfile(water) == nullptr, "");
    emit("built_in_water_phase_agrees",
        runtimeSupportsPhase(water, MatterPhase::Liquid) == supportsPhase(SUBSTANCE_WATER, MatterPhase::Liquid)
            && runtimeSupportsPhase(water, MatterPhase::Gas) == supportsPhase(SUBSTANCE_WATER, MatterPhase::Gas)
            && runtimeSupportsPhase(water, MatterPhase::Solid) == supportsPhase(SUBSTANCE_WATER, MatterPhase::Solid), "");
    emit("substance_none_is_runtime_none",
        runtimeSubstanceIsNone(runtimeBuiltIn(SUBSTANCE_NONE))
            && runtimeBuiltinId(runtimeBuiltIn(SUBSTANCE_NONE)) == SUBSTANCE_NONE, "");
    emit("invalid_builtin_id_rejected",
        runtimeSubstanceIsNone(runtimeBuiltIn(SUBSTANCE_COUNT))
            && runtimeSubstanceIsNone(runtimeBuiltIn(static_cast<SubstanceId>(0xFFFFu))), "");

    SaceMolecularGraph gA{}, gB{};
    saceSyntheticC2H6OGraphA(gA);
    saceSyntheticC2H6OGraphB(gB);
    ElementCount c2h6o[] = {
        {kAtomicCarbon, 2}, {kAtomicHydrogen, 6}, {kAtomicOxygen, 1}
    };
    ChemicalIdentity idA = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-a", c2h6o, 3);
    ChemicalIdentity idB = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-b", c2h6o, 3);
    SaceSubstanceRef saceA = cat.resolve(idA, true);
    SaceSubstanceRef saceB = cat.resolve(idB, true);
    cat.attachMolecularGraph(saceA.generatedId, "synthetic-structure-a", gA);
    cat.attachMolecularGraph(saceB.generatedId, "synthetic-structure-b", gB);
    SaceGeneratedRecord const *recA = cat.record(saceA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(saceB.generatedId);

    char sigBefore[kSaceSignatureCap]{};
    ChemicalIdentity viewBefore = recA ? cat.identityView(*recA) : ChemicalIdentity{};
    writeChemicalSignature(viewBefore, sigBefore, kSaceSignatureCap);
    uint32_t ordinalBefore = recA ? recA->displayOrdinal : 0;
    SaceRecordId recordBefore = recA ? recA->recordId : kSaceRecordNone;
    std::string structureBefore = recA ? recA->structureKey : "";

    RuntimeSubstanceRef genA{};
    RuntimeSubstanceRef genB{};
    bool regA = recA && registerGeneratedRuntimeMaterial(recA->recordId, genA);
    bool regB = recB && registerGeneratedRuntimeMaterial(recB->recordId, genB);
    emit("generated_a_registered",
        regA && runtimeSubstanceValid(genA) && runtimeSubstanceIsGenerated(genA), "");
    emit("generated_b_registered",
        regB && runtimeSubstanceValid(genB) && runtimeSubstanceIsGenerated(genB), "");
    emit("generated_a_roundtrip",
        recA && runtimeGeneratedId(genA) == recA->recordId
            && runtimeGenerated(recA->recordId) == genA
            && runtimeBuiltinId(genA) == SUBSTANCE_NONE, "");
    emit("generated_b_roundtrip",
        recB && runtimeGeneratedId(genB) == recB->recordId
            && runtimeGenerated(recB->recordId) == genB
            && runtimeBuiltinId(genB) == SUBSTANCE_NONE, "");
    emit("generated_a_b_distinct",
        regA && regB && genA != genB
            && runtimeGeneratedId(genA) != runtimeGeneratedId(genB), "");

    std::size_t countAfterAB = generatedRuntimeMaterialCount();
    RuntimeSubstanceRef genA2{};
    bool regA2 = recA && registerGeneratedRuntimeMaterial(recA->recordId, genA2);
    emit("registration_idempotent",
        regA2 && genA2 == genA
            && generatedRuntimeMaterialCount() == countAfterAB
            && generatedRuntimeMaterial(genA) == generatedRuntimeMaterial(genA2), "");

    SaceCompiledSimulationProfile directA{};
    bool compiledA = recA && saceCompileSimulationProfile(*recA, directA);
    GeneratedRuntimeMaterial const *matA = generatedRuntimeMaterial(genA);
    emit("generated_a_liquid_profile_retained",
        compiledA && matA && matA->profile.liquid.valid
            && matA->sourceRecord == recA->recordId
            && near(matA->profile.liquid.densityRelativeToWater,
                directA.liquid.densityRelativeToWater, 1.0e-6f)
            && near(matA->profile.liquid.thermalConductivityWPerMK,
                directA.liquid.thermalConductivityWPerMK, 1.0e-6f), "");
    emit("generated_a_gas_profile_retained",
        compiledA && matA && matA->profile.gas.valid
            && near(matA->profile.gas.thermalConductivityWPerMK,
                directA.gas.thermalConductivityWPerMK, 1.0e-6f)
            && near(matA->profile.gas.specificHeatJPerKgK,
                directA.gas.specificHeatJPerKgK, 1.0e-3f), "");
    emit("generated_a_phase_profile_retained",
        compiledA && matA && matA->profile.phase.liquidGasValid
            && near(matA->profile.phase.normalBoilingPointK,
                directA.phase.normalBoilingPointK, 1.0e-4f), "");

    emit("generated_liquid_supported",
        runtimeSupportsPhase(genA, MatterPhase::Liquid)
            && runtimeSupportsPhase(genB, MatterPhase::Liquid), "");
    emit("generated_gas_supported",
        runtimeSupportsPhase(genA, MatterPhase::Gas)
            && runtimeSupportsPhase(genB, MatterPhase::Gas), "");
    emit("generated_solid_not_supported",
        !runtimeSupportsPhase(genA, MatterPhase::Solid)
            && !runtimeSupportsPhase(genB, MatterPhase::Solid)
            && !runtimeSupportsPhase(genA, MatterPhase::Plasma), "");

    std::size_t countBeforeBad = generatedRuntimeMaterialCount();
    ElementCount c2h4[] = {{kAtomicCarbon, 2}, {kAtomicHydrogen, 4}};
    ChemicalIdentity idE = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H4", "synthetic-ethene", c2h4, 2);
    SaceSubstanceRef saceE = cat.resolve(idE, true);
    cat.attachMolecularGraph(saceE.generatedId, "synthetic-ethene", etheneGraph());
    RuntimeSubstanceRef genE{};
    bool regE = registerGeneratedRuntimeMaterial(saceE.generatedId, genE);
    emit("unsupported_registration_rejected",
        !regE && runtimeSubstanceIsNone(genE)
            && generatedRuntimeMaterialCount() == countBeforeBad
            && generatedRuntimeMaterial(genE) == nullptr, "");

    char displayA[32]{};
    char catalogName[32]{};
    bool namedRuntime = runtimeFormatDisplayName(genA, displayA, 32);
    bool namedCatalog = recA && formatGeneratedElementName(recA->displayOrdinal, catalogName, 32);
    emit("generated_display_matches_catalog",
        namedRuntime && namedCatalog && std::strcmp(displayA, catalogName) == 0
            && std::strncmp(displayA, "Element #", 9) == 0, displayA);
    char waterName[32]{};
    emit("built_in_display_unchanged",
        runtimeFormatDisplayName(water, waterName, 32)
            && std::strcmp(waterName, substanceDef(SUBSTANCE_WATER).displayName) == 0, waterName);

    recA = cat.record(saceA.generatedId);
    char sigAfter[kSaceSignatureCap]{};
    ChemicalIdentity viewAfter = recA ? cat.identityView(*recA) : ChemicalIdentity{};
    writeChemicalSignature(viewAfter, sigAfter, kSaceSignatureCap);
    emit("catalog_identity_unchanged",
        recA && std::strcmp(sigBefore, sigAfter) == 0
            && recA->recordId == recordBefore
            && recA->displayOrdinal == ordinalBefore
            && recA->structureKey == structureBefore, sigAfter);

    RuntimeSubstanceRef stale = genA;
    SaceRecordId oldRecord = runtimeGeneratedId(stale);
    std::size_t countBeforeReset = generatedRuntimeMaterialCount();
    cat.clear();
    emit("registry_reset_invalidates_entry",
        countBeforeReset > 0
            && generatedRuntimeMaterialCount() == 0
            && generatedRuntimeMaterial(stale) == nullptr
            && runtimeGeneratedId(stale) == kSaceRecordNone
            && runtimeGenerated(oldRecord) == runtimeNone(), "");

    SaceSubstanceRef saceA2 = cat.resolve(idA, true);
    cat.attachMolecularGraph(saceA2.generatedId, "synthetic-structure-a", gA);
    RuntimeSubstanceRef genANew{};
    bool regANew = registerGeneratedRuntimeMaterial(saceA2.generatedId, genANew);
    emit("stale_ref_does_not_resolve_after_reset",
        regANew && genANew != stale
            && generatedRuntimeMaterial(stale) == nullptr
            && generatedRuntimeMaterial(genANew) != nullptr
            && runtimeGeneratedId(genANew) == saceA2.generatedId
            && runtimeGeneratedId(stale) != saceA2.generatedId, "");

    recA = cat.record(saceA2.generatedId);
    emit("generated_record_still_unspawnable",
        recA && !recA->spawnable, "");
    emit("substance_count_remains_12", SUBSTANCE_COUNT == 12, std::to_string(SUBSTANCE_COUNT));
    emit("sace_substance_ref_bridge_builtin",
        runtimeRefFromSaceSubstanceRef(SaceSubstanceRef{SaceResolutionKind::BuiltIn, SUBSTANCE_WATER, kSaceRecordNone})
            == water, "");
    emit("runtime_ref_size_is_4",
        sizeof(RuntimeSubstanceRef) == 4, std::to_string(sizeof(RuntimeSubstanceRef)));
    emit("built_ins_not_in_generated_registry",
        generatedRuntimeMaterial(water) == nullptr
            && generatedRuntimeMaterialCount() == 1, "");

    resetGeneratedRuntimeMaterials();
    cat.clear();
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
