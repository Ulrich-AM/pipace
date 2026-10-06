#include "substance/RuntimeComponentDiagnostics.h"

#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceMolecule.h"
#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "substance/GeneratedMaterialRegistry.h"
#include "substance/LiquidMixtureProperties.h"
#include "substance/PhaseTransfer.h"

#include <cmath>
#include <fstream>
#include <string>

namespace {

bool near(double a, double b, double tol) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

} // namespace

void runRuntimeComponentDiagnostics() {
    std::ofstream out(miscFile("runtime_component_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail = "") {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    emit("runtime_ref_is_4_bytes", sizeof(RuntimeSubstanceRef) == 4,
        std::to_string(sizeof(RuntimeSubstanceRef)));
    emit("liquid_component_compact", sizeof(LiquidComponent) == 8,
        std::to_string(sizeof(LiquidComponent)));
    emit("gas_component_compact", sizeof(GasComponent) == 8,
        std::to_string(sizeof(GasComponent)));

    LiquidComponent liquids[kMaxLiquidComponents]{};
    int liquidCount = 0;
    emit("builtin_liquid_add_water",
        addLiquidPayload(liquids, liquidCount, SUBSTANCE_WATER, 0.6f) == 0.0f
            && liquidCount == 1
            && liquids[0].id == SUBSTANCE_WATER
            && near(liquids[0].amount, 0.6, 1.0e-6), "");
    emit("builtin_liquid_merge_water",
        addLiquidPayload(liquids, liquidCount, SUBSTANCE_WATER, 0.2f) == 0.0f
            && liquidCount == 1
            && near(liquids[0].amount, 0.8, 1.0e-6), "");
    emit("builtin_liquid_add_honey",
        addLiquidPayload(liquids, liquidCount, SUBSTANCE_HONEY, 0.2f) == 0.0f
            && liquidCount == 2
            && liquidPayloadAmount(liquids, liquidCount, SUBSTANCE_HONEY) > 0.19f, "");
    emit("air_rejected_as_liquid",
        addLiquidPayload(liquids, liquidCount, SUBSTANCE_AIR, 0.1f) == 0.0f
            && findLiquidComponent(liquids, liquidCount, SUBSTANCE_AIR) < 0, "");

    GasComponent gases[kMaxGasComponents]{};
    int gasCount = 0;
    emit("builtin_gas_add_air",
        addGasPayload(gases, gasCount, SUBSTANCE_AIR, 0.7f) == 0.0f
            && gasCount == 1 && gases[0].id == SUBSTANCE_AIR, "");
    emit("builtin_gas_add_oxygen",
        addGasPayload(gases, gasCount, SUBSTANCE_OXYGEN, 0.3f) == 0.0f
            && gasCount == 2
            && gasPayloadAmount(gases, gasCount, SUBSTANCE_OXYGEN) > 0.29f, "");
    emit("honey_rejected_as_gas",
        addGasPayload(gases, gasCount, SUBSTANCE_HONEY, 0.1f) == 0.0f
            && findGasComponent(gases, gasCount, SUBSTANCE_HONEY) < 0, "");

    resetGeneratedRuntimeMaterials();
    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();

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
    bool attachA = cat.attachMolecularGraph(saceA.generatedId, "synthetic-structure-a", gA);
    bool attachB = cat.attachMolecularGraph(saceB.generatedId, "synthetic-structure-b", gB);
    SaceGeneratedRecord const *recA = cat.record(saceA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(saceB.generatedId);

    RuntimeSubstanceRef genA{}, genB{};
    bool regA = attachA && recA && registerGeneratedRuntimeMaterial(recA->recordId, genA);
    bool regB = attachB && recB && registerGeneratedRuntimeMaterial(recB->recordId, genB);
    emit("generated_refs_registered",
        regA && regB && runtimeSubstanceIsGenerated(genA)
            && runtimeSubstanceIsGenerated(genB) && genA != genB, "");
    emit("generated_refs_phase_valid",
        validLiquidComponentId(genA) && validGasComponentId(genA)
            && validLiquidComponentId(genB) && validGasComponentId(genB), "");

    LiquidComponent generatedLiquid[kMaxLiquidComponents]{};
    int generatedLiquidCount = 0;
    float lUnplaced1 = addLiquidPayload(generatedLiquid, generatedLiquidCount, genA, 0.65f);
    float lUnplaced2 = addLiquidPayload(generatedLiquid, generatedLiquidCount, genA, 0.15f);
    float lUnplaced3 = addLiquidPayload(generatedLiquid, generatedLiquidCount, genB, 0.20f);
    emit("generated_liquid_payload_add_merge",
        lUnplaced1 == 0.0f && lUnplaced2 == 0.0f && lUnplaced3 == 0.0f
            && generatedLiquidCount == 2
            && near(liquidPayloadAmount(generatedLiquid, generatedLiquidCount, genA), 0.80, 1.0e-6)
            && near(liquidPayloadAmount(generatedLiquid, generatedLiquidCount, genB), 0.20, 1.0e-6), "");

    LiquidComponent liquidCopy[kMaxLiquidComponents]{};
    int liquidCopyCount = 0;
    copyLiquidPayload(liquidCopy, liquidCopyCount, generatedLiquid, generatedLiquidCount);
    scaleLiquidPayload(liquidCopy, liquidCopyCount, 0.5f);
    emit("generated_liquid_payload_copy_scale",
        liquidCopyCount == 2
            && liquidCopy[0].id == generatedLiquid[0].id
            && liquidCopy[1].id == generatedLiquid[1].id
            && near(liquidPayloadSum(liquidCopy, liquidCopyCount), 0.50, 1.0e-6), "");

    LiquidComponentView liquidView{};
    copyLiquidPayload(liquidView.items, liquidView.count, generatedLiquid, generatedLiquidCount);
    LiquidMixtureProperties liquidMix = evaluateLiquidMixture(liquidView);
    float liquidMu = evaluateLiquidMixtureViscosity(
        liquidView, kSaceCompilerReferenceTemperatureK);
    emit("generated_liquid_mixture_properties",
        std::isfinite(liquidMix.density) && liquidMix.density > 0.0f
            && std::isfinite(liquidMix.specificHeat) && liquidMix.specificHeat > 0.0f
            && std::isfinite(liquidMix.conductivity) && liquidMix.conductivity > 0.0f
            && std::isfinite(liquidMix.surfaceTension) && liquidMix.surfaceTension > 0.0f
            && std::isfinite(liquidMu) && liquidMu > 0.0f, "");

    GeneratedRuntimeMaterial const *matA = generatedRuntimeMaterial(genA);
    LiquidComponentView pureA{};
    pureA.items[pureA.count++] = {genA, 1.0f};
    LiquidMixtureProperties pureAMix = evaluateLiquidMixture(pureA);
    emit("pure_generated_liquid_matches_compiler",
        matA && near(pureAMix.density, matA->profile.liquid.densityRelativeToWater, 1.0e-6)
            && near(pureAMix.specificHeat, matA->profile.liquid.specificHeatJPerKgK, 1.0e-3)
            && near(pureAMix.conductivity, matA->profile.liquid.thermalConductivityWPerMK, 1.0e-6)
            && near(pureAMix.surfaceTension, matA->profile.liquid.solverSurfaceTension, 1.0e-6), "");

    GasComponent generatedGas[kMaxGasComponents]{};
    int generatedGasCount = 0;
    float gUnplaced1 = addGasPayload(generatedGas, generatedGasCount, genA, 0.55f);
    float gUnplaced2 = addGasPayload(generatedGas, generatedGasCount, genB, 0.45f);
    emit("generated_gas_payload_add",
        gUnplaced1 == 0.0f && gUnplaced2 == 0.0f
            && generatedGasCount == 2
            && near(gasPayloadSum(generatedGas, generatedGasCount), 1.0, 1.0e-6), "");

    GasComponent gasCopy[kMaxGasComponents]{};
    int gasCopyCount = 0;
    copyGasPayload(gasCopy, gasCopyCount, generatedGas, generatedGasCount);
    scaleGasPayload(gasCopy, gasCopyCount, 0.25f);
    emit("generated_gas_payload_copy_scale",
        gasCopyCount == 2
            && gasCopy[0].id == generatedGas[0].id
            && gasCopy[1].id == generatedGas[1].id
            && near(gasPayloadSum(gasCopy, gasCopyCount), 0.25, 1.0e-6), "");

    GasComponentView gasView{};
    copyGasPayload(gasView.items, gasView.count, generatedGas, generatedGasCount);
    float gasCap = gasMixtureThermalCapacity(gasView, 4.0f);
    float gasRho = gasMixtureReferenceDensityKgM3(gasView);
    float gasK = gasMixtureConductivity(gasView);
    emit("generated_gas_mixture_properties",
        std::isfinite(gasCap) && gasCap > 0.0f
            && std::isfinite(gasRho) && gasRho > 0.0f
            && std::isfinite(gasK) && gasK > 0.0f, "");

    FluidEngine fluid;
    GasEngine gas;
    int fi = FluidEngine::ci(GW / 2, GH / 2);
    LiquidComponentView generatedCommit{};
    generatedCommit.items[generatedCommit.count++] = {genA, 0.5f};
    bool fluidCommit = fluid.tryCommitLiquidOccupancy(fi, generatedCommit);
    emit("fluid_engine_generated_commit_deferred",
        !fluidCommit
            && fluid.fill[static_cast<size_t>(fi)] <= kMinLiquidComponent
            && fluid.liquidCompCount[static_cast<size_t>(fi)] == 0, "");

    int gi = GasEngine::ci(GW / 2, GH / 2);
    GasComponentView generatedGasCommit{};
    generatedGasCommit.items[generatedGasCommit.count++] = {genA, 0.5f};
    bool gasCommit = gas.tryCommitGasOccupancy(gi, generatedGasCommit);
    emit("gas_engine_generated_commit_deferred",
        !gasCommit
            && gas.amount[static_cast<size_t>(gi)] <= GAS_MIN_AMOUNT
            && gas.gasCompCount[static_cast<size_t>(gi)] == 0, "");

    emit("generated_records_still_unspawnable",
        recA && recB && !recA->spawnable && !recB->spawnable, "");
    emit("substance_count_remains_12", SUBSTANCE_COUNT == 12,
        std::to_string(SUBSTANCE_COUNT));

    RuntimeSubstanceRef staleA = genA;
    cat.clear();
    emit("stale_generated_liquid_ref_invalid",
        !validLiquidComponentId(staleA), "");
    emit("stale_generated_gas_ref_invalid",
        !validGasComponentId(staleA), "");

    compactLiquidPayload(generatedLiquid, generatedLiquidCount);
    compactGasPayload(generatedGas, generatedGasCount);
    emit("stale_refs_removed_by_compaction",
        generatedLiquidCount == 0 && generatedGasCount == 0, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
