#include "substance/RuntimeComponentDiagnostics.h"

#include "chemistry/SaceCatalog.h"
#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceMolecule.h"
#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "substance/GeneratedMaterialRegistry.h"
#include "substance/LiquidMixtureProperties.h"
#include "substance/PhaseTransfer.h"
#include "substance/RuntimeSubstanceProperties.h"

#include <cmath>
#include <fstream>
#include <string>

namespace {

bool near(double a, double b, double tol = 1.0e-6) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

} // namespace

void runRuntimeComponentDiagnostics() {
    std::ofstream out(miscFile("runtime_component_diag.tsv"));
    int passed = 0;
    int failed = 0;
    auto emit = [&](char const *name, bool ok, std::string detail = "") {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    emit("liquid_component_stays_compact",
        sizeof(LiquidComponent) == sizeof(RuntimeSubstanceRef) + sizeof(float),
        std::to_string(sizeof(LiquidComponent)));
    emit("gas_component_stays_compact",
        sizeof(GasComponent) == sizeof(RuntimeSubstanceRef) + sizeof(float),
        std::to_string(sizeof(GasComponent)));

    FluidEngine fluid;
    GasEngine gas;
    fluid.clearWorld();

    int li = FluidEngine::ci(40, 40);
    LiquidComponentView builtL{};
    builtL.items[builtL.count++] = {SUBSTANCE_WATER, 0.60f};
    builtL.items[builtL.count++] = {SUBSTANCE_HONEY, 0.40f};
    bool builtLCommit = fluid.tryCommitLiquidOccupancy(li, builtL);
    LiquidComponentView builtLRead = fluid.liquidComponents(li);
    emit("builtin_liquid_commit_roundtrip",
        builtLCommit && fluid.liquidCompositionValid(li)
            && builtLRead.count == 2
            && near(fluid.liquidComponentAmount(li, SUBSTANCE_WATER), 0.60)
            && near(fluid.liquidComponentAmount(li, SUBSTANCE_HONEY), 0.40)
            && fluid.dominantLiquidRef(li) == runtimeBuiltIn(SUBSTANCE_WATER),
        "count=" + std::to_string(builtLRead.count));

    int gi = GasEngine::ci(42, 40);
    GasComponentView builtG{};
    builtG.items[builtG.count++] = {SUBSTANCE_AIR, 0.70f};
    builtG.items[builtG.count++] = {SUBSTANCE_OXYGEN, 0.30f};
    bool builtGCommit = gas.tryCommitGasOccupancy(gi, builtG);
    GasComponentView builtGRead = gas.gasComponents(gi);
    emit("builtin_gas_commit_roundtrip",
        builtGCommit && gas.gasCompositionValid(gi)
            && builtGRead.count == 2
            && near(gas.gasComponentAmount(gi, SUBSTANCE_AIR), 0.70)
            && near(gas.gasComponentAmount(gi, SUBSTANCE_OXYGEN), 0.30)
            && gas.dominantGasRef(gi) == runtimeBuiltIn(SUBSTANCE_AIR),
        "count=" + std::to_string(builtGRead.count));

    resetGeneratedRuntimeMaterials();
    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();

    SaceMolecularGraph graph{};
    saceSyntheticC2H6OGraphA(graph);
    ElementCount c2h6o[] = {
        {kAtomicCarbon, 2}, {kAtomicHydrogen, 6}, {kAtomicOxygen, 1}
    };
    ChemicalIdentity identity = saceExactIdentity(
        ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "runtime-component-a", c2h6o, 3);
    SaceSubstanceRef sref = cat.resolve(identity, true);
    bool attached = sref.kind == SaceResolutionKind::Generated
        && cat.attachMolecularGraph(sref.generatedId, "runtime-component-a", graph);

    RuntimeSubstanceRef generated{};
    bool registered = attached
        && registerGeneratedRuntimeMaterial(sref.generatedId, generated);
    emit("generated_runtime_registered",
        registered && runtimeSubstanceIsGenerated(generated)
            && runtimeSupportsPhase(generated, MatterPhase::Liquid)
            && runtimeSupportsPhase(generated, MatterPhase::Gas));

    LiquidComponentView genL{};
    genL.items[genL.count++] = {generated, 0.65f};
    genL.items[genL.count++] = {SUBSTANCE_WATER, 0.35f};
    bool genLCommit = registered && fluid.tryCommitLiquidOccupancy(li, genL);
    emit("generated_liquid_storage_valid",
        genLCommit && fluid.liquidCompositionValid(li)
            && near(fluid.liquidComponentAmount(li, generated), 0.65)
            && near(fluid.liquidComponentAmount(li, SUBSTANCE_WATER), 0.35)
            && fluid.dominantLiquidRef(li) == generated,
        "raw=" + std::to_string(generated.raw));

    emit("builtin_dominant_adapter_does_not_alias_generated",
        genLCommit && fluid.dominantLiquidSubstance(li) == SUBSTANCE_NONE);

    RuntimeLiquidPropertySample genLP{};
    LiquidMixtureProperties mixed = fluid.mixProperties(li);
    bool genLPok = registered
        && sampleRuntimeLiquidProperties(generated,
            kSaceCompilerReferenceTemperatureK, genLP);
    emit("generated_liquid_mix_uses_runtime_properties",
        genLCommit && genLPok
            && std::isfinite(mixed.density) && mixed.density > 0.0f
            && std::isfinite(mixed.specificHeat) && mixed.specificHeat > 0.0f
            && std::isfinite(mixed.conductivity) && mixed.conductivity > 0.0f
            && mixed.density < 1.0f
            && mixed.density > genLP.densityRelativeToWater,
        "mix_rho=" + std::to_string(mixed.density)
            + " gen_rho=" + std::to_string(genLP.densityRelativeToWater));

    float beforeGen = fluid.liquidComponentAmount(li, generated);
    float beforeWater = fluid.liquidComponentAmount(li, SUBSTANCE_WATER);
    LiquidCarry carry = fluid.takeLiquidCarry(li, 0.25f);
    float carryGen = liquidPayloadAmount(carry.comps, carry.compCount, generated);
    float carryWater = liquidPayloadAmount(carry.comps, carry.compCount, SUBSTANCE_WATER);
    emit("generated_liquid_extraction_preserves_identity",
        carryGen > kMinLiquidComponent
            && carryWater > kMinLiquidComponent
            && near(beforeGen,
                fluid.liquidComponentAmount(li, generated) + carryGen, 1.0e-5)
            && near(beforeWater,
                fluid.liquidComponentAmount(li, SUBSTANCE_WATER) + carryWater, 1.0e-5),
        "carry_gen=" + std::to_string(carryGen)
            + " carry_water=" + std::to_string(carryWater));

    LiquidComponent copied[kMaxLiquidComponents]{};
    int copiedCount = 0;
    copyLiquidPayload(copied, copiedCount, carry.comps, carry.compCount);
    scaleLiquidPayload(copied, copiedCount, 0.5f);
    compactLiquidPayload(copied, copiedCount);
    emit("generated_liquid_payload_copy_preserves_ref",
        findLiquidComponent(copied, copiedCount, generated) >= 0
            && liquidPayloadAmount(copied, copiedCount, generated) > 0.0f);

    GasComponentView genG{};
    genG.items[genG.count++] = {generated, 0.75f};
    genG.items[genG.count++] = {SUBSTANCE_AIR, 0.25f};
    bool genGCommit = registered && gas.tryCommitGasOccupancy(gi, genG);
    emit("generated_gas_storage_valid",
        genGCommit && gas.gasCompositionValid(gi)
            && near(gas.gasComponentAmount(gi, generated), 0.75)
            && near(gas.gasComponentAmount(gi, SUBSTANCE_AIR), 0.25)
            && gas.dominantGasRef(gi) == generated,
        "raw=" + std::to_string(generated.raw));

    emit("builtin_gas_dominant_adapter_does_not_alias_generated",
        genGCommit && gas.dominantGasSubstance(gi) == SUBSTANCE_NONE);

    GasComponentView genGRead = gas.gasComponents(gi);
    float gasCap = gasMixtureThermalCapacity(genGRead, fluid.config.cellsPerMeter);
    float gasRho = gasMixtureReferenceDensityKgM3(genGRead);
    float gasK = gasMixtureConductivity(genGRead);
    emit("generated_gas_mixture_properties_finite",
        genGCommit
            && std::isfinite(gasCap) && gasCap > 0.0f
            && std::isfinite(gasRho) && gasRho > 0.0f
            && std::isfinite(gasK) && gasK > 0.0f,
        "cap=" + std::to_string(gasCap)
            + " rho=" + std::to_string(gasRho)
            + " k=" + std::to_string(gasK));

    GasComponent gasCopy[kMaxGasComponents]{};
    int gasCopyCount = 0;
    copyGasPayload(gasCopy, gasCopyCount, genGRead.items, genGRead.count);
    scaleGasPayload(gasCopy, gasCopyCount, 0.5f);
    compactGasPayload(gasCopy, gasCopyCount);
    emit("generated_gas_payload_copy_preserves_ref",
        findGasComponent(gasCopy, gasCopyCount, generated) >= 0
            && gasPayloadAmount(gasCopy, gasCopyCount, generated) > 0.0f);

    SaceGeneratedRecord const *rec = cat.record(sref.generatedId);
    emit("generated_record_still_unspawnable",
        rec && !rec->spawnable);
    emit("substance_count_remains_12",
        SUBSTANCE_COUNT == 12, std::to_string(SUBSTANCE_COUNT));

    RuntimeSubstanceRef stale = generated;
    cat.clear(); // also resets generated runtime registry
    emit("stale_ref_fails_component_validation",
        !validLiquidComponentId(stale)
            && !validGasComponentId(stale)
            && !fluid.liquidCompositionValid(li)
            && !gas.gasCompositionValid(gi));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
