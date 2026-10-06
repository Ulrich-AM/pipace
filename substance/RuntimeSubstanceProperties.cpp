#include "substance/RuntimeSubstanceProperties.h"

#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceMolecule.h"
#include "fluid/DiagOutput.h"
#include "substance/GeneratedMaterialRegistry.h"
#include "substance/SubstanceTypes.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <string>

namespace {

bool finitePositive(double v) {
    return std::isfinite(v) && v > 0.0;
}

bool validTemperature(float temperatureK) {
    return std::isfinite(temperatureK)
        && temperatureK > MIN_SAFE_TEMPERATURE_K
        && temperatureK <= MAX_SAFE_TEMPERATURE_K;
}

bool sampleBuiltInLiquid(
    SubstanceId id,
    float temperatureK,
    RuntimeLiquidPropertySample &out)
{
    if (!validSubstance(id) || id == SUBSTANCE_NONE)
        return false;
    if (!supportsPhase(id, MatterPhase::Liquid) || !hasFluidProperties(id))
        return false;

    FluidProperties const &fluid = fluidForSubstance(id);
    ThermalProperties const &thermal = thermalForSubstance(id);
    if (!fluid.valid || !thermal.valid)
        return false;

    float rho = fluid.densityAtTemperature(temperatureK);
    float cp = thermal.specificHeat;
    float k = thermal.conductivity;
    float mu = fluid.viscosityAtTemperature(temperatureK);
    float sigma = fluid.surfaceTension;
    if (!finitePositive(rho) || !finitePositive(cp) || !finitePositive(k)
        || !finitePositive(mu) || !finitePositive(sigma))
        return false;

    out.valid = true;
    out.temperatureK = temperatureK;
    out.densityRelativeToWater = rho;
    out.specificHeatJPerKgK = cp;
    out.thermalConductivityWPerMK = k;
    out.solverViscosity = mu;
    out.solverSurfaceTension = sigma;
    return true;
}

bool sampleGeneratedLiquid(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeLiquidPropertySample &out)
{
    SaceCompiledLiquidProfile const *liquid = runtimeLiquidProfile(ref);
    if (!liquid || !liquid->valid)
        return false;

    double rhoKgM3 = 0.0;
    if (!saceCostaldSaturatedLiquidDensityKgPerM3(
            liquid->densityModel, temperatureK, rhoKgM3))
        return false;
    float rhoRel = 0.0f;
    if (!sacePhysicalDensityToWaterRelative(rhoKgM3, rhoRel))
        return false;

    double cpMolar = 0.0;
    if (!saceRowlinsonPolingLiquidHeatCapacityJPerMolK(
            liquid->heatCapacityModel, temperatureK, cpMolar))
        return false;
    double cpSpecific = 0.0;
    if (!saceMolarHeatCapacityToSpecificJPerKgK(
            cpMolar, liquid->heatCapacityModel.molarMassKgPerMol, cpSpecific))
        return false;

    double conductivity = 0.0;
    if (!saceSatoRiedelLiquidThermalConductivityWPerMK(
            liquid->conductivityModel, temperatureK, conductivity))
        return false;

    float viscosity = 0.0f;
    if (!saceCompiledSolverViscosityAtTemperature(
            *liquid, temperatureK, viscosity))
        return false;

    double surfacePhysical = 0.0;
    if (!saceSastriRaoSurfaceTensionNPerM(
            liquid->surfaceTensionModel, temperatureK, surfacePhysical))
        return false;
    float surfaceSolver = 0.0f;
    if (!saceCalibrateSolverSurfaceTensionFromPhysicalNPerM(
            surfacePhysical, surfaceSolver))
        return false;

    if (!finitePositive(rhoRel) || !finitePositive(cpSpecific)
        || !finitePositive(conductivity) || !finitePositive(viscosity)
        || !finitePositive(surfaceSolver))
        return false;

    out.valid = true;
    out.temperatureK = temperatureK;
    out.densityRelativeToWater = rhoRel;
    out.specificHeatJPerKgK = static_cast<float>(cpSpecific);
    out.thermalConductivityWPerMK = static_cast<float>(conductivity);
    out.solverViscosity = viscosity;
    out.solverSurfaceTension = surfaceSolver;
    return true;
}

bool sampleBuiltInGas(
    SubstanceId id,
    float temperatureK,
    RuntimeGasPropertySample &out)
{
    if (!validSubstance(id) || id == SUBSTANCE_NONE)
        return false;
    if (!supportsPhase(id, MatterPhase::Gas))
        return false;

    ThermalProperties const &thermal = thermalForSubstance(id);
    if (!thermal.valid)
        return false;

    float molarMass = 0.0f;
    if (id == SUBSTANCE_AIR) {
        molarMass = AIR_MOLAR_MASS_G_MOL;
    } else {
        ChemicalProperties const &chemical = chemicalForSubstance(id);
        if (chemical.valid)
            molarMass = chemical.molarMass;
    }

    float cp = gasPhaseSpecificHeat(id);
    float k = thermal.conductivity;
    if (!finitePositive(molarMass) || !finitePositive(cp) || !finitePositive(k))
        return false;

    out.valid = true;
    out.temperatureK = temperatureK;
    out.molarMassGPerMol = molarMass;
    out.specificHeatJPerKgK = cp;
    out.thermalConductivityWPerMK = k;
    return true;
}

bool sampleGeneratedGas(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeGasPropertySample &out)
{
    SaceCompiledGasProfile const *gas = runtimeGasProfile(ref);
    if (!gas || !gas->valid)
        return false;

    double cpMolar = 0.0;
    if (!saceJobackIdealGasHeatCapacityJPerMolK(
            gas->heatCapacityModel, temperatureK, cpMolar))
        return false;

    double molarMassKg = static_cast<double>(gas->molarMassGPerMol) * 1.0e-3;
    double cpSpecific = 0.0;
    if (!saceMolarHeatCapacityToSpecificJPerKgK(
            cpMolar, molarMassKg, cpSpecific))
        return false;

    double conductivity = 0.0;
    if (!saceGharagheiziGasThermalConductivityWPerMK(
            gas->conductivityModel, temperatureK, conductivity))
        return false;

    if (!finitePositive(gas->molarMassGPerMol)
        || !finitePositive(cpSpecific)
        || !finitePositive(conductivity))
        return false;

    out.valid = true;
    out.temperatureK = temperatureK;
    out.molarMassGPerMol = gas->molarMassGPerMol;
    out.specificHeatJPerKgK = static_cast<float>(cpSpecific);
    out.thermalConductivityWPerMK = static_cast<float>(conductivity);
    return true;
}

bool near(double a, double b, double tol) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

} // namespace

bool sampleRuntimeLiquidProperties(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeLiquidPropertySample &out)
{
    out = {};
    if (!validTemperature(temperatureK))
        return false;

    if (runtimeSubstanceIsBuiltIn(ref))
        return sampleBuiltInLiquid(runtimeBuiltinId(ref), temperatureK, out);
    if (runtimeSubstanceIsGenerated(ref))
        return sampleGeneratedLiquid(ref, temperatureK, out);
    return false;
}

bool sampleRuntimeGasProperties(
    RuntimeSubstanceRef ref,
    float temperatureK,
    RuntimeGasPropertySample &out)
{
    out = {};
    if (!validTemperature(temperatureK))
        return false;

    if (runtimeSubstanceIsBuiltIn(ref))
        return sampleBuiltInGas(runtimeBuiltinId(ref), temperatureK, out);
    if (runtimeSubstanceIsGenerated(ref))
        return sampleGeneratedGas(ref, temperatureK, out);
    return false;
}

void runRuntimeSubstancePropertyDiagnostics() {
    std::ofstream out(miscFile("runtime_property_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail = "") {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    constexpr float tref = kSaceCompilerReferenceTemperatureK;

    RuntimeLiquidPropertySample water{};
    bool waterOk = sampleRuntimeLiquidProperties(
        runtimeBuiltIn(SUBSTANCE_WATER), tref, water);
    emit("builtin_water_liquid_samples", waterOk && water.valid);
    emit("builtin_water_density_matches",
        waterOk && near(water.densityRelativeToWater,
            fluidForSubstance(SUBSTANCE_WATER).densityAtTemperature(tref), 1.0e-6));
    emit("builtin_water_cp_matches",
        waterOk && near(water.specificHeatJPerKgK,
            thermalForSubstance(SUBSTANCE_WATER).specificHeat, 1.0e-4));
    emit("builtin_water_conductivity_matches",
        waterOk && near(water.thermalConductivityWPerMK,
            thermalForSubstance(SUBSTANCE_WATER).conductivity, 1.0e-6));
    emit("builtin_water_viscosity_matches",
        waterOk && near(water.solverViscosity,
            fluidForSubstance(SUBSTANCE_WATER).viscosityAtTemperature(tref), 1.0e-7));
    emit("builtin_water_surface_matches",
        waterOk && near(water.solverSurfaceTension,
            fluidForSubstance(SUBSTANCE_WATER).surfaceTension, 1.0e-7));

    RuntimeLiquidPropertySample honey{};
    bool honeyOk = sampleRuntimeLiquidProperties(
        runtimeBuiltIn(SUBSTANCE_HONEY), tref, honey);
    emit("builtin_honey_liquid_samples", honeyOk && honey.valid);
    emit("builtin_air_not_liquid",
        !sampleRuntimeLiquidProperties(runtimeBuiltIn(SUBSTANCE_AIR), tref, honey));

    RuntimeGasPropertySample air{};
    bool airOk = sampleRuntimeGasProperties(runtimeBuiltIn(SUBSTANCE_AIR), tref, air);
    emit("builtin_air_gas_samples", airOk && air.valid);
    emit("builtin_air_molar_mass_reference",
        airOk && near(air.molarMassGPerMol, AIR_MOLAR_MASS_G_MOL, 1.0e-5));
    emit("builtin_air_cp_matches",
        airOk && near(air.specificHeatJPerKgK,
            gasPhaseSpecificHeat(SUBSTANCE_AIR), 1.0e-4));
    emit("builtin_air_conductivity_matches",
        airOk && near(air.thermalConductivityWPerMK,
            thermalForSubstance(SUBSTANCE_AIR).conductivity, 1.0e-6));

    RuntimeGasPropertySample hydrogen{};
    bool h2Ok = sampleRuntimeGasProperties(
        runtimeBuiltIn(SUBSTANCE_HYDROGEN), tref, hydrogen);
    emit("builtin_hydrogen_gas_samples", h2Ok && hydrogen.valid);
    emit("builtin_hydrogen_molar_mass_matches",
        h2Ok && near(hydrogen.molarMassGPerMol,
            chemicalForSubstance(SUBSTANCE_HYDROGEN).molarMass, 1.0e-5));
    emit("builtin_honey_not_gas",
        !sampleRuntimeGasProperties(runtimeBuiltIn(SUBSTANCE_HONEY), tref, hydrogen));

    RuntimeLiquidPropertySample invalidLiquid{};
    RuntimeGasPropertySample invalidGas{};
    emit("runtime_none_rejected",
        !sampleRuntimeLiquidProperties(runtimeNone(), tref, invalidLiquid)
            && !sampleRuntimeGasProperties(runtimeNone(), tref, invalidGas));
    emit("invalid_temperature_rejected",
        !sampleRuntimeLiquidProperties(runtimeBuiltIn(SUBSTANCE_WATER), 0.0f, invalidLiquid)
            && !sampleRuntimeGasProperties(runtimeBuiltIn(SUBSTANCE_AIR), -1.0f, invalidGas));

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
    bool attachA = cat.attachMolecularGraph(
        saceA.generatedId, "synthetic-structure-a", gA);
    bool attachB = cat.attachMolecularGraph(
        saceB.generatedId, "synthetic-structure-b", gB);

    SaceGeneratedRecord const *recA = cat.record(saceA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(saceB.generatedId);
    bool spawnA = recA && recA->spawnable;
    bool spawnB = recB && recB->spawnable;

    RuntimeSubstanceRef genA{}, genB{};
    bool regA = attachA && recA
        && registerGeneratedRuntimeMaterial(recA->recordId, genA);
    bool regB = attachB && recB
        && registerGeneratedRuntimeMaterial(recB->recordId, genB);
    emit("generated_a_registered", regA && runtimeSubstanceIsGenerated(genA));
    emit("generated_b_registered", regB && runtimeSubstanceIsGenerated(genB));

    RuntimeLiquidPropertySample aLiquid{}, bLiquid{};
    RuntimeGasPropertySample aGas{}, bGas{};
    bool aL = regA && sampleRuntimeLiquidProperties(genA, tref, aLiquid);
    bool bL = regB && sampleRuntimeLiquidProperties(genB, tref, bLiquid);
    bool aG = regA && sampleRuntimeGasProperties(genA, tref, aGas);
    bool bG = regB && sampleRuntimeGasProperties(genB, tref, bGas);
    emit("generated_a_liquid_samples", aL && aLiquid.valid);
    emit("generated_b_liquid_samples", bL && bLiquid.valid);
    emit("generated_a_gas_samples", aG && aGas.valid);
    emit("generated_b_gas_samples", bG && bGas.valid);

    GeneratedRuntimeMaterial const *matA = generatedRuntimeMaterial(genA);
    emit("generated_a_reference_density_matches_compiler",
        aL && matA && near(aLiquid.densityRelativeToWater,
            matA->profile.liquid.densityRelativeToWater, 1.0e-5));
    emit("generated_a_reference_liquid_cp_matches_compiler",
        aL && matA && near(aLiquid.specificHeatJPerKgK,
            matA->profile.liquid.specificHeatJPerKgK, 1.0e-2));
    emit("generated_a_reference_liquid_k_matches_compiler",
        aL && matA && near(aLiquid.thermalConductivityWPerMK,
            matA->profile.liquid.thermalConductivityWPerMK, 1.0e-6));
    emit("generated_a_reference_viscosity_matches_compiler",
        aL && matA && near(aLiquid.solverViscosity,
            matA->profile.liquid.solverViscosity, 1.0e-6));
    emit("generated_a_reference_surface_matches_compiler",
        aL && matA && near(aLiquid.solverSurfaceTension,
            matA->profile.liquid.solverSurfaceTension, 1.0e-6));
    emit("generated_a_reference_gas_cp_matches_compiler",
        aG && matA && near(aGas.specificHeatJPerKgK,
            matA->profile.gas.specificHeatJPerKgK, 1.0e-2));
    emit("generated_a_reference_gas_k_matches_compiler",
        aG && matA && near(aGas.thermalConductivityWPerMK,
            matA->profile.gas.thermalConductivityWPerMK, 1.0e-6));

    emit("generated_isomers_keep_liquid_differences",
        aL && bL
            && !near(aLiquid.densityRelativeToWater, bLiquid.densityRelativeToWater, 1.0e-4)
            && !near(aLiquid.solverViscosity, bLiquid.solverViscosity, 1.0e-5)
            && !near(aLiquid.solverSurfaceTension, bLiquid.solverSurfaceTension, 1.0e-5));
    emit("generated_isomers_keep_gas_differences",
        aG && bG
            && !near(aGas.thermalConductivityWPerMK,
                bGas.thermalConductivityWPerMK, 1.0e-5));

    constexpr float t2 = 320.0f;
    RuntimeLiquidPropertySample aLiquidT2{};
    RuntimeGasPropertySample aGasT2{};
    bool aLT2 = regA && sampleRuntimeLiquidProperties(genA, t2, aLiquidT2);
    bool aGT2 = regA && sampleRuntimeGasProperties(genA, t2, aGasT2);
    float directMuT2 = 0.0f;
    bool directMuOk = matA
        && saceCompiledSolverViscosityAtTemperature(
            matA->profile.liquid, t2, directMuT2);
    double directKL = 0.0, directKG = 0.0;
    bool directKLOk = matA
        && saceSatoRiedelLiquidThermalConductivityWPerMK(
            matA->profile.liquid.conductivityModel, t2, directKL);
    bool directKGOk = matA
        && saceGharagheiziGasThermalConductivityWPerMK(
            matA->profile.gas.conductivityModel, t2, directKG);
    emit("generated_temperature_liquid_query_uses_models",
        aLT2 && directMuOk && directKLOk
            && near(aLiquidT2.solverViscosity, directMuT2, 1.0e-6)
            && near(aLiquidT2.thermalConductivityWPerMK, directKL, 1.0e-6));
    emit("generated_temperature_gas_query_uses_models",
        aGT2 && directKGOk
            && near(aGasT2.thermalConductivityWPerMK, directKG, 1.0e-6));

    emit("generated_records_still_unspawnable",
        recA && recB && !spawnA && !spawnB && !recA->spawnable && !recB->spawnable);
    emit("substance_count_remains_12", SUBSTANCE_COUNT == 12,
        std::to_string(SUBSTANCE_COUNT));

    RuntimeSubstanceRef staleA = genA;
    cat.clear();
    RuntimeLiquidPropertySample staleLiquid{};
    RuntimeGasPropertySample staleGas{};
    emit("stale_generated_ref_rejected_after_catalog_clear",
        !sampleRuntimeLiquidProperties(staleA, tref, staleLiquid)
            && !sampleRuntimeGasProperties(staleA, tref, staleGas)
            && generatedRuntimeMaterial(staleA) == nullptr);

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
