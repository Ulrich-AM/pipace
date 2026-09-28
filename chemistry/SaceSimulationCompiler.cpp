#include "chemistry/SaceSimulationCompiler.h"

#include "chemistry/SaceIdentity.h"
#include "chemistry/SaceMolecule.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceTypes.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <string>

namespace {

bool finitePositive(double v) {
    return std::isfinite(v) && v > 0.0;
}

bool finitePositiveF(float v) {
    return std::isfinite(v) && v > 0.0f;
}

bool trustedGraphOk(SaceGeneratedRecord const &record) {
    if (!record.hasMolecularGraph)
        return false;
    if (validateMolecularGraph(record.molecularGraph, true) != SaceGraphValidation::Ok)
        return false;
    if (!record.molecularGraphStructureKey.empty()
        && record.molecularGraphStructureKey != record.structureKey)
        return false;
    return true;
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

bool compileGas(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceCompiledGasProfile &out)
{
    out = {};
    out.referenceTemperatureK = kSaceCompilerReferenceTemperatureK;

    SaceJobackIdealGasCpModel gasCp{};
    if (!saceBuildJobackIdealGasCpModel(bundle.groups, gasCp) || !gasCp.valid)
        return false;

    SaceGharagheiziGasConductivityModel gasK{};
    if (!saceBuildGharagheiziGasConductivityModelFromJoback(graph, bundle, gasK) || !gasK.valid)
        return false;

    double Tref = kSaceHeatCapacityReferenceK;
    double cpMolar = 0.0;
    if (!saceJobackIdealGasHeatCapacityJPerMolK(gasCp, Tref, cpMolar))
        return false;
    double mwKg = gasK.molarMassGPerMol * 1.0e-3;
    double cpSpec = 0.0;
    if (!saceMolarHeatCapacityToSpecificJPerKgK(cpMolar, mwKg, cpSpec))
        return false;

    double kG = 0.0;
    if (!saceGharagheiziGasThermalConductivityWPerMK(gasK, Tref, kG))
        return false;

    if (!finitePositive(gasK.molarMassGPerMol) || !finitePositive(cpSpec) || !finitePositive(kG))
        return false;

    out.heatCapacityModel = gasCp;
    out.conductivityModel = gasK;
    out.molarMassGPerMol = static_cast<float>(gasK.molarMassGPerMol);
    out.specificHeatJPerKgK = static_cast<float>(cpSpec);
    out.thermalConductivityWPerMK = static_cast<float>(kG);
    out.valid = true;
    return true;
}

bool compileLiquid(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceCompiledLiquidProfile &out)
{
    out = {};
    out.referenceTemperatureK = kSaceCompilerReferenceTemperatureK;

    SaceCostaldLiquidDensityModel dens{};
    if (!saceBuildCostaldLiquidDensityModelFromJoback(graph, bundle, dens) || !dens.valid)
        return false;
    SaceRowlinsonPolingLiquidCpModel liqCp{};
    if (!saceBuildRowlinsonPolingLiquidCpModelFromJoback(graph, bundle, liqCp) || !liqCp.valid)
        return false;
    SaceJobackLiquidViscosityModel visc{};
    if (!saceBuildJobackLiquidViscosityModel(graph, bundle, visc) || !visc.valid)
        return false;
    SaceSastriRaoSurfaceTensionModel sigma{};
    if (!saceBuildSastriRaoSurfaceTensionModelFromJoback(graph, bundle, sigma) || !sigma.valid)
        return false;
    SaceSatoRiedelLiquidConductivityModel liqK{};
    if (!saceBuildSatoRiedelLiquidConductivityModelFromJoback(graph, bundle, liqK) || !liqK.valid)
        return false;

    double Tref = kSaceHeatCapacityReferenceK;
    double rho = 0.0;
    if (!saceCostaldSaturatedLiquidDensityKgPerM3(dens, Tref, rho))
        return false;
    float rhoRel = 0.0f;
    if (!sacePhysicalDensityToWaterRelative(rho, rhoRel))
        return false;

    double cpMolar = 0.0;
    if (!saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqCp, Tref, cpMolar))
        return false;
    double cpSpec = 0.0;
    if (!saceMolarHeatCapacityToSpecificJPerKgK(cpMolar, liqCp.molarMassKgPerMol, cpSpec))
        return false;

    double kL = 0.0;
    if (!saceSatoRiedelLiquidThermalConductivityWPerMK(liqK, Tref, kL))
        return false;

    double muPaS = 0.0;
    if (!saceJobackLiquidViscosityPaS(visc, Tref, muPaS))
        return false;
    float solverMu = 0.0f;
    if (!saceCalibrateSolverViscosityFromPhysicalPaS(muPaS, solverMu))
        return false;

    double sigmaN = 0.0;
    if (!saceSastriRaoSurfaceTensionNPerM(sigma, Tref, sigmaN))
        return false;
    float solverSigma = 0.0f;
    if (!saceCalibrateSolverSurfaceTensionFromPhysicalNPerM(sigmaN, solverSigma))
        return false;

    if (!std::isfinite(visc.aKelvin))
        return false;
    if (!finitePositiveF(rhoRel) || !finitePositive(cpSpec) || !finitePositive(kL)
        || !finitePositiveF(solverMu) || !finitePositiveF(solverSigma))
        return false;

    out.densityModel = dens;
    out.heatCapacityModel = liqCp;
    out.viscosityModel = visc;
    out.surfaceTensionModel = sigma;
    out.conductivityModel = liqK;
    out.densityRelativeToWater = rhoRel;
    out.specificHeatJPerKgK = static_cast<float>(cpSpec);
    out.thermalConductivityWPerMK = static_cast<float>(kL);
    out.solverViscosity = solverMu;
    out.solverSurfaceTension = solverSigma;
    out.viscArrheniusK = static_cast<float>(visc.aKelvin);
    out.valid = true;
    return true;
}

bool compilePhase(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceCompiledPhaseProfile &out)
{
    out = {};
    out.referencePressurePa = static_cast<float>(kLeeKeslerAtmPa);

    SaceLeeKeslerVaporModel vapor{};
    if (!saceBuildLeeKeslerVaporModelFromJoback(bundle, vapor) || !vapor.valid)
        return false;
    SaceWatsonVaporizationModel hvap{};
    if (!saceBuildWatsonVaporizationModelFromJoback(graph, bundle, hvap) || !hvap.valid)
        return false;
    if (!bundle.normalBoilingPointK.known || !bundle.criticalTemperatureK.known
        || !bundle.criticalPressurePa.known || !bundle.acentricFactor.known)
        return false;
    if (!finitePositive(vapor.normalBoilingPointK) || !finitePositive(vapor.criticalTemperatureK)
        || !finitePositive(vapor.criticalPressurePa) || !std::isfinite(vapor.acentricFactor))
        return false;

    double mwG = 0.0;
    if (hvap.molarMassKgPerMol > 0.0 && std::isfinite(hvap.molarMassKgPerMol))
        mwG = hvap.molarMassKgPerMol * 1000.0;

    out.vaporPressureModel = vapor;
    out.vaporizationModel = hvap;
    out.normalBoilingPointK = static_cast<float>(vapor.normalBoilingPointK);
    out.criticalTemperatureK = static_cast<float>(vapor.criticalTemperatureK);
    out.criticalPressurePa = static_cast<float>(vapor.criticalPressurePa);
    out.acentricFactor = vapor.acentricFactor;
    out.molarMassGPerMol = mwG;
    out.liquidGasValid = true;
    return true;
}

bool finitePositiveFields(SaceCompiledLiquidProfile const &liq) {
    return liq.valid
        && finitePositiveF(liq.densityRelativeToWater)
        && finitePositiveF(liq.specificHeatJPerKgK)
        && finitePositiveF(liq.thermalConductivityWPerMK)
        && finitePositiveF(liq.solverViscosity)
        && finitePositiveF(liq.solverSurfaceTension);
}

bool finitePositiveFields(SaceCompiledGasProfile const &gas) {
    return gas.valid
        && finitePositiveF(gas.molarMassGPerMol)
        && finitePositiveF(gas.specificHeatJPerKgK)
        && finitePositiveF(gas.thermalConductivityWPerMK);
}

bool near(double a, double b, double tol) {
    if (!std::isfinite(a) || !std::isfinite(b))
        return false;
    return std::fabs(a - b) <= tol;
}

} // namespace

void saceResetCompiledSimulationProfile(SaceCompiledSimulationProfile &out) {
    out = {};
}

bool sacePhysicalDensityToWaterRelative(double densityKgPerM3, float &outRelative) {
    outRelative = 0.0f;
    if (!finitePositive(densityKgPerM3) || !(WATER_DENSITY_KG_M3 > 0.0f))
        return false;
    double rel = densityKgPerM3 / static_cast<double>(WATER_DENSITY_KG_M3);
    if (!finitePositive(rel))
        return false;
    outRelative = static_cast<float>(rel);
    return true;
}

bool saceCalibrateSolverViscosityFromPhysicalPaS(double physicalPaS, float &outSolver) {
    outSolver = 0.0f;
    if (!finitePositive(physicalPaS) || !finitePositive(kPhysicalWaterViscosityAt298KPaS))
        return false;
    float refSolver = sandboxReferenceLiquid().viscosityAtTemperature(kSaceCompilerReferenceTemperatureK);
    if (!finitePositiveF(refSolver))
        return false;
    double solver = static_cast<double>(refSolver)
        * physicalPaS / kPhysicalWaterViscosityAt298KPaS;
    if (!std::isfinite(solver) || !(solver > 0.0))
        return false;
    outSolver = static_cast<float>(solver);
    return true;
}

bool saceCalibrateSolverSurfaceTensionFromPhysicalNPerM(double physicalNPerM, float &outSolver) {
    outSolver = 0.0f;
    if (!finitePositive(physicalNPerM) || !finitePositive(kPhysicalWaterSurfaceTensionAt298KNPerM))
        return false;
    float refSolver = sandboxReferenceLiquid().surfaceTension;
    if (!finitePositiveF(refSolver))
        return false;
    double solver = static_cast<double>(refSolver)
        * physicalNPerM / kPhysicalWaterSurfaceTensionAt298KNPerM;
    if (!std::isfinite(solver) || !(solver > 0.0))
        return false;
    outSolver = static_cast<float>(solver);
    return true;
}

bool saceCompiledSolverViscosityAtTemperature(
    SaceCompiledLiquidProfile const &liquid,
    float temperatureK,
    float &outSolver)
{
    outSolver = 0.0f;
    if (!liquid.valid || !finitePositiveF(liquid.solverViscosity))
        return false;
    float T = temperatureK;
    if (!(T > 1.0f) || !std::isfinite(T))
        return false;
    float Tref = liquid.referenceTemperatureK > 1.0f
        ? liquid.referenceTemperatureK
        : kSaceCompilerReferenceTemperatureK;
    float mu = liquid.solverViscosity
        * std::exp(liquid.viscArrheniusK * (1.0f / T - 1.0f / Tref));
    if (!std::isfinite(mu) || !(mu > 0.0f))
        return false;
    outSolver = mu;
    return true;
}

bool saceCompileSimulationProfile(
    SaceGeneratedRecord const &record,
    SaceCompiledSimulationProfile &out)
{
    saceResetCompiledSimulationProfile(out);
    out.sourceRecord = record.recordId;

    SaceSimulationReadiness ready = saceAssessSimulationReadiness(record);
    if (!trustedGraphOk(record))
        return false;

    SaceJobackEstimateBundle bundle{};
    if (!saceEstimateJobackBundle(record.molecularGraph, bundle)) {
        // Unsupported chemistry: no Water/Air fallback.
        return false;
    }

    bool gasOk = false;
    bool liquidOk = false;
    bool phaseOk = false;
    if (ready.liveGasReady)
        gasOk = compileGas(record.molecularGraph, bundle, out.gas);
    if (ready.liveLiquidReady)
        liquidOk = compileLiquid(record.molecularGraph, bundle, out.liquid);
    if (ready.liquidGasEquilibriumReady)
        phaseOk = compilePhase(record.molecularGraph, bundle, out.phase);

    out.gasReady = gasOk && out.gas.valid;
    out.liquidReady = liquidOk && out.liquid.valid;
    out.liquidGasPhaseChangeReady = out.gasReady && out.liquidReady
        && phaseOk && out.phase.liquidGasValid;
    out.valid = out.gasReady || out.liquidReady || out.phase.liquidGasValid;
    return out.valid;
}

void runSaceCompilerDiagnostics() {
    std::ofstream out(miscFile("sace_compiler_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

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
    SaceSubstanceRef refA = cat.resolve(idA, true);
    SaceSubstanceRef refB = cat.resolve(idB, true);
    cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    cat.attachMolecularGraph(refB.generatedId, "synthetic-structure-b", gB);
    SaceGeneratedRecord const *recA = cat.record(refA.generatedId);
    SaceGeneratedRecord const *recB = cat.record(refB.generatedId);

    char sigBefore[kSaceSignatureCap]{};
    ChemicalIdentity viewBefore = recA ? cat.identityView(*recA) : ChemicalIdentity{};
    writeChemicalSignature(viewBefore, sigBefore, kSaceSignatureCap);
    uint32_t ordinalBefore = recA ? recA->displayOrdinal : 0;
    SaceRecordId idBefore = recA ? recA->recordId : kSaceRecordNone;
    std::string structureBefore = recA ? recA->structureKey : "";
    bool spawnBefore = recA && recA->spawnable;

    SaceSimulationReadiness aReady = recA ? saceAssessSimulationReadiness(*recA) : SaceSimulationReadiness{};
    SaceCompiledSimulationProfile a{};
    bool aCompile = recA && saceCompileSimulationProfile(*recA, a);
    emit("a_readiness_succeeds",
        recA && aReady.liveGasReady && aReady.liveLiquidReady
            && aReady.liveLiquidGasPhaseChangeReady, "");
    emit("a_compiler_succeeds", aCompile && a.valid, "");
    emit("a_liquid_valid", a.liquid.valid && a.liquidReady, "");
    emit("a_gas_valid", a.gas.valid && a.gasReady, "");
    emit("a_liquid_gas_phase_valid",
        a.phase.liquidGasValid && a.liquidGasPhaseChangeReady, "");
    emit("a_source_record_preserved",
        recA && a.sourceRecord == recA->recordId && a.sourceRecord != kSaceRecordNone, "");
    emit("a_liquid_fields_finite_positive", finitePositiveFields(a.liquid), "");
    emit("a_gas_fields_finite_positive", finitePositiveFields(a.gas), "");

    SaceSimulationReadiness bReady = recB ? saceAssessSimulationReadiness(*recB) : SaceSimulationReadiness{};
    SaceCompiledSimulationProfile b{};
    bool bCompile = recB && saceCompileSimulationProfile(*recB, b);
    emit("b_readiness_succeeds",
        recB && bReady.liveGasReady && bReady.liveLiquidReady
            && bReady.liveLiquidGasPhaseChangeReady, "");
    emit("b_compiler_succeeds", bCompile && b.valid, "");
    emit("b_liquid_valid", b.liquid.valid && b.liquidReady, "");
    emit("b_gas_valid", b.gas.valid && b.gasReady, "");
    emit("b_liquid_gas_phase_valid",
        b.phase.liquidGasValid && b.liquidGasPhaseChangeReady, "");
    emit("b_source_record_preserved",
        recB && b.sourceRecord == recB->recordId && b.sourceRecord != kSaceRecordNone, "");
    emit("b_liquid_fields_finite_positive", finitePositiveFields(b.liquid), "");
    emit("b_gas_fields_finite_positive", finitePositiveFields(b.gas), "");

    emit("a_b_density_differs",
        a.liquid.valid && b.liquid.valid
            && a.liquid.densityRelativeToWater != b.liquid.densityRelativeToWater, "");
    emit("a_b_viscosity_differs",
        a.liquid.valid && b.liquid.valid
            && a.liquid.solverViscosity != b.liquid.solverViscosity, "");
    emit("a_b_surface_tension_differs",
        a.liquid.valid && b.liquid.valid
            && a.liquid.solverSurfaceTension != b.liquid.solverSurfaceTension, "");
    emit("a_b_liquid_conductivity_differs",
        a.liquid.valid && b.liquid.valid
            && a.liquid.thermalConductivityWPerMK != b.liquid.thermalConductivityWPerMK, "");
    emit("a_b_gas_conductivity_differs",
        a.gas.valid && b.gas.valid
            && a.gas.thermalConductivityWPerMK != b.gas.thermalConductivityWPerMK, "");
    emit("a_b_liquid_cp_differs",
        a.liquid.valid && b.liquid.valid
            && a.liquid.specificHeatJPerKgK != b.liquid.specificHeatJPerKgK, "");
    emit("a_b_boiling_point_differs",
        a.phase.liquidGasValid && b.phase.liquidGasValid
            && a.phase.normalBoilingPointK != b.phase.normalBoilingPointK, "");

    double rhoA = 0.0, rhoB = 0.0;
    bool rhoAOk = a.liquid.valid
        && saceCostaldSaturatedLiquidDensityKgPerM3(
            a.liquid.densityModel, kSaceHeatCapacityReferenceK, rhoA);
    bool rhoBOk = b.liquid.valid
        && saceCostaldSaturatedLiquidDensityKgPerM3(
            b.liquid.densityModel, kSaceHeatCapacityReferenceK, rhoB);
    emit("a_density_relative_roundtrip",
        rhoAOk && near(static_cast<double>(a.liquid.densityRelativeToWater) * WATER_DENSITY_KG_M3,
            rhoA, 1.0e-3),
        std::to_string(a.liquid.densityRelativeToWater));
    emit("b_density_relative_roundtrip",
        rhoBOk && near(static_cast<double>(b.liquid.densityRelativeToWater) * WATER_DENSITY_KG_M3,
            rhoB, 1.0e-3),
        std::to_string(b.liquid.densityRelativeToWater));
    emit("a_density_relative_near_expected",
        a.liquid.valid && near(a.liquid.densityRelativeToWater, 0.822, 0.03),
        std::to_string(a.liquid.densityRelativeToWater));
    emit("b_density_relative_near_expected",
        b.liquid.valid && near(b.liquid.densityRelativeToWater, 0.700, 0.03),
        std::to_string(b.liquid.densityRelativeToWater));

    float waterSolverMu = 0.0f;
    bool waterMuBridge = saceCalibrateSolverViscosityFromPhysicalPaS(
        kPhysicalWaterViscosityAt298KPaS, waterSolverMu);
    emit("viscosity_water_unit_bridge",
        waterMuBridge && near(waterSolverMu,
            sandboxReferenceLiquid().viscosityAtTemperature(kSaceCompilerReferenceTemperatureK),
            1.0e-8),
        std::to_string(waterSolverMu));
    float waterSolverSigma = 0.0f;
    bool waterSigmaBridge = saceCalibrateSolverSurfaceTensionFromPhysicalNPerM(
        kPhysicalWaterSurfaceTensionAt298KNPerM, waterSolverSigma);
    emit("surface_tension_water_unit_bridge",
        waterSigmaBridge && near(waterSolverSigma,
            sandboxReferenceLiquid().surfaceTension, 1.0e-8),
        std::to_string(waterSolverSigma));

    double T2 = 320.0;
    double muPhysRef = 0.0, muPhysT2 = 0.0;
    float muSolRef = 0.0f, muSolT2 = 0.0f;
    bool viscShape = a.liquid.valid
        && saceJobackLiquidViscosityPaS(a.liquid.viscosityModel, kSaceHeatCapacityReferenceK, muPhysRef)
        && saceJobackLiquidViscosityPaS(a.liquid.viscosityModel, T2, muPhysT2)
        && saceCompiledSolverViscosityAtTemperature(a.liquid, kSaceCompilerReferenceTemperatureK, muSolRef)
        && saceCompiledSolverViscosityAtTemperature(a.liquid, static_cast<float>(T2), muSolT2)
        && muPhysRef > 0.0 && muSolRef > 0.0f
        && near(muSolT2 / muSolRef, muPhysT2 / muPhysRef, 1.0e-5);
    emit("joback_viscosity_temperature_shape", viscShape,
        viscShape ? std::to_string(muSolT2 / muSolRef) : "");

    double satoA = 0.0, gharA = 0.0;
    bool satoOk = a.liquid.valid && saceSatoRiedelLiquidThermalConductivityWPerMK(
        a.liquid.conductivityModel, kSaceHeatCapacityReferenceK, satoA);
    bool gharOk = a.gas.valid && saceGharagheiziGasThermalConductivityWPerMK(
        a.gas.conductivityModel, kSaceHeatCapacityReferenceK, gharA);
    emit("liquid_conductivity_from_sato_riedel",
        satoOk && near(a.liquid.thermalConductivityWPerMK, satoA, 1.0e-6),
        std::to_string(a.liquid.thermalConductivityWPerMK));
    emit("gas_conductivity_from_gharagheizi",
        gharOk && near(a.gas.thermalConductivityWPerMK, gharA, 1.0e-6),
        std::to_string(a.gas.thermalConductivityWPerMK));
    emit("liquid_and_gas_conductivity_not_collapsed",
        a.liquid.valid && a.gas.valid
            && a.liquid.conductivityModel.valid
            && a.gas.conductivityModel.valid
            && &a.liquid.thermalConductivityWPerMK != &a.gas.thermalConductivityWPerMK
            && satoOk && gharOk
            && near(a.liquid.thermalConductivityWPerMK, satoA, 1.0e-6)
            && near(a.gas.thermalConductivityWPerMK, gharA, 1.0e-6)
            && !near(a.liquid.thermalConductivityWPerMK, gharA, 1.0e-4)
            && !near(a.gas.thermalConductivityWPerMK, satoA, 1.0e-4),
        "");

    emit("physical_models_retained",
        a.liquid.valid && a.gas.valid && a.phase.liquidGasValid
            && a.liquid.densityModel.valid
            && a.liquid.heatCapacityModel.valid
            && a.liquid.viscosityModel.valid
            && a.liquid.surfaceTensionModel.valid
            && a.liquid.conductivityModel.valid
            && a.gas.heatCapacityModel.valid
            && a.gas.conductivityModel.valid
            && a.phase.vaporPressureModel.valid
            && a.phase.vaporizationModel.valid, "");

    ElementCount c2h4[] = {{kAtomicCarbon, 2}, {kAtomicHydrogen, 4}};
    ChemicalIdentity idE = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H4", "synthetic-ethene", c2h4, 2);
    SaceSubstanceRef refE = cat.resolve(idE, true);
    cat.attachMolecularGraph(refE.generatedId, "synthetic-ethene", etheneGraph());
    SaceGeneratedRecord const *recE = cat.record(refE.generatedId);
    SaceSimulationReadiness eReady = recE ? saceAssessSimulationReadiness(*recE) : SaceSimulationReadiness{};
    SaceCompiledSimulationProfile e{};
    bool eCompile = recE && saceCompileSimulationProfile(*recE, e);
    emit("unsupported_readiness_false",
        recE && !eReady.liveGasReady && !eReady.liveLiquidReady
            && !eReady.liveLiquidGasPhaseChangeReady, "");
    emit("unsupported_compiler_refuses",
        recE && !eCompile && !e.valid && !e.liquid.valid && !e.gas.valid
            && !e.phase.liquidGasValid, "");
    emit("unsupported_no_nan_fallback",
        recE && !e.liquid.valid && !e.gas.valid
            && e.liquid.densityRelativeToWater == 0.0f
            && e.liquid.specificHeatJPerKgK == 0.0f
            && e.gas.thermalConductivityWPerMK == 0.0f
            && std::isfinite(e.liquid.solverViscosity)
            && std::isfinite(e.gas.specificHeatJPerKgK), "");

    recA = cat.record(refA.generatedId);
    char sigAfter[kSaceSignatureCap]{};
    ChemicalIdentity viewAfter = recA ? cat.identityView(*recA) : ChemicalIdentity{};
    writeChemicalSignature(viewAfter, sigAfter, kSaceSignatureCap);
    emit("compiler_does_not_alter_signature",
        recA && std::strcmp(sigBefore, sigAfter) == 0, sigAfter);
    emit("compiler_does_not_alter_record_id",
        recA && recA->recordId == idBefore, "");
    emit("compiler_does_not_alter_display_ordinal",
        recA && recA->displayOrdinal == ordinalBefore, "");
    emit("compiler_does_not_alter_structure_key",
        recA && recA->structureKey == structureBefore, "");
    emit("generated_record_unspawnable",
        recA && !recA->spawnable && recB && !recB->spawnable && !spawnBefore, "");
    emit("substance_count_remains_12", SUBSTANCE_COUNT == 12, std::to_string(SUBSTANCE_COUNT));
    emit("profile_valid_means_any_subprofile",
        a.valid && (a.gas.valid || a.liquid.valid || a.phase.liquidGasValid)
            && !e.valid, "");
    emit("a_compiled_reference_values", a.valid,
        std::string("rhoRel=") + std::to_string(a.liquid.densityRelativeToWater)
            + " liqCp=" + std::to_string(a.liquid.specificHeatJPerKgK)
            + " gasCp=" + std::to_string(a.gas.specificHeatJPerKgK)
            + " kL=" + std::to_string(a.liquid.thermalConductivityWPerMK)
            + " kG=" + std::to_string(a.gas.thermalConductivityWPerMK)
            + " muSol=" + std::to_string(a.liquid.solverViscosity)
            + " sigSol=" + std::to_string(a.liquid.solverSurfaceTension)
            + " Tb=" + std::to_string(a.phase.normalBoilingPointK)
            + " MW=" + std::to_string(a.gas.molarMassGPerMol));
    emit("b_compiled_reference_values", b.valid,
        std::string("rhoRel=") + std::to_string(b.liquid.densityRelativeToWater)
            + " liqCp=" + std::to_string(b.liquid.specificHeatJPerKgK)
            + " gasCp=" + std::to_string(b.gas.specificHeatJPerKgK)
            + " kL=" + std::to_string(b.liquid.thermalConductivityWPerMK)
            + " kG=" + std::to_string(b.gas.thermalConductivityWPerMK)
            + " muSol=" + std::to_string(b.liquid.solverViscosity)
            + " sigSol=" + std::to_string(b.liquid.solverSurfaceTension)
            + " Tb=" + std::to_string(b.phase.normalBoilingPointK)
            + " MW=" + std::to_string(b.gas.molarMassGPerMol));

    cat.clear();
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
