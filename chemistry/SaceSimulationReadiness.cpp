#include "chemistry/SaceSimulationReadiness.h"

#include "chemistry/SaceEstimation.h"
#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceTypes.h"

#include <cstring>
#include <fstream>
#include <string>

namespace {

void markMissing(SaceSimulationReadiness &r, SaceSimulationRequirement req) {
    r.missingMask |= static_cast<uint64_t>(req);
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

bool molarMassOk(SaceGeneratedRecord const &record, bool graphOk) {
    if (record.properties.molarMassGPerMol.known
        && record.properties.molarMassGPerMol.value > 0.0f)
        return true;
    if (!graphOk)
        return false;
    ElementalComposition elemental{};
    if (!elementalCompositionFromGraph(record.molecularGraph, elemental))
        return false;
    SaceScalarProperty mm{};
    return saceDeriveMolarMass(elemental, mm) && mm.known && mm.value > 0.0f;
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

bool classMatches(SaceSimulationReadiness const &a, SaceSimulationReadiness const &b) {
    return a.gasThermoReady == b.gasThermoReady
        && a.liquidThermoReady == b.liquidThermoReady
        && a.liquidGasEquilibriumReady == b.liquidGasEquilibriumReady
        && a.liveGasReady == b.liveGasReady
        && a.liveLiquidReady == b.liveLiquidReady
        && a.liveLiquidGasPhaseChangeReady == b.liveLiquidGasPhaseChangeReady
        && saceReadinessMissing(a, SaceSimulationRequirement::LiquidViscosity)
            == saceReadinessMissing(b, SaceSimulationRequirement::LiquidViscosity)
        && saceReadinessMissing(a, SaceSimulationRequirement::LiquidSurfaceTension)
            == saceReadinessMissing(b, SaceSimulationRequirement::LiquidSurfaceTension)
        && saceReadinessMissing(a, SaceSimulationRequirement::LiquidThermalConductivity)
            == saceReadinessMissing(b, SaceSimulationRequirement::LiquidThermalConductivity)
        && saceReadinessMissing(a, SaceSimulationRequirement::GasThermalConductivity)
            == saceReadinessMissing(b, SaceSimulationRequirement::GasThermalConductivity);
}

} // namespace

char const *saceSimulationRequirementKey(SaceSimulationRequirement requirement) {
    switch (requirement) {
        case SaceSimulationRequirement::MolecularGraph: return "molecular_graph";
        case SaceSimulationRequirement::MolarMass: return "molar_mass";
        case SaceSimulationRequirement::NormalBoilingPoint: return "normal_boiling_point";
        case SaceSimulationRequirement::CriticalTemperature: return "critical_temperature";
        case SaceSimulationRequirement::CriticalPressure: return "critical_pressure";
        case SaceSimulationRequirement::AcentricFactor: return "acentric_factor";
        case SaceSimulationRequirement::VaporPressureModel: return "vapor_pressure_model";
        case SaceSimulationRequirement::VaporizationEnthalpyModel: return "vaporization_enthalpy_model";
        case SaceSimulationRequirement::GasHeatCapacityModel: return "gas_heat_capacity_model";
        case SaceSimulationRequirement::LiquidHeatCapacityModel: return "liquid_heat_capacity_model";
        case SaceSimulationRequirement::LiquidDensityModel: return "liquid_density_model";
        case SaceSimulationRequirement::LiquidViscosity: return "liquid_viscosity";
        case SaceSimulationRequirement::LiquidSurfaceTension: return "liquid_surface_tension";
        case SaceSimulationRequirement::LiquidThermalConductivity: return "liquid_thermal_conductivity";
        case SaceSimulationRequirement::GasThermalConductivity: return "gas_thermal_conductivity";
    }
    return "";
}

SaceSimulationReadiness saceAssessSimulationReadiness(SaceGeneratedRecord const &record) {
    SaceSimulationReadiness r{};

    bool graphOk = trustedGraphOk(record);
    if (!graphOk)
        markMissing(r, SaceSimulationRequirement::MolecularGraph);

    bool massOk = molarMassOk(record, graphOk);
    if (!massOk)
        markMissing(r, SaceSimulationRequirement::MolarMass);

    SaceJobackEstimateBundle bundle{};
    bool bundleOk = graphOk && saceEstimateJobackBundle(record.molecularGraph, bundle);

    bool tbOk = bundleOk && bundle.normalBoilingPointK.known;
    if (!tbOk)
        markMissing(r, SaceSimulationRequirement::NormalBoilingPoint);
    bool tcOk = bundleOk && bundle.criticalTemperatureK.known;
    if (!tcOk)
        markMissing(r, SaceSimulationRequirement::CriticalTemperature);
    bool pcOk = bundleOk && bundle.criticalPressurePa.known;
    if (!pcOk)
        markMissing(r, SaceSimulationRequirement::CriticalPressure);
    bool omegaOk = bundleOk && bundle.acentricFactor.known;
    if (!omegaOk)
        markMissing(r, SaceSimulationRequirement::AcentricFactor);

    SaceLeeKeslerVaporModel vapor{};
    bool vaporOk = bundleOk && saceBuildLeeKeslerVaporModelFromJoback(bundle, vapor) && vapor.valid;
    if (!vaporOk)
        markMissing(r, SaceSimulationRequirement::VaporPressureModel);

    SaceWatsonVaporizationModel hvap{};
    bool hvapOk = bundleOk
        && saceBuildWatsonVaporizationModelFromJoback(record.molecularGraph, bundle, hvap)
        && hvap.valid;
    if (!hvapOk)
        markMissing(r, SaceSimulationRequirement::VaporizationEnthalpyModel);

    SaceJobackIdealGasCpModel gasCp{};
    bool gasCpOk = bundleOk && saceBuildJobackIdealGasCpModel(bundle.groups, gasCp) && gasCp.valid;
    if (!gasCpOk)
        markMissing(r, SaceSimulationRequirement::GasHeatCapacityModel);

    SaceRowlinsonPolingLiquidCpModel liqCp{};
    bool liqCpOk = bundleOk
        && saceBuildRowlinsonPolingLiquidCpModelFromJoback(record.molecularGraph, bundle, liqCp)
        && liqCp.valid;
    if (!liqCpOk)
        markMissing(r, SaceSimulationRequirement::LiquidHeatCapacityModel);

    SaceCostaldLiquidDensityModel dens{};
    bool densOk = bundleOk
        && saceBuildCostaldLiquidDensityModelFromJoback(record.molecularGraph, bundle, dens)
        && dens.valid;
    if (!densOk)
        markMissing(r, SaceSimulationRequirement::LiquidDensityModel);

    SaceJobackLiquidViscosityModel visc{};
    bool viscOk = bundleOk
        && saceBuildJobackLiquidViscosityModel(record.molecularGraph, bundle, visc)
        && visc.valid;
    if (!viscOk)
        markMissing(r, SaceSimulationRequirement::LiquidViscosity);

    SaceSastriRaoSurfaceTensionModel sigma{};
    bool sigmaOk = bundleOk
        && saceBuildSastriRaoSurfaceTensionModelFromJoback(record.molecularGraph, bundle, sigma)
        && sigma.valid;
    if (!sigmaOk)
        markMissing(r, SaceSimulationRequirement::LiquidSurfaceTension);

    SaceSatoRiedelLiquidConductivityModel liqK{};
    bool liquidKOk = bundleOk
        && saceBuildSatoRiedelLiquidConductivityModelFromJoback(record.molecularGraph, bundle, liqK)
        && liqK.valid;
    if (!liquidKOk)
        markMissing(r, SaceSimulationRequirement::LiquidThermalConductivity);

    SaceGharagheiziGasConductivityModel gasK{};
    bool gasKOk = bundleOk
        && saceBuildGharagheiziGasConductivityModelFromJoback(record.molecularGraph, bundle, gasK)
        && gasK.valid;
    if (!gasKOk)
        markMissing(r, SaceSimulationRequirement::GasThermalConductivity);

    r.gasThermoReady = graphOk && massOk && gasCpOk;
    r.liquidThermoReady = graphOk && massOk && liqCpOk && densOk;
    r.liquidGasEquilibriumReady = graphOk && massOk && tbOk && tcOk && pcOk && omegaOk
        && vaporOk && hvapOk && gasCpOk && liqCpOk && densOk;
    r.liveGasReady = r.gasThermoReady
        && !saceReadinessMissing(r, SaceSimulationRequirement::GasThermalConductivity);
    r.liveLiquidReady = r.liquidThermoReady
        && !saceReadinessMissing(r, SaceSimulationRequirement::LiquidViscosity)
        && !saceReadinessMissing(r, SaceSimulationRequirement::LiquidSurfaceTension)
        && !saceReadinessMissing(r, SaceSimulationRequirement::LiquidThermalConductivity);
    r.liveLiquidGasPhaseChangeReady = r.liveGasReady && r.liveLiquidReady
        && r.liquidGasEquilibriumReady;
    return r;
}

void runSaceReadinessDiagnostics() {
    std::ofstream out(miscFile("sace_readiness_diag.tsv"));
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

    SaceSimulationReadiness a = recA ? saceAssessSimulationReadiness(*recA) : SaceSimulationReadiness{};
    emit("a_assessment_succeeds", recA != nullptr, "");
    emit("a_gas_thermo_ready", a.gasThermoReady, "");
    emit("a_liquid_thermo_ready", a.liquidThermoReady, "");
    emit("a_liquid_gas_equilibrium_ready", a.liquidGasEquilibriumReady, "");
    emit("a_live_gas_ready", a.liveGasReady, "");
    emit("a_live_liquid_ready", a.liveLiquidReady, "");
    emit("a_live_phase_change_ready", a.liveLiquidGasPhaseChangeReady, "");
    emit("a_gas_conductivity_not_missing",
        !saceReadinessMissing(a, SaceSimulationRequirement::GasThermalConductivity), "");
    emit("a_liquid_conductivity_not_missing",
        !saceReadinessMissing(a, SaceSimulationRequirement::LiquidThermalConductivity), "");
    emit("a_liquid_viscosity_not_missing",
        !saceReadinessMissing(a, SaceSimulationRequirement::LiquidViscosity), "");
    emit("a_liquid_surface_tension_not_missing",
        !saceReadinessMissing(a, SaceSimulationRequirement::LiquidSurfaceTension), "");

    SaceSimulationReadiness b = recB ? saceAssessSimulationReadiness(*recB) : SaceSimulationReadiness{};
    emit("b_gas_thermo_ready", b.gasThermoReady, "");
    emit("b_liquid_thermo_ready", b.liquidThermoReady, "");
    emit("b_liquid_gas_equilibrium_ready", b.liquidGasEquilibriumReady, "");
    emit("b_live_gas_ready", b.liveGasReady, "");
    emit("b_live_liquid_ready", b.liveLiquidReady, "");
    emit("b_live_phase_change_ready", b.liveLiquidGasPhaseChangeReady, "");
    emit("b_gas_conductivity_not_missing",
        !saceReadinessMissing(b, SaceSimulationRequirement::GasThermalConductivity), "");
    emit("b_liquid_conductivity_not_missing",
        !saceReadinessMissing(b, SaceSimulationRequirement::LiquidThermalConductivity), "");
    emit("b_liquid_viscosity_not_missing",
        !saceReadinessMissing(b, SaceSimulationRequirement::LiquidViscosity), "");
    emit("b_liquid_surface_tension_not_missing",
        !saceReadinessMissing(b, SaceSimulationRequirement::LiquidSurfaceTension), "");
    emit("a_b_readiness_class_matches", recA && recB && classMatches(a, b), "");
    emit("a_b_numerical_properties_differ",
        recA && recB
            && recA->properties.normalBoilingPointK.known
            && recB->properties.normalBoilingPointK.known
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && recB->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && recA->properties.normalBoilingPointK.value
                != recB->properties.normalBoilingPointK.value
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.value
                != recB->properties.saturatedLiquidDensityAt298KKgPerM3.value, "");

    ElementCount c2h4[] = {{kAtomicCarbon, 2}, {kAtomicHydrogen, 4}};
    ChemicalIdentity idE = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H4", "synthetic-ethene", c2h4, 2);
    SaceSubstanceRef refE = cat.resolve(idE, true);
    cat.attachMolecularGraph(refE.generatedId, "synthetic-ethene", etheneGraph());
    SaceGeneratedRecord const *recE = cat.record(refE.generatedId);
    SaceSimulationReadiness e = recE ? saceAssessSimulationReadiness(*recE) : SaceSimulationReadiness{};
    emit("unsupported_not_thermo_ready",
        recE && !e.gasThermoReady && !e.liquidThermoReady && !e.liquidGasEquilibriumReady, "");
    emit("unsupported_not_live_ready",
        recE && !e.liveGasReady && !e.liveLiquidReady && !e.liveLiquidGasPhaseChangeReady, "");

    emit("missing_mask_helper_works",
        !saceReadinessMissing(a, SaceSimulationRequirement::LiquidThermalConductivity)
            && !saceReadinessMissing(a, SaceSimulationRequirement::GasThermalConductivity)
            && !saceReadinessMissing(a, SaceSimulationRequirement::LiquidSurfaceTension)
            && !saceReadinessMissing(a, SaceSimulationRequirement::LiquidViscosity)
            && saceReadinessMissing(e, SaceSimulationRequirement::LiquidThermalConductivity)
            && saceReadinessMissing(e, SaceSimulationRequirement::GasThermalConductivity)
            && !saceReadinessMissing(a, SaceSimulationRequirement::MolarMass)
            && !saceReadinessMissing(a, SaceSimulationRequirement::MolecularGraph), "");
    char const *viscKey = saceSimulationRequirementKey(SaceSimulationRequirement::LiquidViscosity);
    char const *gasKKey = saceSimulationRequirementKey(SaceSimulationRequirement::GasThermalConductivity);
    emit("requirement_keys_stable",
        viscKey && std::strcmp(viscKey, "liquid_viscosity") == 0
            && gasKKey && std::strcmp(gasKKey, "gas_thermal_conductivity") == 0
            && saceSimulationRequirementKey(SaceSimulationRequirement::MolecularGraph)[0] != '\0', "");

    recA = cat.record(refA.generatedId);
    char sigAfter[kSaceSignatureCap]{};
    ChemicalIdentity viewAfter = recA ? cat.identityView(*recA) : ChemicalIdentity{};
    writeChemicalSignature(viewAfter, sigAfter, kSaceSignatureCap);
    emit("readiness_does_not_alter_signature",
        recA && std::strcmp(sigBefore, sigAfter) == 0, sigAfter);
    emit("readiness_does_not_alter_record_id",
        recA && recA->recordId == idBefore, "");
    emit("readiness_does_not_alter_display_ordinal",
        recA && recA->displayOrdinal == ordinalBefore, "");
    emit("generated_record_unspawnable", recA && !recA->spawnable && recB && !recB->spawnable, "");
    emit("substance_count_remains_12", SUBSTANCE_COUNT == 12, std::to_string(SUBSTANCE_COUNT));

    cat.clear();
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
