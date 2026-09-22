#include "chemistry/SaceProperties.h"

#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceTypes.h"

#include <cmath>
#include <fstream>
#include <string>

bool saceAtomicMassGPerMol(AtomicNumber z, double &out) {
    out = 0.0;
    if (z == kAtomicHydrogen) { out = 1.008; return true; }
    if (z == kAtomicCarbon) { out = 12.011; return true; }
    if (z == kAtomicOxygen) { out = 15.999; return true; }
    return false;
}

bool saceDeriveMolarMass(ElementalComposition const &elemental, SaceScalarProperty &out) {
    out = saceUnknownScalarProperty();
    if (elemental.count == 0 || elemental.count > kMaxElementalSpecies)
        return false;
    ElementalComposition norm = elemental;
    if (normalizeElementalComposition(norm) != CompositionNormalizeResult::Ok)
        return false;
    if (norm.count == 0)
        return false;
    double sum = 0.0;
    for (int i = 0; i < norm.count; ++i) {
        AtomicNumber z = norm.entries[i].atomicNumber;
        uint16_t n = norm.entries[i].count;
        if (n == 0) continue;
        double mass = 0.0;
        if (!saceAtomicMassGPerMol(z, mass))
            return false;
        sum += mass * static_cast<double>(n);
    }
    if (!(sum > 0.0) || !std::isfinite(sum))
        return false;
    out.value = static_cast<float>(sum);
    out.known = true;
    out.source = SacePropertySource::IdentityDerived;
    out.confidence = SaceConfidence::High;
    return true;
}

bool saceAssignScalarProperty(SaceScalarProperty &dst, SaceScalarProperty const &src) {
    if (!src.known)
        return false;
    if (src.source == SacePropertySource::Unknown || src.confidence == SaceConfidence::Unknown)
        return false;
    if (!dst.known) {
        dst = src;
        return true;
    }
    int dstSource = sacePropertySourceRank(dst.source);
    int srcSource = sacePropertySourceRank(src.source);
    if (srcSource > dstSource) {
        dst = src;
        return true;
    }
    if (srcSource < dstSource)
        return false;
    int dstConf = saceConfidenceRank(dst.confidence);
    int srcConf = saceConfidenceRank(src.confidence);
    if (srcConf > dstConf) {
        dst = src;
        return true;
    }
    return false;
}

void initializeGeneratedProperties(SaceGeneratedProperties &props, ElementalComposition const &elemental) {
    props = {};
    SaceScalarProperty mm{};
    if (saceDeriveMolarMass(elemental, mm))
        saceAssignScalarProperty(props.molarMassGPerMol, mm);
}

namespace {

bool nearMass(double a, double b, double tol = 1.0e-3) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

bool deriveNear(ChemicalIdentity const &id, double expected, SaceScalarProperty &out) {
    return saceDeriveMolarMass(id.elemental, out)
        && out.known
        && out.source == SacePropertySource::IdentityDerived
        && out.confidence == SaceConfidence::High
        && nearMass(out.value, expected);
}

} // namespace

void runSacePropertyDiagnostics() {
    std::ofstream out(miscFile("sace_property_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    double hMass = 0, cMass = 0, oMass = 0, nMass = 0;
    bool hOk = saceAtomicMassGPerMol(kAtomicHydrogen, hMass);
    bool cOk = saceAtomicMassGPerMol(kAtomicCarbon, cMass);
    bool oOk = saceAtomicMassGPerMol(kAtomicOxygen, oMass);
    bool nOk = saceAtomicMassGPerMol(7, nMass);
    emit("h_atomic_mass_available", hOk && nearMass(hMass, 1.008), std::to_string(hMass));
    emit("c_atomic_mass_available", cOk && nearMass(cMass, 12.011), std::to_string(cMass));
    emit("o_atomic_mass_available", oOk && nearMass(oMass, 15.999), std::to_string(oMass));
    emit("unsupported_element_mass_unavailable", !nOk && nMass == 0.0, "");

    SaceScalarProperty mmH2{}, mmO2{}, mmH2O{}, mmC{}, mmCO2{}, mmC2{};
    bool h2Ok = deriveNear(saceBuiltinHydrogenIdentity(),
        substanceDef(SUBSTANCE_HYDROGEN).chemical.molarMass, mmH2);
    bool o2Ok = deriveNear(saceBuiltinOxygenIdentity(),
        substanceDef(SUBSTANCE_OXYGEN).chemical.molarMass, mmO2);
    bool h2oOk = deriveNear(saceBuiltinWaterIdentity(),
        substanceDef(SUBSTANCE_WATER).chemical.molarMass, mmH2O);
    bool cOkMass = deriveNear(saceBuiltinCarbonIdentity(),
        substanceDef(SUBSTANCE_CARBON).chemical.molarMass, mmC);
    bool co2Ok = deriveNear(saceBuiltinCarbonDioxideIdentity(),
        substanceDef(SUBSTANCE_CARBON_DIOXIDE).chemical.molarMass, mmCO2);
    emit("h2_matches_registry", h2Ok, std::to_string(mmH2.value));
    emit("o2_matches_registry", o2Ok, std::to_string(mmO2.value));
    emit("h2o_matches_registry", h2oOk, std::to_string(mmH2O.value));
    emit("c_matches_registry", cOkMass, std::to_string(mmC.value));
    emit("co2_matches_registry", co2Ok, std::to_string(mmCO2.value));

    ElementCount c2h6o[] = {
        {kAtomicCarbon, 2}, {kAtomicHydrogen, 6}, {kAtomicOxygen, 1}
    };
    ChemicalIdentity ethanolish = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-a", c2h6o, 3);
    bool c2Ok = deriveNear(ethanolish, 46.069, mmC2);
    emit("c2h6o_approximately_46_069", c2Ok, std::to_string(mmC2.value));

    ElementCount withN[] = {{kAtomicCarbon, 1}, {7, 1}};
    ElementalComposition unknownEl{};
    unknownEl.entries[0] = withN[0];
    unknownEl.entries[1] = withN[1];
    unknownEl.count = 2;
    SaceScalarProperty unknownMm{};
    bool derivedUnknown = saceDeriveMolarMass(unknownEl, unknownMm);
    emit("unknown_element_composition_unknown_mass",
        !derivedUnknown && !unknownMm.known
            && unknownMm.source == SacePropertySource::Unknown
            && unknownMm.confidence == SaceConfidence::Unknown, "");

    ElementalComposition shuffled{};
    shuffled.entries[0] = {kAtomicOxygen, 1};
    shuffled.entries[1] = {kAtomicHydrogen, 6};
    shuffled.entries[2] = {kAtomicCarbon, 2};
    shuffled.count = 3;
    SaceScalarProperty mmOrder{};
    bool orderOk = saceDeriveMolarMass(ethanolish.elemental, mmC2)
        && saceDeriveMolarMass(shuffled, mmOrder)
        && mmC2.known && mmOrder.known && nearMass(mmC2.value, mmOrder.value);
    emit("identity_order_does_not_affect_mass", orderOk, std::to_string(mmOrder.value));

    ChemicalIdentity formulaVariant = ethanolish;
    formulaVariant.formula = "EtOH";
    SaceScalarProperty mmFormula{};
    emit("formula_text_does_not_affect_mass",
        saceDeriveMolarMass(formulaVariant.elemental, mmFormula)
            && mmC2.known && mmFormula.known && nearMass(mmC2.value, mmFormula.value), "");

    ChemicalIdentity structureB = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H6O", "synthetic-structure-b", c2h6o, 3);
    SaceScalarProperty mmB{};
    emit("structure_key_does_not_alter_composition_mass",
        saceDeriveMolarMass(structureB.elemental, mmB)
            && mmC2.known && mmB.known && nearMass(mmC2.value, mmB.value)
            && !chemicalIdentitiesEquivalent(ethanolish, structureB), "");

    emit("builtin_molar_mass_not_overwritten",
        nearMass(substanceDef(SUBSTANCE_WATER).chemical.molarMass, 18.015)
            && nearMass(substanceDef(SUBSTANCE_CARBON_DIOXIDE).chemical.molarMass, 44.0095), "");

    auto makeKnown = [](float value, SacePropertySource source, SaceConfidence conf) {
        SaceScalarProperty p{};
        p.value = value;
        p.known = true;
        p.source = source;
        p.confidence = conf;
        return p;
    };

    SaceScalarProperty unknown{};
    SaceScalarProperty fallbackLow = makeKnown(1.0f, SacePropertySource::Fallback, SaceConfidence::Low);
    bool acceptedFallback = saceAssignScalarProperty(unknown, fallbackLow);
    emit("unknown_dest_accepts_fallback_low",
        acceptedFallback && unknown.known
            && unknown.source == SacePropertySource::Fallback
            && unknown.confidence == SaceConfidence::Low
            && unknown.value == 1.0f, "");

    SaceScalarProperty referenceHigh = makeKnown(18.015f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty fallbackHigh = makeKnown(99.0f, SacePropertySource::Fallback, SaceConfidence::High);
    bool rejectedFallback = !saceAssignScalarProperty(referenceHigh, fallbackHigh)
        && referenceHigh.source == SacePropertySource::Reference
        && referenceHigh.value == 18.015f;
    emit("reference_high_not_overwritten_by_fallback", rejectedFallback, "");

    SaceScalarProperty structuralHigh = makeKnown(10.0f, SacePropertySource::StructuralEstimate, SaceConfidence::High);
    SaceScalarProperty structuralLow = makeKnown(11.0f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_high_not_overwritten_by_low",
        !saceAssignScalarProperty(structuralHigh, structuralLow)
            && structuralHigh.confidence == SaceConfidence::High
            && structuralHigh.value == 10.0f, "");

    SaceScalarProperty structuralLowDst = makeKnown(10.0f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty structuralHighSrc = makeKnown(12.0f, SacePropertySource::StructuralEstimate, SaceConfidence::High);
    emit("structural_low_replaced_by_high",
        saceAssignScalarProperty(structuralLowDst, structuralHighSrc)
            && structuralLowDst.confidence == SaceConfidence::High
            && structuralLowDst.value == 12.0f, "");

    SaceScalarProperty empiricalMed = makeKnown(5.0f, SacePropertySource::EmpiricalEstimate, SaceConfidence::Medium);
    SaceScalarProperty structuralMed = makeKnown(6.0f, SacePropertySource::StructuralEstimate, SaceConfidence::Medium);
    emit("equal_rank_medium_tie_preserves_existing",
        !saceAssignScalarProperty(empiricalMed, structuralMed)
            && empiricalMed.source == SacePropertySource::EmpiricalEstimate
            && empiricalMed.value == 5.0f, "");

    SaceScalarProperty dstOk = makeKnown(1.0f, SacePropertySource::Fallback, SaceConfidence::Low);
    SaceScalarProperty srcUnknownSource = makeKnown(2.0f, SacePropertySource::Unknown, SaceConfidence::High);
    emit("known_unknown_source_rejected",
        !saceAssignScalarProperty(dstOk, srcUnknownSource) && dstOk.value == 1.0f, "");

    SaceScalarProperty srcUnknownConf = makeKnown(3.0f, SacePropertySource::Fallback, SaceConfidence::Unknown);
    emit("known_unknown_confidence_rejected",
        !saceAssignScalarProperty(dstOk, srcUnknownConf) && dstOk.value == 1.0f, "");

    SaceScalarProperty waterMm{};
    bool waterDerived = saceDeriveMolarMass(saceBuiltinWaterIdentity().elemental, waterMm);
    SaceScalarProperty waterFallback = makeKnown(1.0f, SacePropertySource::Fallback, SaceConfidence::High);
    emit("existing_molar_mass_remains_identity_derived_high",
        waterDerived
            && waterMm.source == SacePropertySource::IdentityDerived
            && waterMm.confidence == SaceConfidence::High
            && !saceAssignScalarProperty(waterMm, waterFallback)
            && waterMm.source == SacePropertySource::IdentityDerived
            && waterMm.confidence == SaceConfidence::High, "");

    SaceGeneratedProperties blankProps{};
    emit("normal_boiling_point_defaults_unknown",
        !blankProps.normalBoilingPointK.known
            && blankProps.normalBoilingPointK.source == SacePropertySource::Unknown
            && blankProps.normalBoilingPointK.confidence == SaceConfidence::Unknown, "");

    SaceScalarProperty refTb = makeKnown(351.44f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackTb = makeKnown(337.54f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_tb",
        !saceAssignScalarProperty(refTb, jobackTb)
            && refTb.source == SacePropertySource::Reference
            && refTb.confidence == SaceConfidence::High
            && nearMass(refTb.value, 351.44), "");

    SaceScalarProperty jobackDst = makeKnown(337.54f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refSrc = makeKnown(351.44f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_tb",
        saceAssignScalarProperty(jobackDst, refSrc)
            && jobackDst.source == SacePropertySource::Reference
            && jobackDst.confidence == SaceConfidence::High
            && nearMass(jobackDst.value, 351.44), "");

    SaceGeneratedProperties blankCrit{};
    emit("critical_temperature_defaults_unknown",
        !blankCrit.criticalTemperatureK.known
            && blankCrit.criticalTemperatureK.source == SacePropertySource::Unknown, "");
    emit("critical_pressure_defaults_unknown",
        !blankCrit.criticalPressurePa.known
            && blankCrit.criticalPressurePa.source == SacePropertySource::Unknown, "");
    emit("critical_molar_volume_defaults_unknown",
        !blankCrit.criticalMolarVolumeM3PerMol.known
            && blankCrit.criticalMolarVolumeM3PerMol.source == SacePropertySource::Unknown, "");

    SaceScalarProperty refTc = makeKnown(513.9f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackTc = makeKnown(499.407f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_tc",
        !saceAssignScalarProperty(refTc, jobackTc)
            && refTc.source == SacePropertySource::Reference
            && refTc.confidence == SaceConfidence::High
            && nearMass(refTc.value, 513.9), "");

    SaceScalarProperty jobackTcDst = makeKnown(499.407f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refTcSrc = makeKnown(513.9f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_tc",
        saceAssignScalarProperty(jobackTcDst, refTcSrc)
            && jobackTcDst.source == SacePropertySource::Reference
            && jobackTcDst.confidence == SaceConfidence::High
            && nearMass(jobackTcDst.value, 513.9), "");

    SaceGeneratedProperties blankOmega{};
    emit("acentric_factor_defaults_unknown",
        !blankOmega.acentricFactor.known
            && blankOmega.acentricFactor.source == SacePropertySource::Unknown, "");
    SaceScalarProperty refOmega = makeKnown(0.644f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackOmega = makeKnown(0.556081f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_omega",
        !saceAssignScalarProperty(refOmega, jobackOmega)
            && refOmega.source == SacePropertySource::Reference
            && nearMass(refOmega.value, 0.644), "");
    SaceScalarProperty jobackOmegaDst = makeKnown(0.556081f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refOmegaSrc = makeKnown(0.644f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_omega",
        saceAssignScalarProperty(jobackOmegaDst, refOmegaSrc)
            && jobackOmegaDst.source == SacePropertySource::Reference
            && nearMass(jobackOmegaDst.value, 0.644), "");

    SaceGeneratedProperties blankHvap{};
    emit("hvap_tb_defaults_unknown",
        !blankHvap.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && blankHvap.enthalpyVaporizationAtNormalBoilingJPerMol.source
                == SacePropertySource::Unknown, "");
    SaceScalarProperty refHvap = makeKnown(38600.0f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackHvap = makeKnown(36725.0f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_hvap",
        !saceAssignScalarProperty(refHvap, jobackHvap)
            && refHvap.source == SacePropertySource::Reference
            && nearMass(refHvap.value, 38600.0), "");
    SaceScalarProperty jobackHvapDst = makeKnown(36725.0f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refHvapSrc = makeKnown(38600.0f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_hvap",
        saceAssignScalarProperty(jobackHvapDst, refHvapSrc)
            && jobackHvapDst.source == SacePropertySource::Reference
            && nearMass(jobackHvapDst.value, 38600.0), "");

    SaceGeneratedProperties blankCp{};
    emit("ideal_gas_cp_298_defaults_unknown",
        !blankCp.idealGasHeatCapacityAt298KJPerMolK.known
            && blankCp.idealGasHeatCapacityAt298KJPerMolK.source == SacePropertySource::Unknown, "");
    emit("liquid_cp_298_defaults_unknown",
        !blankCp.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && blankCp.saturatedLiquidHeatCapacityAt298KJPerMolK.source == SacePropertySource::Unknown, "");

    SaceScalarProperty refGasCp = makeKnown(65.0f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackGasCp = makeKnown(64.6209f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_gas_cp",
        !saceAssignScalarProperty(refGasCp, jobackGasCp)
            && refGasCp.source == SacePropertySource::Reference
            && refGasCp.confidence == SaceConfidence::High
            && nearMass(refGasCp.value, 65.0), "");
    SaceScalarProperty jobackGasCpDst = makeKnown(64.6209f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refGasCpSrc = makeKnown(65.0f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_gas_cp",
        saceAssignScalarProperty(jobackGasCpDst, refGasCpSrc)
            && jobackGasCpDst.source == SacePropertySource::Reference
            && jobackGasCpDst.confidence == SaceConfidence::High
            && nearMass(jobackGasCpDst.value, 65.0), "");

    SaceScalarProperty refLiqCp = makeKnown(112.0f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackLiqCp = makeKnown(148.729f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_liquid_cp",
        !saceAssignScalarProperty(refLiqCp, jobackLiqCp)
            && refLiqCp.source == SacePropertySource::Reference
            && refLiqCp.confidence == SaceConfidence::High
            && nearMass(refLiqCp.value, 112.0), "");
    SaceScalarProperty jobackLiqCpDst = makeKnown(148.729f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refLiqCpSrc = makeKnown(112.0f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_liquid_cp",
        saceAssignScalarProperty(jobackLiqCpDst, refLiqCpSrc)
            && jobackLiqCpDst.source == SacePropertySource::Reference
            && jobackLiqCpDst.confidence == SaceConfidence::High
            && nearMass(jobackLiqCpDst.value, 112.0), "");

    SaceGeneratedProperties blankRho{};
    emit("saturated_liquid_density_298_defaults_unknown",
        !blankRho.saturatedLiquidDensityAt298KKgPerM3.known
            && blankRho.saturatedLiquidDensityAt298KKgPerM3.source == SacePropertySource::Unknown, "");
    SaceScalarProperty refRho = makeKnown(789.0f, SacePropertySource::Reference, SaceConfidence::High);
    SaceScalarProperty jobackRho = makeKnown(821.662f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    emit("structural_estimate_low_does_not_overwrite_reference_high_density",
        !saceAssignScalarProperty(refRho, jobackRho)
            && refRho.source == SacePropertySource::Reference
            && refRho.confidence == SaceConfidence::High
            && nearMass(refRho.value, 789.0), "");
    SaceScalarProperty jobackRhoDst = makeKnown(821.662f, SacePropertySource::StructuralEstimate, SaceConfidence::Low);
    SaceScalarProperty refRhoSrc = makeKnown(789.0f, SacePropertySource::Reference, SaceConfidence::High);
    emit("reference_high_overwrites_structural_estimate_low_density",
        saceAssignScalarProperty(jobackRhoDst, refRhoSrc)
            && jobackRhoDst.source == SacePropertySource::Reference
            && jobackRhoDst.confidence == SaceConfidence::High
            && nearMass(jobackRhoDst.value, 789.0), "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
