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
    if (dst.known && sacePropertySourceRank(dst.source) > sacePropertySourceRank(src.source))
        return false;
    dst = src;
    return true;
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

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
