#include "chemistry/SaceEstimation.h"

#include "chemistry/SaceCatalog.h"
#include "chemistry/SaceFunctional.h"
#include "chemistry/SaceIdentity.h"
#include "fluid/DiagOutput.h"
#include "substance/SubstanceProperties.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace {

bool isH(AtomicNumber z) { return z == kAtomicHydrogen; }
bool isC(AtomicNumber z) { return z == kAtomicCarbon; }
bool isO(AtomicNumber z) { return z == kAtomicOxygen; }
bool isHCO(AtomicNumber z) { return isH(z) || isC(z) || isO(z); }

void addAtom(SaceMolecularGraph &g, AtomicNumber z, int8_t charge = 0) {
    g.atoms.push_back({z, charge});
}
void addBond(SaceMolecularGraph &g, uint16_t a, uint16_t b, SaceBondOrder order = SaceBondOrder::Single) {
    g.bonds.push_back({a, b, order});
}
void addHs(SaceMolecularGraph &g, uint16_t center, int nH) {
    for (int i = 0; i < nH; ++i) {
        uint16_t h = static_cast<uint16_t>(g.atoms.size());
        addAtom(g, kAtomicHydrogen);
        addBond(g, center, h);
    }
}

SaceMolecularGraph methaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addHs(g, 0, 4);
    return g;
}

SaceMolecularGraph isobutaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon); // 0 CH
    addAtom(g, kAtomicCarbon); // 1 CH3
    addAtom(g, kAtomicCarbon); // 2 CH3
    addAtom(g, kAtomicCarbon); // 3 CH3
    addBond(g, 0, 1);
    addBond(g, 0, 2);
    addBond(g, 0, 3);
    addHs(g, 0, 1);
    addHs(g, 1, 3);
    addHs(g, 2, 3);
    addHs(g, 3, 3);
    return g;
}

SaceMolecularGraph neopentaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon); // 0 C
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1);
    addBond(g, 0, 2);
    addBond(g, 0, 3);
    addBond(g, 0, 4);
    addHs(g, 1, 3);
    addHs(g, 2, 3);
    addHs(g, 3, 3);
    addHs(g, 4, 3);
    return g;
}

SaceMolecularGraph cyclopropaneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1);
    addBond(g, 1, 2);
    addBond(g, 2, 0);
    addHs(g, 0, 2);
    addHs(g, 1, 2);
    addHs(g, 2, 2);
    return g;
}

SaceMolecularGraph etheneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1, SaceBondOrder::Double);
    addHs(g, 0, 2);
    addHs(g, 1, 2);
    return g;
}

SaceMolecularGraph ethyneGraph() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicCarbon);
    addAtom(g, kAtomicCarbon);
    addBond(g, 0, 1, SaceBondOrder::Triple);
    addHs(g, 0, 1);
    addHs(g, 1, 1);
    return g;
}

SaceMolecularGraph chargedWaterLike() {
    SaceMolecularGraph g{};
    addAtom(g, kAtomicOxygen, 1);
    addHs(g, 0, 2);
    return g;
}

SaceMolecularGraph nitrogenGraph() {
    SaceMolecularGraph g{};
    addAtom(g, 7);
    addHs(g, 0, 3);
    return g;
}

SaceMolecularGraph permuteGraph(SaceMolecularGraph const &src) {
    SaceMolecularGraph g{};
    int n = atomCount(src);
    if (n <= 0) return g;
    std::vector<int> map(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        map[static_cast<size_t>(i)] = n - 1 - i;
    g.atoms.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        g.atoms[static_cast<size_t>(map[static_cast<size_t>(i)])] = src.atoms[static_cast<size_t>(i)];
    for (SaceBond const &b : src.bonds) {
        SaceBond nb = b;
        nb.atomA = static_cast<uint16_t>(map[b.atomA]);
        nb.atomB = static_cast<uint16_t>(map[b.atomB]);
        g.bonds.push_back(nb);
    }
    return g;
}

SaceMolecularGraph reorderBonds(SaceMolecularGraph g) {
    if (g.bonds.size() >= 2)
        std::swap(g.bonds[0], g.bonds[g.bonds.size() - 1]);
    return g;
}

SaceMolecularGraph reverseEnds(SaceMolecularGraph g) {
    for (SaceBond &b : g.bonds) {
        uint16_t t = b.atomA;
        b.atomA = b.atomB;
        b.atomB = t;
    }
    return g;
}

bool nearK(double a, double b, double tol = 0.02) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tol;
}

} // namespace

bool saceJobackGroupCountsEqual(SaceJobackGroupCounts const &a, SaceJobackGroupCounts const &b) {
    return a.carbonCH3 == b.carbonCH3 && a.carbonCH2 == b.carbonCH2
        && a.carbonCH == b.carbonCH && a.carbonC == b.carbonC
        && a.hydroxylAlcohol == b.hydroxylAlcohol && a.etherNonRing == b.etherNonRing;
}

SaceJobackFragmentationResult saceFragmentJobackGroups(
    SaceMolecularGraph const &graph, SaceJobackGroupCounts &counts)
{
    counts = {};
    if (validateMolecularGraph(graph, true) != SaceGraphValidation::Ok)
        return SaceJobackFragmentationResult::InvalidGraph;

    for (SaceAtom const &atom : graph.atoms) {
        if (!isHCO(atom.atomicNumber))
            return SaceJobackFragmentationResult::UnsupportedElement;
    }
    if (molecularGraphNetFormalCharge(graph) != 0)
        return SaceJobackFragmentationResult::Charged;

    SaceMolecularDescriptors desc{};
    if (!deriveMolecularDescriptors(graph, desc))
        return SaceJobackFragmentationResult::InvalidGraph;
    if (desc.cycleRank != 0)
        return SaceJobackFragmentationResult::Cyclic;

    for (SaceBond const &b : graph.bonds) {
        if (b.order != SaceBondOrder::Single)
            return SaceJobackFragmentationResult::Unsaturated;
    }

    int n = atomCount(graph);
    for (int i = 0; i < n; ++i) {
        if (!isH(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 1)
            return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
        uint16_t idx = static_cast<uint16_t>(i);
        bool okH = false;
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            if (b.order != SaceBondOrder::Single)
                return SaceJobackFragmentationResult::Unsaturated;
            AtomicNumber z = graph.atoms[static_cast<size_t>(other)].atomicNumber;
            if (!isC(z) && !isO(z))
                return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
            okH = true;
        }
        if (!okH)
            return SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment;
    }

    std::vector<char> covered(static_cast<size_t>(n), 0);

    auto cover = [&](int idx) {
        if (idx < 0 || idx >= n) return false;
        if (covered[static_cast<size_t>(idx)]) return false;
        covered[static_cast<size_t>(idx)] = 1;
        return true;
    };

    auto neighborZs = [&](int atom, int &hN, int &cN, int &oN, int &otherN) {
        hN = cN = oN = otherN = 0;
        uint16_t idx = static_cast<uint16_t>(atom);
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            AtomicNumber z = graph.atoms[static_cast<size_t>(other)].atomicNumber;
            if (isH(z)) ++hN;
            else if (isC(z)) ++cN;
            else if (isO(z)) ++oN;
            else ++otherN;
        }
    };

    for (int i = 0; i < n; ++i) {
        if (!isC(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        int hN = 0, cN = 0, oN = 0, otherN = 0;
        neighborZs(i, hN, cN, oN, otherN);
        if (otherN != 0 || hN + cN + oN != 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (hN == 4)
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (hN == 3) ++counts.carbonCH3;
        else if (hN == 2) ++counts.carbonCH2;
        else if (hN == 1) ++counts.carbonCH;
        else if (hN == 0) ++counts.carbonC;
        else
            return SaceJobackFragmentationResult::UnsupportedCarbonEnvironment;
        if (!cover(i))
            return SaceJobackFragmentationResult::IncompleteCoverage;
        uint16_t idx = static_cast<uint16_t>(i);
        for (SaceBond const &b : graph.bonds) {
            int other = -1;
            if (b.atomA == idx) other = b.atomB;
            else if (b.atomB == idx) other = b.atomA;
            if (other < 0) continue;
            if (!isH(graph.atoms[static_cast<size_t>(other)].atomicNumber))
                continue;
            if (!cover(other))
                return SaceJobackFragmentationResult::IncompleteCoverage;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!isO(graph.atoms[static_cast<size_t>(i)].atomicNumber))
            continue;
        if (graph.atoms[static_cast<size_t>(i)].formalCharge != 0)
            return SaceJobackFragmentationResult::Charged;
        if (atomDegree(graph, i) != 2)
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        int hN = 0, cN = 0, oN = 0, otherN = 0;
        neighborZs(i, hN, cN, oN, otherN);
        if (otherN != 0 || oN != 0)
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        if (hN == 1 && cN == 1) {
            ++counts.hydroxylAlcohol;
            if (!cover(i))
                return SaceJobackFragmentationResult::IncompleteCoverage;
            uint16_t idx = static_cast<uint16_t>(i);
            for (SaceBond const &b : graph.bonds) {
                int other = -1;
                if (b.atomA == idx) other = b.atomB;
                else if (b.atomB == idx) other = b.atomA;
                if (other < 0) continue;
                if (!isH(graph.atoms[static_cast<size_t>(other)].atomicNumber))
                    continue;
                if (!cover(other))
                    return SaceJobackFragmentationResult::IncompleteCoverage;
            }
        } else if (hN == 0 && cN == 2) {
            ++counts.etherNonRing;
            if (!cover(i))
                return SaceJobackFragmentationResult::IncompleteCoverage;
        } else {
            return SaceJobackFragmentationResult::UnsupportedOxygenEnvironment;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (!covered[static_cast<size_t>(i)])
            return SaceJobackFragmentationResult::IncompleteCoverage;
    }
    return SaceJobackFragmentationResult::Ok;
}

void saceResetJobackEstimateBundle(SaceJobackEstimateBundle &out) {
    out = {};
}

namespace {

SaceScalarProperty jobackLowScalar(double value) {
    SaceScalarProperty p = saceUnknownScalarProperty();
    if (!(value > 0.0) || !std::isfinite(value))
        return p;
    p.value = static_cast<float>(value);
    p.known = true;
    p.source = SacePropertySource::StructuralEstimate;
    p.confidence = SaceConfidence::Low;
    return p;
}

SaceScalarProperty jobackLowFiniteScalar(double value) {
    SaceScalarProperty p = saceUnknownScalarProperty();
    if (!std::isfinite(value))
        return p;
    p.value = static_cast<float>(value);
    p.known = true;
    p.source = SacePropertySource::StructuralEstimate;
    p.confidence = SaceConfidence::Low;
    return p;
}

bool estimateLeeKeslerOmega(double tb, double tc, double pc, double &omega) {
    omega = 0.0;
    if (!(tb > 0.0) || !(tc > 0.0) || !(pc > 0.0))
        return false;
    if (!std::isfinite(tb) || !std::isfinite(tc) || !std::isfinite(pc))
        return false;
    double tbr = tb / tc;
    if (!(tbr > 0.0) || !(tbr < 1.0) || !std::isfinite(tbr))
        return false;
    double lnTbr = std::log(tbr);
    double tbr6 = tbr * tbr * tbr * tbr * tbr * tbr;
    double num = std::log(kLeeKeslerAtmPa / pc)
        - 5.92714
        + 6.09648 / tbr
        + 1.28862 * lnTbr
        - 0.169347 * tbr6;
    double den = 15.2518
        - 15.6875 / tbr
        - 13.4721 * lnTbr
        + 0.43577 * tbr6;
    if (!std::isfinite(num) || !std::isfinite(den) || std::fabs(den) < 1.0e-14)
        return false;
    omega = num / den;
    return std::isfinite(omega);
}

} // namespace

bool saceLeeKeslerSaturationPressurePa(
    double temperatureK,
    double criticalTemperatureK,
    double criticalPressurePa,
    double acentricFactor,
    double &outPressurePa)
{
    outPressurePa = 0.0;
    if (!std::isfinite(temperatureK) || !std::isfinite(criticalTemperatureK)
        || !std::isfinite(criticalPressurePa) || !std::isfinite(acentricFactor))
        return false;
    if (!(temperatureK > 0.0) || !(criticalTemperatureK > 0.0) || !(criticalPressurePa > 0.0))
        return false;
    if (temperatureK > criticalTemperatureK)
        return false;
    if (temperatureK == criticalTemperatureK) {
        outPressurePa = criticalPressurePa;
        return true;
    }
    double tr = temperatureK / criticalTemperatureK;
    if (!(tr > 0.0) || !(tr < 1.0) || !std::isfinite(tr))
        return false;
    double lnTr = std::log(tr);
    double tr6 = tr * tr * tr * tr * tr * tr;
    double f0 = 5.92714 - 6.09648 / tr - 1.28862 * lnTr + 0.169347 * tr6;
    double f1 = 15.2518 - 15.6875 / tr - 13.4721 * lnTr + 0.43577 * tr6;
    if (!std::isfinite(f0) || !std::isfinite(f1))
        return false;
    double lnPr = f0 + acentricFactor * f1;
    if (!std::isfinite(lnPr))
        return false;
    double pr = std::exp(lnPr);
    if (!std::isfinite(pr) || !(pr > 0.0))
        return false;
    double psat = pr * criticalPressurePa;
    if (!std::isfinite(psat) || !(psat > 0.0))
        return false;
    if (psat > criticalPressurePa)
        return false;
    outPressurePa = psat;
    return true;
}

bool saceBuildLeeKeslerVaporModelFromJoback(
    SaceJobackEstimateBundle const &bundle,
    SaceLeeKeslerVaporModel &out)
{
    out = {};
    if (!bundle.normalBoilingPointK.known || !bundle.criticalTemperatureK.known
        || !bundle.criticalPressurePa.known || !bundle.acentricFactor.known)
        return false;
    if (bundle.normalBoilingPointK.source == SacePropertySource::Unknown
        || bundle.criticalTemperatureK.source == SacePropertySource::Unknown
        || bundle.criticalPressurePa.source == SacePropertySource::Unknown
        || bundle.acentricFactor.source == SacePropertySource::Unknown)
        return false;
    out.normalBoilingPointK = bundle.normalBoilingPointK.value;
    out.criticalTemperatureK = bundle.criticalTemperatureK.value;
    out.criticalPressurePa = bundle.criticalPressurePa.value;
    out.acentricFactor = bundle.acentricFactor.value;
    if (!(out.normalBoilingPointK > 0.0) || !(out.criticalTemperatureK > 0.0)
        || !(out.criticalPressurePa > 0.0)
        || !std::isfinite(out.normalBoilingPointK)
        || !std::isfinite(out.criticalTemperatureK)
        || !std::isfinite(out.criticalPressurePa)
        || !std::isfinite(out.acentricFactor))
        return false;
    if (!(out.normalBoilingPointK < out.criticalTemperatureK))
        return false;
    out.valid = true;
    return true;
}

bool saceBuildWatsonVaporizationModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceWatsonVaporizationModel &out)
{
    out = {};
    if (!bundle.normalBoilingPointK.known || !bundle.criticalTemperatureK.known
        || !bundle.enthalpyVaporizationAtNormalBoilingJPerMol.known)
        return false;
    if (bundle.normalBoilingPointK.source == SacePropertySource::Unknown
        || bundle.criticalTemperatureK.source == SacePropertySource::Unknown
        || bundle.enthalpyVaporizationAtNormalBoilingJPerMol.source == SacePropertySource::Unknown)
        return false;
    ElementalComposition elemental{};
    if (!elementalCompositionFromGraph(graph, elemental))
        return false;
    SaceScalarProperty mm{};
    if (!saceDeriveMolarMass(elemental, mm) || !mm.known || !(mm.value > 0.0f))
        return false;
    out.referenceTemperatureK = bundle.normalBoilingPointK.value;
    out.criticalTemperatureK = bundle.criticalTemperatureK.value;
    out.referenceEnthalpyJPerMol = bundle.enthalpyVaporizationAtNormalBoilingJPerMol.value;
    out.molarMassKgPerMol = static_cast<double>(mm.value) * 1.0e-3;
    if (!(out.referenceTemperatureK > 0.0) || !(out.criticalTemperatureK > 0.0)
        || !(out.referenceEnthalpyJPerMol > 0.0) || !(out.molarMassKgPerMol > 0.0)
        || !std::isfinite(out.referenceTemperatureK)
        || !std::isfinite(out.criticalTemperatureK)
        || !std::isfinite(out.referenceEnthalpyJPerMol)
        || !std::isfinite(out.molarMassKgPerMol))
        return false;
    if (!(out.referenceTemperatureK < out.criticalTemperatureK))
        return false;
    out.valid = true;
    return true;
}

bool saceWatsonEnthalpyVaporizationJPerMol(
    SaceWatsonVaporizationModel const &model,
    double temperatureK,
    double &outJPerMol)
{
    outJPerMol = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !std::isfinite(model.criticalTemperatureK)
        || !std::isfinite(model.referenceTemperatureK)
        || !std::isfinite(model.referenceEnthalpyJPerMol))
        return false;
    if (!(temperatureK > 0.0) || !(model.criticalTemperatureK > 0.0)
        || !(model.referenceTemperatureK > 0.0)
        || !(model.referenceEnthalpyJPerMol > 0.0))
        return false;
    if (!(model.referenceTemperatureK < model.criticalTemperatureK))
        return false;
    if (temperatureK > model.criticalTemperatureK)
        return false;
    if (temperatureK == model.criticalTemperatureK) {
        outJPerMol = 0.0;
        return true;
    }
    double oneMinusTr = 1.0 - temperatureK / model.criticalTemperatureK;
    double oneMinusTbr = 1.0 - model.referenceTemperatureK / model.criticalTemperatureK;
    if (!(oneMinusTr > 0.0) || !(oneMinusTbr > 0.0) || !std::isfinite(oneMinusTr)
        || !std::isfinite(oneMinusTbr))
        return false;
    double ratio = oneMinusTr / oneMinusTbr;
    if (!(ratio >= 0.0) || !std::isfinite(ratio))
        return false;
    double hvap = model.referenceEnthalpyJPerMol * std::pow(ratio, kWatsonHvapExponent);
    if (!std::isfinite(hvap) || hvap < 0.0)
        return false;
    outJPerMol = hvap;
    return true;
}

bool saceWatsonLatentHeatVaporizationJPerKg(
    SaceWatsonVaporizationModel const &model,
    double temperatureK,
    double &outJPerKg)
{
    outJPerKg = 0.0;
    if (!(model.molarMassKgPerMol > 0.0) || !std::isfinite(model.molarMassKgPerMol))
        return false;
    double jPerMol = 0.0;
    if (!saceWatsonEnthalpyVaporizationJPerMol(model, temperatureK, jPerMol))
        return false;
    double jPerKg = jPerMol / model.molarMassKgPerMol;
    if (!std::isfinite(jPerKg) || jPerKg < 0.0)
        return false;
    outJPerKg = jPerKg;
    return true;
}

bool saceBuildJobackIdealGasCpModel(
    SaceJobackGroupCounts const &groups,
    SaceJobackIdealGasCpModel &out)
{
    out = {};
    double a = kJobackCpACH3 * groups.carbonCH3
        + kJobackCpACH2 * groups.carbonCH2
        + kJobackCpACH * groups.carbonCH
        + kJobackCpAC * groups.carbonC
        + kJobackCpAAlcoholOH * groups.hydroxylAlcohol
        + kJobackCpAEtherO * groups.etherNonRing;
    double b = kJobackCpBCH3 * groups.carbonCH3
        + kJobackCpBCH2 * groups.carbonCH2
        + kJobackCpBCH * groups.carbonCH
        + kJobackCpBC * groups.carbonC
        + kJobackCpBAlcoholOH * groups.hydroxylAlcohol
        + kJobackCpBEtherO * groups.etherNonRing;
    double c = kJobackCpCCH3 * groups.carbonCH3
        + kJobackCpCCH2 * groups.carbonCH2
        + kJobackCpCCH * groups.carbonCH
        + kJobackCpCC * groups.carbonC
        + kJobackCpCAlcoholOH * groups.hydroxylAlcohol
        + kJobackCpCEtherO * groups.etherNonRing;
    double d = kJobackCpDCH3 * groups.carbonCH3
        + kJobackCpDCH2 * groups.carbonCH2
        + kJobackCpDCH * groups.carbonCH
        + kJobackCpDC * groups.carbonC
        + kJobackCpDAlcoholOH * groups.hydroxylAlcohol
        + kJobackCpDEtherO * groups.etherNonRing;
    out.A = a + kJobackCpAIntercept;
    out.B = b + kJobackCpBIntercept;
    out.C = c + kJobackCpCIntercept;
    out.D = d + kJobackCpDIntercept;
    if (!std::isfinite(out.A) || !std::isfinite(out.B)
        || !std::isfinite(out.C) || !std::isfinite(out.D)) {
        out = {};
        return false;
    }
    out.valid = true;
    return true;
}

bool saceJobackIdealGasHeatCapacityJPerMolK(
    SaceJobackIdealGasCpModel const &model,
    double temperatureK,
    double &outJPerMolK)
{
    outJPerMolK = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !std::isfinite(model.A) || !std::isfinite(model.B)
        || !std::isfinite(model.C) || !std::isfinite(model.D))
        return false;
    if (temperatureK < kSaceHeatCapacityTMinK || temperatureK > kSaceHeatCapacityTMaxK)
        return false;
    double t = temperatureK;
    double t2 = t * t;
    double t3 = t2 * t;
    double cp = model.A + model.B * t + model.C * t2 + model.D * t3;
    if (!std::isfinite(cp) || !(cp > 0.0))
        return false;
    outJPerMolK = cp;
    return true;
}

bool saceMolarHeatCapacityToSpecificJPerKgK(
    double cpJPerMolK,
    double molarMassKgPerMol,
    double &outJPerKgK)
{
    outJPerKgK = 0.0;
    if (!std::isfinite(cpJPerMolK) || !std::isfinite(molarMassKgPerMol))
        return false;
    if (!(cpJPerMolK > 0.0) || !(molarMassKgPerMol > 0.0))
        return false;
    double spec = cpJPerMolK / molarMassKgPerMol;
    if (!std::isfinite(spec) || !(spec > 0.0))
        return false;
    outJPerKgK = spec;
    return true;
}

bool saceBuildRowlinsonPolingLiquidCpModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceRowlinsonPolingLiquidCpModel &out)
{
    out = {};
    if (!bundle.criticalTemperatureK.known || !bundle.acentricFactor.known)
        return false;
    if (!saceBuildJobackIdealGasCpModel(bundle.groups, out.idealGasCp))
        return false;
    ElementalComposition elemental{};
    if (!elementalCompositionFromGraph(graph, elemental))
        return false;
    SaceScalarProperty mm{};
    if (!saceDeriveMolarMass(elemental, mm) || !mm.known || !(mm.value > 0.0f))
        return false;
    out.criticalTemperatureK = bundle.criticalTemperatureK.value;
    out.acentricFactor = bundle.acentricFactor.value;
    out.molarMassKgPerMol = static_cast<double>(mm.value) * 1.0e-3;
    if (!(out.criticalTemperatureK > 0.0) || !std::isfinite(out.criticalTemperatureK)
        || !std::isfinite(out.acentricFactor) || !(out.molarMassKgPerMol > 0.0)
        || !std::isfinite(out.molarMassKgPerMol)) {
        out = {};
        return false;
    }
    out.valid = true;
    return true;
}

bool saceRowlinsonPolingLiquidHeatCapacityJPerMolK(
    SaceRowlinsonPolingLiquidCpModel const &model,
    double temperatureK,
    double &outJPerMolK)
{
    outJPerMolK = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !std::isfinite(model.criticalTemperatureK)
        || !std::isfinite(model.acentricFactor))
        return false;
    if (temperatureK < kSaceHeatCapacityTMinK || temperatureK > kSaceHeatCapacityTMaxK)
        return false;
    if (!(temperatureK < model.criticalTemperatureK) || !(model.criticalTemperatureK > 0.0))
        return false;
    double tr = temperatureK / model.criticalTemperatureK;
    if (!(tr > 0.0) || !std::isfinite(tr) || tr >= kSaceLiquidCpTrReject)
        return false;
    double x = 1.0 - tr;
    if (!(x > 0.0) || !std::isfinite(x))
        return false;
    double cpig = 0.0;
    if (!saceJobackIdealGasHeatCapacityJPerMolK(model.idealGasCp, temperatureK, cpig))
        return false;
    double r = static_cast<double>(UNIVERSAL_GAS_R_J_MOL_K);
    double term = 1.586 + 0.49 / x
        + model.acentricFactor * (4.2775 + 6.3 * std::cbrt(x) / tr + 0.4355 / x);
    double cpl = cpig + r * term;
    if (!std::isfinite(cpl) || !(cpl > 0.0))
        return false;
    outJPerMolK = cpl;
    return true;
}

bool saceBuildCostaldLiquidDensityModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceCostaldLiquidDensityModel &out)
{
    out = {};
    if (!bundle.criticalTemperatureK.known
        || !bundle.criticalMolarVolumeM3PerMol.known
        || !bundle.acentricFactor.known)
        return false;
    ElementalComposition elemental{};
    if (!elementalCompositionFromGraph(graph, elemental))
        return false;
    SaceScalarProperty mm{};
    if (!saceDeriveMolarMass(elemental, mm) || !mm.known || !(mm.value > 0.0f))
        return false;
    out.criticalTemperatureK = bundle.criticalTemperatureK.value;
    out.characteristicVolumeM3PerMol = bundle.criticalMolarVolumeM3PerMol.value;
    out.acentricFactor = bundle.acentricFactor.value;
    out.molarMassKgPerMol = static_cast<double>(mm.value) * 1.0e-3;
    if (!(out.criticalTemperatureK > 0.0) || !std::isfinite(out.criticalTemperatureK)
        || !(out.characteristicVolumeM3PerMol > 0.0)
        || !std::isfinite(out.characteristicVolumeM3PerMol)
        || !std::isfinite(out.acentricFactor)
        || !(out.molarMassKgPerMol > 0.0) || !std::isfinite(out.molarMassKgPerMol)) {
        out = {};
        return false;
    }
    out.valid = true;
    return true;
}

bool saceCostaldSaturatedLiquidMolarVolumeM3PerMol(
    SaceCostaldLiquidDensityModel const &model,
    double temperatureK,
    double &outM3PerMol)
{
    outM3PerMol = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !(temperatureK > 0.0))
        return false;
    if (!std::isfinite(model.criticalTemperatureK) || !(model.criticalTemperatureK > 0.0)
        || !std::isfinite(model.characteristicVolumeM3PerMol)
        || !(model.characteristicVolumeM3PerMol > 0.0)
        || !std::isfinite(model.acentricFactor))
        return false;
    double tr = temperatureK / model.criticalTemperatureK;
    if (!std::isfinite(tr) || !(tr > kSaceCostaldTrMinExclusive)
        || !(tr < kSaceCostaldTrMaxExclusive))
        return false;
    double tau = 1.0 - tr;
    if (!(tau > 0.0) || !std::isfinite(tau))
        return false;
    double tau13 = std::cbrt(tau);
    double tau23 = tau13 * tau13;
    double tau43 = tau * tau13;
    double v0 = 1.0 - 1.52816 * tau13 + 1.43907 * tau23 - 0.81446 * tau + 0.190454 * tau43;
    if (!std::isfinite(v0) || !(v0 > 0.0))
        return false;
    double tr2 = tr * tr;
    double tr3 = tr2 * tr;
    double vdeltaNum = -0.296123 + 0.386914 * tr - 0.0427258 * tr2 - 0.0480645 * tr3;
    double vdeltaDen = tr - 1.00001;
    if (!std::isfinite(vdeltaNum) || !std::isfinite(vdeltaDen) || vdeltaDen == 0.0)
        return false;
    double vdelta = vdeltaNum / vdeltaDen;
    if (!std::isfinite(vdelta))
        return false;
    double corr = 1.0 - model.acentricFactor * vdelta;
    if (!std::isfinite(corr) || !(corr > 0.0))
        return false;
    double vs = model.characteristicVolumeM3PerMol * v0 * corr;
    if (!std::isfinite(vs) || !(vs > 0.0))
        return false;
    outM3PerMol = vs;
    return true;
}

bool saceCostaldSaturatedLiquidDensityKgPerM3(
    SaceCostaldLiquidDensityModel const &model,
    double temperatureK,
    double &outKgPerM3)
{
    outKgPerM3 = 0.0;
    if (!std::isfinite(model.molarMassKgPerMol) || !(model.molarMassKgPerMol > 0.0))
        return false;
    double vs = 0.0;
    if (!saceCostaldSaturatedLiquidMolarVolumeM3PerMol(model, temperatureK, vs))
        return false;
    double rho = model.molarMassKgPerMol / vs;
    if (!std::isfinite(rho) || !(rho > 0.0))
        return false;
    outKgPerM3 = rho;
    return true;
}

bool saceBuildJobackLiquidViscosityModel(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceJobackLiquidViscosityModel &out)
{
    out = {};
    SaceJobackGroupCounts const &g = bundle.groups;
    int nGroups = g.carbonCH3 + g.carbonCH2 + g.carbonCH + g.carbonC
        + g.hydroxylAlcohol + g.etherNonRing;
    if (nGroups <= 0)
        return false;
    ElementalComposition elemental{};
    if (!elementalCompositionFromGraph(graph, elemental))
        return false;
    SaceScalarProperty mm{};
    if (!saceDeriveMolarMass(elemental, mm) || !mm.known || !(mm.value > 0.0f))
        return false;
    double aKelvin = saceJobackViscosityMuASum(g) - kJobackViscAOffset;
    double bDimensionless = saceJobackViscosityMuBSum(g) - kJobackViscBOffset;
    double molarMassGPerMol = static_cast<double>(mm.value);
    if (!std::isfinite(aKelvin) || !std::isfinite(bDimensionless)
        || !std::isfinite(molarMassGPerMol) || !(molarMassGPerMol > 0.0)) {
        out = {};
        return false;
    }
    out.aKelvin = aKelvin;
    out.bDimensionless = bDimensionless;
    out.molarMassGPerMol = molarMassGPerMol;
    out.valid = true;
    return true;
}

bool saceJobackLiquidViscosityPaS(
    SaceJobackLiquidViscosityModel const &model,
    double temperatureK,
    double &outPaS)
{
    outPaS = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !(temperatureK > 0.0))
        return false;
    if (!std::isfinite(model.aKelvin) || !std::isfinite(model.bDimensionless)
        || !std::isfinite(model.molarMassGPerMol) || !(model.molarMassGPerMol > 0.0))
        return false;
    double exponent = model.aKelvin / temperatureK + model.bDimensionless;
    if (!std::isfinite(exponent))
        return false;
    double e = std::exp(exponent);
    if (!std::isfinite(e))
        return false;
    double mu = model.molarMassGPerMol * e;
    if (!std::isfinite(mu) || !(mu > 0.0))
        return false;
    outPaS = mu;
    return true;
}

bool saceBuildSastriRaoSurfaceTensionModelFromJoback(
    SaceMolecularGraph const &graph,
    SaceJobackEstimateBundle const &bundle,
    SaceSastriRaoSurfaceTensionModel &out)
{
    out = {};
    if (!bundle.normalBoilingPointK.known || !bundle.criticalTemperatureK.known
        || !bundle.criticalPressurePa.known)
        return false;
    double tb = bundle.normalBoilingPointK.value;
    double tc = bundle.criticalTemperatureK.value;
    double pc = bundle.criticalPressurePa.value;
    if (!std::isfinite(tb) || !std::isfinite(tc) || !std::isfinite(pc))
        return false;
    if (!(tb > 0.0) || !(tc > 0.0) || !(pc > 0.0) || !(tb < tc))
        return false;
    SaceFunctionalProfile profile{};
    if (!deriveFunctionalProfile(graph, profile))
        return false;
    out.normalBoilingPointK = tb;
    out.criticalTemperatureK = tc;
    out.criticalPressurePa = pc;
    out.chemicalClass = (profile.supported && profile.hydroxylCount > 0)
        ? SaceSurfaceTensionClass::Alcohol
        : SaceSurfaceTensionClass::GeneralOrganic;
    out.valid = true;
    return true;
}

bool saceSastriRaoSurfaceTensionNPerM(
    SaceSastriRaoSurfaceTensionModel const &model,
    double temperatureK,
    double &outNPerM)
{
    outNPerM = 0.0;
    if (!model.valid)
        return false;
    if (!std::isfinite(temperatureK) || !(temperatureK > 0.0))
        return false;
    double tb = model.normalBoilingPointK;
    double tc = model.criticalTemperatureK;
    double pc = model.criticalPressurePa;
    if (!std::isfinite(tb) || !std::isfinite(tc) || !std::isfinite(pc))
        return false;
    if (!(tb > 0.0) || !(tc > 0.0) || !(pc > 0.0) || !(tb < tc))
        return false;
    if (temperatureK > tc)
        return false;
    if (!(temperatureK < tc)) {
        outNPerM = 0.0;
        return true;
    }
    double tbr = tb / tc;
    double tr = temperatureK / tc;
    if (!(tbr > 0.0) || !(tbr < 1.0) || !(tr > 0.0) || !(tr < 1.0))
        return false;
    double oneMinusTbr = 1.0 - tbr;
    double oneMinusTr = 1.0 - tr;
    if (!(oneMinusTbr > 0.0) || !(oneMinusTr > 0.0))
        return false;
    double k = kSastriRaoGeneralK;
    double x = kSastriRaoGeneralX;
    double y = kSastriRaoGeneralY;
    double z = kSastriRaoGeneralZ;
    double m = kSastriRaoGeneralM;
    if (model.chemicalClass == SaceSurfaceTensionClass::Alcohol) {
        k = kSastriRaoAlcoholK;
        x = kSastriRaoAlcoholX;
        y = kSastriRaoAlcoholY;
        z = kSastriRaoAlcoholZ;
        m = kSastriRaoAlcoholM;
    }
    double pcBar = pc * kSastriRaoPcPaToBar;
    if (!std::isfinite(pcBar) || !(pcBar > 0.0))
        return false;
    double pcTerm = std::pow(pcBar, x);
    double tbTerm = std::pow(tb, y);
    double tcTerm = std::pow(tc, z);
    double ratio = oneMinusTr / oneMinusTbr;
    double ratioTerm = std::pow(ratio, m);
    if (!std::isfinite(pcTerm) || !std::isfinite(tbTerm) || !std::isfinite(tcTerm)
        || !std::isfinite(ratio) || !std::isfinite(ratioTerm))
        return false;
    double sigmaMNPerM = k * pcTerm * tbTerm * tcTerm * ratioTerm;
    double sigmaNPerM = sigmaMNPerM * kSastriRaoMNPerMToNPerM;
    if (!std::isfinite(sigmaNPerM) || !(sigmaNPerM >= 0.0))
        return false;
    outNPerM = sigmaNPerM;
    return true;
}

bool saceEstimateJobackBundle(SaceMolecularGraph const &graph, SaceJobackEstimateBundle &out) {
    saceResetJobackEstimateBundle(out);
    SaceJobackFragmentationResult r = saceFragmentJobackGroups(graph, out.groups);
    if (r != SaceJobackFragmentationResult::Ok) {
        out.groups = {};
        return false;
    }

    double tb = kJobackTbInterceptK
        + kJobackTbCH3K * out.groups.carbonCH3
        + kJobackTbCH2K * out.groups.carbonCH2
        + kJobackTbCHK * out.groups.carbonCH
        + kJobackTbCK * out.groups.carbonC
        + kJobackTbAlcoholOHK * out.groups.hydroxylAlcohol
        + kJobackTbEtherOK * out.groups.etherNonRing;
    if (tb > 0.0 && std::isfinite(tb))
        out.normalBoilingPointK = jobackLowScalar(tb);

    double sTc = saceJobackTcContributionSum(out.groups);
    double denom = 0.584 + 0.965 * sTc - sTc * sTc;
    if (tb > 0.0 && std::isfinite(tb) && denom > 0.0 && std::isfinite(denom)) {
        double tc = tb / denom;
        if (tc > 0.0 && std::isfinite(tc))
            out.criticalTemperatureK = jobackLowScalar(tc);
    }

    double sPc = saceJobackPcContributionSum(out.groups);
    double nAtoms = static_cast<double>(graph.atoms.size());
    double base = 0.113 + 0.0032 * nAtoms - sPc;
    if (base > 0.0 && std::isfinite(base)) {
        double pcBar = 1.0 / (base * base);
        if (pcBar > 0.0 && std::isfinite(pcBar))
            out.criticalPressurePa = jobackLowScalar(pcBar * kJobackBarToPa);
    }

    double vcCm3 = kJobackVcInterceptCm3 + saceJobackVcContributionSumCm3PerMol(out.groups);
    if (vcCm3 > 0.0 && std::isfinite(vcCm3))
        out.criticalMolarVolumeM3PerMol = jobackLowScalar(vcCm3 * kJobackCm3ToM3);

    if (tb > 0.0 && std::isfinite(tb)
        && out.criticalTemperatureK.known && out.criticalPressurePa.known) {
        double tc = tb / denom;
        double pc = 0.0;
        if (base > 0.0 && std::isfinite(base)) {
            double pcBar = 1.0 / (base * base);
            if (pcBar > 0.0 && std::isfinite(pcBar))
                pc = pcBar * kJobackBarToPa;
        }
        double omega = 0.0;
        if (pc > 0.0 && estimateLeeKeslerOmega(tb, tc, pc, omega))
            out.acentricFactor = jobackLowFiniteScalar(omega);
    }

    double sHvapKJ = saceJobackHvapContributionSumKJPerMol(out.groups);
    double hvapKJ = kJobackHvapInterceptKJPerMol + sHvapKJ;
    if (hvapKJ > 0.0 && std::isfinite(hvapKJ))
        out.enthalpyVaporizationAtNormalBoilingJPerMol = jobackLowScalar(hvapKJ * kJobackKJToJ);

    SaceJobackIdealGasCpModel igCp{};
    if (saceBuildJobackIdealGasCpModel(out.groups, igCp)) {
        double cpig = 0.0;
        if (saceJobackIdealGasHeatCapacityJPerMolK(igCp, kSaceHeatCapacityReferenceK, cpig))
            out.idealGasHeatCapacityAt298KJPerMolK = jobackLowScalar(cpig);
    }
    SaceRowlinsonPolingLiquidCpModel liqCp{};
    if (saceBuildRowlinsonPolingLiquidCpModelFromJoback(graph, out, liqCp)) {
        double cpl = 0.0;
        if (saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqCp, kSaceHeatCapacityReferenceK, cpl))
            out.saturatedLiquidHeatCapacityAt298KJPerMolK = jobackLowScalar(cpl);
    }
    SaceCostaldLiquidDensityModel dens{};
    bool densOk = saceBuildCostaldLiquidDensityModelFromJoback(graph, out, dens);
    if (densOk) {
        double rho = 0.0;
        if (saceCostaldSaturatedLiquidDensityKgPerM3(dens, kSaceHeatCapacityReferenceK, rho))
            out.saturatedLiquidDensityAt298KKgPerM3 = jobackLowScalar(rho);
    }
    SaceJobackLiquidViscosityModel visc{};
    if (saceBuildJobackLiquidViscosityModel(graph, out, visc)) {
        double mu = 0.0;
        double rhoGate = 0.0;
        if (saceJobackLiquidViscosityPaS(visc, kSaceHeatCapacityReferenceK, mu)
            && densOk
            && saceCostaldSaturatedLiquidDensityKgPerM3(dens, kSaceHeatCapacityReferenceK, rhoGate))
            out.liquidDynamicViscosityAt298KPaS = jobackLowScalar(mu);
    }
    SaceSastriRaoSurfaceTensionModel sigma{};
    if (saceBuildSastriRaoSurfaceTensionModelFromJoback(graph, out, sigma)
        && kSaceHeatCapacityReferenceK < sigma.criticalTemperatureK) {
        double st = 0.0;
        if (saceSastriRaoSurfaceTensionNPerM(sigma, kSaceHeatCapacityReferenceK, st))
            out.liquidSurfaceTensionAt298KNPerM = jobackLowScalar(st);
    }

    return true;
}

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out, SaceJobackGroupCounts &counts)
{
    SaceJobackEstimateBundle bundle{};
    bool ok = saceEstimateJobackBundle(graph, bundle);
    counts = bundle.groups;
    out = bundle.normalBoilingPointK;
    return ok && out.known;
}

bool saceEstimateJobackNormalBoilingPoint(
    SaceMolecularGraph const &graph, SaceScalarProperty &out)
{
    SaceJobackGroupCounts counts{};
    return saceEstimateJobackNormalBoilingPoint(graph, out, counts);
}

void runSaceEstimationDiagnostics() {
    std::ofstream out(miscFile("sace_estimation_diag.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    auto unsupported = [&](char const *name, SaceMolecularGraph const &g,
        SaceJobackFragmentationResult expect)
    {
        SaceJobackGroupCounts c{};
        SaceScalarProperty p = {};
        p.known = true;
        p.value = 99.0f;
        SaceJobackFragmentationResult r = saceFragmentJobackGroups(g, c);
        bool est = saceEstimateJobackNormalBoilingPoint(g, p);
        emit(name, r == expect && !est && !p.known
            && p.source == SacePropertySource::Unknown
            && p.confidence == SaceConfidence::Unknown, saceJobackFragmentationResultKey(r));
    };

    SaceMolecularGraph empty{};
    SaceScalarProperty junk{};
    junk.known = true;
    junk.value = 400.0f;
    junk.source = SacePropertySource::StructuralEstimate;
    emit("invalid_graph_unknown",
        !saceEstimateJobackNormalBoilingPoint(empty, junk) && !junk.known, "");

    SaceMolecularGraph water{}, h2{}, o2{}, co2{};
    saceBuiltinMolecularGraph(SUBSTANCE_WATER, water);
    saceBuiltinMolecularGraph(SUBSTANCE_HYDROGEN, h2);
    saceBuiltinMolecularGraph(SUBSTANCE_OXYGEN, o2);
    saceBuiltinMolecularGraph(SUBSTANCE_CARBON_DIOXIDE, co2);
    unsupported("water_unsupported", water, SaceJobackFragmentationResult::UnsupportedOxygenEnvironment);
    unsupported("h2_unsupported", h2, SaceJobackFragmentationResult::UnsupportedHydrogenEnvironment);
    unsupported("o2_unsupported", o2, SaceJobackFragmentationResult::Unsaturated);
    unsupported("co2_unsupported", co2, SaceJobackFragmentationResult::Unsaturated);
    unsupported("methane_unsupported", methaneGraph(), SaceJobackFragmentationResult::UnsupportedCarbonEnvironment);
    unsupported("charged_unsupported", chargedWaterLike(), SaceJobackFragmentationResult::Charged);
    unsupported("ring_unsupported", cyclopropaneGraph(), SaceJobackFragmentationResult::Cyclic);
    unsupported("nitrogen_unsupported", nitrogenGraph(), SaceJobackFragmentationResult::UnsupportedElement);
    unsupported("cc_double_unsupported", etheneGraph(), SaceJobackFragmentationResult::Unsaturated);
    unsupported("cc_triple_unsupported", ethyneGraph(), SaceJobackFragmentationResult::Unsaturated);

    SaceMolecularGraph gA{}, gB{};
    saceSyntheticC2H6OGraphA(gA);
    saceSyntheticC2H6OGraphB(gB);
    SaceJobackGroupCounts cA{}, cB{};
    SaceScalarProperty tbA{}, tbB{};
    emit("c2h6o_a_fragmentation_ok",
        saceEstimateJobackNormalBoilingPoint(gA, tbA, cA)
            && saceFragmentJobackGroups(gA, cA) == SaceJobackFragmentationResult::Ok, "");
    emit("a_ch3_count_1", cA.carbonCH3 == 1, std::to_string(cA.carbonCH3));
    emit("a_ch2_count_1", cA.carbonCH2 == 1, std::to_string(cA.carbonCH2));
    emit("a_oh_count_1", cA.hydroxylAlcohol == 1, std::to_string(cA.hydroxylAlcohol));
    emit("a_ether_count_0", cA.etherNonRing == 0, std::to_string(cA.etherNonRing));
    emit("a_tb_337_54", tbA.known && nearK(tbA.value, 337.54), std::to_string(tbA.value));
    emit("a_source_structural_estimate",
        tbA.source == SacePropertySource::StructuralEstimate, sacePropertySourceKey(tbA.source));
    emit("a_confidence_low",
        tbA.confidence == SaceConfidence::Low, saceConfidenceKey(tbA.confidence));

    emit("c2h6o_b_fragmentation_ok",
        saceEstimateJobackNormalBoilingPoint(gB, tbB, cB), "");
    emit("b_ch3_count_2", cB.carbonCH3 == 2, std::to_string(cB.carbonCH3));
    emit("b_oh_count_0", cB.hydroxylAlcohol == 0, std::to_string(cB.hydroxylAlcohol));
    emit("b_ether_count_1", cB.etherNonRing == 1, std::to_string(cB.etherNonRing));
    emit("b_tb_267_78", tbB.known && nearK(tbB.value, 267.78), std::to_string(tbB.value));
    emit("b_source_structural_estimate",
        tbB.source == SacePropertySource::StructuralEstimate, "");
    emit("b_confidence_low", tbB.confidence == SaceConfidence::Low, "");

    SaceScalarProperty mmA{}, mmB{};
    ElementalComposition elA{}, elB{};
    elementalCompositionFromGraph(gA, elA);
    elementalCompositionFromGraph(gB, elB);
    saceDeriveMolarMass(elA, mmA);
    saceDeriveMolarMass(elB, mmB);
    emit("a_b_molar_mass_equal",
        mmA.known && mmB.known && nearK(mmA.value, mmB.value, 1e-4), "");
    emit("a_b_boiling_points_differ",
        tbA.known && tbB.known && !nearK(tbA.value, tbB.value, 0.5), "");
    emit("a_tb_greater_than_b", tbA.known && tbB.known && tbA.value > tbB.value, "");
    emit("tb_difference_69_76",
        tbA.known && tbB.known && nearK(tbA.value - tbB.value, 69.76),
        std::to_string(tbA.value - tbB.value));

    emit("ethanol_topology_vs_literature", true,
        "joback=" + std::to_string(tbA.value) + " ref=351.44 err="
            + std::to_string(tbA.value - 351.44f));
    emit("dme_topology_vs_literature", true,
        "joback=" + std::to_string(tbB.value) + " ref=248.3 err="
            + std::to_string(tbB.value - 248.3f));

    SaceJobackGroupCounts isoC{};
    SaceScalarProperty isoTb{};
    bool isoOk = saceEstimateJobackNormalBoilingPoint(isobutaneGraph(), isoTb, isoC)
        && isoC.carbonCH == 1 && isoC.carbonCH3 == 3;
    emit("isobutane_exercises_ch", isoOk, std::to_string(isoC.carbonCH));
    SaceJobackGroupCounts neoC{};
    SaceScalarProperty neoTb{};
    bool neoOk = saceEstimateJobackNormalBoilingPoint(neopentaneGraph(), neoTb, neoC)
        && neoC.carbonC == 1 && neoC.carbonCH3 == 4;
    emit("neopentane_exercises_quaternary_c", neoOk, std::to_string(neoC.carbonC));

    SaceJobackGroupCounts permC{}, bondC{}, endC{};
    SaceScalarProperty permTb{}, bondTb{}, endTb{};
    emit("atom_permutation_preserves_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(permuteGraph(gA), permTb, permC)
            && saceJobackGroupCountsEqual(cA, permC) && nearK(permTb.value, tbA.value), "");
    emit("bond_reorder_preserves_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(reorderBonds(gA), bondTb, bondC)
            && saceJobackGroupCountsEqual(cA, bondC) && nearK(bondTb.value, tbA.value), "");
    emit("reversed_endpoints_preserve_groups_and_tb",
        saceEstimateJobackNormalBoilingPoint(reverseEnds(gA), endTb, endC)
            && saceJobackGroupCountsEqual(cA, endC) && nearK(endTb.value, tbA.value), "");

    SaceScalarProperty resetP{};
    resetP.known = true;
    resetP.value = 12.0f;
    resetP.source = SacePropertySource::StructuralEstimate;
    resetP.confidence = SaceConfidence::Low;
    emit("failed_estimation_resets_unknown",
        !saceEstimateJobackNormalBoilingPoint(empty, resetP) && !resetP.known
            && resetP.source == SacePropertySource::Unknown, "");

    char sig[kSaceSignatureCap]{};
    writeChemicalSignature(saceBuiltinWaterIdentity(), sig, kSaceSignatureCap);
    emit("joback_not_in_canonical_signature",
        std::strstr(sig, "337") == nullptr && std::strstr(sig, "Joback") == nullptr
            && std::strstr(sig, "boiling") == nullptr, sig);

    SaceScalarProperty refHigh{};
    refHigh.known = true;
    refHigh.value = 351.44f;
    refHigh.source = SacePropertySource::Reference;
    refHigh.confidence = SaceConfidence::High;
    SaceScalarProperty jobackLow = tbA;
    emit("higher_quality_not_overwritten_by_joback",
        !saceAssignScalarProperty(refHigh, jobackLow)
            && refHigh.source == SacePropertySource::Reference
            && nearK(refHigh.value, 351.44), "");

    SaceCatalog &cat = saceGeneratedCatalog();
    cat.clear();
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
    emit("catalog_a_joback_tb",
        recA && recA->hasMolecularGraph && recA->hasMolecularDescriptors
            && recA->hasFunctionalProfile
            && recA->properties.normalBoilingPointK.known
            && recA->properties.normalBoilingPointK.source == SacePropertySource::StructuralEstimate
            && recA->properties.normalBoilingPointK.confidence == SaceConfidence::Low
            && nearK(recA->properties.normalBoilingPointK.value, 337.54),
        recA ? std::to_string(recA->properties.normalBoilingPointK.value) : "");
    emit("catalog_b_joback_tb",
        recB && recB->properties.normalBoilingPointK.known
            && nearK(recB->properties.normalBoilingPointK.value, 267.78),
        recB ? std::to_string(recB->properties.normalBoilingPointK.value) : "");

    ElementCount c2h4[] = {{kAtomicCarbon, 2}, {kAtomicHydrogen, 4}};
    ChemicalIdentity idE = saceExactIdentity(ChemicalRepresentationKind::SmallMolecule,
        "C2H4", "synthetic-ethene", c2h4, 2);
    SaceSubstanceRef refE = cat.resolve(idE, true);
    SaceMolecularGraph ethene = etheneGraph();
    bool attached = cat.attachMolecularGraph(refE.generatedId, "synthetic-ethene", ethene);
    SaceGeneratedRecord const *recE = cat.record(refE.generatedId);
    emit("unsupported_joback_still_attaches_graph",
        attached && recE && recE->hasMolecularGraph && recE->hasFunctionalProfile
            && recE->hasMolecularDescriptors
            && !recE->properties.normalBoilingPointK.known
            && !recE->properties.criticalTemperatureK.known
            && !recE->properties.criticalPressurePa.known
            && !recE->properties.criticalMolarVolumeM3PerMol.known
            && !recE->properties.acentricFactor.known
            && !recE->properties.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && !recE->properties.idealGasHeatCapacityAt298KJPerMolK.known
            && !recE->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && !recE->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && !recE->properties.liquidDynamicViscosityAt298KPaS.known
            && !recE->properties.liquidSurfaceTensionAt298KNPerM.known, "");

    SaceJobackEstimateBundle bunA{}, bunB{};
    emit("c2h6o_a_bundle_ok", saceEstimateJobackBundle(gA, bunA), "");
    double aTcSum = saceJobackTcContributionSum(bunA.groups);
    double aPcSum = saceJobackPcContributionSum(bunA.groups);
    double aVcSum = saceJobackVcContributionSumCm3PerMol(bunA.groups);
    emit("a_tc_contribution_sum_0_1071", nearK(aTcSum, 0.1071, 1e-6), std::to_string(aTcSum));
    emit("a_pc_contribution_sum_0_0100", nearK(aPcSum, 0.0100, 1e-6), std::to_string(aPcSum));
    emit("a_vc_contribution_sum_149", nearK(aVcSum, 149.0, 1e-6), std::to_string(aVcSum));
    emit("a_tc_499_407", bunA.criticalTemperatureK.known
        && nearK(bunA.criticalTemperatureK.value, 499.407, 0.05),
        std::to_string(bunA.criticalTemperatureK.value));
    emit("a_pc_5_75664e6", bunA.criticalPressurePa.known
        && nearK(bunA.criticalPressurePa.value, 5.75664e6, 80.0),
        std::to_string(bunA.criticalPressurePa.value));
    emit("a_vc_1_665e-4", bunA.criticalMolarVolumeM3PerMol.known
        && nearK(bunA.criticalMolarVolumeM3PerMol.value, 1.665e-4, 1e-7),
        std::to_string(bunA.criticalMolarVolumeM3PerMol.value));
    emit("a_critical_structural_estimate_low",
        bunA.criticalTemperatureK.source == SacePropertySource::StructuralEstimate
            && bunA.criticalTemperatureK.confidence == SaceConfidence::Low
            && bunA.criticalPressurePa.source == SacePropertySource::StructuralEstimate
            && bunA.criticalPressurePa.confidence == SaceConfidence::Low
            && bunA.criticalMolarVolumeM3PerMol.source == SacePropertySource::StructuralEstimate
            && bunA.criticalMolarVolumeM3PerMol.confidence == SaceConfidence::Low, "");

    emit("c2h6o_b_bundle_ok", saceEstimateJobackBundle(gB, bunB), "");
    double bTcSum = saceJobackTcContributionSum(bunB.groups);
    double bPcSum = saceJobackPcContributionSum(bunB.groups);
    double bVcSum = saceJobackVcContributionSumCm3PerMol(bunB.groups);
    emit("b_tc_contribution_sum_0_0450", nearK(bTcSum, 0.0450, 1e-6), std::to_string(bTcSum));
    emit("b_pc_contribution_sum_minus_0_0009", nearK(bPcSum, -0.0009, 1e-6), std::to_string(bPcSum));
    emit("b_vc_contribution_sum_148", nearK(bVcSum, 148.0, 1e-6), std::to_string(bVcSum));
    emit("b_tc_428_174", bunB.criticalTemperatureK.known
        && nearK(bunB.criticalTemperatureK.value, 428.174, 0.05),
        std::to_string(bunB.criticalTemperatureK.value));
    emit("b_pc_4_91080e6", bunB.criticalPressurePa.known
        && nearK(bunB.criticalPressurePa.value, 4.91080e6, 80.0),
        std::to_string(bunB.criticalPressurePa.value));
    emit("b_vc_1_655e-4", bunB.criticalMolarVolumeM3PerMol.known
        && nearK(bunB.criticalMolarVolumeM3PerMol.value, 1.655e-4, 1e-7),
        std::to_string(bunB.criticalMolarVolumeM3PerMol.value));
    emit("b_critical_structural_estimate_low",
        bunB.criticalTemperatureK.source == SacePropertySource::StructuralEstimate
            && bunB.criticalTemperatureK.confidence == SaceConfidence::Low
            && bunB.criticalPressurePa.source == SacePropertySource::StructuralEstimate
            && bunB.criticalPressurePa.confidence == SaceConfidence::Low
            && bunB.criticalMolarVolumeM3PerMol.source == SacePropertySource::StructuralEstimate
            && bunB.criticalMolarVolumeM3PerMol.confidence == SaceConfidence::Low, "");

    emit("a_b_tc_differ", bunA.criticalTemperatureK.known && bunB.criticalTemperatureK.known
        && !nearK(bunA.criticalTemperatureK.value, bunB.criticalTemperatureK.value, 0.5), "");
    emit("a_b_pc_differ", bunA.criticalPressurePa.known && bunB.criticalPressurePa.known
        && !nearK(bunA.criticalPressurePa.value, bunB.criticalPressurePa.value, 100.0), "");
    emit("a_b_vc_differ", bunA.criticalMolarVolumeM3PerMol.known && bunB.criticalMolarVolumeM3PerMol.known
        && !nearK(bunA.criticalMolarVolumeM3PerMol.value, bunB.criticalMolarVolumeM3PerMol.value, 1e-8), "");

    auto bundleMatches = [&](SaceJobackEstimateBundle const &x, SaceJobackEstimateBundle const &y) {
        return saceJobackGroupCountsEqual(x.groups, y.groups)
            && x.normalBoilingPointK.known && y.normalBoilingPointK.known
            && x.criticalTemperatureK.known && y.criticalTemperatureK.known
            && x.criticalPressurePa.known && y.criticalPressurePa.known
            && x.criticalMolarVolumeM3PerMol.known && y.criticalMolarVolumeM3PerMol.known
            && x.acentricFactor.known && y.acentricFactor.known
            && x.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && y.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && nearK(x.normalBoilingPointK.value, y.normalBoilingPointK.value)
            && nearK(x.criticalTemperatureK.value, y.criticalTemperatureK.value, 0.05)
            && nearK(x.criticalPressurePa.value, y.criticalPressurePa.value, 80.0)
            && nearK(x.criticalMolarVolumeM3PerMol.value, y.criticalMolarVolumeM3PerMol.value, 1e-7)
            && nearK(x.acentricFactor.value, y.acentricFactor.value, 1e-5)
            && nearK(x.enthalpyVaporizationAtNormalBoilingJPerMol.value,
                y.enthalpyVaporizationAtNormalBoilingJPerMol.value, 1.0)
            && x.idealGasHeatCapacityAt298KJPerMolK.known
            && y.idealGasHeatCapacityAt298KJPerMolK.known
            && x.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && y.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && nearK(x.idealGasHeatCapacityAt298KJPerMolK.value,
                y.idealGasHeatCapacityAt298KJPerMolK.value, 0.02)
            && nearK(x.saturatedLiquidHeatCapacityAt298KJPerMolK.value,
                y.saturatedLiquidHeatCapacityAt298KJPerMolK.value, 0.05)
            && x.saturatedLiquidDensityAt298KKgPerM3.known
            && y.saturatedLiquidDensityAt298KKgPerM3.known
            && nearK(x.saturatedLiquidDensityAt298KKgPerM3.value,
                y.saturatedLiquidDensityAt298KKgPerM3.value, 0.05)
            && x.liquidDynamicViscosityAt298KPaS.known
            && y.liquidDynamicViscosityAt298KPaS.known
            && nearK(x.liquidDynamicViscosityAt298KPaS.value,
                y.liquidDynamicViscosityAt298KPaS.value, 1e-8)
            && x.liquidSurfaceTensionAt298KNPerM.known
            && y.liquidSurfaceTensionAt298KNPerM.known
            && nearK(x.liquidSurfaceTensionAt298KNPerM.value,
                y.liquidSurfaceTensionAt298KNPerM.value, 1e-8);
    };
    SaceJobackEstimateBundle permB{}, bondB{}, endB{};
    emit("atom_permutation_preserves_joback_bundle",
        saceEstimateJobackBundle(permuteGraph(gA), permB) && bundleMatches(bunA, permB), "");
    emit("bond_reorder_preserves_joback_bundle",
        saceEstimateJobackBundle(reorderBonds(gA), bondB) && bundleMatches(bunA, bondB), "");
    emit("reversed_endpoints_preserve_joback_bundle",
        saceEstimateJobackBundle(reverseEnds(gA), endB) && bundleMatches(bunA, endB), "");

    auto allJobackUnknown = [](SaceJobackEstimateBundle const &b) {
        return !b.normalBoilingPointK.known && !b.criticalTemperatureK.known
            && !b.criticalPressurePa.known && !b.criticalMolarVolumeM3PerMol.known
            && !b.acentricFactor.known
            && !b.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && !b.idealGasHeatCapacityAt298KJPerMolK.known
            && !b.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && !b.saturatedLiquidDensityAt298KKgPerM3.known
            && !b.liquidDynamicViscosityAt298KPaS.known
            && !b.liquidSurfaceTensionAt298KNPerM.known;
    };
    SaceJobackEstimateBundle waterB{}, methaneB{}, unsB{};
    emit("water_all_joback_unknown",
        !saceEstimateJobackBundle(water, waterB) && allJobackUnknown(waterB), "");
    emit("methane_all_joback_unknown",
        !saceEstimateJobackBundle(methaneGraph(), methaneB) && allJobackUnknown(methaneB), "");
    emit("unsaturated_all_joback_unknown",
        !saceEstimateJobackBundle(ethene, unsB) && allJobackUnknown(unsB), "");

    emit("critical_properties_not_in_canonical_signature",
        std::strstr(sig, "499") == nullptr && std::strstr(sig, "critical") == nullptr
            && std::strstr(sig, "Joback") == nullptr, sig);

    emit("catalog_a_joback_critical_and_unspawnable",
        recA && !recA->spawnable
            && recA->properties.molarMassGPerMol.known
            && recA->properties.molarMassGPerMol.source == SacePropertySource::IdentityDerived
            && recA->properties.molarMassGPerMol.confidence == SaceConfidence::High
            && recA->properties.criticalTemperatureK.known
            && recA->properties.criticalTemperatureK.source == SacePropertySource::StructuralEstimate
            && recA->properties.criticalTemperatureK.confidence == SaceConfidence::Low
            && recA->properties.criticalPressurePa.known
            && recA->properties.criticalPressurePa.source == SacePropertySource::StructuralEstimate
            && recA->properties.criticalMolarVolumeM3PerMol.known
            && recA->properties.criticalMolarVolumeM3PerMol.source == SacePropertySource::StructuralEstimate
            && recA->properties.acentricFactor.known
            && recA->properties.acentricFactor.source == SacePropertySource::StructuralEstimate
            && recA->properties.acentricFactor.confidence == SaceConfidence::Low
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.source
                == SacePropertySource::StructuralEstimate
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.confidence
                == SaceConfidence::Low
            && nearK(recA->properties.criticalTemperatureK.value, 499.407, 0.05), "");

    SaceGeneratedRecord *mutA = cat.recordMutable(refA.generatedId);
    if (mutA) {
        mutA->properties.criticalTemperatureK = saceUnknownScalarProperty();
        mutA->properties.criticalPressurePa = saceUnknownScalarProperty();
        mutA->properties.criticalMolarVolumeM3PerMol = saceUnknownScalarProperty();
        mutA->properties.acentricFactor = saceUnknownScalarProperty();
        mutA->properties.enthalpyVaporizationAtNormalBoilingJPerMol = saceUnknownScalarProperty();
    }
    bool reattachFill = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_critical",
        reattachFill && recA && recA->properties.criticalTemperatureK.known
            && recA->properties.criticalPressurePa.known
            && recA->properties.criticalMolarVolumeM3PerMol.known
            && recA->properties.criticalTemperatureK.source == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.criticalTemperatureK.value, 499.407, 0.05), "");

    SaceGeneratedRecord *mutA2 = cat.recordMutable(refA.generatedId);
    if (mutA2) {
        mutA2->properties.criticalTemperatureK.known = true;
        mutA2->properties.criticalTemperatureK.value = 513.9f;
        mutA2->properties.criticalTemperatureK.source = SacePropertySource::Reference;
        mutA2->properties.criticalTemperatureK.confidence = SaceConfidence::High;
    }
    bool reattachKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high",
        reattachKeep && recA
            && recA->properties.criticalTemperatureK.source == SacePropertySource::Reference
            && recA->properties.criticalTemperatureK.confidence == SaceConfidence::High
            && nearK(recA->properties.criticalTemperatureK.value, 513.9, 0.05), "");

    auto omegaOk = [](SaceScalarProperty const &p) {
        return p.known
            && p.source == SacePropertySource::StructuralEstimate
            && p.confidence == SaceConfidence::Low;
    };
    emit("a_omega_0_556081", bunA.acentricFactor.known
        && nearK(bunA.acentricFactor.value, 0.556081, 2e-5),
        std::to_string(bunA.acentricFactor.value));
    emit("b_omega_0_193764", bunB.acentricFactor.known
        && nearK(bunB.acentricFactor.value, 0.193764, 2e-5),
        std::to_string(bunB.acentricFactor.value));
    emit("a_b_omega_structural_estimate_low",
        omegaOk(bunA.acentricFactor) && omegaOk(bunB.acentricFactor), "");

    SaceLeeKeslerVaporModel modelA{}, modelB{};
    emit("a_vapor_model_from_joback_bundle",
        saceBuildLeeKeslerVaporModelFromJoback(bunA, modelA) && modelA.valid, "");
    emit("b_vapor_model_from_joback_bundle",
        saceBuildLeeKeslerVaporModelFromJoback(bunB, modelB) && modelB.valid, "");

    auto psatAt = [](SaceLeeKeslerVaporModel const &m, double t, double &p) {
        return m.valid && saceLeeKeslerSaturationPressurePa(
            t, m.criticalTemperatureK, m.criticalPressurePa, m.acentricFactor, p);
    };
    double aTbPsat = 0, bTbPsat = 0, aTcPsat = 0, bTcPsat = 0;
    bool aTbOk = psatAt(modelA, modelA.normalBoilingPointK, aTbPsat)
        && nearK(aTbPsat, kLeeKeslerAtmPa, 50.0);
    emit("a_psat_at_tb_approx_1_atm", aTbOk, std::to_string(aTbPsat));
    bool bTbOk = psatAt(modelB, modelB.normalBoilingPointK, bTbPsat)
        && nearK(bTbPsat, kLeeKeslerAtmPa, 50.0);
    emit("b_psat_at_tb_approx_1_atm", bTbOk, std::to_string(bTbPsat));
    bool aTcOk = psatAt(modelA, modelA.criticalTemperatureK, aTcPsat)
        && aTcPsat == modelA.criticalPressurePa;
    emit("a_psat_at_tc_equals_pc", aTcOk, std::to_string(aTcPsat));
    bool bTcOk = psatAt(modelB, modelB.criticalTemperatureK, bTcPsat)
        && bTcPsat == modelB.criticalPressurePa;
    emit("b_psat_at_tc_equals_pc", bTcOk, std::to_string(bTcPsat));

    auto omegaDefOk = [&](SaceLeeKeslerVaporModel const &m, double expect) {
        double t = 0.7 * m.criticalTemperatureK;
        double p = 0.0;
        if (!psatAt(m, t, p) || !(p > 0.0) || !(m.criticalPressurePa > 0.0))
            return false;
        double pr = p / m.criticalPressurePa;
        if (!(pr > 0.0) || !std::isfinite(pr))
            return false;
        double omegaDef = -std::log10(pr) - 1.0;
        return nearK(omegaDef, expect, 2e-4) && nearK(omegaDef, m.acentricFactor, 2e-4);
    };
    emit("a_omega_definition_at_tr_0_7", omegaDefOk(modelA, bunA.acentricFactor.value), "");
    emit("b_omega_definition_at_tr_0_7", omegaDefOk(modelB, bunB.acentricFactor.value), "");

    auto checkTr = [&](char const *name, SaceLeeKeslerVaporModel const &m, double tr, double expect, double tol) {
        double p = 0.0;
        bool ok = psatAt(m, tr * m.criticalTemperatureK, p) && nearK(p, expect, tol);
        emit(name, ok, std::to_string(p));
        return p;
    };
    double a06 = checkTr("a_psat_tr_0_6", modelA, 0.6, 17651.9, 20.0);
    double a07 = checkTr("a_psat_tr_0_7", modelA, 0.7, 159999.0, 80.0);
    double a08 = checkTr("a_psat_tr_0_8", modelA, 0.8, 741479.0, 200.0);
    double a09 = checkTr("a_psat_tr_0_9", modelA, 0.9, 2.30695e6, 400.0);
    double b06 = checkTr("b_psat_tr_0_6", modelB, 0.6, 63957.8, 40.0);
    double b07 = checkTr("b_psat_tr_0_7", modelB, 0.7, 314341.0, 150.0);
    double b08 = checkTr("b_psat_tr_0_8", modelB, 0.8, 990254.0, 300.0);
    double b09 = checkTr("b_psat_tr_0_9", modelB, 0.9, 2.38267e6, 400.0);
    emit("a_psat_monotonic_sampled_tr",
        a06 < a07 && a07 < a08 && a08 < a09 && a09 < modelA.criticalPressurePa, "");
    emit("b_psat_monotonic_sampled_tr",
        b06 < b07 && b07 < b08 && b08 < b09 && b09 < modelB.criticalPressurePa, "");

    double rejected = 0.0;
    emit("psat_rejects_t_above_tc",
        !saceLeeKeslerSaturationPressurePa(modelA.criticalTemperatureK + 1.0,
            modelA.criticalTemperatureK, modelA.criticalPressurePa, modelA.acentricFactor, rejected), "");
    emit("psat_rejects_nonpositive_t",
        !saceLeeKeslerSaturationPressurePa(0.0, modelA.criticalTemperatureK,
            modelA.criticalPressurePa, modelA.acentricFactor, rejected)
            && !saceLeeKeslerSaturationPressurePa(-1.0, modelA.criticalTemperatureK,
                modelA.criticalPressurePa, modelA.acentricFactor, rejected), "");
    emit("psat_rejects_invalid_tc",
        !saceLeeKeslerSaturationPressurePa(300.0, 0.0, modelA.criticalPressurePa,
            modelA.acentricFactor, rejected), "");
    emit("psat_rejects_invalid_pc",
        !saceLeeKeslerSaturationPressurePa(300.0, modelA.criticalTemperatureK, 0.0,
            modelA.acentricFactor, rejected), "");
    double nan = std::numeric_limits<double>::quiet_NaN();
    double inf = std::numeric_limits<double>::infinity();
    emit("psat_rejects_nan_inf_inputs",
        !saceLeeKeslerSaturationPressurePa(nan, modelA.criticalTemperatureK,
            modelA.criticalPressurePa, modelA.acentricFactor, rejected)
            && !saceLeeKeslerSaturationPressurePa(300.0, inf, modelA.criticalPressurePa,
                modelA.acentricFactor, rejected)
            && !saceLeeKeslerSaturationPressurePa(300.0, modelA.criticalTemperatureK, nan,
                modelA.acentricFactor, rejected)
            && !saceLeeKeslerSaturationPressurePa(300.0, modelA.criticalTemperatureK,
                modelA.criticalPressurePa, inf, rejected), "");

    emit("water_no_joback_omega", !waterB.acentricFactor.known, "");
    emit("methane_no_joback_omega", !methaneB.acentricFactor.known, "");
    emit("omega_not_in_canonical_signature",
        std::strstr(sig, "omega") == nullptr && std::strstr(sig, "0.556") == nullptr
            && std::strstr(sig, "Lee") == nullptr, sig);
    emit("vapor_model_not_in_canonical_identity",
        std::strstr(sig, "vapor") == nullptr && std::strstr(sig, "Psat") == nullptr
            && std::strstr(sig, "acentric") == nullptr, sig);

    SaceGeneratedRecord *mutOmega = cat.recordMutable(refA.generatedId);
    if (mutOmega)
        mutOmega->properties.acentricFactor = saceUnknownScalarProperty();
    bool reattachOmega = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_omega",
        reattachOmega && recA && recA->properties.acentricFactor.known
            && recA->properties.acentricFactor.source == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.acentricFactor.value, 0.556081, 2e-5), "");

    SaceGeneratedRecord *mutOmegaRef = cat.recordMutable(refA.generatedId);
    if (mutOmegaRef) {
        mutOmegaRef->properties.acentricFactor.known = true;
        mutOmegaRef->properties.acentricFactor.value = 0.644f;
        mutOmegaRef->properties.acentricFactor.source = SacePropertySource::Reference;
        mutOmegaRef->properties.acentricFactor.confidence = SaceConfidence::High;
    }
    bool reattachOmegaKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_omega",
        reattachOmegaKeep && recA
            && recA->properties.acentricFactor.source == SacePropertySource::Reference
            && recA->properties.acentricFactor.confidence == SaceConfidence::High
            && nearK(recA->properties.acentricFactor.value, 0.644, 1e-4), "");
    emit("generated_record_still_unspawnable", recA && !recA->spawnable, "");

    double aHvapSum = saceJobackHvapContributionSumKJPerMol(bunA.groups);
    double bHvapSum = saceJobackHvapContributionSumKJPerMol(bunB.groups);
    emit("a_hvap_contribution_sum_21_425", nearK(aHvapSum, 21.425, 1e-6), std::to_string(aHvapSum));
    emit("a_hvap_tb_36725", bunA.enthalpyVaporizationAtNormalBoilingJPerMol.known
        && nearK(bunA.enthalpyVaporizationAtNormalBoilingJPerMol.value, 36725.0, 1.0),
        std::to_string(bunA.enthalpyVaporizationAtNormalBoilingJPerMol.value));
    emit("b_hvap_contribution_sum_7_156", nearK(bHvapSum, 7.156, 1e-6), std::to_string(bHvapSum));
    emit("b_hvap_tb_22456", bunB.enthalpyVaporizationAtNormalBoilingJPerMol.known
        && nearK(bunB.enthalpyVaporizationAtNormalBoilingJPerMol.value, 22456.0, 1.0),
        std::to_string(bunB.enthalpyVaporizationAtNormalBoilingJPerMol.value));
    emit("a_b_hvap_differ", bunA.enthalpyVaporizationAtNormalBoilingJPerMol.known
        && bunB.enthalpyVaporizationAtNormalBoilingJPerMol.known
        && bunA.enthalpyVaporizationAtNormalBoilingJPerMol.value
            > bunB.enthalpyVaporizationAtNormalBoilingJPerMol.value, "");
    emit("a_b_hvap_structural_estimate_low",
        bunA.enthalpyVaporizationAtNormalBoilingJPerMol.source == SacePropertySource::StructuralEstimate
            && bunA.enthalpyVaporizationAtNormalBoilingJPerMol.confidence == SaceConfidence::Low
            && bunB.enthalpyVaporizationAtNormalBoilingJPerMol.source == SacePropertySource::StructuralEstimate
            && bunB.enthalpyVaporizationAtNormalBoilingJPerMol.confidence == SaceConfidence::Low, "");

    SaceWatsonVaporizationModel watA{}, watB{};
    emit("a_watson_model_builds",
        saceBuildWatsonVaporizationModelFromJoback(gA, bunA, watA) && watA.valid, "");
    emit("b_watson_model_builds",
        saceBuildWatsonVaporizationModelFromJoback(gB, bunB, watB) && watB.valid, "");
    emit("a_watson_molar_mass_0_046069",
        watA.valid && nearK(watA.molarMassKgPerMol, 0.046069, 1e-6),
        std::to_string(watA.molarMassKgPerMol));
    emit("b_watson_molar_mass_0_046069",
        watB.valid && nearK(watB.molarMassKgPerMol, 0.046069, 1e-6),
        std::to_string(watB.molarMassKgPerMol));

    auto watsonAt = [](SaceWatsonVaporizationModel const &m, double t, double &h) {
        return saceWatsonEnthalpyVaporizationJPerMol(m, t, h);
    };
    double aHref = 0, bHref = 0, aHtc = 0, bHtc = 0;
    bool aHrefOk = watsonAt(watA, watA.referenceTemperatureK, aHref) && nearK(aHref, 36725.0, 1.0);
    emit("a_hvap_at_tref_36725", aHrefOk, std::to_string(aHref));
    bool bHrefOk = watsonAt(watB, watB.referenceTemperatureK, bHref) && nearK(bHref, 22456.0, 1.0);
    emit("b_hvap_at_tref_22456", bHrefOk, std::to_string(bHref));
    bool aHtcOk = watsonAt(watA, watA.criticalTemperatureK, aHtc) && aHtc == 0.0;
    emit("a_hvap_at_tc_zero", aHtcOk, std::to_string(aHtc));
    bool bHtcOk = watsonAt(watB, watB.criticalTemperatureK, bHtc) && bHtc == 0.0;
    emit("b_hvap_at_tc_zero", bHtcOk, std::to_string(bHtc));

    auto checkWatsonTr = [&](char const *name, SaceWatsonVaporizationModel const &m,
        double tr, double expect, double tol)
    {
        double h = 0.0;
        bool ok = watsonAt(m, tr * m.criticalTemperatureK, h) && nearK(h, expect, tol);
        emit(name, ok, std::to_string(h));
        return h;
    };
    double aW06 = checkWatsonTr("a_watson_tr_0_6", watA, 0.6, 39781.1, 2.0);
    double aW07 = checkWatsonTr("a_watson_tr_0_7", watA, 0.7, 35661.6, 2.0);
    double aW08 = checkWatsonTr("a_watson_tr_0_8", watA, 0.8, 30569.3, 2.0);
    double aW09 = checkWatsonTr("a_watson_tr_0_9", watA, 0.9, 23490.6, 2.0);
    double bW06 = checkWatsonTr("b_watson_tr_0_6", watB, 0.6, 23022.9, 2.0);
    double bW07 = checkWatsonTr("b_watson_tr_0_7", watB, 0.7, 20638.7, 2.0);
    double bW08 = checkWatsonTr("b_watson_tr_0_8", watB, 0.8, 17691.6, 2.0);
    double bW09 = checkWatsonTr("b_watson_tr_0_9", watB, 0.9, 13594.9, 2.0);
    (void)aW06;
    (void)bW06;

    double aKg = 0, bKg = 0;
    bool aKgOk = saceWatsonLatentHeatVaporizationJPerKg(watA, watA.referenceTemperatureK, aKg)
        && nearK(aKg, 797174.0, 50.0);
    emit("a_specific_latent_tb_797174", aKgOk, std::to_string(aKg));
    bool bKgOk = saceWatsonLatentHeatVaporizationJPerKg(watB, watB.referenceTemperatureK, bKg)
        && nearK(bKg, 487443.0, 50.0);
    emit("b_specific_latent_tb_487443", bKgOk, std::to_string(bKg));

    double rejectedH = 0.0;
    emit("watson_rejects_t_above_tc",
        !saceWatsonEnthalpyVaporizationJPerMol(watA, watA.criticalTemperatureK + 1.0, rejectedH), "");
    emit("watson_rejects_nonpositive_t",
        !saceWatsonEnthalpyVaporizationJPerMol(watA, 0.0, rejectedH)
            && !saceWatsonEnthalpyVaporizationJPerMol(watA, -1.0, rejectedH), "");
    double nanH = std::numeric_limits<double>::quiet_NaN();
    double infH = std::numeric_limits<double>::infinity();
    emit("watson_rejects_nan_inf",
        !saceWatsonEnthalpyVaporizationJPerMol(watA, nanH, rejectedH)
            && !saceWatsonEnthalpyVaporizationJPerMol(watA, infH, rejectedH), "");
    SaceWatsonVaporizationModel invalidW{};
    emit("watson_rejects_invalid_model",
        !saceWatsonEnthalpyVaporizationJPerMol(invalidW, 300.0, rejectedH), "");
    emit("sampled_hvap_decreases_toward_tc",
        aW07 > aW08 && aW08 > aW09 && aW09 > 0.0
            && bW07 > bW08 && bW08 > bW09 && bW09 > 0.0, "");
    emit("water_no_joback_hvap", !waterB.enthalpyVaporizationAtNormalBoilingJPerMol.known, "");
    emit("methane_no_joback_hvap", !methaneB.enthalpyVaporizationAtNormalBoilingJPerMol.known, "");
    emit("hvap_not_in_canonical_signature",
        std::strstr(sig, "36725") == nullptr && std::strstr(sig, "Hvap") == nullptr
            && std::strstr(sig, "Watson") == nullptr, sig);
    emit("watson_model_not_in_canonical_identity",
        std::strstr(sig, "watson") == nullptr && std::strstr(sig, "latent") == nullptr, sig);

    SaceGeneratedRecord *mutHvap = cat.recordMutable(refA.generatedId);
    if (mutHvap)
        mutHvap->properties.enthalpyVaporizationAtNormalBoilingJPerMol = saceUnknownScalarProperty();
    bool reattachHvap = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_hvap",
        reattachHvap && recA && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.known
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.source
                == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.value, 36725.0, 1.0), "");

    SaceGeneratedRecord *mutHvapRef = cat.recordMutable(refA.generatedId);
    if (mutHvapRef) {
        mutHvapRef->properties.enthalpyVaporizationAtNormalBoilingJPerMol.known = true;
        mutHvapRef->properties.enthalpyVaporizationAtNormalBoilingJPerMol.value = 38600.0f;
        mutHvapRef->properties.enthalpyVaporizationAtNormalBoilingJPerMol.source
            = SacePropertySource::Reference;
        mutHvapRef->properties.enthalpyVaporizationAtNormalBoilingJPerMol.confidence
            = SaceConfidence::High;
    }
    bool reattachHvapKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_hvap",
        reattachHvapKeep && recA
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.source
                == SacePropertySource::Reference
            && recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.confidence
                == SaceConfidence::High
            && nearK(recA->properties.enthalpyVaporizationAtNormalBoilingJPerMol.value, 38600.0, 1.0), "");
    emit("generated_record_still_unspawnable_after_hvap", recA && !recA->spawnable, "");

    SaceJobackIdealGasCpModel cpA{}, cpB{};
    emit("a_cp_polynomial_builds", saceBuildJobackIdealGasCpModel(bunA.groups, cpA) && cpA.valid, "");
    emit("a_cp_coefficients",
        cpA.valid
            && nearK(cpA.A, 6.361, 1e-6)
            && nearK(cpA.B, 0.22782, 1e-7)
            && nearK(cpA.C, -1.154e-4, 1e-9)
            && nearK(cpA.D, 2.24e-8, 1e-11),
        std::to_string(cpA.A) + "," + std::to_string(cpA.B) + ","
            + std::to_string(cpA.C) + "," + std::to_string(cpA.D));
    emit("b_cp_polynomial_builds", saceBuildJobackIdealGasCpModel(bunB.groups, cpB) && cpB.valid, "");
    emit("b_cp_coefficients",
        cpB.valid
            && nearK(cpB.A, 26.57, 1e-6)
            && nearK(cpB.B, 0.13064, 1e-7)
            && nearK(cpB.C, 2.60e-5, 1e-9)
            && nearK(cpB.D, -4.22e-8, 1e-11),
        std::to_string(cpB.A) + "," + std::to_string(cpB.B) + ","
            + std::to_string(cpB.C) + "," + std::to_string(cpB.D));

    double aIg298 = 0, bIg298 = 0, aIg300 = 0, bIg300 = 0, aIg1000 = 0, bIg1000 = 0;
    bool aIg298Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpA, 298.15, aIg298);
    bool bIg298Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpB, 298.15, bIg298);
    bool aIg300Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpA, 300.0, aIg300);
    bool bIg300Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpB, 300.0, bIg300);
    bool aIg1000Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpA, 1000.0, aIg1000);
    bool bIg1000Ok = saceJobackIdealGasHeatCapacityJPerMolK(cpB, 1000.0, bIg1000);
    emit("a_cpig_298_15", aIg298Ok && nearK(aIg298, 64.6209, 0.005), std::to_string(aIg298));
    emit("b_cpig_298_15", bIg298Ok && nearK(bIg298, 66.7131, 0.005), std::to_string(bIg298));
    emit("a_cpig_300", aIg300Ok && nearK(aIg300, 64.9258, 0.005), std::to_string(aIg300));
    emit("b_cpig_300", bIg300Ok && nearK(bIg300, 66.9626, 0.005), std::to_string(bIg300));
    emit("a_cpig_1000", aIg1000Ok && nearK(aIg1000, 141.181, 0.02), std::to_string(aIg1000));
    emit("b_cpig_1000", bIg1000Ok && nearK(bIg1000, 141.010, 0.02), std::to_string(bIg1000));

    double rejectedCp = 0.0;
    emit("cpig_rejects_below_298",
        !saceJobackIdealGasHeatCapacityJPerMolK(cpA, 297.0, rejectedCp) && rejectedCp == 0.0, "");
    emit("cpig_rejects_above_1000",
        !saceJobackIdealGasHeatCapacityJPerMolK(cpA, 1000.01, rejectedCp) && rejectedCp == 0.0, "");
    emit("cpig_rejects_nan_inf",
        !saceJobackIdealGasHeatCapacityJPerMolK(cpA, nanH, rejectedCp)
            && !saceJobackIdealGasHeatCapacityJPerMolK(cpA, infH, rejectedCp), "");

    SaceRowlinsonPolingLiquidCpModel liqA{}, liqB{};
    emit("a_rowlinson_poling_model_builds",
        saceBuildRowlinsonPolingLiquidCpModelFromJoback(gA, bunA, liqA) && liqA.valid, "");
    emit("b_rowlinson_poling_model_builds",
        saceBuildRowlinsonPolingLiquidCpModelFromJoback(gB, bunB, liqB) && liqB.valid, "");

    double aL298 = 0, bL298 = 0;
    bool aL298Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA, 298.15, aL298);
    bool bL298Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqB, 298.15, bL298);
    emit("a_cpl_298_15", aL298Ok && nearK(aL298, 148.729, 0.02), std::to_string(aL298));
    emit("b_cpl_298_15", bL298Ok && nearK(bL298, 112.315, 0.02), std::to_string(bL298));

    double aT08 = 0.8 * liqA.criticalTemperatureK;
    double bT08 = 0.8 * liqB.criticalTemperatureK;
    double aL08 = 0, bL08 = 0, aL09 = 0, bL09 = 0;
    bool aL08Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA, aT08, aL08);
    bool bL08Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqB, bT08, bL08);
    bool aL09Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA, 0.9 * liqA.criticalTemperatureK, aL09);
    bool bL09Ok = saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqB, 0.9 * liqB.criticalTemperatureK, bL09);
    emit("a_cpl_tr_0_8",
        nearK(aT08, 399.526, 0.05) && aL08Ok && nearK(aL08, 165.084, 0.05),
        std::to_string(aL08));
    emit("b_cpl_tr_0_8",
        nearK(bT08, 342.539, 0.05) && bL08Ok && nearK(bL08, 124.050, 0.05),
        std::to_string(bL08));
    emit("a_cpl_tr_0_9", aL09Ok && nearK(aL09, 196.342, 0.05), std::to_string(aL09));
    emit("b_cpl_tr_0_9", bL09Ok && nearK(bL09, 151.428, 0.05), std::to_string(bL09));

    emit("liquid_cp_rejects_tr_ge_0_98",
        !saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA,
            kSaceLiquidCpTrReject * liqA.criticalTemperatureK, rejectedCp), "");
    emit("liquid_cp_rejects_t_ge_tc",
        !saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA, liqA.criticalTemperatureK, rejectedCp)
            && !saceRowlinsonPolingLiquidHeatCapacityJPerMolK(liqA,
                liqA.criticalTemperatureK + 1.0, rejectedCp), "");
    SaceRowlinsonPolingLiquidCpModel invalidLiq{};
    emit("liquid_cp_rejects_invalid_model",
        !saceRowlinsonPolingLiquidHeatCapacityJPerMolK(invalidLiq, 298.15, rejectedCp), "");

    double aGasSpec = 0, bGasSpec = 0, aLiqSpec = 0, bLiqSpec = 0;
    bool aGasSpecOk = saceMolarHeatCapacityToSpecificJPerKgK(aIg298, 0.046069, aGasSpec);
    bool bGasSpecOk = saceMolarHeatCapacityToSpecificJPerKgK(bIg298, 0.046069, bGasSpec);
    bool aLiqSpecOk = saceMolarHeatCapacityToSpecificJPerKgK(aL298, 0.046069, aLiqSpec);
    bool bLiqSpecOk = saceMolarHeatCapacityToSpecificJPerKgK(bL298, 0.046069, bLiqSpec);
    emit("a_gas_specific_cp_298_15",
        aGasSpecOk && nearK(aGasSpec, 1402.70, 0.5), std::to_string(aGasSpec));
    emit("b_gas_specific_cp_298_15",
        bGasSpecOk && nearK(bGasSpec, 1448.11, 0.5), std::to_string(bGasSpec));
    emit("a_liquid_specific_cp_298_15",
        aLiqSpecOk && nearK(aLiqSpec, 3228.40, 1.0), std::to_string(aLiqSpec));
    emit("b_liquid_specific_cp_298_15",
        bLiqSpecOk && nearK(bLiqSpec, 2437.97, 1.0), std::to_string(bLiqSpec));

    emit("cp_reference_structural_estimate_low",
        bunA.idealGasHeatCapacityAt298KJPerMolK.source == SacePropertySource::StructuralEstimate
            && bunA.idealGasHeatCapacityAt298KJPerMolK.confidence == SaceConfidence::Low
            && bunA.saturatedLiquidHeatCapacityAt298KJPerMolK.source == SacePropertySource::StructuralEstimate
            && bunA.saturatedLiquidHeatCapacityAt298KJPerMolK.confidence == SaceConfidence::Low
            && bunB.idealGasHeatCapacityAt298KJPerMolK.source == SacePropertySource::StructuralEstimate
            && bunB.saturatedLiquidHeatCapacityAt298KJPerMolK.source == SacePropertySource::StructuralEstimate
            && recA && recA->properties.idealGasHeatCapacityAt298KJPerMolK.source
                == SacePropertySource::StructuralEstimate
            && recA->properties.idealGasHeatCapacityAt298KJPerMolK.confidence == SaceConfidence::Low
            && recB && recB->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.source
                == SacePropertySource::StructuralEstimate, "");

    auto cpModelMatches = [&](SaceJobackIdealGasCpModel const &x, SaceJobackIdealGasCpModel const &y) {
        return x.valid && y.valid
            && nearK(x.A, y.A, 1e-9) && nearK(x.B, y.B, 1e-9)
            && nearK(x.C, y.C, 1e-12) && nearK(x.D, y.D, 1e-14);
    };
    SaceJobackIdealGasCpModel permCp{}, bondCp{}, endCp{};
    SaceRowlinsonPolingLiquidCpModel permLiq{}, bondLiq{}, endLiq{};
    double permIg = 0, permL = 0, bondIg = 0, bondL = 0, endIg = 0, endL = 0;
    emit("atom_permutation_preserves_cp",
        saceBuildJobackIdealGasCpModel(permB.groups, permCp)
            && cpModelMatches(cpA, permCp)
            && saceBuildRowlinsonPolingLiquidCpModelFromJoback(permuteGraph(gA), permB, permLiq)
            && saceJobackIdealGasHeatCapacityJPerMolK(permCp, 298.15, permIg)
            && saceRowlinsonPolingLiquidHeatCapacityJPerMolK(permLiq, 298.15, permL)
            && nearK(permIg, aIg298, 1e-9) && nearK(permL, aL298, 1e-6), "");
    emit("bond_reorder_preserves_cp",
        saceBuildJobackIdealGasCpModel(bondB.groups, bondCp)
            && cpModelMatches(cpA, bondCp)
            && saceBuildRowlinsonPolingLiquidCpModelFromJoback(reorderBonds(gA), bondB, bondLiq)
            && saceJobackIdealGasHeatCapacityJPerMolK(bondCp, 298.15, bondIg)
            && saceRowlinsonPolingLiquidHeatCapacityJPerMolK(bondLiq, 298.15, bondL)
            && nearK(bondIg, aIg298, 1e-9) && nearK(bondL, aL298, 1e-6), "");
    emit("reversed_endpoints_preserve_cp",
        saceBuildJobackIdealGasCpModel(endB.groups, endCp)
            && cpModelMatches(cpA, endCp)
            && saceBuildRowlinsonPolingLiquidCpModelFromJoback(reverseEnds(gA), endB, endLiq)
            && saceJobackIdealGasHeatCapacityJPerMolK(endCp, 298.15, endIg)
            && saceRowlinsonPolingLiquidHeatCapacityJPerMolK(endLiq, 298.15, endL)
            && nearK(endIg, aIg298, 1e-9) && nearK(endL, aL298, 1e-6), "");

    emit("water_no_joback_cp",
        !waterB.idealGasHeatCapacityAt298KJPerMolK.known
            && !waterB.saturatedLiquidHeatCapacityAt298KJPerMolK.known, "");
    emit("methane_no_joback_cp",
        !methaneB.idealGasHeatCapacityAt298KJPerMolK.known
            && !methaneB.saturatedLiquidHeatCapacityAt298KJPerMolK.known, "");
    emit("cp_not_in_canonical_identity",
        std::strstr(sig, "64.62") == nullptr && std::strstr(sig, "Cp") == nullptr
            && std::strstr(sig, "Rowlinson") == nullptr && std::strstr(sig, "heat") == nullptr, sig);

    emit("catalog_cp_matches_direct_evaluation",
        recA && recB
            && nearK(recA->properties.idealGasHeatCapacityAt298KJPerMolK.value, aIg298, 0.05)
            && nearK(recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.value, aL298, 0.05)
            && nearK(recB->properties.idealGasHeatCapacityAt298KJPerMolK.value, bIg298, 0.05)
            && nearK(recB->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.value, bL298, 0.05)
            && nearK(bunA.idealGasHeatCapacityAt298KJPerMolK.value, aIg298, 0.05)
            && nearK(bunA.saturatedLiquidHeatCapacityAt298KJPerMolK.value, aL298, 0.05), "");

    SaceGeneratedRecord *mutCp = cat.recordMutable(refA.generatedId);
    if (mutCp) {
        mutCp->properties.idealGasHeatCapacityAt298KJPerMolK = saceUnknownScalarProperty();
        mutCp->properties.saturatedLiquidHeatCapacityAt298KJPerMolK = saceUnknownScalarProperty();
    }
    bool reattachCp = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_cp",
        reattachCp && recA
            && recA->properties.idealGasHeatCapacityAt298KJPerMolK.known
            && recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.known
            && recA->properties.idealGasHeatCapacityAt298KJPerMolK.source
                == SacePropertySource::StructuralEstimate
            && recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.source
                == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.idealGasHeatCapacityAt298KJPerMolK.value, aIg298, 0.05)
            && nearK(recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.value, aL298, 0.05), "");

    SaceGeneratedRecord *mutCpRef = cat.recordMutable(refA.generatedId);
    if (mutCpRef) {
        mutCpRef->properties.idealGasHeatCapacityAt298KJPerMolK.known = true;
        mutCpRef->properties.idealGasHeatCapacityAt298KJPerMolK.value = 65.0f;
        mutCpRef->properties.idealGasHeatCapacityAt298KJPerMolK.source = SacePropertySource::Reference;
        mutCpRef->properties.idealGasHeatCapacityAt298KJPerMolK.confidence = SaceConfidence::High;
        mutCpRef->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.known = true;
        mutCpRef->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.value = 112.0f;
        mutCpRef->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.source = SacePropertySource::Reference;
        mutCpRef->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.confidence = SaceConfidence::High;
    }
    bool reattachCpKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_cp",
        reattachCpKeep && recA
            && recA->properties.idealGasHeatCapacityAt298KJPerMolK.source == SacePropertySource::Reference
            && recA->properties.idealGasHeatCapacityAt298KJPerMolK.confidence == SaceConfidence::High
            && recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.source
                == SacePropertySource::Reference
            && recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.confidence
                == SaceConfidence::High
            && nearK(recA->properties.idealGasHeatCapacityAt298KJPerMolK.value, 65.0, 0.02)
            && nearK(recA->properties.saturatedLiquidHeatCapacityAt298KJPerMolK.value, 112.0, 0.02), "");
    emit("generated_record_still_unspawnable_after_cp", recA && !recA->spawnable, "");

    SaceCostaldLiquidDensityModel densA{}, densB{};
    emit("a_costald_model_builds",
        saceBuildCostaldLiquidDensityModelFromJoback(gA, bunA, densA) && densA.valid, "");
    emit("b_costald_model_builds",
        saceBuildCostaldLiquidDensityModelFromJoback(gB, bunB, densB) && densB.valid, "");
    emit("a_costald_tc", densA.valid && nearK(densA.criticalTemperatureK, 499.407379, 0.05),
        std::to_string(densA.criticalTemperatureK));
    emit("a_costald_vstar", densA.valid && nearK(densA.characteristicVolumeM3PerMol, 1.665e-4, 1e-7),
        std::to_string(densA.characteristicVolumeM3PerMol));
    emit("a_costald_omega", densA.valid && nearK(densA.acentricFactor, 0.556081, 2e-5),
        std::to_string(densA.acentricFactor));
    emit("a_costald_molar_mass", densA.valid && nearK(densA.molarMassKgPerMol, 0.046069, 1e-6),
        std::to_string(densA.molarMassKgPerMol));
    emit("b_costald_inputs",
        densB.valid
            && nearK(densB.criticalTemperatureK, 428.173981, 0.05)
            && nearK(densB.characteristicVolumeM3PerMol, 1.655e-4, 1e-7)
            && nearK(densB.acentricFactor, 0.193764, 2e-5)
            && nearK(densB.molarMassKgPerMol, 0.046069, 1e-6),
        std::to_string(densB.criticalTemperatureK));

    double aTr298 = 298.15 / densA.criticalTemperatureK;
    double bTr298 = 298.15 / densB.criticalTemperatureK;
    emit("a_costald_tr_298_15", nearK(aTr298, 0.597008, 2e-6), std::to_string(aTr298));
    emit("b_costald_tr_298_15", nearK(bTr298, 0.696329, 2e-6), std::to_string(bTr298));

    double aVs298 = 0, aRho298 = 0, bVs298 = 0, bRho298 = 0;
    bool aVs298Ok = saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA, 298.15, aVs298);
    bool aRho298Ok = saceCostaldSaturatedLiquidDensityKgPerM3(densA, 298.15, aRho298);
    bool bVs298Ok = saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densB, 298.15, bVs298);
    bool bRho298Ok = saceCostaldSaturatedLiquidDensityKgPerM3(densB, 298.15, bRho298);
    emit("a_vs_298_15", aVs298Ok && nearK(aVs298, 5.60681e-5, 2e-9), std::to_string(aVs298));
    emit("a_rho_298_15", aRho298Ok && nearK(aRho298, 821.662, 0.05), std::to_string(aRho298));
    emit("b_vs_298_15", bVs298Ok && nearK(bVs298, 6.58211e-5, 2e-9), std::to_string(bVs298));
    emit("b_rho_298_15", bRho298Ok && nearK(bRho298, 699.912, 0.05), std::to_string(bRho298));
    emit("a_density_gt_b_at_298_15", aRho298Ok && bRho298Ok && aRho298 > bRho298, "");

    auto costaldAtTr = [&](SaceCostaldLiquidDensityModel const &m, double tr,
        double &vs, double &rho)
    {
        double t = tr * m.criticalTemperatureK;
        return saceCostaldSaturatedLiquidMolarVolumeM3PerMol(m, t, vs)
            && saceCostaldSaturatedLiquidDensityKgPerM3(m, t, rho);
    };
    double aVs06 = 0, aRho06 = 0, aVs07 = 0, aRho07 = 0, aVs08 = 0, aRho08 = 0, aVs09 = 0, aRho09 = 0;
    double bVs06 = 0, bRho06 = 0, bVs07 = 0, bRho07 = 0, bVs08 = 0, bRho08 = 0, bVs09 = 0, bRho09 = 0;
    bool a06Ok = costaldAtTr(densA, 0.6, aVs06, aRho06);
    bool a07Ok = costaldAtTr(densA, 0.7, aVs07, aRho07);
    bool a08Ok = costaldAtTr(densA, 0.8, aVs08, aRho08);
    bool a09Ok = costaldAtTr(densA, 0.9, aVs09, aRho09);
    bool b06Ok = costaldAtTr(densB, 0.6, bVs06, bRho06);
    bool b07Ok = costaldAtTr(densB, 0.7, bVs07, bRho07);
    bool b08Ok = costaldAtTr(densB, 0.8, bVs08, bRho08);
    bool b09Ok = costaldAtTr(densB, 0.9, bVs09, bRho09);
    emit("a_costald_tr_0_6",
        a06Ok && nearK(aVs06, 5.61984e-5, 2e-9) && nearK(aRho06, 819.756, 0.05),
        std::to_string(aRho06));
    emit("a_costald_tr_0_7",
        a07Ok && nearK(aVs07, 6.11866e-5, 2e-9) && nearK(aRho07, 752.926, 0.05),
        std::to_string(aRho07));
    emit("a_costald_tr_0_8",
        a08Ok && nearK(aVs08, 6.80648e-5, 2e-9) && nearK(aRho08, 676.840, 0.05),
        std::to_string(aRho08));
    emit("a_costald_tr_0_9",
        a09Ok && nearK(aVs09, 7.93532e-5, 2e-9) && nearK(aRho09, 580.557, 0.05),
        std::to_string(aRho09));
    emit("b_costald_tr_0_6",
        b06Ok && nearK(bVs06, 6.10485e-5, 2e-9) && nearK(bRho06, 754.629, 0.05),
        std::to_string(bRho06));
    emit("b_costald_tr_0_7",
        b07Ok && nearK(bVs07, 6.60305e-5, 2e-9) && nearK(bRho07, 697.693, 0.05),
        std::to_string(bRho07));
    emit("b_costald_tr_0_8",
        b08Ok && nearK(bVs08, 7.29471e-5, 2e-9) && nearK(bRho08, 631.540, 0.05),
        std::to_string(bRho08));
    emit("b_costald_tr_0_9",
        b09Ok && nearK(bVs09, 8.44329e-5, 2e-9) && nearK(bRho09, 545.628, 0.05),
        std::to_string(bRho09));
    emit("a_sampled_density_decreases_with_t",
        a06Ok && a07Ok && a08Ok && a09Ok
            && aRho06 > aRho07 && aRho07 > aRho08 && aRho08 > aRho09
            && aVs06 < aVs07 && aVs07 < aVs08 && aVs08 < aVs09, "");
    emit("b_sampled_density_decreases_with_t",
        b06Ok && b07Ok && b08Ok && b09Ok
            && bRho06 > bRho07 && bRho07 > bRho08 && bRho08 > bRho09
            && bVs06 < bVs07 && bVs07 < bVs08 && bVs08 < bVs09, "");

    double rejectedVs = 0.0, rejectedRho = 0.0;
    emit("costald_rejects_tr_le_0_25",
        !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA,
            kSaceCostaldTrMinExclusive * densA.criticalTemperatureK, rejectedVs), "");
    emit("costald_rejects_tr_ge_0_95",
        !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA,
            kSaceCostaldTrMaxExclusive * densA.criticalTemperatureK, rejectedVs), "");
    emit("costald_rejects_nonpositive_t",
        !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA, 0.0, rejectedVs)
            && !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA, -1.0, rejectedVs), "");
    emit("costald_rejects_nan_inf",
        !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA, nanH, rejectedVs)
            && !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(densA, infH, rejectedVs), "");
    SaceCostaldLiquidDensityModel invalidDens{};
    emit("costald_rejects_invalid_model",
        !saceCostaldSaturatedLiquidMolarVolumeM3PerMol(invalidDens, 298.15, rejectedVs), "");
    SaceCostaldLiquidDensityModel badMass = densA;
    badMass.molarMassKgPerMol = 0.0;
    emit("costald_density_rejects_invalid_molar_mass",
        !saceCostaldSaturatedLiquidDensityKgPerM3(badMass, 298.15, rejectedRho)
            && rejectedRho == 0.0, "");

    SaceCostaldLiquidDensityModel permDens{}, bondDens{}, endDens{};
    double permRho = 0, bondRho = 0, endRho = 0;
    emit("atom_permutation_preserves_costald",
        saceBuildCostaldLiquidDensityModelFromJoback(permuteGraph(gA), permB, permDens)
            && saceCostaldSaturatedLiquidDensityKgPerM3(permDens, 298.15, permRho)
            && nearK(permRho, aRho298, 1e-8), "");
    emit("bond_reorder_preserves_costald",
        saceBuildCostaldLiquidDensityModelFromJoback(reorderBonds(gA), bondB, bondDens)
            && saceCostaldSaturatedLiquidDensityKgPerM3(bondDens, 298.15, bondRho)
            && nearK(bondRho, aRho298, 1e-8), "");
    emit("reversed_endpoints_preserve_costald",
        saceBuildCostaldLiquidDensityModelFromJoback(reverseEnds(gA), endB, endDens)
            && saceCostaldSaturatedLiquidDensityKgPerM3(endDens, 298.15, endRho)
            && nearK(endRho, aRho298, 1e-8), "");

    emit("catalog_a_density_matches_direct",
        recA && recA->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && nearK(recA->properties.saturatedLiquidDensityAt298KKgPerM3.value, aRho298, 0.05)
            && nearK(bunA.saturatedLiquidDensityAt298KKgPerM3.value, aRho298, 0.05), "");
    emit("catalog_b_density_matches_direct",
        recB && recB->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && nearK(recB->properties.saturatedLiquidDensityAt298KKgPerM3.value, bRho298, 0.05)
            && nearK(bunB.saturatedLiquidDensityAt298KKgPerM3.value, bRho298, 0.05), "");
    emit("density_structural_estimate_low",
        bunA.saturatedLiquidDensityAt298KKgPerM3.source == SacePropertySource::StructuralEstimate
            && bunA.saturatedLiquidDensityAt298KKgPerM3.confidence == SaceConfidence::Low
            && recA && recA->properties.saturatedLiquidDensityAt298KKgPerM3.source
                == SacePropertySource::StructuralEstimate
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.confidence == SaceConfidence::Low, "");
    emit("water_no_costald_density", !waterB.saturatedLiquidDensityAt298KKgPerM3.known, "");
    emit("methane_no_costald_density", !methaneB.saturatedLiquidDensityAt298KKgPerM3.known, "");
    emit("density_not_in_canonical_identity",
        std::strstr(sig, "821") == nullptr && std::strstr(sig, "density") == nullptr, sig);
    emit("costald_not_in_canonical_identity",
        std::strstr(sig, "COSTALD") == nullptr && std::strstr(sig, "Vstar") == nullptr, sig);

    SaceGeneratedRecord *mutRho = cat.recordMutable(refA.generatedId);
    if (mutRho)
        mutRho->properties.saturatedLiquidDensityAt298KKgPerM3 = saceUnknownScalarProperty();
    bool reattachRho = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_density",
        reattachRho && recA && recA->properties.saturatedLiquidDensityAt298KKgPerM3.known
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.source
                == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.saturatedLiquidDensityAt298KKgPerM3.value, aRho298, 0.05), "");

    SaceGeneratedRecord *mutRhoRef = cat.recordMutable(refA.generatedId);
    if (mutRhoRef) {
        mutRhoRef->properties.saturatedLiquidDensityAt298KKgPerM3.known = true;
        mutRhoRef->properties.saturatedLiquidDensityAt298KKgPerM3.value = 789.0f;
        mutRhoRef->properties.saturatedLiquidDensityAt298KKgPerM3.source = SacePropertySource::Reference;
        mutRhoRef->properties.saturatedLiquidDensityAt298KKgPerM3.confidence = SaceConfidence::High;
    }
    bool reattachRhoKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_density",
        reattachRhoKeep && recA
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.source == SacePropertySource::Reference
            && recA->properties.saturatedLiquidDensityAt298KKgPerM3.confidence == SaceConfidence::High
            && nearK(recA->properties.saturatedLiquidDensityAt298KKgPerM3.value, 789.0, 0.05), "");
    emit("generated_record_still_unspawnable_after_density", recA && !recA->spawnable, "");

    SaceJobackLiquidViscosityModel visA{}, visB{};
    emit("a_viscosity_model_builds",
        saceBuildJobackLiquidViscosityModel(gA, bunA, visA) && visA.valid, "");
    emit("b_viscosity_model_builds",
        saceBuildJobackLiquidViscosityModel(gB, bunB, visB) && visB.valid, "");
    emit("a_viscosity_a_2218_35", nearK(visA.aKelvin, 2218.35, 1e-4), std::to_string(visA.aKelvin));
    emit("a_viscosity_b_-18_177", nearK(visA.bDimensionless, -18.177, 1e-4),
        std::to_string(visA.bDimensionless));
    emit("b_viscosity_a_620_85", nearK(visB.aKelvin, 620.85, 1e-4), std::to_string(visB.aKelvin));
    emit("b_viscosity_b_-15_026", nearK(visB.bDimensionless, -15.026, 1e-4),
        std::to_string(visB.bDimensionless));
    emit("a_viscosity_mw_46_069", nearK(visA.molarMassGPerMol, 46.069, 1e-4),
        std::to_string(visA.molarMassGPerMol));
    emit("b_viscosity_mw_46_069", nearK(visB.molarMassGPerMol, 46.069, 1e-4),
        std::to_string(visB.molarMassGPerMol));

    double aMu298 = 0.0, bMu298 = 0.0, aMu300 = 0.0, bMu300 = 0.0;
    double aMuTb = 0.0, bMuTb = 0.0, aMu350 = 0.0, aMu400 = 0.0;
    double bMu350 = 0.0, bMu400 = 0.0;
    bool a298Ok = saceJobackLiquidViscosityPaS(visA, 298.15, aMu298);
    bool b298Ok = saceJobackLiquidViscosityPaS(visB, 298.15, bMu298);
    bool a300Ok = saceJobackLiquidViscosityPaS(visA, 300.0, aMu300);
    bool b300Ok = saceJobackLiquidViscosityPaS(visB, 300.0, bMu300);
    bool aMuTbOk = bunA.normalBoilingPointK.known
        && saceJobackLiquidViscosityPaS(visA, bunA.normalBoilingPointK.value, aMuTb);
    bool bMuTbOk = bunB.normalBoilingPointK.known
        && saceJobackLiquidViscosityPaS(visB, bunB.normalBoilingPointK.value, bMuTb);
    bool a350Ok = saceJobackLiquidViscosityPaS(visA, 350.0, aMu350);
    bool a400Ok = saceJobackLiquidViscosityPaS(visA, 400.0, aMu400);
    bool b350Ok = saceJobackLiquidViscosityPaS(visB, 350.0, bMu350);
    bool b400Ok = saceJobackLiquidViscosityPaS(visB, 400.0, bMu400);
    emit("a_mu_298_15", a298Ok && nearK(aMu298, 0.001001279, 2e-9), std::to_string(aMu298));
    emit("b_mu_298_15", b298Ok && nearK(bMu298, 0.0001101664, 2e-9), std::to_string(bMu298));
    emit("a_mu_300", a300Ok && nearK(aMu300, 0.000956376, 2e-9), std::to_string(aMu300));
    emit("b_mu_300", b300Ok && nearK(bMu300, 0.000108761, 2e-9), std::to_string(bMu300));
    emit("a_mu_tb", aMuTbOk && nearK(aMuTb, 0.000420213, 2e-9), std::to_string(aMuTb));
    emit("b_mu_tb", bMuTbOk && nearK(bMuTb, 0.000139513, 2e-9), std::to_string(bMuTb));
    emit("a_sampled_viscosity_decreases_with_t",
        a298Ok && a350Ok && a400Ok
            && aMu298 > aMu350 && aMu350 > aMu400
            && nearK(aMu350, 3.32551e-4, 2e-9) && nearK(aMu400, 1.50585e-4, 2e-9), "");
    emit("b_sampled_viscosity_decreases_with_t",
        b298Ok && b350Ok && b400Ok
            && bMu298 > bMu350 && bMu350 > bMu400
            && nearK(bMu350, 8.09238e-5, 2e-9) && nearK(bMu400, 6.48305e-5, 2e-9), "");

    SaceJobackEstimateBundle isoBun{}, neoBun{};
    SaceJobackLiquidViscosityModel isoVis{}, neoVis{};
    bool isoVisOk = saceEstimateJobackBundle(isobutaneGraph(), isoBun)
        && saceBuildJobackLiquidViscosityModel(isobutaneGraph(), isoBun, isoVis);
    bool neoVisOk = saceEstimateJobackBundle(neopentaneGraph(), neoBun)
        && saceBuildJobackLiquidViscosityModel(neopentaneGraph(), neoBun, neoVis);
    emit("isobutane_viscosity_a_724_90",
        isoVisOk && nearK(isoVis.aKelvin, 724.90, 1e-2), std::to_string(isoVis.aKelvin));
    emit("isobutane_viscosity_b_-15_172",
        isoVisOk && nearK(isoVis.bDimensionless, -15.172, 1e-3),
        std::to_string(isoVis.bDimensionless));
    emit("neopentane_viscosity_a_1021_78",
        neoVisOk && nearK(neoVis.aKelvin, 1021.78, 1e-2), std::to_string(neoVis.aKelvin));
    emit("neopentane_viscosity_b_-15_771",
        neoVisOk && nearK(neoVis.bDimensionless, -15.771, 1e-3),
        std::to_string(neoVis.bDimensionless));

    double rejectedMu = 0.0;
    emit("viscosity_rejects_nonpositive_t",
        !saceJobackLiquidViscosityPaS(visA, 0.0, rejectedMu)
            && !saceJobackLiquidViscosityPaS(visA, -1.0, rejectedMu)
            && rejectedMu == 0.0, "");
    emit("viscosity_rejects_nan_inf",
        !saceJobackLiquidViscosityPaS(visA, nanH, rejectedMu)
            && !saceJobackLiquidViscosityPaS(visA, infH, rejectedMu), "");
    SaceJobackLiquidViscosityModel invalidVis{};
    emit("viscosity_rejects_invalid_model",
        !saceJobackLiquidViscosityPaS(invalidVis, 298.15, rejectedMu), "");
    SaceJobackLiquidViscosityModel badMw = visA;
    badMw.molarMassGPerMol = 0.0;
    emit("viscosity_rejects_invalid_mw",
        !saceJobackLiquidViscosityPaS(badMw, 298.15, rejectedMu) && rejectedMu == 0.0, "");

    emit("catalog_a_viscosity_matches_direct",
        recA && recA->properties.liquidDynamicViscosityAt298KPaS.known
            && nearK(recA->properties.liquidDynamicViscosityAt298KPaS.value, aMu298, 2e-8)
            && nearK(bunA.liquidDynamicViscosityAt298KPaS.value, aMu298, 2e-8), "");
    emit("catalog_b_viscosity_matches_direct",
        recB && recB->properties.liquidDynamicViscosityAt298KPaS.known
            && nearK(recB->properties.liquidDynamicViscosityAt298KPaS.value, bMu298, 2e-8)
            && nearK(bunB.liquidDynamicViscosityAt298KPaS.value, bMu298, 2e-8), "");
    emit("viscosity_structural_estimate_low",
        bunA.liquidDynamicViscosityAt298KPaS.source == SacePropertySource::StructuralEstimate
            && bunA.liquidDynamicViscosityAt298KPaS.confidence == SaceConfidence::Low
            && recA && recA->properties.liquidDynamicViscosityAt298KPaS.source
                == SacePropertySource::StructuralEstimate
            && recA->properties.liquidDynamicViscosityAt298KPaS.confidence == SaceConfidence::Low, "");

    SaceJobackLiquidViscosityModel permVis{}, bondVis{}, endVis{};
    double permMu = 0, bondMu = 0, endMu = 0;
    emit("atom_permutation_preserves_viscosity",
        saceBuildJobackLiquidViscosityModel(permuteGraph(gA), permB, permVis)
            && saceJobackLiquidViscosityPaS(permVis, 298.15, permMu)
            && nearK(permMu, aMu298, 1e-12), "");
    emit("bond_reorder_preserves_viscosity",
        saceBuildJobackLiquidViscosityModel(reorderBonds(gA), bondB, bondVis)
            && saceJobackLiquidViscosityPaS(bondVis, 298.15, bondMu)
            && nearK(bondMu, aMu298, 1e-12), "");
    emit("reversed_endpoints_preserve_viscosity",
        saceBuildJobackLiquidViscosityModel(reverseEnds(gA), endB, endVis)
            && saceJobackLiquidViscosityPaS(endVis, 298.15, endMu)
            && nearK(endMu, aMu298, 1e-12), "");

    emit("water_no_generated_joback_viscosity", !waterB.liquidDynamicViscosityAt298KPaS.known, "");
    emit("methane_no_generated_joback_viscosity", !methaneB.liquidDynamicViscosityAt298KPaS.known, "");
    emit("viscosity_not_in_canonical_identity",
        std::strstr(sig, "viscosity") == nullptr && std::strstr(sig, "0.00100") == nullptr, sig);
    emit("viscosity_model_not_in_canonical_identity",
        std::strstr(sig, "JobackVisc") == nullptr && std::strstr(sig, "mu_liq") == nullptr, sig);

    SaceGeneratedRecord *mutVisc = cat.recordMutable(refA.generatedId);
    if (mutVisc)
        mutVisc->properties.liquidDynamicViscosityAt298KPaS = saceUnknownScalarProperty();
    bool reattachVisc = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_viscosity",
        reattachVisc && recA && recA->properties.liquidDynamicViscosityAt298KPaS.known
            && recA->properties.liquidDynamicViscosityAt298KPaS.source
                == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.liquidDynamicViscosityAt298KPaS.value, aMu298, 2e-8), "");

    SaceGeneratedRecord *mutViscRef = cat.recordMutable(refA.generatedId);
    if (mutViscRef) {
        mutViscRef->properties.liquidDynamicViscosityAt298KPaS.known = true;
        mutViscRef->properties.liquidDynamicViscosityAt298KPaS.value = 0.0012f;
        mutViscRef->properties.liquidDynamicViscosityAt298KPaS.source = SacePropertySource::Reference;
        mutViscRef->properties.liquidDynamicViscosityAt298KPaS.confidence = SaceConfidence::High;
    }
    bool reattachViscKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_viscosity",
        reattachViscKeep && recA
            && recA->properties.liquidDynamicViscosityAt298KPaS.source == SacePropertySource::Reference
            && recA->properties.liquidDynamicViscosityAt298KPaS.confidence == SaceConfidence::High
            && nearK(recA->properties.liquidDynamicViscosityAt298KPaS.value, 0.0012, 1e-8), "");
    emit("generated_record_still_unspawnable_after_viscosity", recA && !recA->spawnable, "");

    SaceSastriRaoSurfaceTensionModel stA{}, stB{};
    bool stAOk = saceBuildSastriRaoSurfaceTensionModelFromJoback(gA, bunA, stA);
    bool stBOk = saceBuildSastriRaoSurfaceTensionModelFromJoback(gB, bunB, stB);
    emit("a_sastri_rao_model_builds", stAOk && stA.valid, "");
    emit("b_sastri_rao_model_builds", stBOk && stB.valid, "");
    emit("a_surface_tension_class_alcohol",
        stA.chemicalClass == SaceSurfaceTensionClass::Alcohol, "");
    emit("b_surface_tension_class_general_organic",
        stB.chemicalClass == SaceSurfaceTensionClass::GeneralOrganic, "");
    emit("a_sastri_rao_tb_tc_pc",
        stAOk && nearK(stA.normalBoilingPointK, 337.540009, 0.02)
            && nearK(stA.criticalTemperatureK, 499.407379, 0.05)
            && nearK(stA.criticalPressurePa, 5756641.5, 80.0),
        std::to_string(stA.criticalPressurePa));
    emit("b_sastri_rao_tb_tc_pc",
        stBOk && nearK(stB.normalBoilingPointK, 267.779999, 0.02)
            && nearK(stB.criticalTemperatureK, 428.173981, 0.05)
            && nearK(stB.criticalPressurePa, 4910798.0, 80.0),
        std::to_string(stB.criticalPressurePa));

    double aSig298 = 0.0, bSig298 = 0.0;
    bool aSig298Ok = saceSastriRaoSurfaceTensionNPerM(stA, 298.15, aSig298);
    bool bSig298Ok = saceSastriRaoSurfaceTensionNPerM(stB, 298.15, bSig298);
    emit("a_sigma_298_15", aSig298Ok && nearK(aSig298, 0.0207066, 2e-7), std::to_string(aSig298));
    emit("b_sigma_298_15", bSig298Ok && nearK(bSig298, 0.0144424, 2e-7), std::to_string(bSig298));
    emit("a_sigma_greater_than_b_at_298_15",
        aSig298Ok && bSig298Ok && aSig298 > bSig298, "");

    double aTr06 = 0, aTr07 = 0, aTr08 = 0, aTr09 = 0, aTr10 = 0;
    double bTr06 = 0, bTr07 = 0, bTr08 = 0, bTr09 = 0, bTr10 = 0;
    bool aSig06Ok = saceSastriRaoSurfaceTensionNPerM(stA, 0.6 * stA.criticalTemperatureK, aTr06);
    bool aSig07Ok = saceSastriRaoSurfaceTensionNPerM(stA, 0.7 * stA.criticalTemperatureK, aTr07);
    bool aSig08Ok = saceSastriRaoSurfaceTensionNPerM(stA, 0.8 * stA.criticalTemperatureK, aTr08);
    bool aSig09Ok = saceSastriRaoSurfaceTensionNPerM(stA, 0.9 * stA.criticalTemperatureK, aTr09);
    bool aSig10Ok = saceSastriRaoSurfaceTensionNPerM(stA, stA.criticalTemperatureK, aTr10);
    bool bSig06Ok = saceSastriRaoSurfaceTensionNPerM(stB, 0.6 * stB.criticalTemperatureK, bTr06);
    bool bSig07Ok = saceSastriRaoSurfaceTensionNPerM(stB, 0.7 * stB.criticalTemperatureK, bTr07);
    bool bSig08Ok = saceSastriRaoSurfaceTensionNPerM(stB, 0.8 * stB.criticalTemperatureK, bTr08);
    bool bSig09Ok = saceSastriRaoSurfaceTensionNPerM(stB, 0.9 * stB.criticalTemperatureK, bTr09);
    bool bSig10Ok = saceSastriRaoSurfaceTensionNPerM(stB, stB.criticalTemperatureK, bTr10);
    emit("a_sigma_tr_0_6", aSig06Ok && nearK(aTr06, 0.0205835, 2e-7), std::to_string(aTr06));
    emit("a_sigma_tr_0_7", aSig07Ok && nearK(aTr07, 0.0163519, 2e-7), std::to_string(aTr07));
    emit("a_sigma_tr_0_8", aSig08Ok && nearK(aTr08, 0.0118221, 2e-7), std::to_string(aTr08));
    emit("a_sigma_tr_0_9", aSig09Ok && nearK(aTr09, 0.00679001, 2e-7), std::to_string(aTr09));
    emit("b_sigma_tr_0_6", bSig06Ok && nearK(bTr06, 0.0202249, 2e-7), std::to_string(bTr06));
    emit("b_sigma_tr_0_7", bSig07Ok && nearK(bTr07, 0.0142293, 2e-7), std::to_string(bTr07));
    emit("b_sigma_tr_0_8", bSig08Ok && nearK(bTr08, 0.00866885, 2e-7), std::to_string(bTr08));
    emit("b_sigma_tr_0_9", bSig09Ok && nearK(bTr09, 0.00371566, 2e-7), std::to_string(bTr09));
    emit("a_sigma_tc_exactly_zero", aSig10Ok && aTr10 == 0.0, std::to_string(aTr10));
    emit("b_sigma_tc_exactly_zero", bSig10Ok && bTr10 == 0.0, std::to_string(bTr10));

    double rejectedSig = 0.0;
    emit("surface_tension_rejects_t_gt_tc",
        !saceSastriRaoSurfaceTensionNPerM(stA, stA.criticalTemperatureK + 1.0, rejectedSig), "");
    emit("surface_tension_rejects_nonpositive_t",
        !saceSastriRaoSurfaceTensionNPerM(stA, 0.0, rejectedSig)
            && !saceSastriRaoSurfaceTensionNPerM(stA, -1.0, rejectedSig)
            && rejectedSig == 0.0, "");
    emit("surface_tension_rejects_nan_inf",
        !saceSastriRaoSurfaceTensionNPerM(stA, nanH, rejectedSig)
            && !saceSastriRaoSurfaceTensionNPerM(stA, infH, rejectedSig), "");
    SaceSastriRaoSurfaceTensionModel invalidSt{};
    emit("surface_tension_rejects_invalid_model",
        !saceSastriRaoSurfaceTensionNPerM(invalidSt, 298.15, rejectedSig), "");
    emit("a_sampled_sigma_decreases_toward_tc",
        aSig06Ok && aSig07Ok && aSig08Ok && aSig09Ok && aSig10Ok
            && aTr06 > aTr07 && aTr07 > aTr08 && aTr08 > aTr09 && aTr09 > aTr10
            && aTr10 == 0.0, "");
    emit("b_sampled_sigma_decreases_toward_tc",
        bSig06Ok && bSig07Ok && bSig08Ok && bSig09Ok && bSig10Ok
            && bTr06 > bTr07 && bTr07 > bTr08 && bTr08 > bTr09 && bTr09 > bTr10
            && bTr10 == 0.0, "");

    SaceSastriRaoSurfaceTensionModel permSt{}, bondSt{}, endSt{};
    double permSig = 0, bondSig = 0, endSig = 0;
    emit("atom_permutation_preserves_surface_tension",
        saceBuildSastriRaoSurfaceTensionModelFromJoback(permuteGraph(gA), permB, permSt)
            && permSt.chemicalClass == SaceSurfaceTensionClass::Alcohol
            && saceSastriRaoSurfaceTensionNPerM(permSt, 298.15, permSig)
            && nearK(permSig, aSig298, 1e-12), "");
    emit("bond_reorder_preserves_surface_tension",
        saceBuildSastriRaoSurfaceTensionModelFromJoback(reorderBonds(gA), bondB, bondSt)
            && bondSt.chemicalClass == SaceSurfaceTensionClass::Alcohol
            && saceSastriRaoSurfaceTensionNPerM(bondSt, 298.15, bondSig)
            && nearK(bondSig, aSig298, 1e-12), "");
    emit("reversed_endpoints_preserve_surface_tension",
        saceBuildSastriRaoSurfaceTensionModelFromJoback(reverseEnds(gA), endB, endSt)
            && endSt.chemicalClass == SaceSurfaceTensionClass::Alcohol
            && saceSastriRaoSurfaceTensionNPerM(endSt, 298.15, endSig)
            && nearK(endSig, aSig298, 1e-12), "");

    emit("water_no_generated_sastri_rao_surface_tension",
        !waterB.liquidSurfaceTensionAt298KNPerM.known, "");
    emit("methane_no_generated_surface_tension",
        !methaneB.liquidSurfaceTensionAt298KNPerM.known, "");
    emit("surface_tension_not_in_canonical_identity",
        std::strstr(sig, "surface") == nullptr && std::strstr(sig, "0.0207") == nullptr, sig);
    emit("surface_tension_model_not_in_canonical_identity",
        std::strstr(sig, "Sastri") == nullptr && std::strstr(sig, "sigma") == nullptr, sig);

    emit("catalog_a_surface_tension_matches_direct",
        recA && recA->properties.liquidSurfaceTensionAt298KNPerM.known
            && nearK(recA->properties.liquidSurfaceTensionAt298KNPerM.value, aSig298, 2e-7)
            && nearK(bunA.liquidSurfaceTensionAt298KNPerM.value, aSig298, 2e-7), "");
    emit("catalog_b_surface_tension_matches_direct",
        recB && recB->properties.liquidSurfaceTensionAt298KNPerM.known
            && nearK(recB->properties.liquidSurfaceTensionAt298KNPerM.value, bSig298, 2e-7)
            && nearK(bunB.liquidSurfaceTensionAt298KNPerM.value, bSig298, 2e-7), "");

    SaceGeneratedRecord *mutSt = cat.recordMutable(refA.generatedId);
    if (mutSt)
        mutSt->properties.liquidSurfaceTensionAt298KNPerM = saceUnknownScalarProperty();
    bool reattachSt = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_reattach_backfills_surface_tension",
        reattachSt && recA && recA->properties.liquidSurfaceTensionAt298KNPerM.known
            && recA->properties.liquidSurfaceTensionAt298KNPerM.source
                == SacePropertySource::StructuralEstimate
            && nearK(recA->properties.liquidSurfaceTensionAt298KNPerM.value, aSig298, 2e-7), "");

    SaceGeneratedRecord *mutStRef = cat.recordMutable(refA.generatedId);
    if (mutStRef) {
        mutStRef->properties.liquidSurfaceTensionAt298KNPerM.known = true;
        mutStRef->properties.liquidSurfaceTensionAt298KNPerM.value = 0.022f;
        mutStRef->properties.liquidSurfaceTensionAt298KNPerM.source = SacePropertySource::Reference;
        mutStRef->properties.liquidSurfaceTensionAt298KNPerM.confidence = SaceConfidence::High;
    }
    bool reattachStKeep = cat.attachMolecularGraph(refA.generatedId, "synthetic-structure-a", gA);
    recA = cat.record(refA.generatedId);
    emit("same_graph_backfill_preserves_reference_high_surface_tension",
        reattachStKeep && recA
            && recA->properties.liquidSurfaceTensionAt298KNPerM.source == SacePropertySource::Reference
            && recA->properties.liquidSurfaceTensionAt298KNPerM.confidence == SaceConfidence::High
            && nearK(recA->properties.liquidSurfaceTensionAt298KNPerM.value, 0.022, 1e-8), "");
    emit("generated_record_still_unspawnable_after_surface_tension", recA && !recA->spawnable, "");

    cat.clear();
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
