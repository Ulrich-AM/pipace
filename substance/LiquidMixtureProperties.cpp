#include "substance/LiquidMixtureProperties.h"

#include "substance/GeneratedMaterialRegistry.h"

#include <cmath>

namespace {

float sanitizeTemperatureK(float temperatureK) {
    if (!std::isfinite(temperatureK) || temperatureK <= 1.0f)
        return AMBIENT_TEMPERATURE_K;
    return temperatureK;
}

float safeViscosity(float mu) {
    if (!std::isfinite(mu) || mu < kMinMixtureViscosity) return kMinMixtureViscosity;
    return mu;
}

int collectValidComponents(LiquidComponentView const &composition,
    RuntimeSubstanceRef *ids, float *amts, float &tot)
{
    int n = 0;
    tot = 0.0f;
    for (int i = 0; i < composition.count && i < kMaxLiquidComponents; ++i) {
        RuntimeSubstanceRef id = composition.items[i].id;
        float amt = composition.items[i].amount;
        if (!validLiquidComponentId(id) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        ids[n] = id;
        amts[n] = amt;
        tot += amt;
        ++n;
    }
    return n;
}

void sanitizeIntrinsic(LiquidMixtureProperties &out) {
    if (!std::isfinite(out.density) || out.density <= 0.0f)
        out.density = sandboxReferenceLiquid().density;
    if (!std::isfinite(out.specificHeat) || out.specificHeat <= 0.0f)
        out.specificHeat = thermalForSubstance(SUBSTANCE_WATER).specificHeat;
    if (!std::isfinite(out.conductivity) || out.conductivity < 0.0f)
        out.conductivity = thermalForSubstance(SUBSTANCE_WATER).conductivity;
    if (!std::isfinite(out.surfaceTension) || out.surfaceTension < 0.0f)
        out.surfaceTension = sandboxReferenceLiquid().surfaceTension;
}

bool intrinsicFor(RuntimeSubstanceRef id, LiquidMixtureProperties &out) {
    out = {};
    if (runtimeSubstanceIsBuiltIn(id)) {
        SubstanceId sid = runtimeBuiltinId(id);
        if (!validLiquidComponentId(sid)) return false;
        FluidProperties const &f = fluidForSubstance(sid);
        ThermalProperties const &th = thermalForSubstance(sid);
        out.density = f.density;
        out.specificHeat = th.specificHeat;
        out.conductivity = th.conductivity;
        out.viscosity = 0.0f;
        out.surfaceTension = f.surfaceTension;
        sanitizeIntrinsic(out);
        return true;
    }

    SaceCompiledLiquidProfile const *p = runtimeLiquidProfile(id);
    if (!p || !p->valid) return false;
    out.density = p->densityRelativeToWater;
    out.specificHeat = p->specificHeatJPerKgK;
    out.conductivity = p->thermalConductivityWPerMK;
    out.viscosity = 0.0f;
    out.surfaceTension = p->solverSurfaceTension;
    return std::isfinite(out.density) && out.density > 0.0f
        && std::isfinite(out.specificHeat) && out.specificHeat > 0.0f
        && std::isfinite(out.conductivity) && out.conductivity > 0.0f
        && std::isfinite(out.surfaceTension) && out.surfaceTension > 0.0f;
}

float viscosityFor(RuntimeSubstanceRef id, float temperatureK) {
    if (runtimeSubstanceIsBuiltIn(id))
        return safeViscosity(fluidForSubstance(runtimeBuiltinId(id)).viscosityAtTemperature(temperatureK));
    SaceCompiledLiquidProfile const *p = runtimeLiquidProfile(id);
    float mu = 0.0f;
    if (!p || !saceCompiledSolverViscosityAtTemperature(*p, temperatureK, mu))
        return kMinMixtureViscosity;
    return safeViscosity(mu);
}

LiquidMixtureProperties fromPureIntrinsic(RuntimeSubstanceRef id) {
    LiquidMixtureProperties out;
    if (!intrinsicFor(id, out))
        return LiquidMixtureProperties{};
    return out;
}

} // namespace

LiquidMixtureProperties referenceLiquidMixture() {
    return fromPureIntrinsic(runtimeBuiltIn(SUBSTANCE_WATER));
}

LiquidMixtureProperties evaluateLiquidMixture(LiquidComponentView const &composition) {
    RuntimeSubstanceRef ids[kMaxLiquidComponents]{};
    float amts[kMaxLiquidComponents]{};
    float tot = 0.0f;
    int n = collectValidComponents(composition, ids, amts, tot);
    if (n == 0 || !(tot > kMinLiquidComponent))
        return referenceLiquidMixture();
    if (n == 1)
        return fromPureIntrinsic(ids[0]);

    float rho = 0.0f;
    float k = 0.0f;
    float gamma = 0.0f;
    float massCp = 0.0f;
    float massSum = 0.0f;
    for (int i = 0; i < n; ++i) {
        LiquidMixtureProperties p{};
        if (!intrinsicFor(ids[i], p)) continue;
        float phi = amts[i] / tot;
        rho += phi * p.density;
        k += phi * p.conductivity;
        gamma += phi * p.surfaceTension;
        float m = amts[i] * std::max(0.0f, p.density);
        massCp += m * p.specificHeat;
        massSum += m;
    }

    LiquidMixtureProperties out;
    out.density = rho;
    out.conductivity = k;
    out.surfaceTension = gamma;
    out.viscosity = 0.0f;
    out.specificHeat = (massSum > kMinLiquidComponent) ? (massCp / massSum)
        : thermalForSubstance(SUBSTANCE_WATER).specificHeat;
    sanitizeIntrinsic(out);
    return out;
}

float evaluateLiquidMixtureViscosity(LiquidComponentView const &composition, float temperatureK) {
    float T = sanitizeTemperatureK(temperatureK);
    RuntimeSubstanceRef ids[kMaxLiquidComponents]{};
    float amts[kMaxLiquidComponents]{};
    float tot = 0.0f;
    int n = collectValidComponents(composition, ids, amts, tot);
    if (n == 0 || !(tot > kMinLiquidComponent))
        return safeViscosity(sandboxReferenceLiquid().viscosityAtTemperature(T));
    if (n == 1)
        return viscosityFor(ids[0], T);

    float logMu = 0.0f;
    for (int i = 0; i < n; ++i) {
        float phi = amts[i] / tot;
        logMu += phi * std::log(viscosityFor(ids[i], T));
    }
    return safeViscosity(std::exp(logMu));
}

bool liquidCompositionIsPureWater(LiquidComponentView const &composition) {
    float water = 0.0f;
    float other = 0.0f;
    for (int i = 0; i < composition.count && i < kMaxLiquidComponents; ++i) {
        RuntimeSubstanceRef id = composition.items[i].id;
        float amt = composition.items[i].amount;
        if (!validLiquidComponentId(id) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        if (id == SUBSTANCE_WATER) water += amt;
        else other += amt;
    }
    return water > kMinLiquidComponent && other <= kMeaningfulLiquidComponent;
}
