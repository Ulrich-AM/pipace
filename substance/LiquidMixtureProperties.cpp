#include "substance/LiquidMixtureProperties.h"

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

LiquidMixtureProperties fromPureSubstance(SubstanceId id, float temperatureK) {
    FluidProperties const &f = fluidForSubstance(id);
    ThermalProperties const &th = thermalForSubstance(id);
    LiquidMixtureProperties out;
    out.density = f.density;
    out.specificHeat = th.specificHeat;
    out.conductivity = th.conductivity;
    out.viscosity = safeViscosity(f.viscosityAtTemperature(temperatureK));
    out.surfaceTension = f.surfaceTension;
    if (!std::isfinite(out.density) || out.density <= 0.0f)
        out.density = sandboxReferenceLiquid().density;
    if (!std::isfinite(out.specificHeat) || out.specificHeat <= 0.0f)
        out.specificHeat = thermalForSubstance(SUBSTANCE_WATER).specificHeat;
    if (!std::isfinite(out.conductivity) || out.conductivity < 0.0f)
        out.conductivity = thermalForSubstance(SUBSTANCE_WATER).conductivity;
    if (!std::isfinite(out.surfaceTension) || out.surfaceTension < 0.0f)
        out.surfaceTension = sandboxReferenceLiquid().surfaceTension;
    return out;
}

} // namespace

LiquidMixtureProperties referenceLiquidMixture(float temperatureK) {
    return fromPureSubstance(SUBSTANCE_WATER, sanitizeTemperatureK(temperatureK));
}

LiquidMixtureProperties evaluateLiquidMixture(
    LiquidComponentView const &composition,
    float temperatureK)
{
    float T = sanitizeTemperatureK(temperatureK);
    SubstanceId ids[kMaxLiquidComponents];
    float amts[kMaxLiquidComponents];
    int n = 0;
    float tot = 0.0f;
    for (int i = 0; i < composition.count && i < kMaxLiquidComponents; ++i) {
        SubstanceId id = composition.items[i].id;
        float amt = composition.items[i].amount;
        if (!validLiquidComponentId(id) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        ids[n] = id;
        amts[n] = amt;
        tot += amt;
        ++n;
    }
    if (n == 0 || !(tot > kMinLiquidComponent))
        return referenceLiquidMixture(T);
    if (n == 1)
        return fromPureSubstance(ids[0], T);

    float rho = 0.0f;
    float k = 0.0f;
    float gamma = 0.0f;
    float logMu = 0.0f;
    float massCp = 0.0f;
    float massSum = 0.0f;
    for (int i = 0; i < n; ++i) {
        float phi = amts[i] / tot;
        FluidProperties const &f = fluidForSubstance(ids[i]);
        ThermalProperties const &th = thermalForSubstance(ids[i]);
        rho += phi * f.density;
        k += phi * th.conductivity;
        gamma += phi * f.surfaceTension;
        logMu += phi * std::log(safeViscosity(f.viscosityAtTemperature(T)));
        float m = amts[i] * std::max(0.0f, f.density);
        massCp += m * th.specificHeat;
        massSum += m;
    }

    LiquidMixtureProperties out;
    out.density = rho;
    out.conductivity = k;
    out.surfaceTension = gamma;
    out.specificHeat = (massSum > kMinLiquidComponent) ? (massCp / massSum)
        : thermalForSubstance(SUBSTANCE_WATER).specificHeat;
    float mu = std::exp(logMu);
    out.viscosity = safeViscosity(mu);

    if (!std::isfinite(out.density) || out.density <= 0.0f)
        out.density = sandboxReferenceLiquid().density;
    if (!std::isfinite(out.specificHeat) || out.specificHeat <= 0.0f)
        out.specificHeat = thermalForSubstance(SUBSTANCE_WATER).specificHeat;
    if (!std::isfinite(out.conductivity) || out.conductivity < 0.0f)
        out.conductivity = thermalForSubstance(SUBSTANCE_WATER).conductivity;
    if (!std::isfinite(out.surfaceTension) || out.surfaceTension < 0.0f)
        out.surfaceTension = sandboxReferenceLiquid().surfaceTension;
    return out;
}
