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

int collectValidComponents(LiquidComponentView const &composition,
    SubstanceId *ids, float *amts, float &tot)
{
    int n = 0;
    tot = 0.0f;
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

LiquidMixtureProperties fromPureIntrinsic(SubstanceId id) {
    FluidProperties const &f = fluidForSubstance(id);
    ThermalProperties const &th = thermalForSubstance(id);
    LiquidMixtureProperties out;
    out.density = f.density;
    out.specificHeat = th.specificHeat;
    out.conductivity = th.conductivity;
    out.viscosity = 0.0f;
    out.surfaceTension = f.surfaceTension;
    sanitizeIntrinsic(out);
    return out;
}

} // namespace

LiquidMixtureProperties referenceLiquidMixture() {
    return fromPureIntrinsic(SUBSTANCE_WATER);
}

LiquidMixtureProperties evaluateLiquidMixture(LiquidComponentView const &composition) {
    SubstanceId ids[kMaxLiquidComponents];
    float amts[kMaxLiquidComponents];
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
        float phi = amts[i] / tot;
        FluidProperties const &f = fluidForSubstance(ids[i]);
        ThermalProperties const &th = thermalForSubstance(ids[i]);
        rho += phi * f.density;
        k += phi * th.conductivity;
        gamma += phi * f.surfaceTension;
        float m = amts[i] * std::max(0.0f, f.density);
        massCp += m * th.specificHeat;
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
    SubstanceId ids[kMaxLiquidComponents];
    float amts[kMaxLiquidComponents];
    float tot = 0.0f;
    int n = collectValidComponents(composition, ids, amts, tot);
    if (n == 0 || !(tot > kMinLiquidComponent))
        return safeViscosity(sandboxReferenceLiquid().viscosityAtTemperature(T));
    if (n == 1)
        return safeViscosity(fluidForSubstance(ids[0]).viscosityAtTemperature(T));

    float logMu = 0.0f;
    for (int i = 0; i < n; ++i) {
        float phi = amts[i] / tot;
        float mu = fluidForSubstance(ids[i]).viscosityAtTemperature(T);
        logMu += phi * std::log(safeViscosity(mu));
    }
    return safeViscosity(std::exp(logMu));
}

bool liquidCompositionIsPureWater(LiquidComponentView const &composition) {
    float water = 0.0f;
    float other = 0.0f;
    for (int i = 0; i < composition.count && i < kMaxLiquidComponents; ++i) {
        SubstanceId id = composition.items[i].id;
        float amt = composition.items[i].amount;
        if (!validLiquidComponentId(id) || !std::isfinite(amt) || amt <= kMinLiquidComponent)
            continue;
        if (id == SUBSTANCE_WATER) water += amt;
        else other += amt;
    }
    return water > kMinLiquidComponent && other <= kMeaningfulLiquidComponent;
}
