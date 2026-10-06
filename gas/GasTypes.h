#pragma once

#include "fluid/FluidTypes.h"

// Gas pressure is a reference-state ideal-gas approximation:
//   P_atm = (amount / volume) * (T / Tref)
// amount is conserved cell-atmospheres at Tref = AMBIENT_TEMPERATURE_K.
// At ambient T the old isothermal P = amount/volume is recovered.
//
// This is one-way coupling: thermal energy -> temperature -> pressure.
// There is no PdV work, no adiabatic compression heating, and no acoustic
// shock solver. Pressure is derived from state and does not modify heat.
// Pressure in the solver is stored in atmospheres; convert with referencePressurePa.
constexpr float GAS_MIN_VOLUME = 1.0e-4f;
constexpr float GAS_MIN_AMOUNT = 1.0e-8f;
constexpr float GAS_REFERENCE_PRESSURE_PA = 101325.0f;
constexpr float GAS_PRESSURE_TREF_K = AMBIENT_TEMPERATURE_K;
// Solver safety cap, atmospheres. Comfortably above brush max (8) and ordinary
// sealed-room / hot-combustion overpressure. Not a gameplay limiter.
constexpr float GAS_MAX_PRESSURE_ATM = 64.0f;

inline float gasPressureAtmFromState(float amount, float volume, float temperatureK) {
    if (!(amount > GAS_MIN_AMOUNT) || !(volume >= GAS_MIN_VOLUME)) return 0.0f;
    if (!std::isfinite(amount) || !std::isfinite(volume)) return 0.0f;
    float T = temperatureK;
    if (!(T > 0.0f) || !std::isfinite(T)) T = GAS_PRESSURE_TREF_K;
    T = std::clamp(T, MIN_SAFE_TEMPERATURE_K, MAX_SAFE_TEMPERATURE_K);
    float tRef = GAS_PRESSURE_TREF_K;
    if (!(tRef > 1.0f) || !std::isfinite(tRef)) tRef = 293.15f;
    float p = (amount / volume) * (T / tRef);
    if (!std::isfinite(p) || p < 0.0f) return 0.0f;
    if (p > GAS_MAX_PRESSURE_ATM) return GAS_MAX_PRESSURE_ATM;
    return p;
}

inline float gasAmountFromPressureAtm(float pressureAtm, float volume, float temperatureK) {
    if (!(volume >= GAS_MIN_VOLUME) || !std::isfinite(volume)) return 0.0f;
    if (!(pressureAtm > 0.0f) || !std::isfinite(pressureAtm)) return 0.0f;
    float T = temperatureK;
    if (!(T > 0.0f) || !std::isfinite(T)) T = GAS_PRESSURE_TREF_K;
    T = std::clamp(T, MIN_SAFE_TEMPERATURE_K, MAX_SAFE_TEMPERATURE_K);
    float tRef = GAS_PRESSURE_TREF_K;
    if (!(tRef > 1.0f) || !std::isfinite(tRef)) tRef = 293.15f;
    float p = std::min(pressureAtm, GAS_MAX_PRESSURE_ATM);
    float a = p * volume * (tRef / T);
    if (!std::isfinite(a) || a < 0.0f) return 0.0f;
    return a;
}

// Fixed-capacity gas composition. amount[] is total cell-atmospheres; these
// slots say what that amount is. 4 covers Air + Water vapor plus two near-term
// extra gases without per-cell heap objects.
constexpr int kMaxGasComponents = 4;
constexpr float kMinGasComponent = GAS_MIN_AMOUNT;

struct GasComponent {
    RuntimeSubstanceRef id{};
    float amount = 0.0f;

    GasComponent() = default;
    GasComponent(RuntimeSubstanceRef runtimeId, float a) : id(runtimeId), amount(a) {}
    GasComponent(SubstanceId builtinId, float a) : id(runtimeBuiltIn(builtinId)), amount(a) {}
};

struct GasComponentView {
    GasComponent items[kMaxGasComponents]{};
    int count = 0;
};

inline bool validGasComponentId(RuntimeSubstanceRef id) {
    if (runtimeSubstanceIsBuiltIn(id)) {
        SubstanceId sid = runtimeBuiltinId(id);
        return validSubstance(sid)
            && sid != SUBSTANCE_NONE
            && supportsPhase(sid, MatterPhase::Gas);
    }
    return runtimeSubstanceIsGenerated(id)
        && runtimeSupportsPhase(id, MatterPhase::Gas);
}

inline bool validGasComponentId(SubstanceId id) {
    return validGasComponentId(runtimeBuiltIn(id));
}

inline int findGasComponent(GasComponent const *items, int count, RuntimeSubstanceRef id) {
    for (int n = 0; n < count; ++n)
        if (items[n].id == id) return n;
    return -1;
}

inline int findGasComponent(GasComponent const *items, int count, SubstanceId id) {
    return findGasComponent(items, count, runtimeBuiltIn(id));
}

inline float gasPayloadAmount(GasComponent const *items, int count, RuntimeSubstanceRef id) {
    int n = findGasComponent(items, count, id);
    return n >= 0 ? items[n].amount : 0.0f;
}

inline float gasPayloadAmount(GasComponent const *items, int count, SubstanceId id) {
    return gasPayloadAmount(items, count, runtimeBuiltIn(id));
}

inline float gasPayloadSum(GasComponent const *items, int count) {
    float s = 0.0f;
    for (int n = 0; n < count; ++n) s += items[n].amount;
    return s;
}

// Merge id into a fixed payload. Returns the amount that did not fit (overflow).
inline float addGasPayload(GasComponent *items, int &count, RuntimeSubstanceRef id, float amount) {
    if (!validGasComponentId(id) || !(amount > kMinGasComponent)) return 0.0f;
    int n = findGasComponent(items, count, id);
    if (n >= 0) {
        items[n].amount += amount;
        return 0.0f;
    }
    if (count >= kMaxGasComponents) return amount;
    items[count++] = {id, amount};
    return 0.0f;
}

inline float addGasPayload(GasComponent *items, int &count, SubstanceId id, float amount) {
    return addGasPayload(items, count, runtimeBuiltIn(id), amount);
}

inline void compactGasPayload(GasComponent *items, int &count) {
    int w = 0;
    for (int n = 0; n < count; ++n) {
        if (items[n].amount > kMinGasComponent && validGasComponentId(items[n].id))
            items[w++] = items[n];
    }
    for (int n = w; n < count; ++n) items[n] = {};
    count = w;
}

inline void scaleGasPayload(GasComponent *items, int count, float frac) {
    for (int n = 0; n < count; ++n) items[n].amount *= frac;
}

inline void copyGasPayload(GasComponent *dst, int &dstCount,
    GasComponent const *src, int srcCount)
{
    dstCount = srcCount;
    for (int n = 0; n < srcCount; ++n) dst[n] = src[n];
    for (int n = srcCount; n < kMaxGasComponents; ++n) dst[n] = {};
}

// True if every occupied src species can merge into dst's fixed slots
// (existing id or a free slot). Used to reject a whole transfer rather than
// drop a species or move identity-less amount.
inline bool gasPayloadCanMerge(GasComponent const *dst, int dstCount,
    GasComponent const *src, int srcCount)
{
    if (!dst || !src || dstCount < 0 || srcCount < 0) return false;
    if (dstCount > kMaxGasComponents) return false;
    int extra = 0;
    for (int s = 0; s < srcCount; ++s) {
        if (!(src[s].amount > kMinGasComponent) || src[s].id == SUBSTANCE_NONE) continue;
        if (findGasComponent(dst, dstCount, src[s].id) >= 0) continue;
        ++extra;
    }
    return dstCount + extra <= kMaxGasComponents;
}

// Compatibility labels only. Authoritative storage is SubstanceId slots:
// Air = SUBSTANCE_AIR, water vapor = SUBSTANCE_WATER + MatterPhase::Gas.
enum class GasSpecies : uint8_t { Air = 0, WaterVapor = 1 };
constexpr int GAS_SPECIES_COUNT = 2;

inline SubstanceId substanceForGasSpecies(GasSpecies species = GasSpecies::Air) {
    return species == GasSpecies::WaterVapor ? SUBSTANCE_WATER : SUBSTANCE_AIR;
}

inline MatterIdentity identityForGasSpecies(GasSpecies species = GasSpecies::Air) {
    return makeMatterIdentity(substanceForGasSpecies(species), MatterPhase::Gas);
}

enum class GasBoundary : uint8_t {
    Sealed = 0,
    OpenAmbient = 1 // infinite 1 atm reservoir; flux is accounted as escaped/entered
};

enum class GasSimMode : uint8_t { Off = 0, Half = 1, Full = 2 };

struct GasTimings {
    double occupancy = 0.0;
    double displace = 0.0;
    double transfer = 0.0;
    double chunks = 0.0;
    double physics = 0.0;
};
