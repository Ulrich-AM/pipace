#pragma once

#include "fluid/FluidTypes.h"

// Pressure remains isothermal in this update: 1.0 amount in 1.0 available cell
// volume = 1 atmosphere. Gas cells also store thermal energy (see thermal/).
// Coupling P ∝ T/T_amb is reserved; ThermalConfig::coupleGasPressureToTemperature
// is off so the existing flow solver stays stable.
// Pressure in the solver is stored in atmospheres; convert with referencePressurePa.
constexpr float GAS_MIN_VOLUME = 1.0e-4f;
constexpr float GAS_MIN_AMOUNT = 1.0e-8f;
constexpr float GAS_REFERENCE_PRESSURE_PA = 101325.0f;

// Fixed-capacity gas composition. amount[] is total cell-atmospheres; these
// slots say what that amount is. 4 covers Air + Water vapor plus two near-term
// extra gases without per-cell heap objects.
constexpr int kMaxGasComponents = 4;
constexpr float kMinGasComponent = GAS_MIN_AMOUNT;

struct GasComponent {
    SubstanceId id = SUBSTANCE_NONE;
    float amount = 0.0f;
};

struct GasComponentView {
    GasComponent items[kMaxGasComponents]{};
    int count = 0;
};

inline bool validGasComponentId(SubstanceId id) {
    return validSubstance(id)
        && id != SUBSTANCE_NONE
        && supportsPhase(id, MatterPhase::Gas);
}

inline int findGasComponent(GasComponent const *items, int count, SubstanceId id) {
    for (int n = 0; n < count; ++n)
        if (items[n].id == id) return n;
    return -1;
}

inline float gasPayloadAmount(GasComponent const *items, int count, SubstanceId id) {
    int n = findGasComponent(items, count, id);
    return n >= 0 ? items[n].amount : 0.0f;
}

inline float gasPayloadSum(GasComponent const *items, int count) {
    float s = 0.0f;
    for (int n = 0; n < count; ++n) s += items[n].amount;
    return s;
}

// Merge id into a fixed payload. Returns the amount that did not fit (overflow).
inline float addGasPayload(GasComponent *items, int &count, SubstanceId id, float amount) {
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
