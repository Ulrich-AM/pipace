#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

// Thermal units (internally consistent SI on a 2D grid with implied depth):
//   Temperature: Kelvin
//   Specific heat: J/(kg·K)
//   Conductivity: W/(m·K)  — physical; sandbox scale is ThermalConfig::conductivityScale
//   Energy: Joules
//
// Mapping: fluid.config.cellsPerMeter (default 4) → cell edge 0.25 m.
// Implied depth = one cell, so cell volume = dx³ = 0.015625 m³.
// Mass of a full water cell ≈ 15.625 kg (densityRel 1.0 × 1000 kg/m³ × V).
//
// Real metal would take thousands of seconds to equalize across 0.25 m. Gameplay uses
// conductivityScale (~2.5e4) so relative rates stay physical (metal ≫ stone/glass ≫ wood
// ≫ water ≫ air) while conduction is visible in seconds. Do not treat scaled k as
// laboratory SI in UI copy.
//
// Gas pressure remains isothermal (amount/volume). Path for P∝T/T_amb is documented
// in ThermalEngine; it is not enabled in this update.
//
// Melting/boiling/latent fields are metadata only — no phase change yet.

constexpr float AMBIENT_TEMPERATURE_K = 293.15f;
constexpr float WATER_DENSITY_KG_M3 = 1000.0f;
constexpr float AIR_DENSITY_KG_M3 = 1.204f;
constexpr float MIN_THERMAL_CAPACITY = 1.0e-6f; // J/K
constexpr float MIN_THERMAL_MASS_KG = 1.0e-9f;
constexpr float MIN_SAFE_TEMPERATURE_K = 0.05f;
constexpr float MAX_SAFE_TEMPERATURE_K = 1.0e7f;

struct ThermalProperties {
    float specificHeat = 1000.0f;   // J/(kg·K)
    float conductivity = 0.026f;    // W/(m·K) physical
    // Future-use only (phase-change update). Solver must ignore these.
    float meltingPointK = 0.0f;
    float boilingPointK = 0.0f;
    float latentFusion = 0.0f;      // J/kg
    float latentVapor = 0.0f;       // J/kg
    // Hooks only: thermal expansion does not resize bodies; softening/melting do not run.
    float expansionCoeff = 0.0f;    // 1/K (linear solids, volumetric fluids)
    float softeningTempK = 0.0f;
};

// Approximate 20 °C engineering-table values, not laboratory certificates.
// Relative behavior is the point: metal ≫ stone/glass ≫ wood ≫ water ≫ air.
inline ThermalProperties const &kWaterThermal() {
    // Water at ~20 °C: cp ≈ 4184 J/(kg·K), k ≈ 0.598 W/(m·K).
    // expansionCoeff is volumetric (~2.07e-4 /K near 20 °C). The 4 °C density
    // maximum is NOT modeled — see LiquidProperties::densityAtTemperature.
    static ThermalProperties const p{4184.0f, 0.598f, 273.15f, 373.15f, 3.34e5f, 2.26e6f, 2.07e-4f, 0.0f};
    return p;
}
inline ThermalProperties const &kHoneyThermal() {
    // Kitchen honey ballpark: denser, lower cp than water, similar k. Not a lab curve.
    static ThermalProperties const p{2200.0f, 0.50f, 0.0f, 0.0f, 0.0f, 0.0f, 2.07e-4f, 0.0f};
    return p;
}
inline ThermalProperties const &kAirThermal() {
    // Dry air ~20 °C, 1 atm: cp ≈ 1005 J/(kg·K), k ≈ 0.026 W/(m·K).
    // expansionCoeff ≈ 1/T_amb for an ideal gas; unused by the isothermal pressure solver.
    static ThermalProperties const p{1005.0f, 0.026f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f / AMBIENT_TEMPERATURE_K, 0.0f};
    return p;
}
inline ThermalProperties const &kWoodThermal() {
    // Dry softwood ballpark: cp ~1700, k ~0.12, α ~5e-6 /K. Softening ~450 K is pyrolysis-ish metadata.
    static ThermalProperties const p{1700.0f, 0.12f, 0.0f, 0.0f, 0.0f, 0.0f, 5.0e-6f, 450.0f};
    return p;
}
inline ThermalProperties const &kStoneThermal() {
    // Granite-like: cp ~880, k ~1.7, α ~8e-6 /K, melt ~1200 °C.
    static ThermalProperties const p{880.0f, 1.70f, 1473.0f, 0.0f, 0.0f, 0.0f, 8.0e-6f, 0.0f};
    return p;
}
inline ThermalProperties const &kGlassThermal() {
    // Soda-lime glass: cp ~840, k ~1.0, α ~9e-6 /K, softening ~800 K, melt ~1700 K.
    static ThermalProperties const p{840.0f, 1.00f, 1700.0f, 0.0f, 0.0f, 0.0f, 9.0e-6f, 800.0f};
    return p;
}
inline ThermalProperties const &kMetalThermal() {
    // Carbon-steel ballpark: cp ~490, k ~50 (between stainless and mild steel), α ~12e-6 /K.
    static ThermalProperties const p{490.0f, 50.0f, 1811.0f, 0.0f, 0.0f, 0.0f, 1.2e-5f, 1000.0f};
    return p;
}
inline ThermalProperties const &kWallThermal() {
    return kStoneThermal();
}
inline ThermalProperties const &kEmptyThermal() {
    static ThermalProperties const p{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    return p;
}

// MaterialId values match rigid/RigidBodyTypes.h (0 empty, 1 wood, 2 stone, 3 glass, 4 metal).
inline ThermalProperties const &thermalForMaterial(uint16_t materialId) {
    switch (materialId) {
        case 1: return kWoodThermal();
        case 2: return kStoneThermal();
        case 3: return kGlassThermal();
        case 4: return kMetalThermal();
        default: return kEmptyThermal();
    }
}

inline float cellLengthM(float cellsPerMeter) {
    float cpm = cellsPerMeter > 0.1f ? cellsPerMeter : 4.0f;
    return 1.0f / cpm;
}
inline float cellVolumeM3(float cellsPerMeter) {
    float dx = cellLengthM(cellsPerMeter);
    return dx * dx * dx;
}
inline float massKg(float densityRelToWater, float fillOrFraction, float cellsPerMeter = 4.0f) {
    float m = densityRelToWater * WATER_DENSITY_KG_M3 * cellVolumeM3(cellsPerMeter) * std::max(0.0f, fillOrFraction);
    return m;
}
inline float gasMassKg(float amount, float cellsPerMeter = 4.0f) {
    return std::max(0.0f, amount) * AIR_DENSITY_KG_M3 * cellVolumeM3(cellsPerMeter);
}
inline float thermalCapacity(float mass, float specificHeat) {
    return std::max(0.0f, mass) * std::max(0.0f, specificHeat);
}
inline float energyFromTemp(float capacity, float temperatureK) {
    if (capacity < MIN_THERMAL_CAPACITY) return 0.0f;
    return capacity * temperatureK;
}
inline float tempFromEnergy(float energy, float capacity, float fallbackK = AMBIENT_TEMPERATURE_K) {
    if (!(capacity > MIN_THERMAL_CAPACITY) || !std::isfinite(energy) || !std::isfinite(capacity))
        return fallbackK;
    float t = energy / capacity;
    if (!std::isfinite(t)) return fallbackK;
    if (t < MIN_SAFE_TEMPERATURE_K) return MIN_SAFE_TEMPERATURE_K;
    if (t > MAX_SAFE_TEMPERATURE_K) return MAX_SAFE_TEMPERATURE_K;
    return t;
}
inline float harmonicConductivity(float ka, float kb) {
    ka = std::max(0.0f, ka);
    kb = std::max(0.0f, kb);
    float s = ka + kb;
    if (s <= 1.0e-20f) return 0.0f;
    return 2.0f * ka * kb / s;
}

// Conserves Ea+Eb. Clamps so one step cannot leap past the two-node equilibrium.
inline float exchangeThermalEnergy(float &energyA, float capacityA, float kA,
    float &energyB, float capacityB, float kB,
    float dt, float areaOverDx, float conductivityScale)
{
    if (capacityA < MIN_THERMAL_CAPACITY || capacityB < MIN_THERMAL_CAPACITY) return 0.0f;
    if (dt <= 0.0f || areaOverDx <= 0.0f) return 0.0f;
    float k = harmonicConductivity(kA, kB) * std::max(0.0f, conductivityScale);
    if (k <= 0.0f) return 0.0f;
    float Ta = tempFromEnergy(energyA, capacityA);
    float Tb = tempFromEnergy(energyB, capacityB);
    float dT = Ta - Tb;
    if (std::abs(dT) < 1.0e-6f) return 0.0f;
    float dQ = k * areaOverDx * dT * dt;
    float Csum = capacityA + capacityB;
    float eqT = (energyA + energyB) / Csum;
    float maxQ = (Ta - eqT) * capacityA;
    if (dQ > 0.0f) dQ = std::min(dQ, std::max(0.0f, maxQ));
    else dQ = std::max(dQ, std::min(0.0f, maxQ));
    energyA -= dQ;
    energyB += dQ;
    if (!std::isfinite(energyA)) energyA = eqT * capacityA;
    if (!std::isfinite(energyB)) energyB = eqT * capacityB;
    if (energyA < 0.0f) energyA = 0.0f;
    if (energyB < 0.0f) energyB = 0.0f;
    return dQ;
}
