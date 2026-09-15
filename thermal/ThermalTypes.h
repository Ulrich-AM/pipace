#pragma once

#include "substance/SubstanceTypes.h"

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
// ≫ water ≫ air) while conduction is visible in seconds. Gas uses a separate
// gasConductivityScale (~400): enough that a hot wall warms adjacent air in seconds,
// far below the solid scale. Hot-air transport is mainly buoyancy-driven convection.
//
// Gas pressure remains isothermal (amount/volume). Path for P∝T/T_amb is documented
// in ThermalEngine; it is not enabled in this update.
//
// Melting/boiling/latent live on PhaseProperties (copied from thermal authoring
// at registry build). No phase-change solver yet.
// Authoritative heat tables: thermalForSubstance(id) / substanceDef(id).thermal.

// ThermalProperties and ambient/density constants live in substance/SubstanceProperties.h.

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

// Infinite ambient reservoir at AMBIENT_TEMPERATURE_K. dQ > 0 leaves the cell.
inline float exchangeThermalEnergyWithAmbient(float &energy, float capacity,
    float kCell, float kAmbient, float dt, float areaOverDx, float conductivityScale)
{
    if (capacity < MIN_THERMAL_CAPACITY || dt <= 0.0f || areaOverDx <= 0.0f) return 0.0f;
    float k = harmonicConductivity(kCell, kAmbient) * std::max(0.0f, conductivityScale);
    if (k <= 0.0f) return 0.0f;
    float T = tempFromEnergy(energy, capacity);
    float dT = T - AMBIENT_TEMPERATURE_K;
    if (std::abs(dT) < 1.0e-6f) return 0.0f;
    float dQ = k * areaOverDx * dT * dt;
    float maxQ = dT * capacity;
    if (dQ > 0.0f) dQ = std::min(dQ, std::max(0.0f, maxQ));
    else dQ = std::max(dQ, std::min(0.0f, maxQ));
    energy -= dQ;
    if (!std::isfinite(energy) || energy < 0.0f) energy = 0.0f;
    if (capacity > MIN_THERMAL_CAPACITY) {
        float t = tempFromEnergy(energy, capacity);
        energy = energyFromTemp(capacity, t);
    }
    return dQ;
}
