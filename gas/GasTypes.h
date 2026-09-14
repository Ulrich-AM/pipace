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

// Future mixtures occupy the same cell and are transported together.
// Do not add per-species engines; extend amount into species amounts later.
enum class GasSpecies : uint8_t { Air = 0 }; // engine identity: SUBSTANCE_AIR
constexpr int GAS_SPECIES_COUNT = 1;

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
