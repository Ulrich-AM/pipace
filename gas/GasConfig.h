#pragma once

#include "GasTypes.h"
#include "thermal/ThermalTypes.h"

struct GasConfig {
    float referencePressurePa = GAS_REFERENCE_PRESSURE_PA;
    // Face transfer rate: amount per atmosphere of ΔP, per unit face volume, per second.
    float conductivity = 3.2f;
    float maxFaceTransfer = 0.40f;     // amount per face per substep
    float maxExtractFraction = 0.45f;  // donor limiter
    float velocityDamping = 0.22f;
    float pressureAccel = 18.0f;       // cells/s² per atm/cell (not physical sound speed)
    float maxVelocity = 40.0f;
    float sleepSpeed = 0.08f;
    float sleepPressureDelta = 0.004f; // atmospheres
    int sleepQuietTicks = 22;
    int maxSubsteps = 4;
    float maxTravelPerSubstep = 0.45f;
    // Force on rigid faces uses normalized atmospheres so a uniform field cancels.
    float rigidForceScale = 22.0f;
    GasBoundary boundary = GasBoundary::OpenAmbient;
    float ambientPressureAtm = 1.0f;
    float brushAtmPerSec = 2.0f;   // pressurize/depressurize at Power 1
    float brushMaxAtm = 8.0f;
    GasSimMode simMode = GasSimMode::Full;
    ThermalProperties thermal = kAirThermal();
};
