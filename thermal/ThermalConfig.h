#pragma once

#include "ThermalTypes.h"

struct ThermalConfig {
    bool enabled = true;
    // Conduction/sleep runs every N physics ticks. Liquid/gas still carry heat every tick.
    int intervalTicks = 2;
    float sleepTempEps = 0.08f;       // K, neighbor |ΔT| to stay awake
    float sleepAmbientEps = 0.25f;    // K from ambient to stay awake
    int sleepQuietTicks = 12;
    float conductivityScale = 25000.0f;
    // Cell edge is 0.25 m, so physical air diffusion across a cell is thousands of
    // seconds. This sandbox scale (~400) makes gas-gas / gas-surface conduction
    // visible in a few seconds while staying far below the solid scale (25000).
    // Bulk hot-air mixing is buoyancy (GasConfig::buoyancyScale), not this number.
    // Unwalled map edges use this scale plus air conductivity against an infinite
    // reservoir at AMBIENT_TEMPERATURE_K (FluidConfig::walledBorders == false).
    float gasConductivityScale = 400.0f;
    float heatToolWatts = 2.5e6f;     // sandbox watts at strength 1, per covered cell
    bool coupleGasPressureToTemperature = false; // documented path; keep off
};

inline void applyThermalQualityKnobs(ThermalConfig &c, int level) {
    if (level <= 0) {
        c.intervalTicks = 3;
        c.sleepTempEps = 0.20f;
        c.sleepAmbientEps = 0.50f;
        c.sleepQuietTicks = 8;
    } else if (level == 1) {
        c.intervalTicks = 2;
        c.sleepTempEps = 0.08f;
        c.sleepAmbientEps = 0.25f;
        c.sleepQuietTicks = 12;
    } else {
        c.intervalTicks = 1;
        c.sleepTempEps = 0.05f;
        c.sleepAmbientEps = 0.15f;
        c.sleepQuietTicks = 16;
    }
}
