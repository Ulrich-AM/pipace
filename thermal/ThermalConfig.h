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
    float gasConductivityScale = 40.0f; // air is a poor conductor; do not use the solid sandbox scale
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
