#pragma once

#include "ThermalTypes.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

// TEMP debug visualization only. Does not clamp or rewrite simulation temperatures.
// Ambient-centered, piecewise |ΔT| scale with a display dead band so 0.01 K of
// floating-point noise is not a different color from 293.15 K.

enum class ThermalSampleKind : uint8_t {
    Empty = 0,
    Gas,
    Liquid,
    Wall,
    Rigid
};

struct ThermalCellSample {
    float temperatureK = AMBIENT_TEMPERATURE_K;
    float energyJ = 0.0f;
    float capacityJK = 0.0f;
    ThermalSampleKind kind = ThermalSampleKind::Empty;
    uint16_t materialId = 0;
    bool hasMatter = false;
};

constexpr float TEMP_VIZ_DEADBAND_K = 0.35f;

inline float tempVizSignedUnit(float temperatureK) {
    if (!std::isfinite(temperatureK)) return 0.0f;
    float dT = temperatureK - AMBIENT_TEMPERATURE_K;
    float mag = std::abs(dT) - TEMP_VIZ_DEADBAND_K;
    if (mag <= 0.0f) return 0.0f;
    float u = 0.0f;
    if (mag < 20.0f) u = 0.40f * (mag / 20.0f);                         // ~273 K / 313 K
    else if (mag < 80.0f) u = 0.40f + 0.25f * ((mag - 20.0f) / 60.0f);  // ~373 K
    else if (mag < 280.0f) u = 0.65f + 0.20f * ((mag - 80.0f) / 200.0f); // ~573 K
    else {
        float tail = 1.0f - std::exp(-(mag - 280.0f) / 900.0f);
        u = 0.85f + 0.15f * tail;                                       // 1000 K+
    }
    u = std::clamp(u, 0.0f, 1.0f);
    return (dT < 0.0f) ? -u : u;
}

inline void tempVizLerp(int r0, int g0, int b0, int r1, int g1, int b1, float t,
    int &r, int &g, int &b)
{
    t = std::clamp(t, 0.0f, 1.0f);
    r = static_cast<int>(std::lround(r0 + (r1 - r0) * t));
    g = static_cast<int>(std::lround(g0 + (g1 - g0) * t));
    b = static_cast<int>(std::lround(b0 + (b1 - b0) * t));
}

// Continuous through ambient (same RGB from both sides at unit=0).
// Cool: slate → cyan → blue → deep blue. Warm: slate → yellow → orange → red → pale.
inline void tempVizRgb(float signedUnit, int &r, int &g, int &b) {
    constexpr int nR = 46, nG = 54, nB = 60; // muted slate, not cyan and not green
    float u = std::clamp(signedUnit, -1.0f, 1.0f);
    if (u <= 0.0f) {
        float q = -u;
        if (q < 0.45f) tempVizLerp(nR, nG, nB, 36, 118, 168, q / 0.45f, r, g, b);
        else if (q < 0.80f) tempVizLerp(36, 118, 168, 22, 48, 150, (q - 0.45f) / 0.35f, r, g, b);
        else tempVizLerp(22, 48, 150, 10, 14, 72, (q - 0.80f) / 0.20f, r, g, b);
    } else {
        float q = u;
        if (q < 0.35f) tempVizLerp(nR, nG, nB, 196, 168, 62, q / 0.35f, r, g, b);
        else if (q < 0.65f) tempVizLerp(196, 168, 62, 214, 108, 36, (q - 0.35f) / 0.30f, r, g, b);
        else if (q < 0.88f) tempVizLerp(214, 108, 36, 196, 42, 28, (q - 0.65f) / 0.23f, r, g, b);
        else tempVizLerp(196, 42, 28, 255, 244, 214, (q - 0.88f) / 0.12f, r, g, b);
    }
}

inline void tempVizOverlay(float temperatureK, ThermalSampleKind kind, int &r, int &g, int &b, float &alpha) {
    float unit = tempVizSignedUnit(temperatureK);
    tempVizRgb(unit, r, g, b);
    float mag = std::abs(unit);
    if (kind == ThermalSampleKind::Empty) {
        alpha = 0.0f;
        return;
    }
    if (mag <= 1.0e-6f) {
        // Ambient matter: gas almost invisible; solids/liquid still faintly readable.
        alpha = (kind == ThermalSampleKind::Gas) ? 0.055f : 0.20f;
        return;
    }
    float visible = 0.42f + 0.58f * mag;
    if (kind == ThermalSampleKind::Gas) visible = 0.28f + 0.72f * mag;
    alpha = std::clamp(visible, 0.0f, 1.0f);
}
