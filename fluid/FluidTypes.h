#pragma once

#include "thermal/ThermalTypes.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#ifndef PIPACE_GRID_WIDTH
#define PIPACE_GRID_WIDTH 200
#endif
#ifndef PIPACE_GRID_HEIGHT
#define PIPACE_GRID_HEIGHT 120
#endif

constexpr int GW = PIPACE_GRID_WIDTH;
constexpr int GH = PIPACE_GRID_HEIGHT;
constexpr int CHUNK = 16;
constexpr int CHUNK_W = (GW + CHUNK - 1) / CHUNK;
constexpr int CHUNK_H = (GH + CHUNK - 1) / CHUNK;
constexpr float MIN_RENDER_FILL = 0.045f;
constexpr float MIN_ACTIVE_FILL = 0.010f;
constexpr float MIN_PRESSURE_FILL = 0.030f;
constexpr float MIN_CONSOLIDATE_FILL = 0.55f;
constexpr float MIN_SUBSTANTIAL_FILL = 0.50f;
constexpr float PHYSICS_DT = 1.0f / 30.0f;

enum class Tool : uint8_t {
    Solid = 0, Water = 1, Eraser = 2, Rigid = 3, Grab = 4,
    Heat = 5, Cool = 6, Pressurize = 7, Depressurize = 8
};

// Brush extras for Tool::Water. Dye never creates fill; honey is composition inside fill.
// asHoney remains the paint-path flag. Identity is paint.substance().
struct LiquidPaint {
    bool asHoney = false;
    bool dyeOnly = false;
    float dyeR = 0.0f, dyeG = 0.0f, dyeB = 0.0f; // 0–1 stain color; all zero = clear
    float dyeStrength = 0.60f;
    SubstanceId substance() const { return substanceForLiquidPaint(asHoney); }
};

// Conserved extras that ride with a volume parcel (heat, dye mass, honey volume).
struct LiquidCarry {
    float heat = 0.0f;
    float dyeR = 0.0f, dyeG = 0.0f, dyeB = 0.0f;
    float honey = 0.0f;
};

inline bool isEnergyTool(Tool t) {
    return t == Tool::Heat || t == Tool::Cool || t == Tool::Pressurize || t == Tool::Depressurize;
}
inline bool isThermalEnergyTool(Tool t) {
    return t == Tool::Heat || t == Tool::Cool;
}
inline bool isPressureEnergyTool(Tool t) {
    return t == Tool::Pressurize || t == Tool::Depressurize;
}
enum class DebugView : uint8_t {
    Normal, Fill, Pressure, Velocity, Divergence, Chunks, Rigid,
    GasPressure, GasAmount, GasVelocity, Temperature, Moisture
};
enum class VelocityAdvection : uint8_t {
    None = 0,
    FirstOrderUpwind,
    NearestSemiLagrangian,
    SemiLagrangian,
    MacCormack,
    BFECC
};

enum class QualityPreset : uint8_t { Low = 0, Medium = 1, High = 2, Auto = 3 };

inline VelocityAdvection nextVelocityAdvection(VelocityAdvection mode) {
    switch (mode) {
        case VelocityAdvection::None: return VelocityAdvection::FirstOrderUpwind;
        case VelocityAdvection::FirstOrderUpwind: return VelocityAdvection::NearestSemiLagrangian;
        case VelocityAdvection::NearestSemiLagrangian: return VelocityAdvection::SemiLagrangian;
        case VelocityAdvection::SemiLagrangian: return VelocityAdvection::MacCormack;
        case VelocityAdvection::MacCormack: return VelocityAdvection::BFECC;
        case VelocityAdvection::BFECC: return VelocityAdvection::None;
    }
    return VelocityAdvection::SemiLagrangian;
}

inline char const *velocityAdvectionKey(VelocityAdvection mode) {
    switch (mode) {
        case VelocityAdvection::None: return "btn_adv_none";
        case VelocityAdvection::FirstOrderUpwind: return "btn_adv_fou";
        case VelocityAdvection::NearestSemiLagrangian: return "btn_adv_nsl";
        case VelocityAdvection::SemiLagrangian: return "btn_adv_sl";
        case VelocityAdvection::MacCormack: return "btn_adv_macc";
        case VelocityAdvection::BFECC: return "btn_adv_bfecc";
    }
    return "btn_adv_sl";
}

inline char const *velocityAdvectionName(VelocityAdvection mode) {
    switch (mode) {
        case VelocityAdvection::None: return "none";
        case VelocityAdvection::FirstOrderUpwind: return "fou";
        case VelocityAdvection::NearestSemiLagrangian: return "nsl";
        case VelocityAdvection::SemiLagrangian: return "sl";
        case VelocityAdvection::MacCormack: return "maccormack";
        case VelocityAdvection::BFECC: return "bfecc";
    }
    return "sl";
}

// Compatibility name. Authoritative liquid tables: fluidForSubstance(id).
using LiquidProperties = FluidProperties;

struct SplashParticle {
    float x = 0.0f, y = 0.0f;
    float vx = 0.0f, vy = 0.0f;
    float volume = 0.0f;
    float life = 0.0f;
    int originX = -1, originY = -1;
    float heat = 0.0f; // Joules carried with this droplet
    uint8_t bounceCount = 0;
    float dyeR = 0.0f, dyeG = 0.0f, dyeB = 0.0f;
    float honey = 0.0f;
};

struct TimingAverages {
    double advection = 0.0;
    double viscosity = 0.0;
    double forces = 0.0;
    double surface = 0.0;
    double pressure = 0.0;
    double transport = 0.0;
    double residual = 0.0;
    double drain = 0.0;
    double splashes = 0.0;
    double chunks = 0.0;
    double listsBuild = 0.0;
    double physics = 0.0;
    double render = 0.0;
    double thermal = 0.0;
};

struct WorkCounts {
    int pressureCells = 0;
    int surfaceCells = 0;
    int nonzeroFluxU = 0;
    int nonzeroFluxV = 0;
    int limiterPasses = 0;
    int splashCount = 0;
    int subvisibleCells = 0;
    int residualTransfers = 0;
    int activeUFaces = 0;
    int activeVFaces = 0;
    int thermalCells = 0;
    int thermalChunks = 0;
    int thermalPairs = 0;
};
