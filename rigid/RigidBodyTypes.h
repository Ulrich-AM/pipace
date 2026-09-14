#pragma once

#include <cstdint>
#include <vector>

using MaterialId = uint16_t;
constexpr MaterialId MATERIAL_EMPTY = 0;
constexpr MaterialId MATERIAL_WOOD = 1;
constexpr MaterialId MATERIAL_STONE = 2;
constexpr MaterialId MATERIAL_GLASS = 3;
constexpr MaterialId MATERIAL_METAL = 4;
constexpr MaterialId MATERIAL_COUNT = 5;

enum class StructuralDamageType : uint8_t {
    Impact = 0,
    Shear,
    Compression,
    Tension,
    Chemical,
    Thermal,
    Pressure,
    Generic
};

struct StructuralDamageRequest {
    float wx = 0.0f;
    float wy = 0.0f;
    float amount = 0.0f; // severity in the same units as internal impact (modest values ~0.2-1.0)
    float radius = 2.0f;
    float dirX = 0.0f;
    float dirY = 0.0f;
    StructuralDamageType type = StructuralDamageType::Generic;
    bool canPropagate = true;
};

struct StructuralBond {
    float damage = 0.0f;
    bool broken = false;
};

struct MaterialDefinition {
    char const *name;
    float density;   // relative to water density 1.0
    int colorR, colorG, colorB;
    float friction;
    float restitution;
    float hardness;
    float toughness;
    float brittleness;
    float tensileStrength;
    float compressiveStrength;
    float shearStrength;
    float fractureToughness;
    float porosity;
    float permeability;
    float moistureCapacity;
    float absorptionRate;
    float dryingRate;
    float wetStrengthMultiplier;
    float wetFrictionMultiplier;
    float wetFractureToughnessMultiplier;
    float wetDensityContribution;
};

inline MaterialDefinition const &materialDef(MaterialId id) {
    static MaterialDefinition const defs[MATERIAL_COUNT] = {
        {"empty", 0.00f, 0, 0, 0,     0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f, 1.00f, 1.00f, 0.00f},
        {"wood",  0.45f, 158, 112, 62, 0.48f, 0.08f, 0.42f, 0.52f, 0.34f, 0.55f, 0.70f, 0.48f, 0.50f, 0.55f, 0.68f, 0.62f, 0.26f, 0.0010f, 0.62f, 1.30f, 0.70f, 0.00f},
        {"stone", 2.20f, 118, 122, 128, 0.55f, 0.06f, 1.35f, 1.05f, 0.28f, 1.35f, 1.80f, 1.10f, 0.88f, 0.08f, 0.08f, 0.12f, 0.040f, 0.0006f, 0.92f, 1.05f, 0.95f, 0.00f},
        {"glass", 2.50f, 168, 204, 214, 0.22f, 0.05f, 1.75f, 0.16f, 0.92f, 0.32f, 1.40f, 0.28f, 0.12f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f, 1.00f, 1.00f, 0.00f},
        {"metal", 7.80f, 148, 152, 158, 0.38f, 0.12f, 2.10f, 1.85f, 0.12f, 2.40f, 2.80f, 1.90f, 1.65f, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f, 1.00f, 1.00f, 0.00f},
    };
    if (id >= MATERIAL_COUNT) return defs[0];
    return defs[id];
}

inline bool materialIsAbsorbent(MaterialId id) {
    return materialDef(id).moistureCapacity > 1.0e-5f && materialDef(id).porosity > 1.0e-5f;
}

struct PixelRun {
    int localY = 0;
    int x0 = 0;
    int x1 = 0;
};

struct RigidContact {
    int bodyA = -1;
    int bodyB = -1; // -1 = static world
    float x = 0.0f, y = 0.0f;
    float nx = 0.0f, ny = 0.0f;
    float rAx = 0.0f, rAy = 0.0f;
    float rBx = 0.0f, rBy = 0.0f;
    float penetration = 0.0f;
    float jn = 0.0f; // accumulated normal impulse (debug)
    float jt = 0.0f; // accumulated friction impulse (debug)
    float impactSpeed = 0.0f; // approaching relative speed recorded on first hit
};

struct OccupancyConflict {
    int bodyA = 0;
    int bodyB = 0;
    int x = 0;
    int y = 0;
};

struct GrabState {
    bool active = false;
    bool strong = false;
    uint32_t bodyId = 0;
    float localX = 0.0f, localY = 0.0f;
    float targetX = 0.0f, targetY = 0.0f;
    float worldX = 0.0f, worldY = 0.0f;
    float lastFx = 0.0f, lastFy = 0.0f;
};

struct RigidBody {
    float x = 0.0f, y = 0.0f;
    float theta = 0.0f;
    float vx = 0.0f, vy = 0.0f, omega = 0.0f;
    float fx = 0.0f, fy = 0.0f, torque = 0.0f;
    float mass = 1.0f, invMass = 1.0f;
    float inertia = 1.0f, invInertia = 1.0f;
    float comLocalX = 0.0f, comLocalY = 0.0f;
    int maskW = 0, maskH = 0;
    std::vector<MaterialId> mask;
    std::vector<float> materialDamage; // 0 intact, 1 degraded substance (pixels remain)
    std::vector<float> moisture;       // absorbed liquid mass per pixel
    std::vector<float> heat;           // Joules per local source pixel
    std::vector<StructuralBond> bondsRight; // (x,y) -> (x+1,y)
    std::vector<StructuralBond> bondsDown;  // (x,y) -> (x,y+1)
    std::vector<PixelRun> runs;
    std::vector<int> occupiedLocal;
    float aabbX0 = 0.0f, aabbY0 = 0.0f, aabbX1 = 0.0f, aabbY1 = 0.0f;
    float prevAabbX0 = 0.0f, prevAabbY0 = 0.0f, prevAabbX1 = 0.0f, prevAabbY1 = 0.0f;
    int quietTicks = 0;
    bool sleeping = false;
    bool supported = false;
    bool anchored = false;
    bool dormant = false; // UI "Sleeping": pinned until disturbed
    bool structureDirty = false;
    bool moistureActive = false;
    float maxPenetration = 0.0f;
    float cachedFriction = 0.42f;
    float cachedRestitution = 0.12f;
    float cachedWetness = 0.0f;
    float absorbedLiquid = 0.0f;
    float debugJn = 0.0f;
    float debugJt = 0.0f;
    float debugPosCorrX = 0.0f;
    float debugPosCorrY = 0.0f;
    float maxDamage = 0.0f;
    float maxBondDamage = 0.0f;
    int brokenBondCount = 0;
    uint32_t id = 0;

    bool immobile() const { return anchored || dormant; }
};
