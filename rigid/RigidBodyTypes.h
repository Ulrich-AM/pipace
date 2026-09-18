#pragma once

#include "substance/SubstanceTypes.h"

#include <cstdint>
#include <vector>

using MaterialId = uint16_t;
constexpr MaterialId MATERIAL_EMPTY = 0;
constexpr MaterialId MATERIAL_WOOD = 1;
constexpr MaterialId MATERIAL_STONE = 2;
constexpr MaterialId MATERIAL_GLASS = 3;
constexpr MaterialId MATERIAL_METAL = 4;
constexpr MaterialId MATERIAL_WATER_SOLID = 5; // compatibility mask for SUBSTANCE_WATER + Solid
constexpr MaterialId MATERIAL_COUNT = 6;
// Rigid source masks still store MaterialId. Map to engine identity with
// substanceForMaterial(). MATERIAL_WATER_SOLID is not a SubstanceId.

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

// Compatibility adapter for code that still keys off MaterialId (masks,
// fracture, moisture, mass, RGB). Values are copied from SubstanceDefinition
// at first use; do not add independent physical constants here.
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

inline MaterialDefinition materialDefinitionFromSubstance(SubstanceDefinition const &s) {
    MaterialDefinition m{};
    m.name = s.internalName;
    m.density = s.mechanical.densityRel;
    m.colorR = s.visual.colorR;
    m.colorG = s.visual.colorG;
    m.colorB = s.visual.colorB;
    m.friction = s.mechanical.friction;
    m.restitution = s.mechanical.restitution;
    m.hardness = s.mechanical.hardness;
    m.toughness = s.mechanical.toughness;
    m.brittleness = s.mechanical.brittleness;
    m.tensileStrength = s.mechanical.tensileStrength;
    m.compressiveStrength = s.mechanical.compressiveStrength;
    m.shearStrength = s.mechanical.shearStrength;
    m.fractureToughness = s.mechanical.fractureToughness;
    m.porosity = s.porous.porosity;
    m.permeability = s.porous.permeability;
    m.moistureCapacity = s.porous.moistureCapacity;
    m.absorptionRate = s.porous.absorptionRate;
    m.dryingRate = s.porous.dryingRate;
    m.wetStrengthMultiplier = s.porous.wetStrengthMultiplier;
    m.wetFrictionMultiplier = s.porous.wetFrictionMultiplier;
    m.wetFractureToughnessMultiplier = s.porous.wetFractureToughnessMultiplier;
    m.wetDensityContribution = s.porous.wetDensityContribution;
    return m;
}

inline MaterialDefinition const &materialDef(MaterialId id) {
    static MaterialDefinition const defs[MATERIAL_COUNT] = {
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(0))),
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(1))),
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(2))),
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(3))),
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(4))),
        materialDefinitionFromSubstance(substanceDef(substanceForMaterialId(5))),
    };
    if (id >= MATERIAL_COUNT) return defs[0];
    return defs[id];
}

inline bool materialIsAbsorbent(MaterialId id) {
    PorousProperties const &p = porousForSubstance(substanceForMaterialId(id));
    return p.valid && p.moistureCapacity > 1.0e-5f && p.porosity > 1.0e-5f;
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

struct GrabMember {
    uint32_t bodyId = 0;
    float localX = 0.0f, localY = 0.0f;
    float relX = 0.0f, relY = 0.0f;
};

struct GrabState {
    bool active = false;
    bool strong = false;
    bool phantom = false;
    float strength = 1.0f;
    uint32_t bodyId = 0;
    float localX = 0.0f, localY = 0.0f;
    float targetX = 0.0f, targetY = 0.0f;
    float worldX = 0.0f, worldY = 0.0f;
    float lastFx = 0.0f, lastFy = 0.0f;
    std::vector<GrabMember> members;
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
    std::vector<float> solidRemain;    // 0..1 remaining fraction of the solid source pixel
                                       // (ice melt and chemistry share this; no second array)
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
    bool massDirty = false; // remain changed without mask topology change
    bool moistureActive = false;
    float maxPenetration = 0.0f;
    float cachedFriction = 0.42f;
    float cachedRestitution = 0.12f;
    float cachedWetness = 0.0f;
    float absorbedLiquid = 0.0f;
    float pendingDrip = 0.0f; // conserved exudation waiting for a visible quantum
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
