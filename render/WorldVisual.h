#pragma once

#include "rigid/RigidBodyTypes.h"
#include "substance/SubstanceTypes.h"

#include <cstdint>

enum class WorldRenderStyle : uint8_t {
    FlatColor = 0,
    NoisyFlat,
    Detailed,
    Realistic,
    AlphaFlat,
    Legacy
};

struct WorldLook {
    WorldRenderStyle style = WorldRenderStyle::Detailed;
    bool glowingLiquids = false;
    bool outlines = false;
    float noiseAmount = 1.0f; // 0..2, 1 = 100%
};

// Cosmetic only. Physics lives on SubstanceDefinition grouped properties.
// Physics must not include this header (no MaterialVisual in engines).
struct MaterialVisual {
    int r = 128;
    int g = 128;
    int b = 128;
    float noiseStrength = 0.12f;
    float depthStrength = 0.45f;
    float realisticDetail = 0.55f;
    bool glowEligible = false;
    float glowStrength = 0.0f;
    float outlineDarken = 0.32f;
    float alphaResponse = 1.0f; // Alpha Flat: 1 = fill/amount drives opacity
};

inline MaterialVisual const kVisualVoid{10, 17, 28, 0.0f, 0.0f, 0.0f, false, 0.0f, 0.0f, 0.0f};
inline MaterialVisual const kVisualWall{96, 100, 108, 0.18f, 0.28f, 0.40f, false, 0.0f, 0.28f, 0.0f};
inline MaterialVisual const kVisualWater{22, 126, 214, 0.10f, 0.55f, 0.72f, true, 0.85f, 0.30f, 1.0f};
inline MaterialVisual const kVisualHoney{176, 110, 22, 0.09f, 0.50f, 0.65f, true, 0.80f, 0.28f, 1.0f};

// Compatibility (Prompt 1): rigid MaterialId visual. Prefer visualForSubstance
// when the caller already has a SubstanceId.
inline MaterialVisual visualForSolid(MaterialId id) {
    SubstanceVisualMetadata const &vis = substanceDef(substanceForMaterialId(id)).visual;
    MaterialVisual v;
    v.r = vis.colorR;
    v.g = vis.colorG;
    v.b = vis.colorB;
    v.glowEligible = false;
    v.glowStrength = 0.0f;
    v.alphaResponse = 0.0f;
    v.outlineDarken = 0.30f;
    v.noiseStrength = 0.14f;
    v.depthStrength = 0.42f;
    v.realisticDetail = 0.50f;
    if (id == MATERIAL_WOOD) {
        v.noiseStrength = 0.16f;
        v.depthStrength = 0.40f;
        v.realisticDetail = 0.52f;
        v.outlineDarken = 0.34f;
    } else if (id == MATERIAL_STONE) {
        v.noiseStrength = 0.13f;
        v.depthStrength = 0.44f;
        v.realisticDetail = 0.48f;
    } else if (id == MATERIAL_GLASS) {
        v.noiseStrength = 0.05f;
        v.depthStrength = 0.32f;
        v.realisticDetail = 0.58f;
        v.alphaResponse = 0.55f;
        v.outlineDarken = 0.22f;
    } else if (id == MATERIAL_METAL) {
        v.noiseStrength = 0.08f;
        v.depthStrength = 0.46f;
        v.realisticDetail = 0.70f;
        v.outlineDarken = 0.26f;
    }
    return v;
}

inline MaterialVisual visualForSubstance(SubstanceId id) {
    switch (id) {
        case SUBSTANCE_WATER: return kVisualWater;
        case SUBSTANCE_HONEY: return kVisualHoney;
        case SUBSTANCE_AIR: return kVisualVoid;
        case SUBSTANCE_WOOD: return visualForSolid(MATERIAL_WOOD);
        case SUBSTANCE_STONE: return visualForSolid(MATERIAL_STONE);
        case SUBSTANCE_GLASS: return visualForSolid(MATERIAL_GLASS);
        case SUBSTANCE_METAL: return visualForSolid(MATERIAL_METAL);
        default: return visualForSolid(MATERIAL_EMPTY);
    }
}
