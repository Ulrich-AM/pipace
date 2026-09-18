#pragma once

#include "render/WorldVisual.h"

#include <cstdint>
#include <vector>

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;

struct WorldRenderer {
    void paintNormal(FluidEngine &fluid, RigidBodyEngine const &rigid, WorldLook const &look,
        float const *reactionActivity = nullptr, GasEngine const *gas = nullptr);
    void ensureSize();

    std::vector<uint8_t> visualLiquid;
    std::vector<int> liquidDepth;
    std::vector<int> rigidDepth;
    std::vector<int> workQueue;
    std::vector<float> glowDist;
    std::vector<uint8_t> glowR;
    std::vector<uint8_t> glowG;
    std::vector<uint8_t> glowB;
};

void runWorldLookBenchmark(FluidEngine &fluid, RigidBodyEngine &rigid, WorldRenderer &renderer);
