#pragma once

#include "ThermalConfig.h"
#include "fluid/BrushGeom.h"
#include "ThermalViz.h"

#include <cstdint>
#include <vector>

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct RigidBody;

struct ThermalWorkCounts {
    int activeCells = 0;
    int activeChunks = 0;
    int conductionPairs = 0;
    int activeBodies = 0;
};

struct ThermalWorldStats {
    double energy = 0.0;
    float tMin = 0.0f;
    float tMax = 0.0f;
    float tMean = 0.0f;
    float tMinGas = 0.0f, tMaxGas = 0.0f;
    float tMinLiquid = 0.0f, tMaxLiquid = 0.0f;
    float tMinSolid = 0.0f, tMaxSolid = 0.0f;
    float tMinRigid = 0.0f, tMaxRigid = 0.0f;
    int nSampled = 0;
    int nGas = 0, nLiquid = 0, nWall = 0, nRigid = 0, nEmpty = 0;
    int nNan = 0, nInf = 0, nNegK = 0, nZeroK = 0;
    int leftoverLiquidCells = 0, leftoverSolidCells = 0;
    double leftoverLiquidEnergy = 0.0, leftoverSolidEnergy = 0.0;
    int nAwayFromAmbient = 0; // |ΔT| > TEMP_VIZ_DEADBAND_K
    float maxAbsDeltaAmbient = 0.0f;
    float tMinLiquidFill = 0.0f, tMaxLiquidFill = 0.0f;
};

struct ThermalEngine {
    ThermalConfig config;
    std::vector<uint8_t> chunkActivity;
    std::vector<uint8_t> chunkQuietTicks;

    ThermalWorkCounts work{};
    double lastStepMs = 0.0;
    double timingAccum = 0.0;
    double timingAverage = 0.0;
    int timingTicks = 0;
    uint32_t tickNo = 0;

    ThermalEngine();

    void clear();
    void seedAmbient(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas);
    void wakeCell(int x, int y);
    void wakeRect(int x0, int y0, int x1, int y1);
    void ingestEngineWakes(FluidEngine &fluid, GasEngine &gas);

    void applyBrush(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
        int cx, int cy, int brushRadius, float signedStrength, float dt,
        BrushShape shape = BrushShape::Circle);

    void simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas, float dt);

    static float sampleTemperatureK(FluidEngine const &fluid, RigidBodyEngine const &rigid,
        GasEngine const &gas, int x, int y);
    static ThermalCellSample sampleCell(FluidEngine const &fluid, RigidBodyEngine const &rigid,
        GasEngine const &gas, int x, int y);
    ThermalWorldStats collectStats(FluidEngine const &fluid, RigidBodyEngine const &rigid,
        GasEngine const &gas) const;
    bool isChunkActive(int x, int y) const;
    double totalThermalEnergy(FluidEngine const &fluid, RigidBodyEngine const &rigid,
        GasEngine const &gas) const;

    static float liquidCapacity(FluidEngine const &fluid, int index);
    static float wallCapacity(FluidEngine const &fluid, int index);
    static float gasCapacity(GasEngine const &gas, int index);
    static float rigidPixelCapacity(RigidBody const &b, int localIndex);
    static float liquidTempK(FluidEngine const &fluid, int index);
    static float wallTempK(FluidEngine const &fluid, int index);
    static float gasTempK(GasEngine const &gas, int index);
    static float rigidPixelTempK(RigidBody const &b, int localIndex);

    void runDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas);
    void runSpreadDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas);

private:
    void conductActive(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas, float dt);
    void conductRigidBodies(RigidBodyEngine &rigid, float dt, float areaOverDx);
    void sleepChunks(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas);
    void wakeMovingBodies(RigidBodyEngine const &rigid);
    void publishTiming(FluidEngine &fluid);
    void scrubMasslessHeat(FluidEngine &fluid);
};
