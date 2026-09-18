#pragma once

#include "GasConfig.h"
#include "fluid/BrushGeom.h"

#include <chrono>
#include <cstdint>
#include <vector>

struct FluidEngine;
struct RigidBodyEngine;

struct GasEngine {
    using Clock = std::chrono::steady_clock;

    GasConfig config;

    std::vector<float> amount;     // total cell-atmospheres; must match sum of gas components
    std::vector<float> heat;       // Joules associated with amount
    std::vector<float> heatNext;
    std::vector<float> volume;     // available gas volume in the cell, 0..1
    std::vector<float> pressure;   // cached atmospheres from amount, volume, T
    std::vector<float> u;          // MAC horizontal, cells/s
    std::vector<float> v;          // MAC vertical, cells/s
    std::vector<float> fluxU;
    std::vector<float> fluxV;
    std::vector<float> outgoing;
    std::vector<uint8_t> chunkActivity;
    std::vector<uint8_t> chunkQuietTicks;
    std::vector<uint8_t> chunkSolveMask;
    std::vector<uint8_t> thermalChunkWake;

    // Fixed-capacity SoA composition. amount[] is occupancy; these slots say what it is.
    std::vector<SubstanceId> gasCompId; // cell * kMaxGasComponents + slot
    std::vector<float> gasCompAmt;
    std::vector<uint8_t> gasCompCount;
    std::vector<SubstanceId> nextGasCompId;
    std::vector<float> nextGasCompAmt;
    std::vector<uint8_t> nextGasCompCount;

    double currentAmount = 0.0;
    double expectedAmount = 0.0;
    double amountError = 0.0;
    double escapedAmount = 0.0; // net leaving through OpenAmbient edges (entered is negative)
    double escapedHeat = 0.0;   // Joules carried by escapedAmount; inflow is negative
    double currentWaterVapor = 0.0;
    double expectedWaterVapor = 0.0;
    double escapedWaterVapor = 0.0;
    int activeGasCells = 0;
    int activeChunks = 0;
    int lastSubsteps = 1;
    bool lastWasSleeping = false;
    double lastStepMs = 0.0;
    GasTimings timing{};

    GasEngine();

    static int ci(int x, int y) { return y * GW + x; }
    static int ui(int x, int y) { return y * (GW + 1) + x; }
    static int vi(int x, int y) { return y * GW + x; }

    float pressurePa(int index) const;
    float pressureAtm(int index) const { return index >= 0 && index < GW * GH ? pressure[static_cast<size_t>(index)] : 0.0f; }
    float cellTemperatureK(int index) const;
    void recomputePressure();
    float cellU(int x, int y) const;
    float cellV(int x, int y) const;
    bool isAccessible(FluidEngine const &fluid, int x, int y) const;
    float availableVolume(FluidEngine const &fluid, int x, int y) const;

    void resetAmbient(FluidEngine &fluid);
    void handleWorldEdit(FluidEngine &fluid);
    void applyPressureBrush(FluidEngine &fluid, int cx, int cy, int brushRadius,
        float signedAtmPerSec, float dt, BrushShape shape = BrushShape::Circle);
    void applyGasBrush(FluidEngine &fluid, int cx, int cy, int brushRadius,
        SubstanceId gas, float amountPerSec, float dt, BrushShape shape = BrushShape::Circle);
    void eraseAmountBrush(FluidEngine &fluid, int cx, int cy, int brushRadius,
        BrushShape shape = BrushShape::Circle);
    void simulationTick(FluidEngine &fluid);
    void applyPressureForces(RigidBodyEngine &rigid, FluidEngine const &fluid) const;
    void loadTestScene(FluidEngine &fluid, RigidBodyEngine &rigid, int scene);
    void runDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid);
    void runCompositionSanityCheck(FluidEngine &fluid, RigidBodyEngine &rigid);

    float gasComponentAmount(int index, SubstanceId id) const;
    float gasComponentFraction(int index, SubstanceId id) const;
    SubstanceId dominantGasSubstance(int index) const;
    bool gasCompositionValid(int index) const;
    GasComponentView gasComponents(int index) const;
    void setGasComponentAmount(int index, SubstanceId id, float componentAmount);
    void addGasComponentAmount(int index, SubstanceId id, float delta);
    float takeGasComponentAmount(int index, SubstanceId id, float da);
    bool tryCommitGasOccupancy(int index, GasComponentView const &view);
    void clearGasComposition(int index);
    float gasPartialPressurePa(int index, SubstanceId id) const;

    // Compatibility wrappers over generic composition. Water vapor is
    // SUBSTANCE_WATER + MatterPhase::Gas; Air is SUBSTANCE_AIR. Not a second store.
    float vaporAmount(int index) const;
    float airAmount(int index) const;
    float vaporFraction(int index) const;
    void clampSpecies(int index);
    void addWaterVapor(int index, float da);
    float takeWaterVapor(int index, float da);
    void wakeAt(int x, int y);
    double sumWaterVapor() const;

private:
    std::vector<uint32_t> relocateStamp;
    std::vector<int> relocateQueue;
    uint32_t relocateEpoch = 1;
    std::vector<float> prevVolume;
    std::vector<int> lastOccupancyScan;
    int volumeScanTick = 0;

    static int compositionSlot(int cell, int slot);
    void compactGasComposition(int index);
    void copyGasCompToNext(int index);
    void commitGasCompFromNext(int index);
    float addGasComponentUntracked(int index, SubstanceId id, float componentAmount);
    float addNextGasComponentUntracked(int index, SubstanceId id, float componentAmount);
    void writePureGas(int index, SubstanceId id, float componentAmount);
    void scaleGasComposition(int index, float frac);
    void syncAmountFromComposition(int index);
    void rebuildVolumes(FluidEngine const &fluid, bool &volumeChanged, bool fullGrid = true);
    float relocateAmount(FluidEngine const &fluid, int x, int y, GasComponentView parcel,
        float leftoverHeatPerAmount = 0.0f, float maxAtm = 0.0f);
    void wakeThermalAt(int x, int y);
    void displaceBlocked(FluidEngine const &fluid);
    void displaceLiquidOverflow(FluidEngine const &fluid);
    void applyMovingBoundary(FluidEngine const &fluid);
    void accumulateTransfer(FluidEngine const &fluid, float dt);
    void applyBoundaryFlux(float dt);
    void applyFluxes();
    void integrateVelocity(FluidEngine const &fluid, float dt);
    void wakeChunkAtCell(int x, int y, bool resetQuiet = true);
    void wakeAll();
    void sleepAll();
    void rebuildActivity(bool advanceSleep);
    void wallRect(FluidEngine &fluid, int x0, int y0, int x1, int y1);
    void fillRectAmount(FluidEngine const &fluid, int x0, int y0, int x1, int y1, float atm);
    void vacuumAll();
    void commitExpected();
    double sumAmount() const;
    void transferSpecies(int donor, int receiver, float q);
};
