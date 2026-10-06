#pragma once

#include "FluidConfig.h"
#include "WorkerPool.h"
#include "BrushGeom.h"
#include "substance/LiquidMixtureProperties.h"

#include <chrono>
#include <cstdint>
#include <vector>

struct PressureStencil {
    int index = 0;
    int pressureNeighbors[4]{};
    int pressureNeighborCount = 0;
    int freeDx[4]{};
    int freeDy[4]{};
    int freeNeighborCount = 0;
};

struct FluidEngine {
    using Clock = std::chrono::steady_clock;

    FluidConfig config;

    std::vector<uint8_t> solid;
    std::vector<float> solidHeat;
    std::vector<uint8_t> dynamicSolid;
    std::vector<float> dynamicVelX;
    std::vector<float> dynamicVelY;
    std::vector<int> dynamicOccupiedCells;
    std::vector<float> fill;
    std::vector<float> nextFill;
    std::vector<float> liquidHeat;
    std::vector<float> nextHeat;
    std::vector<float> dyeR, dyeG, dyeB; // dye mass; concentration = dye / fill
    std::vector<RuntimeSubstanceRef> liquidCompId; // SoA runtime identity: cell * kMaxLiquidComponents + slot
    std::vector<float> liquidCompAmt;
    std::vector<uint8_t> liquidCompCount;
    std::vector<SubstanceId> solidifyPendingId; // one pending solid SubstanceId per cell
    std::vector<float> solidifyPendingKg;       // sub-pixel solid mass waiting for a rigid pixel
    std::vector<float> solidifyPendingHeatJ;    // sensible energy of that pending solid mass
    std::vector<float> nextDyeR, nextDyeG, nextDyeB;
    std::vector<RuntimeSubstanceRef> nextCompId;
    std::vector<float> nextCompAmt;
    std::vector<uint8_t> nextCompCount;
    std::vector<float> pressure;
    std::vector<float> divergenceField;
    std::vector<float> outgoing;
    std::vector<float> incoming;
    std::vector<float> donorScale;
    std::vector<float> receiverScale;
    std::vector<float> cellForceX;
    std::vector<float> cellForceY;
    std::vector<float> curlField;
    std::vector<float> smoothedFill;
    std::vector<float> surfaceNormalX;
    std::vector<float> surfaceNormalY;
    std::vector<float> surfaceCurvature;
    std::vector<uint8_t> surfaceMask;
    std::vector<float> pressureBefore;
    std::vector<float> previousFill;
    std::vector<int> residualTarget;
    std::vector<float> residualTransfer;
    std::vector<int8_t> waterShade;

    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> uTemp;
    std::vector<float> vTemp;
    std::vector<float> uScratch;
    std::vector<float> vScratch;
    std::vector<float> fluxU;
    std::vector<float> fluxV;

    std::vector<SplashParticle> splashes;
    std::vector<uint8_t> chunkActivity;
    std::vector<uint8_t> chunkQuietTicks;
    std::vector<uint8_t> chunkHasFluid;
    std::vector<uint8_t> chunkSolveMask;
    std::vector<uint8_t> chunkHaloSource;
    std::vector<uint8_t> thermalChunkWake;
    std::vector<uint32_t> pixels;

    // Sparse work lists / caches
    std::vector<int> pressureRed;
    std::vector<int> pressureBlack;
    std::vector<PressureStencil> pressureStencils;
    std::vector<int> surfaceCells;
    std::vector<int> nonzeroFluxUFaces;
    std::vector<int> nonzeroFluxVFaces;
    std::vector<int> fluxTouchedCells;

    uint32_t tickNo = 1;
    double expectedVolume = 0.0;
    double currentVolume = 0.0;
    double volumeError = 0.0;
    double volumeLostRigid = 0.0;
    double volumeDisplacedRigid = 0.0;
    double escapedHeat = 0.0; // Joules carried out with liquid/splash through open rims
    float measuredMaxVelocity = 0.0f;
    int activeFluidCells = 0;
    int activeChunks = 0;
    int lastFluidSubsteps = 1;
    int lastPressureIterations = 0;
    bool substepCapReached = false;
    int thinCellCount = 0;
    double thinVolume = 0.0;
    double momentumX = 0.0, momentumY = 0.0, kineticEnergy = 0.0;
    int solveX0 = 0, solveY0 = 0, solveX1 = GW - 1, solveY1 = GH - 1;
    bool hasActiveSolveRegion = true;
    bool residualConsolidationEnabled = true;
    bool paintDirty = false;
    bool useSparsePressure = true;
    int lastResolvedWorkers = 1;
    bool lastPressureParallel = false;
    // -1 = choose from flux density; 0 = dense limiter; 1 = sparse limiter.
    int forceSparseFlux = -1;
    double lastAdvectOverflowVol = 0.0;
    double accumAdvectOverflowVol = 0.0;
    int pendingSolidCellCount = 0;

    TimingAverages timingAccum{};
    TimingAverages timingAverage{};
    WorkCounts workCounts{};
    int timingTicks = 0;
    int renderSamples = 0;

    WorkerPool workerPool;
    std::vector<uint32_t> relocateStamp;
    std::vector<int> relocateQueue;
    uint32_t relocateEpoch = 1;
    std::vector<int> advectOverflowIndex;
    std::vector<float> advectOverflowVol;
    std::vector<LiquidCarry> advectOverflowCarry;

    FluidEngine();

    static int ci(int x, int y) { return y * GW + x; }
    static int ui(int x, int y) { return y * (GW + 1) + x; }
    static int vi(int x, int y) { return y * GW + x; }
    static bool inside(int x, int y) { return x >= 0 && x < GW && y >= 0 && y < GH; }
    static float liquidFaceFraction(float a, float b);
    static uint32_t hashCell(uint32_t x, uint32_t y, uint32_t tick);
    static double elapsedMs(Clock::time_point start);

    bool isSolid(int x, int y) const;
    bool isStaticSolid(int x, int y) const;
    bool isPaintedSolid(int x, int y) const;
    bool isFluid(int x, int y) const;
    bool isPressureFluid(int x, int y) const;
    float gridGravity() const;
    bool openUFace(int x, int y) const;
    bool openVFace(int x, int y) const;
    int activeX0(int halo = 0) const;
    int activeY0(int halo = 0) const;
    int activeX1(int halo = 0) const;
    int activeY1(int halo = 0) const;
    int8_t makeShade(int x, int y) const;
    float cellU(int x, int y) const;
    float cellV(int x, int y) const;
    int splashAtCell(int x, int y) const;
    void syncWorkerPool();
    int resolvedWorkerCount() const;
    int autoWorkerCount() const;
    int maxSelectableWorkers() const;
    bool useParallelPressure(int cellCount) const;

    void enforceSolidBoundaries();
    void enforceActiveBoundaries();
    void advectVelocity(float dt);
    void diffuseVelocity(float dt);
    void applyGravity(float dt);
    bool surfaceCell(int x, int y) const;
    void updateSurfaceField();
    void applySurfaceTension(float dt);
    void applyVorticityConfinement(float dt);
    void clampVelocity();
    void buildPressureWorkLists();
    void projectVelocity(float dt);
    void advectLiquidVolume(float dt);
    void consolidateResidualVolume();
    void repairEnclosedUndersaturatedCells();
    void promoteUnsupportedIsolatedLiquid();
    float depositVolume(float x, float y, float amount, float momentumX, float momentumY, LiquidCarry *carry = nullptr);
    // fromOccupiedCell: start cell is now solid (rigid squeeze). Never use this for
    // splash/paint deposit — that would seed both sides of a one-pixel wall.
    float relocateVolumeTopologySafe(int x, int y, float amount, float momentumX, float momentumY, LiquidCarry *carry = nullptr, bool fromOccupiedCell = false);
    void wakeThermalAt(int x, int y);
    void seedAmbientHeat();
    void addLiquidFill(int index, float dFill, float dHeat);
    void clearSolidifyPending(int index);
    bool addSolidifyPendingKg(int index, SubstanceId id, float kg, float heatJ = 0.0f);
    bool takeSolidifyPendingKg(int index, SubstanceId id, float kg, float *outHeatJ = nullptr);
    bool hasPendingSolid() const { return pendingSolidCellCount > 0; }
    SubstanceId solidifyPendingSubstance(int index) const;
    float solidifyPendingMassKg(int index) const;
    float solidifyPendingSensibleJ(int index) const;
    // Water diagnostics compatibility: pending kg only when the slot is Water.
    float waterFrozenPendingKg(int index) const;
    LiquidCarry takeLiquidCarry(int index, float amount);
    float takeLiquidVolume(int index, float amount);
    float honeyFraction(int index) const; // transport convenience; = liquidComponentFraction(..., HONEY)
    float liquidComponentAmount(int index, RuntimeSubstanceRef id) const;
    float liquidComponentFraction(int index, RuntimeSubstanceRef id) const;
    void setLiquidComponentAmount(int index, RuntimeSubstanceRef id, float amount);
    void addLiquidComponentAmount(int index, RuntimeSubstanceRef id, float delta);
    RuntimeSubstanceRef dominantLiquidRef(int index) const;

    // Built-in compatibility wrappers. Storage is RuntimeSubstanceRef.
    float liquidComponentAmount(int index, SubstanceId id) const {
        return liquidComponentAmount(index, runtimeBuiltIn(id));
    }
    float liquidComponentFraction(int index, SubstanceId id) const {
        return liquidComponentFraction(index, runtimeBuiltIn(id));
    }
    void setLiquidComponentAmount(int index, SubstanceId id, float amount) {
        setLiquidComponentAmount(index, runtimeBuiltIn(id), amount);
    }
    void addLiquidComponentAmount(int index, SubstanceId id, float delta) {
        addLiquidComponentAmount(index, runtimeBuiltIn(id), delta);
    }
    SubstanceId dominantLiquidSubstance(int index) const {
        return runtimeBuiltinId(dominantLiquidRef(index));
    }
    bool liquidCompositionValid(int index) const;
    LiquidComponentView liquidComponents(int index) const;
    template<typename Fn>
    void forEachLiquidComponent(int index, Fn &&fn) const {
        LiquidComponentView view = liquidComponents(index);
        for (int n = 0; n < view.count; ++n)
            fn(view.items[n].id, view.items[n].amount);
    }
    LiquidMixtureProperties mixProperties(int index) const; // density, Cp, k, gamma; viscosity unused
    float mixDensity(int index) const;
    float mixSpecificHeat(int index) const;
    float mixConductivity(int index) const;
    // Occupancy write for chemistry: sets slots + fill from the view, does not
    // rescale heat/dye or invent Water. false = no mutation.
    bool tryCommitLiquidOccupancy(int index, LiquidComponentView const &view);
    void clearEmptyLiquidCell(int index);
    void updateSplashParticles(float dt);
    void spawnSurfaceSpray();
    void drainWaterIntoVoid();
    void wakeRegion(int x0, int y0, int x1, int y1);
    void clearDynamicOccupancy();
    void setDynamicOccupancy(int x, int y, float velX, float velY);
    void applyMovingBoundaryVelocity();
    void displaceFluidFromDynamicSolids();
    void computeSolveRegion();
    void wakeAllFluidChunks();
    void wakeUntrackedLiquid();
    void wakeChunkAtCell(int x, int y, bool resetQuiet = true);
    void rebuildActivityAndMetrics(bool advanceSleep = false);
    void flushPaintDirty();
    void simulationTick();
    void zeroFluidState();
    void clearWorld();
    void resetWorld();
    void loadTestScene(int scene);
    void addSloshImpulse();
    void paintDisc(int cx, int cy, Tool tool, int brushRadius, LiquidPaint paint = {},
        BrushShape shape = BrushShape::Circle, bool eraseWalls = true, bool eraseLiquids = true);
    void paintLine(int x0, int y0, int x1, int y1, Tool tool, int brushRadius, LiquidPaint paint = {},
        BrushShape shape = BrushShape::Circle, bool eraseWalls = true, bool eraseLiquids = true);
    void finalizePaint();
    void runHeadlessBenchmark();
    void runScaleBenchmark();
    void runLiquidBugDiagnostics();
    void runLiquidCompositionDiagnostics();
    void runThreadBenchmark();
    void runAdvectionBenchmark();

private:
    static float bilinear(std::vector<float> const &a, int width, int height, float x, float y);
    static float nearestSample(std::vector<float> const &a, int width, int height, float x, float y);
    static float clampToSourceExtrema(std::vector<float> const &field, int width, int height, float x, float y, float value);
    static float finiteOrZero(float value);
    float sampleU(float worldX, float worldY) const;
    float sampleV(float worldX, float worldY) const;
    float sampleUField(std::vector<float> const &field, float worldX, float worldY) const;
    float sampleVField(std::vector<float> const &field, float worldX, float worldY) const;
    float sampleUNearest(float worldX, float worldY) const;
    float sampleVNearest(float worldX, float worldY) const;
    bool liveUFace(int x, int y) const;
    bool liveVFace(int x, int y) const;
    float uFaceOrWall(int x, int y) const;
    float vFaceOrWall(int x, int y) const;
    void advectVelocitySemiLagrangian(float dt, bool nearest);
    void advectVelocityFirstOrderUpwind(float dt);
    void advectVelocityMacCormack(float dt);
    void advectVelocityBfecc(float dt);
    void commitAdvectedVelocity();
    void sorPressureColor(std::vector<int> const &cells, float invDt, float sor);
    void sorPressureColorParallel(std::vector<int> const &cells, float invDt, float sor);
    float pressureResidual(float invDt, float dt) const;
    float mixViscosity(int index) const;
    float mixSurfaceTension(int index) const;
    void applyCarry(int index, LiquidCarry const &c);
    LiquidCarry extractVolume(int index, float amount);
    static int compositionSlot(int cell, int slot);
    void clearComposition(int index);
    void compactComposition(int index);
    void copyCompositionToNext(int index);
    void commitNextComposition(int index);
    float addComponentUntracked(int index, RuntimeSubstanceRef id, float amount);
    float addNextComponentUntracked(int index, RuntimeSubstanceRef id, float amount);
    bool cellHasDuplicateComponents(int index) const;
    void returnAdvectOverflow(int origin, float vol, LiquidCarry &carry, int x0, int y0, int x1, int y1);
};
