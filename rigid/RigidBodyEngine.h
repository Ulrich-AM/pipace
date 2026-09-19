#pragma once

#include "RigidBodyTypes.h"
#include "fluid/BrushGeom.h"

#include <cstdint>
#include <vector>

struct FluidEngine;

struct RigidBodyEngine {
    std::vector<RigidBody> bodies;
    std::vector<MaterialId> pending;
    std::vector<int> occupant;          // world cell -> body index, or -1
    std::vector<MaterialId> occupantMat;
    std::vector<float> occupantDamage;
    std::vector<float> occupantMoisture;
    std::vector<float> occupantCrack;
    std::vector<RigidContact> lastContacts;
    GrabState grab;
    MaterialId drawMaterial = MATERIAL_WOOD;
    bool debugOverlay = false;
    bool placeAnchored = false;
    bool placeSleeping = false;
    bool placePowder = false;
    int powderParticleSize = 2;
    float gravityScale = 1.0f;
    int velocityContactIters = 8;
    int positionalIters = 1;
    int sleepQuietTicks = 18;
    float sleepLin = 0.16f;
    uint32_t nextId = 1;
    double lastStepMs = 0.0;
    int lastFractureSplits = 0;
    double lastAbsorbed = 0.0;
    double lastDried = 0.0;
    int lastSaturatedSurface = 0;
    int lastDripSites = 0;
    int lastSourceSites = 0;
    int lastOutletSites = 0;
    double lastPendingQueued = 0.0;
    std::vector<float> moistureFlux;

    RigidBodyEngine();

    void clear();
    void paintPendingDisc(int cx, int cy, int brushRadius);
    void paintPendingLine(int x0, int y0, int x1, int y1, int brushRadius);
    void clearPending();
    int commitPending(FluidEngine &fluid);
    // Per requested world cell: whether addSameMaterialWorldCells created a unique
    // source pixel (bodyId + localIndex) that was empty before this call.
    struct SpawnedSourcePixel {
        int worldIndex = -1;
        bool uniqueCreated = false;
        uint32_t bodyId = 0;
        int localIndex = -1;
    };
    int addSameMaterialWorldCells(FluidEngine &fluid, std::vector<int> const &worldCells,
        MaterialId material, float temperatureK,
        std::vector<SpawnedSourcePixel> *uniqueCreated = nullptr);
    void finalizeMaskEdits(FluidEngine &fluid);
    void refreshMassProperties(int bodyIndex);
    // After chemistry: refresh remain-only mass, then connectivity for mask clears.
    void finalizeChemistryEdits(FluidEngine &fluid);

    // World cell -> rigid source pixel (center sample, rotated bodies OK).
    // Static FluidEngine walls are never chemistry solids.
    struct SourcePixel {
        bool valid = false;
        uint32_t bodyId = 0;
        int bodyIndex = -1;
        int localIndex = -1;
        int worldX = 0;
        int worldY = 0;
        MaterialId material = MATERIAL_EMPTY;
        SubstanceId substance = SUBSTANCE_NONE;
        float fraction = 0.0f;
        float heatJ = 0.0f;
    };
    SourcePixel resolveSourcePixel(FluidEngine const &fluid, int worldX, int worldY) const;
    bool sourcePixelStillValid(SourcePixel const &site) const;
    // Transactional remain/heat write. bodyId is the stable identity.
    // newFraction below the removal epsilon clears the mask pixel so
    // ensurePixelState cannot restore it. Does not rebuild mass/occupancy.
    bool commitSourcePixelState(uint32_t bodyId, int localIndex, float newFraction, float newHeatJ);
    void eraseDisc(int cx, int cy, int brushRadius, FluidEngine &fluid,
        BrushShape shape = BrushShape::Circle, bool strictBodies = false,
        std::vector<uint32_t> *strokeSeen = nullptr);
    void eraseLine(int x0, int y0, int x1, int y1, int brushRadius, FluidEngine &fluid,
        BrushShape shape = BrushShape::Circle, bool strictBodies = false,
        std::vector<uint32_t> *strokeSeen = nullptr);
    int bodyAtCell(int x, int y) const;
    void removeBody(int index, FluidEngine &fluid);

    void step(FluidEngine &fluid, float dt);
    void gatherFluidForces(FluidEngine const &fluid);
    void syncOccupancy(FluidEngine &fluid);

    bool beginGrab(float wx, float wy, float strength = 1.0f, bool group = false,
        float groupRadius = 8.0f, bool phantom = false);
    void updateGrabTarget(float wx, float wy, bool strong);
    void endGrab();
    int indexOfId(uint32_t id) const;
    bool grabContainsId(uint32_t id) const;
    bool grabContainsIndex(int index) const;
    void pruneGrab();
    void collectBodiesInRadius(float wx, float wy, float radius, std::vector<int> &out, bool includeAnchored = false) const;
    void applyTouch(float wx, float wy, bool group, float groupRadius, bool toggleAnchor);

    void depositImpactDamage(int bodyIndex, float wx, float wy, float speed, float impulse, float nx = 0.0f, float ny = 0.0f);
    void applyStructuralDamage(int bodyIndex, StructuralDamageRequest const &request);
    void applyMaterialDamage(int bodyIndex, float wx, float wy, float amount, float radius,
        StructuralDamageType type = StructuralDamageType::Generic);
    void applyBondDamage(int bodyIndex, float wx, float wy, float amount, float radius,
        float dirX = 0.0f, float dirY = 0.0f, StructuralDamageType type = StructuralDamageType::Generic);
    void applyLocalStructuralStress(int bodyIndex, float wx, float wy, float amount, float radius,
        float dirX, float dirY, StructuralDamageType type = StructuralDamageType::Generic);
    void modifyLocalStrength(int bodyIndex, float wx, float wy, float radius, float factor);
    bool breakBond(int bodyIndex, int lx, int ly, int nx, int ny);
    bool bondConnects(RigidBody const &b, int lx, int ly, int nx, int ny) const;
    float getLocalCrackFactor(RigidBody const &b, int lx, int ly) const;
    float getLocalCrackFactorWorld(int x, int y) const;
    float effectiveStrength(RigidBody const &b, int localIndex, StructuralDamageType type = StructuralDamageType::Generic) const;
    int brokenBondCountOf(RigidBody const &b) const;
    bool inspectLocalStructure(int wx, int wy, float &materialDamage, float &bondDamage, int &brokenNeighbors,
        float &strength, float &crack, float &wetness) const;
    void processMoisture(FluidEngine &fluid, float dt);
    int totalSolidPixels() const;
    float totalSolidMass() const;
    float totalAbsorbedLiquid() const;
    float totalPendingDrip() const;
    void runMoistureDripDiagnostics(FluidEngine &fluid);

    bool worldCellOccupied(int x, int y) const;
    MaterialId worldCellMaterial(int x, int y) const;
    MatterIdentity worldCellIdentity(int x, int y) const;
    void loadTestScene(FluidEngine &fluid, int scene);
    void runConservationBenchmark(FluidEngine &fluid);
    void runContactDiagnostics(FluidEngine &fluid);
    void runMoistureDiagnostics(FluidEngine &fluid);
    void runSolidDiagnostics(FluidEngine &fluid);

    static void worldToLocal(RigidBody const &b, float wx, float wy, float &lx, float &ly);
    static void localToWorld(RigidBody const &b, float lx, float ly, float &wx, float &wy);
    static bool maskOccupied(RigidBody const &b, int lx, int ly);
    static MaterialId maskMaterial(RigidBody const &b, int lx, int ly);
    static bool occupiesWorldCell(RigidBody const &b, int gx, int gy);
    static bool occupiesWorldCellCenter(RigidBody const &b, int gx, int gy);

private:
    std::vector<OccupancyConflict> occupancyConflicts;
    std::vector<uint8_t> moistureFlags;
    struct MoistureDripSite {
        int localIndex = 0;
        int destX = 0;
        int destY = 0;
        float excess = 0.0f;
        float score = 0.0f;
    };
    std::vector<MoistureDripSite> dripSites;

    void rebuildDerived(RigidBody &b);
    void ensurePixelState(RigidBody &b);
    void computeMassProperties(RigidBody &b);
    void refreshMaterialCache(RigidBody &b);
    void updateAabb(RigidBody &b);
    void applyComShift(RigidBody &b, float oldComX, float oldComY, int padL, int padT);
    void addWorldCellsToBody(int bodyIndex, std::vector<int> const &worldCells, MaterialId material,
        std::vector<SpawnedSourcePixel> *createdPixels = nullptr);
    void shrinkMask(RigidBody &b);
    void refreshBodyAfterMaskEdit(int index, std::vector<RigidBody> &spawned);
    void carveWorldCells(std::vector<int> const &cells, FluidEngine &fluid);
    int attachedBodyForComponent(std::vector<int> const &component, MaterialId requiredMat = MATERIAL_EMPTY) const;
    RigidBody makeBodyFromCells(std::vector<int> const &cells, MaterialId material,
        std::vector<SpawnedSourcePixel> *createdPixels = nullptr);
    void integrate(RigidBody &b, float dt);
    void collectStaticContacts(FluidEngine const &fluid, RigidBody const &b, int index, std::vector<RigidContact> &out) const;
    void collectBodyContacts(std::vector<RigidContact> &out);
    void reduceStaticManifold(std::vector<RigidContact> &contacts) const;
    bool estimateTerrainNormal(FluidEngine const &fluid, int solidX, int solidY,
        float faceNx, float faceNy, float &nx, float &ny) const;
    void solveVelocityContacts(std::vector<RigidContact> &contacts);
    void solvePositionalContacts(std::vector<RigidContact> const &contacts);
    void rasterizeBodyOccupancy(RigidBody const &b, int index, bool writeFluid, bool conservative, FluidEngine *fluid);
    bool bodySupported(FluidEngine const &fluid, RigidBody const &b) const;
    float measureFacePenetration(RigidBody const &b, float faceX, float faceY, float nx, float ny) const;
    void emitContact(std::vector<RigidContact> &out, int bodyA, int bodyB, float x, float y, float nx, float ny, float penetration) const;
    void markSupport(std::vector<RigidContact> const &contacts);
    bool sweptThroughStaticSolid(FluidEngine const &fluid, RigidBody const &b, int index,
        float x0, float y0, float theta0) const;
    void cullBodiesLeftIntoVoid();
    int spawnPattern(FluidEngine &fluid, int ox, int oy, char const *const *rows, int rowCount, MaterialId material);
    void applyGrabForces(float gravity);
    void applyContactDamage(std::vector<RigidContact> const &contacts);
    void splitDirtyBodies();
    void splitBodyByConnectivity(int index, std::vector<RigidBody> &spawned);
    RigidBody packFragment(RigidBody const &parent, std::vector<int> const &localCells, bool keepAnchored);
    void setFragmentVelocity(RigidBody &child, float px, float py, float pvx, float pvy, float pomega);
    float pixelMass(RigidBody const &b, int localIndex) const;
    float pixelWetness(RigidBody const &b, int localIndex) const;
    StructuralBond *bondPtr(RigidBody &b, int lx, int ly, int nx, int ny);
    StructuralBond const *bondPtr(RigidBody const &b, int lx, int ly, int nx, int ny) const;
    bool damageBond(RigidBody &b, int lx, int ly, int nx, int ny, float add, bool *newlyBroken);
    void propagateCracks(RigidBody &b, std::vector<int> const &seedBonds);
    void recountStructure(RigidBody &b);
    float bodyContactFriction(RigidBody const &b) const;
    float bodyContactRestitution(RigidBody const &b) const;
};
