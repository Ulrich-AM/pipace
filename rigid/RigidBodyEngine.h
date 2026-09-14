#pragma once

#include "RigidBodyTypes.h"

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
    uint32_t nextId = 1;
    double lastStepMs = 0.0;
    int lastFractureSplits = 0;
    double lastAbsorbed = 0.0;
    double lastDried = 0.0;

    RigidBodyEngine();

    void clear();
    void paintPendingDisc(int cx, int cy, int brushRadius);
    void paintPendingLine(int x0, int y0, int x1, int y1, int brushRadius);
    void clearPending();
    int commitPending(FluidEngine &fluid);
    void eraseDisc(int cx, int cy, int brushRadius, FluidEngine &fluid);
    void eraseLine(int x0, int y0, int x1, int y1, int brushRadius, FluidEngine &fluid);
    int bodyAtCell(int x, int y) const;
    void removeBody(int index, FluidEngine &fluid);

    void step(FluidEngine &fluid, float dt);
    void gatherFluidForces(FluidEngine const &fluid);
    void syncOccupancy(FluidEngine &fluid);

    bool beginGrab(float wx, float wy);
    void updateGrabTarget(float wx, float wy, bool strong);
    void endGrab();
    int indexOfId(uint32_t id) const;

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

    bool worldCellOccupied(int x, int y) const;
    MaterialId worldCellMaterial(int x, int y) const;
    void loadTestScene(FluidEngine &fluid, int scene);
    void runConservationBenchmark(FluidEngine &fluid);
    void runContactDiagnostics(FluidEngine &fluid);
    void runSolidDiagnostics(FluidEngine &fluid);

    static void worldToLocal(RigidBody const &b, float wx, float wy, float &lx, float &ly);
    static void localToWorld(RigidBody const &b, float lx, float ly, float &wx, float &wy);
    static bool maskOccupied(RigidBody const &b, int lx, int ly);
    static MaterialId maskMaterial(RigidBody const &b, int lx, int ly);
    static bool occupiesWorldCell(RigidBody const &b, int gx, int gy);
    static bool occupiesWorldCellCenter(RigidBody const &b, int gx, int gy);

private:
    std::vector<OccupancyConflict> occupancyConflicts;

    void rebuildDerived(RigidBody &b);
    void ensurePixelState(RigidBody &b);
    void computeMassProperties(RigidBody &b);
    void refreshMaterialCache(RigidBody &b);
    void updateAabb(RigidBody &b);
    void applyComShift(RigidBody &b, float oldComX, float oldComY, int padL, int padT);
    void addWorldCellsToBody(int bodyIndex, std::vector<int> const &worldCells, MaterialId material);
    void shrinkMask(RigidBody &b);
    void refreshBodyAfterMaskEdit(int index, std::vector<RigidBody> &spawned);
    void carveWorldCells(std::vector<int> const &cells, FluidEngine &fluid);
    int attachedBodyForComponent(std::vector<int> const &component) const;
    RigidBody makeBodyFromCells(std::vector<int> const &cells, MaterialId material);
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
