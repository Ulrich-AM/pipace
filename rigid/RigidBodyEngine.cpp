#include "RigidBodyEngine.h"

#include "fluid/FluidEngine.h"
#include "fluid/DiagOutput.h"
#include "thermal/ThermalTypes.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr float kSleepLin = 0.16f;
constexpr float kSleepAng = 0.16f;
constexpr int kSleepTicks = 18;
constexpr float kFriction = 0.42f;
constexpr float kSlop = 0.12f;
constexpr float kBaumgarte = 0.22f;
constexpr float kMaxPen = 1.0f;
constexpr float kMaxLin = 70.0f;
constexpr float kMaxAng = 18.0f;
constexpr float kRestitution = 0.12f;
constexpr float kRestitutionThreshold = 1.5f;
constexpr float kManifoldDot = 0.75f;
constexpr int kMaxManifoldGroups = 3;
constexpr int kMaxManifoldPoints = 2;
constexpr float kImpactSpeedMin = 4.2f;
constexpr float kImpactEnergyScale = 720.0f;
constexpr float kFailDamage = 0.999f;

constexpr int kN4x[4] = {-1, 1, 0, 0};
constexpr int kN4y[4] = {0, 0, -1, 1};

float clampf(float v, float a, float b) { return std::max(a, std::min(b, v)); }

void invTerms(RigidBody const &b, float &invM, float &invI) {
    if (b.immobile()) { invM = 0.0f; invI = 0.0f; }
    else { invM = b.invMass; invI = b.invInertia; }
}

void wakeDormant(RigidBody &b) {
    if (b.anchored) return;
    b.dormant = false;
    b.sleeping = false;
    b.quietTicks = 0;
}

} // namespace

RigidBodyEngine::RigidBodyEngine()
    : pending(static_cast<size_t>(GW * GH), MATERIAL_EMPTY)
    , occupant(static_cast<size_t>(GW * GH), -1)
    , occupantMat(static_cast<size_t>(GW * GH), MATERIAL_EMPTY)
    , occupantDamage(static_cast<size_t>(GW * GH), 0.0f)
    , occupantMoisture(static_cast<size_t>(GW * GH), 0.0f)
    , occupantCrack(static_cast<size_t>(GW * GH), 0.0f)
{
    bodies.reserve(64);
    lastContacts.reserve(128);
}

void RigidBodyEngine::clear() {
    bodies.clear();
    clearPending();
    std::fill(occupant.begin(), occupant.end(), -1);
    std::fill(occupantMat.begin(), occupantMat.end(), MATERIAL_EMPTY);
    std::fill(occupantDamage.begin(), occupantDamage.end(), 0.0f);
    std::fill(occupantMoisture.begin(), occupantMoisture.end(), 0.0f);
    std::fill(occupantCrack.begin(), occupantCrack.end(), 0.0f);
    lastContacts.clear();
    grab = GrabState{};
    nextId = 1;
    lastFractureSplits = 0;
    lastAbsorbed = 0.0;
    lastDried = 0.0;
}

void RigidBodyEngine::clearPending() {
    std::fill(pending.begin(), pending.end(), MATERIAL_EMPTY);
}

void RigidBodyEngine::paintPendingDisc(int cx, int cy, int brushRadius) {
    int r2 = brushRadius * brushRadius;
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) > r2) continue;
            pending[static_cast<size_t>(FluidEngine::ci(x, y))] = drawMaterial;
        }
}

void RigidBodyEngine::paintPendingLine(int x0, int y0, int x1, int y1, int brushRadius) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1, dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1, error = dx + dy;
    for (;;) {
        paintPendingDisc(x0, y0, brushRadius);
        if (x0 == x1 && y0 == y1) break;
        int twice = 2 * error;
        if (twice >= dy) { error += dy; x0 += sx; }
        if (twice <= dx) { error += dx; y0 += sy; }
    }
}

void RigidBodyEngine::worldToLocal(RigidBody const &b, float wx, float wy, float &lx, float &ly) {
    float dx = wx - b.x, dy = wy - b.y;
    float c = std::cos(b.theta), s = std::sin(b.theta);
    lx = c * dx + s * dy + b.comLocalX;
    ly = -s * dx + c * dy + b.comLocalY;
}

void RigidBodyEngine::localToWorld(RigidBody const &b, float lx, float ly, float &wx, float &wy) {
    float c = std::cos(b.theta), s = std::sin(b.theta);
    float rx = lx - b.comLocalX, ry = ly - b.comLocalY;
    wx = b.x + c * rx - s * ry;
    wy = b.y + s * rx + c * ry;
}

bool RigidBodyEngine::maskOccupied(RigidBody const &b, int lx, int ly) {
    if (lx < 0 || ly < 0 || lx >= b.maskW || ly >= b.maskH) return false;
    return b.mask[static_cast<size_t>(ly * b.maskW + lx)] != MATERIAL_EMPTY;
}

MaterialId RigidBodyEngine::maskMaterial(RigidBody const &b, int lx, int ly) {
    if (lx < 0 || ly < 0 || lx >= b.maskW || ly >= b.maskH) return MATERIAL_EMPTY;
    return b.mask[static_cast<size_t>(ly * b.maskW + lx)];
}

bool RigidBodyEngine::occupiesWorldCell(RigidBody const &b, int gx, int gy) {
    // Conservative 5-tap for fluid exclusion / render so water cannot leak through.
    constexpr float ox[5] = {0.5f, 0.12f, 0.88f, 0.12f, 0.88f};
    constexpr float oy[5] = {0.5f, 0.12f, 0.12f, 0.88f, 0.88f};
    for (int s = 0; s < 5; ++s) {
        float lx, ly;
        worldToLocal(b, static_cast<float>(gx) + ox[s], static_cast<float>(gy) + oy[s], lx, ly);
        int ix = static_cast<int>(std::floor(lx));
        int iy = static_cast<int>(std::floor(ly));
        if (maskOccupied(b, ix, iy)) return true;
    }
    return false;
}

bool RigidBodyEngine::occupiesWorldCellCenter(RigidBody const &b, int gx, int gy) {
    // Stricter center sample for collision so grazing floor cells do not become contacts.
    float lx, ly;
    worldToLocal(b, static_cast<float>(gx) + 0.5f, static_cast<float>(gy) + 0.5f, lx, ly);
    return maskOccupied(b, static_cast<int>(std::floor(lx)), static_cast<int>(std::floor(ly)));
}

bool RigidBodyEngine::estimateTerrainNormal(FluidEngine const &fluid, int solidX, int solidY,
    float faceNx, float faceNy, float &nx, float &ny) const {
    float gx = 0.0f, gy = 0.0f;
    int hits = 0;
    for (int oy = -2; oy <= 2; ++oy) {
        for (int ox = -2; ox <= 2; ++ox) {
            if (!fluid.isPaintedSolid(solidX + ox, solidY + oy)) continue;
            float w = (ox == 0 && oy == 0) ? 1.5f : 1.0f;
            gx += w * static_cast<float>(ox);
            gy += w * static_cast<float>(oy);
            ++hits;
        }
    }
    if (hits < 3) {
        nx = faceNx; ny = faceNy;
        return false;
    }
    float len = std::sqrt(gx * gx + gy * gy);
    if (len < 1e-3f) {
        nx = faceNx; ny = faceNy;
        return false;
    }
    nx = -gx / len;
    ny = -gy / len;
    if (nx * faceNx + ny * faceNy < 0.20f) {
        nx = faceNx; ny = faceNy;
        return false;
    }
    return true;
}

float RigidBodyEngine::measureFacePenetration(RigidBody const &b, float faceX, float faceY, float nx, float ny) const {
    float depth = 0.0f;
    for (int i = 1; i <= 10; ++i) {
        float t = 0.10f * static_cast<float>(i);
        float lx, ly;
        worldToLocal(b, faceX - nx * t, faceY - ny * t, lx, ly);
        if (maskOccupied(b, static_cast<int>(std::floor(lx)), static_cast<int>(std::floor(ly))))
            depth = t;
        else
            break;
    }
    return std::min(depth, kMaxPen);
}

void RigidBodyEngine::emitContact(std::vector<RigidContact> &out, int bodyA, int bodyB,
    float x, float y, float nx, float ny, float penetration) const {
    float len = std::sqrt(nx * nx + ny * ny);
    if (len < 1e-5f) return;
    nx /= len; ny /= len;
    RigidContact c;
    c.bodyA = bodyA; c.bodyB = bodyB;
    c.x = x; c.y = y; c.nx = nx; c.ny = ny;
    c.rAx = x - bodies[static_cast<size_t>(bodyA)].x;
    c.rAy = y - bodies[static_cast<size_t>(bodyA)].y;
    if (bodyB >= 0) {
        c.rBx = x - bodies[static_cast<size_t>(bodyB)].x;
        c.rBy = y - bodies[static_cast<size_t>(bodyB)].y;
    }
    c.penetration = std::min(std::max(penetration, 0.0f), kMaxPen);
    out.push_back(c);
}

void RigidBodyEngine::rebuildDerived(RigidBody &b) {
    b.runs.clear();
    b.occupiedLocal.clear();
    for (int y = 0; y < b.maskH; ++y) {
        int x = 0;
        while (x < b.maskW) {
            if (!maskOccupied(b, x, y)) { ++x; continue; }
            int x0 = x;
            while (x < b.maskW && maskOccupied(b, x, y)) ++x;
            b.runs.push_back({y, x0, x - 1});
            for (int xx = x0; xx < x; ++xx) b.occupiedLocal.push_back(y * b.maskW + xx);
        }
    }
    ensurePixelState(b);
}

void RigidBodyEngine::ensurePixelState(RigidBody &b) {
    size_t n = b.mask.size();
    auto fitFloat = [&](std::vector<float> &v) {
        if (v.size() == n) return;
        std::vector<float> next(n, 0.0f);
        size_t copy = std::min(n, v.size());
        for (size_t i = 0; i < copy; ++i) next[i] = v[i];
        v.swap(next);
    };
    auto fitBond = [&](std::vector<StructuralBond> &v) {
        if (v.size() == n) return;
        std::vector<StructuralBond> next(n);
        size_t copy = std::min(n, v.size());
        for (size_t i = 0; i < copy; ++i) next[i] = v[i];
        v.swap(next);
    };
    fitFloat(b.materialDamage);
    fitFloat(b.moisture);
    fitFloat(b.heat);
    fitBond(b.bondsRight);
    fitBond(b.bondsDown);
}

float RigidBodyEngine::pixelMass(RigidBody const &b, int localIndex) const {
    if (localIndex < 0 || localIndex >= static_cast<int>(b.mask.size())) return 0.0f;
    MaterialId id = b.mask[static_cast<size_t>(localIndex)];
    if (id == MATERIAL_EMPTY) return 0.0f;
    MaterialDefinition const &mat = materialDef(id);
    float water = 0.0f;
    if (localIndex < static_cast<int>(b.moisture.size())) water = std::max(0.0f, b.moisture[static_cast<size_t>(localIndex)]);
    return mat.density + water * (1.0f + mat.wetDensityContribution);
}

float RigidBodyEngine::pixelWetness(RigidBody const &b, int localIndex) const {
    if (localIndex < 0 || localIndex >= static_cast<int>(b.mask.size())) return 0.0f;
    MaterialDefinition const &mat = materialDef(b.mask[static_cast<size_t>(localIndex)]);
    if (mat.moistureCapacity <= 1.0e-8f) return 0.0f;
    float mst = 0.0f;
    if (localIndex < static_cast<int>(b.moisture.size())) mst = std::max(0.0f, b.moisture[static_cast<size_t>(localIndex)]);
    return clampf(mst / mat.moistureCapacity, 0.0f, 1.0f);
}

float RigidBodyEngine::effectiveStrength(RigidBody const &b, int localIndex, StructuralDamageType type) const {
    if (localIndex < 0 || localIndex >= static_cast<int>(b.mask.size())) return 0.0f;
    MaterialId id = b.mask[static_cast<size_t>(localIndex)];
    if (id == MATERIAL_EMPTY) return 0.0f;
    MaterialDefinition const &mat = materialDef(id);
    float base = mat.tensileStrength;
    if (type == StructuralDamageType::Compression) base = mat.compressiveStrength;
    else if (type == StructuralDamageType::Shear) base = mat.shearStrength;
    else if (type == StructuralDamageType::Impact)
        base = 0.55f * mat.tensileStrength + 0.45f * mat.compressiveStrength;
    float wet = pixelWetness(b, localIndex);
    base *= (1.0f + (mat.wetStrengthMultiplier - 1.0f) * wet);
    float md = 0.0f;
    if (localIndex < static_cast<int>(b.materialDamage.size()))
        md = clampf(b.materialDamage[static_cast<size_t>(localIndex)], 0.0f, 1.0f);
    base *= (1.0f - md * (0.50f + 0.40f * mat.brittleness));
    float tough = mat.fractureToughness * (1.0f + (mat.wetFractureToughnessMultiplier - 1.0f) * wet);
    tough *= (1.0f - 0.35f * md);
    base *= (0.55f + 0.45f * clampf(tough, 0.05f, 2.0f));
    return std::max(base, 0.04f);
}

void RigidBodyEngine::computeMassProperties(RigidBody &b) {
    ensurePixelState(b);
    double mass = 0.0, cx = 0.0, cy = 0.0;
    for (int i : b.occupiedLocal) {
        int lx = i % b.maskW, ly = i / b.maskW;
        float m = pixelMass(b, i);
        mass += m;
        cx += m * (lx + 0.5);
        cy += m * (ly + 0.5);
    }
    if (mass < 1e-8) {
        b.mass = 1.0f; b.invMass = 1.0f; b.inertia = 1.0f; b.invInertia = 1.0f;
        b.comLocalX = 0.5f * b.maskW; b.comLocalY = 0.5f * b.maskH;
        refreshMaterialCache(b);
        return;
    }
    b.comLocalX = static_cast<float>(cx / mass);
    b.comLocalY = static_cast<float>(cy / mass);
    double inertia = 0.0;
    for (int i : b.occupiedLocal) {
        int lx = i % b.maskW, ly = i / b.maskW;
        float m = pixelMass(b, i);
        double dx = (lx + 0.5) - b.comLocalX;
        double dy = (ly + 0.5) - b.comLocalY;
        inertia += m * (dx * dx + dy * dy + 1.0 / 6.0);
    }
    b.mass = static_cast<float>(mass);
    b.inertia = std::max(static_cast<float>(inertia), b.mass * 0.12f);
    b.invMass = 1.0f / b.mass;
    b.invInertia = 1.0f / b.inertia;
    refreshMaterialCache(b);
}

void RigidBodyEngine::refreshMaterialCache(RigidBody &b) {
    ensurePixelState(b);
    double fr = 0.0, rest = 0.0, wet = 0.0, absorbed = 0.0, w = 0.0, maxD = 0.0;
    bool absorbent = false;
    for (int i : b.occupiedLocal) {
        MaterialId id = b.mask[static_cast<size_t>(i)];
        MaterialDefinition const &mat = materialDef(id);
        float wetness = pixelWetness(b, i);
        float mu = mat.friction * (1.0f + (mat.wetFrictionMultiplier - 1.0f) * wetness);
        fr += mu;
        rest += mat.restitution;
        wet += wetness;
        absorbed += b.moisture[static_cast<size_t>(i)];
        maxD = std::max(maxD, static_cast<double>(b.materialDamage[static_cast<size_t>(i)]));
        if (materialIsAbsorbent(id)) absorbent = true;
        w += 1.0;
    }
    if (w < 1.0e-8) {
        b.cachedFriction = 0.42f;
        b.cachedRestitution = 0.12f;
        b.cachedWetness = 0.0f;
        b.absorbedLiquid = 0.0f;
        b.maxDamage = 0.0f;
        b.maxBondDamage = 0.0f;
        b.brokenBondCount = 0;
        b.moistureActive = false;
        return;
    }
    b.cachedFriction = static_cast<float>(fr / w);
    b.cachedRestitution = static_cast<float>(rest / w);
    b.cachedWetness = static_cast<float>(wet / w);
    b.absorbedLiquid = static_cast<float>(absorbed);
    b.maxDamage = static_cast<float>(maxD);
    b.moistureActive = absorbent;
    recountStructure(b);
}

void RigidBodyEngine::updateAabb(RigidBody &b) {
    float x0 = 1e9f, y0 = 1e9f, x1 = -1e9f, y1 = -1e9f;
    float corners[4][2] = {
        {0.0f, 0.0f}, {static_cast<float>(b.maskW), 0.0f},
        {0.0f, static_cast<float>(b.maskH)}, {static_cast<float>(b.maskW), static_cast<float>(b.maskH)}
    };
    for (auto &c : corners) {
        float wx, wy; localToWorld(b, c[0], c[1], wx, wy);
        x0 = std::min(x0, wx); y0 = std::min(y0, wy);
        x1 = std::max(x1, wx); y1 = std::max(y1, wy);
    }
    b.prevAabbX0 = b.aabbX0; b.prevAabbY0 = b.aabbY0; b.prevAabbX1 = b.aabbX1; b.prevAabbY1 = b.aabbY1;
    b.aabbX0 = x0 - 1.0f; b.aabbY0 = y0 - 1.0f; b.aabbX1 = x1 + 1.0f; b.aabbY1 = y1 + 1.0f;
}

RigidBody RigidBodyEngine::makeBodyFromCells(std::vector<int> const &cells, MaterialId material) {
    RigidBody b;
    int minX = GW, minY = GH, maxX = -1, maxY = -1;
    for (int index : cells) {
        int x = index % GW, y = index / GW;
        minX = std::min(minX, x); maxX = std::max(maxX, x);
        minY = std::min(minY, y); maxY = std::max(maxY, y);
    }
    b.maskW = maxX - minX + 1;
    b.maskH = maxY - minY + 1;
    size_t nMask = static_cast<size_t>(b.maskW * b.maskH);
    b.mask.assign(nMask, MATERIAL_EMPTY);
    b.materialDamage.assign(nMask, 0.0f);
    b.moisture.assign(nMask, 0.0f);
    b.heat.assign(nMask, 0.0f);
    b.bondsRight.assign(nMask, StructuralBond{});
    b.bondsDown.assign(nMask, StructuralBond{});
    for (int index : cells) {
        int x = index % GW, y = index / GW;
        b.mask[static_cast<size_t>((y - minY) * b.maskW + (x - minX))] = material;
    }
    rebuildDerived(b);
    computeMassProperties(b);
    if (b.heat.size() != b.mask.size()) b.heat.assign(b.mask.size(), 0.0f);
    for (int li : b.occupiedLocal) {
        MaterialId id = b.mask[static_cast<size_t>(li)];
        float cap = thermalCapacity(massKg(materialDef(id).density, 1.0f), thermalForMaterial(id).specificHeat);
        b.heat[static_cast<size_t>(li)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
    }
    b.x = static_cast<float>(minX) + b.comLocalX;
    b.y = static_cast<float>(minY) + b.comLocalY;
    b.theta = 0.0f;
    b.id = nextId++;
    b.anchored = placeAnchored;
    b.dormant = placeSleeping && !placeAnchored;
    if (b.immobile()) {
        b.vx = b.vy = b.omega = 0.0f;
        b.sleeping = true;
    }
    updateAabb(b);
    b.prevAabbX0 = b.aabbX0; b.prevAabbY0 = b.aabbY0; b.prevAabbX1 = b.aabbX1; b.prevAabbY1 = b.aabbY1;
    return b;
}

void RigidBodyEngine::applyComShift(RigidBody &b, float oldComX, float oldComY, int padL, int padT) {
    float c = std::cos(b.theta), s = std::sin(b.theta);
    float dlx = b.comLocalX - (oldComX + static_cast<float>(padL));
    float dly = b.comLocalY - (oldComY + static_cast<float>(padT));
    b.x += c * dlx - s * dly;
    b.y += s * dlx + c * dly;
}

void RigidBodyEngine::shrinkMask(RigidBody &b) {
    if (b.occupiedLocal.empty()) return;
    int minX = b.maskW, minY = b.maskH, maxX = -1, maxY = -1;
    for (int i : b.occupiedLocal) {
        int lx = i % b.maskW, ly = i / b.maskW;
        minX = std::min(minX, lx); maxX = std::max(maxX, lx);
        minY = std::min(minY, ly); maxY = std::max(maxY, ly);
    }
    if (minX == 0 && minY == 0 && maxX == b.maskW - 1 && maxY == b.maskH - 1) return;
    int newW = maxX - minX + 1, newH = maxY - minY + 1;
    std::vector<MaterialId> packed(static_cast<size_t>(newW * newH), MATERIAL_EMPTY);
    std::vector<float> packedD(static_cast<size_t>(newW * newH), 0.0f);
    std::vector<float> packedM(static_cast<size_t>(newW * newH), 0.0f);
    std::vector<float> packedH(static_cast<size_t>(newW * newH), 0.0f);
    std::vector<StructuralBond> packedR(static_cast<size_t>(newW * newH));
    std::vector<StructuralBond> packedDn(static_cast<size_t>(newW * newH));
    ensurePixelState(b);
    for (int y = minY; y <= maxY; ++y) for (int x = minX; x <= maxX; ++x) {
        size_t src = static_cast<size_t>(y * b.maskW + x);
        size_t dst = static_cast<size_t>((y - minY) * newW + (x - minX));
        packed[dst] = b.mask[src];
        packedD[dst] = b.materialDamage[src];
        packedM[dst] = b.moisture[src];
        packedH[dst] = (src < b.heat.size()) ? b.heat[src] : 0.0f;
        if (x < maxX) packedR[dst] = b.bondsRight[src];
        if (y < maxY) packedDn[dst] = b.bondsDown[src];
    }
    float oldComX = b.comLocalX, oldComY = b.comLocalY;
    b.mask.swap(packed);
    b.materialDamage.swap(packedD);
    b.moisture.swap(packedM);
    b.heat.swap(packedH);
    b.bondsRight.swap(packedR);
    b.bondsDown.swap(packedDn);
    b.maskW = newW;
    b.maskH = newH;
    rebuildDerived(b);
    computeMassProperties(b);
    applyComShift(b, oldComX, oldComY, -minX, -minY);
}

void RigidBodyEngine::addWorldCellsToBody(int bodyIndex, std::vector<int> const &worldCells, MaterialId material) {
    if (bodyIndex < 0 || bodyIndex >= static_cast<int>(bodies.size()) || worldCells.empty()) return;
    RigidBody &b = bodies[static_cast<size_t>(bodyIndex)];
    struct Local { int x, y; };
    std::vector<Local> locals;
    locals.reserve(worldCells.size());
    int minLx = 0, minLy = 0, maxLx = std::max(0, b.maskW - 1), maxLy = std::max(0, b.maskH - 1);
    bool any = b.maskW > 0 && b.maskH > 0;
    for (int index : worldCells) {
        int gx = index % GW, gy = index / GW;
        float lx, ly;
        worldToLocal(b, static_cast<float>(gx) + 0.5f, static_cast<float>(gy) + 0.5f, lx, ly);
        int ix = static_cast<int>(std::floor(lx));
        int iy = static_cast<int>(std::floor(ly));
        locals.push_back({ix, iy});
        if (!any) { minLx = maxLx = ix; minLy = maxLy = iy; any = true; }
        else {
            minLx = std::min(minLx, ix); maxLx = std::max(maxLx, ix);
            minLy = std::min(minLy, iy); maxLy = std::max(maxLy, iy);
        }
    }
    int padL = std::max(0, -minLx);
    int padT = std::max(0, -minLy);
    int padR = std::max(0, maxLx - (b.maskW - 1));
    int padB = std::max(0, maxLy - (b.maskH - 1));
    float oldComX = b.comLocalX, oldComY = b.comLocalY;
    if (padL || padT || padR || padB) {
        int newW = b.maskW + padL + padR;
        int newH = b.maskH + padT + padB;
        std::vector<MaterialId> grown(static_cast<size_t>(std::max(1, newW * newH)), MATERIAL_EMPTY);
        std::vector<float> grownD(static_cast<size_t>(std::max(1, newW * newH)), 0.0f);
        std::vector<float> grownM(static_cast<size_t>(std::max(1, newW * newH)), 0.0f);
        std::vector<float> grownH(static_cast<size_t>(std::max(1, newW * newH)), 0.0f);
        std::vector<StructuralBond> grownR(static_cast<size_t>(std::max(1, newW * newH)));
        std::vector<StructuralBond> grownDn(static_cast<size_t>(std::max(1, newW * newH)));
        ensurePixelState(b);
        for (int y = 0; y < b.maskH; ++y) for (int x = 0; x < b.maskW; ++x) {
            size_t src = static_cast<size_t>(y * b.maskW + x);
            size_t dst = static_cast<size_t>((y + padT) * newW + (x + padL));
            grown[dst] = b.mask[src];
            grownD[dst] = b.materialDamage[src];
            grownM[dst] = b.moisture[src];
            grownH[dst] = (src < b.heat.size()) ? b.heat[src] : 0.0f;
            grownR[dst] = b.bondsRight[src];
            grownDn[dst] = b.bondsDown[src];
        }
        b.mask.swap(grown);
        b.materialDamage.swap(grownD);
        b.moisture.swap(grownM);
        b.heat.swap(grownH);
        b.bondsRight.swap(grownR);
        b.bondsDown.swap(grownDn);
        b.maskW = newW;
        b.maskH = newH;
    }
    for (Local const &p : locals) {
        int x = p.x + padL, y = p.y + padT;
        if (x < 0 || y < 0 || x >= b.maskW || y >= b.maskH) continue;
        size_t i = static_cast<size_t>(y * b.maskW + x);
        if (b.mask[i] == MATERIAL_EMPTY) b.mask[i] = material;
        if (b.heat.size() == b.mask.size() && b.heat[i] == 0.0f && b.mask[i] != MATERIAL_EMPTY) {
            float cap = thermalCapacity(massKg(materialDef(b.mask[i]).density, 1.0f),
                thermalForMaterial(b.mask[i]).specificHeat);
            b.heat[i] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
        }
    }
    rebuildDerived(b);
    if (b.occupiedLocal.empty()) return;
    computeMassProperties(b);
    applyComShift(b, oldComX, oldComY, padL, padT);
    if (b.dormant) wakeDormant(b);
    else {
        b.sleeping = false;
        b.quietTicks = 0;
    }
    updateAabb(b);
}

int RigidBodyEngine::attachedBodyForComponent(std::vector<int> const &component) const {
    constexpr int dx[4] = {-1, 1, 0, 0};
    constexpr int dy[4] = {0, 0, -1, 1};
    int best = -1, bestCount = 0;
    std::vector<int> votes(bodies.size(), 0);
    for (int index : component) {
        int x = index % GW, y = index / GW;
        auto consider = [&](int ox, int oy) {
            if (!FluidEngine::inside(ox, oy)) return;
            int body = occupant[static_cast<size_t>(FluidEngine::ci(ox, oy))];
            if (body < 0 || body >= static_cast<int>(bodies.size())) return;
            votes[static_cast<size_t>(body)] += 1;
            if (votes[static_cast<size_t>(body)] > bestCount) {
                bestCount = votes[static_cast<size_t>(body)];
                best = body;
            }
        };
        consider(x, y);
        for (int n = 0; n < 4; ++n) consider(x + dx[n], y + dy[n]);
    }
    return best;
}

void RigidBodyEngine::refreshBodyAfterMaskEdit(int index, std::vector<RigidBody> &spawned) {
    RigidBody &b = bodies[static_cast<size_t>(index)];
    rebuildDerived(b);
    if (b.occupiedLocal.empty()) return;
    std::vector<uint8_t> seen(static_cast<size_t>(b.maskW * b.maskH), 0);
    constexpr int dx[4] = {-1, 1, 0, 0};
    constexpr int dy[4] = {0, 0, -1, 1};
    std::vector<std::vector<int>> comps;
    for (int start : b.occupiedLocal) {
        if (seen[static_cast<size_t>(start)]) continue;
        std::vector<int> comp;
        std::vector<int> stack{start};
        seen[static_cast<size_t>(start)] = 1;
        while (!stack.empty()) {
            int i = stack.back(); stack.pop_back();
            comp.push_back(i);
            int lx = i % b.maskW, ly = i / b.maskW;
            for (int n = 0; n < 4; ++n) {
                int nx = lx + dx[n], ny = ly + dy[n];
                if (nx < 0 || ny < 0 || nx >= b.maskW || ny >= b.maskH) continue;
                int ni = ny * b.maskW + nx;
                if (seen[static_cast<size_t>(ni)] || !maskOccupied(b, nx, ny)) continue;
                seen[static_cast<size_t>(ni)] = 1;
                stack.push_back(ni);
            }
        }
        comps.push_back(std::move(comp));
    }
    if (comps.size() == 1) {
        float oldComX = b.comLocalX, oldComY = b.comLocalY;
        computeMassProperties(b);
        applyComShift(b, oldComX, oldComY, 0, 0);
        shrinkMask(b);
        if (b.anchored) {
            b.vx = b.vy = b.omega = 0.0f;
            b.sleeping = true;
        } else if (b.dormant) {
            wakeDormant(b);
        } else {
            b.sleeping = false;
            b.quietTicks = 0;
        }
        updateAabb(b);
        return;
    }
    size_t keep = 0;
    size_t keepCount = comps[0].size();
    for (size_t c = 1; c < comps.size(); ++c) {
        if (comps[c].size() > keepCount) {
            keep = c;
            keepCount = comps[c].size();
        }
    }
    uint32_t keepId = b.id;
    bool parentAnchored = b.anchored;
    RigidBody kept = packFragment(b, comps[keep], parentAnchored);
    kept.id = keepId;
    for (size_t c = 0; c < comps.size(); ++c) {
        if (c == keep) continue;
        spawned.push_back(packFragment(b, comps[c], false));
    }
    b = std::move(kept);
}

int RigidBodyEngine::commitPending(FluidEngine &fluid) {
    if (placePowder) {
        int s = std::clamp(powderParticleSize, 1, 8);
        std::map<uint64_t, std::pair<MaterialId, std::vector<int>>> grains;
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int i = FluidEngine::ci(x, y);
            MaterialId mat = pending[static_cast<size_t>(i)];
            if (mat == MATERIAL_EMPTY) continue;
            pending[static_cast<size_t>(i)] = MATERIAL_EMPTY;
            if (fluid.isStaticSolid(x, y)) continue;
            if (occupant[static_cast<size_t>(i)] >= 0) continue;
            int gx = x / s;
            int gy = y / s;
            uint64_t key = (static_cast<uint64_t>(mat) << 32)
                | (static_cast<uint64_t>(static_cast<uint32_t>(gx) & 0xffffu) << 16)
                | static_cast<uint64_t>(static_cast<uint32_t>(gy) & 0xffffu);
            auto &grain = grains[key];
            grain.first = mat;
            grain.second.push_back(i);
        }
        int spawned = 0;
        for (auto &entry : grains) {
            if (entry.second.second.empty()) continue;
            bodies.push_back(makeBodyFromCells(entry.second.second, entry.second.first));
            ++spawned;
        }
        if (spawned > 0) syncOccupancy(fluid);
        return spawned;
    }
    std::vector<uint8_t> seen(static_cast<size_t>(GW * GH), 0);
    int spawned = 0;
    constexpr int dx[4] = {-1, 1, 0, 0};
    constexpr int dy[4] = {0, 0, -1, 1};
    bool changed = false;
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int start = FluidEngine::ci(x, y);
        if (pending[static_cast<size_t>(start)] == MATERIAL_EMPTY || seen[static_cast<size_t>(start)]) continue;
        if (fluid.isStaticSolid(x, y)) { pending[static_cast<size_t>(start)] = MATERIAL_EMPTY; continue; }
        MaterialId mat = pending[static_cast<size_t>(start)];
        std::vector<int> component;
        std::vector<int> stack{start};
        seen[static_cast<size_t>(start)] = 1;
        while (!stack.empty()) {
            int i = stack.back(); stack.pop_back();
            int cx = i % GW, cy = i / GW;
            if (fluid.isStaticSolid(cx, cy)) continue;
            component.push_back(i);
            for (int n = 0; n < 4; ++n) {
                int nx = cx + dx[n], ny = cy + dy[n];
                if (!FluidEngine::inside(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (seen[static_cast<size_t>(ni)] || pending[static_cast<size_t>(ni)] == MATERIAL_EMPTY) continue;
                seen[static_cast<size_t>(ni)] = 1;
                stack.push_back(ni);
            }
        }
        if (component.empty()) continue;
        std::vector<int> fresh;
        fresh.reserve(component.size());
        int attach = attachedBodyForComponent(component);
        for (int index : component) {
            if (occupant[static_cast<size_t>(index)] >= 0) continue;
            fresh.push_back(index);
        }
        if (fresh.empty()) continue;
        if (attach >= 0) {
            addWorldCellsToBody(attach, fresh, mat);
            changed = true;
        } else {
            bodies.push_back(makeBodyFromCells(fresh, mat));
            ++spawned;
            changed = true;
        }
    }
    clearPending();
    if (changed) syncOccupancy(fluid);
    return spawned;
}

int RigidBodyEngine::bodyAtCell(int x, int y) const {
    if (!FluidEngine::inside(x, y)) return -1;
    return occupant[static_cast<size_t>(FluidEngine::ci(x, y))];
}

void RigidBodyEngine::removeBody(int index, FluidEngine &fluid) {
    if (index < 0 || index >= static_cast<int>(bodies.size())) return;
    if (grab.active && bodies[static_cast<size_t>(index)].id == grab.bodyId) grab = GrabState{};
    bodies.erase(bodies.begin() + index);
    lastContacts.clear();
    syncOccupancy(fluid);
}

void RigidBodyEngine::carveWorldCells(std::vector<int> const &cells, FluidEngine &fluid) {
    if (cells.empty() || bodies.empty()) return;
    std::vector<uint8_t> hit(bodies.size(), 0);
    bool any = false;
    constexpr float ox[5] = {0.5f, 0.12f, 0.88f, 0.12f, 0.88f};
    constexpr float oy[5] = {0.5f, 0.12f, 0.12f, 0.88f, 0.88f};
    auto clearLocal = [&](RigidBody &b, int ix, int iy) {
        if (!maskOccupied(b, ix, iy)) return;
        b.mask[static_cast<size_t>(iy * b.maskW + ix)] = MATERIAL_EMPTY;
        any = true;
    };
    for (int index : cells) {
        if (index < 0 || index >= GW * GH) continue;
        int x = index % GW, y = index / GW;
        int body = occupant[static_cast<size_t>(index)];
        if (body < 0 || body >= static_cast<int>(bodies.size())) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                if (occupiesWorldCell(bodies[i], x, y)) { body = static_cast<int>(i); break; }
            }
        }
        if (body < 0) continue;
        RigidBody &b = bodies[static_cast<size_t>(body)];
        for (int s = 0; s < 5; ++s) {
            float lx, ly;
            worldToLocal(b, static_cast<float>(x) + ox[s], static_cast<float>(y) + oy[s], lx, ly);
            clearLocal(b, static_cast<int>(std::floor(lx)), static_cast<int>(std::floor(ly)));
        }
        hit[static_cast<size_t>(body)] = 1;
    }
    if (!any) return;
    std::vector<RigidBody> spawned;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (!hit[i]) continue;
        refreshBodyAfterMaskEdit(static_cast<int>(i), spawned);
        if (bodies[i].occupiedLocal.empty()) continue;
        if (bodies[i].anchored) {
            bodies[i].vx = bodies[i].vy = bodies[i].omega = 0.0f;
            bodies[i].sleeping = true;
        } else if (bodies[i].dormant) {
            wakeDormant(bodies[i]);
        } else {
            bodies[i].sleeping = false;
            bodies[i].quietTicks = 0;
        }
    }
    size_t write = 0;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].occupiedLocal.empty()) continue;
        if (write != i) bodies[write] = std::move(bodies[i]);
        ++write;
    }
    bodies.resize(write);
    for (RigidBody &extra : spawned) bodies.push_back(std::move(extra));
    lastContacts.clear();
    syncOccupancy(fluid);
}

void RigidBodyEngine::eraseDisc(int cx, int cy, int brushRadius, FluidEngine &fluid) {
    int r2 = brushRadius * brushRadius;
    std::vector<int> cells;
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) > r2) continue;
            cells.push_back(FluidEngine::ci(x, y));
        }
    carveWorldCells(cells, fluid);
}

void RigidBodyEngine::eraseLine(int x0, int y0, int x1, int y1, int brushRadius, FluidEngine &fluid) {
    int r2 = brushRadius * brushRadius;
    std::vector<int> cells;
    auto stamp = [&](int cx, int cy) {
        for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
            for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
                if (!FluidEngine::inside(x, y)) continue;
                if ((x - cx) * (x - cx) + (y - cy) * (y - cy) > r2) continue;
                cells.push_back(FluidEngine::ci(x, y));
            }
    };
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1, dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1, error = dx + dy;
    for (;;) {
        stamp(x0, y0);
        if (x0 == x1 && y0 == y1) break;
        int twice = 2 * error;
        if (twice >= dy) { error += dy; x0 += sx; }
        if (twice <= dx) { error += dx; y0 += sy; }
    }
    carveWorldCells(cells, fluid);
}

void RigidBodyEngine::integrate(RigidBody &b, float dt) {
    b.x += b.vx * dt;
    b.y += b.vy * dt;
    b.theta += b.omega * dt;
    float twoPi = 6.2831853f;
    if (b.theta > twoPi || b.theta < -twoPi) b.theta = std::fmod(b.theta, twoPi);
}

bool RigidBodyEngine::bodySupported(FluidEngine const &fluid, RigidBody const &b) const {
    int x0 = std::max(0, static_cast<int>(std::floor(b.aabbX0)));
    int y0 = std::max(0, static_cast<int>(std::floor(b.aabbY0)));
    int x1 = std::min(GW - 1, static_cast<int>(std::ceil(b.aabbX1)));
    int y1 = std::min(GH - 1, static_cast<int>(std::ceil(b.aabbY1)));
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        if (!occupiesWorldCellCenter(b, x, y)) continue;
        if (fluid.isPaintedSolid(x, y + 1)) return true;
        if (FluidEngine::inside(x, y + 1)) {
            int below = occupant[static_cast<size_t>(FluidEngine::ci(x, y + 1))];
            if (below >= 0 && below < static_cast<int>(bodies.size()) && bodies[static_cast<size_t>(below)].immobile())
                return true;
        }
    }
    return false;
}

void RigidBodyEngine::markSupport(std::vector<RigidContact> const &contacts) {
    for (RigidBody &b : bodies) { b.supported = false; b.maxPenetration = 0.0f; }
    for (RigidContact const &c : contacts) {
        auto apply = [&](int index, float ny) {
            if (index < 0) return;
            RigidBody &b = bodies[static_cast<size_t>(index)];
            b.maxPenetration = std::max(b.maxPenetration, c.penetration);
            if (ny < -0.35f) b.supported = true;
        };
        apply(c.bodyA, c.ny);
        apply(c.bodyB, c.bodyB >= 0 ? -c.ny : 0.0f);
    }
}

bool RigidBodyEngine::sweptThroughStaticSolid(FluidEngine const &fluid, RigidBody const &b, int index,
    float x0, float y0, float theta0) const {
    float c0 = std::cos(theta0), s0 = std::sin(theta0);
    float c1 = std::cos(b.theta), s1 = std::sin(b.theta);
    auto worldOf = [&](float x, float y, float c, float s, float lx, float ly, float &wx, float &wy) {
        float rx = lx - b.comLocalX, ry = ly - b.comLocalY;
        wx = x + c * rx - s * ry;
        wy = y + s * rx + c * ry;
    };
    for (int idx : b.occupiedLocal) {
        int lx = idx % b.maskW, ly = idx / b.maskW;
        float px0, py0, px1, py1;
        worldOf(x0, y0, c0, s0, static_cast<float>(lx) + 0.5f, static_cast<float>(ly) + 0.5f, px0, py0);
        worldOf(b.x, b.y, c1, s1, static_cast<float>(lx) + 0.5f, static_cast<float>(ly) + 0.5f, px1, py1);
        int gx0 = static_cast<int>(std::floor(px0)), gy0 = static_cast<int>(std::floor(py0));
        int gx1 = static_cast<int>(std::floor(px1)), gy1 = static_cast<int>(std::floor(py1));
        int steps = std::max(1, std::max(std::abs(gx1 - gx0), std::abs(gy1 - gy0)));
        for (int s = 1; s < steps; ++s) {
            float t = static_cast<float>(s) / static_cast<float>(steps);
            int gx = static_cast<int>(std::floor(px0 + (px1 - px0) * t));
            int gy = static_cast<int>(std::floor(py0 + (py1 - py0) * t));
            if ((gx == gx0 && gy == gy0) || (gx == gx1 && gy == gy1)) continue;
            if (fluid.isPaintedSolid(gx, gy)) return true;
            if (FluidEngine::inside(gx, gy)) {
                int occ = occupant[static_cast<size_t>(FluidEngine::ci(gx, gy))];
                if (occ >= 0 && occ != index && occ < static_cast<int>(bodies.size())
                    && bodies[static_cast<size_t>(occ)].immobile())
                    return true;
            }
        }
    }
    return false;
}

void RigidBodyEngine::collectStaticContacts(FluidEngine const &fluid, RigidBody const &b, int index, std::vector<RigidContact> &out) const {
    int x0 = std::max(0, static_cast<int>(std::floor(b.aabbX0)));
    int y0 = std::max(0, static_cast<int>(std::floor(b.aabbY0)));
    int x1 = std::min(GW - 1, static_cast<int>(std::ceil(b.aabbX1)));
    int y1 = std::min(GH - 1, static_cast<int>(std::ceil(b.aabbY1)));
    constexpr int ndx[4] = {-1, 1, 0, 0};
    constexpr int ndy[4] = {0, 0, -1, 1};
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        if (!occupiesWorldCellCenter(b, x, y)) continue;
        if (fluid.isPaintedSolid(x, y)) {
            float faceNx = 0.0f, faceNy = 0.0f;
            for (int n = 0; n < 4; ++n) {
                if (!fluid.isPaintedSolid(x + ndx[n], y + ndy[n])) {
                    faceNx += static_cast<float>(ndx[n]);
                    faceNy += static_cast<float>(ndy[n]);
                }
            }
            if (faceNx == 0.0f && faceNy == 0.0f) {
                faceNx = b.x - (x + 0.5f);
                faceNy = b.y - (y + 0.5f);
            }
            float len = std::sqrt(faceNx * faceNx + faceNy * faceNy);
            if (len < 1e-5f) { faceNx = 0.0f; faceNy = -1.0f; len = 1.0f; }
            faceNx /= len; faceNy /= len;
            float nx = faceNx, ny = faceNy;
            estimateTerrainNormal(fluid, x, y, faceNx, faceNy, nx, ny);
            float faceX = x + 0.5f + 0.5f * nx;
            float faceY = y + 0.5f + 0.5f * ny;
            float pen = measureFacePenetration(b, faceX, faceY, nx, ny);
            emitContact(out, index, -1, faceX, faceY, nx, ny, pen);
            continue;
        }
        for (int n = 0; n < 4; ++n) {
            int sx = x + ndx[n], sy = y + ndy[n];
            bool wall = fluid.isPaintedSolid(sx, sy)
                || (fluid.config.walledBorders && !FluidEngine::inside(sx, sy));
            if (!wall) continue;
            float faceNx = static_cast<float>(x - sx);
            float faceNy = static_cast<float>(y - sy);
            float nx = faceNx, ny = faceNy;
            estimateTerrainNormal(fluid, sx, sy, faceNx, faceNy, nx, ny);
            float faceX = 0.5f * static_cast<float>(x + sx) + 0.5f;
            float faceY = 0.5f * static_cast<float>(y + sy) + 0.5f;
            float pen = measureFacePenetration(b, faceX, faceY, nx, ny);
            emitContact(out, index, -1, faceX, faceY, nx, ny, pen);
        }
    }
}

void RigidBodyEngine::collectBodyContacts(std::vector<RigidContact> &out) {
    std::vector<int> pairCount(static_cast<size_t>(std::max(1, static_cast<int>(bodies.size()) * static_cast<int>(bodies.size()))), 0);
    auto pairKey = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        return a * static_cast<int>(bodies.size()) + b;
    };
    auto tryEmit = [&](int a, int b, float x, float y, float nx, float ny, float pen) {
        if (a == b || a < 0 || b < 0) return;
        int key = pairKey(a, b);
        if (pairCount[static_cast<size_t>(key)] >= 12) return;
        ++pairCount[static_cast<size_t>(key)];
        emitContact(out, a, b, x, y, nx, ny, pen);
    };

    for (OccupancyConflict const &ov : occupancyConflicts) {
        RigidBody const &ba = bodies[static_cast<size_t>(ov.bodyA)];
        RigidBody const &bb = bodies[static_cast<size_t>(ov.bodyB)];
        float nx = ba.x - bb.x, ny = ba.y - bb.y;
        if (std::abs(nx) + std::abs(ny) < 1e-5f) { nx = 0.0f; ny = -1.0f; }
        tryEmit(ov.bodyA, ov.bodyB, ov.x + 0.5f, ov.y + 0.5f, nx, ny, 0.42f);
    }

    constexpr int ndx[2] = {1, 0};
    constexpr int ndy[2] = {0, 1};
    for (size_t bi = 0; bi < bodies.size(); ++bi) {
        RigidBody const &body = bodies[bi];
        int x0 = std::max(0, static_cast<int>(std::floor(body.aabbX0)));
        int y0 = std::max(0, static_cast<int>(std::floor(body.aabbY0)));
        int x1 = std::min(GW - 1, static_cast<int>(std::ceil(body.aabbX1)));
        int y1 = std::min(GH - 1, static_cast<int>(std::ceil(body.aabbY1)));
        if (x1 < x0 || y1 < y0) continue;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            int i = FluidEngine::ci(x, y);
            int a = occupant[static_cast<size_t>(i)];
            if (a != static_cast<int>(bi)) continue;
            for (int n = 0; n < 2; ++n) {
                int nx = x + ndx[n], ny = y + ndy[n];
                if (!FluidEngine::inside(nx, ny)) continue;
                int b = occupant[static_cast<size_t>(FluidEngine::ci(nx, ny))];
                if (b < 0 || b == a) continue;
                float nnx = static_cast<float>(x - nx);
                float nny = static_cast<float>(y - ny);
                float faceX = x + 0.5f + 0.5f * ndx[n];
                float faceY = y + 0.5f + 0.5f * ndy[n];
                int mobile = bodies[static_cast<size_t>(a)].immobile() ? b : a;
                float nnLen = std::sqrt(nnx * nnx + nny * nny);
                float mnx = nnx, mny = nny;
                if (nnLen > 1e-5f) { mnx /= nnLen; mny /= nnLen; }
                if (mobile == b) { mnx = -mnx; mny = -mny; }
                float pen = measureFacePenetration(bodies[static_cast<size_t>(mobile)], faceX, faceY, mnx, mny);
                tryEmit(a, b, faceX, faceY, nnx, nny, std::max(0.28f, pen));
            }
        }
    }
}

void RigidBodyEngine::reduceStaticManifold(std::vector<RigidContact> &contacts) const {
    if (contacts.empty() || bodies.empty()) return;
    std::vector<RigidContact> bodyBody;
    std::vector<std::vector<int>> perBody(bodies.size());
    bodyBody.reserve(contacts.size());
    for (int i = 0; i < static_cast<int>(contacts.size()); ++i) {
        RigidContact const &c = contacts[static_cast<size_t>(i)];
        if (c.bodyB >= 0) {
            bodyBody.push_back(c);
            continue;
        }
        if (c.bodyA < 0 || c.bodyA >= static_cast<int>(bodies.size())) continue;
        perBody[static_cast<size_t>(c.bodyA)].push_back(i);
    }

    std::vector<RigidContact> kept;
    kept.reserve(contacts.size());
    for (size_t bi = 0; bi < perBody.size(); ++bi) {
        std::vector<int> const &idxs = perBody[bi];
        if (idxs.empty()) continue;
        std::vector<char> used(idxs.size(), 0);
        struct Group {
            std::vector<int> members;
            float maxPen = 0.0f;
        };
        std::vector<Group> groups;
        for (size_t a = 0; a < idxs.size(); ++a) {
            if (used[a]) continue;
            RigidContact const &ca = contacts[static_cast<size_t>(idxs[a])];
            Group g;
            g.members.push_back(idxs[a]);
            g.maxPen = ca.penetration;
            used[a] = 1;
            for (size_t b = a + 1; b < idxs.size(); ++b) {
                if (used[b]) continue;
                RigidContact const &cb = contacts[static_cast<size_t>(idxs[b])];
                if (ca.nx * cb.nx + ca.ny * cb.ny < kManifoldDot) continue;
                g.members.push_back(idxs[b]);
                g.maxPen = std::max(g.maxPen, cb.penetration);
                used[b] = 1;
            }
            groups.push_back(std::move(g));
        }
        if (static_cast<int>(groups.size()) > kMaxManifoldGroups) {
            std::sort(groups.begin(), groups.end(), [](Group const &l, Group const &r) {
                return l.maxPen > r.maxPen;
            });
            groups.resize(static_cast<size_t>(kMaxManifoldGroups));
        }
        for (Group const &g : groups) {
            if (g.members.size() <= static_cast<size_t>(kMaxManifoldPoints)) {
                for (int idx : g.members) kept.push_back(contacts[static_cast<size_t>(idx)]);
                continue;
            }
            RigidContact const &seed = contacts[static_cast<size_t>(g.members[0])];
            float tx = -seed.ny, ty = seed.nx;
            float tlen = std::sqrt(tx * tx + ty * ty);
            if (tlen > 1e-5f) { tx /= tlen; ty /= tlen; }
            int iMin = g.members[0], iMax = g.members[0];
            float pMin = 1e9f, pMax = -1e9f;
            for (int idx : g.members) {
                RigidContact const &c = contacts[static_cast<size_t>(idx)];
                float p = c.x * tx + c.y * ty;
                if (p < pMin) { pMin = p; iMin = idx; }
                if (p > pMax) { pMax = p; iMax = idx; }
            }
            kept.push_back(contacts[static_cast<size_t>(iMin)]);
            if (iMax != iMin) kept.push_back(contacts[static_cast<size_t>(iMax)]);
        }
    }
    kept.insert(kept.end(), bodyBody.begin(), bodyBody.end());
    contacts.swap(kept);
}

void RigidBodyEngine::solveVelocityContacts(std::vector<RigidContact> &contacts) {
    for (RigidContact &c : contacts) {
        RigidBody &a = bodies[static_cast<size_t>(c.bodyA)];
        RigidBody *b = (c.bodyB >= 0) ? &bodies[static_cast<size_t>(c.bodyB)] : nullptr;
        if (b) {
            float speedA = std::abs(a.vx) + std::abs(a.vy) + std::abs(a.omega);
            float speedB = std::abs(b->vx) + std::abs(b->vy) + std::abs(b->omega);
            if (a.dormant && !b->immobile() && speedB > 0.80f) wakeDormant(a);
            if (b->dormant && !a.immobile() && speedA > 0.80f) wakeDormant(*b);
        }
        if (a.sleeping && !a.immobile() && (!b || (b->sleeping && !b->immobile()))) continue;

        float vax = a.vx - a.omega * c.rAy;
        float vay = a.vy + a.omega * c.rAx;
        float vbx = 0.0f, vby = 0.0f;
        if (b) { vbx = b->vx - b->omega * c.rBy; vby = b->vy + b->omega * c.rBx; }
        float rvx = vax - vbx, rvy = vay - vby;
        float vn = rvx * c.nx + rvy * c.ny;
        if (vn < -kImpactSpeedMin) c.impactSpeed = std::max(c.impactSpeed, -vn);
        if (vn > 0.0f) continue;

        float rAn = c.rAx * c.ny - c.rAy * c.nx;
        float rBn = b ? (c.rBx * c.ny - c.rBy * c.nx) : 0.0f;
        float invMa, invIa, invMb = 0.0f, invIb = 0.0f;
        invTerms(a, invMa, invIa);
        if (b) invTerms(*b, invMb, invIb);
        float keff = invMa + invIa * rAn * rAn + (b ? (invMb + invIb * rBn * rBn) : 0.0f);
        if (keff < 1e-8f) continue;

        float eA = bodyContactRestitution(a);
        float eB = b ? bodyContactRestitution(*b) : 0.08f;
        float eMat = 0.5f * (eA + eB);
        float e = (vn < -kRestitutionThreshold) ? std::max(kRestitution, eMat) : 0.0f;
        float j = -(1.0f + e) * vn / keff;
        if (j < 0.0f) j = 0.0f;
        float jx = j * c.nx, jy = j * c.ny;
        a.vx += jx * invMa; a.vy += jy * invMa;
        a.omega += (c.rAx * jy - c.rAy * jx) * invIa;
        if (b) {
            b->vx -= jx * invMb; b->vy -= jy * invMb;
            b->omega -= (c.rBx * jy - c.rBy * jx) * invIb;
        }
        c.jn += j;
        a.debugJn += j;
        if (b) b->debugJn += j;

        vax = a.vx - a.omega * c.rAy;
        vay = a.vy + a.omega * c.rAx;
        if (b) { vbx = b->vx - b->omega * c.rBy; vby = b->vy + b->omega * c.rBx; }
        else { vbx = 0.0f; vby = 0.0f; }
        rvx = vax - vbx; rvy = vay - vby;
        float tx = -c.ny, ty = c.nx;
        float vt = rvx * tx + rvy * ty;
        float rAt = c.rAx * ty - c.rAy * tx;
        float rBt = b ? (c.rBx * ty - c.rBy * tx) : 0.0f;
        float kefft = invMa + invIa * rAt * rAt + (b ? (invMb + invIb * rBt * rBt) : 0.0f);
        if (kefft > 1e-8f) {
            float jt = -vt / kefft;
            float muA = bodyContactFriction(a);
            float muB = b ? bodyContactFriction(*b) : kFriction;
            float mu = 0.5f * (muA + muB);
            float maxF = mu * std::abs(j);
            jt = clampf(jt, -maxF, maxF);
            float fx = jt * tx, fy = jt * ty;
            a.vx += fx * invMa; a.vy += fy * invMa;
            a.omega += (c.rAx * fy - c.rAy * fx) * invIa;
            if (b) {
                b->vx -= fx * invMb; b->vy -= fy * invMb;
                b->omega -= (c.rBx * fy - c.rBy * fx) * invIb;
            }
            c.jt += jt;
            a.debugJt += jt;
            if (b) b->debugJt += jt;
        }
        if (b && std::abs(vn) >= kSleepLin) {
            if (a.sleeping) { a.sleeping = false; a.quietTicks = 0; }
            if (b->sleeping) { b->sleeping = false; b->quietTicks = 0; }
        }
    }
}

void RigidBodyEngine::solvePositionalContacts(std::vector<RigidContact> const &contacts) {
    struct Group {
        int bodyA = -1, bodyB = -1;
        float nx = 0.0f, ny = 0.0f, pen = 0.0f;
        int count = 0;
    };
    std::vector<Group> groups;
    groups.reserve(16);
    for (RigidContact const &c : contacts) {
        Group *found = nullptr;
        for (Group &g : groups) {
            if (g.bodyA != c.bodyA || g.bodyB != c.bodyB) continue;
            if (g.count == 0) continue;
            float inv = 1.0f / static_cast<float>(g.count);
            float gnx = g.nx * inv, gny = g.ny * inv;
            float glen = std::sqrt(gnx * gnx + gny * gny);
            if (glen > 1e-5f) { gnx /= glen; gny /= glen; }
            if (gnx * c.nx + gny * c.ny < kManifoldDot) continue;
            found = &g;
            break;
        }
        if (!found) {
            groups.push_back({c.bodyA, c.bodyB, c.nx, c.ny, c.penetration, 1});
        } else {
            found->nx += c.nx; found->ny += c.ny;
            found->pen = std::max(found->pen, c.penetration);
            ++found->count;
        }
    }
    for (Group const &g : groups) {
        RigidBody &a = bodies[static_cast<size_t>(g.bodyA)];
        RigidBody *b = (g.bodyB >= 0) ? &bodies[static_cast<size_t>(g.bodyB)] : nullptr;
        if (a.sleeping && !a.immobile() && (!b || (b->sleeping && !b->immobile()))) continue;
        float nx = g.nx, ny = g.ny;
        float len = std::sqrt(nx * nx + ny * ny);
        if (len < 1e-5f) continue;
        nx /= len; ny /= len;
        float depth = g.pen - kSlop;
        if (depth <= 0.0f) continue;
        float corr = depth * kBaumgarte;
        float invMa, invIa, invMb = 0.0f, invIb = 0.0f;
        invTerms(a, invMa, invIa);
        if (b) invTerms(*b, invMb, invIb);
        float keff = invMa + invMb;
        if (keff < 1e-8f) continue;
        float lambda = corr / keff;
        float dx = lambda * nx, dy = lambda * ny;
        a.x += dx * invMa; a.y += dy * invMa;
        a.debugPosCorrX += dx * invMa;
        a.debugPosCorrY += dy * invMa;
        if (b) {
            b->x -= dx * invMb; b->y -= dy * invMb;
            b->debugPosCorrX -= dx * invMb;
            b->debugPosCorrY -= dy * invMb;
        }
        (void)invIa; (void)invIb;
    }
}

void RigidBodyEngine::rasterizeBodyOccupancy(RigidBody const &b, int index, bool writeFluid, bool conservative, FluidEngine *fluid) {
    int x0 = std::max(0, static_cast<int>(std::floor(b.aabbX0)));
    int y0 = std::max(0, static_cast<int>(std::floor(b.aabbY0)));
    int x1 = std::min(GW - 1, static_cast<int>(std::ceil(b.aabbX1)));
    int y1 = std::min(GH - 1, static_cast<int>(std::ceil(b.aabbY1)));
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        bool hit = conservative ? occupiesWorldCell(b, x, y) : occupiesWorldCellCenter(b, x, y);
        if (!hit) continue;
        int i = FluidEngine::ci(x, y);
        int prev = occupant[static_cast<size_t>(i)];
        if (!writeFluid && prev >= 0 && prev != index)
            occupancyConflicts.push_back({prev, index, x, y});
        occupant[static_cast<size_t>(i)] = index;
        float lx, ly;
        worldToLocal(b, x + 0.5f, y + 0.5f, lx, ly);
        occupantMat[static_cast<size_t>(i)] = maskMaterial(b, static_cast<int>(std::floor(lx)), static_cast<int>(std::floor(ly)));
        int ix = static_cast<int>(std::floor(lx));
        int iy = static_cast<int>(std::floor(ly));
        if (ix >= 0 && iy >= 0 && ix < b.maskW && iy < b.maskH) {
            size_t li = static_cast<size_t>(iy * b.maskW + ix);
            if (li < b.materialDamage.size()) occupantDamage[static_cast<size_t>(i)] = b.materialDamage[li];
            if (li < b.moisture.size()) occupantMoisture[static_cast<size_t>(i)] = b.moisture[li];
            occupantCrack[static_cast<size_t>(i)] = getLocalCrackFactor(b, ix, iy);
        }
        if (occupantMat[static_cast<size_t>(i)] == MATERIAL_EMPTY) {
            occupantMat[static_cast<size_t>(i)] = MATERIAL_WOOD;
            for (int iOcc : b.occupiedLocal) {
                occupantMat[static_cast<size_t>(i)] = b.mask[static_cast<size_t>(iOcc)];
                break;
            }
        }
        if (writeFluid && fluid) {
            float rx = (x + 0.5f) - b.x, ry = (y + 0.5f) - b.y;
            float svx = b.vx - b.omega * ry;
            float svy = b.vy + b.omega * rx;
            fluid->setDynamicOccupancy(x, y, svx, svy);
        }
    }
}

void RigidBodyEngine::syncOccupancy(FluidEngine &fluid) {
    std::fill(occupant.begin(), occupant.end(), -1);
    std::fill(occupantMat.begin(), occupantMat.end(), MATERIAL_EMPTY);
    std::fill(occupantDamage.begin(), occupantDamage.end(), 0.0f);
    std::fill(occupantMoisture.begin(), occupantMoisture.end(), 0.0f);
    std::fill(occupantCrack.begin(), occupantCrack.end(), 0.0f);
    fluid.clearDynamicOccupancy();
    auto wakeBody = [&](size_t i) {
        int wx0 = static_cast<int>(std::floor(std::min(bodies[i].prevAabbX0, bodies[i].aabbX0))) - 1;
        int wy0 = static_cast<int>(std::floor(std::min(bodies[i].prevAabbY0, bodies[i].aabbY0))) - 1;
        int wx1 = static_cast<int>(std::ceil(std::max(bodies[i].prevAabbX1, bodies[i].aabbX1))) + 1;
        int wy1 = static_cast<int>(std::ceil(std::max(bodies[i].prevAabbY1, bodies[i].aabbY1))) + 1;
        fluid.wakeRegion(wx0, wy0, wx1, wy1);
    };
    for (size_t i = 0; i < bodies.size(); ++i) {
        updateAabb(bodies[i]);
        wakeBody(i);
        if (!bodies[i].immobile())
            rasterizeBodyOccupancy(bodies[i], static_cast<int>(i), true, true, &fluid);
    }
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].immobile())
            rasterizeBodyOccupancy(bodies[i], static_cast<int>(i), true, true, &fluid);
    }
    fluid.displaceFluidFromDynamicSolids();
}

void RigidBodyEngine::gatherFluidForces(FluidEngine const &fluid) {
    // Future pass should share substeps: rigid integrate → rasterize/displace/moving BC
    // → fluid project → forces inside the same substep. This tick still uses the
    // split worldTick (rigid then fluid then gather) with clamped hydrostatic-dominant forces.
    constexpr int ndx[4] = {-1, 1, 0, 0};
    constexpr int ndy[4] = {0, 0, -1, 1};
    float gravity = fluid.gridGravity();
    for (size_t bi = 0; bi < bodies.size(); ++bi) {
        RigidBody &b = bodies[bi];
        b.fx = b.fy = b.torque = 0.0f;
        if (b.anchored) continue;
        int x0 = std::max(1, static_cast<int>(std::floor(b.aabbX0)) - 1);
        int y0 = std::max(1, static_cast<int>(std::floor(b.aabbY0)) - 1);
        int x1 = std::min(GW - 2, static_cast<int>(std::ceil(b.aabbX1)) + 1);
        int y1 = std::min(GH - 2, static_cast<int>(std::ceil(b.aabbY1)) + 1);
        float vxSamp = 0.0f, vySamp = 0.0f, wet = 0.0f;
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
            if (occupant[static_cast<size_t>(FluidEngine::ci(x, y))] != static_cast<int>(bi)) continue;
            for (int n = 0; n < 4; ++n) {
                int nx = x + ndx[n], ny = y + ndy[n];
                if (!FluidEngine::inside(nx, ny) || fluid.isSolid(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                float fill = fluid.fill[static_cast<size_t>(ni)];
                if (fill < MIN_ACTIVE_FILL) continue;
                float persist = std::clamp(fluid.previousFill[static_cast<size_t>(ni)] / 0.25f, 0.12f, 1.0f);
                float p = fluid.pressure[static_cast<size_t>(ni)] * std::min(fill, 1.0f) * persist;
                float fx = -p * static_cast<float>(ndx[n]);
                float fy = -p * static_cast<float>(ndy[n]);
                float rx = (x + 0.5f) - b.x, ry = (y + 0.5f) - b.y;
                b.fx += fx; b.fy += fy;
                b.torque += rx * fy - ry * fx;
                vxSamp += fluid.cellU(nx, ny);
                vySamp += fluid.cellV(nx, ny);
                wet += 1.0f;
            }
        }
        if (wet > 0.0f) {
            vxSamp /= wet; vySamp /= wet;
            float drag = 1.15f * b.mass;
            b.fx -= (b.vx - vxSamp) * drag;
            b.fy -= (b.vy - vySamp) * drag;
        }
        float accel = std::sqrt(b.fx * b.fx + b.fy * b.fy) * b.invMass;
        float aMax = 3.0f * gravity;
        if (accel > aMax && accel > 1e-5f) {
            float s = aMax / accel;
            b.fx *= s; b.fy *= s; b.torque *= s;
            accel = aMax;
        }
        float ang = std::abs(b.torque) * b.invInertia;
        float angMax = 12.0f;
        if (ang > angMax && ang > 1e-5f) b.torque *= angMax / ang;
        if (b.dormant) {
            if (accel > 2.5f * gravity) wakeDormant(b);
            else { b.fx = b.fy = b.torque = 0.0f; }
            continue;
        }
        if (b.sleeping && accel > 2.0f * gravity) { b.sleeping = false; b.quietTicks = 0; }
    }
}

void RigidBodyEngine::step(FluidEngine &fluid, float dt) {
    auto start = FluidEngine::Clock::now();
    float gravity = fluid.gridGravity() * gravityScale;
    applyGrabForces(gravity);
    float maxMotion = 0.0f;
    for (RigidBody &b : bodies) {
        if (b.immobile()) continue;
        if (b.sleeping) continue;
        float radius = 0.5f * std::sqrt(static_cast<float>(b.maskW * b.maskW + b.maskH * b.maskH));
        maxMotion = std::max(maxMotion, std::abs(b.vx) + std::abs(b.vy) + std::abs(b.omega) * radius);
    }
    int sub = std::clamp(static_cast<int>(std::ceil(maxMotion * dt / 0.28f)), 1, 16);
    float h = dt / static_cast<float>(sub);
    for (RigidBody &b : bodies) {
        if (b.immobile()) {
            b.vx = b.vy = b.omega = 0.0f;
            b.fx = b.fy = b.torque = 0.0f;
            b.sleeping = true;
            continue;
        }
        if (b.sleeping) {
            if (!bodySupported(fluid, b)) { b.sleeping = false; b.quietTicks = 0; }
            else continue;
        }
        b.debugJn = b.debugJt = 0.0f;
        b.debugPosCorrX = b.debugPosCorrY = 0.0f;
        b.vy += gravity * dt;
        b.vx += b.fx * b.invMass * dt;
        b.vy += b.fy * b.invMass * dt;
        b.omega += b.torque * b.invInertia * dt;
        b.fx = b.fy = b.torque = 0.0f;
        b.vx = clampf(b.vx, -kMaxLin, kMaxLin);
        b.vy = clampf(b.vy, -kMaxLin, kMaxLin);
        b.omega = clampf(b.omega, -kMaxAng, kMaxAng);
    }
    for (int s = 0; s < sub; ++s) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            RigidBody &b = bodies[i];
            if (b.sleeping || b.immobile()) continue;
            float ox = b.x, oy = b.y, oth = b.theta;
            integrate(b, h);
            if (sweptThroughStaticSolid(fluid, b, static_cast<int>(i), ox, oy, oth)) {
                b.x = ox; b.y = oy; b.theta = oth;
                b.vx *= 0.35f; b.vy *= 0.35f; b.omega *= 0.35f;
            }
            updateAabb(b);
        }
        std::vector<RigidContact> contacts;
        contacts.reserve(64);
        for (size_t i = 0; i < bodies.size(); ++i) collectStaticContacts(fluid, bodies[i], static_cast<int>(i), contacts);
        occupancyConflicts.clear();
        std::fill(occupant.begin(), occupant.end(), -1);
        std::fill(occupantDamage.begin(), occupantDamage.end(), 0.0f);
        std::fill(occupantMoisture.begin(), occupantMoisture.end(), 0.0f);
        std::fill(occupantCrack.begin(), occupantCrack.end(), 0.0f);
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (bodies[i].immobile()) continue;
            rasterizeBodyOccupancy(bodies[i], static_cast<int>(i), false, false, nullptr);
        }
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (!bodies[i].immobile()) continue;
            rasterizeBodyOccupancy(bodies[i], static_cast<int>(i), false, false, nullptr);
        }
        collectBodyContacts(contacts);
        reduceStaticManifold(contacts);
        markSupport(contacts);
        for (int it = 0; it < 8; ++it) solveVelocityContacts(contacts);
        applyContactDamage(contacts);
        solvePositionalContacts(contacts);
        if (s == sub - 1) lastContacts = contacts;
    }
    splitDirtyBodies();
    if (grab.active && indexOfId(grab.bodyId) < 0) grab = GrabState{};
    for (RigidBody &b : bodies) {
        if (b.immobile()) {
            b.vx = b.vy = b.omega = 0.0f;
            updateAabb(b);
            continue;
        }
        b.vx *= 0.997f; b.vy *= 0.997f; b.omega *= 0.995f;
        updateAabb(b);
        if (b.sleeping) continue;
        if (grab.active && b.id == grab.bodyId) { b.quietTicks = 0; continue; }
        bool still = std::abs(b.vx) < kSleepLin && std::abs(b.vy) < kSleepLin && std::abs(b.omega) < kSleepAng;
        if (still && b.supported) {
            b.quietTicks++;
            if (b.quietTicks >= kSleepTicks) {
                b.sleeping = true; b.vx = b.vy = b.omega = 0.0f;
            }
        } else b.quietTicks = 0;
    }
    cullBodiesLeftIntoVoid();
    if (grab.active && indexOfId(grab.bodyId) < 0) grab = GrabState{};
    syncOccupancy(fluid);
    processMoisture(fluid, dt);
    lastStepMs = FluidEngine::elapsedMs(start);
}

void RigidBodyEngine::cullBodiesLeftIntoVoid() {
    size_t write = 0;
    constexpr float margin = 2.0f;
    for (size_t i = 0; i < bodies.size(); ++i) {
        RigidBody const &b = bodies[i];
        bool gone = b.aabbX1 < -margin || b.aabbX0 > static_cast<float>(GW) + margin
            || b.aabbY1 < -margin || b.aabbY0 > static_cast<float>(GH) + margin;
        if (gone) continue;
        if (write != i) bodies[write] = std::move(bodies[i]);
        ++write;
    }
    if (write != bodies.size()) {
        bodies.resize(write);
        lastContacts.clear();
    }
}

bool RigidBodyEngine::worldCellOccupied(int x, int y) const {
    if (!FluidEngine::inside(x, y)) return false;
    return occupant[static_cast<size_t>(FluidEngine::ci(x, y))] >= 0;
}

MaterialId RigidBodyEngine::worldCellMaterial(int x, int y) const {
    if (!FluidEngine::inside(x, y)) return MATERIAL_EMPTY;
    return occupantMat[static_cast<size_t>(FluidEngine::ci(x, y))];
}

int RigidBodyEngine::indexOfId(uint32_t id) const {
    for (size_t i = 0; i < bodies.size(); ++i)
        if (bodies[i].id == id) return static_cast<int>(i);
    return -1;
}

float RigidBodyEngine::bodyContactFriction(RigidBody const &b) const {
    return std::max(0.05f, b.cachedFriction);
}

float RigidBodyEngine::bodyContactRestitution(RigidBody const &b) const {
    return clampf(b.cachedRestitution, 0.0f, 0.6f);
}

bool RigidBodyEngine::beginGrab(float wx, float wy) {
    int gx = static_cast<int>(std::floor(wx));
    int gy = static_cast<int>(std::floor(wy));
    int index = bodyAtCell(gx, gy);
    if (index < 0) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            float lx, ly;
            worldToLocal(bodies[i], wx, wy, lx, ly);
            int ix = static_cast<int>(std::floor(lx));
            int iy = static_cast<int>(std::floor(ly));
            if (maskOccupied(bodies[i], ix, iy)) { index = static_cast<int>(i); break; }
        }
    }
    if (index < 0 || index >= static_cast<int>(bodies.size())) return false;
    RigidBody &b = bodies[static_cast<size_t>(index)];
    if (b.anchored) return false;
    float lx, ly;
    worldToLocal(b, wx, wy, lx, ly);
    grab.active = true;
    grab.bodyId = b.id;
    grab.localX = lx;
    grab.localY = ly;
    grab.targetX = wx;
    grab.targetY = wy;
    grab.worldX = wx;
    grab.worldY = wy;
    grab.lastFx = grab.lastFy = 0.0f;
    wakeDormant(b);
    b.sleeping = false;
    b.quietTicks = 0;
    return true;
}

void RigidBodyEngine::updateGrabTarget(float wx, float wy, bool strong) {
    if (!grab.active) return;
    grab.targetX = wx;
    grab.targetY = wy;
    grab.strong = strong;
}

void RigidBodyEngine::endGrab() {
    grab = GrabState{};
}

void RigidBodyEngine::applyGrabForces(float gravity) {
    if (!grab.active) return;
    int index = indexOfId(grab.bodyId);
    if (index < 0) { grab = GrabState{}; return; }
    RigidBody &b = bodies[static_cast<size_t>(index)];
    if (b.anchored) { grab = GrabState{}; return; }
    wakeDormant(b);
    b.sleeping = false;
    b.quietTicks = 0;
    float wx, wy;
    localToWorld(b, grab.localX, grab.localY, wx, wy);
    grab.worldX = wx;
    grab.worldY = wy;
    float rx = wx - b.x, ry = wy - b.y;
    float pvx = b.vx - b.omega * ry;
    float pvy = b.vy + b.omega * rx;
    float ex = grab.targetX - wx;
    float ey = grab.targetY - wy;
    float stiff = (grab.strong ? 210.0f : 90.0f) * b.mass;
    float damp = (grab.strong ? 26.0f : 13.0f) * b.mass;
    float fx = ex * stiff - pvx * damp;
    float fy = ey * stiff - pvy * damp;
    float g = std::max(gravity, 8.0f);
    float fmax = (grab.strong ? 36.0f : 16.0f) * b.mass * g;
    float mag = std::sqrt(fx * fx + fy * fy);
    if (mag > fmax && mag > 1.0e-5f) {
        float s = fmax / mag;
        fx *= s; fy *= s;
    }
    b.fx += fx;
    b.fy += fy;
    b.torque += rx * fy - ry * fx;
    grab.lastFx = fx;
    grab.lastFy = fy;
}

void RigidBodyEngine::setFragmentVelocity(RigidBody &child, float px, float py, float pvx, float pvy, float pomega) {
    float rx = child.x - px;
    float ry = child.y - py;
    child.vx = pvx - pomega * ry;
    child.vy = pvy + pomega * rx;
    child.omega = pomega;
    if (child.anchored) child.vx = child.vy = child.omega = 0.0f;
}

RigidBody RigidBodyEngine::packFragment(RigidBody const &parent, std::vector<int> const &localCells, bool keepAnchored) {
    RigidBody child;
    child.theta = parent.theta;
    int minX = parent.maskW, minY = parent.maskH, maxX = -1, maxY = -1;
    for (int i : localCells) {
        int lx = i % parent.maskW, ly = i / parent.maskW;
        minX = std::min(minX, lx); maxX = std::max(maxX, lx);
        minY = std::min(minY, ly); maxY = std::max(maxY, ly);
    }
    child.maskW = maxX - minX + 1;
    child.maskH = maxY - minY + 1;
    size_t n = static_cast<size_t>(child.maskW * child.maskH);
    child.mask.assign(n, MATERIAL_EMPTY);
    child.materialDamage.assign(n, 0.0f);
    child.moisture.assign(n, 0.0f);
    child.heat.assign(n, 0.0f);
    child.bondsRight.assign(n, StructuralBond{});
    child.bondsDown.assign(n, StructuralBond{});
    std::vector<uint8_t> keep(static_cast<size_t>(parent.maskW * parent.maskH), 0);
    for (int i : localCells) keep[static_cast<size_t>(i)] = 1;
    for (int i : localCells) {
        int lx = i % parent.maskW, ly = i / parent.maskW;
        size_t dst = static_cast<size_t>((ly - minY) * child.maskW + (lx - minX));
        child.mask[dst] = parent.mask[static_cast<size_t>(i)];
        if (i < static_cast<int>(parent.materialDamage.size())) child.materialDamage[dst] = parent.materialDamage[static_cast<size_t>(i)];
        if (i < static_cast<int>(parent.moisture.size())) child.moisture[dst] = parent.moisture[static_cast<size_t>(i)];
        if (i < static_cast<int>(parent.heat.size())) child.heat[dst] = parent.heat[static_cast<size_t>(i)];
        int rx = lx + 1, ry = ly;
        if (rx < parent.maskW && keep[static_cast<size_t>(ry * parent.maskW + rx)])
            child.bondsRight[dst] = parent.bondsRight[static_cast<size_t>(i)];
        int dx = lx, dy = ly + 1;
        if (dy < parent.maskH && keep[static_cast<size_t>(dy * parent.maskW + dx)])
            child.bondsDown[dst] = parent.bondsDown[static_cast<size_t>(i)];
    }
    rebuildDerived(child);
    computeMassProperties(child);
    float wx, wy;
    localToWorld(parent, static_cast<float>(minX) + child.comLocalX, static_cast<float>(minY) + child.comLocalY, wx, wy);
    child.x = wx;
    child.y = wy;
    setFragmentVelocity(child, parent.x, parent.y, parent.vx, parent.vy, parent.omega);
    child.id = nextId++;
    child.anchored = keepAnchored;
    child.dormant = false;
    child.sleeping = child.anchored;
    child.structureDirty = false;
    if (child.anchored) child.vx = child.vy = child.omega = 0.0f;
    updateAabb(child);
    child.prevAabbX0 = child.aabbX0; child.prevAabbY0 = child.aabbY0;
    child.prevAabbX1 = child.aabbX1; child.prevAabbY1 = child.aabbY1;
    return child;
}

StructuralBond *RigidBodyEngine::bondPtr(RigidBody &b, int lx, int ly, int nx, int ny) {
    return const_cast<StructuralBond *>(static_cast<RigidBodyEngine const *>(this)->bondPtr(b, lx, ly, nx, ny));
}

StructuralBond const *RigidBodyEngine::bondPtr(RigidBody const &b, int lx, int ly, int nx, int ny) const {
    if (!maskOccupied(b, lx, ly) || !maskOccupied(b, nx, ny)) return nullptr;
    int dx = nx - lx, dy = ny - ly;
    if (std::abs(dx) + std::abs(dy) != 1) return nullptr;
    int x = lx, y = ly;
    bool down = false;
    if (dx == 1) { /* right from lx,ly */ }
    else if (dx == -1) { x = nx; y = ny; }
    else if (dy == 1) { down = true; }
    else { x = nx; y = ny; down = true; }
    size_t i = static_cast<size_t>(y * b.maskW + x);
    if (down) {
        if (i >= b.bondsDown.size()) return nullptr;
        return &b.bondsDown[i];
    }
    if (i >= b.bondsRight.size()) return nullptr;
    return &b.bondsRight[i];
}

bool RigidBodyEngine::bondConnects(RigidBody const &b, int lx, int ly, int nx, int ny) const {
    StructuralBond const *bond = bondPtr(b, lx, ly, nx, ny);
    return bond && !bond->broken;
}

float RigidBodyEngine::getLocalCrackFactor(RigidBody const &b, int lx, int ly) const {
    if (!maskOccupied(b, lx, ly)) return 0.0f;
    float f = 0.0f;
    for (int n = 0; n < 4; ++n) {
        StructuralBond const *bond = bondPtr(b, lx, ly, lx + kN4x[n], ly + kN4y[n]);
        if (!bond) continue;
        if (bond->broken) return 1.0f;
        f = std::max(f, bond->damage);
    }
    return f;
}

float RigidBodyEngine::getLocalCrackFactorWorld(int x, int y) const {
    if (!FluidEngine::inside(x, y)) return 0.0f;
    size_t i = static_cast<size_t>(FluidEngine::ci(x, y));
    if (i >= occupantCrack.size()) return 0.0f;
    return occupantCrack[i];
}

int RigidBodyEngine::brokenBondCountOf(RigidBody const &b) const {
    return b.brokenBondCount;
}

void RigidBodyEngine::recountStructure(RigidBody &b) {
    ensurePixelState(b);
    int broken = 0;
    float maxB = 0.0f;
    auto tally = [&](std::vector<StructuralBond> const &bonds, int axis) {
        for (int y = 0; y < b.maskH; ++y) for (int x = 0; x < b.maskW; ++x) {
            int nx = x + (axis == 0 ? 1 : 0);
            int ny = y + (axis == 1 ? 1 : 0);
            if (!maskOccupied(b, x, y) || !maskOccupied(b, nx, ny)) continue;
            StructuralBond const &bond = bonds[static_cast<size_t>(y * b.maskW + x)];
            maxB = std::max(maxB, bond.damage);
            if (bond.broken) ++broken;
        }
    };
    tally(b.bondsRight, 0);
    tally(b.bondsDown, 1);
    b.brokenBondCount = broken;
    b.maxBondDamage = maxB;
}

bool RigidBodyEngine::damageBond(RigidBody &b, int lx, int ly, int nx, int ny, float add, bool *newlyBroken) {
    if (newlyBroken) *newlyBroken = false;
    StructuralBond *bond = bondPtr(b, lx, ly, nx, ny);
    if (!bond || add <= 1.0e-6f) return false;
    if (bond->broken) {
        bond->damage = 1.0f;
        return false;
    }
    bond->damage = std::min(1.0f, bond->damage + add);
    if (bond->damage >= kFailDamage) {
        bond->broken = true;
        bond->damage = 1.0f;
        if (newlyBroken) *newlyBroken = true;
        return true;
    }
    return false;
}

bool RigidBodyEngine::breakBond(int bodyIndex, int lx, int ly, int nx, int ny) {
    if (bodyIndex < 0 || bodyIndex >= static_cast<int>(bodies.size())) return false;
    RigidBody &b = bodies[static_cast<size_t>(bodyIndex)];
    ensurePixelState(b);
    StructuralBond *bond = bondPtr(b, lx, ly, nx, ny);
    if (!bond) return false;
    bool was = bond->broken;
    bond->broken = true;
    bond->damage = 1.0f;
    if (!was) {
        b.structureDirty = true;
        recountStructure(b);
    }
    return true;
}

void RigidBodyEngine::propagateCracks(RigidBody &b, std::vector<int> const &seedBonds) {
    if (seedBonds.empty()) return;
    std::vector<int> queue = seedBonds;
    std::vector<uint8_t> seen(static_cast<size_t>(std::max(1, b.maskW * b.maskH * 2)), 0);
    auto mark = [&](int packed) {
        int id = packed & 0x0fffffff;
        int axis = packed >> 28;
        size_t key = static_cast<size_t>(axis * b.maskW * b.maskH + id);
        if (key >= seen.size()) return false;
        if (seen[key]) return false;
        seen[key] = 1;
        return true;
    };
    for (int s : queue) mark(s);
    int guard = 0;
    size_t head = 0;
    while (head < queue.size() && guard++ < 48) {
        int packed = queue[head++];
        int id = packed & 0x0fffffff;
        int axis = packed >> 28;
        int x = id % b.maskW, y = id / b.maskW;
        MaterialDefinition const &mat = materialDef(maskMaterial(b, x, y));
        auto consider = [&](int ax, int ay, int naxis) {
            if (!maskOccupied(b, ax, ay)) return;
            int ni = ay * b.maskW + ax;
            if (naxis == 0 && ax + 1 >= b.maskW) return;
            if (naxis == 1 && ay + 1 >= b.maskH) return;
            int nx = ax + (naxis == 0 ? 1 : 0);
            int ny = ay + (naxis == 1 ? 1 : 0);
            if (!maskOccupied(b, nx, ny)) return;
            MaterialDefinition const &nm = materialDef(maskMaterial(b, ax, ay));
            float wetN = pixelWetness(b, ni);
            float toughN = nm.fractureToughness * (1.0f + (nm.wetFractureToughnessMultiplier - 1.0f) * wetN);
            toughN = std::max(toughN, 0.08f);
            float add = clampf(0.08f + 0.95f * nm.brittleness - 0.55f * toughN, 0.0f, 0.95f);
            if (add < 0.04f) return;
            bool fresh = false;
            damageBond(b, ax, ay, nx, ny, add, &fresh);
            if (fresh) {
                int packedN = (naxis << 28) | ni;
                if (mark(packedN)) queue.push_back(packedN);
            }
        };
        if (axis == 0) {
            consider(x, y - 1, 0);
            consider(x, y + 1, 0);
            consider(x + 1, y, 0);
            consider(x - 1, y, 0);
        } else {
            consider(x - 1, y, 1);
            consider(x + 1, y, 1);
            consider(x, y + 1, 1);
            consider(x, y - 1, 1);
        }
        if (mat.brittleness > 0.55f) {
            consider(x, y, axis == 0 ? 1 : 0);
        }
    }
}

void RigidBodyEngine::applyStructuralDamage(int bodyIndex, StructuralDamageRequest const &request) {
    if (bodyIndex < 0 || bodyIndex >= static_cast<int>(bodies.size())) return;
    if (request.amount <= 1.0e-6f) return;
    RigidBody &b = bodies[static_cast<size_t>(bodyIndex)];
    ensurePixelState(b);
    int brokenBefore = b.brokenBondCount;
    float lx, ly;
    worldToLocal(b, request.wx, request.wy, lx, ly);
    float dirX = request.dirX, dirY = request.dirY;
    float dlen = std::sqrt(dirX * dirX + dirY * dirY);
    if (dlen > 1.0e-5f) { dirX /= dlen; dirY /= dlen; }
    else { dirX = 0.0f; dirY = 0.0f; }
    std::vector<int> newlyBroken;
    bool anyMat = false;
    float radius = std::max(0.75f, request.radius);
    for (int i : b.occupiedLocal) {
        int px = i % b.maskW, py = i / b.maskW;
        MaterialDefinition const &mat = materialDef(b.mask[static_cast<size_t>(i)]);
        float dx = (px + 0.5f) - lx;
        float dy = (py + 0.5f) - ly;
        float dist = std::sqrt(dx * dx + dy * dy);
        float rad = radius * clampf(1.05f + 2.8f * (1.0f - mat.brittleness), 0.85f, 3.8f);
        if (dist > rad + 0.35f) continue;
        float falloff = 1.0f - dist / (rad + 0.35f);
        falloff = falloff * falloff;
        float resist = effectiveStrength(b, i, request.type);
        float focus = 0.30f + 1.55f * mat.brittleness;
        float add = (request.amount * falloff * focus) / (resist * 1.35f + 0.12f);
        if (add < 1.0e-5f) continue;
        float &md = b.materialDamage[static_cast<size_t>(i)];
        md = std::min(1.0f, md + add * 0.55f);
        b.maxDamage = std::max(b.maxDamage, md);
        anyMat = true;
        for (int n = 0; n < 4; ++n) {
            int nx = px + kN4x[n], ny = py + kN4y[n];
            if (!maskOccupied(b, nx, ny)) continue;
            float bx = 0.5f * static_cast<float>(kN4x[n]);
            float by = 0.5f * static_cast<float>(kN4y[n]);
            float align = 1.0f;
            if (dlen > 1.0e-5f) {
                float along = std::abs(bx * dirX + by * dirY);
                align = 0.45f + 0.90f * along;
                if (request.type == StructuralDamageType::Tension) align = 0.25f + 1.20f * along;
                if (request.type == StructuralDamageType::Compression) align = 1.10f - 0.55f * along;
                if (request.type == StructuralDamageType::Shear) align = 0.70f + 0.50f * (1.0f - along);
            }
            bool fresh = false;
            damageBond(b, px, py, nx, ny, add * align * (0.85f + 0.80f * mat.brittleness), &fresh);
            if (fresh) {
                int axis = (kN4y[n] != 0) ? 1 : 0;
                int ox = std::min(px, nx), oy = std::min(py, ny);
                newlyBroken.push_back((axis << 28) | (oy * b.maskW + ox));
            }
        }
    }
    if (request.canPropagate && !newlyBroken.empty())
        propagateCracks(b, newlyBroken);
    recountStructure(b);
    if (b.brokenBondCount > brokenBefore) b.structureDirty = true;
    if ((anyMat || !newlyBroken.empty()) && !b.anchored) {
        b.sleeping = false;
        b.quietTicks = 0;
        b.dormant = false;
    }
}

void RigidBodyEngine::applyMaterialDamage(int bodyIndex, float wx, float wy, float amount, float radius, StructuralDamageType type) {
    (void)type;
    if (bodyIndex < 0 || bodyIndex >= static_cast<int>(bodies.size())) return;
    if (amount <= 1.0e-6f) return;
    RigidBody &b = bodies[static_cast<size_t>(bodyIndex)];
    ensurePixelState(b);
    float lx, ly;
    worldToLocal(b, wx, wy, lx, ly);
    float rad = std::max(0.75f, radius);
    bool any = false;
    for (int i : b.occupiedLocal) {
        int px = i % b.maskW, py = i / b.maskW;
        float dx = (px + 0.5f) - lx, dy = (py + 0.5f) - ly;
        float dist = std::sqrt(dx * dx + dy * dy);
        if (dist > rad + 0.35f) continue;
        float falloff = 1.0f - dist / (rad + 0.35f);
        falloff = falloff * falloff;
        float add = amount * falloff;
        if (add < 1.0e-5f) continue;
        float &md = b.materialDamage[static_cast<size_t>(i)];
        md = std::min(1.0f, md + add);
        b.maxDamage = std::max(b.maxDamage, md);
        any = true;
    }
    if (any && !b.anchored) {
        b.sleeping = false;
        b.quietTicks = 0;
        b.dormant = false;
    }
}

void RigidBodyEngine::applyBondDamage(int bodyIndex, float wx, float wy, float amount, float radius,
    float dirX, float dirY, StructuralDamageType type) {
    StructuralDamageRequest req;
    req.wx = wx; req.wy = wy; req.amount = amount; req.radius = radius;
    req.dirX = dirX; req.dirY = dirY; req.type = type; req.canPropagate = true;
    applyStructuralDamage(bodyIndex, req);
}

void RigidBodyEngine::applyLocalStructuralStress(int bodyIndex, float wx, float wy, float amount, float radius,
    float dirX, float dirY, StructuralDamageType type) {
    applyBondDamage(bodyIndex, wx, wy, amount, radius, dirX, dirY, type);
}

void RigidBodyEngine::modifyLocalStrength(int bodyIndex, float wx, float wy, float radius, float factor) {
    if (bodyIndex < 0 || bodyIndex >= static_cast<int>(bodies.size())) return;
    RigidBody &b = bodies[static_cast<size_t>(bodyIndex)];
    ensurePixelState(b);
    float lx, ly;
    worldToLocal(b, wx, wy, lx, ly);
    float add = clampf(1.0f - factor, 0.0f, 1.0f);
    if (add <= 1.0e-6f) return;
    for (int i : b.occupiedLocal) {
        int px = i % b.maskW, py = i / b.maskW;
        float dx = (px + 0.5f) - lx, dy = (py + 0.5f) - ly;
        float dist = std::sqrt(dx * dx + dy * dy);
        if (dist > radius + 0.35f) continue;
        float falloff = 1.0f - dist / (radius + 0.35f);
        b.materialDamage[static_cast<size_t>(i)] = std::min(1.0f,
            b.materialDamage[static_cast<size_t>(i)] + add * falloff * falloff);
        b.maxDamage = std::max(b.maxDamage, b.materialDamage[static_cast<size_t>(i)]);
    }
}

bool RigidBodyEngine::inspectLocalStructure(int wx, int wy, float &materialDamage, float &bondDamage, int &brokenNeighbors,
    float &strength, float &crack, float &wetness) const {
    materialDamage = bondDamage = strength = crack = wetness = 0.0f;
    brokenNeighbors = 0;
    int gx = wx, gy = wy;
    int bi = bodyAtCell(gx, gy);
    if (bi < 0 || bi >= static_cast<int>(bodies.size())) return false;
    RigidBody const &b = bodies[static_cast<size_t>(bi)];
    float lx, ly;
    worldToLocal(b, static_cast<float>(gx) + 0.5f, static_cast<float>(gy) + 0.5f, lx, ly);
    int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
    if (!maskOccupied(b, ix, iy)) return false;
    int i = iy * b.maskW + ix;
    materialDamage = (i < static_cast<int>(b.materialDamage.size())) ? b.materialDamage[static_cast<size_t>(i)] : 0.0f;
    wetness = pixelWetness(b, i);
    strength = effectiveStrength(b, i, StructuralDamageType::Generic);
    crack = getLocalCrackFactor(b, ix, iy);
    float maxBond = 0.0f;
    for (int n = 0; n < 4; ++n) {
        StructuralBond const *bond = bondPtr(b, ix, iy, ix + kN4x[n], iy + kN4y[n]);
        if (!bond) continue;
        maxBond = std::max(maxBond, bond->damage);
        if (bond->broken) ++brokenNeighbors;
    }
    bondDamage = maxBond;
    return true;
}

void RigidBodyEngine::depositImpactDamage(int bodyIndex, float wx, float wy, float speed, float impulse, float nx, float ny) {
    if (speed < kImpactSpeedMin) return;
    StructuralDamageRequest req;
    req.wx = wx;
    req.wy = wy;
    req.type = StructuralDamageType::Impact;
    req.dirX = nx;
    req.dirY = ny;
    req.canPropagate = true;
    float core = std::max(0.0f, speed - kImpactSpeedMin);
    float energyLike = 0.5f * core * core;
    float mass = 1.0f;
    if (bodyIndex >= 0 && bodyIndex < static_cast<int>(bodies.size()))
        mass = std::max(bodies[static_cast<size_t>(bodyIndex)].mass, 0.2f);
    if (impulse > 0.0f) energyLike = std::max(energyLike, 0.35f * impulse * core / mass);
    req.amount = energyLike / (kImpactEnergyScale * 0.67f);
    req.radius = 2.4f;
    applyStructuralDamage(bodyIndex, req);
}

void RigidBodyEngine::applyContactDamage(std::vector<RigidContact> const &contacts) {
    for (RigidContact const &c : contacts) {
        if (c.impactSpeed < kImpactSpeedMin) continue;
        if (c.jn < 0.08f && c.impactSpeed < kImpactSpeedMin + 1.5f) continue;
        depositImpactDamage(c.bodyA, c.x, c.y, c.impactSpeed, c.jn, c.nx, c.ny);
        if (c.bodyB >= 0) depositImpactDamage(c.bodyB, c.x, c.y, c.impactSpeed, c.jn, -c.nx, -c.ny);
    }
}

void RigidBodyEngine::splitBodyByConnectivity(int index, std::vector<RigidBody> &spawned) {
    RigidBody &b = bodies[static_cast<size_t>(index)];
    rebuildDerived(b);
    if (b.occupiedLocal.empty()) return;
    ensurePixelState(b);
    std::vector<uint8_t> seen(static_cast<size_t>(b.maskW * b.maskH), 0);
    std::vector<std::vector<int>> comps;
    for (int start : b.occupiedLocal) {
        if (seen[static_cast<size_t>(start)]) continue;
        std::vector<int> comp;
        std::vector<int> stack{start};
        seen[static_cast<size_t>(start)] = 1;
        while (!stack.empty()) {
            int i = stack.back(); stack.pop_back();
            comp.push_back(i);
            int lx = i % b.maskW, ly = i / b.maskW;
            for (int n = 0; n < 4; ++n) {
                int nx = lx + kN4x[n], ny = ly + kN4y[n];
                if (!bondConnects(b, lx, ly, nx, ny)) continue;
                int ni = ny * b.maskW + nx;
                if (seen[static_cast<size_t>(ni)]) continue;
                seen[static_cast<size_t>(ni)] = 1;
                stack.push_back(ni);
            }
        }
        comps.push_back(std::move(comp));
    }
    b.structureDirty = false;
    if (comps.size() <= 1) return;
    size_t keep = 0;
    size_t keepCount = comps[0].size();
    bool parentAnchored = b.anchored;
    for (size_t c = 1; c < comps.size(); ++c) {
        if (comps[c].size() > keepCount) {
            keep = c;
            keepCount = comps[c].size();
        }
    }
    uint32_t keepId = b.id;
    RigidBody kept = packFragment(b, comps[keep], parentAnchored);
    kept.id = keepId;
    for (size_t c = 0; c < comps.size(); ++c) {
        if (c == keep) continue;
        spawned.push_back(packFragment(b, comps[c], false));
    }
    b = std::move(kept);
}

void RigidBodyEngine::splitDirtyBodies() {
    lastFractureSplits = 0;
    std::vector<RigidBody> spawned;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (!bodies[i].structureDirty) continue;
        splitBodyByConnectivity(static_cast<int>(i), spawned);
    }
    if (spawned.empty()) return;
    lastFractureSplits = static_cast<int>(spawned.size());
    for (RigidBody &child : spawned) bodies.push_back(std::move(child));
    size_t write = 0;
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].occupiedLocal.empty()) continue;
        if (write != i) bodies[write] = std::move(bodies[i]);
        ++write;
    }
    bodies.resize(write);
    if (grab.active) {
        int gi = indexOfId(grab.bodyId);
        if (gi < 0) grab = GrabState{};
        else {
            RigidBody const &b = bodies[static_cast<size_t>(gi)];
            int lx = static_cast<int>(std::floor(grab.localX));
            int ly = static_cast<int>(std::floor(grab.localY));
            if (!maskOccupied(b, lx, ly)) grab = GrabState{};
        }
    }
}

void RigidBodyEngine::processMoisture(FluidEngine &fluid, float dt) {
    lastAbsorbed = 0.0;
    lastDried = 0.0;
    if (dt <= 0.0f) return;
    constexpr int ndx[4] = {-1, 1, 0, 0};
    constexpr int ndy[4] = {0, 0, -1, 1};
    for (RigidBody &b : bodies) {
        ensurePixelState(b);
        refreshMaterialCache(b);
        if (!b.moistureActive) continue;
        bool massDirty = false;
        bool touchingLiquid = false;
        for (int i : b.occupiedLocal) {
            int lx = i % b.maskW, ly = i / b.maskW;
            MaterialId id = b.mask[static_cast<size_t>(i)];
            if (!materialIsAbsorbent(id)) continue;
            MaterialDefinition const &mat = materialDef(id);
            float wx, wy;
            localToWorld(b, lx + 0.5f, ly + 0.5f, wx, wy);
            int gx = static_cast<int>(std::floor(wx));
            int gy = static_cast<int>(std::floor(wy));
            float &stored = b.moisture[static_cast<size_t>(i)];
            float room = mat.moistureCapacity - stored;
            if (room < 0.0f) room = 0.0f;
            float faceBudget = mat.absorptionRate * mat.porosity * dt;
            for (int n = -1; n < 4; ++n) {
                int nx = (n < 0) ? gx : gx + ndx[n];
                int ny = (n < 0) ? gy : gy + ndy[n];
                if (!FluidEngine::inside(nx, ny) || fluid.isSolid(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (occupant[static_cast<size_t>(ni)] >= 0) continue;
                float avail = fluid.fill[static_cast<size_t>(ni)];
                if (avail <= 1.0e-8f) continue;
                touchingLiquid = true;
                if (room <= 1.0e-8f || faceBudget <= 1.0e-8f) continue;
                float take = std::min(room, std::min(faceBudget, avail));
                if (take <= 1.0e-8f) continue;
                float actual = fluid.takeLiquidVolume(ni, take);
                if (actual <= 1.0e-12f) continue;
                stored += actual;
                room -= actual;
                fluid.expectedVolume -= actual;
                lastAbsorbed += actual;
                massDirty = true;
                fluid.wakeRegion(nx - 1, ny - 1, nx + 1, ny + 1);
            }
        }
        refreshMaterialCache(b);
        if (b.absorbedLiquid > 1.0e-6f) {
            if (moistureFlux.size() < b.moisture.size()) moistureFlux.assign(b.moisture.size(), 0.0f);
            else std::fill(moistureFlux.begin(), moistureFlux.begin() + static_cast<std::ptrdiff_t>(b.moisture.size()), 0.0f);
            double sum0 = 0.0;
            for (int i : b.occupiedLocal) sum0 += b.moisture[static_cast<size_t>(i)];
            for (int i : b.occupiedLocal) {
                int lx = i % b.maskW, ly = i / b.maskW;
                MaterialId id = b.mask[static_cast<size_t>(i)];
                MaterialDefinition const &mat = materialDef(id);
                if (mat.moistureCapacity <= 1.0e-8f || mat.permeability <= 1.0e-8f) continue;
                float ci = b.moisture[static_cast<size_t>(i)] / mat.moistureCapacity;
                for (int n = 0; n < 4; ++n) {
                    int nx = lx + ndx[n], ny = ly + ndy[n];
                    if (!maskOccupied(b, nx, ny)) continue;
                    int ni = ny * b.maskW + nx;
                    if (ni <= i) continue;
                    MaterialId oid = b.mask[static_cast<size_t>(ni)];
                    if (!materialIsAbsorbent(oid)) continue;
                    MaterialDefinition const &om = materialDef(oid);
                    if (om.moistureCapacity <= 1.0e-8f || om.permeability <= 1.0e-8f) continue;
                    float cj = b.moisture[static_cast<size_t>(ni)] / om.moistureCapacity;
                    float diff = ci - cj;
                    if (std::abs(diff) <= 1.0e-6f) continue;
                    float cond = 0.5f * (mat.permeability + om.permeability) * dt;
                    float move = diff * cond * std::min(mat.moistureCapacity, om.moistureCapacity);
                    if (diff > 0.0f) {
                        move = std::min(move, b.moisture[static_cast<size_t>(i)] * 0.22f);
                        float room = om.moistureCapacity - b.moisture[static_cast<size_t>(ni)];
                        move = std::min(move, std::max(0.0f, room));
                        moistureFlux[static_cast<size_t>(i)] -= move;
                        moistureFlux[static_cast<size_t>(ni)] += move;
                    } else {
                        move = std::min(-move, b.moisture[static_cast<size_t>(ni)] * 0.22f);
                        float room = mat.moistureCapacity - b.moisture[static_cast<size_t>(i)];
                        move = std::min(move, std::max(0.0f, room));
                        moistureFlux[static_cast<size_t>(ni)] -= move;
                        moistureFlux[static_cast<size_t>(i)] += move;
                    }
                }
            }
            for (int i : b.occupiedLocal) {
                MaterialDefinition const &mat = materialDef(b.mask[static_cast<size_t>(i)]);
                float next = b.moisture[static_cast<size_t>(i)] + moistureFlux[static_cast<size_t>(i)];
                if (next < 0.0f) next = 0.0f;
                if (mat.moistureCapacity > 0.0f && next > mat.moistureCapacity) next = mat.moistureCapacity;
                if (std::abs(next - b.moisture[static_cast<size_t>(i)]) > 1.0e-8f) massDirty = true;
                b.moisture[static_cast<size_t>(i)] = next;
            }
            double sum1 = 0.0;
            for (int i : b.occupiedLocal) sum1 += b.moisture[static_cast<size_t>(i)];
            double err = sum0 - sum1;
            if (std::abs(err) > 1.0e-7) {
                if (err > 0.0) {
                    for (int i : b.occupiedLocal) {
                        if (err <= 1.0e-10) break;
                        MaterialDefinition const &mat = materialDef(b.mask[static_cast<size_t>(i)]);
                        float room = mat.moistureCapacity - b.moisture[static_cast<size_t>(i)];
                        if (room <= 1.0e-8f) continue;
                        float add = static_cast<float>(std::min(static_cast<double>(room), err));
                        b.moisture[static_cast<size_t>(i)] += add;
                        err -= add;
                    }
                } else {
                    err = -err;
                    for (int i : b.occupiedLocal) {
                        if (err <= 1.0e-10) break;
                        float have = b.moisture[static_cast<size_t>(i)];
                        if (have <= 1.0e-8f) continue;
                        float sub = static_cast<float>(std::min(static_cast<double>(have), err));
                        b.moisture[static_cast<size_t>(i)] -= sub;
                        err -= sub;
                    }
                }
            }
        }
        for (int i : b.occupiedLocal) {
            int lx = i % b.maskW, ly = i / b.maskW;
            MaterialId id = b.mask[static_cast<size_t>(i)];
            MaterialDefinition const &mat = materialDef(id);
            if (mat.dryingRate <= 1.0e-8f || b.moisture[static_cast<size_t>(i)] <= 1.0e-8f) continue;
            float wx, wy;
            localToWorld(b, lx + 0.5f, ly + 0.5f, wx, wy);
            int gx = static_cast<int>(std::floor(wx));
            int gy = static_cast<int>(std::floor(wy));
            int tx = -1, ty = -1;
            float bestRoom = 0.0f;
            for (int n = -1; n < 4; ++n) {
                int nx = (n < 0) ? gx : gx + ndx[n];
                int ny = (n < 0) ? gy : gy + ndy[n];
                if (!FluidEngine::inside(nx, ny) || fluid.isSolid(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (occupant[static_cast<size_t>(ni)] >= 0) continue;
                float room = std::max(0.0f, 1.0f - fluid.fill[static_cast<size_t>(ni)]);
                if (room > bestRoom) { bestRoom = room; tx = nx; ty = ny; }
            }
            if (tx < 0 || bestRoom <= 1.0e-8f) continue;
            int ti = FluidEngine::ci(tx, ty);
            float roomNow = std::max(0.0f, 1.0f - fluid.fill[static_cast<size_t>(ti)]);
            float give = std::min(b.moisture[static_cast<size_t>(i)], mat.dryingRate * dt);
            give = std::min(give, roomNow);
            if (give <= 1.0e-8f) continue;
            b.moisture[static_cast<size_t>(i)] -= give;
            float cap = thermalCapacity(massKg(fluid.config.water.density, give, fluid.config.cellsPerMeter),
                fluid.config.water.thermal.specificHeat);
            fluid.addLiquidFill(ti, give, energyFromTemp(cap, AMBIENT_TEMPERATURE_K));
            fluid.expectedVolume += give;
            lastDried += give;
            massDirty = true;
            fluid.wakeRegion(tx - 1, ty - 1, tx + 1, ty + 1);
        }
        if (massDirty || touchingLiquid) {
            float oldCx = b.comLocalX, oldCy = b.comLocalY;
            computeMassProperties(b);
            applyComShift(b, oldCx, oldCy, 0, 0);
            if (!b.anchored && touchingLiquid && b.absorbedLiquid > 0.02f) {
                b.sleeping = false;
                b.quietTicks = 0;
            }
        } else refreshMaterialCache(b);
    }
}

int RigidBodyEngine::totalSolidPixels() const {
    int n = 0;
    for (RigidBody const &b : bodies) n += static_cast<int>(b.occupiedLocal.size());
    return n;
}

float RigidBodyEngine::totalSolidMass() const {
    float s = 0.0f;
    for (RigidBody const &b : bodies) s += b.mass;
    return s;
}

float RigidBodyEngine::totalAbsorbedLiquid() const {
    float s = 0.0f;
    for (RigidBody const &b : bodies) s += b.absorbedLiquid;
    return s;
}

int RigidBodyEngine::spawnPattern(FluidEngine &fluid, int ox, int oy, char const *const *rows, int rowCount, MaterialId material) {
    std::vector<int> cells;
    for (int y = 0; y < rowCount; ++y) {
        int w = static_cast<int>(std::strlen(rows[y]));
        for (int x = 0; x < w; ++x) {
            if (rows[y][x] != '#') continue;
            int gx = ox + x, gy = oy + y;
            if (!FluidEngine::inside(gx, gy) || fluid.isStaticSolid(gx, gy)) continue;
            cells.push_back(FluidEngine::ci(gx, gy));
        }
    }
    if (cells.empty()) return 0;
    bodies.push_back(makeBodyFromCells(cells, material));
    return 1;
}

void RigidBodyEngine::loadTestScene(FluidEngine &fluid, int scene) {
    clear();
    fluid.clearWorld();
    auto wall = [&](int x, int y) { if (FluidEngine::inside(x, y)) fluid.solid[FluidEngine::ci(x, y)] = 1; };
    auto waterRect = [&](int x0, int y0, int x1, int y1) {
        for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x)
            if (FluidEngine::inside(x, y) && !fluid.solid[FluidEngine::ci(x, y)]) {
                fluid.fill[FluidEngine::ci(x, y)] = 1.0f;
                fluid.waterShade[FluidEngine::ci(x, y)] = fluid.makeShade(x, y);
            }
    };
    for (int x = 20; x <= 180; ++x) wall(x, 112);
    for (int y = 20; y <= 112; ++y) { wall(20, y); wall(180, y); }

    if (scene == 1) {
        char const *rows[] = {"####", "####", "####"};
        spawnPattern(fluid, 96, 40, rows, 3, MATERIAL_STONE);
    } else if (scene == 2) {
        char const *rows[] = {"##...", "####.", ".###.", "..#.."};
        spawnPattern(fluid, 94, 36, rows, 4, MATERIAL_STONE);
    } else if (scene == 3) {
        char const *rows[] = {"########", "#......#", "#......#", "########"};
        spawnPattern(fluid, 96, 40, rows, 4, MATERIAL_STONE);
    } else if (scene == 4) {
        for (int y = 88; y <= 111; ++y) wall(108, y);
        char const *rows[] = {"#......", "#######", "#......"};
        spawnPattern(fluid, 90, 40, rows, 3, MATERIAL_STONE);
    } else if (scene == 5) {
        waterRect(21, 70, 179, 111);
        char const *rows[] = {"######", "######", "######"};
        spawnPattern(fluid, 97, 30, rows, 3, MATERIAL_STONE);
    } else if (scene == 6) {
        waterRect(21, 70, 179, 111);
        char const *rows[] = {"######", "######", "######", "######"};
        spawnPattern(fluid, 97, 28, rows, 4, MATERIAL_WOOD);
    } else if (scene == 7) {
        waterRect(21, 70, 179, 111);
        char const *rows[] = {"######", "######", "######", "######"};
        spawnPattern(fluid, 97, 28, rows, 4, MATERIAL_STONE);
    } else if (scene == 8) {
        waterRect(70, 88, 130, 111);
        char const *rows[] = {"#......#", "#......#", "#......#", "########"};
        spawnPattern(fluid, 96, 48, rows, 4, MATERIAL_WOOD);
    } else if (scene == 9) {
        char const *rows[] = {
            "....#....",
            "....#....",
            "....#....",
            "#########",
            "#.......#",
            "#########"
        };
        spawnPattern(fluid, 95, 40, rows, 6, MATERIAL_STONE);
    } else if (scene == 10) {
        waterRect(21, 70, 179, 111);
        char const *rows[] = {"####", "####"};
        spawnPattern(fluid, 40, 28, rows, 2, MATERIAL_STONE);
        spawnPattern(fluid, 90, 24, rows, 2, MATERIAL_WOOD);
        spawnPattern(fluid, 140, 32, rows, 2, MATERIAL_STONE);
    } else if (scene == 11) {
        char const *rows[] = {
            "..####..",
            ".######.",
            "########",
            "########",
            "########",
            "########",
            ".######.",
            "..####.."
        };
        spawnPattern(fluid, 96, 102, rows, 8, MATERIAL_STONE);
    } else if (scene == 12) {
        for (int x = 40; x <= 160; ++x) {
            int yTop = 70 + (x - 40) / 4;
            for (int y = yTop; y <= 112; ++y) wall(x, y);
        }
        char const *rows[] = {
            "..####..",
            ".######.",
            "########",
            "########",
            "########",
            "########",
            ".######.",
            "..####.."
        };
        spawnPattern(fluid, 48, 54, rows, 8, MATERIAL_STONE);
    } else if (scene == 13) {
        for (int x = 40; x <= 160; ++x) {
            int yTop = 70 + (160 - x) / 4;
            for (int y = yTop; y <= 112; ++y) wall(x, y);
        }
        char const *rows[] = {
            "..####..",
            ".######.",
            "########",
            "########",
            "########",
            "########",
            ".######.",
            "..####.."
        };
        spawnPattern(fluid, 136, 54, rows, 8, MATERIAL_STONE);
    } else if (scene == 14) {
        char const *rows[] = {
            "..####........",
            ".######.......",
            "########......",
            ".##########...",
            "..############",
            "...##########.",
            "......########",
            ".......######.",
            "........####.."
        };
        spawnPattern(fluid, 90, 100, rows, 9, MATERIAL_STONE);
    } else if (scene == 15) {
        char const *rows[] = {
            "..####..",
            ".######.",
            "########",
            "########",
            "########",
            "########",
            ".######.",
            "..####.."
        };
        spawnPattern(fluid, 96, 48, rows, 8, MATERIAL_STONE);
        if (!bodies.empty()) {
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
            bodies[0].sleeping = false;
            bodies[0].quietTicks = 0;
        }
    } else if (scene == 16) {
        char const *rows[] = {
            "..####..",
            ".######.",
            "########",
            "########",
            "########",
            "########",
            ".######.",
            "..####.."
        };
        spawnPattern(fluid, 96, 104, rows, 8, MATERIAL_STONE);
        if (!bodies.empty()) {
            bodies[0].y += 0.28f;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
            updateAabb(bodies[0]);
        }
    } else if (scene == 17) {
        char const *rows[] = {
            "######",
            "#.....",
            "#.....",
            "#....."
        };
        spawnPattern(fluid, 92, 86, rows, 4, MATERIAL_STONE);
    } else if (scene == 18) {
        char const *rows[] = {"######", "######", "######", "######"};
        spawnPattern(fluid, 97, 18, rows, 4, MATERIAL_GLASS);
    } else if (scene == 19) {
        char const *rows[] = {"######", "######", "######", "######"};
        spawnPattern(fluid, 97, 18, rows, 4, MATERIAL_METAL);
    } else if (scene == 20) {
        char const *rows[] = {"#####", "#####", "#####"};
        spawnPattern(fluid, 50, 18, rows, 3, MATERIAL_GLASS);
        spawnPattern(fluid, 140, 18, rows, 3, MATERIAL_METAL);
    } else if (scene == 21) {
        char const *rows[] = {"####################", "####################"};
        spawnPattern(fluid, 88, 36, rows, 2, MATERIAL_GLASS);
    } else if (scene == 22) {
        char const *rows[] = {"##########", "##########", "##########", "##########"};
        spawnPattern(fluid, 95, 88, rows, 4, MATERIAL_WOOD);
        if (!bodies.empty()) {
            bodies[0].anchored = true;
            bodies[0].sleeping = true;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
        }
        syncOccupancy(fluid);
        for (int y = 78; y <= 111; ++y) for (int x = 70; x <= 130; ++x)
            if (FluidEngine::inside(x, y) && !fluid.isSolid(x, y)) {
                fluid.fill[FluidEngine::ci(x, y)] = 1.0f;
                fluid.waterShade[FluidEngine::ci(x, y)] = fluid.makeShade(x, y);
            }
    } else if (scene == 23) {
        char const *wallRows[] = {
            "################",
            "################",
            "################",
            "################",
            "################",
            "################"
        };
        spawnPattern(fluid, 92, 70, wallRows, 6, MATERIAL_STONE);
        if (!bodies.empty()) {
            bodies[0].anchored = true;
            bodies[0].sleeping = true;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
        }
        char const *shot[] = {"####", "####", "####"};
        spawnPattern(fluid, 96, 16, shot, 3, MATERIAL_GLASS);
    } else if (scene == 24) {
        char const *rows[] = {"############", "############", "############"};
        spawnPattern(fluid, 94, 50, rows, 3, MATERIAL_GLASS);
        if (!bodies.empty()) {
            bodies[0].omega = 9.0f;
            bodies[0].sleeping = false;
            gravityScale = 0.0f;
        }
    } else if (scene == 25) {
        char const *rows[] = {
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########"
        };
        spawnPattern(fluid, 95, 84, rows, 8, MATERIAL_WOOD);
        if (!bodies.empty()) {
            bodies[0].anchored = true;
            bodies[0].sleeping = true;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
            gravityScale = 1.0f;
        }
        syncOccupancy(fluid);
        for (int y = 90; y <= 111; ++y) for (int x = 70; x <= 130; ++x)
            if (FluidEngine::inside(x, y) && !fluid.isSolid(x, y)) {
                fluid.fill[FluidEngine::ci(x, y)] = 1.0f;
                fluid.waterShade[FluidEngine::ci(x, y)] = fluid.makeShade(x, y);
            }
    }

    fluid.expectedVolume = 0.0;
    for (float amount : fluid.fill) fluid.expectedVolume += amount;
    fluid.wakeAllFluidChunks();
    fluid.enforceSolidBoundaries();
    fluid.rebuildActivityAndMetrics();
    syncOccupancy(fluid);
}

void RigidBodyEngine::runConservationBenchmark(FluidEngine &fluid) {
    loadTestScene(fluid, 10);
    double before = fluid.expectedVolume;
    int ticks = 240;
    auto start = FluidEngine::Clock::now();
    for (int i = 0; i < ticks; ++i) {
        step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gatherFluidForces(fluid);
    }
    double ms = FluidEngine::elapsedMs(start) / ticks;
    std::ofstream out(miscFile("rigid_benchmark_" + std::to_string(GW) + "x" + std::to_string(GH) + ".tsv"));
    out << "ticks\tbodies\tvolume\texpected\terror\tlost_rigid\tms_per_tick\trigid_ms\n";
    out << ticks << '\t' << bodies.size() << '\t' << fluid.currentVolume << '\t' << fluid.expectedVolume << '\t'
        << fluid.volumeError << '\t' << fluid.volumeLostRigid << '\t' << ms << '\t' << lastStepMs << '\n';
    (void)before;

    loadTestScene(fluid, 1);
    int settleTicks = 150;
    for (int i = 0; i < settleTicks; ++i) {
        step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        gatherFluidForces(fluid);
    }
    int sleeping = 0;
    float maxSpeed = 0.0f;
    for (RigidBody const &b : bodies) {
        if (b.sleeping) ++sleeping;
        maxSpeed = std::max(maxSpeed, std::sqrt(b.vx * b.vx + b.vy * b.vy) + std::abs(b.omega));
    }
    out << "settle_f1\t" << settleTicks << '\t' << bodies.size() << '\t' << sleeping << '\t' << maxSpeed << '\t'
        << lastContacts.size() << '\t' << lastStepMs;
    if (!bodies.empty())
        out << '\t' << (bodies[0].supported ? 1 : 0) << '\t' << bodies[0].maxPenetration << '\t' << bodies[0].y;
    out << '\n';
}

void RigidBodyEngine::runContactDiagnostics(FluidEngine &fluid) {
    std::ofstream out(miscFile("rigid_contact_diag_" + std::to_string(GW) + "x" + std::to_string(GH) + ".tsv"));
    out << "test\tresult\tdetail\tvx\tvy\tomega\tx\ty\tsleep\tsupported\tpen\tjn\tjt\tposx\tposy\n";
    int passed = 0, failed = 0;
    auto dump = [&](char const *name, bool ok, std::string const &detail) {
        RigidBody const *b = bodies.empty() ? nullptr : &bodies[0];
        out << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail;
        if (b) {
            out << '\t' << b->vx << '\t' << b->vy << '\t' << b->omega << '\t' << b->x << '\t' << b->y
                << '\t' << (b->sleeping ? 1 : 0) << '\t' << (b->supported ? 1 : 0) << '\t' << b->maxPenetration
                << '\t' << b->debugJn << '\t' << b->debugJt << '\t' << b->debugPosCorrX << '\t' << b->debugPosCorrY;
        }
        out << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto ticks = [&](int n) {
        for (int i = 0; i < n; ++i) step(fluid, PHYSICS_DT);
    };

    gravityScale = 1.0f;
    loadTestScene(fluid, 11);
    ticks(180);
    if (bodies.empty()) dump("flat_round", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        bool ok = std::abs(b.vx) < 0.30f && std::abs(b.omega) < 0.30f && std::abs(b.vy) < 0.50f;
        dump("flat_round", ok, ok ? "settled without spontaneous roll" : "still moving on flat floor");
    }

    loadTestScene(fluid, 12);
    ticks(90);
    if (bodies.empty()) dump("shallow_slope", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        bool ok = b.vx > 0.35f;
        dump("shallow_slope", ok, ok ? "moved downhill +x" : "did not roll/slide downhill");
    }

    loadTestScene(fluid, 13);
    ticks(90);
    if (bodies.empty()) dump("reverse_slope", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        bool ok = b.vx < -0.35f;
        dump("reverse_slope", ok, ok ? "moved downhill -x" : "did not reverse with slope");
    }

    loadTestScene(fluid, 14);
    ticks(80);
    float ke80 = 0.0f;
    if (!bodies.empty()) {
        RigidBody const &b = bodies[0];
        ke80 = 0.5f * b.mass * (b.vx * b.vx + b.vy * b.vy) + 0.5f * b.inertia * b.omega * b.omega;
    }
    ticks(120);
    if (bodies.empty()) dump("peanut_flat", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        float ke = 0.5f * b.mass * (b.vx * b.vx + b.vy * b.vy) + 0.5f * b.inertia * b.omega * b.omega;
        bool still = std::abs(b.vx) < 0.35f && std::abs(b.omega) < 0.35f && std::abs(b.vy) < 0.50f;
        bool noGrowth = ke <= ke80 * 1.6f + 0.08f;
        bool ok = still && noGrowth;
        dump("peanut_flat", ok, ok ? "settled without crawl/energy growth" : "asymmetric body still crawling");
    }

    gravityScale = 0.0f;
    loadTestScene(fluid, 15);
    float x0 = bodies.empty() ? 0.0f : bodies[0].x;
    float y0 = bodies.empty() ? 0.0f : bodies[0].y;
    if (!bodies.empty()) {
        bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
        bodies[0].sleeping = false;
    }
    ticks(120);
    if (bodies.empty()) dump("no_gravity", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        bool ok = std::abs(b.x - x0) < 0.08f && std::abs(b.y - y0) < 0.08f
            && std::abs(b.vx) < 0.05f && std::abs(b.vy) < 0.05f && std::abs(b.omega) < 0.05f;
        dump("no_gravity", ok, ok ? "remained stationary" : "moved with gravity off");
    }
    gravityScale = 1.0f;

    loadTestScene(fluid, 16);
    ticks(30);
    if (bodies.empty()) dump("slight_pen", false, "no body");
    else {
        RigidBody const &b = bodies[0];
        float speed = std::sqrt(b.vx * b.vx + b.vy * b.vy);
        bool ok = speed < 8.0f && std::abs(b.omega) < 4.0f && b.maxPenetration < 0.45f;
        dump("slight_pen", ok, ok ? "resolved overlap without launch" : "overlap launched or spun the body");
    }

    loadTestScene(fluid, 17);
    float peakOmega = 0.0f;
    for (int i = 0; i < 90; ++i) {
        step(fluid, PHYSICS_DT);
        if (!bodies.empty()) peakOmega = std::max(peakOmega, std::abs(bodies[0].omega));
    }
    if (bodies.empty()) dump("off_center", false, "no body");
    else {
        bool ok = peakOmega > 0.80f;
        dump("off_center", ok, ok
            ? ("contact impulse produced torque, peak|omega|=" + std::to_string(peakOmega))
            : ("no meaningful rotation on off-center hit, peak|omega|=" + std::to_string(peakOmega)));
    }

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
    gravityScale = 1.0f;
}

void RigidBodyEngine::runSolidDiagnostics(FluidEngine &fluid) {
    std::ofstream out(miscFile("solid_diag.tsv"));
    out << "test\tresult\tdetail\tbodies\tpixels\tmass\tmax_damage\tabsorbed\tvolume\n";
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\t'
            << bodies.size() << '\t' << totalSolidPixels() << '\t' << totalSolidMass() << '\t'
            << (bodies.empty() ? 0.0f : bodies[0].maxDamage) << '\t' << totalAbsorbedLiquid() << '\t'
            << fluid.currentVolume << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto ticks = [&](int n) {
        for (int i = 0; i < n; ++i) {
            step(fluid, PHYSICS_DT);
            fluid.simulationTick();
        }
        fluid.rebuildActivityAndMetrics();
    };

    gravityScale = 0.0f;
    loadTestScene(fluid, 1);
    gravityScale = 0.0f;
    if (bodies.empty()) emit("A_grab", false, "no body");
    else {
        bodies[0].sleeping = false;
        bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
        float x0 = bodies[0].x;
        bool grabbed = beginGrab(bodies[0].x, bodies[0].y);
        updateGrabTarget(x0 + 18.0f, bodies[0].y, false);
        ticks(24);
        endGrab();
        float x1 = bodies.empty() ? x0 : bodies[0].x;
        bool moved = x1 > x0 + 2.0f;
        bool notTeleported = bodies.empty() ? false : std::abs(bodies[0].x - (x0 + 18.0f)) > 0.35f;
        bool ok = grabbed && moved && notTeleported;
        emit("A_grab", ok, ok ? "followed spring constraint" :
            ("grab=" + std::to_string(grabbed) + " x0=" + std::to_string(x0) + " x1=" + std::to_string(x1)));
    }
    gravityScale = 1.0f;

    clear();
    fluid.clearWorld();
    for (int x = 20; x <= 180; ++x) if (FluidEngine::inside(x, 112)) fluid.solid[FluidEngine::ci(x, 112)] = 1;
    {
        char const *rows[] = {"####", "####", "####"};
        spawnPattern(fluid, 96, 100, rows, 3, MATERIAL_METAL);
    }
    fluid.wakeAllFluidChunks();
    syncOccupancy(fluid);
    ticks(50);
    float restDamage = bodies.empty() ? 1.0f : bodies[0].maxDamage;
    int restBodies = static_cast<int>(bodies.size());
    ticks(40);
    float laterDamage = bodies.empty() ? 1.0f : bodies[0].maxDamage;
    bool okB = restBodies == 1 && restDamage < 0.12f && laterDamage <= restDamage + 0.02f;
    emit("B_safe_drop", okB, okB ? "low metal drop did not fracture or idle-damage" : "low drop damaged or split");

    loadTestScene(fluid, 18);
    int glassPixels = totalSolidPixels();
    ticks(70);
    bool fractured = static_cast<int>(bodies.size()) > 1;
    int broken = 0;
    float maxBond = 0.0f;
    for (RigidBody const &b : bodies) {
        broken += b.brokenBondCount;
        maxBond = std::max(maxBond, b.maxBondDamage);
        fractured = fractured || b.maxDamage > 0.55f || b.brokenBondCount > 0;
    }
    bool keptMass = totalSolidPixels() == glassPixels;
    emit("C_hard_drop", fractured && keptMass, fractured
        ? "glass took impact damage/cracks"
        : ("glass survived high drop broken=" + std::to_string(broken) + " bond=" + std::to_string(maxBond)));

    loadTestScene(fluid, 20);
    ticks(70);
    int glassBits = 0, metalBits = 0;
    float glassD = 0.0f, metalD = 0.0f;
    int glassBroken = 0, metalBroken = 0;
    for (RigidBody const &b : bodies) {
        bool isGlass = false, isMetal = false;
        for (int i : b.occupiedLocal) {
            if (b.mask[static_cast<size_t>(i)] == MATERIAL_GLASS) isGlass = true;
            if (b.mask[static_cast<size_t>(i)] == MATERIAL_METAL) isMetal = true;
        }
        if (isGlass) { ++glassBits; glassD = std::max(glassD, b.maxDamage); glassBroken += b.brokenBondCount; }
        if (isMetal) { ++metalBits; metalD = std::max(metalD, b.maxDamage); metalBroken += b.brokenBondCount; }
    }
    bool okD = (glassBits > 1 || glassBroken > metalBroken || glassD > metalD + 0.12f) && metalBits == 1;
    emit("D_material_diff", okD, okD ? "glass failed more than metal" : "materials behaved too similarly");

    loadTestScene(fluid, 1);
    if (bodies.empty()) emit("E_fatigue", false, "no body");
    else {
        int id = 0;
        RigidBody &b = bodies[0];
        float wx, wy;
        localToWorld(b, b.comLocalX, static_cast<float>(b.maskH) - 0.4f, wx, wy);
        for (int n = 0; n < 8; ++n) depositImpactDamage(id, wx, wy, 12.0f, 4.0f);
        splitDirtyBodies();
        int midBodies = static_cast<int>(bodies.size());
        float midD = bodies.empty() ? 0.0f : bodies[0].maxDamage;
        for (int n = 0; n < 12; ++n) {
            if (bodies.empty()) break;
            localToWorld(bodies[0], bodies[0].comLocalX, static_cast<float>(bodies[0].maskH) - 0.4f, wx, wy);
            depositImpactDamage(0, wx, wy, 12.0f, 4.0f);
        }
        splitDirtyBodies();
        bool grew = midD > 0.04f;
        bool later = static_cast<int>(bodies.size()) > midBodies || (!bodies.empty() && bodies[0].maxDamage > midD + 0.03f);
        emit("E_fatigue", grew && later, grew ? "subcritical hits accumulated" : "no accumulation");
    }

    loadTestScene(fluid, 21);
    if (bodies.empty()) emit("F_asymmetric", false, "no body");
    else {
        RigidBody &b = bodies[0];
        float leftD = 0.0f, rightD = 0.0f;
        depositImpactDamage(0, b.aabbX0 + 0.6f, b.y, 38.0f, 12.0f);
        ensurePixelState(b);
        for (int i : b.occupiedLocal) {
            int lx = i % b.maskW;
            if (lx < b.maskW / 3) leftD = std::max(leftD, b.materialDamage[static_cast<size_t>(i)]);
            if (lx > (2 * b.maskW) / 3) rightD = std::max(rightD, b.materialDamage[static_cast<size_t>(i)]);
        }
        bool okF = leftD > rightD + 0.15f;
        emit("F_asymmetric", okF, okF ? "damage stayed near struck end" : "damage spread across beam");
    }

    loadTestScene(fluid, 22);
    ticks(1);
    auto totalLiquid = [&]() {
        double s = 0.0;
        for (float amount : fluid.fill) s += amount;
        for (SplashParticle const &p : fluid.splashes) s += p.volume;
        s += totalAbsorbedLiquid();
        return s;
    };
    auto freeLiquid = [&]() {
        double s = 0.0;
        for (float amount : fluid.fill) s += amount;
        return s;
    };
    double water0 = freeLiquid();
    float abs0 = totalAbsorbedLiquid();
    double tot0 = totalLiquid();
    ticks(90);
    double water1 = freeLiquid();
    float abs1 = totalAbsorbedLiquid();
    double tot1 = totalLiquid();
    bool took = abs1 > abs0 + 0.3f;
    bool waterDown = water1 < water0 - 0.25f;
    bool cons = std::abs(tot1 - tot0) < 0.12;
    emit("G_absorb", took && waterDown && cons, took
        ? ("wood absorbed free water abs=" + std::to_string(abs1) + " dTotal=" + std::to_string(tot1 - tot0))
        : "no conservative absorption");

    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"##########", "##########"};
        spawnPattern(fluid, 90, 40, rows, 2, MATERIAL_WOOD);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("H_spread", false, "no body");
    else {
        RigidBody &b = bodies[0];
        ensurePixelState(b);
        for (int i : b.occupiedLocal) {
            int lx = i % b.maskW;
            if (lx <= 1) b.moisture[static_cast<size_t>(i)] = materialDef(MATERIAL_WOOD).moistureCapacity * 0.95f;
        }
        computeMassProperties(b);
        float right0 = 0.0f, right1 = 0.0f, left1 = 0.0f;
        auto sample = [&](float &left, float &right) {
            left = right = 0.0f;
            for (int i : b.occupiedLocal) {
                int lx = i % b.maskW;
                if (lx <= 1) left += b.moisture[static_cast<size_t>(i)];
                if (lx >= b.maskW - 2) right += b.moisture[static_cast<size_t>(i)];
            }
        };
        sample(left1, right0);
        for (int i = 0; i < 180; ++i) processMoisture(fluid, PHYSICS_DT);
        sample(left1, right1);
        float mid = 0.0f;
        for (int i : b.occupiedLocal) {
            int lx = i % b.maskW;
            if (lx >= 3 && lx <= 5) mid += b.moisture[static_cast<size_t>(i)];
        }
        bool okH = mid > 0.05f && right1 < left1;
        emit("H_spread", okH, okH ? "moisture diffused along wood" : "no gradual spread");
    }

    auto smash = [&](MaterialId mat, bool saturate) {
        clear();
        fluid.clearWorld();
        for (int x = 20; x <= 180; ++x) if (FluidEngine::inside(x, 112)) fluid.solid[FluidEngine::ci(x, 112)] = 1;
        char const *rows[] = {"########", "########"};
        spawnPattern(fluid, 96, 80, rows, 2, mat);
        if (!bodies.empty()) {
            RigidBody &b = bodies[0];
            ensurePixelState(b);
            if (saturate) {
                for (int i : b.occupiedLocal) b.moisture[static_cast<size_t>(i)] = materialDef(mat).moistureCapacity;
            }
            computeMassProperties(b);
            float wx, wy;
            localToWorld(b, b.comLocalX, static_cast<float>(b.maskH) - 0.4f, wx, wy);
            depositImpactDamage(0, wx, wy, 16.0f, 6.0f);
            splitDirtyBodies();
        }
        int n = static_cast<int>(bodies.size());
        float dmg = 0.0f;
        for (RigidBody const &b : bodies) dmg = std::max(dmg, b.maxDamage);
        return std::pair<int, float>(n, dmg);
    };
    auto dry = smash(MATERIAL_WOOD, false);
    auto wet = smash(MATERIAL_WOOD, true);
    bool okI = wet.second > dry.second + 0.04f || wet.first > dry.first;
    emit("I_wet_fracture", okI, okI
        ? "wet wood weaker than dry"
        : ("dry n/d=" + std::to_string(dry.first) + "/" + std::to_string(dry.second)
            + " wet n/d=" + std::to_string(wet.first) + "/" + std::to_string(wet.second)));

    loadTestScene(fluid, 18);
    int pix0 = totalSolidPixels();
    float mass0 = totalSolidMass();
    ticks(70);
    int pix1 = totalSolidPixels();
    float mass1 = totalSolidMass();
    bool okJ = pix1 == pix0 && std::abs(mass1 - mass0) < 0.08f;
    emit("J_fragment_mass", okJ, okJ ? "pixels and mass conserved" : "fracture lost or created mass");

    gravityScale = 0.0f;
    loadTestScene(fluid, 24);
    if (bodies.size() != 1) emit("K_spin_split", false, "expected one spinning body");
    else {
        RigidBody &b = bodies[0];
        b.omega = 8.0f;
        depositImpactDamage(0, b.aabbX0 + 0.5f, b.y, 50.0f, 20.0f);
        splitDirtyBodies();
        bool okK = false;
        if (bodies.size() >= 2) {
            float dvx = bodies[0].vx - bodies[1].vx;
            float dvy = bodies[0].vy - bodies[1].vy;
            okK = std::sqrt(dvx * dvx + dvy * dvy) > 0.40f;
        }
        emit("K_spin_split", okK, okK ? "fragments inherited spin velocities" : "fragments shared one velocity");
    }
    gravityScale = 1.0f;

    loadTestScene(fluid, 23);
    ticks(70);
    int anchored = 0, loose = 0;
    for (RigidBody const &b : bodies) {
        if (b.anchored) ++anchored;
        else ++loose;
    }
    bool okL = anchored >= 1 && (loose >= 1 || (!bodies.empty() && (bodies[0].maxDamage > 0.4f || bodies[0].brokenBondCount > 0)));
    emit("L_anchored", okL, okL ? "anchored mass remained after impact" : "anchored fracture failed");

    gravityScale = 0.0f;
    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"########", "########", "########", "########"};
        spawnPattern(fluid, 90, 40, rows, 4, MATERIAL_STONE);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("M_crack_persist", false, "no body");
    else {
        int h = bodies[0].maskH;
        for (int y = 0; y < h - 1; ++y) breakBond(0, 3, y, 4, y);
        splitDirtyBodies();
        int broken0 = bodies.empty() ? 0 : bodies[0].brokenBondCount;
        int n0 = static_cast<int>(bodies.size());
        if (!bodies.empty()) {
            bodies[0].sleeping = true;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
        }
        ticks(36);
        if (!bodies.empty()) {
            bodies[0].sleeping = false;
            bodies[0].quietTicks = 0;
        }
        ticks(8);
        int broken1 = 0;
        for (RigidBody const &b : bodies) broken1 += b.brokenBondCount;
        bool okM = n0 == 1 && static_cast<int>(bodies.size()) == 1 && broken0 >= 3 && broken1 == broken0;
        emit("M_crack_persist", okM, okM ? "crack survived sleep/wake" :
            ("n0=" + std::to_string(n0) + " n=" + std::to_string(bodies.size())
                + " broken0=" + std::to_string(broken0) + " broken1=" + std::to_string(broken1)));
    }

    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"########", "########", "########", "########"};
        spawnPattern(fluid, 90, 40, rows, 4, MATERIAL_STONE);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("N_crack_no_split", false, "no body");
    else {
        int h = bodies[0].maskH;
        for (int y = 0; y < h - 1; ++y) breakBond(0, 3, y, 4, y);
        splitDirtyBodies();
        bool one = bodies.size() == 1;
        int brokenN = bodies.empty() ? 0 : bodies[0].brokenBondCount;
        bool okN = one && brokenN >= 3;
        emit("N_crack_no_split", okN, okN ? "crack kept a single rigid body" :
            ("bodies=" + std::to_string(bodies.size()) + " broken=" + std::to_string(brokenN)));
    }

    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"########", "########", "########", "########"};
        spawnPattern(fluid, 90, 40, rows, 4, MATERIAL_STONE);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("O_full_split", false, "no body");
    else {
        int h = bodies[0].maskH;
        int pix = totalSolidPixels();
        for (int y = 0; y < h; ++y) breakBond(0, 3, y, 4, y);
        splitDirtyBodies();
        bool okO = bodies.size() == 2 && totalSolidPixels() == pix;
        emit("O_full_split", okO, okO ? "through-crack split into two bodies" :
            ("bodies=" + std::to_string(bodies.size()) + " pix=" + std::to_string(totalSolidPixels())));
    }

    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"##########", "##########", "##########", "##########"};
        spawnPattern(fluid, 88, 40, rows, 4, MATERIAL_STONE);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("P_fragment_crack", false, "no body");
    else {
        ensurePixelState(bodies[0]);
        for (int i : bodies[0].occupiedLocal) {
            int lx = i % bodies[0].maskW;
            if (lx <= 2) bodies[0].moisture[static_cast<size_t>(i)] = 0.35f;
        }
        computeMassProperties(bodies[0]);
        int h = bodies[0].maskH;
        for (int y = 0; y < h - 1; ++y) breakBond(0, 1, y, 2, y);
        splitDirtyBodies();
        int preBodies = static_cast<int>(bodies.size());
        int preBroken = bodies.empty() ? 0 : bodies[0].brokenBondCount;
        if (!bodies.empty()) {
            for (int y = 0; y < bodies[0].maskH; ++y) breakBond(0, 6, y, 7, y);
        }
        splitDirtyBodies();
        bool retained = false;
        bool wetKept = false;
        for (RigidBody const &b : bodies) {
            if (b.brokenBondCount > 0 && static_cast<int>(b.occupiedLocal.size()) >= 8) retained = true;
            if (b.absorbedLiquid > 0.10f && b.brokenBondCount > 0) wetKept = true;
        }
        bool okP = preBodies == 1 && preBroken >= 3 && bodies.size() >= 2 && retained && wetKept;
        emit("P_fragment_crack", okP, okP ? "fragment kept internal crack" :
            ("preN=" + std::to_string(preBodies) + " n=" + std::to_string(bodies.size())
                + " preBroken=" + std::to_string(preBroken) + " retained=" + std::to_string(retained)
                + " wet=" + std::to_string(wetKept)));
    }

    clear();
    fluid.clearWorld();
    {
        char const *rows[] = {"######", "######", "######"};
        spawnPattern(fluid, 94, 40, rows, 3, MATERIAL_WOOD);
    }
    syncOccupancy(fluid);
    if (bodies.empty()) emit("Q_generic_api", false, "no body");
    else {
        float wx, wy;
        localToWorld(bodies[0], bodies[0].comLocalX, bodies[0].comLocalY, wx, wy);
        StructuralDamageRequest req;
        req.wx = wx;
        req.wy = wy;
        req.amount = 0.55f;
        req.radius = 1.8f;
        req.dirX = 1.0f;
        req.dirY = 0.0f;
        req.type = StructuralDamageType::Chemical;
        req.canPropagate = true;
        applyStructuralDamage(0, req);
        splitDirtyBodies();
        float md = 0.0f, bd = 0.0f;
        int brokenQ = 0;
        for (RigidBody const &b : bodies) {
            md = std::max(md, b.maxDamage);
            bd = std::max(bd, b.maxBondDamage);
            brokenQ += b.brokenBondCount;
        }
        bool okQ = md > 0.08f || bd > 0.20f || brokenQ > 0;
        emit("Q_generic_api", okQ, okQ ? "external structural API damaged bonds" :
            ("md=" + std::to_string(md) + " bd=" + std::to_string(bd) + " broken=" + std::to_string(brokenQ)));
    }

    gravityScale = 1.0f;
    clear();
    fluid.clearWorld();
    for (int x = 20; x <= 180; ++x) if (FluidEngine::inside(x, 112)) fluid.solid[FluidEngine::ci(x, 112)] = 1;
    {
        char const *rows[] = {"####", "####", "####"};
        spawnPattern(fluid, 96, 102, rows, 3, MATERIAL_STONE);
    }
    fluid.wakeAllFluidChunks();
    syncOccupancy(fluid);
    ticks(40);
    float rest0 = 0.0f;
    int brokenRest0 = 0;
    for (RigidBody const &b : bodies) {
        rest0 = std::max(rest0, b.maxDamage);
        brokenRest0 += b.brokenBondCount;
    }
    ticks(160);
    float rest1 = 0.0f;
    int brokenRest1 = 0;
    for (RigidBody const &b : bodies) {
        rest1 = std::max(rest1, b.maxDamage);
        brokenRest1 += b.brokenBondCount;
    }
    bool okR = bodies.size() == 1 && rest1 <= rest0 + 0.02f && brokenRest1 == brokenRest0 && rest1 < 0.12f;
    emit("R_resting", okR, okR ? "resting body did not self-damage" :
        ("d0=" + std::to_string(rest0) + " d1=" + std::to_string(rest1)
            + " b0=" + std::to_string(brokenRest0) + " b1=" + std::to_string(brokenRest1)));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
    gravityScale = 1.0f;
    endGrab();
}

void RigidBodyEngine::runMoistureDiagnostics(FluidEngine &fluid) {
    std::ofstream out(miscFile("moisture_diag.tsv"));
    out << std::setprecision(8);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto splashVol = [&]() {
        double s = 0.0;
        for (SplashParticle const &p : fluid.splashes) s += p.volume;
        return s;
    };
    auto freeVol = [&]() {
        double s = 0.0;
        for (float amount : fluid.fill) s += amount;
        return s;
    };
    auto sample = [&](int &wetPixels, float &avgWet, float &maxWet, float &mass, float &cx, float &cy) {
        wetPixels = 0;
        avgWet = 0.0f;
        maxWet = 0.0f;
        mass = totalSolidMass();
        cx = 0.0f;
        cy = 0.0f;
        int n = 0;
        double wetSum = 0.0;
        for (RigidBody const &b : bodies) {
            cx = b.x;
            cy = b.y;
            for (int i : b.occupiedLocal) {
                float w = pixelWetness(b, i);
                wetSum += w;
                maxWet = std::max(maxWet, w);
                if (b.moisture[static_cast<size_t>(i)] > 1.0e-4f) ++wetPixels;
                ++n;
            }
        }
        avgWet = (n > 0) ? static_cast<float>(wetSum / n) : 0.0f;
    };
    auto tickOnce = [&]() {
        step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        fluid.rebuildActivityAndMetrics();
    };
    auto soak = [&](MaterialId mat, bool rotate) {
        clear();
        fluid.clearWorld();
        gravityScale = 1.0f;
        auto wall = [&](int x, int y) { if (FluidEngine::inside(x, y)) fluid.solid[FluidEngine::ci(x, y)] = 1; };
        for (int x = 20; x <= 180; ++x) wall(x, 112);
        for (int y = 20; y <= 112; ++y) { wall(20, y); wall(180, y); }
        char const *rows[] = {
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########",
            "##########"
        };
        spawnPattern(fluid, 95, 84, rows, 8, mat);
        if (!bodies.empty()) {
            bodies[0].anchored = true;
            bodies[0].sleeping = true;
            bodies[0].vx = bodies[0].vy = bodies[0].omega = 0.0f;
            if (rotate) {
                bodies[0].theta = 0.55f;
                updateAabb(bodies[0]);
            }
        }
        syncOccupancy(fluid);
        for (int y = 90; y <= 111; ++y) for (int x = 70; x <= 130; ++x)
            if (FluidEngine::inside(x, y) && !fluid.isSolid(x, y)) {
                fluid.fill[FluidEngine::ci(x, y)] = 1.0f;
                fluid.waterShade[FluidEngine::ci(x, y)] = fluid.makeShade(x, y);
            }
        fluid.expectedVolume = 0.0;
        for (float amount : fluid.fill) fluid.expectedVolume += amount;
        fluid.wakeAllFluidChunks();
        fluid.enforceSolidBoundaries();
        fluid.rebuildActivityAndMetrics();
        syncOccupancy(fluid);
    };

    out << "section\twood_soak\n";
    out << "tick\tfree_liquid\tabsorbed_liquid\tsplash_liquid\ttotal_liquid\taverage_wetness\tmax_wetness\twet_pixel_count\trigid_body_mass\tcenter_of_mass_x\tcenter_of_mass_y\n";
    soak(MATERIAL_WOOD, false);
    double total0 = 0.0;
    double free0 = 0.0;
    float mass0 = 0.0f;
    float abs0 = 0.0f;
    auto start = FluidEngine::Clock::now();
    constexpr int kWoodTicks = 180;
    for (int tick = 0; tick <= kWoodTicks; ++tick) {
        if (tick > 0) tickOnce();
        double freeL = freeVol();
        double splashL = splashVol();
        double absL = totalAbsorbedLiquid();
        double totalL = freeL + absL + splashL;
        int wetPixels = 0;
        float avgWet = 0.0f, maxWet = 0.0f, mass = 0.0f, cx = 0.0f, cy = 0.0f;
        sample(wetPixels, avgWet, maxWet, mass, cx, cy);
        if (tick == 0) {
            total0 = totalL;
            free0 = freeL;
            mass0 = mass;
            abs0 = static_cast<float>(absL);
        }
        out << tick << '\t' << freeL << '\t' << absL << '\t' << splashL << '\t' << totalL << '\t'
            << avgWet << '\t' << maxWet << '\t' << wetPixels << '\t' << mass << '\t' << cx << '\t' << cy << '\n';
    }
    double moistureMs = FluidEngine::elapsedMs(start) / kWoodTicks;
    double free1 = freeVol();
    float abs1 = totalAbsorbedLiquid();
    double splash1 = splashVol();
    double total1 = free1 + abs1 + splash1;
    int wetPixels = 0;
    float avgWet = 0.0f, maxWet = 0.0f, mass1 = 0.0f, cx = 0.0f, cy = 0.0f;
    sample(wetPixels, avgWet, maxWet, mass1, cx, cy);
    float interior = 0.0f, surface = 0.0f;
    if (!bodies.empty()) {
        RigidBody const &b = bodies[0];
        for (int i : b.occupiedLocal) {
            int ly = i / b.maskW;
            if (ly >= b.maskH - 2) surface += b.moisture[static_cast<size_t>(i)];
            if (ly <= 2) interior += b.moisture[static_cast<size_t>(i)];
        }
    }
    bool woodAbsorbed = abs1 > abs0 + 1.2f;
    bool freeDropped = free1 < free0 - 1.0f;
    bool conserved = std::abs(total1 - total0) < 0.12;
    bool spreadIn = interior > 0.05f && surface > interior;
    bool heavier = mass1 > mass0 + 1.0f;
    emit("A_wood", woodAbsorbed && freeDropped && conserved && spreadIn && heavier,
        "abs=" + std::to_string(abs1) + " dFree=" + std::to_string(free0 - free1)
            + " dTotal=" + std::to_string(total1 - total0) + " interior=" + std::to_string(interior)
            + " surface=" + std::to_string(surface) + " dMass=" + std::to_string(mass1 - mass0)
            + " splash=" + std::to_string(splash1));

    auto soakMeasure = [&](MaterialId mat, bool rotate) {
        soak(mat, rotate);
        double t0 = freeVol() + totalAbsorbedLiquid() + splashVol();
        float a0 = totalAbsorbedLiquid();
        for (int i = 0; i < 180; ++i) tickOnce();
        float a1 = totalAbsorbedLiquid();
        double t1 = freeVol() + a1 + splashVol();
        return std::tuple<float, float, double, double>(a0, a1, t0, t1);
    };

    auto stoneR = soakMeasure(MATERIAL_STONE, false);
    bool stoneOk = std::get<1>(stoneR) > 0.008f && std::get<1>(stoneR) < abs1 * 0.45f
        && std::abs(std::get<3>(stoneR) - std::get<2>(stoneR)) < 0.12;
    emit("B_stone", stoneOk, stoneOk ? "stone absorbed slowly and weakly"
        : ("stone_abs=" + std::to_string(std::get<1>(stoneR)) + " wood_abs=" + std::to_string(abs1)));

    auto glassR = soakMeasure(MATERIAL_GLASS, false);
    emit("C_glass", std::get<1>(glassR) < 1.0e-4f && std::abs(std::get<3>(glassR) - std::get<2>(glassR)) < 0.12,
        std::get<1>(glassR) < 1.0e-4f ? "glass absorbed ~0" : ("glass_abs=" + std::to_string(std::get<1>(glassR))));

    auto metalR = soakMeasure(MATERIAL_METAL, false);
    emit("D_metal", std::get<1>(metalR) < 1.0e-4f && std::abs(std::get<3>(metalR) - std::get<2>(metalR)) < 0.12,
        std::get<1>(metalR) < 1.0e-4f ? "metal absorbed ~0" : ("metal_abs=" + std::to_string(std::get<1>(metalR))));

    auto rotR = soakMeasure(MATERIAL_WOOD, true);
    emit("E_rotated_wood", std::get<1>(rotR) > 0.6f && std::abs(std::get<3>(rotR) - std::get<2>(rotR)) < 0.12,
        std::get<1>(rotR) > 0.8f ? "rotated wood still absorbed"
            : ("rot_abs=" + std::to_string(std::get<1>(rotR))));

    soak(MATERIAL_WOOD, false);
    if (bodies.empty()) emit("F_fracture", false, "no body");
    else {
        RigidBody &b = bodies[0];
        ensurePixelState(b);
        for (int i : b.occupiedLocal)
            b.moisture[static_cast<size_t>(i)] = materialDef(MATERIAL_WOOD).moistureCapacity * 0.80f;
        computeMassProperties(b);
        float parentAbs = b.absorbedLiquid;
        int mid = b.maskW / 2;
        for (int ly = 0; ly < b.maskH; ++ly)
            breakBond(0, mid - 1, ly, mid, ly);
        splitDirtyBodies();
        float childAbs = totalAbsorbedLiquid();
        bool split = static_cast<int>(bodies.size()) > 1;
        bool kept = std::abs(childAbs - parentAbs) < 0.08f;
        emit("F_fracture", split && kept,
            split && kept ? ("fragments conserved moisture, bodies=" + std::to_string(bodies.size()))
                : ("split=" + std::to_string(split) + " parent=" + std::to_string(parentAbs)
                    + " children=" + std::to_string(childAbs)));
    }

    out << "perf_ms_per_tick\t" << moistureMs << '\n';
    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
    gravityScale = 1.0f;
}
