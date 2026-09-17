#include "render/WorldRenderer.h"

#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "rigid/RigidBodyEngine.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

namespace {

constexpr int kN4x[4] = {-1, 1, 0, 0};
constexpr int kN4y[4] = {0, 0, -1, 1};

int clampByte(int v) {
    return std::clamp(v, 0, 255);
}

uint32_t packRgb(int r, int g, int b) {
    return static_cast<uint32_t>(clampByte(b))
        | (static_cast<uint32_t>(clampByte(g)) << 8)
        | (static_cast<uint32_t>(clampByte(r)) << 16);
}

void tintToward(int &r, int &g, int &b, int r1, int g1, int b1, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    r = static_cast<int>(std::lround(r + (r1 - r) * t));
    g = static_cast<int>(std::lround(g + (g1 - g) * t));
    b = static_cast<int>(std::lround(b + (b1 - b) * t));
}

uint32_t scaleRgb(int r, int g, int b, float s) {
    return packRgb(static_cast<int>(std::lround(static_cast<float>(r) * s)),
        static_cast<int>(std::lround(static_cast<float>(g) * s)),
        static_cast<int>(std::lround(static_cast<float>(b) * s)));
}

uint32_t lerpRgb(int r0, int g0, int b0, int r1, int g1, int b1, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    return packRgb(static_cast<int>(std::lround(r0 + (r1 - r0) * t)),
        static_cast<int>(std::lround(g0 + (g1 - g0) * t)),
        static_cast<int>(std::lround(b0 + (b1 - b0) * t)));
}

uint32_t mixToward(uint32_t dst, int r, int g, int b, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    int dr = static_cast<int>((dst >> 16) & 255);
    int dg = static_cast<int>((dst >> 8) & 255);
    int db = static_cast<int>(dst & 255);
    return lerpRgb(dr, dg, db, r, g, b, t);
}

// Rendering only. Full cells stay opaque; tiny fill stays visible.
float liquidVisualAlpha(float fill) {
    float f = std::clamp(fill, 0.0f, 1.0f);
    return 0.15f + 0.85f * std::sqrt(f);
}

uint32_t applyLiquidFillAlpha(uint32_t color, uint32_t bg, float fill) {
    int cr = static_cast<int>((color >> 16) & 255);
    int cg = static_cast<int>((color >> 8) & 255);
    int cb = static_cast<int>(color & 255);
    return mixToward(bg, cr, cg, cb, liquidVisualAlpha(fill));
}

float stableNoise(int x, int y, uint32_t salt) {
    uint32_t h = FluidEngine::hashCell(static_cast<uint32_t>(x), static_cast<uint32_t>(y), salt);
    return (static_cast<int>(h & 255u) / 127.5f) - 1.0f;
}

void lerpVisual(MaterialVisual &a, MaterialVisual const &b, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    a.r = static_cast<int>(std::lround(a.r + (b.r - a.r) * t));
    a.g = static_cast<int>(std::lround(a.g + (b.g - a.g) * t));
    a.b = static_cast<int>(std::lround(a.b + (b.b - a.b) * t));
    a.noiseStrength += (b.noiseStrength - a.noiseStrength) * t;
    a.depthStrength += (b.depthStrength - a.depthStrength) * t;
    a.realisticDetail += (b.realisticDetail - a.realisticDetail) * t;
    a.glowStrength += (b.glowStrength - a.glowStrength) * t;
    a.outlineDarken += (b.outlineDarken - a.outlineDarken) * t;
    a.alphaResponse += (b.alphaResponse - a.alphaResponse) * t;
    a.glowEligible = a.glowEligible || b.glowEligible;
}

void applyDye(MaterialVisual &vis, float ir, float ig, float ib) {
    float dyeI = std::clamp(std::max(ir, std::max(ig, ib)), 0.0f, 1.0f);
    if (dyeI > 0.02f)
        tintToward(vis.r, vis.g, vis.b, static_cast<int>(ir * 255.0f), static_cast<int>(ig * 255.0f),
            static_cast<int>(ib * 255.0f), dyeI);
}

MaterialVisual blendLiquidVisual(LiquidComponent const *items, int count) {
    float tot = 0.0f;
    for (int n = 0; n < count; ++n) {
        if (validLiquidComponentId(items[n].id) && items[n].amount > kMinLiquidComponent)
            tot += items[n].amount;
    }
    if (!(tot > kMinLiquidComponent))
        return visualForSubstance(SUBSTANCE_WATER);
    MaterialVisual vis{};
    float acc = 0.0f;
    bool first = true;
    for (int n = 0; n < count; ++n) {
        if (!validLiquidComponentId(items[n].id) || items[n].amount <= kMinLiquidComponent)
            continue;
        if (first) {
            vis = visualForSubstance(items[n].id);
            acc = items[n].amount;
            first = false;
            continue;
        }
        float w = items[n].amount / (acc + items[n].amount);
        lerpVisual(vis, visualForSubstance(items[n].id), w);
        acc += items[n].amount;
    }
    return vis;
}

MaterialVisual sampleLiquidCell(FluidEngine const &fluid, int index) {
    LiquidComponentView view = fluid.liquidComponents(index);
    MaterialVisual vis = blendLiquidVisual(view.items, view.count);
    float f = std::max(fluid.fill[static_cast<size_t>(index)], 1.0e-8f);
    applyDye(vis,
        std::clamp(fluid.dyeR[static_cast<size_t>(index)] / f, 0.0f, 1.0f),
        std::clamp(fluid.dyeG[static_cast<size_t>(index)] / f, 0.0f, 1.0f),
        std::clamp(fluid.dyeB[static_cast<size_t>(index)] / f, 0.0f, 1.0f));
    return vis;
}

MaterialVisual sampleSplash(SplashParticle const &p) {
    MaterialVisual vis = blendLiquidVisual(p.comps, p.compCount);
    float vol = std::max(p.volume, 1.0e-8f);
    applyDye(vis,
        std::clamp(p.dyeR / vol, 0.0f, 1.0f),
        std::clamp(p.dyeG / vol, 0.0f, 1.0f),
        std::clamp(p.dyeB / vol, 0.0f, 1.0f));
    return vis;
}

bool needsDepth(WorldRenderStyle style) {
    return style == WorldRenderStyle::Detailed || style == WorldRenderStyle::Realistic;
}

bool needsHysteresis(WorldRenderStyle style) {
    return needsDepth(style) || style == WorldRenderStyle::Legacy;
}

bool needsNoise(WorldRenderStyle style) {
    return style == WorldRenderStyle::NoisyFlat
        || style == WorldRenderStyle::Realistic
        || style == WorldRenderStyle::Legacy;
}

float noiseScale(WorldLook const &look) {
    return std::clamp(look.noiseAmount, 0.0f, 2.0f);
}

float depthShade(int depth, float strength) {
    if (depth <= 0) return 1.0f + 0.16f * strength;
    float u = 1.0f - std::min(1.0f, static_cast<float>(depth) / 8.0f);
    return (0.70f + 0.30f * u) * (0.82f + 0.18f * (1.0f - strength)) + 0.08f * strength * u;
}

uint32_t shadeLegacy(MaterialVisual const &vis, WorldLook const &look, int x, int y, int depth,
    float amount, float speed, uint32_t salt, bool liquid, bool topSurface, bool solidEdge) {
    int r = vis.r, g = vis.g, b = vis.b;
    float n = vis.noiseStrength * noiseScale(look) * stableNoise(x, y, salt);
    if (!liquid) {
        float s = 1.0f + 0.18f * n;
        if (solidEdge) s *= 0.86f;
        s = std::clamp(s, 0.55f, 1.25f);
        return scaleRgb(r, g, b, s);
    }
    float amountLift = 0.48f + 0.52f * std::clamp(amount, 0.0f, 1.0f);
    float depthMul = 1.0f - 0.046f * static_cast<float>(std::min(depth, 10));
    depthMul = std::max(depthMul, 0.54f);
    float flow = std::tanh(speed / 24.0f);
    float motion = 1.0f + 0.13f * flow * vis.realisticDetail;
    float s = amountLift * depthMul * motion * (1.0f + 0.16f * n);
    if (topSurface) {
        s *= 1.24f;
        int lr = std::min(255, r + 72);
        int lg = std::min(255, g + 84);
        int lb = std::min(255, b + 48);
        tintToward(r, g, b, lr, lg, lb, 0.42f);
    }
    s = std::clamp(s, 0.42f, 1.70f);
    return scaleRgb(r, g, b, s);
}

uint32_t shadeMaterial(MaterialVisual const &vis, WorldLook const &look, int x, int y, int depth,
    float amount, float speed, uint32_t bg, uint32_t salt, bool liquid = false, bool topSurface = false,
    bool solidEdge = false) {
    if (look.style == WorldRenderStyle::Legacy)
        return shadeLegacy(vis, look, x, y, depth, amount, speed, salt, liquid, topSurface, solidEdge);
    int r = vis.r, g = vis.g, b = vis.b;
    float s = 1.0f;
    WorldRenderStyle const st = look.style;
    if (st == WorldRenderStyle::Detailed || st == WorldRenderStyle::Realistic)
        s *= depthShade(depth, vis.depthStrength);
    if (needsNoise(st))
        s *= 1.0f + vis.noiseStrength * 0.12f * noiseScale(look) * stableNoise(x, y, salt);
    if (st == WorldRenderStyle::Realistic) {
        float flow = std::tanh(speed / 28.0f);
        s *= 1.0f + vis.realisticDetail * 0.10f * flow;
        if (depth == 0) {
            float hi = vis.realisticDetail * (0.10f + 0.08f * flow);
            r = clampByte(r + static_cast<int>(18.0f * hi));
            g = clampByte(g + static_cast<int>(16.0f * hi));
            b = clampByte(b + static_cast<int>(12.0f * hi));
        }
    }
    uint32_t c = scaleRgb(r, g, b, s);
    if (st == WorldRenderStyle::AlphaFlat) {
        float a = vis.alphaResponse > 0.01f ? std::clamp(amount, 0.10f, 1.0f) * vis.alphaResponse
                                            : std::max(vis.alphaResponse, 1.0f);
        if (vis.alphaResponse <= 0.01f) a = 1.0f;
        int cr = static_cast<int>((c >> 16) & 255);
        int cg = static_cast<int>((c >> 8) & 255);
        int cb = static_cast<int>(c & 255);
        c = mixToward(bg, cr, cg, cb, a);
    }
    return c;
}

float wetDarkenMul(WorldRenderStyle style, float wetness) {
    wetness = std::clamp(wetness, 0.0f, 1.0f);
    if (wetness <= 1.0e-5f) return 1.0f;
    float sat = 0.20f;
    switch (style) {
        case WorldRenderStyle::FlatColor: sat = 0.13f; break;
        case WorldRenderStyle::NoisyFlat: sat = 0.20f; break;
        case WorldRenderStyle::Detailed: sat = 0.28f; break;
        case WorldRenderStyle::Realistic: sat = 0.34f; break;
        case WorldRenderStyle::AlphaFlat: sat = 0.10f; break;
        case WorldRenderStyle::Legacy: sat = 0.32f; break;
    }
    return 1.0f - sat * wetness;
}

void applyWetLook(int &r, int &g, int &b, WorldRenderStyle style, float wetness) {
    wetness = std::clamp(wetness, 0.0f, 1.0f);
    if (wetness <= 1.0e-5f) return;
    if (style == WorldRenderStyle::Realistic) {
        float mean = (static_cast<float>(r) + static_cast<float>(g) + static_cast<float>(b)) / 3.0f;
        float chroma = 1.0f + 0.18f * wetness;
        float contrast = 1.0f + 0.12f * wetness;
        auto push = [&](int ch) {
            float v = static_cast<float>(ch);
            v = mean + (v - mean) * chroma;
            v = mean + (v - mean) * contrast;
            return clampByte(static_cast<int>(std::lround(v)));
        };
        r = push(r); g = push(g); b = push(b);
    }
}

} // namespace

void WorldRenderer::ensureSize() {
    size_t n = static_cast<size_t>(GW * GH);
    if (visualLiquid.size() != n) {
        visualLiquid.assign(n, 0);
        liquidDepth.assign(n, -1);
        rigidDepth.assign(n, -1);
        glowDist.assign(n, 1.0e9f);
        glowR.assign(n, 0);
        glowG.assign(n, 0);
        glowB.assign(n, 0);
    }
}

void WorldRenderer::paintNormal(FluidEngine &fluid, RigidBodyEngine const &rigid, WorldLook const &look) {
    ensureSize();
    WorldRenderStyle const style = look.style;
    bool const useHysteresis = needsHysteresis(style);
    uint32_t const bg = packRgb(kVisualVoid.r, kVisualVoid.g, kVisualVoid.b);

    static std::vector<uint8_t> prevLiquid;
    if (prevLiquid.size() != visualLiquid.size()) prevLiquid.assign(visualLiquid.size(), 0);
    std::fill(visualLiquid.begin(), visualLiquid.end(), 0);
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = FluidEngine::ci(x, y);
        if (fluid.solid[i] || fluid.dynamicSolid[i]) continue;
        float amount = fluid.fill[static_cast<size_t>(i)];
        if (amount >= MIN_RENDER_FILL) visualLiquid[static_cast<size_t>(i)] = 1;
        else if (useHysteresis && amount >= MIN_ACTIVE_FILL && prevLiquid[static_cast<size_t>(i)])
            visualLiquid[static_cast<size_t>(i)] = 1;
    }
    if (useHysteresis) prevLiquid = visualLiquid;
    else std::fill(prevLiquid.begin(), prevLiquid.end(), 0);

    auto visibleLiquid = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return false;
        return visualLiquid[static_cast<size_t>(FluidEngine::ci(x, y))] != 0;
    };

    std::fill(liquidDepth.begin(), liquidDepth.end(), -1);
    if (needsDepth(style)) {
        workQueue.clear();
        if (workQueue.capacity() < static_cast<size_t>(GW * GH) / 4)
            workQueue.reserve(static_cast<size_t>(GW * GH) / 4);
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            if (!visibleLiquid(x, y)) continue;
            bool surface = false;
            for (int k = 0; k < 4; ++k) {
                if (!visibleLiquid(x + kN4x[k], y + kN4y[k])) { surface = true; break; }
            }
            if (!surface) continue;
            int i = FluidEngine::ci(x, y);
            liquidDepth[static_cast<size_t>(i)] = 0;
            workQueue.push_back(i);
        }
        for (size_t head = 0; head < workQueue.size(); ++head) {
            int i = workQueue[head];
            int x = i % GW, y = i / GW;
            int d = liquidDepth[static_cast<size_t>(i)];
            for (int k = 0; k < 4; ++k) {
                int nx = x + kN4x[k], ny = y + kN4y[k];
                if (!visibleLiquid(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (liquidDepth[static_cast<size_t>(ni)] >= 0) continue;
                liquidDepth[static_cast<size_t>(ni)] = d + 1;
                workQueue.push_back(ni);
            }
        }
    }

    auto occupantAt = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return -1;
        return rigid.occupant[static_cast<size_t>(FluidEngine::ci(x, y))];
    };
    auto sameBody = [&](int x, int y, int body) {
        return body >= 0 && occupantAt(x, y) == body;
    };

    std::fill(rigidDepth.begin(), rigidDepth.end(), -1);
    if (needsDepth(style)) {
        workQueue.clear();
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int body = occupantAt(x, y);
            if (body < 0) continue;
            bool edge = false;
            for (int k = 0; k < 4; ++k) {
                if (!sameBody(x + kN4x[k], y + kN4y[k], body)) { edge = true; break; }
            }
            if (!edge) continue;
            int i = FluidEngine::ci(x, y);
            rigidDepth[static_cast<size_t>(i)] = 0;
            workQueue.push_back(i);
        }
        for (size_t head = 0; head < workQueue.size(); ++head) {
            int i = workQueue[head];
            int x = i % GW, y = i / GW;
            int body = occupantAt(x, y);
            int d = rigidDepth[static_cast<size_t>(i)];
            for (int k = 0; k < 4; ++k) {
                int nx = x + kN4x[k], ny = y + kN4y[k];
                if (!sameBody(nx, ny, body)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (rigidDepth[static_cast<size_t>(ni)] >= 0) continue;
                rigidDepth[static_cast<size_t>(ni)] = d + 1;
                workQueue.push_back(ni);
            }
        }
    }

    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int index = FluidEngine::ci(x, y);
        uint32_t color = bg;
        if (fluid.solid[index]) {
            MaterialVisual vis = kVisualWall;
            int depth = 0;
            bool edge = false;
            if (needsDepth(style) || style == WorldRenderStyle::Legacy) {
                for (int k = 0; k < 4; ++k) {
                    int nx = x + kN4x[k], ny = y + kN4y[k];
                    if (!FluidEngine::inside(nx, ny) || !fluid.solid[FluidEngine::ci(nx, ny)]) {
                        edge = true;
                        break;
                    }
                }
                depth = edge ? 0 : 3;
            }
            color = shadeMaterial(vis, look, x, y, depth, 1.0f, 0.0f, bg, 17u, false, false, edge);
        } else if (visualLiquid[static_cast<size_t>(index)]) {
            MaterialVisual vis = sampleLiquidCell(fluid, index);
            int depth = 0;
            bool topSurface = false;
            if (needsDepth(style)) {
                depth = std::max(0, liquidDepth[static_cast<size_t>(index)]);
                topSurface = depth == 0;
            } else if (style == WorldRenderStyle::Legacy) {
                topSurface = !visibleLiquid(x, y - 1);
                for (int k = 1; k <= 10; ++k) {
                    if (!visibleLiquid(x, y - k)) break;
                    ++depth;
                }
            }
            float speed = 0.0f;
            if (style == WorldRenderStyle::Realistic || style == WorldRenderStyle::Legacy)
                speed = std::sqrt(fluid.cellU(x, y) * fluid.cellU(x, y) + fluid.cellV(x, y) * fluid.cellV(x, y));
            float amount = std::clamp(fluid.fill[static_cast<size_t>(index)], 0.0f, 1.0f);
            color = shadeMaterial(vis, look, x, y, depth, amount, speed, bg, 41u, true, topSurface, false);
            if (style != WorldRenderStyle::AlphaFlat)
                color = applyLiquidFillAlpha(color, bg, amount);
        }
        fluid.pixels[static_cast<size_t>(index)] = color;
    }

    if (look.glowingLiquids) {
        constexpr float glowRadius = 5.0f;
        std::fill(glowDist.begin(), glowDist.end(), 1.0e9f);
        workQueue.clear();
        auto seedGlow = [&](int x, int y, int r, int g, int b) {
            if (!FluidEngine::inside(x, y)) return;
            int i = FluidEngine::ci(x, y);
            if (fluid.solid[i] || fluid.dynamicSolid[i]) return;
            if (glowDist[static_cast<size_t>(i)] <= 0.0f) return;
            glowDist[static_cast<size_t>(i)] = 0.0f;
            glowR[static_cast<size_t>(i)] = static_cast<uint8_t>(clampByte(r));
            glowG[static_cast<size_t>(i)] = static_cast<uint8_t>(clampByte(g));
            glowB[static_cast<size_t>(i)] = static_cast<uint8_t>(clampByte(b));
            workQueue.push_back(i);
        };
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            if (!visibleLiquid(x, y)) continue;
            MaterialVisual vis = sampleLiquidCell(fluid, FluidEngine::ci(x, y));
            if (!vis.glowEligible) continue;
            seedGlow(x, y, vis.r + 28, vis.g + 32, vis.b + 12);
        }
        for (SplashParticle const &p : fluid.splashes) {
            MaterialVisual vis = sampleSplash(p);
            if (!vis.glowEligible) continue;
            seedGlow(static_cast<int>(p.x), static_cast<int>(p.y), vis.r + 28, vis.g + 32, vis.b + 12);
        }
        for (size_t head = 0; head < workQueue.size(); ++head) {
            int i = workQueue[head];
            int x = i % GW, y = i / GW;
            float d = glowDist[static_cast<size_t>(i)];
            uint8_t sr = glowR[static_cast<size_t>(i)];
            uint8_t sg = glowG[static_cast<size_t>(i)];
            uint8_t sb = glowB[static_cast<size_t>(i)];
            for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
                if (ox == 0 && oy == 0) continue;
                int nx = x + ox, ny = y + oy;
                if (!FluidEngine::inside(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (fluid.solid[ni] || fluid.dynamicSolid[ni]) continue;
                float nd = d + ((ox != 0 && oy != 0) ? 1.41421356f : 1.0f);
                if (nd > glowRadius || nd >= glowDist[static_cast<size_t>(ni)]) continue;
                glowDist[static_cast<size_t>(ni)] = nd;
                glowR[static_cast<size_t>(ni)] = sr;
                glowG[static_cast<size_t>(ni)] = sg;
                glowB[static_cast<size_t>(ni)] = sb;
                workQueue.push_back(ni);
            }
        }
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int index = FluidEngine::ci(x, y);
            float d = glowDist[static_cast<size_t>(index)];
            if (d > glowRadius || d <= 0.0f) continue;
            if (fluid.solid[index] || fluid.dynamicSolid[index]) continue;
            float u = 1.0f - d / glowRadius;
            fluid.pixels[static_cast<size_t>(index)] = mixToward(fluid.pixels[static_cast<size_t>(index)],
                glowR[static_cast<size_t>(index)], glowG[static_cast<size_t>(index)], glowB[static_cast<size_t>(index)],
                std::min(1.0f, u * u * (1.20f + 0.45f * u)));
        }
    }

    for (SplashParticle const &p : fluid.splashes) {
        int x = static_cast<int>(p.x), y = static_cast<int>(p.y);
        if (!FluidEngine::inside(x, y) || fluid.solid[FluidEngine::ci(x, y)] || fluid.dynamicSolid[FluidEngine::ci(x, y)])
            continue;
        MaterialVisual vis = sampleSplash(p);
        uint32_t dst = fluid.pixels[static_cast<size_t>(FluidEngine::ci(x, y))];
        fluid.pixels[static_cast<size_t>(FluidEngine::ci(x, y))] =
            mixToward(dst, vis.r, vis.g, vis.b, liquidVisualAlpha(p.volume));
    }

    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int index = FluidEngine::ci(x, y);
        MaterialId pendingMat = rigid.pending[static_cast<size_t>(index)];
        MaterialId occMat = rigid.worldCellMaterial(x, y);
        if (occMat != MATERIAL_EMPTY) {
            MaterialVisual vis = visualForSolid(occMat);
            float wetness = 0.0f;
            MaterialDefinition const &mat = materialDef(occMat);
            if (!fluid.solid[index] && mat.moistureCapacity > 1.0e-8f)
                wetness = std::clamp(rigid.occupantMoisture[static_cast<size_t>(index)] / mat.moistureCapacity, 0.0f, 1.0f);
            float dmg = std::clamp(rigid.occupantDamage[static_cast<size_t>(index)], 0.0f, 1.0f);
            float crack = 0.0f;
            if (index < static_cast<int>(rigid.occupantCrack.size()))
                crack = std::clamp(rigid.occupantCrack[static_cast<size_t>(index)], 0.0f, 1.0f);
            int bodyI = rigid.occupant[static_cast<size_t>(index)];
            if (bodyI >= 0 && bodyI < static_cast<int>(rigid.bodies.size())) {
                RigidBody const &rb = rigid.bodies[static_cast<size_t>(bodyI)];
                float lx, ly;
                RigidBodyEngine::worldToLocal(rb, x + 0.5f, y + 0.5f, lx, ly);
                int ix = static_cast<int>(std::floor(lx)), iy = static_cast<int>(std::floor(ly));
                for (int k = 0; k < 4; ++k) {
                    int nx = x + kN4x[k], ny = y + kN4y[k];
                    if (!FluidEngine::inside(nx, ny)) continue;
                    if (rigid.occupant[static_cast<size_t>(FluidEngine::ci(nx, ny))] != bodyI) continue;
                    float nlx, nly;
                    RigidBodyEngine::worldToLocal(rb, nx + 0.5f, ny + 0.5f, nlx, nly);
                    int nix = static_cast<int>(std::floor(nlx)), niy = static_cast<int>(std::floor(nly));
                    if (std::abs(nix - ix) + std::abs(niy - iy) != 1) continue;
                    if (!rigid.bondConnects(rb, ix, iy, nix, niy) &&
                        RigidBodyEngine::maskOccupied(rb, ix, iy) &&
                        RigidBodyEngine::maskOccupied(rb, nix, niy)) {
                        crack = 1.0f;
                        break;
                    }
                }
            }
            int depth = needsDepth(style) ? std::max(0, rigidDepth[static_cast<size_t>(index)]) : 0;
            bool edge = false;
            if (style == WorldRenderStyle::Legacy) {
                int body = occupantAt(x, y);
                edge = body < 0;
                if (body >= 0) {
                    for (int k = 0; k < 4; ++k) {
                        if (!sameBody(x + kN4x[k], y + kN4y[k], body)) { edge = true; break; }
                    }
                }
            }
            uint32_t c = shadeMaterial(vis, look, x, y, depth, 1.0f, 0.0f, bg,
                73u + static_cast<uint32_t>(occMat), false, false, edge);
            int cr = static_cast<int>((c >> 16) & 255);
            int cg = static_cast<int>((c >> 8) & 255);
            int cb = static_cast<int>(c & 255);
            applyWetLook(cr, cg, cb, style, wetness);
            float stain = wetDarkenMul(style, wetness) * (1.0f - 0.38f * dmg) * (1.0f - 0.62f * crack);
            fluid.pixels[static_cast<size_t>(index)] = scaleRgb(cr, cg, cb, stain);
        } else if (pendingMat != MATERIAL_EMPTY) {
            MaterialVisual vis = visualForSolid(pendingMat);
            fluid.pixels[static_cast<size_t>(index)] = packRgb(std::min(255, vis.r + 40),
                std::min(255, vis.g + 40), std::min(255, vis.b + 20));
        }
    }

    if (look.outlines) {
        auto isLiq = [&](int x, int y) {
            return FluidEngine::inside(x, y) && visualLiquid[static_cast<size_t>(FluidEngine::ci(x, y))] != 0;
        };
        auto isWall = [&](int x, int y) {
            return FluidEngine::inside(x, y) && fluid.solid[FluidEngine::ci(x, y)] != 0;
        };
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int index = FluidEngine::ci(x, y);
            bool edge = false;
            float dark = 0.30f;
            if (isLiq(x, y)) {
                dark = sampleLiquidCell(fluid, index).outlineDarken;
                for (int k = 0; k < 4; ++k) {
                    if (!isLiq(x + kN4x[k], y + kN4y[k])) { edge = true; break; }
                }
            } else if (isWall(x, y)) {
                dark = kVisualWall.outlineDarken;
                for (int k = 0; k < 4; ++k) {
                    if (!isWall(x + kN4x[k], y + kN4y[k])) { edge = true; break; }
                }
            } else {
                int body = occupantAt(x, y);
                if (body < 0) continue;
                dark = visualForSolid(rigid.worldCellMaterial(x, y)).outlineDarken;
                for (int k = 0; k < 4; ++k) {
                    if (!sameBody(x + kN4x[k], y + kN4y[k], body)) { edge = true; break; }
                }
            }
            if (!edge) continue;
            uint32_t c = fluid.pixels[static_cast<size_t>(index)];
            int cr = static_cast<int>((c >> 16) & 255);
            int cg = static_cast<int>((c >> 8) & 255);
            int cb = static_cast<int>(c & 255);
            float s = std::clamp(1.0f - dark, 0.35f, 1.0f);
            fluid.pixels[static_cast<size_t>(index)] = scaleRgb(cr, cg, cb, s);
        }
    }
}

void runWorldLookBenchmark(FluidEngine &fluid, RigidBodyEngine &rigid, WorldRenderer &renderer) {
    rigid.loadTestScene(fluid, 10);
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int i = FluidEngine::ci(x, y);
        if (fluid.fill[static_cast<size_t>(i)] < 0.2f) continue;
        if (x < 70) fluid.setLiquidComponentAmount(i, SUBSTANCE_HONEY, fluid.fill[static_cast<size_t>(i)]);
        else if (x < 110) {
            float f = fluid.fill[static_cast<size_t>(i)];
            fluid.dyeR[static_cast<size_t>(i)] = f * 0.88f;
            fluid.dyeG[static_cast<size_t>(i)] = f * 0.12f;
            fluid.dyeB[static_cast<size_t>(i)] = f * 0.16f;
        }
    }
    for (int t = 0; t < 8; ++t) {
        rigid.step(fluid, PHYSICS_DT);
        fluid.simulationTick();
        rigid.gatherFluidForces(fluid);
    }

    struct Case { char const *name; WorldLook look; };
    Case cases[] = {
        {"FlatColor", {WorldRenderStyle::FlatColor, false, false}},
        {"AlphaFlat", {WorldRenderStyle::AlphaFlat, false, false}},
        {"NoisyFlat", {WorldRenderStyle::NoisyFlat, false, false}},
        {"Detailed", {WorldRenderStyle::Detailed, false, false}},
        {"Realistic", {WorldRenderStyle::Realistic, false, false}},
        {"Legacy", {WorldRenderStyle::Legacy, false, false, 1.0f}},
        {"Realistic+Glow+Outlines", {WorldRenderStyle::Realistic, true, true, 1.0f}},
    };

    std::ofstream out(miscFile("look_bench.tsv"));
    out << "case\tms_per_paint\n";
    constexpr int kIters = 80;
    for (Case const &c : cases) {
        renderer.paintNormal(fluid, rigid, c.look);
        auto start = FluidEngine::Clock::now();
        for (int i = 0; i < kIters; ++i)
            renderer.paintNormal(fluid, rigid, c.look);
        double ms = FluidEngine::elapsedMs(start) / kIters;
        out << c.name << '\t' << ms << '\n';
    }
}
