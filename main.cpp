#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <windowsx.h>

#include "fluid/FluidEngine.h"
#include "fluid/DiagOutput.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "thermal/ThermalEngine.h"
#include "ui/UiShell.h"
#include "ui/UiTheme.h"
#include "ui/UiLanguage.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cwchar>
#include <fstream>
#include <string>
#include <vector>

namespace {

FluidEngine engine;
RigidBodyEngine rigid;
GasEngine gas;
ThermalEngine thermal;
ui::ShellState shell;

constexpr std::array<float, 6> speedScales{{0.10f, 0.25f, 0.50f, 1.0f, 2.0f, 4.0f}};
size_t speedScaleIndex = 3;
Tool activeTool = Tool::Water;
DebugView debugView = DebugView::Normal;
int brushRadius = 3;
bool paused = false;
bool mousePainting = false;
bool mouseGrabbing = false;
bool linePainting = false;
int lastPaintX = -1, lastPaintY = -1;
int lineStartX = -1, lineStartY = -1;
int lineEndX = -1, lineEndY = -1;
int hoverX = -1, hoverY = -1;
bool liquidGlow = false;
bool flatOutline = false;
bool settingsOpen = false;
std::vector<uint8_t> prevVisualLiquid(static_cast<size_t>(GW * GH), 0);
HWND mainWindow = nullptr;
HDC backDc = nullptr;
HBITMAP backBitmap = nullptr;
HGDIOBJ backOld = nullptr;
int backWidth = 0, backHeight = 0;

void releaseBackbuffer() {
    if (backDc) {
        if (backOld) SelectObject(backDc, backOld);
        if (backBitmap) DeleteObject(backBitmap);
        DeleteDC(backDc);
    }
    backDc = nullptr; backBitmap = nullptr; backOld = nullptr; backWidth = backHeight = 0;
}

void ensureBackbuffer(HDC dc, int width, int height) {
    width = std::max(1, width); height = std::max(1, height);
    if (backDc && backWidth == width && backHeight == height) return;
    releaseBackbuffer();
    backDc = CreateCompatibleDC(dc);
    backBitmap = CreateCompatibleBitmap(dc, width, height);
    backOld = SelectObject(backDc, backBitmap);
    backWidth = width; backHeight = height;
}

void changeTickRate(int direction) {
    int next = std::clamp(static_cast<int>(speedScaleIndex) + direction, 0, static_cast<int>(speedScales.size()) - 1);
    speedScaleIndex = static_cast<size_t>(next);
}

void refreshLayout(HWND hwnd) {
    RECT rc{}; GetClientRect(hwnd, &rc);
    ui::computeShellLayout(shell, rc.right, rc.bottom);
}

void clientToGrid(int mx, int my, int &gx, int &gy) {
    RECT const &c = shell.layout.canvasInner;
    int w = std::max(1, static_cast<int>(c.right - c.left));
    int h = std::max(1, static_cast<int>(c.bottom - c.top));
    gx = std::clamp(static_cast<int>((mx - c.left) * GW / w), 0, GW - 1);
    gy = std::clamp(static_cast<int>((my - c.top) * GH / h), 0, GH - 1);
}

void clientToGridF(int mx, int my, float &gx, float &gy) {
    RECT const &c = shell.layout.canvasInner;
    float w = std::max(1.0f, static_cast<float>(c.right - c.left));
    float h = std::max(1.0f, static_cast<float>(c.bottom - c.top));
    gx = std::clamp((static_cast<float>(mx - c.left) + 0.5f) * static_cast<float>(GW) / w, 0.0f, static_cast<float>(GW) - 0.001f);
    gy = std::clamp((static_cast<float>(my - c.top) + 0.5f) * static_cast<float>(GH) / h, 0.0f, static_cast<float>(GH) - 0.001f);
}

bool inCanvas(int x, int y) {
    return ui::ptIn(shell.layout.canvasInner, x, y);
}

ui::View makeView() {
    ui::View v;
    v.engine = &engine;
    v.rigid = &rigid;
    v.gas = &gas;
    v.tool = activeTool;
    v.debugView = debugView;
    v.brushRadius = brushRadius;
    v.paused = paused;
    v.settingsOpen = settingsOpen;
    v.liquidGlow = liquidGlow;
    v.flatOutline = flatOutline;
    v.speedIndex = speedScaleIndex;
    v.speedValue = speedScales[speedScaleIndex];
    v.hoverX = hoverX;
    v.hoverY = hoverY;
    v.heatPower = shell.heatPower;
    v.thermalEnabled = thermal.config.enabled;
    v.thermal = &thermal;
    return v;
}

void worldTick();
void invalidate();
void applySimQuality(int level, bool touchLook);
void setQualityPreset(QualityPreset preset);
void adaptAutoQuality(double worldMs);
void paintEnergyDisc(int cx, int cy);

void applySimQuality(int level, bool touchLook) {
    applyFluidQualityKnobs(engine.config, level);
    applyThermalQualityKnobs(thermal.config, level);
    engine.residualConsolidationEnabled = true;
    if (level <= 0) {
        gas.config.simMode = GasSimMode::Off;
        if (touchLook) {
            flatOutline = true;
            liquidGlow = false;
        }
    } else {
        gas.config.simMode = GasSimMode::Full;
        if (touchLook && level == 1) flatOutline = false;
    }
}

void setQualityPreset(QualityPreset preset) {
    engine.config.quality = preset;
    if (preset == QualityPreset::Auto) {
        engine.config.autoQualityLevel = 1;
        applySimQuality(1, true);
    } else if (preset == QualityPreset::Low) {
        applySimQuality(0, true);
    } else if (preset == QualityPreset::High) {
        applySimQuality(2, true);
    } else {
        applySimQuality(1, true);
    }
}

void adaptAutoQuality(double worldMs) {
    if (engine.config.quality != QualityPreset::Auto) return;
    static int slowStreak = 0;
    static int fastStreak = 0;
    if (worldMs > 24.0) { ++slowStreak; fastStreak = 0; }
    else if (worldMs < 8.0) { ++fastStreak; slowStreak = 0; }
    else { slowStreak = 0; fastStreak = 0; }
    int level = engine.config.autoQualityLevel;
    if (slowStreak >= 8 && level > 0) {
        engine.config.autoQualityLevel = 0;
        applySimQuality(0, false);
        slowStreak = 0;
    } else if (fastStreak >= 24 && level < 1) {
        engine.config.autoQualityLevel = 1;
        applySimQuality(1, false);
        fastStreak = 0;
    }
}

void handleMenuCommand(ui::MenuCmd cmd) {
    using ui::MenuCmd;
    switch (cmd) {
        case MenuCmd::None: break;
        case MenuCmd::Close:
            settingsOpen = false;
            shell.settingsDragging = false;
            shell.settingsScrollDragging = false;
            shell.settingsFlyout = ui::SettingsFlyout::None;
            break;
        case MenuCmd::Pause: paused = !paused; break;
        case MenuCmd::Step: if (paused) worldTick(); break;
        case MenuCmd::Clear: rigid.clear(); engine.clearWorld(); gas.resetAmbient(engine); thermal.seedAmbient(engine, rigid, gas); shell.log(ui::tr("log_cleared")); break;
        case MenuCmd::Reset: rigid.clear(); engine.resetWorld(); gas.resetAmbient(engine); thermal.seedAmbient(engine, rigid, gas); shell.log(ui::tr("log_reset")); break;
        case MenuCmd::Slosh: engine.addSloshImpulse(); break;
        case MenuCmd::WalledBorders:
            engine.config.walledBorders = !engine.config.walledBorders;
            gas.config.boundary = engine.config.walledBorders ? GasBoundary::Sealed : GasBoundary::OpenAmbient;
            break;
        case MenuCmd::ViewNormal: debugView = DebugView::Normal; break;
        case MenuCmd::ViewFill: debugView = DebugView::Fill; break;
        case MenuCmd::ViewPressure: debugView = DebugView::Pressure; break;
        case MenuCmd::ViewVelocity: debugView = DebugView::Velocity; break;
        case MenuCmd::ViewDivergence: debugView = DebugView::Divergence; break;
        case MenuCmd::ViewChunks: debugView = DebugView::Chunks; break;
        case MenuCmd::ViewRigid: debugView = DebugView::Rigid; break;
        case MenuCmd::ViewGasPressure: debugView = DebugView::GasPressure; break;
        case MenuCmd::ViewGasAmount: debugView = DebugView::GasAmount; break;
        case MenuCmd::ViewGasVelocity: debugView = DebugView::GasVelocity; break;
        case MenuCmd::ThermalOn: thermal.config.enabled = true; break;
        case MenuCmd::ThermalOff: thermal.config.enabled = false; break;
        case MenuCmd::Glow: liquidGlow = !liquidGlow; break;
        case MenuCmd::FlatOutline: flatOutline = !flatOutline; break;
        case MenuCmd::Overlay: rigid.debugOverlay = !rigid.debugOverlay; break;
        case MenuCmd::MatWood: shell.applyPalette(ui::PaletteId::Wood, activeTool, rigid.drawMaterial); break;
        case MenuCmd::MatStone: shell.applyPalette(ui::PaletteId::Stone, activeTool, rigid.drawMaterial); break;
        case MenuCmd::MatGlass: shell.applyPalette(ui::PaletteId::Glass, activeTool, rigid.drawMaterial); break;
        case MenuCmd::MatMetal: shell.applyPalette(ui::PaletteId::Metal, activeTool, rigid.drawMaterial); break;
        case MenuCmd::Vorticity: engine.config.vorticityEnabled = !engine.config.vorticityEnabled; break;
        case MenuCmd::Advection:
            engine.config.velocityAdvection = nextVelocityAdvection(engine.config.velocityAdvection);
            break;
        case MenuCmd::AdvectNone: engine.config.velocityAdvection = VelocityAdvection::None; break;
        case MenuCmd::AdvectFou: engine.config.velocityAdvection = VelocityAdvection::FirstOrderUpwind; break;
        case MenuCmd::AdvectNsl: engine.config.velocityAdvection = VelocityAdvection::NearestSemiLagrangian; break;
        case MenuCmd::AdvectSl: engine.config.velocityAdvection = VelocityAdvection::SemiLagrangian; break;
        case MenuCmd::AdvectMacc: engine.config.velocityAdvection = VelocityAdvection::MacCormack; break;
        case MenuCmd::AdvectBfecc: engine.config.velocityAdvection = VelocityAdvection::BFECC; break;
        case MenuCmd::VorticityMinus:
            engine.config.vorticityStrength = std::max(0.005f, engine.config.vorticityStrength - 0.005f);
            break;
        case MenuCmd::VorticityPlus:
            engine.config.vorticityStrength = std::min(0.12f, engine.config.vorticityStrength + 0.005f);
            break;
        case MenuCmd::CatchUpMinus:
            engine.config.catchUpTicks = std::max(1, engine.config.catchUpTicks - 1);
            break;
        case MenuCmd::CatchUpPlus:
            engine.config.catchUpTicks = std::min(4, engine.config.catchUpTicks + 1);
            break;
        case MenuCmd::ThermalIntervalMinus:
            thermal.config.intervalTicks = std::max(1, thermal.config.intervalTicks - 1);
            break;
        case MenuCmd::ThermalIntervalPlus:
            thermal.config.intervalTicks = std::min(6, thermal.config.intervalTicks + 1);
            break;
        case MenuCmd::RigidGravity:
            rigid.gravityScale = (rigid.gravityScale > 0.5f) ? 0.0f : 1.0f;
            break;
        case MenuCmd::Residual: engine.residualConsolidationEnabled = !engine.residualConsolidationEnabled; break;
        case MenuCmd::QualityLow: setQualityPreset(QualityPreset::Low); break;
        case MenuCmd::QualityMed: setQualityPreset(QualityPreset::Medium); break;
        case MenuCmd::QualityHigh: setQualityPreset(QualityPreset::High); break;
        case MenuCmd::QualityAuto: setQualityPreset(QualityPreset::Auto); break;
        case MenuCmd::Tension: engine.config.surfaceTensionEnabled = !engine.config.surfaceTensionEnabled; break;
        case MenuCmd::Spray: engine.config.sprayEnabled = !engine.config.sprayEnabled; break;
        case MenuCmd::SubstepsMinus: engine.config.maxSubsteps = std::max(1, engine.config.maxSubsteps - 1); break;
        case MenuCmd::SubstepsPlus: engine.config.maxSubsteps = std::min(6, engine.config.maxSubsteps + 1); break;
        case MenuCmd::LimiterMinus: engine.config.maxLimiterPasses = std::max(4, engine.config.maxLimiterPasses - 4); break;
        case MenuCmd::LimiterPlus: engine.config.maxLimiterPasses = std::min(16, engine.config.maxLimiterPasses + 4); break;
        case MenuCmd::GasOff: gas.config.simMode = GasSimMode::Off; break;
        case MenuCmd::GasHalf: gas.config.simMode = GasSimMode::Half; break;
        case MenuCmd::GasFull: gas.config.simMode = GasSimMode::Full; break;
        case MenuCmd::Hz20: engine.config.physicsHz = 20; engine.config.catchUpTicks = 1; break;
        case MenuCmd::Hz30: engine.config.physicsHz = 30; engine.config.catchUpTicks = 2; break;
        case MenuCmd::ThreadsAuto: engine.config.workerCount = 0; engine.syncWorkerPool(); break;
        case MenuCmd::Threads1: engine.config.workerCount = 1; engine.syncWorkerPool(); break;
        case MenuCmd::Threads2: engine.config.workerCount = 2; engine.syncWorkerPool(); break;
        case MenuCmd::Threads4: engine.config.workerCount = 4; engine.syncWorkerPool(); break;
        case MenuCmd::Threads6: engine.config.workerCount = 6; engine.syncWorkerPool(); break;
        case MenuCmd::Threads8: engine.config.workerCount = 8; engine.syncWorkerPool(); break;
        case MenuCmd::PressureMinus: engine.config.maxPressureIterations = std::max(8, engine.config.maxPressureIterations - 2); break;
        case MenuCmd::PressurePlus: engine.config.maxPressureIterations = std::min(30, engine.config.maxPressureIterations + 2); break;
        case MenuCmd::BrushMinus: brushRadius = std::max(1, brushRadius - 1); break;
        case MenuCmd::BrushPlus: brushRadius = std::min(14, brushRadius + 1); break;
        case MenuCmd::Speed0: case MenuCmd::Speed1: case MenuCmd::Speed2:
        case MenuCmd::Speed3: case MenuCmd::Speed4: case MenuCmd::Speed5:
            speedScaleIndex = static_cast<size_t>(static_cast<int>(cmd) - static_cast<int>(MenuCmd::Speed0));
            break;
        case MenuCmd::Scene1: case MenuCmd::Scene2: case MenuCmd::Scene3: case MenuCmd::Scene4: case MenuCmd::Scene5:
        case MenuCmd::Scene6: case MenuCmd::Scene7: case MenuCmd::Scene8: case MenuCmd::Scene9: case MenuCmd::Scene10:
            rigid.clear();
            engine.loadTestScene(static_cast<int>(cmd) - static_cast<int>(MenuCmd::Scene1) + 1);
            gas.resetAmbient(engine);
            thermal.seedAmbient(engine, rigid, gas);
            shell.log(ui::tr("log_fluid_scene"));
            break;
        case MenuCmd::Rigid1: case MenuCmd::Rigid2: case MenuCmd::Rigid3: case MenuCmd::Rigid4: case MenuCmd::Rigid5:
        case MenuCmd::Rigid6: case MenuCmd::Rigid7: case MenuCmd::Rigid8: case MenuCmd::Rigid9: case MenuCmd::Rigid10:
        case MenuCmd::Rigid11: case MenuCmd::Rigid12: case MenuCmd::Rigid13: case MenuCmd::Rigid14: case MenuCmd::Rigid15:
        case MenuCmd::Rigid16: case MenuCmd::Rigid17: case MenuCmd::Rigid18: case MenuCmd::Rigid19: case MenuCmd::Rigid20:
        case MenuCmd::Rigid21: case MenuCmd::Rigid22: case MenuCmd::Rigid23: case MenuCmd::Rigid24:
            rigid.gravityScale = (static_cast<int>(cmd) - static_cast<int>(MenuCmd::Rigid1) + 1 == 24) ? 0.0f : 1.0f;
            rigid.loadTestScene(engine, static_cast<int>(cmd) - static_cast<int>(MenuCmd::Rigid1) + 1);
            gas.resetAmbient(engine);
            thermal.seedAmbient(engine, rigid, gas);
            shell.log(ui::tr("log_rigid_scene"));
            break;
        case MenuCmd::Gas1: case MenuCmd::Gas2: case MenuCmd::Gas3: case MenuCmd::Gas4:
        case MenuCmd::Gas5: case MenuCmd::Gas6: case MenuCmd::Gas7: case MenuCmd::Gas8:
            if (gas.config.simMode == GasSimMode::Off) gas.config.simMode = GasSimMode::Full;
            gas.loadTestScene(engine, rigid, static_cast<int>(cmd) - static_cast<int>(MenuCmd::Gas1) + 1);
            thermal.seedAmbient(engine, rigid, gas);
            shell.log(ui::tr("log_gas_scene"));
            break;
    }
}

void handleHit(ui::HitId id) {
    using ui::HitId;
    switch (id) {
        case HitId::SpeedMinus: changeTickRate(-1); break;
        case HitId::SpeedPlus: changeTickRate(1); break;
        case HitId::Settings:
            settingsOpen = !settingsOpen;
            if (settingsOpen) ui::ensureSettingsPlacement(shell);
            else {
                shell.settingsDragging = false;
                shell.settingsScrollDragging = false;
                shell.settingsFlyout = ui::SettingsFlyout::None;
            }
            break;
        case HitId::Pause: paused = !paused; break;
        case HitId::Clear: handleMenuCommand(ui::MenuCmd::Clear); break;
        case HitId::ViewNorm: debugView = DebugView::Normal; break;
        case HitId::ViewChnk: debugView = DebugView::Chunks; break;
        case HitId::ViewFill: debugView = DebugView::Fill; break;
        case HitId::ViewLiqp: debugView = DebugView::Pressure; break;
        case HitId::ViewLvel: debugView = DebugView::Velocity; break;
        case HitId::ViewLdiv: debugView = DebugView::Divergence; break;
        case HitId::ViewRgdn: debugView = DebugView::Rigid; break;
        case HitId::ViewRgdo: rigid.debugOverlay = !rigid.debugOverlay; break;
        case HitId::ViewGasp: debugView = DebugView::GasPressure; break;
        case HitId::ViewGasa: debugView = DebugView::GasAmount; break;
        case HitId::ViewGasv: debugView = DebugView::GasVelocity; break;
        case HitId::ViewTemp: debugView = DebugView::Temperature; break;
        case HitId::Dye0: case HitId::Dye1: case HitId::Dye2: case HitId::Dye3:
        case HitId::Dye4: case HitId::Dye5: case HitId::Dye6: case HitId::Dye7:
            shell.dyeSwatch = static_cast<int>(id) - static_cast<int>(HitId::Dye0);
            break;
        case HitId::DyeOnly: shell.dyeOnly = !shell.dyeOnly; break;
        case HitId::DyeStrMinus: shell.dyeStrength = std::max(0.15f, shell.dyeStrength - 0.15f); break;
        case HitId::DyeStrPlus: shell.dyeStrength = std::min(1.0f, shell.dyeStrength + 0.15f); break;
        case HitId::Pal0: case HitId::Pal1: case HitId::Pal2: case HitId::Pal3:
        case HitId::Pal4: case HitId::Pal5: case HitId::Pal6: case HitId::Pal7: {
            int slot = static_cast<int>(id) - static_cast<int>(HitId::Pal0);
            if (slot < shell.elementCount()) shell.applyPalette(shell.elementAt(slot), activeTool, rigid.drawMaterial);
            break;
        }
        case HitId::Cat0: case HitId::Cat1: case HitId::Cat2: case HitId::Cat3:
        case HitId::Cat4: case HitId::Cat5: case HitId::Cat6: {
            int slot = static_cast<int>(id) - static_cast<int>(HitId::Cat0);
            shell.applyCategory(static_cast<ui::Category>(slot), activeTool, rigid.drawMaterial);
            break;
        }
        case HitId::BrushMinus: brushRadius = std::max(1, brushRadius - 1); break;
        case HitId::BrushPlus: brushRadius = std::min(14, brushRadius + 1); break;
        case HitId::PowerMinus: shell.heatPower = std::max(0.25f, shell.heatPower - 0.25f); break;
        case HitId::PowerPlus: shell.heatPower = std::min(8.0f, shell.heatPower + 0.25f); break;
        case HitId::Anchored:
            if (shell.category == ui::Category::Solids) {
                shell.placeAnchored = !shell.placeAnchored;
                if (shell.placeAnchored) {
                    shell.placeSleeping = false;
                    shell.placePowder = false;
                }
            }
            break;
        case HitId::Sleeping:
            if (shell.category == ui::Category::Solids) {
                shell.placeSleeping = !shell.placeSleeping;
                if (shell.placeSleeping) {
                    shell.placeAnchored = false;
                    shell.placePowder = false;
                }
            }
            break;
        case HitId::Powder:
            if (shell.category == ui::Category::Solids) {
                shell.placePowder = !shell.placePowder;
                if (shell.placePowder) {
                    shell.placeAnchored = false;
                    shell.placeSleeping = false;
                }
            }
            break;
        case HitId::GrainMinus:
            if (shell.category == ui::Category::Solids)
                shell.powderParticleSize = std::max(1, shell.powderParticleSize - 1);
            break;
        case HitId::GrainPlus:
            if (shell.category == ui::Category::Solids)
                shell.powderParticleSize = std::min(8, shell.powderParticleSize + 1);
            break;
        case HitId::Search:
            shell.searchFocused = true;
            shell.consoleFocused = false;
            break;
        case HitId::Console:
            shell.consoleFocused = true;
            shell.searchFocused = false;
            break;
        default: break;
    }
}

uint32_t rgb(int r, int g, int b) {
    return static_cast<uint32_t>(std::clamp(b, 0, 255)) | (static_cast<uint32_t>(std::clamp(g, 0, 255)) << 8)
        | (static_cast<uint32_t>(std::clamp(r, 0, 255)) << 16);
}

constexpr int kWaterR = 22, kWaterG = 126, kWaterB = 214;
constexpr int kWaterDarkR = 14, kWaterDarkG = 86, kWaterDarkB = 164;
constexpr int kWaterLightR = 54, kWaterLightG = 164, kWaterLightB = 226;
constexpr int kWaterRimR = 72, kWaterRimG = 176, kWaterRimB = 230;
constexpr int kHoneyR = 176, kHoneyG = 110, kHoneyB = 22;
constexpr int kHoneyDarkR = 120, kHoneyDarkG = 70, kHoneyDarkB = 12;
constexpr int kHoneyLightR = 210, kHoneyLightG = 150, kHoneyLightB = 40;
constexpr int kHoneyRimR = 198, kHoneyRimG = 140, kHoneyRimB = 48;

void tintChannels(int &r, int &g, int &b, int r1, int g1, int b1, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    r = static_cast<int>(std::lround(r + (r1 - r) * t));
    g = static_cast<int>(std::lround(g + (g1 - g) * t));
    b = static_cast<int>(std::lround(b + (b1 - b) * t));
}

void liquidAppearance(FluidEngine const &eng, int index, int wr, int wg, int wb, int &r, int &g, int &b) {
    r = wr; g = wg; b = wb;
    float f = std::max(eng.fill[static_cast<size_t>(index)], 1.0e-8f);
    float h = eng.honeyFraction(index);
    if (h > 0.001f) {
        int hr = kHoneyR, hg = kHoneyG, hb = kHoneyB;
        if (wr == kWaterRimR) { hr = kHoneyRimR; hg = kHoneyRimG; hb = kHoneyRimB; }
        else if (wr == kWaterDarkR) { hr = kHoneyDarkR; hg = kHoneyDarkG; hb = kHoneyDarkB; }
        else if (wr == kWaterLightR) { hr = kHoneyLightR; hg = kHoneyLightG; hb = kHoneyLightB; }
        tintChannels(r, g, b, hr, hg, hb, h);
    }
    float ir = std::clamp(eng.dyeR[static_cast<size_t>(index)] / f, 0.0f, 1.0f);
    float ig = std::clamp(eng.dyeG[static_cast<size_t>(index)] / f, 0.0f, 1.0f);
    float ib = std::clamp(eng.dyeB[static_cast<size_t>(index)] / f, 0.0f, 1.0f);
    float dyeI = std::clamp(std::max(ir, std::max(ig, ib)), 0.0f, 1.0f);
    if (dyeI > 0.02f)
        tintChannels(r, g, b, static_cast<int>(ir * 255.0f), static_cast<int>(ig * 255.0f), static_cast<int>(ib * 255.0f), dyeI);
}

uint32_t scaleRgb(int r, int g, int b, float s) {
    return rgb(static_cast<int>(std::lround(static_cast<float>(r) * s)),
        static_cast<int>(std::lround(static_cast<float>(g) * s)),
        static_cast<int>(std::lround(static_cast<float>(b) * s)));
}

uint32_t lerpRgb(int r0, int g0, int b0, int r1, int g1, int b1, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    return rgb(static_cast<int>(std::lround(r0 + (r1 - r0) * t)),
        static_cast<int>(std::lround(g0 + (g1 - g0) * t)),
        static_cast<int>(std::lround(b0 + (b1 - b0) * t)));
}

uint32_t mixToward(uint32_t dst, int r, int g, int b, float t);

template <typename Fn>
void forEachBrushCell(int cx, int cy, Fn &&fn) {
    int r2 = brushRadius * brushRadius;
    for (int y = cy - brushRadius; y <= cy + brushRadius; ++y)
        for (int x = cx - brushRadius; x <= cx + brushRadius; ++x) {
            if (!FluidEngine::inside(x, y)) continue;
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) > r2) continue;
            fn(x, y);
        }
}

template <typename Fn>
void forEachBrushLine(int x0, int y0, int x1, int y1, Fn &&fn) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;
    for (;;) {
        forEachBrushCell(x0, y0, fn);
        if (x0 == x1 && y0 == y1) break;
        int twice = 2 * error;
        if (twice >= dy) { error += dy; x0 += sx; }
        if (twice <= dx) { error += dx; y0 += sy; }
    }
}

void ghostTint(int &r, int &g, int &b) {
    r = 214; g = 208; b = 190;
    if (activeTool == Tool::Water) {
        r = 70; g = 168; b = 214;
        if (shell.palette == ui::PaletteId::Honey) { r = 210; g = 150; b = 48; }
        if (shell.dyeOnly || shell.dyeSwatch > 0) {
            float dr, dg, db;
            shell.dyeChannels(dr, dg, db);
            if (dr + dg + db > 1.0e-6f) {
                r = static_cast<int>(dr * 255.0f);
                g = static_cast<int>(dg * 255.0f);
                b = static_cast<int>(db * 255.0f);
            }
        }
    }
    else if (activeTool == Tool::Solid) { r = 92; g = 94; b = 100; }
    else if (activeTool == Tool::Eraser) { r = 196; g = 92; b = 72; }
    else if (activeTool == Tool::Heat) { r = 210; g = 110; b = 60; }
    else if (activeTool == Tool::Cool) { r = 80; g = 130; b = 200; }
    else if (activeTool == Tool::Pressurize) { r = 200; g = 160; b = 70; }
    else if (activeTool == Tool::Depressurize) { r = 110; g = 112; b = 140; }
    else if (activeTool == Tool::Rigid) {
        MaterialDefinition const &mat = materialDef(rigid.drawMaterial);
        r = std::min(255, mat.colorR + 36);
        g = std::min(255, mat.colorG + 36);
        b = std::min(255, mat.colorB + 24);
    }
}

void overlayLineGhost() {
    if (!linePainting) return;
    int gr, gg, gb;
    ghostTint(gr, gg, gb);
    static std::vector<uint8_t> mark;
    mark.assign(static_cast<size_t>(GW * GH), 0);
    forEachBrushLine(lineStartX, lineStartY, lineEndX, lineEndY, [&](int x, int y) {
        size_t i = static_cast<size_t>(FluidEngine::ci(x, y));
        if (mark[i]) return;
        mark[i] = 1;
        engine.pixels[i] = mixToward(engine.pixels[i], gr, gg, gb, 0.46f);
    });
}

void gridToClientCenter(int gx, int gy, int &sx, int &sy) {
    RECT const &c = shell.layout.canvasInner;
    float w = std::max(1.0f, static_cast<float>(c.right - c.left));
    float h = std::max(1.0f, static_cast<float>(c.bottom - c.top));
    sx = c.left + static_cast<int>(std::lround((static_cast<float>(gx) + 0.5f) * w / static_cast<float>(GW)));
    sy = c.top + static_cast<int>(std::lround((static_cast<float>(gy) + 0.5f) * h / static_cast<float>(GH)));
}

void drawRing(HDC dc, int cx, int cy, int rx, int ry) {
    auto stroke = [&](int prx, int pry, COLORREF color) {
        if (prx < 1 || pry < 1) return;
        HPEN pen = CreatePen(PS_SOLID, 1, color);
        HGDIOBJ oldPen = SelectObject(dc, pen);
        HGDIOBJ oldBr = SelectObject(dc, GetStockObject(HOLLOW_BRUSH));
        Ellipse(dc, cx - prx, cy - pry, cx + prx + 1, cy + pry + 1);
        SelectObject(dc, oldPen);
        SelectObject(dc, oldBr);
        DeleteObject(pen);
    };
    stroke(rx, ry, RGB(28, 26, 22));
    if (rx > 2 && ry > 2) stroke(rx - 1, ry - 1, RGB(236, 228, 210));
}

void drawBrushOverlay(HDC dc) {
    if (settingsOpen || activeTool == Tool::Grab) return;
    if (!shell.hasPlacement() && !linePainting) return;
    bool over = inCanvas(shell.mouseX, shell.mouseY);
    if (!over && !linePainting) return;
    RECT const &c = shell.layout.canvasInner;
    float cellW = std::max(1.0f, static_cast<float>(c.right - c.left) / static_cast<float>(GW));
    float cellH = std::max(1.0f, static_cast<float>(c.bottom - c.top) / static_cast<float>(GH));
    int rx = std::max(2, static_cast<int>(std::lround(static_cast<float>(brushRadius) * cellW)));
    int ry = std::max(2, static_cast<int>(std::lround(static_cast<float>(brushRadius) * cellH)));
    HRGN clip = CreateRectRgn(c.left, c.top, c.right, c.bottom);
    SelectClipRgn(dc, clip);
    if (over) drawRing(dc, shell.mouseX, shell.mouseY, rx, ry);
    if (linePainting) {
        int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
        gridToClientCenter(lineStartX, lineStartY, x0, y0);
        gridToClientCenter(lineEndX, lineEndY, x1, y1);
        drawRing(dc, x0, y0, rx, ry);
        HPEN axis = CreatePen(PS_SOLID, 1, RGB(28, 26, 22));
        HGDIOBJ oldPen = SelectObject(dc, axis);
        MoveToEx(dc, x0, y0, nullptr);
        LineTo(dc, x1, y1);
        HPEN axisLite = CreatePen(PS_SOLID, 1, RGB(236, 228, 210));
        SelectObject(dc, axisLite);
        MoveToEx(dc, x0, y0 - 1, nullptr);
        LineTo(dc, x1, y1 - 1);
        SelectObject(dc, oldPen);
        DeleteObject(axis);
        DeleteObject(axisLite);
        if (!over) drawRing(dc, x1, y1, rx, ry);
    }
    SelectClipRgn(dc, nullptr);
    DeleteObject(clip);
}

uint32_t mixToward(uint32_t dst, int r, int g, int b, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    int dr = static_cast<int>((dst >> 16) & 255);
    int dg = static_cast<int>((dst >> 8) & 255);
    int db = static_cast<int>(dst & 255);
    return lerpRgb(dr, dg, db, r, g, b, t);
}

uint32_t fillDebugColor(float amount) {
    float f = std::clamp(amount, 0.0f, 1.0f); if (f <= 0.0f) return rgb(0,0,0);
    if (f < 0.25f) return rgb(0, static_cast<int>(40 + f * 240), static_cast<int>(85 + f * 500));
    if (f < 0.50f) return rgb(0, static_cast<int>(100 + (f-0.25f)*400), 220);
    if (f < 0.75f) return rgb(static_cast<int>((f-0.5f)*320), 210, 245);
    return rgb(static_cast<int>(80 + (f-0.75f)*700), static_cast<int>(210 + (f-0.75f)*180), 255);
}

uint32_t gasPressureColor(float atm) {
    if (atm <= 1.0e-4f) return rgb(6, 6, 8);
    if (atm < 1.0f) {
        float t = std::clamp(atm, 0.0f, 1.0f);
        return lerpRgb(12, 10, 18, 72, 92, 108, t);
    }
    float t = std::clamp((atm - 1.0f) / 1.5f, 0.0f, 1.0f);
    return lerpRgb(72, 92, 108, 245, 220, 70, t);
}

uint32_t gasAmountColor(float amount) {
    float t = std::clamp(amount / 2.0f, 0.0f, 1.0f);
    if (t <= 0.0f) return rgb(6, 6, 8);
    return lerpRgb(20, 28, 36, 90, 200, 170, t);
}

void fillWorldPixels() {
    float pressureScale = 0.0f, divergenceScale = 0.0f;
    bool const flatPaint = flatOutline && debugView == DebugView::Normal;
    bool const needFieldScale = debugView == DebugView::Pressure || debugView == DebugView::Divergence;
    constexpr int n4x[4] = {-1, 1, 0, 0};
    constexpr int n4y[4] = {0, 0, -1, 1};
    if (needFieldScale) {
        for (int i = 0; i < GW * GH; ++i) {
            pressureScale = std::max(pressureScale, std::abs(engine.pressure[i]));
            divergenceScale = std::max(divergenceScale, std::abs(engine.divergenceField[i]));
        }
        pressureScale = std::max(pressureScale, 1e-4f); divergenceScale = std::max(divergenceScale, 1e-4f);
    }
    static std::vector<uint8_t> visualLiquid(static_cast<size_t>(GW * GH), 0);
    static std::vector<int> liquidDepth(static_cast<size_t>(GW * GH), -1);
    static std::vector<int> depthQueue;
    if (debugView == DebugView::Normal) {
        std::fill(visualLiquid.begin(), visualLiquid.end(), 0);
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int i = FluidEngine::ci(x, y);
            if (engine.solid[i] || engine.dynamicSolid[i]) continue;
            float amount = engine.fill[i];
            if (amount >= MIN_RENDER_FILL) visualLiquid[static_cast<size_t>(i)] = 1;
            else if (!flatPaint && amount >= MIN_ACTIVE_FILL && prevVisualLiquid[static_cast<size_t>(i)]) visualLiquid[static_cast<size_t>(i)] = 1;
        }
        prevVisualLiquid = visualLiquid;
    }
    auto visibleLiquid = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return false;
        int i = FluidEngine::ci(x, y);
        if (debugView == DebugView::Normal) return visualLiquid[static_cast<size_t>(i)] != 0;
        return !engine.solid[i] && !engine.dynamicSolid[i] && engine.fill[i] >= MIN_RENDER_FILL;
    };
    std::fill(liquidDepth.begin(), liquidDepth.end(), -1);
    if (debugView == DebugView::Normal && !flatPaint) {
        depthQueue.clear();
        if (depthQueue.capacity() < static_cast<size_t>(GW * GH) / 4)
            depthQueue.reserve(static_cast<size_t>(GW * GH) / 4);
        constexpr int ox[4] = {-1, 1, 0, 0};
        constexpr int oy[4] = {0, 0, -1, 1};
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            if (!visibleLiquid(x, y)) continue;
            bool surface = false;
            for (int k = 0; k < 4; ++k) {
                if (!visibleLiquid(x + ox[k], y + oy[k])) { surface = true; break; }
            }
            if (!surface) continue;
            int i = FluidEngine::ci(x, y);
            liquidDepth[static_cast<size_t>(i)] = 0;
            depthQueue.push_back(i);
        }
        for (size_t head = 0; head < depthQueue.size(); ++head) {
            int i = depthQueue[head];
            int x = i % GW, y = i / GW;
            int d = liquidDepth[static_cast<size_t>(i)];
            for (int k = 0; k < 4; ++k) {
                int nx = x + ox[k], ny = y + oy[k];
                if (!visibleLiquid(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (liquidDepth[static_cast<size_t>(ni)] >= 0) continue;
                liquidDepth[static_cast<size_t>(ni)] = d + 1;
                depthQueue.push_back(ni);
            }
        }
    }
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int index = FluidEngine::ci(x, y); uint32_t color = rgb(10, 17, 28);
        if (engine.solid[index]) {
            if (flatPaint) {
                bool edge = false;
                for (int k = 0; k < 4; ++k) {
                    int nx = x + n4x[k], ny = y + n4y[k];
                    if (!FluidEngine::inside(nx, ny) || !engine.solid[FluidEngine::ci(nx, ny)]) { edge = true; break; }
                }
                color = edge ? rgb(58, 62, 70) : rgb(96, 100, 108);
            } else {
                int shade = 88 + static_cast<int>(FluidEngine::hashCell(x, y, 1) & 15u);
                color = rgb(shade, shade + 5, shade + 12);
            }
        }
        else if (debugView == DebugView::GasPressure && !engine.dynamicSolid[index]) {
            color = gasPressureColor(gas.pressure[static_cast<size_t>(index)]);
        } else if (debugView == DebugView::GasAmount && !engine.dynamicSolid[index]) {
            color = gasAmountColor(gas.amount[static_cast<size_t>(index)]);
        } else if (debugView == DebugView::GasVelocity && !engine.dynamicSolid[index] && gas.volume[static_cast<size_t>(index)] >= GAS_MIN_VOLUME) {
            float speed = std::sqrt(gas.cellU(x, y) * gas.cellU(x, y) + gas.cellV(x, y) * gas.cellV(x, y));
            float q = std::clamp(speed / 20.0f, 0.0f, 1.0f);
            color = rgb(int(40 + 215 * q), int(80 + 100 * (1.0f - q)), int(255 * (1.0f - q)));
        } else if (debugView == DebugView::Pressure && engine.fill[index] >= MIN_RENDER_FILL) {
            float q = std::clamp(engine.pressure[index] / pressureScale, -1.0f, 1.0f); color = q >= 0 ? rgb(35 + int(q * 210), 35, 55) : rgb(35, 55, 35 + int(-q * 210));
        } else if (debugView == DebugView::Velocity && engine.fill[index] >= MIN_RENDER_FILL) {
            float speed = std::sqrt(engine.cellU(x, y) * engine.cellU(x, y) + engine.cellV(x, y) * engine.cellV(x, y)); float q = std::clamp(speed / 30.0f, 0.0f, 1.0f); color = rgb(int(255 * q), int(190 * (1.0f - std::abs(q - 0.5f) * 2.0f)), int(255 * (1.0f - q)));
        } else if (debugView == DebugView::Divergence && engine.fill[index] >= MIN_RENDER_FILL) {
            float q = std::clamp(engine.divergenceField[index] / divergenceScale, -1.0f, 1.0f); color = q >= 0 ? rgb(245, int(80 * (1.0f - q)), 50) : rgb(40, int(80 * (1.0f + q)), 245);
        } else if (debugView == DebugView::Normal ? visualLiquid[static_cast<size_t>(index)] != 0 : engine.fill[index] >= MIN_RENDER_FILL) {
            if (flatPaint) {
                bool edge = false;
                for (int k = 0; k < 4; ++k) {
                    if (!visibleLiquid(x + n4x[k], y + n4y[k])) { edge = true; break; }
                }
                int lr, lg, lb;
                if (edge) liquidAppearance(engine, index, kWaterRimR, kWaterRimG, kWaterRimB, lr, lg, lb);
                else liquidAppearance(engine, index, kWaterR, kWaterG, kWaterB, lr, lg, lb);
                color = rgb(lr, lg, lb);
            } else {
                int depth = liquidDepth[static_cast<size_t>(index)];
                int lr, lg, lb;
                if (debugView == DebugView::Normal && depth >= 0) {
                    if (depth == 0) liquidAppearance(engine, index, kWaterRimR, kWaterRimG, kWaterRimB, lr, lg, lb);
                    else {
                        float u = 1.0f - std::min(1.0f, static_cast<float>(depth) / 8.0f);
                        uint32_t base = lerpRgb(kWaterDarkR, kWaterDarkG, kWaterDarkB, kWaterLightR, kWaterLightG, kWaterLightB, u);
                        int wr = (base >> 16) & 255, wg = (base >> 8) & 255, wb = base & 255;
                        liquidAppearance(engine, index, wr, wg, wb, lr, lg, lb);
                    }
                } else {
                    liquidAppearance(engine, index, kWaterR, kWaterG, kWaterB, lr, lg, lb);
                }
                color = rgb(lr, lg, lb);
            }
        }
        if (debugView == DebugView::Fill && !engine.solid[index]) color = fillDebugColor(engine.fill[index]);
        if (debugView == DebugView::Chunks && engine.chunkActivity[(y / CHUNK) * CHUNK_W + x / CHUNK] && !engine.solid[index]) color = engine.fill[index] >= MIN_RENDER_FILL ? rgb(35, 180, 225) : rgb(25, 55, 47);
        if (debugView == DebugView::Chunks && (x % CHUNK == 0 || y % CHUNK == 0)) color = rgb(80, 110, 90);
        engine.pixels[index] = color;
    }

    if (debugView == DebugView::Normal && liquidGlow && !flatPaint) {
        constexpr float glowRadius = 5.0f;
        std::vector<float> glowDist(static_cast<size_t>(GW * GH), 1.0e9f);
        std::vector<int> glowQueue;
        glowQueue.reserve(static_cast<size_t>(GW * GH) / 4);
        auto seedGlow = [&](int x, int y) {
            if (!FluidEngine::inside(x, y) || engine.solid[FluidEngine::ci(x, y)] || engine.dynamicSolid[FluidEngine::ci(x, y)]) return;
            int i = FluidEngine::ci(x, y);
            if (glowDist[static_cast<size_t>(i)] <= 0.0f) return;
            glowDist[static_cast<size_t>(i)] = 0.0f;
            glowQueue.push_back(i);
        };
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x)
            if (visibleLiquid(x, y)) seedGlow(x, y);
        for (SplashParticle const &p : engine.splashes) seedGlow(static_cast<int>(p.x), static_cast<int>(p.y));
        for (size_t head = 0; head < glowQueue.size(); ++head) {
            int i = glowQueue[head];
            int x = i % GW, y = i / GW;
            float d = glowDist[static_cast<size_t>(i)];
            for (int oy = -1; oy <= 1; ++oy) for (int ox = -1; ox <= 1; ++ox) {
                if (ox == 0 && oy == 0) continue;
                int nx = x + ox, ny = y + oy;
                if (!FluidEngine::inside(nx, ny)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (engine.solid[ni] || engine.dynamicSolid[ni]) continue;
                float nd = d + ((ox != 0 && oy != 0) ? 1.41421356f : 1.0f);
                if (nd > glowRadius || nd >= glowDist[static_cast<size_t>(ni)]) continue;
                glowDist[static_cast<size_t>(ni)] = nd;
                glowQueue.push_back(ni);
            }
        }
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int index = FluidEngine::ci(x, y);
            float d = glowDist[static_cast<size_t>(index)];
            if (d > glowRadius) continue;
            if (d <= 0.0f) continue;
            if (engine.solid[index] || engine.dynamicSolid[index]) continue;
            float u = 1.0f - d / glowRadius;
            engine.pixels[index] = mixToward(engine.pixels[index], kWaterLightR, kWaterLightG, kWaterLightB,
                std::min(1.0f, u * u * (1.20f + 0.45f * u)));
        }
    }
    for (SplashParticle const &p : engine.splashes) {
        int x = static_cast<int>(p.x), y = static_cast<int>(p.y);
        if (!FluidEngine::inside(x, y) || engine.solid[FluidEngine::ci(x, y)] || engine.dynamicSolid[FluidEngine::ci(x, y)]) continue;
        if (debugView == DebugView::Fill) continue;
        int lr = kWaterRimR, lg = kWaterRimG, lb = kWaterRimB;
        float vol = std::max(p.volume, 1.0e-8f);
        float h = std::clamp(p.honey / vol, 0.0f, 1.0f);
        if (h > 0.001f) tintChannels(lr, lg, lb, kHoneyRimR, kHoneyRimG, kHoneyRimB, h);
        float ir = std::clamp(p.dyeR / vol, 0.0f, 1.0f);
        float ig = std::clamp(p.dyeG / vol, 0.0f, 1.0f);
        float ib = std::clamp(p.dyeB / vol, 0.0f, 1.0f);
        float dyeI = std::clamp(std::max(ir, std::max(ig, ib)), 0.0f, 1.0f);
        if (dyeI > 0.02f) tintChannels(lr, lg, lb, int(ir * 255.0f), int(ig * 255.0f), int(ib * 255.0f), dyeI);
        engine.pixels[FluidEngine::ci(x, y)] = debugView == DebugView::Normal ? rgb(lr, lg, lb) : rgb(115, 220, 255);
    }
    auto occupantAt = [&](int x, int y) {
        if (!FluidEngine::inside(x, y)) return -1;
        return rigid.occupant[static_cast<size_t>(FluidEngine::ci(x, y))];
    };
    auto sameBody = [&](int x, int y, int body) {
        return body >= 0 && occupantAt(x, y) == body;
    };
    std::vector<int> rigidDepth(static_cast<size_t>(GW * GH), -1);
    if (debugView != DebugView::Rigid && !flatPaint) {
        std::vector<int> depthQueue;
        depthQueue.reserve(static_cast<size_t>(GW * GH) / 8);
        constexpr int ox[4] = {-1, 1, 0, 0};
        constexpr int oy[4] = {0, 0, -1, 1};
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int body = occupantAt(x, y);
            if (body < 0) continue;
            bool edge = false;
            for (int k = 0; k < 4; ++k) {
                if (!sameBody(x + ox[k], y + oy[k], body)) { edge = true; break; }
            }
            if (!edge) continue;
            int i = FluidEngine::ci(x, y);
            rigidDepth[static_cast<size_t>(i)] = 0;
            depthQueue.push_back(i);
        }
        for (size_t head = 0; head < depthQueue.size(); ++head) {
            int i = depthQueue[head];
            int x = i % GW, y = i / GW;
            int body = occupantAt(x, y);
            int d = rigidDepth[static_cast<size_t>(i)];
            for (int k = 0; k < 4; ++k) {
                int nx = x + ox[k], ny = y + oy[k];
                if (!sameBody(nx, ny, body)) continue;
                int ni = FluidEngine::ci(nx, ny);
                if (rigidDepth[static_cast<size_t>(ni)] >= 0) continue;
                rigidDepth[static_cast<size_t>(ni)] = d + 1;
                depthQueue.push_back(ni);
            }
        }
    }
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int index = FluidEngine::ci(x, y);
        MaterialId pendingMat = rigid.pending[static_cast<size_t>(index)];
        MaterialId occMat = rigid.worldCellMaterial(x, y);
        if (debugView == DebugView::Rigid) {
            if (occMat != MATERIAL_EMPTY) engine.pixels[index] = rgb(210, 70, 90);
            else if (pendingMat != MATERIAL_EMPTY) engine.pixels[index] = rgb(230, 180, 70);
            continue;
        }
        if (occMat != MATERIAL_EMPTY) {
            MaterialDefinition const &mat = materialDef(occMat);
            float wetness = 0.0f;
            if (mat.moistureCapacity > 1.0e-8f)
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
                    int nx = x + n4x[k], ny = y + n4y[k];
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
            // Wet, degraded, and structurally cracked pixels darken from real state.
            float stain = (1.0f - 0.46f * wetness) * (1.0f - 0.38f * dmg) * (1.0f - 0.62f * crack);
            if (flatPaint) {
                int body = occupantAt(x, y);
                bool edge = body < 0;
                if (body >= 0) {
                    for (int k = 0; k < 4; ++k) {
                        if (!sameBody(x + n4x[k], y + n4y[k], body)) { edge = true; break; }
                    }
                }
                float s = (edge ? 0.58f : 1.0f) * stain;
                engine.pixels[index] = scaleRgb(mat.colorR, mat.colorG, mat.colorB, s);
            } else {
                int depth = std::max(0, rigidDepth[static_cast<size_t>(index)]);
                float s = 1.22f;
                if (depth > 0) {
                    float u = 1.0f - std::min(1.0f, static_cast<float>(depth) / 7.0f);
                    s = 0.62f + 0.38f * u;
                }
                int cr = std::max(0, static_cast<int>(mat.colorR * (1.0f - 0.16f * wetness)));
                int cg = std::max(0, static_cast<int>(mat.colorG * (1.0f - 0.22f * wetness)));
                int cb = std::max(0, static_cast<int>(mat.colorB * (1.0f - 0.10f * wetness)));
                engine.pixels[index] = scaleRgb(cr, cg, cb, s * stain);
            }
        } else if (pendingMat != MATERIAL_EMPTY) {
            MaterialDefinition const &mat = materialDef(pendingMat);
            engine.pixels[index] = rgb(std::min(255, mat.colorR + 40), std::min(255, mat.colorG + 40), std::min(255, mat.colorB + 20));
        }
    }
    if (rigid.debugOverlay || debugView == DebugView::Rigid) {
        auto plot = [&](int x, int y, uint32_t color) {
            if (FluidEngine::inside(x, y)) engine.pixels[FluidEngine::ci(x, y)] = color;
        };
        for (RigidBody const &b : rigid.bodies) {
            int cx = static_cast<int>(std::floor(b.x)), cy = static_cast<int>(std::floor(b.y));
            uint32_t com = b.sleeping ? rgb(80, 200, 90) : (b.supported ? rgb(240, 220, 70) : rgb(255, 160, 60));
            plot(cx, cy, com); plot(cx + 1, cy, com); plot(cx, cy + 1, com);
            int vx = cx + static_cast<int>(std::round(b.vx * 0.12f));
            int vy = cy + static_cast<int>(std::round(b.vy * 0.12f));
            plot(vx, vy, rgb(255, 120, 40));
            int ox = cx + static_cast<int>(std::round(-b.omega * 0.35f));
            int oy = cy + static_cast<int>(std::round(b.omega * 0.12f));
            plot(ox, oy, rgb(220, 90, 220));
            int px = cx + static_cast<int>(std::round(b.debugPosCorrX * 8.0f));
            int py = cy + static_cast<int>(std::round(b.debugPosCorrY * 8.0f));
            plot(px, py, rgb(180, 255, 80));
        }
        for (RigidContact const &c : rigid.lastContacts) {
            int cx = static_cast<int>(std::floor(c.x));
            int cy = static_cast<int>(std::floor(c.y));
            plot(cx, cy, rgb(255, 40, 40));
            plot(cx + static_cast<int>(std::round(c.nx * 2.0f)), cy + static_cast<int>(std::round(c.ny * 2.0f)), rgb(80, 220, 255));
            float tx = -c.ny, ty = c.nx;
            float js = std::max(-2.0f, std::min(2.0f, c.jt * 0.35f));
            plot(cx + static_cast<int>(std::round(tx * js)), cy + static_cast<int>(std::round(ty * js)), rgb(255, 180, 40));
            int pd = std::max(1, static_cast<int>(std::round(c.penetration * 4.0f)));
            plot(cx - static_cast<int>(std::round(c.nx * static_cast<float>(pd))),
                cy - static_cast<int>(std::round(c.ny * static_cast<float>(pd))), rgb(255, 90, 160));
        }
    }
    if (rigid.grab.active) {
        auto plot = [&](int x, int y, uint32_t color) {
            if (FluidEngine::inside(x, y)) engine.pixels[FluidEngine::ci(x, y)] = color;
        };
        int x0 = static_cast<int>(std::floor(rigid.grab.worldX));
        int y0 = static_cast<int>(std::floor(rigid.grab.worldY));
        int x1 = static_cast<int>(std::floor(rigid.grab.targetX));
        int y1 = static_cast<int>(std::floor(rigid.grab.targetY));
        int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1, dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1, err = dx + dy;
        for (;;) {
            plot(x0, y0, rgb(70, 210, 110));
            if (x0 == x1 && y0 == y1) break;
            int e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
        plot(static_cast<int>(std::floor(rigid.grab.worldX)), static_cast<int>(std::floor(rigid.grab.worldY)), rgb(255, 240, 80));
    }
    if (debugView == DebugView::Temperature) {
        for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
            int index = FluidEngine::ci(x, y);
            ThermalCellSample s = ThermalEngine::sampleCell(engine, rigid, gas, x, y);
            int tr = 0, tg = 0, tb = 0;
            float a = 0.0f;
            tempVizOverlay(s.temperatureK, s.kind, tr, tg, tb, a);
            if (a <= 0.001f) continue;
            engine.pixels[index] = mixToward(engine.pixels[index], tr, tg, tb, a);
        }
    }
    overlayLineGhost();
}

void blitCanvas(HDC dc) {
    RECT const &c = shell.layout.canvasInner;
    int dw = std::max(1, static_cast<int>(c.right - c.left));
    int dh = std::max(1, static_cast<int>(c.bottom - c.top));
    BITMAPINFO bmi{}; bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER); bmi.bmiHeader.biWidth = GW; bmi.bmiHeader.biHeight = -GH;
    bmi.bmiHeader.biPlanes = 1; bmi.bmiHeader.biBitCount = 32; bmi.bmiHeader.biCompression = BI_RGB;
    SetStretchBltMode(dc, COLORONCOLOR);
    StretchDIBits(dc, c.left, c.top, dw, dh, 0, 0, GW, GH, engine.pixels.data(), &bmi, DIB_RGB_COLORS, SRCCOPY);
}

void render(HWND hwnd, HDC dc) {
    auto renderStart = FluidEngine::Clock::now();
    refreshLayout(hwnd);
    fillWorldPixels();
    ui::View view = makeView();
    ui::drawShell(dc, shell, view);
    blitCanvas(dc);
    drawBrushOverlay(dc);
    ui::drawSettingsPanel(dc, shell, view, shell.layout.client);
    ui::drawHoverTip(dc, shell, view);
    engine.timingAccum.render += FluidEngine::elapsedMs(renderStart);
    if (++engine.renderSamples >= 60) {
        engine.timingAverage.render = engine.timingAccum.render / engine.renderSamples;
        engine.timingAccum.render = 0.0;
        engine.renderSamples = 0;
    }
}

void worldTick() {
    auto tickStart = FluidEngine::Clock::now();
    engine.syncWorkerPool();
    rigid.step(engine, PHYSICS_DT);
    engine.simulationTick();
    gas.simulationTick(engine);
    rigid.gatherFluidForces(engine);
    gas.applyPressureForces(rigid, engine);
    thermal.simulationTick(engine, rigid, gas, PHYSICS_DT);
    adaptAutoQuality(FluidEngine::elapsedMs(tickStart));
}

void invalidate() { if (mainWindow) InvalidateRect(mainWindow, nullptr, FALSE); }

void preparePlacement() {
    rigid.placePowder = shell.placePowder && shell.category == ui::Category::Solids;
    rigid.powderParticleSize = std::clamp(shell.powderParticleSize, 1, 8);
    rigid.placeAnchored = shell.placeAnchored && shell.category == ui::Category::Solids && !rigid.placePowder;
    rigid.placeSleeping = shell.placeSleeping && shell.category == ui::Category::Solids
        && !shell.placeAnchored && !rigid.placePowder;
}

void cancelLineStroke() {
    linePainting = false;
    lineStartX = lineStartY = lineEndX = lineEndY = -1;
}

LiquidPaint liquidPaintFromShell() {
    LiquidPaint p;
    p.asHoney = shell.palette == ui::PaletteId::Honey;
    p.dyeOnly = shell.dyeOnly;
    p.dyeStrength = shell.dyeStrength;
    shell.dyeChannels(p.dyeR, p.dyeG, p.dyeB);
    return p;
}

void commitLineStroke() {
    if (!linePainting) return;
    preparePlacement();
    if (activeTool == Tool::Rigid) {
        rigid.clearPending();
        rigid.paintPendingLine(lineStartX, lineStartY, lineEndX, lineEndY, brushRadius);
        rigid.commitPending(engine);
        rigid.syncOccupancy(engine);
        gas.handleWorldEdit(engine);
    } else if (isEnergyTool(activeTool)) {
        int dx = std::abs(lineEndX - lineStartX), sx = lineStartX < lineEndX ? 1 : -1;
        int dy = -std::abs(lineEndY - lineStartY), sy = lineStartY < lineEndY ? 1 : -1, error = dx + dy;
        int x0 = lineStartX, y0 = lineStartY;
        for (;;) {
            paintEnergyDisc(x0, y0);
            if (x0 == lineEndX && y0 == lineEndY) break;
            int twice = 2 * error;
            if (twice >= dy) { error += dy; x0 += sx; }
            if (twice <= dx) { error += dx; y0 += sy; }
        }
    } else {
        engine.paintLine(lineStartX, lineStartY, lineEndX, lineEndY, activeTool, brushRadius, liquidPaintFromShell());
        if (activeTool == Tool::Eraser) rigid.eraseLine(lineStartX, lineStartY, lineEndX, lineEndY, brushRadius, engine);
        engine.finalizePaint();
        engine.flushPaintDirty();
        gas.handleWorldEdit(engine);
    }
    cancelLineStroke();
}

void paintEnergyDisc(int cx, int cy) {
    if (isThermalEnergyTool(activeTool)) {
        float sign = activeTool == Tool::Heat ? 1.0f : -1.0f;
        thermal.applyBrush(engine, rigid, gas, cx, cy, brushRadius, sign * shell.heatPower, PHYSICS_DT);
        return;
    }
    float sign = activeTool == Tool::Pressurize ? 1.0f : -1.0f;
    gas.applyPressureBrush(engine, cx, cy, brushRadius,
        sign * gas.config.brushAtmPerSec * shell.heatPower, PHYSICS_DT);
}

void beginPaintStroke(int x, int y) {
    if (!shell.hasPlacement()) return;
    preparePlacement();
    clientToGrid(x, y, lastPaintX, lastPaintY);
    if (isEnergyTool(activeTool)) {
        paintEnergyDisc(lastPaintX, lastPaintY);
        return;
    }
    if (activeTool == Tool::Rigid) rigid.paintPendingDisc(lastPaintX, lastPaintY, brushRadius);
    else {
        engine.paintDisc(lastPaintX, lastPaintY, activeTool, brushRadius, liquidPaintFromShell());
        if (activeTool == Tool::Eraser) rigid.eraseDisc(lastPaintX, lastPaintY, brushRadius, engine);
        engine.finalizePaint();
        gas.handleWorldEdit(engine);
    }
}

void continuePaintStroke(int x, int y) {
    if (!shell.hasPlacement() || lastPaintX < 0) return;
    int gx, gy; clientToGrid(x, y, gx, gy);
    if (isEnergyTool(activeTool)) {
        paintEnergyDisc(gx, gy);
        lastPaintX = gx; lastPaintY = gy;
        return;
    }
    if (activeTool == Tool::Rigid) rigid.paintPendingLine(lastPaintX, lastPaintY, gx, gy, brushRadius);
    else {
        engine.paintLine(lastPaintX, lastPaintY, gx, gy, activeTool, brushRadius, liquidPaintFromShell());
        if (activeTool == Tool::Eraser) rigid.eraseLine(lastPaintX, lastPaintY, gx, gy, brushRadius, engine);
        engine.finalizePaint();
        gas.handleWorldEdit(engine);
    }
    lastPaintX = gx; lastPaintY = gy;
}

LRESULT CALLBACK wndProc(HWND hwnd, UINT message, WPARAM wp, LPARAM lp) {
    switch (message) {
        case WM_LBUTTONDOWN: {
            SetCapture(hwnd);
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            refreshLayout(hwnd);
            shell.mouseX = x; shell.mouseY = y; shell.mouseDown = true;
            ui::View view = makeView();
            ui::MenuCmd settingsCmd = ui::MenuCmd::None;
            ui::SettingsMouseResult settingsHit = ui::handleSettingsMouseDown(shell, view, x, y, settingsCmd);
            if (settingsHit == ui::SettingsMouseResult::Drag) {
                mousePainting = false;
                cancelLineStroke();
                invalidate();
                return 0;
            }
            if (settingsHit != ui::SettingsMouseResult::Miss) {
                if (settingsHit == ui::SettingsMouseResult::Command) handleMenuCommand(settingsCmd);
                mousePainting = false;
                cancelLineStroke();
                ReleaseCapture();
                invalidate();
                return 0;
            }
            ui::HitId hit = ui::hitTest(shell, x, y, settingsOpen);
            shell.hoverId = static_cast<int>(hit);
            if (hit != ui::HitId::None) {
                if (hit != ui::HitId::Search && hit != ui::HitId::Console) {
                    shell.searchFocused = false;
                    shell.consoleFocused = false;
                }
                handleHit(hit);
                mousePainting = false;
                cancelLineStroke();
                ReleaseCapture();
                invalidate();
                return 0;
            }
            shell.searchFocused = false;
            shell.consoleFocused = false;
            if (inCanvas(x, y) && activeTool == Tool::Grab) {
                float gx, gy;
                clientToGridF(x, y, gx, gy);
                mouseGrabbing = rigid.beginGrab(gx, gy);
                mousePainting = false;
                cancelLineStroke();
                if (mouseGrabbing) {
                    rigid.updateGrabTarget(gx, gy, (wp & MK_SHIFT) != 0);
                    shell.log(ui::tr("log_grab"));
                }
            } else if (inCanvas(x, y) && shell.hasPlacement()) {
                if (wp & MK_SHIFT) {
                    mousePainting = false;
                    linePainting = true;
                    preparePlacement();
                    clientToGrid(x, y, lineStartX, lineStartY);
                    lineEndX = lineStartX;
                    lineEndY = lineStartY;
                } else {
                    cancelLineStroke();
                    mousePainting = true;
                    beginPaintStroke(x, y);
                }
            } else {
                mousePainting = false;
                cancelLineStroke();
                ReleaseCapture();
            }
            invalidate(); return 0;
        }
        case WM_MOUSEMOVE: {
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            refreshLayout(hwnd);
            shell.mouseX = x; shell.mouseY = y;
            shell.hoverId = static_cast<int>(ui::hitTest(shell, x, y, settingsOpen));
            TRACKMOUSEEVENT tme{};
            tme.cbSize = sizeof(tme);
            tme.dwFlags = TME_LEAVE;
            tme.hwndTrack = hwnd;
            TrackMouseEvent(&tme);
            if (shell.settingsScrollDragging && (wp & MK_LBUTTON)) {
                ui::dragSettingsScroll(shell, y);
                invalidate();
                return 0;
            }
            if (shell.settingsDragging && (wp & MK_LBUTTON)) {
                ui::dragSettingsWindow(shell, x, y);
                invalidate();
                return 0;
            }
            if (inCanvas(x, y)) clientToGrid(x, y, hoverX, hoverY);
            else hoverX = hoverY = -1;
            if (mouseGrabbing && (wp & MK_LBUTTON) && inCanvas(x, y)) {
                float gx, gy;
                clientToGridF(x, y, gx, gy);
                rigid.updateGrabTarget(gx, gy, (wp & MK_SHIFT) != 0);
            } else if (linePainting && (wp & MK_LBUTTON)) {
                clientToGrid(x, y, lineEndX, lineEndY);
            } else if (!shell.settingsDragging && !shell.settingsScrollDragging && mousePainting && (wp & MK_LBUTTON) && inCanvas(x, y))
                continuePaintStroke(x, y);
            invalidate();
            return 0;
        }
        case WM_MOUSELEAVE:
            shell.hoverId = 0;
            shell.mouseX = -1;
            shell.mouseY = -1;
            invalidate();
            return 0;
        case WM_LBUTTONUP:
            shell.mouseDown = false;
            shell.settingsDragging = false;
            shell.settingsScrollDragging = false;
            if (mouseGrabbing) {
                rigid.endGrab();
                mouseGrabbing = false;
            }
            if (linePainting) {
                commitLineStroke();
            } else if (activeTool == Tool::Rigid && mousePainting) {
                rigid.commitPending(engine);
                rigid.syncOccupancy(engine);
                gas.handleWorldEdit(engine);
            }
            else if (mousePainting) engine.flushPaintDirty();
            mousePainting = false; lastPaintX = lastPaintY = -1; ReleaseCapture();
            invalidate();
            return 0;
        case WM_RBUTTONDOWN: {
            refreshLayout(hwnd);
            int x = GET_X_LPARAM(lp), y = GET_Y_LPARAM(lp);
            if (settingsOpen && ui::settingsCoversPoint(shell, x, y)) return 0;
            if (!inCanvas(x, y)) return 0;
            int gx, gy; clientToGrid(x, y, gx, gy);
            int body = rigid.bodyAtCell(gx, gy);
            if (body >= 0) { rigid.removeBody(body, engine); shell.log(ui::tr("log_rigid_deleted")); }
            invalidate();
            return 0;
        }
        case WM_MOUSEWHEEL: {
            POINT p{};
            p.x = GET_X_LPARAM(lp);
            p.y = GET_Y_LPARAM(lp);
            ScreenToClient(hwnd, &p);
            refreshLayout(hwnd);
            if (settingsOpen && ui::handleSettingsWheel(shell, p.x, p.y, GET_WHEEL_DELTA_WPARAM(wp))) {
                invalidate();
                return 0;
            }
            brushRadius = std::clamp(brushRadius + (GET_WHEEL_DELTA_WPARAM(wp) > 0 ? 1 : -1), 1, 14);
            invalidate(); return 0;
        }
        case WM_CHAR:
            if (shell.searchFocused) {
                if (wp == 8) { if (!shell.searchText.empty()) shell.searchText.pop_back(); }
                else if (wp >= 32 && wp != 127) { if (shell.searchText.size() < 48) shell.searchText.push_back(static_cast<wchar_t>(wp)); }
                invalidate(); return 0;
            }
            if (shell.consoleFocused) {
                if (wp == 8) { if (!shell.consoleDraft.empty()) shell.consoleDraft.pop_back(); }
                else if (wp == L'\r') {
                    if (!shell.consoleDraft.empty()) {
                        shell.log((L"> " + shell.consoleDraft).c_str());
                        shell.log(ui::tr("log_no_command"));
                        shell.consoleDraft.clear();
                    }
                } else if (wp >= 32 && wp != 127) { if (shell.consoleDraft.size() < 80) shell.consoleDraft.push_back(static_cast<wchar_t>(wp)); }
                invalidate(); return 0;
            }
            return 0;
        case WM_KEYDOWN:
            if (shell.searchFocused || shell.consoleFocused) {
                if (wp == VK_ESCAPE) { shell.searchFocused = false; shell.consoleFocused = false; settingsOpen = false; invalidate(); }
                return 0;
            }
            if (wp == VK_ESCAPE) {
                if (linePainting) {
                    cancelLineStroke();
                    shell.mouseDown = false;
                    ReleaseCapture();
                    invalidate();
                    return 0;
                }
                settingsOpen = !settingsOpen;
                if (settingsOpen) ui::ensureSettingsPlacement(shell);
                else {
                    shell.settingsDragging = false;
                    shell.settingsScrollDragging = false;
                    shell.settingsFlyout = ui::SettingsFlyout::None;
                }
                invalidate();
                return 0;
            }
            if (wp == '1') shell.applyCategory(ui::Category::Tools, activeTool, rigid.drawMaterial);
            else if (wp == '2') shell.applyCategory(ui::Category::Fluids, activeTool, rigid.drawMaterial);
            else if (wp == '3') shell.applyCategory(ui::Category::Solids, activeTool, rigid.drawMaterial);
            else if (wp == '4') shell.applyCategory(ui::Category::Misc, activeTool, rigid.drawMaterial);
            else if (wp == VK_SPACE) paused = !paused; else if (wp == 'S' && paused) worldTick();
            else if (wp == 'A') engine.config.velocityAdvection = nextVelocityAdvection(engine.config.velocityAdvection);
            else if (wp == VK_OEM_4) brushRadius = std::max(1, brushRadius - 1); else if (wp == VK_OEM_6) brushRadius = std::min(14, brushRadius + 1);
            else if (wp == VK_OEM_MINUS || wp == VK_SUBTRACT) changeTickRate(-1); else if (wp == VK_OEM_PLUS || wp == VK_ADD) changeTickRate(1);
            invalidate(); return 0;
        case WM_SETCURSOR: {
            POINT p{}; GetCursorPos(&p); ScreenToClient(hwnd, &p);
            refreshLayout(hwnd);
            if (shell.settingsDragging || (settingsOpen && ui::ptIn(shell.settingsTitleBar, p.x, p.y))) {
                SetCursor(LoadCursor(nullptr, IDC_SIZEALL));
                return TRUE;
            }
            if (settingsOpen && (ui::ptIn(shell.settingsThumb, p.x, p.y) || ui::ptIn(shell.settingsScrollTrack, p.x, p.y))) {
                SetCursor(LoadCursor(nullptr, IDC_ARROW));
                return TRUE;
            }
            if (inCanvas(p.x, p.y) && !(settingsOpen && ui::settingsCoversPoint(shell, p.x, p.y))) {
                SetCursor(LoadCursor(nullptr, activeTool == Tool::Grab ? IDC_HAND : IDC_CROSS));
                return TRUE;
            }
            SetCursor(LoadCursor(nullptr, IDC_ARROW));
            return TRUE;
        }
        case WM_ERASEBKGND: return 1;
        case WM_PAINT: {
            PAINTSTRUCT ps{}; HDC dc = BeginPaint(hwnd, &ps); RECT rc{}; GetClientRect(hwnd, &rc);
            ensureBackbuffer(dc, rc.right, rc.bottom);
            render(hwnd, backDc);
            BitBlt(dc, 0, 0, rc.right, rc.bottom, backDc, 0, 0, SRCCOPY);
            EndPaint(hwnd, &ps); return 0;
        }
        case WM_SIZE:
            releaseBackbuffer();
            return 0;
        case WM_DESTROY:
            releaseBackbuffer();
            shell.releaseFonts();
            PostQuitMessage(0); return 0;
    }
    return DefWindowProcW(hwnd, message, wp, lp);
}

} // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, LPWSTR commandLine, int show) {
    char const *kUiLanguage = "english";
    ui::loadLanguage(kUiLanguage);
    if (commandLine && wcsstr(commandLine, L"--scale-benchmark")) { engine.runScaleBenchmark(); return 0; }
    if (commandLine && wcsstr(commandLine, L"--thread-benchmark")) {
        engine.runThreadBenchmark();
        std::ofstream out(miscFile("thread_benchmark_" + std::to_string(GW) + "x" + std::to_string(GH) + ".tsv"), std::ios::app);
        int counts[] = {1, 2, 4, 0};
        char const *modes[] = {"1", "2", "4", "auto"};
        int ticks = 90;
        for (int m = 0; m < 4; ++m) {
            engine.config.workerCount = counts[m];
            engine.syncWorkerPool();
            engine.residualConsolidationEnabled = true;
            rigid.loadTestScene(engine, 5);
            for (int i = 0; i < 20; ++i) {
                rigid.step(engine, PHYSICS_DT);
                engine.simulationTick();
                rigid.gatherFluidForces(engine);
            }
            engine.timingAccum = TimingAverages{};
            engine.timingAverage = TimingAverages{};
            engine.timingTicks = 0;
            auto start = FluidEngine::Clock::now();
            for (int i = 0; i < ticks; ++i) {
                rigid.step(engine, PHYSICS_DT);
                engine.simulationTick();
                rigid.gatherFluidForces(engine);
            }
            engine.rebuildActivityAndMetrics();
            double perTick = FluidEngine::elapsedMs(start) / ticks;
            out << engine.lastResolvedWorkers << '\t' << modes[m] << '\t' << "E_rigid_pool" << '\t' << ticks << '\t'
                << perTick << '\t' << engine.timingAverage.pressure << '\t' << engine.timingAverage.advection << '\t'
                << engine.timingAverage.transport << '\t' << engine.timingAverage.surface << '\t' << engine.timingAverage.residual << '\t'
                << rigid.lastStepMs << '\t' << engine.activeFluidCells << '\t' << engine.workCounts.pressureCells << '\t'
                << engine.currentVolume << '\t' << engine.expectedVolume << '\t' << engine.volumeError << '\t'
                << (engine.lastPressureParallel ? 1 : 0) << '\n';
        }
        engine.config.workerCount = 1;
        engine.syncWorkerPool();
        return 0;
    }
    if (commandLine && wcsstr(commandLine, L"--solid-diag")) { rigid.runSolidDiagnostics(engine); return 0; }
    if (commandLine && wcsstr(commandLine, L"--thermal-diag")) { thermal.runDiagnostics(engine, rigid, gas); return 0; }
    if (commandLine && wcsstr(commandLine, L"--gas-diag")) { gas.runDiagnostics(engine, rigid); return 0; }
    if (commandLine && wcsstr(commandLine, L"--benchmark")) { engine.runHeadlessBenchmark(); return 0; }
    if (commandLine && wcsstr(commandLine, L"--advection-benchmark")) { engine.runAdvectionBenchmark(); return 0; }
    if (commandLine && wcsstr(commandLine, L"--liquid-diag")) { engine.runLiquidBugDiagnostics(); return 0; }
    if (commandLine && wcsstr(commandLine, L"--rigid-benchmark")) { rigid.runConservationBenchmark(engine); return 0; }
    if (commandLine && wcsstr(commandLine, L"--rigid-contact-diag")) { rigid.runContactDiagnostics(engine); return 0; }
    SetProcessDPIAware();
    WNDCLASSW wc{}; wc.style = CS_HREDRAW | CS_VREDRAW; wc.lpfnWndProc = wndProc; wc.hInstance = instance; wc.hCursor = LoadCursor(nullptr, IDC_ARROW); wc.lpszClassName = L"PipaceFluidWindow";
    if (!RegisterClassW(&wc)) return 1;
    mainWindow = CreateWindowExW(0, wc.lpszClassName, ui::tr("window_title"), WS_OVERLAPPEDWINDOW | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, 1440, 860, nullptr, nullptr, instance, nullptr);
    if (!mainWindow) return 2;
    engine.resetWorld();
    gas.resetAmbient(engine);
    thermal.seedAmbient(engine, rigid, gas);
    shell.applyPalette(ui::PaletteId::Water, activeTool, rigid.drawMaterial);
    shell.log(ui::tr("log_ready"));
    ShowWindow(mainWindow, show); UpdateWindow(mainWindow);
    LARGE_INTEGER frequency{}, previous{}; QueryPerformanceFrequency(&frequency); QueryPerformanceCounter(&previous); double accumulator = 0.0;
    MSG msg{}; bool running = true;
    while (running) {
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) { if (msg.message == WM_QUIT) { running = false; break; } TranslateMessage(&msg); DispatchMessageW(&msg); }
        if (!running) break;
        LARGE_INTEGER now{}; QueryPerformanceCounter(&now); double elapsed = static_cast<double>(now.QuadPart - previous.QuadPart) / frequency.QuadPart; previous = now;
        bool updated = false;
        if (paused) accumulator = 0.0;
        else {
            accumulator += std::min(elapsed, 0.1) * speedScales[speedScaleIndex];
            double tickStep = 1.0 / static_cast<double>(std::max(10, engine.config.physicsHz));
            int catchUp = 0;
            int catchUpMax = std::max(1, engine.config.catchUpTicks);
            while (accumulator >= tickStep && catchUp < catchUpMax) { worldTick(); accumulator -= tickStep; updated = true; ++catchUp; }
            if (catchUp == catchUpMax) accumulator = 0.0;
        }
        if (updated) invalidate();
        Sleep(1);
    }
    return 0;
}
