#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include "ui/UiCommands.h"
#include "ui/ToolWindow.h"
#include "render/WorldVisual.h"
#include "ui/UiLayout.h"

#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"

#include <string>
#include <vector>

struct ThermalEngine;

namespace ui {

struct ShellState {
    Layout layout{};
    Category category = Category::Fluids;
    PaletteId palette = PaletteId::Water;
    bool placeAnchored = false;
    bool placeSleeping = false;
    bool placePowder = false;
    int powderParticleSize = 2;
    float heatPower = 1.0f;
    int dyeSwatch = 0;
    bool dyeOnly = false;
    float dyeStrength = 0.60f;
    bool searchFocused = false;
    bool consoleFocused = false;
    bool consoleExpanded = false;
    std::wstring searchText;
    std::wstring consoleDraft;
    std::vector<std::wstring> consoleLines;
    ToolSettings tools{};
    ToolWindowState toolWin{};
    int hoverId = 0;
    int mouseX = 0;
    int mouseY = 0;
    bool mouseDown = false;
    std::vector<MenuHit> menuHits;
    int settingsX = -1;
    int settingsY = -1;
    int settingsW = 472;
    int settingsH = 560;
    bool settingsDragging = false;
    int settingsDragOX = 0;
    int settingsDragOY = 0;
    RECT settingsTitleBar{};
    RECT settingsBody{};
    RECT settingsFlyoutRc{};
    RECT settingsScrollTrack{};
    RECT settingsThumb{};
    SettingsFlyout settingsFlyout = SettingsFlyout::None;
    int settingsScrollY = 0;
    int settingsContentH = 0;
    bool settingsScrollDragging = false;
    int settingsScrollDragOY = 0;
    HFONT uiFont = nullptr;
    HFONT smallFont = nullptr;
    HFONT consoleFont = nullptr;

    void ensureFonts();
    void releaseFonts();
    void log(wchar_t const *line);
    void applyPalette(PaletteId id, Tool &tool, MaterialId &drawMaterial, bool openToolWindow = true);
    void applyCategory(Category cat, Tool &tool, MaterialId &drawMaterial);
    PaletteId elementAt(int slot) const;
    PaletteId elementAt(Category cat, int slot) const;
    int elementCount() const;
    int elementCount(Category cat) const;
    bool hasPlacement() const;
    wchar_t const *categoryName(Category cat) const;
    wchar_t const *paletteName(PaletteId id) const;
    wchar_t const *paletteHint(PaletteId id) const;
    void dyeChannels(float &r, float &g, float &b) const;
};

struct View {
    FluidEngine const *engine = nullptr;
    RigidBodyEngine const *rigid = nullptr;
    GasEngine const *gas = nullptr;
    ThermalEngine const *thermal = nullptr;
    Tool tool = Tool::Water;
    DebugView debugView = DebugView::Normal;
    int brushRadius = 3;
    bool paused = false;
    bool settingsOpen = false;
    WorldLook worldLook;
    size_t speedIndex = 3;
    float speedValue = 1.0f;
    int hoverX = -1;
    int hoverY = -1;
    float heatPower = 1.0f;
    bool thermalEnabled = true;
};

void computeShellLayout(ShellState &shell, int clientW, int clientH);
HitId hitTest(ShellState const &shell, int x, int y, bool settingsOpen);
void drawShell(HDC dc, ShellState &shell, View const &view);
void drawSettingsPanel(HDC dc, ShellState &shell, View const &view, RECT const &client);
void drawHoverTip(HDC dc, ShellState &shell, View const &view);
SettingsMouseResult handleSettingsMouseDown(ShellState &shell, View const &view, int x, int y, MenuCmd &outCmd);
void dragSettingsWindow(ShellState &shell, int x, int y);
void dragSettingsScroll(ShellState &shell, int y);
bool handleSettingsWheel(ShellState &shell, int x, int y, int wheelDelta);
void clampSettingsWindow(ShellState &shell);
void ensureSettingsPlacement(ShellState &shell);
bool settingsCoversPoint(ShellState const &shell, int x, int y);
RECT settingsPanelRect(ShellState const &shell);

} // namespace ui
