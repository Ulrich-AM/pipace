#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include "fluid/BrushGeom.h"
#include "ui/UiCommands.h"

#include <cstdint>

namespace ui {

struct ShellState;

struct EraseToolSettings {
    int brushSize = 3;
    BrushShape shape = BrushShape::Circle;
    bool deleteSolids = true;
    bool strictSolids = false;
    bool deleteLiquids = true;
    bool deleteGases = false;
};

struct GrabToolSettings {
    float strength = 1.0f;
    bool groupGrab = false;
    float groupRadius = 8.0f;
    bool phantom = false;
};

struct BrushToolSettings {
    int brushSize = 3;
    BrushShape shape = BrushShape::Circle;
    bool colorMode = false;
};

struct TouchToolSettings {
    bool groupTouch = false;
    float groupRadius = 8.0f;
    bool toggleAnchor = false;
};

struct HeatToolSettings {
    float power = 1.0f; // signed; negative cools. Existing add/remove-heat tool.
    int brushSize = 3;
    BrushShape shape = BrushShape::Circle;
};

struct ToolSettings {
    EraseToolSettings erase;
    GrabToolSettings grab;
    BrushToolSettings brush;
    TouchToolSettings touch;
    HeatToolSettings heat;
};

struct ToolWindowState {
    bool open = false;
    bool userClosed = false;
    int x = -1;
    int y = -1;
    int w = 252;
    int h = 220;
    bool dragging = false;
    int dragOX = 0;
    int dragOY = 0;
    RECT panel{};
    RECT titleBar{};
    RECT closeBtn{};
    RECT hits[28]{};
    int hitCount = 0;
    int hitId[28]{};
    char const *hitTip[28]{};
};

bool paletteHasToolWindow(PaletteId id);
void syncToolWindow(ShellState &shell, PaletteId id, bool openIfTool);
void clampToolWindow(ShellState &shell);
void ensureToolWindowPlacement(ShellState &shell);
bool toolWindowCoversPoint(ShellState const &shell, int x, int y);
char const *toolWindowTipKey(ShellState const &shell);
void drawToolWindow(HDC dc, ShellState &shell);
SettingsMouseResult handleToolWindowMouseDown(ShellState &shell, int x, int y);
void dragToolWindow(ShellState &shell, int x, int y);
bool handleToolWindowWheel(ShellState &shell, int x, int y, int wheelDelta);

BrushShape activeBrushShape(ShellState const &shell);
int activeToolBrushSize(ShellState const &shell);
void setActiveToolBrushSize(ShellState &shell, int size);
void writeBackActiveBrushSize(ShellState &shell, int size);

} // namespace ui
