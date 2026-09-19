#include "ui/ToolWindow.h"

#include "ui/UiLanguage.h"
#include "ui/UiShell.h"
#include "ui/UiTheme.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cwchar>
#include <string>

namespace ui {
namespace {

constexpr int kTwHitCap = 28;

enum TwHit : int {
    TwNone = 0,
    TwTitle,
    TwClose,
    TwSpin0M,
    TwSpin0P,
    TwSpin1M,
    TwSpin1P,
    TwShape0,
    TwShape1,
    TwShape2,
    TwSolids,
    TwStrictSolids,
    TwLiquids,
    TwGases,
    TwGroup,
    TwPhantom,
    TwColorMode,
    TwToggleAnchor,
    TwIdle,
    TwDye0,
    TwDye1,
    TwDye2,
    TwDye3,
    TwDye4,
    TwDye5,
    TwDye6,
    TwDye7
};

COLORREF const kDyeRgb[8] = {
    RGB(210, 210, 210), RGB(220, 40, 42), RGB(230, 118, 28), RGB(230, 200, 36),
    RGB(40, 170, 55), RGB(36, 72, 220), RGB(140, 50, 190), RGB(28, 28, 32)
};

void clearHits(ToolWindowState &w) {
    w.hitCount = 0;
    for (int i = 0; i < kTwHitCap; ++i) {
        w.hits[i] = RECT{};
        w.hitId[i] = TwNone;
        w.hitTip[i] = nullptr;
    }
}

void addHit(ToolWindowState &w, RECT const &rc, int id, char const *tip = nullptr) {
    if (w.hitCount >= kTwHitCap) return;
    int i = w.hitCount++;
    w.hits[i] = rc;
    w.hitId[i] = id;
    w.hitTip[i] = tip;
}

int hitAt(ToolWindowState const &w, int x, int y) {
    for (int i = w.hitCount - 1; i >= 0; --i) {
        if (ptIn(w.hits[i], x, y)) return w.hitId[i];
    }
    return TwNone;
}

void drawCheck(HDC dc, RECT box, bool on, bool enabled) {
    fillRect(dc, box, kBlack);
    frameRect(dc, box, enabled ? kBorder : kBorderDim);
    if (on) {
        RECT inn = insetRect(box, 3);
        fillRect(dc, inn, enabled ? kText : kDimText);
    }
}

} // namespace

bool paletteHasToolWindow(PaletteId id) {
    return id == PaletteId::Erase || id == PaletteId::Grab || id == PaletteId::Brush
        || id == PaletteId::Touch || id == PaletteId::Heat;
}

void syncToolWindow(ShellState &shell, PaletteId id, bool openIfTool) {
    if (!paletteHasToolWindow(id))
        return;
    if (!openIfTool)
        return;
    shell.toolWin.userClosed = false;
    shell.toolWin.open = true;
}

char const *toolWindowTipKey(ShellState const &shell) {
    if (!shell.toolWin.open) return nullptr;
    for (int i = shell.toolWin.hitCount - 1; i >= 0; --i) {
        if (!ptIn(shell.toolWin.hits[i], shell.mouseX, shell.mouseY)) continue;
        if (shell.toolWin.hitTip[i] && shell.toolWin.hitTip[i][0]) return shell.toolWin.hitTip[i];
        return nullptr;
    }
    return nullptr;
}

void clampToolWindow(ShellState &shell) {
    RECT const &c = shell.layout.client;
    int w = std::max(180, shell.toolWin.w);
    int h = std::max(80, shell.toolWin.h);
    int minX = static_cast<int>(c.left + 4);
    int maxX = static_cast<int>(std::max(c.left + 4, c.right - w - 4));
    int minY = static_cast<int>(c.top + 4);
    int maxY = static_cast<int>(std::max(c.top + 4, c.bottom - h - 4));
    shell.toolWin.x = std::clamp(shell.toolWin.x, minX, maxX);
    shell.toolWin.y = std::clamp(shell.toolWin.y, minY, maxY);
}

void ensureToolWindowPlacement(ShellState &shell) {
    RECT world = shell.layout.canvasInner;
    if (world.right <= world.left) world = shell.layout.canvas;
    if (shell.toolWin.x < 0 || shell.toolWin.y < 0) {
        shell.toolWin.x = world.left + 10;
        shell.toolWin.y = world.top + 10;
    }
    clampToolWindow(shell);
}

bool toolWindowCoversPoint(ShellState const &shell, int x, int y) {
    if (!shell.toolWin.open) return false;
    return ptIn(shell.toolWin.panel, x, y) != 0;
}

BrushShape activeBrushShape(ShellState const &shell) {
    switch (shell.editTool) {
        case EditTool::Erase: return shell.tools.erase.shape;
        case EditTool::Heat: return shell.tools.heat.shape;
        case EditTool::Brush:
        default: return shell.tools.brush.shape;
    }
}

int activeToolBrushSize(ShellState const &shell) {
    switch (shell.editTool) {
        case EditTool::Erase: return shell.tools.erase.brushSize;
        case EditTool::Heat: return shell.tools.heat.brushSize;
        case EditTool::Brush:
        default: return shell.tools.brush.brushSize;
    }
}

void setActiveToolBrushSize(ShellState &shell, int size) {
    size = std::clamp(size, 1, 14);
    switch (shell.editTool) {
        case EditTool::Erase: shell.tools.erase.brushSize = size; break;
        case EditTool::Heat: shell.tools.heat.brushSize = size; break;
        case EditTool::Brush:
        default: shell.tools.brush.brushSize = size; break;
    }
}

void writeBackActiveBrushSize(ShellState &shell, int size) {
    setActiveToolBrushSize(shell, size);
}

void drawToolWindow(HDC dc, ShellState &shell) {
    if (!shell.toolWin.open) return;
    PaletteId pal = paletteForEditTool(shell.editTool);
    if (!paletteHasToolWindow(pal)) return;
    shell.ensureFonts();
    ensureToolWindowPlacement(shell);
    ToolWindowState &w = shell.toolWin;
    clearHits(w);

    int x = w.x, y0 = w.y;
    int width = w.w;
    int headH = 30;
    int pad = 8;

    auto measureEnd = [&]() {
        switch (shell.editTool) {
            case EditTool::Erase: return 228;
            case EditTool::Grab: return 176;
            case EditTool::Brush: return 218;
            case EditTool::Touch: return 158;
            case EditTool::Heat: return 176;
            default: return 128;
        }
    };
    w.h = measureEnd();
    clampToolWindow(shell);
    x = w.x; y0 = w.y;
    RECT panel{x, y0, x + width, y0 + w.h};
    w.panel = panel;
    fillRect(dc, panel, kPanel);
    frameRect(dc, panel, kBorder);
    RECT inner = insetRect(panel, 1);
    frameRect(dc, inner, kBorderDim);

    RECT closeRc{panel.right - 28, panel.top + 6, panel.right - 8, panel.top + 28};
    RECT title{panel.left + 1, panel.top + 1, closeRc.left - 4, panel.top + headH};
    w.titleBar = title;
    w.closeBtn = closeRc;
    addHit(w, panel, TwIdle, nullptr);
    bool titleOver = ptIn(title, shell.mouseX, shell.mouseY);
    if (titleOver || w.dragging) fillRect(dc, title, RGB(52, 52, 56));
    HFONT old = static_cast<HFONT>(SelectObject(dc, shell.smallFont ? shell.smallFont : shell.uiFont));
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, kText);
    RECT titlePad{title.left + 12, title.top + 6, title.right, title.bottom};
    drawLabel(dc, titlePad, shell.paletteName(pal), kText, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
    bool closeOver = ptIn(closeRc, shell.mouseX, shell.mouseY);
    drawButton(dc, closeRc, tr("tw_close"), closeOver ? BtnState::Hover : BtnState::Normal);
    addHit(w, title, TwTitle, "tip_tw_title");
    addHit(w, closeRc, TwClose, "tip_tw_close");

    int left = x + pad;
    int right = x + width - pad;
    int y = y0 + headH + 8;

    auto section = [&](char const *key) {
        RECT r{left, y, right, y + 16};
        drawLabel(dc, r, tr(key), kDimText, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        y += 18;
    };
    auto spinRow = [&](char const *labelKey, std::wstring const &value, int idM, int idP, bool enabled, char const *tip) {
        RECT lab{left, y, right - 70, y + 18};
        drawLabel(dc, lab, tr(labelKey), enabled ? kText : kDimText, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        RECT minus{right - 66, y, right - 46, y + 18};
        RECT val{right - 46, y, right - 20, y + 18};
        RECT plus{right - 20, y, right, y + 18};
        bool overM = enabled && ptIn(minus, shell.mouseX, shell.mouseY);
        bool overP = enabled && ptIn(plus, shell.mouseX, shell.mouseY);
        drawButton(dc, minus, L"-", enabled ? (overM ? BtnState::Hover : BtnState::Normal) : BtnState::Disabled);
        drawLabel(dc, val, value.c_str(), enabled ? kText : kDimText, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
        drawButton(dc, plus, L"+", enabled ? (overP ? BtnState::Hover : BtnState::Normal) : BtnState::Disabled);
        RECT row{left, y, right, y + 18};
        addHit(w, row, TwIdle, tip);
        addHit(w, minus, idM, tip);
        addHit(w, plus, idP, tip);
        y += 22;
    };
    auto checkRow = [&](char const *key, bool on, int id, bool enabled, int indent, char const *tip) {
        RECT box{left + indent, y + 2, left + indent + 14, y + 16};
        drawCheck(dc, box, on, enabled);
        RECT lab{box.right + 6, y, right, y + 18};
        drawLabel(dc, lab, tr(key), enabled ? kText : kDimText, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        RECT hit{left + indent, y, right, y + 18};
        addHit(w, hit, enabled ? id : TwIdle, tip);
        y += 20;
    };
    auto shapeRow = [&](BrushShape cur) {
        section("tw_brush_shape");
        int bw = (right - left - 8) / 3;
        char const *keys[3] = {"tw_shape_circle", "tw_shape_square", "tw_shape_triangle"};
        int ids[3] = {TwShape0, TwShape1, TwShape2};
        for (int i = 0; i < 3; ++i) {
            RECT rc{left + i * (bw + 4), y, left + i * (bw + 4) + bw, y + 20};
            bool on = static_cast<int>(cur) == i;
            bool over = ptIn(rc, shell.mouseX, shell.mouseY);
            fillRect(dc, rc, on ? kBtnSel : (over ? kBtnHover : kBtn));
            frameRect(dc, rc, on || over ? kBorder : kBorderDim);
            drawLabel(dc, rc, tr(keys[i]), kText, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
            addHit(w, rc, ids[i], "tip_tw_brush_shape");
        }
        y += 26;
    };

    PaletteId palUi = pal;
    if (palUi == PaletteId::Erase) {
        EraseToolSettings &s = shell.tools.erase;
        spinRow("tw_brush_size", std::to_wstring(s.brushSize), TwSpin0M, TwSpin0P, true, "tip_tw_brush_size");
        section("tw_delete");
        checkRow("tw_solids", s.deleteSolids, TwSolids, true, 0, "tip_tw_solids");
        checkRow("tw_strict_delete", s.strictSolids, TwStrictSolids, s.deleteSolids, 16, "tip_tw_strict_delete");
        checkRow("tw_liquids", s.deleteLiquids, TwLiquids, true, 0, "tip_tw_liquids");
        checkRow("tw_gases", s.deleteGases, TwGases, true, 0, "tip_tw_gases");
        shapeRow(s.shape);
    } else if (pal == PaletteId::Grab) {
        GrabToolSettings &s = shell.tools.grab;
        spinRow("tw_grab_strength", shortFloat(s.strength, 2), TwSpin0M, TwSpin0P, true, "tip_tw_grab_strength");
        checkRow("tw_group_grab", s.groupGrab, TwGroup, true, 0, "tip_tw_group_grab");
        spinRow("tw_group_grab_radius", shortFloat(s.groupRadius, 0), TwSpin1M, TwSpin1P, s.groupGrab, "tip_tw_group_grab_radius");
        checkRow("tw_phantom_grab", s.phantom, TwPhantom, true, 0, "tip_tw_phantom_grab");
    } else if (pal == PaletteId::Brush) {
        BrushToolSettings &s = shell.tools.brush;
        spinRow("tw_brush_size", std::to_wstring(s.brushSize), TwSpin0M, TwSpin0P, true, "tip_tw_brush_size");
        shapeRow(s.shape);
        checkRow("tw_color_mode", s.colorMode, TwColorMode, true, 0, "tip_tw_color_mode");
        RECT dLab{left, y, right, y + 16};
        drawLabel(dc, dLab, tr("tw_dye_color"), s.colorMode ? kText : kDimText, DT_LEFT | DT_VCENTER | DT_SINGLELINE);
        addHit(w, dLab, TwIdle, "tip_tw_dye_color");
        y += 18;
        int sw = (right - left - 14) / 8;
        for (int i = 0; i < 8; ++i) {
            RECT rc{left + i * (sw + 2), y, left + i * (sw + 2) + sw, y + 16};
            fillRect(dc, rc, kDyeRgb[i]);
            bool over = ptIn(rc, shell.mouseX, shell.mouseY);
            frameRect(dc, rc, (s.colorMode && (shell.dyeSwatch == i || over)) ? kBorder : kBorderDim);
            addHit(w, rc, TwDye0 + i, "tip_tw_dye_color");
        }
        y += 22;
    } else if (pal == PaletteId::Touch) {
        TouchToolSettings &s = shell.tools.touch;
        checkRow("tw_group_touch", s.groupTouch, TwGroup, true, 0, "tip_tw_group_touch");
        spinRow("tw_group_touch_radius", shortFloat(s.groupRadius, 0), TwSpin1M, TwSpin1P, s.groupTouch, "tip_tw_group_touch_radius");
        checkRow("tw_toggle_anchor", s.toggleAnchor, TwToggleAnchor, true, 0, "tip_tw_toggle_anchor");
    } else if (pal == PaletteId::Heat) {
        HeatToolSettings &s = shell.tools.heat;
        wchar_t pbuf[24]{};
        swprintf_s(pbuf, L"%+.2f", static_cast<double>(s.power));
        spinRow("tw_heat_power", pbuf, TwSpin0M, TwSpin0P, true, "tip_tw_heat_power");
        spinRow("tw_brush_size", std::to_wstring(s.brushSize), TwSpin1M, TwSpin1P, true, "tip_tw_brush_size");
        shapeRow(s.shape);
    }

    int contentBottom = y + 6;
    if (contentBottom - y0 > w.h) {
        w.h = contentBottom - y0;
        panel.bottom = y0 + w.h;
        w.panel = panel;
        frameRect(dc, panel, kBorder);
        frameRect(dc, insetRect(panel, 1), kBorderDim);
    }
    SelectObject(dc, old);
}

SettingsMouseResult handleToolWindowMouseDown(ShellState &shell, int x, int y) {
    if (!shell.toolWin.open) return SettingsMouseResult::Miss;
    if (!ptIn(shell.toolWin.panel, x, y)) return SettingsMouseResult::Miss;

    int hit = hitAt(shell.toolWin, x, y);
    if (hit == TwClose) {
        shell.toolWin.open = false;
        shell.toolWin.userClosed = true;
        shell.toolWin.dragging = false;
        if (shell.editTool == EditTool::Brush)
            shell.brushEnabled = false;
        return SettingsMouseResult::Consume;
    }
    if (hit == TwTitle) {
        shell.toolWin.dragging = true;
        shell.toolWin.dragOX = x - shell.toolWin.x;
        shell.toolWin.dragOY = y - shell.toolWin.y;
        return SettingsMouseResult::Drag;
    }
    if (hit == TwIdle || hit == TwNone) return SettingsMouseResult::Consume;

    auto nudgeInt = [](int &v, int d, int lo, int hi) { v = std::clamp(v + d, lo, hi); };
    auto nudgeF = [](float &v, float d, float lo, float hi) { v = std::clamp(v + d, lo, hi); };

    switch (shell.editTool) {
        case EditTool::Erase: {
            EraseToolSettings &s = shell.tools.erase;
            if (hit == TwSpin0M) nudgeInt(s.brushSize, -1, 1, 14);
            if (hit == TwSpin0P) nudgeInt(s.brushSize, 1, 1, 14);
            if (hit == TwShape0) s.shape = BrushShape::Circle;
            if (hit == TwShape1) s.shape = BrushShape::Square;
            if (hit == TwShape2) s.shape = BrushShape::Triangle;
            if (hit == TwSolids) s.deleteSolids = !s.deleteSolids;
            if (hit == TwStrictSolids && s.deleteSolids) s.strictSolids = !s.strictSolids;
            if (hit == TwLiquids) s.deleteLiquids = !s.deleteLiquids;
            if (hit == TwGases) s.deleteGases = !s.deleteGases;
            break;
        }
        case EditTool::Grab: {
            GrabToolSettings &s = shell.tools.grab;
            if (hit == TwSpin0M) nudgeF(s.strength, -0.25f, 0.25f, 2.0f);
            if (hit == TwSpin0P) nudgeF(s.strength, 0.25f, 0.25f, 2.0f);
            if (hit == TwGroup) s.groupGrab = !s.groupGrab;
            if (hit == TwSpin1M && s.groupGrab) nudgeF(s.groupRadius, -1.0f, 2.0f, 24.0f);
            if (hit == TwSpin1P && s.groupGrab) nudgeF(s.groupRadius, 1.0f, 2.0f, 24.0f);
            if (hit == TwPhantom) s.phantom = !s.phantom;
            break;
        }
        case EditTool::Brush: {
            BrushToolSettings &s = shell.tools.brush;
            if (hit == TwSpin0M) nudgeInt(s.brushSize, -1, 1, 14);
            if (hit == TwSpin0P) nudgeInt(s.brushSize, 1, 1, 14);
            if (hit == TwShape0) s.shape = BrushShape::Circle;
            if (hit == TwShape1) s.shape = BrushShape::Square;
            if (hit == TwShape2) s.shape = BrushShape::Triangle;
            if (hit == TwColorMode) s.colorMode = !s.colorMode;
            if (hit >= TwDye0 && hit <= TwDye7 && s.colorMode)
                shell.dyeSwatch = hit - TwDye0;
            break;
        }
        case EditTool::Touch: {
            TouchToolSettings &s = shell.tools.touch;
            if (hit == TwGroup) s.groupTouch = !s.groupTouch;
            if (hit == TwSpin1M && s.groupTouch) nudgeF(s.groupRadius, -1.0f, 2.0f, 24.0f);
            if (hit == TwSpin1P && s.groupTouch) nudgeF(s.groupRadius, 1.0f, 2.0f, 24.0f);
            if (hit == TwToggleAnchor) s.toggleAnchor = !s.toggleAnchor;
            break;
        }
        case EditTool::Heat: {
            HeatToolSettings &s = shell.tools.heat;
            if (hit == TwSpin0M) nudgeF(s.power, -0.25f, -8.0f, 8.0f);
            if (hit == TwSpin0P) nudgeF(s.power, 0.25f, -8.0f, 8.0f);
            if (hit == TwSpin1M) nudgeInt(s.brushSize, -1, 1, 14);
            if (hit == TwSpin1P) nudgeInt(s.brushSize, 1, 1, 14);
            if (hit == TwShape0) s.shape = BrushShape::Circle;
            if (hit == TwShape1) s.shape = BrushShape::Square;
            if (hit == TwShape2) s.shape = BrushShape::Triangle;
            break;
        }
        default: break;
    }
    return SettingsMouseResult::Consume;
}

void dragToolWindow(ShellState &shell, int x, int y) {
    if (!shell.toolWin.dragging) return;
    shell.toolWin.x = x - shell.toolWin.dragOX;
    shell.toolWin.y = y - shell.toolWin.dragOY;
    clampToolWindow(shell);
}

bool handleToolWindowWheel(ShellState &shell, int x, int y, int wheelDelta) {
    if (!toolWindowCoversPoint(shell, x, y)) return false;
    int d = wheelDelta > 0 ? 1 : -1;
    switch (shell.editTool) {
        case EditTool::Erase:
            shell.tools.erase.brushSize = std::clamp(shell.tools.erase.brushSize + d, 1, 14);
            return true;
        case EditTool::Brush:
            shell.tools.brush.brushSize = std::clamp(shell.tools.brush.brushSize + d, 1, 14);
            return true;
        case EditTool::Heat:
            shell.tools.heat.brushSize = std::clamp(shell.tools.heat.brushSize + d, 1, 14);
            return true;
        case EditTool::Grab:
            shell.tools.grab.strength = std::clamp(shell.tools.grab.strength + d * 0.25f, 0.25f, 2.0f);
            return true;
        case EditTool::Touch:
            if (shell.tools.touch.groupTouch)
                shell.tools.touch.groupRadius = std::clamp(shell.tools.touch.groupRadius + static_cast<float>(d), 2.0f, 24.0f);
            return true;
        default:
            return true;
    }
}

} // namespace ui
