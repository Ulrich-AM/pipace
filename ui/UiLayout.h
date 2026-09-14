#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <algorithm>

namespace ui {

inline RECT insetLike(RECT rc, int m) {
    rc.left += m; rc.top += m; rc.right -= m; rc.bottom -= m;
    return rc;
}

struct Layout {
    RECT client{};
    RECT topBar{};
    RECT speedMinus{};
    RECT speedLabel{};
    RECT speedPlus{};
    RECT viewBtn[13]{};
    int viewCount = 13;
    RECT settingsBtn{};
    RECT pauseBtn{};
    RECT clearBtn{};
    RECT canvas{};
    RECT canvasInner{};
    RECT infoBar{};
    RECT tray{};
    RECT traySlot[8]{};
    int catCount = 0;
    RECT console{};
    RECT sidebar{};
    RECT search{};
    RECT catCol{};
    RECT catBody{};
    RECT propsCol{};
    RECT propsHeader{};
    RECT propsBody{};
    RECT inspectCol{};
    RECT inspectHeader{};
    RECT inspectBody{};
    RECT propBrushMinus{};
    RECT propBrushPlus{};
    RECT propAnchored{};
    RECT propSleeping{};
    RECT propPowder{};
    RECT propGrainMinus{};
    RECT propGrainPlus{};
    RECT propPowerMinus{};
    RECT propPowerPlus{};
    RECT propDyeSwatch[8]{};
    RECT propDyeOnly{};
    RECT propDyeStrMinus{};
    RECT propDyeStrPlus{};

    int palCount = 0;
    RECT palSlot[8]{};
};

inline Layout computeLayout(int clientW, int clientH) {
    Layout L{};
    L.client = RECT{0, 0, clientW, clientH};
    int const pad = 8;
    int const gap = 6;
    int const topH = 32;
    int const infoH = 34;
    int const trayH = 38;
    int const consH = 22;
    int sidebarW = std::clamp(clientW * 38 / 100, 360, 520);
    int leftW = std::max(280, clientW - sidebarW - pad);

    int x = pad;
    int y = pad;
    L.topBar = RECT{x, y, x + leftW - pad, y + topH};

    int bx = L.topBar.left;
    int by = L.topBar.top;
    int bh = topH;
    L.speedMinus = RECT{bx, by, bx + 32, by + bh};
    L.speedLabel = RECT{bx + 36, by, bx + 92, by + bh};
    L.speedPlus  = RECT{bx + 96, by, bx + 128, by + bh};
    int right = L.topBar.right;
    L.clearBtn    = RECT{right - 80, by, right, by + bh};
    L.pauseBtn    = RECT{right - 168, by, right - 86, by + bh};
    L.settingsBtn = RECT{right - 270, by, right - 174, by + bh};

    L.viewCount = 13;
    {
        int viewLeft = L.speedPlus.right + 8;
        int viewRight = L.settingsBtn.left - 8;
        int avail = std::max(0, viewRight - viewLeft);
        int groupGap = 5;
        int innerGap = 2;
        int extra = 3 * groupGap;
        int innerGaps = 9;
        int btnW = avail > 0 ? (avail - extra - innerGaps * innerGap) / L.viewCount : 26;
        btnW = std::clamp(btnW, 20, 44);
        int groupOf[13] = {0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3};
        int x = viewLeft;
        for (int i = 0; i < L.viewCount; ++i) {
            if (i > 0) x += (groupOf[i] != groupOf[i - 1]) ? groupGap : innerGap;
            L.viewBtn[i] = RECT{x, by, x + btnW, by + bh};
            x = L.viewBtn[i].right;
        }
    }

    y = L.topBar.bottom + gap;
    int bottomStack = infoH + gap + trayH + gap + consH + pad;
    int canvasBottom = std::max(y + 80, clientH - bottomStack);
    L.canvas = RECT{x, y, x + leftW - pad, canvasBottom};
    L.canvasInner = insetLike(L.canvas, 2);

    y = L.canvas.bottom + gap;
    L.infoBar = RECT{x, y, x + leftW - pad, y + infoH};
    y = L.infoBar.bottom + gap;
    L.tray = RECT{x, y, x + leftW - pad, y + trayH};
    L.catCount = 7;
    {
        int inner = static_cast<int>(L.tray.right - L.tray.left);
        int slotW = std::max(1, inner / L.catCount);
        for (int i = 0; i < L.catCount; ++i) {
            L.traySlot[i] = RECT{L.tray.left + i * slotW, L.tray.top, (i + 1 == L.catCount) ? L.tray.right : L.tray.left + (i + 1) * slotW, L.tray.bottom};
        }
    }
    y = L.tray.bottom + gap;
    L.console = RECT{x, y, x + leftW - pad, y + consH};

    int sx = clientW - sidebarW;
    L.sidebar = RECT{sx, pad, clientW - pad, clientH - pad};
    L.search = RECT{L.sidebar.left + 6, L.sidebar.top + 6, L.sidebar.right - 6, L.sidebar.top + 32};

    int colsTop = L.search.bottom + 8;
    int catW = 78;
    int rest = L.sidebar.right - 6 - (L.sidebar.left + 6 + catW);
    int propW = rest / 2;
    L.catCol = RECT{L.sidebar.left + 6, colsTop, L.sidebar.left + 6 + catW, L.sidebar.bottom - 6};
    L.catBody = L.catCol;
    L.propsCol = RECT{L.catCol.right + 6, colsTop, L.catCol.right + 6 + propW, L.sidebar.bottom - 6};
    L.inspectCol = RECT{L.propsCol.right + 6, colsTop, L.sidebar.right - 6, L.sidebar.bottom - 6};

    int headH = 28;
    L.propsHeader   = RECT{L.propsCol.left, L.propsCol.top, L.propsCol.right, L.propsCol.top + headH};
    L.propsBody     = RECT{L.propsCol.left, L.propsHeader.bottom + 4, L.propsCol.right, L.propsCol.bottom};
    L.inspectHeader = RECT{L.inspectCol.left, L.inspectCol.top, L.inspectCol.right, L.inspectCol.top + headH};
    L.inspectBody   = RECT{L.inspectCol.left, L.inspectHeader.bottom + 4, L.inspectCol.right, L.inspectCol.bottom};

    int px = L.propsBody.left + 8;
    int py = L.propsBody.top + 54;
    L.propBrushMinus = RECT{px, py, px + 28, py + 22};
    L.propBrushPlus  = RECT{px + 34, py, px + 62, py + 22};
    L.propPowerMinus = RECT{px, py + 52, px + 28, py + 74};
    L.propPowerPlus  = RECT{px + 34, py + 52, px + 62, py + 74};
    L.propAnchored   = RECT{px, py + 52, px + 18, py + 70};
    L.propSleeping   = RECT{px, py + 100, px + 18, py + 118};
    L.propPowder     = RECT{px, py + 148, px + 18, py + 166};
    L.propGrainMinus = RECT{px, py + 190, px + 28, py + 212};
    L.propGrainPlus  = RECT{px + 34, py + 190, px + 62, py + 212};
    int dpy = py + 88;
    for (int i = 0; i < 8; ++i) {
        int col = i % 4, row = i / 4;
        int dx = px + col * 22;
        int dy = dpy + row * 22;
        L.propDyeSwatch[i] = RECT{dx, dy, dx + 18, dy + 18};
    }
    L.propDyeOnly = RECT{px, dpy + 70, px + 18, dpy + 88};
    L.propDyeStrMinus = RECT{px, dpy + 114, px + 28, dpy + 136};
    L.propDyeStrPlus  = RECT{px + 34, dpy + 114, px + 62, dpy + 136};

    return L;
}

} // namespace ui
