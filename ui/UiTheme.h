#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <algorithm>
#include <string>

namespace ui {

constexpr COLORREF kBg       = RGB(50, 50, 50);
constexpr COLORREF kPanel    = RGB(42, 42, 42);
constexpr COLORREF kBtn      = RGB(90, 90, 90);
constexpr COLORREF kBtnHover = RGB(112, 112, 112);
constexpr COLORREF kBtnDown  = RGB(70, 70, 70);
constexpr COLORREF kBtnSel   = RGB(120, 120, 120);
constexpr COLORREF kBorder   = RGB(230, 230, 230);
constexpr COLORREF kBorderDim= RGB(160, 160, 160);
constexpr COLORREF kText     = RGB(245, 245, 245);
constexpr COLORREF kDimText  = RGB(168, 168, 168);
constexpr COLORREF kBlack    = RGB(0, 0, 0);
constexpr COLORREF kSlot     = RGB(72, 72, 72);
constexpr COLORREF kPauseOn  = RGB(120, 90, 48);
constexpr COLORREF kClear    = RGB(120, 56, 60);
constexpr COLORREF kSettingsOn = RGB(70, 90, 110);

inline void fillRect(HDC dc, RECT const &rc, COLORREF color) {
    HBRUSH b = CreateSolidBrush(color);
    FillRect(dc, &rc, b);
    DeleteObject(b);
}

inline void frameRect(HDC dc, RECT const &rc, COLORREF color) {
    HBRUSH b = CreateSolidBrush(color);
    FrameRect(dc, &rc, b);
    DeleteObject(b);
}

inline RECT insetRect(RECT rc, int m) {
    rc.left += m; rc.top += m; rc.right -= m; rc.bottom -= m;
    return rc;
}

inline bool ptIn(RECT const &rc, int x, int y) {
    POINT p{x, y};
    return PtInRect(&rc, p) != 0;
}

enum class BtnState { Normal, Hover, Selected, Disabled };

constexpr COLORREF kViewGeneral = RGB(72, 72, 72);
constexpr COLORREF kViewGeneralSel = RGB(108, 108, 108);
constexpr COLORREF kViewLiquid = RGB(48, 92, 128);
constexpr COLORREF kViewLiquidSel = RGB(62, 122, 168);
constexpr COLORREF kViewSolid = RGB(158, 96, 52);
constexpr COLORREF kViewSolidSel = RGB(198, 64, 58);
constexpr COLORREF kViewGas = RGB(132, 136, 142);
constexpr COLORREF kViewGasSel = RGB(168, 172, 178);

inline void drawButton(HDC dc, RECT const &rc, wchar_t const *label, BtnState state,
    COLORREF fill = kBtn, COLORREF selFill = kBtnSel) {
    COLORREF bg = fill;
    if (state == BtnState::Hover) {
        int r = GetRValue(fill), g = GetGValue(fill), b = GetBValue(fill);
        bg = RGB(std::min(255, r + 30), std::min(255, g + 30), std::min(255, b + 30));
    }
    if (state == BtnState::Selected) bg = selFill;
    if (state == BtnState::Disabled) bg = RGB(64, 64, 64);
    fillRect(dc, rc, bg);
    frameRect(dc, rc, state == BtnState::Selected ? kBorder : kBorderDim);
    COLORREF old = SetTextColor(dc, state == BtnState::Disabled ? kDimText : kText);
    int oldBk = SetBkMode(dc, TRANSPARENT);
    DrawTextW(dc, label, -1, const_cast<RECT *>(&rc), DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
    SetTextColor(dc, old);
    SetBkMode(dc, oldBk);
}

inline void drawPanel(HDC dc, RECT const &rc, COLORREF fill = kPanel) {
    fillRect(dc, rc, fill);
    frameRect(dc, rc, kBorder);
}

inline void drawLabel(HDC dc, RECT rc, wchar_t const *text, COLORREF color = kText, UINT fmt = DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS) {
    COLORREF old = SetTextColor(dc, color);
    int oldBk = SetBkMode(dc, TRANSPARENT);
    DrawTextW(dc, text, -1, &rc, fmt);
    SetTextColor(dc, old);
    SetBkMode(dc, oldBk);
}

inline std::wstring shortFloat(double value, int decimals = 2) {
    wchar_t buffer[64]{};
    swprintf_s(buffer, L"%.*f", decimals, value);
    return buffer;
}

} // namespace ui
