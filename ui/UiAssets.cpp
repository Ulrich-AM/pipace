#include "ui/UiAssets.h"

#include <objidl.h>
#include <gdiplus.h>

#include <algorithm>
#include <string>

namespace ui {
namespace {

ULONG_PTR gGdiplusToken = 0;
bool gGdiplusOk = false;
Gdiplus::Bitmap *gCatIcon[7]{};

wchar_t const *kCatFiles[7] = {
    L"tools.png",
    L"liquid.png",
    L"solid.png",
    L"gas.png",
    L"plasma.png",
    L"energy.png",
    L"misc.png",
};

std::wstring exeDir() {
    wchar_t buf[MAX_PATH]{};
    DWORD n = GetModuleFileNameW(nullptr, buf, MAX_PATH);
    std::wstring path(buf, n ? n : 0);
    size_t slash = path.find_last_of(L"\\/");
    if (slash != std::wstring::npos) path.resize(slash);
    return path;
}

bool fileExists(std::wstring const &path) {
    DWORD attr = GetFileAttributesW(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY) == 0;
}

std::wstring findAsset(wchar_t const *file) {
    std::wstring dir = exeDir();
    std::wstring candidates[] = {
        dir + L"\\assets\\" + file,
        dir + L"\\..\\assets\\" + file,
        std::wstring(L"assets\\") + file,
    };
    for (std::wstring const &p : candidates)
        if (fileExists(p)) return p;
    return {};
}

} // namespace

void initAssets() {
    if (gGdiplusOk) return;
    Gdiplus::GdiplusStartupInput input;
    if (Gdiplus::GdiplusStartup(&gGdiplusToken, &input, nullptr) != Gdiplus::Ok) return;
    gGdiplusOk = true;
    for (int i = 0; i < 7; ++i) {
        std::wstring path = findAsset(kCatFiles[i]);
        if (path.empty()) continue;
        Gdiplus::Bitmap *bmp = Gdiplus::Bitmap::FromFile(path.c_str(), FALSE);
        if (!bmp || bmp->GetLastStatus() != Gdiplus::Ok) {
            delete bmp;
            continue;
        }
        gCatIcon[i] = bmp;
    }
}

void shutdownAssets() {
    for (int i = 0; i < 7; ++i) {
        delete gCatIcon[i];
        gCatIcon[i] = nullptr;
    }
    if (gGdiplusOk) {
        Gdiplus::GdiplusShutdown(gGdiplusToken);
        gGdiplusToken = 0;
        gGdiplusOk = false;
    }
}

int drawCategoryIcon(HDC dc, RECT const &slot, int categoryIndex, bool dimmed) {
    if (!gGdiplusOk || categoryIndex < 0 || categoryIndex >= 7) return 0;
    Gdiplus::Bitmap *bmp = gCatIcon[categoryIndex];
    if (!bmp) return 0;
    int slotH = slot.bottom - slot.top;
    int slotW = slot.right - slot.left;
    int icon = std::clamp(slotH - 12, 14, 22);
    if (slotW < icon + 10) icon = std::max(12, slotW - 8);
    int x = slot.left + 6;
    int y = slot.top + (slotH - icon) / 2;
    Gdiplus::Graphics g(dc);
    g.SetInterpolationMode(Gdiplus::InterpolationModeNearestNeighbor);
    g.SetPixelOffsetMode(Gdiplus::PixelOffsetModeHalf);
    g.SetSmoothingMode(Gdiplus::SmoothingModeNone);
    Gdiplus::Rect dest(x, y, icon, icon);
    if (dimmed) {
        Gdiplus::ColorMatrix cm{};
        cm.m[0][0] = cm.m[1][1] = cm.m[2][2] = cm.m[4][4] = 1.0f;
        cm.m[3][3] = 0.38f;
        Gdiplus::ImageAttributes attr;
        attr.SetColorMatrix(&cm);
        g.DrawImage(bmp, dest, 0, 0, static_cast<INT>(bmp->GetWidth()), static_cast<INT>(bmp->GetHeight()),
            Gdiplus::UnitPixel, &attr);
    } else {
        g.DrawImage(bmp, dest, 0, 0, static_cast<INT>(bmp->GetWidth()), static_cast<INT>(bmp->GetHeight()),
            Gdiplus::UnitPixel);
    }
    return icon + 10;
}

} // namespace ui
