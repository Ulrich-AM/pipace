#include "ui/UiShell.h"
#include "thermal/ThermalEngine.h"

#include "ui/UiAssets.h"
#include "ui/UiLanguage.h"
#include "ui/UiTheme.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cwchar>

namespace ui {
namespace {

BtnState btnState(ShellState const &shell, HitId id, bool selected) {
    if (selected) return BtnState::Selected;
    if (shell.hoverId == static_cast<int>(id) && shell.mouseDown) return BtnState::Normal;
    if (shell.hoverId == static_cast<int>(id)) return BtnState::Hover;
    return BtnState::Normal;
}

void layoutPalette(ShellState &shell) {
    RECT body = shell.layout.catBody;
    int n = shell.elementCount();
    shell.layout.palCount = n;
    int slotH = 40;
    for (int i = 0; i < 8; ++i) shell.layout.palSlot[i] = RECT{};
    for (int i = 0; i < n; ++i) {
        int y = body.top + i * (slotH + 4);
        if (y + slotH > body.bottom) break;
        shell.layout.palSlot[i] = RECT{body.left, y, body.right, y + slotH};
    }
}

struct DyeDef { float r, g, b; COLORREF vis; };
DyeDef const kDyeSwatches[8] = {
    {0.00f, 0.00f, 0.00f, RGB(210, 210, 210)},
    {0.92f, 0.12f, 0.14f, RGB(220, 40, 42)},
    {0.95f, 0.45f, 0.08f, RGB(230, 118, 28)},
    {0.95f, 0.82f, 0.12f, RGB(230, 200, 36)},
    {0.12f, 0.72f, 0.22f, RGB(40, 170, 55)},
    {0.10f, 0.28f, 0.92f, RGB(36, 72, 220)},
    {0.55f, 0.18f, 0.78f, RGB(140, 50, 190)},
    {0.08f, 0.08f, 0.10f, RGB(28, 28, 32)}
};

const wchar_t *debugName(DebugView v) {
    switch (v) {
        case DebugView::Fill: return tr("bar_fill");
        case DebugView::Pressure: return tr("bar_liqp");
        case DebugView::Velocity: return tr("bar_lvel");
        case DebugView::Divergence: return tr("bar_ldiv");
        case DebugView::Chunks: return tr("bar_chnk");
        case DebugView::Rigid: return tr("bar_rgdn");
        case DebugView::GasPressure: return tr("bar_gasp");
        case DebugView::GasAmount: return tr("bar_gasa");
        case DebugView::GasVelocity: return tr("bar_gasv");
        case DebugView::Temperature: return tr("bar_temp");
        default: return tr("bar_norm");
    }
}

struct ViewBarDef {
    HitId hit;
    char const *labelKey;
    char const *titleKey;
    char const *helpKey;
    COLORREF fill;
    COLORREF sel;
};

ViewBarDef const kViewBar[] = {
    {HitId::ViewNorm, "bar_norm", "view_title_norm", "view_help_norm", kViewGeneral, kViewGeneralSel},
    {HitId::ViewChnk, "bar_chnk", "view_title_chnk", "view_help_chnk", kViewGeneral, kViewGeneralSel},
    {HitId::ViewFill, "bar_fill", "view_title_fill", "view_help_fill", kViewLiquid, kViewLiquidSel},
    {HitId::ViewLiqp, "bar_liqp", "view_title_liqp", "view_help_liqp", kViewLiquid, kViewLiquidSel},
    {HitId::ViewLvel, "bar_lvel", "view_title_lvel", "view_help_lvel", kViewLiquid, kViewLiquidSel},
        {HitId::ViewLdiv, "bar_ldiv", "view_title_ldiv", "view_help_ldiv", kViewLiquid, kViewLiquidSel},
        {HitId::ViewTemp, "bar_temp", "view_title_temp", "view_help_temp", kViewLiquid, kViewLiquidSel},
        {HitId::ViewRgdn, "bar_rgdn", "view_title_rgdn", "view_help_rgdn", kViewSolid, kViewSolidSel},
    {HitId::ViewRgdo, "bar_rgdo", "view_title_rgdo", "view_help_rgdo", kViewSolid, kViewSolidSel},
    {HitId::ViewGasp, "bar_gasp", "view_title_gasp", "view_help_gasp", kViewGas, kViewGasSel},
    {HitId::ViewGasa, "bar_gasa", "view_title_gasa", "view_help_gasa", kViewGas, kViewGasSel},
    {HitId::ViewGasv, "bar_gasv", "view_title_gasv", "view_help_gasv", kViewGas, kViewGasSel},
};

constexpr int kViewBarCount = static_cast<int>(sizeof(kViewBar) / sizeof(kViewBar[0]));

ViewBarDef const *viewBarByHit(int hoverId) {
    for (int i = 0; i < kViewBarCount; ++i)
        if (static_cast<int>(kViewBar[i].hit) == hoverId) return &kViewBar[i];
    return nullptr;
}

bool viewBarSelected(ViewBarDef const &item, View const &view) {
    switch (item.hit) {
        case HitId::ViewNorm: return view.debugView == DebugView::Normal;
        case HitId::ViewChnk: return view.debugView == DebugView::Chunks;
        case HitId::ViewFill: return view.debugView == DebugView::Fill;
        case HitId::ViewLiqp: return view.debugView == DebugView::Pressure;
        case HitId::ViewLvel: return view.debugView == DebugView::Velocity;
        case HitId::ViewLdiv: return view.debugView == DebugView::Divergence;
        case HitId::ViewTemp: return view.debugView == DebugView::Temperature;
        case HitId::ViewRgdn: return view.debugView == DebugView::Rigid;
        case HitId::ViewRgdo: return view.rigid && view.rigid->debugOverlay;
        case HitId::ViewGasp: return view.debugView == DebugView::GasPressure;
        case HitId::ViewGasa: return view.debugView == DebugView::GasAmount;
        case HitId::ViewGasv: return view.debugView == DebugView::GasVelocity;
        default: return false;
    }
}

} // namespace

void ShellState::ensureFonts() {
    initAssets();
    if (!uiFont)
        uiFont = CreateFontW(16, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, NONANTIALIASED_QUALITY, DEFAULT_PITCH, L"Tahoma");
    if (!smallFont)
        smallFont = CreateFontW(13, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, NONANTIALIASED_QUALITY, DEFAULT_PITCH, L"Tahoma");
    if (!consoleFont)
        consoleFont = CreateFontW(14, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, NONANTIALIASED_QUALITY, FIXED_PITCH, L"Consolas");
}

void ShellState::releaseFonts() {
    if (uiFont) { DeleteObject(uiFont); uiFont = nullptr; }
    if (smallFont) { DeleteObject(smallFont); smallFont = nullptr; }
    if (consoleFont) { DeleteObject(consoleFont); consoleFont = nullptr; }
    shutdownAssets();
}

void ShellState::log(wchar_t const *line) {
    consoleLines.emplace_back(line);
    if (consoleLines.size() > 80) consoleLines.erase(consoleLines.begin(), consoleLines.begin() + 20);
}

void ShellState::applyPalette(PaletteId id, Tool &tool, MaterialId &drawMaterial) {
    palette = id;
    switch (id) {
        case PaletteId::None: break;
        case PaletteId::Water: category = Category::Fluids; tool = Tool::Water; break;
        case PaletteId::Honey: category = Category::Fluids; tool = Tool::Water; break;
        case PaletteId::Wood:  category = Category::Solids; tool = Tool::Rigid; drawMaterial = MATERIAL_WOOD; break;
        case PaletteId::Stone: category = Category::Solids; tool = Tool::Rigid; drawMaterial = MATERIAL_STONE; break;
        case PaletteId::Glass: category = Category::Solids; tool = Tool::Rigid; drawMaterial = MATERIAL_GLASS; break;
        case PaletteId::Metal: category = Category::Solids; tool = Tool::Rigid; drawMaterial = MATERIAL_METAL; break;
        case PaletteId::Erase: category = Category::Tools; tool = Tool::Eraser; break;
        case PaletteId::Grab:  category = Category::Tools; tool = Tool::Grab; break;
        case PaletteId::Wall:  category = Category::Misc; tool = Tool::Solid; break;
        case PaletteId::Heat:  category = Category::Energy; tool = Tool::Heat; break;
        case PaletteId::Cool:  category = Category::Energy; tool = Tool::Cool; break;
        case PaletteId::Pressurize: category = Category::Energy; tool = Tool::Pressurize; break;
        case PaletteId::Depressurize: category = Category::Energy; tool = Tool::Depressurize; break;
    }
}

void ShellState::applyCategory(Category cat, Tool &tool, MaterialId &drawMaterial) {
    category = cat;
    if (elementCount() <= 0) {
        palette = PaletteId::None;
        return;
    }
    applyPalette(elementAt(0), tool, drawMaterial);
}

int ShellState::elementCount() const {
    switch (category) {
        case Category::Tools: return 2;
        case Category::Fluids: return 2;
        case Category::Solids: return 4;
        case Category::Misc: return 1;
        case Category::Energy: return 4;
        default: return 0;
    }
}

PaletteId ShellState::elementAt(int slot) const {
    switch (category) {
        case Category::Tools: return slot == 0 ? PaletteId::Erase : PaletteId::Grab;
        case Category::Fluids: return slot == 1 ? PaletteId::Honey : PaletteId::Water;
        case Category::Solids:
            if (slot == 1) return PaletteId::Stone;
            if (slot == 2) return PaletteId::Glass;
            if (slot == 3) return PaletteId::Metal;
            return PaletteId::Wood;
        case Category::Misc: return PaletteId::Wall;
        case Category::Energy:
            if (slot == 1) return PaletteId::Cool;
            if (slot == 2) return PaletteId::Pressurize;
            if (slot == 3) return PaletteId::Depressurize;
            return PaletteId::Heat;
        default: return PaletteId::None;
    }
}

bool ShellState::hasPlacement() const {
    return palette != PaletteId::None && palette != PaletteId::Grab;
}

wchar_t const *ShellState::categoryName(Category cat) const {
    switch (cat) {
        case Category::Tools: return tr("cat_tools");
        case Category::Fluids: return tr("cat_liquids");
        case Category::Solids: return tr("cat_solids");
        case Category::Gases: return tr("cat_gases");
        case Category::Plasma: return tr("cat_plasma");
        case Category::Energy: return tr("cat_energy");
        case Category::Misc: return tr("cat_misc");
        default: return L"";
    }
}

wchar_t const *ShellState::paletteName(PaletteId id) const {
    switch (id) {
        case PaletteId::Water: return tr("el_water");
        case PaletteId::Honey: return tr("el_honey");
        case PaletteId::Wood:  return tr("el_wood");
        case PaletteId::Stone: return tr("el_stone");
        case PaletteId::Glass: return tr("el_glass");
        case PaletteId::Metal: return tr("el_metal");
        case PaletteId::Erase: return tr("el_erase");
        case PaletteId::Grab:  return tr("el_grab");
        case PaletteId::Wall:  return tr("el_wall");
        case PaletteId::Heat:  return tr("el_heat");
        case PaletteId::Cool:  return tr("el_cool");
        case PaletteId::Pressurize: return tr("el_pressurize");
        case PaletteId::Depressurize: return tr("el_depressurize");
        case PaletteId::None:  return L"";
    }
    return L"";
}

wchar_t const *ShellState::paletteHint(PaletteId id) const {
    switch (id) {
        case PaletteId::Water: return tr("hint_water");
        case PaletteId::Honey: return tr("hint_honey");
        case PaletteId::Wood:  return tr("hint_wood");
        case PaletteId::Stone: return tr("hint_stone");
        case PaletteId::Glass: return tr("hint_glass");
        case PaletteId::Metal: return tr("hint_metal");
        case PaletteId::Erase: return tr("hint_erase");
        case PaletteId::Grab:  return tr("hint_grab");
        case PaletteId::Wall:  return tr("hint_wall");
        case PaletteId::Heat:  return tr("hint_heat");
        case PaletteId::Cool:  return tr("hint_cool");
        case PaletteId::Pressurize: return tr("hint_pressurize");
        case PaletteId::Depressurize: return tr("hint_depressurize");
        case PaletteId::None:  return tr("hint_none");
    }
    return L"";
}

void ShellState::dyeChannels(float &r, float &g, float &b) const {
    int i = std::clamp(dyeSwatch, 0, 7);
    r = kDyeSwatches[i].r;
    g = kDyeSwatches[i].g;
    b = kDyeSwatches[i].b;
}

void computeShellLayout(ShellState &shell, int clientW, int clientH) {
    shell.layout = computeLayout(clientW, clientH);
    layoutPalette(shell);
}

HitId hitTest(ShellState const &shell, int x, int y, bool settingsOpen) {
    Layout const &L = shell.layout;
    auto check = [&](RECT const &rc, HitId id) {
        return ptIn(rc, x, y) ? id : HitId::None;
    };
    HitId id = HitId::None;
    if ((id = check(L.speedMinus, HitId::SpeedMinus)) != HitId::None) return id;
    if ((id = check(L.speedPlus, HitId::SpeedPlus)) != HitId::None) return id;
    for (int i = 0; i < L.viewCount; ++i) {
        if ((id = check(L.viewBtn[i], kViewBar[i].hit)) != HitId::None) return id;
    }
    if ((id = check(L.settingsBtn, HitId::Settings)) != HitId::None) return id;
    if ((id = check(L.pauseBtn, HitId::Pause)) != HitId::None) return id;
    if ((id = check(L.clearBtn, HitId::Clear)) != HitId::None) return id;
    if ((id = check(L.search, HitId::Search)) != HitId::None) return id;
    if ((id = check(L.console, HitId::Console)) != HitId::None) return id;
    if (!settingsOpen && (id = check(L.propBrushMinus, HitId::BrushMinus)) != HitId::None) return id;
    if (!settingsOpen && (id = check(L.propBrushPlus, HitId::BrushPlus)) != HitId::None) return id;
    if (!settingsOpen && shell.category == Category::Energy) {
        if ((id = check(L.propPowerMinus, HitId::PowerMinus)) != HitId::None) return id;
        if ((id = check(L.propPowerPlus, HitId::PowerPlus)) != HitId::None) return id;
    }
    if (!settingsOpen && shell.category == Category::Fluids) {
        for (int i = 0; i < 8; ++i) {
            if ((id = check(L.propDyeSwatch[i], static_cast<HitId>(static_cast<int>(HitId::Dye0) + i))) != HitId::None)
                return id;
        }
        RECT dyeHit = L.propDyeOnly;
        dyeHit.left = L.propsCol.left + 4;
        dyeHit.right = L.propsCol.right - 4;
        dyeHit.top -= 18;
        if ((id = check(dyeHit, HitId::DyeOnly)) != HitId::None) return id;
        if ((id = check(L.propDyeStrMinus, HitId::DyeStrMinus)) != HitId::None) return id;
        if ((id = check(L.propDyeStrPlus, HitId::DyeStrPlus)) != HitId::None) return id;
    }
    if (!settingsOpen && shell.category == Category::Solids) {
        RECT aHit = L.propAnchored;
        aHit.left = L.propsCol.left + 4;
        aHit.right = L.propsCol.right - 4;
        aHit.top -= 18;
        if ((id = check(aHit, HitId::Anchored)) != HitId::None) return id;
        RECT sHit = L.propSleeping;
        sHit.left = L.propsCol.left + 4;
        sHit.right = L.propsCol.right - 4;
        sHit.top -= 18;
        if ((id = check(sHit, HitId::Sleeping)) != HitId::None) return id;
        RECT pHit = L.propPowder;
        pHit.left = L.propsCol.left + 4;
        pHit.right = L.propsCol.right - 4;
        pHit.top -= 18;
        if ((id = check(pHit, HitId::Powder)) != HitId::None) return id;
        if ((id = check(L.propGrainMinus, HitId::GrainMinus)) != HitId::None) return id;
        if ((id = check(L.propGrainPlus, HitId::GrainPlus)) != HitId::None) return id;
    }
    for (int i = 0; i < shell.layout.catCount; ++i) {
        if (ptIn(shell.layout.traySlot[i], x, y))
            return static_cast<HitId>(static_cast<int>(HitId::Cat0) + i);
    }
    for (int i = 0; i < shell.layout.palCount; ++i) {
        if (ptIn(shell.layout.palSlot[i], x, y))
            return static_cast<HitId>(static_cast<int>(HitId::Pal0) + i);
    }
    return HitId::None;
}

void drawShell(HDC dc, ShellState &shell, View const &view) {
    shell.ensureFonts();
    Layout const &L = shell.layout;
    fillRect(dc, L.client, kBg);

    HFONT old = static_cast<HFONT>(SelectObject(dc, shell.uiFont));
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, kText);

    drawButton(dc, L.speedMinus, L"-", btnState(shell, HitId::SpeedMinus, false));
    fillRect(dc, L.speedLabel, kPanel);
    frameRect(dc, L.speedLabel, kBorderDim);
    std::wstring speed = shortFloat(view.speedValue, view.speedValue < 1.0f ? 2 : 2) + L"x";
    drawLabel(dc, L.speedLabel, speed.c_str(), kText, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    drawButton(dc, L.speedPlus, L"+", btnState(shell, HitId::SpeedPlus, false));

    SelectObject(dc, shell.smallFont);
    for (int i = 0; i < kViewBarCount && i < L.viewCount; ++i) {
        ViewBarDef const &item = kViewBar[i];
        drawButton(dc, L.viewBtn[i], tr(item.labelKey),
            btnState(shell, item.hit, viewBarSelected(item, view)),
            item.fill, item.sel);
    }
    SelectObject(dc, shell.uiFont);

    drawButton(dc, L.settingsBtn, tr("settings"), btnState(shell, HitId::Settings, view.settingsOpen), kBtn, kSettingsOn);
    drawButton(dc, L.pauseBtn, view.paused ? tr("resume") : tr("pause"), btnState(shell, HitId::Pause, view.paused), kBtn, kPauseOn);
    drawButton(dc, L.clearBtn, tr("clear"), btnState(shell, HitId::Clear, false), kClear, kClear);

    fillRect(dc, L.canvas, kBlack);
    frameRect(dc, L.canvas, kBorder);

    fillRect(dc, L.infoBar, kPanel);
    frameRect(dc, L.infoBar, kBorder);
    std::wstring info = trf("info_brush", std::to_wstring(view.brushRadius));
    if (view.engine && FluidEngine::inside(view.hoverX, view.hoverY)) {
        float fill = view.engine->fill[static_cast<size_t>(FluidEngine::ci(view.hoverX, view.hoverY))];
        info += L"    " + trf("info_fill", shortFloat(fill * 100.0, 0));
        info += L"    " + trf("info_cell", std::to_wstring(view.hoverX), std::to_wstring(view.hoverY));
    }
    if (shell.hasPlacement()) {
        info += L"    ";
        info += tr("info_line");
    }
    info += L"    " + trf("info_view", debugName(view.debugView));
    RECT infoPad = insetLike(L.infoBar, 8);
    drawLabel(dc, infoPad, info.c_str(), kText);

    fillRect(dc, L.tray, kSlot);
    frameRect(dc, L.tray, kBorder);
    SelectObject(dc, shell.smallFont);
    for (int i = 0; i < L.catCount; ++i) {
        auto cat = static_cast<Category>(i);
        HitId hid = static_cast<HitId>(static_cast<int>(HitId::Cat0) + i);
        RECT slot = L.traySlot[i];
        if (i > 0) slot.left += 1;
        bool empty = (cat == Category::Gases || cat == Category::Plasma);
        BtnState st = btnState(shell, hid, shell.category == cat);
        drawButton(dc, slot, L"", st, empty ? RGB(70, 70, 70) : kBtn, kBtnSel);
        int iconW = drawCategoryIcon(dc, slot, i, empty && shell.category != cat);
        RECT label = slot;
        UINT fmt = DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS;
        if (iconW > 0) {
            label.left += iconW;
            label.right -= 4;
            fmt |= DT_LEFT;
        } else {
            fmt |= DT_CENTER;
        }
        COLORREF textCol = empty ? kDimText : kText;
        drawLabel(dc, label, shell.categoryName(cat), textCol, fmt);
    }
    SelectObject(dc, shell.uiFont);

    fillRect(dc, L.console, kBlack);
    frameRect(dc, L.console, kBorder);
    SelectObject(dc, shell.consoleFont);
    std::wstring cons = L"> ";
    if (shell.consoleFocused) {
        cons += shell.consoleDraft;
        cons += L"_";
    } else if (!shell.consoleLines.empty()) {
        cons += shell.consoleLines.back();
    } else {
        cons += tr("console_hint");
    }
    drawLabel(dc, insetLike(L.console, 4), cons.c_str(), RGB(200, 200, 200));
    SelectObject(dc, shell.uiFont);

    fillRect(dc, L.sidebar, RGB(40, 40, 40));
    frameRect(dc, L.sidebar, kBorder);

    fillRect(dc, L.search, kBlack);
    frameRect(dc, L.search, kBorder);
    SelectObject(dc, shell.consoleFont);
    std::wstring search = L"> ";
    if (shell.searchText.empty() && !shell.searchFocused) search += tr("search_hint");
    else {
        search += shell.searchText;
        if (shell.searchFocused) search += L"_";
    }
    drawLabel(dc, insetLike(L.search, 4), search.c_str(), shell.searchFocused || !shell.searchText.empty() ? kText : kDimText);
    SelectObject(dc, shell.uiFont);

    fillRect(dc, L.catCol, RGB(58, 58, 58));
    frameRect(dc, L.catCol, kBorder);
    SelectObject(dc, shell.smallFont);
    if (shell.layout.palCount == 0) {
        drawLabel(dc, insetLike(L.catCol, 6), tr("none_yet"), kDimText, DT_CENTER | DT_TOP | DT_WORDBREAK);
    }
    for (int i = 0; i < shell.layout.palCount; ++i) {
        PaletteId pid = shell.elementAt(i);
        HitId hid = static_cast<HitId>(static_cast<int>(HitId::Pal0) + i);
        COLORREF fill = kBtn;
        if (pid == PaletteId::Water) fill = RGB(22, 90, 150);
        if (pid == PaletteId::Honey) fill = RGB(176, 110, 22);
        if (pid == PaletteId::Wood) fill = RGB(120, 88, 50);
        if (pid == PaletteId::Stone) fill = RGB(96, 96, 100);
        if (pid == PaletteId::Glass) fill = RGB(70, 120, 140);
        if (pid == PaletteId::Metal) fill = RGB(110, 112, 118);
        if (pid == PaletteId::Erase) fill = RGB(130, 55, 58);
        if (pid == PaletteId::Grab) fill = RGB(48, 96, 62);
        if (pid == PaletteId::Wall) fill = RGB(86, 91, 102);
        if (pid == PaletteId::Heat) fill = RGB(160, 70, 40);
        if (pid == PaletteId::Cool) fill = RGB(50, 90, 150);
        if (pid == PaletteId::Pressurize) fill = RGB(150, 110, 50);
        if (pid == PaletteId::Depressurize) fill = RGB(70, 72, 92);
        drawButton(dc, shell.layout.palSlot[i], shell.paletteName(pid),
            btnState(shell, hid, shell.palette == pid), fill, fill);
    }

    fillRect(dc, L.propsCol, kPanel);
    frameRect(dc, L.propsCol, kBorder);
    fillRect(dc, L.propsHeader, kBtn);
    frameRect(dc, L.propsHeader, kBorder);
    SelectObject(dc, shell.uiFont);
    drawLabel(dc, L.propsHeader, tr("properties"), kText, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    SelectObject(dc, shell.smallFont);
    RECT pb = insetLike(L.propsBody, 8);
    int y = pb.top;
    auto propLine = [&](wchar_t const *t) {
        RECT r{pb.left, y, pb.right, y + 18};
        drawLabel(dc, r, t, kText, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        y += 20;
    };
    propLine(tr("prop_intro"));
    y += 4;
    RECT bLab{pb.left, y, pb.right, y + 18};
    std::wstring br = trf("prop_brush", std::to_wstring(view.brushRadius));
    drawLabel(dc, bLab, br.c_str(), kText);
    drawButton(dc, L.propBrushMinus, L"-", btnState(shell, HitId::BrushMinus, false));
    drawButton(dc, L.propBrushPlus, L"+", btnState(shell, HitId::BrushPlus, false));
    y = L.propBrushMinus.bottom + 10;
    if (shell.category == Category::Solids) {
        RECT aLab{pb.left, y, pb.right, y + 18};
        drawLabel(dc, aLab, tr("prop_anchored"), kText);
        RECT box = L.propAnchored;
        fillRect(dc, box, kBlack);
        frameRect(dc, box, kBorder);
        if (shell.placeAnchored) {
            RECT inn = insetLike(box, 4);
            fillRect(dc, inn, kText);
        }
        RECT al{box.right + 8, box.top, pb.right, box.bottom};
        drawLabel(dc, al, shell.placeAnchored ? tr("prop_on_pinned") : tr("prop_off"), kDimText);
        y = box.bottom + 8;
        RECT sLab{pb.left, y, pb.right, y + 18};
        drawLabel(dc, sLab, tr("prop_sleeping"), kText);
        RECT sbox = L.propSleeping;
        fillRect(dc, sbox, kBlack);
        frameRect(dc, sbox, kBorder);
        if (shell.placeSleeping) {
            RECT inn = insetLike(sbox, 4);
            fillRect(dc, inn, kText);
        }
        RECT sl{sbox.right + 8, sbox.top, pb.right, sbox.bottom};
        drawLabel(dc, sl, shell.placeSleeping ? tr("prop_on_wakes") : tr("prop_off"), kDimText);
        y = sbox.bottom + 8;
        RECT pLab{pb.left, y, pb.right, y + 18};
        drawLabel(dc, pLab, tr("prop_powder"), kText);
        RECT pbox = L.propPowder;
        fillRect(dc, pbox, kBlack);
        frameRect(dc, pbox, kBorder);
        if (shell.placePowder) {
            RECT inn = insetLike(pbox, 4);
            fillRect(dc, inn, kText);
        }
        RECT pl{pbox.right + 8, pbox.top, pb.right, pbox.bottom};
        drawLabel(dc, pl, shell.placePowder ? tr("prop_on") : tr("prop_off"), kDimText);
        y = pbox.bottom + 8;
        RECT gLab{pb.left, y, pb.right, y + 18};
        drawLabel(dc, gLab, trf("prop_grain", std::to_wstring(shell.powderParticleSize)).c_str(),
            shell.placePowder ? kText : kDimText);
        drawButton(dc, L.propGrainMinus, L"-", btnState(shell, HitId::GrainMinus, false));
        drawButton(dc, L.propGrainPlus, L"+", btnState(shell, HitId::GrainPlus, false));
        y = L.propGrainMinus.bottom + 8;
        RECT note{pb.left, y, pb.right, y + 72};
        drawLabel(dc, note, tr("prop_solids_note"), kDimText, DT_LEFT | DT_TOP | DT_WORDBREAK);
        y += 76;
        if (shell.palette == PaletteId::Stone) propLine(tr("prop_material_stone"));
        else if (shell.palette == PaletteId::Glass) propLine(tr("prop_material_glass"));
        else if (shell.palette == PaletteId::Metal) propLine(tr("prop_material_metal"));
        else propLine(tr("prop_material_wood"));
    } else if (shell.category == Category::Fluids) {
        propLine(shell.palette == PaletteId::Honey ? tr("prop_substance_honey") : tr("prop_substance_water"));
        propLine(shell.dyeOnly ? tr("prop_dye_mode") : tr("prop_fill_paint"));
        RECT dLab{pb.left, y, pb.right, y + 16};
        drawLabel(dc, dLab, tr("prop_dye"), kText);
        for (int i = 0; i < 8; ++i) {
            RECT rc = L.propDyeSwatch[i];
            fillRect(dc, rc, kDyeSwatches[i].vis);
            frameRect(dc, rc, shell.dyeSwatch == i ? kText : kBorder);
            if (i == 0) {
                RECT inn = insetLike(rc, 4);
                frameRect(dc, inn, kDimText);
            }
            if (shell.dyeSwatch == i) {
                RECT inn = insetLike(rc, 2);
                frameRect(dc, inn, kText);
            }
        }
        RECT dyeLab{pb.left, L.propDyeOnly.top - 18, pb.right, L.propDyeOnly.top};
        drawLabel(dc, dyeLab, tr("prop_dye_only"), kText);
        RECT dbox = L.propDyeOnly;
        fillRect(dc, dbox, kBlack);
        frameRect(dc, dbox, kBorder);
        if (shell.dyeOnly) {
            RECT inn = insetLike(dbox, 4);
            fillRect(dc, inn, kText);
        }
        RECT dl{dbox.right + 8, dbox.top, pb.right, dbox.bottom};
        drawLabel(dc, dl, shell.dyeOnly ? tr("prop_on") : tr("prop_off"), kDimText);
        RECT strLab{pb.left, L.propDyeStrMinus.top - 18, pb.right, L.propDyeStrMinus.top};
        drawLabel(dc, strLab, trf("prop_dye_strength", shortFloat(shell.dyeStrength, 2)).c_str(), kText);
        drawButton(dc, L.propDyeStrMinus, L"-", btnState(shell, HitId::DyeStrMinus, false));
        drawButton(dc, L.propDyeStrPlus, L"+", btnState(shell, HitId::DyeStrPlus, false));
        y = L.propDyeStrMinus.bottom + 8;
        RECT note{pb.left, y, pb.right, y + 56};
        drawLabel(dc, note, tr("prop_dye_note"), kDimText, DT_LEFT | DT_TOP | DT_WORDBREAK);
    } else if (shell.category == Category::Tools) {
        if (shell.palette == PaletteId::Grab) propLine(tr("prop_tool_grab"));
        else propLine(tr("prop_tool_eraser"));
    } else if (shell.category == Category::Energy) {
        if (shell.palette == PaletteId::Cool) propLine(tr("prop_tool_cool"));
        else if (shell.palette == PaletteId::Pressurize) propLine(tr("prop_tool_pressurize"));
        else if (shell.palette == PaletteId::Depressurize) propLine(tr("prop_tool_depressurize"));
        else propLine(tr("prop_tool_heat"));
        std::wstring pw = trf("prop_power", shortFloat(view.heatPower, 2));
        RECT pLab{pb.left, y, pb.right, y + 18};
        drawLabel(dc, pLab, pw.c_str(), kText);
        drawButton(dc, L.propPowerMinus, L"-", btnState(shell, HitId::PowerMinus, false));
        drawButton(dc, L.propPowerPlus, L"+", btnState(shell, HitId::PowerPlus, false));
        y = L.propPowerPlus.bottom + 8;
    } else if (shell.category == Category::Misc) {
        propLine(tr("prop_tool_wall"));
    } else {
        propLine(tr("prop_empty_tab"));
    }

    fillRect(dc, L.inspectCol, kPanel);
    frameRect(dc, L.inspectCol, kBorder);
    fillRect(dc, L.inspectHeader, kBtn);
    frameRect(dc, L.inspectHeader, kBorder);
    SelectObject(dc, shell.uiFont);
    drawLabel(dc, L.inspectHeader, tr("inspector"), kText, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
    SelectObject(dc, shell.smallFont);
    RECT ib = insetLike(L.inspectBody, 8);
    y = ib.top;
    auto ins = [&](std::wstring const &t) {
        RECT r{ib.left, y, ib.right, y + 16};
        drawLabel(dc, r, t.c_str(), kText, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_END_ELLIPSIS);
        y += 17;
    };
    ins(shell.paletteName(shell.palette));
    {
        RECT desc{ib.left, y, ib.right, y + 56};
        drawLabel(dc, desc, shell.paletteHint(shell.palette), kDimText, DT_LEFT | DT_TOP | DT_WORDBREAK);
        y += 60;
    }
    ins(L"");
    if (ViewBarDef const *hoverView = viewBarByHit(shell.hoverId)) {
        ins(tr(hoverView->titleKey));
        RECT help{ib.left, y, ib.right, ib.bottom};
        drawLabel(dc, help, tr(hoverView->helpKey), kDimText, DT_LEFT | DT_TOP | DT_WORDBREAK);
        SelectObject(dc, old);
        return;
    }
    FluidEngine const *e = view.engine;
    RigidBodyEngine const *rg = view.rigid;
    if (!e || !FluidEngine::inside(view.hoverX, view.hoverY)) {
        ins(tr("ins_nothing"));
        ins(tr("ins_hover_1"));
        ins(tr("ins_hover_2"));
    } else {
        int hx = view.hoverX, hy = view.hoverY;
        int hi = FluidEngine::ci(hx, hy);
        ins(trf("info_cell", std::to_wstring(hx), std::to_wstring(hy)));
        if (rg && view.gas) {
            ThermalCellSample s = ThermalEngine::sampleCell(*e, *rg, *view.gas, hx, hy);
            if (view.debugView == DebugView::Temperature) {
                ins(tr("ins_temp_header"));
                ins(trf("ins_temp", shortFloat(s.temperatureK, 2)));
                wchar_t dTbuf[32]{};
                swprintf_s(dTbuf, L"%+.2f", static_cast<double>(s.temperatureK - AMBIENT_TEMPERATURE_K));
                ins(trf("ins_dT", dTbuf));
                wchar_t const *matName = tr("ins_mat_none");
                if (s.kind == ThermalSampleKind::Liquid)
                    matName = (e->honeyFraction(hi) > 0.5f) ? tr("ins_mat_honey") : tr("ins_mat_water");
                else if (s.kind == ThermalSampleKind::Gas) matName = tr("ins_mat_air");
                else if (s.kind == ThermalSampleKind::Wall) matName = tr("stone");
                else if (s.kind == ThermalSampleKind::Rigid) {
                    if (s.materialId == MATERIAL_STONE) matName = tr("stone");
                    else if (s.materialId == MATERIAL_GLASS) matName = tr("glass");
                    else if (s.materialId == MATERIAL_METAL) matName = tr("metal");
                    else if (s.materialId == MATERIAL_WOOD) matName = tr("wood");
                    else matName = tr("ins_mat_none");
                }
                ins(trf("ins_material", matName));
                if (s.capacityJK > MIN_THERMAL_CAPACITY)
                    ins(trf("ins_thermal_energy", shortFloat(s.energyJ, 0)));
                bool active = view.thermal && view.thermal->isChunkActive(hx, hy);
                ins(trf("ins_thermal_active", active ? tr("yes") : tr("no")));
            } else {
                ins(trf("ins_temp", shortFloat(s.temperatureK, 1)));
            }
        }
        if (e->solid[static_cast<size_t>(hi)]) ins(tr("ins_kind_wall"));
        int body = rg ? rg->occupant[static_cast<size_t>(hi)] : -1;
        if (body >= 0 && rg && body < static_cast<int>(rg->bodies.size())) {
            RigidBody const &b = rg->bodies[static_cast<size_t>(body)];
            MaterialId mat = rg->worldCellMaterial(hx, hy);
            wchar_t const *matName = tr("wood");
            if (mat == MATERIAL_STONE) matName = tr("stone");
            else if (mat == MATERIAL_GLASS) matName = tr("glass");
            else if (mat == MATERIAL_METAL) matName = tr("metal");
            ins(tr("ins_kind_rigid"));
            ins(trf("ins_material", matName));
            ins(trf("ins_id", std::to_wstring(b.id)));
            ins(trf("ins_component", std::to_wstring(b.id)));
            ins(trf("ins_mass", shortFloat(b.mass, 2)));
            ins(trf("ins_pos", shortFloat(b.x, 1), shortFloat(b.y, 1)));
            ins(trf("ins_vel", shortFloat(b.vx, 2), shortFloat(b.vy, 2)));
            ins(trf("ins_omega", shortFloat(b.omega, 2)));
            ins(trf("ins_supported", b.supported ? tr("yes") : tr("no")));
            ins(trf("ins_anchored", b.anchored ? tr("yes") : tr("no")));
            ins(trf("ins_sleeping", b.dormant ? tr("yes") : tr("no")));
            ins(trf("ins_at_rest", (b.sleeping && !b.anchored && !b.dormant) ? tr("yes") : tr("no")));
            ins(trf("ins_pen", shortFloat(b.maxPenetration, 3)));
            ins(trf("ins_jnjt", shortFloat(b.debugJn, 2), shortFloat(b.debugJt, 2)));
            ins(trf("ins_poscorr", shortFloat(b.debugPosCorrX, 3), shortFloat(b.debugPosCorrY, 3)));
            ins(trf("ins_damage", shortFloat(b.maxDamage, 2)));
            ins(trf("ins_bond_body", shortFloat(b.maxBondDamage, 2), std::to_wstring(b.brokenBondCount)));
            ins(trf("ins_moisture", shortFloat(b.absorbedLiquid, 2), shortFloat(b.cachedWetness, 2)));
            float matD = 0.0f, bondD = 0.0f, strength = 0.0f, crack = 0.0f, wet = 0.0f;
            int brokenN = 0;
            if (rg->inspectLocalStructure(hx, hy, matD, bondD, brokenN, strength, crack, wet)) {
                ins(trf("ins_mat_damage", shortFloat(matD, 2)));
                ins(trf("ins_bond_damage", shortFloat(bondD, 2)));
                ins(trf("ins_broken", std::to_wstring(brokenN)));
                ins(trf("ins_strength", shortFloat(strength, 2)));
                ins(trf("ins_crack", shortFloat(crack, 2)));
                ins(trf("ins_saturation", shortFloat(wet, 2)));
                if (brokenN > 0) ins(tr("ins_bond_broken"));
            }
        } else if (!e->solid[static_cast<size_t>(hi)]) {
            float fill = e->fill[static_cast<size_t>(hi)];
            ins(tr("ins_kind_liquid"));
            ins(trf("ins_fill", shortFloat(fill, 3)));
            float dr = e->dyeR[static_cast<size_t>(hi)];
            float dg = e->dyeG[static_cast<size_t>(hi)];
            float db = e->dyeB[static_cast<size_t>(hi)];
            float dyeAmt = (fill > 1.0e-8f) ? std::max(dr, std::max(dg, db)) / fill : 0.0f;
            if (dyeAmt > 0.02f) ins(trf("ins_dye", shortFloat(dyeAmt * 100.0, 0)));
            ins(trf("ins_pressure", shortFloat(e->pressure[static_cast<size_t>(hi)], 2)));
            ins(trf("ins_vel", shortFloat(e->cellU(hx, hy), 2), shortFloat(e->cellV(hx, hy), 2)));
            ins(trf("ins_surface", e->surfaceMask[static_cast<size_t>(hi)] ? tr("yes") : tr("no")));
            if (fill > 1.0e-8f) {
                ins(L"");
                int honeyPct = std::clamp(static_cast<int>(std::lround(static_cast<double>(e->honeyFraction(hi) * 100.0f))), 0, 100);
                int waterPct = 100 - honeyPct;
                ins(trf("ins_comp_water", std::to_wstring(waterPct)));
                ins(trf("ins_comp_honey", std::to_wstring(honeyPct)));
            }
        }
        ins(trf("ins_volume", shortFloat(e->currentVolume, 1), shortFloat(e->expectedVolume, 1)));
        GasEngine const *g = view.gas;
        if (g && !e->solid[static_cast<size_t>(hi)] && body < 0) {
            int gi = hi;
            ins(tr("ins_gas_material"));
            ins(trf("ins_gas_amount", shortFloat(g->amount[static_cast<size_t>(gi)], 3)));
            ins(trf("ins_gas_vol", shortFloat(g->volume[static_cast<size_t>(gi)], 3)));
            ins(trf("ins_gas_pa", shortFloat(g->pressurePa(gi), 0)));
            ins(trf("ins_gas_atm", shortFloat(g->pressureAtm(gi), 3)));
            ins(trf("ins_gas_vel", shortFloat(g->cellU(hx, hy), 2), shortFloat(g->cellV(hx, hy), 2)));
            ins(tr("ins_gas_comp"));
            ins(trf("ins_gas_world", shortFloat(g->currentAmount, 1), shortFloat(g->expectedAmount, 1)));
        }
    }

    SelectObject(dc, old);
}

RECT settingsPanelRect(ShellState const &shell) {
    return RECT{shell.settingsX, shell.settingsY, shell.settingsX + shell.settingsW, shell.settingsY + shell.settingsH};
}

void clampSettingsWindow(ShellState &shell) {
    int clientW = std::max(1, static_cast<int>(shell.layout.client.right - shell.layout.client.left));
    int clientH = std::max(1, static_cast<int>(shell.layout.client.bottom - shell.layout.client.top));
    int minVisible = 36;
    shell.settingsH = std::clamp(shell.settingsH, 220, std::max(220, clientH - 16));
    shell.settingsX = std::clamp(shell.settingsX, 8 - shell.settingsW + 80, clientW - minVisible);
    shell.settingsY = std::clamp(shell.settingsY, 8, std::max(8, clientH - minVisible));
    if (shell.settingsY + shell.settingsH > clientH - 8)
        shell.settingsY = std::max(8, clientH - 8 - shell.settingsH);
}

void clampSettingsScroll(ShellState &shell) {
    int viewH = std::max(1, static_cast<int>(shell.settingsBody.bottom - shell.settingsBody.top));
    int maxScroll = std::max(0, shell.settingsContentH - viewH);
    shell.settingsScrollY = std::clamp(shell.settingsScrollY, 0, maxScroll);
}

void dragSettingsWindow(ShellState &shell, int x, int y) {
    shell.settingsX = x - shell.settingsDragOX;
    shell.settingsY = y - shell.settingsDragOY;
    clampSettingsWindow(shell);
}

void dragSettingsScroll(ShellState &shell, int y) {
    int trackH = shell.settingsScrollTrack.bottom - shell.settingsScrollTrack.top;
    int thumbH = shell.settingsThumb.bottom - shell.settingsThumb.top;
    int viewH = std::max(1, static_cast<int>(shell.settingsBody.bottom - shell.settingsBody.top));
    int maxScroll = std::max(0, shell.settingsContentH - viewH);
    int travel = std::max(1, trackH - thumbH);
    int top = y - shell.settingsScrollDragOY - shell.settingsScrollTrack.top;
    float t = std::clamp(static_cast<float>(top) / static_cast<float>(travel), 0.0f, 1.0f);
    shell.settingsScrollY = static_cast<int>(t * static_cast<float>(maxScroll) + 0.5f);
    clampSettingsScroll(shell);
}

bool handleSettingsWheel(ShellState &shell, int x, int y, int wheelDelta) {
    if (!ptIn(settingsPanelRect(shell), x, y)) return false;
    int viewH = std::max(1, static_cast<int>(shell.settingsBody.bottom - shell.settingsBody.top));
    int step = std::max(24, viewH / 8);
    int notches = wheelDelta / 120;
    if (notches == 0) notches = wheelDelta > 0 ? 1 : -1;
    shell.settingsScrollY -= notches * step;
    shell.settingsFlyout = SettingsFlyout::None;
    clampSettingsScroll(shell);
    return true;
}

void ensureSettingsPlacement(ShellState &shell) {
    int clientW = std::max(1, static_cast<int>(shell.layout.client.right - shell.layout.client.left));
    int clientH = std::max(1, static_cast<int>(shell.layout.client.bottom - shell.layout.client.top));
    if (shell.settingsX < 0 || shell.settingsY < 0) {
        shell.settingsX = std::max(8, (clientW - shell.settingsW) / 2);
        shell.settingsY = std::max(8, (clientH - shell.settingsH) / 4);
    }
    clampSettingsWindow(shell);
}

bool settingsCoversPoint(ShellState const &shell, int x, int y) {
    RECT panel = settingsPanelRect(shell);
    if (ptIn(panel, x, y)) return true;
    if (shell.settingsFlyout != SettingsFlyout::None && ptIn(shell.settingsFlyoutRc, x, y)) return true;
    return false;
}

namespace {

struct FlyItem {
    wchar_t const *label = nullptr;
    MenuCmd cmd = MenuCmd::None;
    bool on = false;
};

RECT flyoutRectFor(RECT const &anchor, int cols, int rows, int colW, int rowH, RECT const &client) {
    int pad = 4;
    int fw = cols * colW + pad * 2;
    int fh = rows * rowH + pad * 2;
    int x = anchor.right + 4;
    int y = anchor.top;
    if (x + fw > client.right - 8) x = anchor.left - fw - 4;
    if (x < client.left + 8) x = client.left + 8;
    if (y + fh > client.bottom - 8) y = std::max(client.top + 8, client.bottom - 8 - fh);
    if (y < client.top + 8) y = client.top + 8;
    return RECT{x, y, x + fw, y + fh};
}

void drawFlyoutList(HDC dc, ShellState &shell, RECT const &box, FlyItem const *items, int count, int cols, int colW, int rowH) {
    fillRect(dc, box, RGB(36, 36, 36));
    frameRect(dc, box, kBorder);
    int pad = 4;
    int rows = (count + cols - 1) / std::max(1, cols);
    for (int i = 0; i < count; ++i) {
        int col = i / rows;
        int row = i % rows;
        RECT rc{box.left + pad + col * colW, box.top + pad + row * rowH,
            box.left + pad + (col + 1) * colW - 2, box.top + pad + (row + 1) * rowH - 2};
        bool over = ptIn(rc, shell.mouseX, shell.mouseY);
        drawButton(dc, rc, items[i].label,
            items[i].on ? BtnState::Selected : (over ? BtnState::Hover : BtnState::Normal),
            RGB(70, 70, 70), RGB(110, 110, 110));
        shell.menuHits.push_back({rc, items[i].cmd, nullptr});
    }
}

} // namespace

SettingsMouseResult handleSettingsMouseDown(ShellState &shell, View const &view, int x, int y, MenuCmd &outCmd) {
    outCmd = MenuCmd::None;
    if (!view.settingsOpen) return SettingsMouseResult::Miss;
    if (ptIn(shell.settingsThumb, x, y) && shell.settingsContentH > (shell.settingsBody.bottom - shell.settingsBody.top)) {
        shell.settingsScrollDragging = true;
        shell.settingsScrollDragOY = y - shell.settingsThumb.top;
        return SettingsMouseResult::Drag;
    }
    if (ptIn(shell.settingsScrollTrack, x, y)) {
        int viewH = std::max(1, static_cast<int>(shell.settingsBody.bottom - shell.settingsBody.top));
        int page = std::max(24, viewH - 24);
        if (y < shell.settingsThumb.top) shell.settingsScrollY -= page;
        else shell.settingsScrollY += page;
        clampSettingsScroll(shell);
        return SettingsMouseResult::Consume;
    }
    for (int i = static_cast<int>(shell.menuHits.size()) - 1; i >= 0; --i) {
        MenuHit const &hit = shell.menuHits[static_cast<size_t>(i)];
        if (!ptIn(hit.rc, x, y)) continue;
        if (hit.cmd == MenuCmd::None) return SettingsMouseResult::Consume;
        outCmd = hit.cmd;
        return SettingsMouseResult::Command;
    }
    if (ptIn(shell.settingsTitleBar, x, y)) {
        shell.settingsDragging = true;
        shell.settingsDragOX = x - shell.settingsX;
        shell.settingsDragOY = y - shell.settingsY;
        return SettingsMouseResult::Drag;
    }
    if (settingsCoversPoint(shell, x, y)) return SettingsMouseResult::Consume;
    return SettingsMouseResult::Miss;
}

void drawSettingsPanel(HDC dc, ShellState &shell, View const &view, RECT const &client) {
    shell.menuHits.clear();
    shell.settingsFlyoutRc = RECT{};
    if (!view.settingsOpen || !view.engine || !view.rigid) {
        shell.settingsFlyout = SettingsFlyout::None;
        shell.settingsDragging = false;
        return;
    }
    ensureSettingsPlacement(shell);
    FluidEngine const &engine = *view.engine;
    RigidBodyEngine const &rigid = *view.rigid;
    RECT panel = settingsPanelRect(shell);

    fillRect(dc, panel, RGB(42, 42, 42));
    frameRect(dc, panel, kBorder);
    RECT inner = insetLike(panel, 1);
    frameRect(dc, inner, kBorderDim);

    HFONT old = static_cast<HFONT>(SelectObject(dc, shell.smallFont ? shell.smallFont : GetStockObject(DEFAULT_GUI_FONT)));
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, kText);

    RECT closeRc{panel.right - 72, panel.top + 6, panel.right - 8, panel.top + 28};
    shell.settingsTitleBar = RECT{panel.left + 1, panel.top + 1, closeRc.left - 4, panel.top + 30};
    bool titleOver = ptIn(shell.settingsTitleBar, shell.mouseX, shell.mouseY);
    if (titleOver || shell.settingsDragging) fillRect(dc, shell.settingsTitleBar, RGB(52, 52, 56));
    RECT title{panel.left + 12, panel.top + 6, closeRc.left - 6, panel.top + 28};
    drawLabel(dc, title, tr("settings_title"), kText);
    bool closeOver = ptIn(closeRc, shell.mouseX, shell.mouseY);
    drawButton(dc, closeRc, tr("close"), closeOver ? BtnState::Hover : BtnState::Normal);
    shell.menuHits.push_back({closeRc, MenuCmd::Close, "tip_close"});

    shell.settingsBody = RECT{panel.left + 1, panel.top + 31, panel.right - 1, panel.bottom - 1};
    int const scrollW = 10;
    RECT body = shell.settingsBody;
    body.right -= scrollW + 2;
    clampSettingsScroll(shell);

    int saved = SaveDC(dc);
    IntersectClipRect(dc, body.left, body.top, body.right, body.bottom);

    int mid = (body.left + body.right) / 2;
    int left0 = body.left + 10, left1 = mid - 8;
    int right0 = mid + 8, right1 = body.right - 10;
    int rowH = 24;
    int yL = body.top + 6 - shell.settingsScrollY;
    int yR = body.top + 6 - shell.settingsScrollY;
    int x = left0;
    int const originY = body.top + 6 - shell.settingsScrollY;

    auto visible = [&](RECT const &rc) -> bool {
        RECT hit{};
        return IntersectRect(&hit, &rc, &body) != 0;
    };
    auto section = [&](int x0, int x1, int &y, wchar_t const *name) {
        RECT lab{x0, y, x1, y + 16};
        drawLabel(dc, lab, name, kDimText, DT_LEFT | DT_BOTTOM | DT_SINGLELINE);
        y += 18;
    };
    auto button = [&](int x0, int x1, int &x, int &y, int w, wchar_t const *label, MenuCmd cmd, bool on, COLORREF color = RGB(80, 80, 80)) -> RECT {
        if (x + w > x1) { y += rowH + 4; x = x0; }
        RECT rc{x, y, std::min(x + w, x1), y + rowH};
        bool over = visible(rc) && ptIn(rc, shell.mouseX, shell.mouseY);
        drawButton(dc, rc, label, on ? BtnState::Selected : (over ? BtnState::Hover : BtnState::Normal), color, RGB(110, 110, 110));
        if (visible(rc)) shell.menuHits.push_back({rc, cmd, nullptr});
        x += w + 6;
        return rc;
    };

    section(left0, left1, yL, tr("sec_quality"));
    x = left0;
    button(left0, left1, x, yL, 52, tr("quality_low"), MenuCmd::QualityLow, engine.config.quality == QualityPreset::Low);
    button(left0, left1, x, yL, 68, tr("quality_med"), MenuCmd::QualityMed, engine.config.quality == QualityPreset::Medium);
    button(left0, left1, x, yL, 52, tr("quality_high"), MenuCmd::QualityHigh, engine.config.quality == QualityPreset::High);
    button(left0, left1, x, yL, 52, tr("quality_auto"), MenuCmd::QualityAuto, engine.config.quality == QualityPreset::Auto);
    yL += rowH + 4;
    x = left0;
    button(left0, left1, x, yL, 84, tr("thermal_on"), MenuCmd::ThermalOn, view.thermalEnabled);
    button(left0, left1, x, yL, 88, tr("thermal_off"), MenuCmd::ThermalOff, !view.thermalEnabled);
    yL += rowH + 4;
    x = left0;
    {
        int interval = view.thermal ? view.thermal->config.intervalTicks : 2;
        std::wstring heatEvery = trf("heat_every", std::to_wstring(interval));
        button(left0, left1, x, yL, 24, L"-", MenuCmd::ThermalIntervalMinus, false);
        button(left0, left1, x, yL, 110, heatEvery.c_str(), MenuCmd::None, false);
        button(left0, left1, x, yL, 24, L"+", MenuCmd::ThermalIntervalPlus, false);
    }
    yL += rowH + 4;
    {
        int shown = engine.config.quality == QualityPreset::Auto ? engine.config.autoQualityLevel
            : (engine.config.quality == QualityPreset::High ? 2
                : (engine.config.quality == QualityPreset::Low ? 0 : 1));
        std::wstring usingName = shown <= 0 ? tr("quality_low") : (shown == 1 ? tr("quality_med") : tr("quality_high"));
        std::wstring qStatus = engine.config.quality == QualityPreset::Auto
            ? trf("quality_auto_using", usingName)
            : trf("quality_using", usingName);
        RECT qRc{left0, yL, left1, yL + 16};
        drawLabel(dc, qRc, qStatus.c_str(), kDimText);
        yL += 18; x = left0;
    }

    section(left0, left1, yL, tr("sec_world"));
    x = left0;
    button(left0, left1, x, yL, 68, view.paused ? tr("resume_word") : tr("pause_word"), MenuCmd::Pause, view.paused, kPauseOn);
    button(left0, left1, x, yL, 56, tr("step"), MenuCmd::Step, false);
    button(left0, left1, x, yL, 56, tr("clear_word"), MenuCmd::Clear, false, kClear);
    button(left0, left1, x, yL, 56, tr("reset"), MenuCmd::Reset, false);
    button(left0, left1, x, yL, 56, tr("slosh"), MenuCmd::Slosh, false);
    button(left0, left1, x, yL, 108, tr("btn_walled"), MenuCmd::WalledBorders, engine.config.walledBorders);
    button(left0, left1, x, yL, 72, tr("btn_gravity"), MenuCmd::RigidGravity, rigid.gravityScale > 0.5f);
    yL += rowH + 8; x = left0;

    section(left0, left1, yL, tr("sec_look"));
    x = left0;
    {
        RECT sLab{left0, yL, left1, yL + 16};
        drawLabel(dc, sLab, tr("lbl_style"), kDimText, DT_LEFT | DT_BOTTOM | DT_SINGLELINE);
        yL += 18; x = left0;
    }
    button(left0, left1, x, yL, 48, tr("style_flat"), MenuCmd::StyleFlat, view.worldLook.style == WorldRenderStyle::FlatColor);
    button(left0, left1, x, yL, 56, tr("style_noisy"), MenuCmd::StyleNoisy, view.worldLook.style == WorldRenderStyle::NoisyFlat);
    button(left0, left1, x, yL, 60, tr("style_detailed"), MenuCmd::StyleDetailed, view.worldLook.style == WorldRenderStyle::Detailed);
    button(left0, left1, x, yL, 48, tr("style_realistic"), MenuCmd::StyleRealistic, view.worldLook.style == WorldRenderStyle::Realistic);
    button(left0, left1, x, yL, 52, tr("style_alpha"), MenuCmd::StyleAlpha, view.worldLook.style == WorldRenderStyle::AlphaFlat);
    yL += rowH + 4; x = left0;
    {
        RECT eLab{left0, yL, left1, yL + 16};
        drawLabel(dc, eLab, tr("lbl_effects"), kDimText, DT_LEFT | DT_BOTTOM | DT_SINGLELINE);
        yL += 18; x = left0;
    }
    button(left0, left1, x, yL, 120, tr("btn_glow_liquids"), MenuCmd::Glow, view.worldLook.glowingLiquids);
    button(left0, left1, x, yL, 80, tr("btn_outlines"), MenuCmd::Outlines, view.worldLook.outlines);
    button(left0, left1, x, yL, 96, tr("btn_rigid_overlay"), MenuCmd::Overlay, rigid.debugOverlay);
    yL += rowH + 8; x = left0;

    section(left0, left1, yL, tr("sec_solver"));
    x = left0;
    button(left0, left1, x, yL, 80, tr("btn_vorticity"), MenuCmd::Vorticity, engine.config.vorticityEnabled);
    button(left0, left1, x, yL, 80, tr("btn_res_merge"), MenuCmd::Residual, engine.residualConsolidationEnabled);
    button(left0, left1, x, yL, 72, tr("btn_tension"), MenuCmd::Tension, engine.config.surfaceTensionEnabled);
    button(left0, left1, x, yL, 56, tr("btn_spray"), MenuCmd::Spray, engine.config.sprayEnabled);
    yL += rowH + 4;
    x = left0;
    RECT advTrigger{};
    {
        char const *advKey = "adv_menu_sl";
        switch (engine.config.velocityAdvection) {
            case VelocityAdvection::None: advKey = "adv_menu_none"; break;
            case VelocityAdvection::FirstOrderUpwind: advKey = "adv_menu_fou"; break;
            case VelocityAdvection::NearestSemiLagrangian: advKey = "adv_menu_nsl"; break;
            case VelocityAdvection::SemiLagrangian: advKey = "adv_menu_sl"; break;
            case VelocityAdvection::MacCormack: advKey = "adv_menu_macc"; break;
            case VelocityAdvection::BFECC: advKey = "adv_menu_bfecc"; break;
        }
        std::wstring advLabel = trf("adv_trigger", tr(advKey));
        advTrigger = button(left0, left1, x, yL, std::min(left1 - left0, 220), advLabel.c_str(), MenuCmd::None,
            shell.settingsFlyout == SettingsFlyout::Advection);
    }
    yL += rowH + 4;
    x = left0;
    std::wstring vort = trf("vort_str", shortFloat(engine.config.vorticityStrength, 3));
    button(left0, left1, x, yL, 24, L"-", MenuCmd::VorticityMinus, false);
    button(left0, left1, x, yL, 96, vort.c_str(), MenuCmd::None, false);
    button(left0, left1, x, yL, 24, L"+", MenuCmd::VorticityPlus, false);
    yL += rowH + 4;
    x = left0;
    std::wstring pressure = trf("p_iter", std::to_wstring(engine.config.maxPressureIterations));
    std::wstring brush = trf("brush", std::to_wstring(view.brushRadius));
    button(left0, left1, x, yL, 24, L"-", MenuCmd::PressureMinus, false);
    button(left0, left1, x, yL, 80, pressure.c_str(), MenuCmd::None, false);
    button(left0, left1, x, yL, 24, L"+", MenuCmd::PressurePlus, false);
    button(left0, left1, x, yL, 24, L"-", MenuCmd::BrushMinus, false);
    button(left0, left1, x, yL, 72, brush.c_str(), MenuCmd::None, false);
    button(left0, left1, x, yL, 24, L"+", MenuCmd::BrushPlus, false);
    yL += rowH + 4;
    x = left0;
    std::wstring substeps = trf("substeps", std::to_wstring(engine.config.maxSubsteps));
    std::wstring limiter = trf("limiter", std::to_wstring(engine.config.maxLimiterPasses));
    button(left0, left1, x, yL, 24, L"-", MenuCmd::SubstepsMinus, false);
    button(left0, left1, x, yL, 88, substeps.c_str(), MenuCmd::None, false);
    button(left0, left1, x, yL, 24, L"+", MenuCmd::SubstepsPlus, false);
    button(left0, left1, x, yL, 24, L"-", MenuCmd::LimiterMinus, false);
    button(left0, left1, x, yL, 88, limiter.c_str(), MenuCmd::None, false);
    button(left0, left1, x, yL, 24, L"+", MenuCmd::LimiterPlus, false);
    yL += rowH + 8;
    x = left0;

    section(right0, right1, yR, tr("sec_gas_sim"));
    x = right0;
    button(right0, right1, x, yR, 52, tr("gas_off"), MenuCmd::GasOff, view.gas && view.gas->config.simMode == GasSimMode::Off);
    button(right0, right1, x, yR, 52, tr("gas_half"), MenuCmd::GasHalf, view.gas && view.gas->config.simMode == GasSimMode::Half);
    button(right0, right1, x, yR, 52, tr("gas_full"), MenuCmd::GasFull, view.gas && view.gas->config.simMode == GasSimMode::Full);
    yR += rowH + 8;
    x = right0;

    section(right0, right1, yR, tr("sec_sim_rate"));
    x = right0;
    button(right0, right1, x, yR, 60, tr("hz_20"), MenuCmd::Hz20, engine.config.physicsHz == 20);
    button(right0, right1, x, yR, 60, tr("hz_30"), MenuCmd::Hz30, engine.config.physicsHz == 30);
    yR += rowH + 4;
    x = right0;
    std::wstring catchUp = trf("catch_up", std::to_wstring(engine.config.catchUpTicks));
    button(right0, right1, x, yR, 24, L"-", MenuCmd::CatchUpMinus, false);
    button(right0, right1, x, yR, 100, catchUp.c_str(), MenuCmd::None, false);
    button(right0, right1, x, yR, 24, L"+", MenuCmd::CatchUpPlus, false);
    yR += rowH + 8;
    x = right0;

    section(right0, right1, yR, tr("sec_threads"));
    x = right0;
    int maxW = engine.maxSelectableWorkers();
    button(right0, right1, x, yR, 52, tr("threads_auto"), MenuCmd::ThreadsAuto, engine.config.workerCount <= 0);
    button(right0, right1, x, yR, 32, L"1", MenuCmd::Threads1, engine.config.workerCount == 1);
    if (maxW >= 2) button(right0, right1, x, yR, 32, L"2", MenuCmd::Threads2, engine.config.workerCount == 2);
    if (maxW >= 4) button(right0, right1, x, yR, 32, L"4", MenuCmd::Threads4, engine.config.workerCount == 4);
    if (maxW >= 6) button(right0, right1, x, yR, 32, L"6", MenuCmd::Threads6, engine.config.workerCount == 6);
    if (maxW >= 8) button(right0, right1, x, yR, 32, L"8", MenuCmd::Threads8, engine.config.workerCount == 8);
    yR += rowH + 4;
    std::wstring thrStatus = engine.config.workerCount <= 0
        ? trf("threads_auto_using", std::to_wstring(engine.lastResolvedWorkers))
        : (trf("threads_using", std::to_wstring(engine.lastResolvedWorkers))
            + (engine.lastPressureParallel ? tr("threads_parallel") : tr("threads_serial")));
    RECT thrRc{right0, yR, right1, yR + 16};
    drawLabel(dc, thrRc, thrStatus.c_str(), kDimText);
    yR += 18;
    x = right0;

    section(right0, right1, yR, tr("sec_speed"));
    x = right0;
    wchar_t const *speedNames[6] = {L"0.1x", L"0.25x", L"0.5x", L"1x", L"2x", L"4x"};
    MenuCmd speedCmds[6] = {MenuCmd::Speed0, MenuCmd::Speed1, MenuCmd::Speed2, MenuCmd::Speed3, MenuCmd::Speed4, MenuCmd::Speed5};
    for (int i = 0; i < 6; ++i)
        button(right0, right1, x, yR, 44, speedNames[i], speedCmds[i], view.speedIndex == static_cast<size_t>(i));
    yR += rowH + 8;
    x = right0;

    section(right0, right1, yR, tr("sec_scenes"));
    x = right0;
    RECT fluidTrigger = button(right0, right1, x, yR, 108, tr("menu_fluid_scenes"), MenuCmd::None, shell.settingsFlyout == SettingsFlyout::FluidScenes);
    RECT rigidTrigger = button(right0, right1, x, yR, 108, tr("menu_rigid_scenes"), MenuCmd::None, shell.settingsFlyout == SettingsFlyout::RigidScenes);
    RECT gasTrigger = button(right0, right1, x, yR, 100, tr("menu_gas_scenes"), MenuCmd::None, shell.settingsFlyout == SettingsFlyout::GasScenes);

    shell.settingsContentH = std::max(yL, yR) + rowH + 8 - originY;
    RestoreDC(dc, saved);
    clampSettingsScroll(shell);

    int viewH = std::max(1, static_cast<int>(shell.settingsBody.bottom - shell.settingsBody.top));
    bool needScroll = shell.settingsContentH > viewH;
    shell.settingsScrollTrack = RECT{};
    shell.settingsThumb = RECT{};
    if (needScroll) {
        RECT track{panel.right - 12, shell.settingsBody.top + 2, panel.right - 4, shell.settingsBody.bottom - 2};
        shell.settingsScrollTrack = track;
        fillRect(dc, track, RGB(32, 32, 32));
        int maxScroll = std::max(1, shell.settingsContentH - viewH);
        int trackH = static_cast<int>(track.bottom - track.top);
        int thumbH = std::max(18, viewH * viewH / std::max(1, shell.settingsContentH));
        thumbH = std::min(thumbH, trackH);
        int travel = std::max(0, trackH - thumbH);
        int thumbY = track.top + (travel * shell.settingsScrollY) / maxScroll;
        RECT thumb{track.left + 1, thumbY, track.right - 1, thumbY + thumbH};
        shell.settingsThumb = thumb;
        bool thumbOver = ptIn(thumb, shell.mouseX, shell.mouseY) || shell.settingsScrollDragging;
        fillRect(dc, thumb, thumbOver ? RGB(150, 150, 150) : RGB(118, 118, 118));
    }

    RECT advBox = flyoutRectFor(advTrigger, 1, 6, 168, 22, client);
    RECT fluidBox = flyoutRectFor(fluidTrigger, 1, 10, 168, 22, client);
    RECT rigidBox = flyoutRectFor(rigidTrigger, 2, 12, 148, 22, client);
    RECT gasBox = flyoutRectFor(gasTrigger, 1, 8, 168, 22, client);

    SettingsFlyout want = SettingsFlyout::None;
    if (!shell.settingsDragging && !shell.settingsScrollDragging) {
        bool advVis = visible(advTrigger);
        bool fluidVis = visible(fluidTrigger);
        bool rigidVis = visible(rigidTrigger);
        bool gasVis = visible(gasTrigger);
        if ((advVis && ptIn(advTrigger, shell.mouseX, shell.mouseY)) || (shell.settingsFlyout == SettingsFlyout::Advection && ptIn(advBox, shell.mouseX, shell.mouseY)))
            want = SettingsFlyout::Advection;
        else if ((fluidVis && ptIn(fluidTrigger, shell.mouseX, shell.mouseY)) || (shell.settingsFlyout == SettingsFlyout::FluidScenes && ptIn(fluidBox, shell.mouseX, shell.mouseY)))
            want = SettingsFlyout::FluidScenes;
        else if ((rigidVis && ptIn(rigidTrigger, shell.mouseX, shell.mouseY)) || (shell.settingsFlyout == SettingsFlyout::RigidScenes && ptIn(rigidBox, shell.mouseX, shell.mouseY)))
            want = SettingsFlyout::RigidScenes;
        else if ((gasVis && ptIn(gasTrigger, shell.mouseX, shell.mouseY)) || (shell.settingsFlyout == SettingsFlyout::GasScenes && ptIn(gasBox, shell.mouseX, shell.mouseY)))
            want = SettingsFlyout::GasScenes;
    }
    shell.settingsFlyout = want;

    if (want == SettingsFlyout::Advection) {
        VelocityAdvection mode = engine.config.velocityAdvection;
        FlyItem items[6] = {
            {tr("adv_menu_none"), MenuCmd::AdvectNone, mode == VelocityAdvection::None},
            {tr("adv_menu_fou"), MenuCmd::AdvectFou, mode == VelocityAdvection::FirstOrderUpwind},
            {tr("adv_menu_nsl"), MenuCmd::AdvectNsl, mode == VelocityAdvection::NearestSemiLagrangian},
            {tr("adv_menu_sl"), MenuCmd::AdvectSl, mode == VelocityAdvection::SemiLagrangian},
            {tr("adv_menu_macc"), MenuCmd::AdvectMacc, mode == VelocityAdvection::MacCormack},
            {tr("adv_menu_bfecc"), MenuCmd::AdvectBfecc, mode == VelocityAdvection::BFECC},
        };
        shell.settingsFlyoutRc = advBox;
        drawFlyoutList(dc, shell, advBox, items, 6, 1, 168, 22);
    } else if (want == SettingsFlyout::FluidScenes) {
        FlyItem items[10];
        MenuCmd cmds[10] = {
            MenuCmd::Scene1, MenuCmd::Scene2, MenuCmd::Scene3, MenuCmd::Scene4, MenuCmd::Scene5,
            MenuCmd::Scene6, MenuCmd::Scene7, MenuCmd::Scene8, MenuCmd::Scene9, MenuCmd::Scene10
        };
        char key[24];
        for (int i = 0; i < 10; ++i) {
            std::snprintf(key, sizeof(key), "fluid_scene_%d", i + 1);
            items[i] = {tr(key), cmds[i], false};
        }
        shell.settingsFlyoutRc = fluidBox;
        drawFlyoutList(dc, shell, fluidBox, items, 10, 1, 168, 22);
    } else if (want == SettingsFlyout::RigidScenes) {
        FlyItem items[24];
        MenuCmd cmds[24] = {
            MenuCmd::Rigid1, MenuCmd::Rigid2, MenuCmd::Rigid3, MenuCmd::Rigid4, MenuCmd::Rigid5,
            MenuCmd::Rigid6, MenuCmd::Rigid7, MenuCmd::Rigid8, MenuCmd::Rigid9, MenuCmd::Rigid10,
            MenuCmd::Rigid11, MenuCmd::Rigid12, MenuCmd::Rigid13, MenuCmd::Rigid14, MenuCmd::Rigid15,
            MenuCmd::Rigid16, MenuCmd::Rigid17, MenuCmd::Rigid18, MenuCmd::Rigid19, MenuCmd::Rigid20,
            MenuCmd::Rigid21, MenuCmd::Rigid22, MenuCmd::Rigid23, MenuCmd::Rigid24
        };
        char key[24];
        for (int i = 0; i < 24; ++i) {
            std::snprintf(key, sizeof(key), "rigid_scene_%d", i + 1);
            items[i] = {tr(key), cmds[i], false};
        }
        shell.settingsFlyoutRc = rigidBox;
        drawFlyoutList(dc, shell, rigidBox, items, 24, 2, 148, 22);
    } else if (want == SettingsFlyout::GasScenes) {
        FlyItem items[8];
        MenuCmd cmds[8] = {
            MenuCmd::Gas1, MenuCmd::Gas2, MenuCmd::Gas3, MenuCmd::Gas4,
            MenuCmd::Gas5, MenuCmd::Gas6, MenuCmd::Gas7, MenuCmd::Gas8
        };
        char key[24];
        for (int i = 0; i < 8; ++i) {
            std::snprintf(key, sizeof(key), "gas_scene_%d", i + 1);
            items[i] = {tr(key), cmds[i], false};
        }
        shell.settingsFlyoutRc = gasBox;
        drawFlyoutList(dc, shell, gasBox, items, 8, 1, 168, 22);
    }

    SelectObject(dc, old);
}

char const *tipKeyForCmd(MenuCmd cmd) {
    if (cmd >= MenuCmd::Scene1 && cmd <= MenuCmd::Scene10) return "tip_fluid_scene";
    if (cmd >= MenuCmd::Rigid1 && cmd <= MenuCmd::Rigid24) return "tip_rigid_scene";
    if (cmd >= MenuCmd::Gas1 && cmd <= MenuCmd::Gas8) return "tip_gas_scene";
    if (cmd >= MenuCmd::Speed0 && cmd <= MenuCmd::Speed5) return "tip_speed";
    if (cmd >= MenuCmd::Threads1 && cmd <= MenuCmd::Threads8) return "tip_threads";
    switch (cmd) {
        case MenuCmd::Close: return "tip_close";
        case MenuCmd::Pause: return "tip_pause";
        case MenuCmd::Step: return "tip_step";
        case MenuCmd::Clear: return "tip_clear";
        case MenuCmd::Reset: return "tip_reset";
        case MenuCmd::Slosh: return "tip_slosh";
        case MenuCmd::WalledBorders: return "tip_walled";
        case MenuCmd::Glow: return "tip_glow";
        case MenuCmd::Outlines: return "tip_outlines";
        case MenuCmd::StyleFlat: return "tip_style_flat";
        case MenuCmd::StyleNoisy: return "tip_style_noisy";
        case MenuCmd::StyleDetailed: return "tip_style_detailed";
        case MenuCmd::StyleRealistic: return "tip_style_realistic";
        case MenuCmd::StyleAlpha: return "tip_style_alpha";
        case MenuCmd::Overlay: return "tip_overlay";
        case MenuCmd::MatWood: return "tip_mat_wood";
        case MenuCmd::MatStone: return "tip_mat_stone";
        case MenuCmd::MatGlass: return "tip_mat_glass";
        case MenuCmd::MatMetal: return "tip_mat_metal";
        case MenuCmd::Vorticity: return "tip_vorticity";
        case MenuCmd::Advection:
        case MenuCmd::AdvectNone: case MenuCmd::AdvectFou: case MenuCmd::AdvectNsl:
        case MenuCmd::AdvectSl: case MenuCmd::AdvectMacc: case MenuCmd::AdvectBfecc:
            return "tip_advection";
        case MenuCmd::VorticityMinus: case MenuCmd::VorticityPlus: return "tip_vort_str";
        case MenuCmd::CatchUpMinus: case MenuCmd::CatchUpPlus: return "tip_catch_up";
        case MenuCmd::ThermalIntervalMinus: case MenuCmd::ThermalIntervalPlus: return "tip_heat_every";
        case MenuCmd::RigidGravity: return "tip_gravity";
        case MenuCmd::Residual: return "tip_residual";
        case MenuCmd::QualityLow: return "tip_quality_low";
        case MenuCmd::QualityMed: return "tip_quality_med";
        case MenuCmd::QualityHigh: return "tip_quality_high";
        case MenuCmd::QualityAuto: return "tip_quality_auto";
        case MenuCmd::ThermalOn: return "tip_thermal_on";
        case MenuCmd::ThermalOff: return "tip_thermal_off";
        case MenuCmd::Tension: return "tip_tension";
        case MenuCmd::Spray: return "tip_spray";
        case MenuCmd::SubstepsMinus: case MenuCmd::SubstepsPlus: return "tip_substeps";
        case MenuCmd::LimiterMinus: case MenuCmd::LimiterPlus: return "tip_limiter";
        case MenuCmd::GasOff: return "tip_gas_off";
        case MenuCmd::GasHalf: return "tip_gas_half";
        case MenuCmd::GasFull: return "tip_gas_full";
        case MenuCmd::Hz20: return "tip_hz_20";
        case MenuCmd::Hz30: return "tip_hz_30";
        case MenuCmd::ThreadsAuto: return "tip_threads_auto";
        case MenuCmd::PressureMinus: case MenuCmd::PressurePlus: return "tip_pressure";
        case MenuCmd::BrushMinus: case MenuCmd::BrushPlus: return "tip_brush";
        default: return nullptr;
    }
}

bool hoverTipBody(ShellState const &shell, View const &view, wchar_t const *&body) {
    body = nullptr;
    if (view.settingsOpen && shell.settingsFlyout != SettingsFlyout::None) return false;
    if (view.settingsOpen) {
        for (MenuHit const &hit : shell.menuHits) {
            if (hit.cmd == MenuCmd::None) continue;
            if (!ptIn(hit.rc, shell.mouseX, shell.mouseY)) continue;
            char const *key = hit.tip ? hit.tip : tipKeyForCmd(hit.cmd);
            if (!key) return false;
            body = tr(key);
            return body && body[0] != 0;
        }
        return false;
    }
    if (ViewBarDef const *hoverView = viewBarByHit(shell.hoverId)) {
        body = tr(hoverView->helpKey);
        return body && body[0] != 0;
    }
    HitId id = static_cast<HitId>(shell.hoverId);
    if (id >= HitId::Pal0 && id <= HitId::Pal7) {
        int slot = static_cast<int>(id) - static_cast<int>(HitId::Pal0);
        if (slot >= 0 && slot < shell.elementCount()) {
            body = shell.paletteHint(shell.elementAt(slot));
            return body && body[0] != 0;
        }
    }
    if (id >= HitId::Cat0 && id <= HitId::Cat6) {
        int slot = static_cast<int>(id) - static_cast<int>(HitId::Cat0);
        char const *keys[] = {"tip_cat_tools", "tip_cat_fluids", "tip_cat_solids", "tip_cat_gases",
            "tip_cat_plasma", "tip_cat_energy", "tip_cat_misc"};
        if (slot >= 0 && slot < 7) {
            body = tr(keys[slot]);
            return body && body[0] != 0;
        }
    }
    char const *key = nullptr;
    switch (id) {
        case HitId::SpeedMinus: case HitId::SpeedPlus: key = "tip_speed"; break;
        case HitId::Settings: key = "tip_settings"; break;
        case HitId::Pause: key = "tip_pause"; break;
        case HitId::Clear: key = "tip_clear"; break;
        case HitId::BrushMinus: case HitId::BrushPlus: key = "tip_brush"; break;
        case HitId::Anchored: key = "tip_anchored"; break;
        case HitId::Sleeping: key = "tip_sleeping"; break;
        case HitId::Powder: key = "tip_powder"; break;
        case HitId::GrainMinus: case HitId::GrainPlus: key = "tip_grain"; break;
        default: break;
    }
    if (!key) return false;
    body = tr(key);
    return body && body[0] != 0;
}

void drawHoverTip(HDC dc, ShellState &shell, View const &view) {
    wchar_t const *body = nullptr;
    if (!hoverTipBody(shell, view, body)) return;
    shell.ensureFonts();
    HFONT font = shell.smallFont ? shell.smallFont : static_cast<HFONT>(GetStockObject(DEFAULT_GUI_FONT));
    HFONT old = static_cast<HFONT>(SelectObject(dc, font));
    int const pad = 8;
    int const maxW = 300;
    RECT text{0, 0, maxW, 0};
    DrawTextW(dc, body, -1, &text, DT_CALCRECT | DT_WORDBREAK | DT_LEFT | DT_TOP);
    int w = (text.right - text.left) + pad * 2;
    int h = (text.bottom - text.top) + pad * 2;
    int x = shell.mouseX;
    int y = shell.mouseY;
    RECT box{x, y, x + w, y + h};
    fillRect(dc, box, RGB(32, 32, 32));
    frameRect(dc, box, kBorder);
    RECT inner = insetLike(box, pad);
    drawLabel(dc, inner, body, kText, DT_LEFT | DT_TOP | DT_WORDBREAK);
    SelectObject(dc, old);
}

} // namespace ui
