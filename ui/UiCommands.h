#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <cstdint>

namespace ui {

enum class MenuCmd : int {
    None = 0,
    Close,
    Pause,
    Step,
    Clear,
    Reset,
    Slosh,
    ViewNormal, ViewFill, ViewPressure, ViewVelocity, ViewDivergence, ViewChunks, ViewRigid, ViewMoisture,
    ViewGasPressure, ViewGasAmount, ViewGasVelocity,
    Glow,
    Outlines,
    StyleFlat,
    StyleNoisy,
    StyleDetailed,
    StyleRealistic,
    StyleAlpha,
    StyleLegacy,
    Overlay,
    MatWood,
    MatStone,
    MatGlass,
    MatMetal,
    Vorticity,
    Advection,
    Residual,
    QualityLow, QualityMed, QualityHigh, QualityAuto, QualityCustom,
    Tension,
    Spray,
    SubstepsMinus,
    SubstepsPlus,
    LimiterMinus,
    LimiterPlus,
    GasOff, GasHalf, GasFull,
    Hz20, Hz30,
    WalledBorders,
    ThermalOn,
    ThermalOff,
    ThreadsAuto, Threads1, Threads2, Threads4, Threads6, Threads8,
    PressureMinus,
    PressurePlus,
    BrushMinus,
    BrushPlus,
    Speed0, Speed1, Speed2, Speed3, Speed4, Speed5,
    Scene1, Scene2, Scene3, Scene4, Scene5, Scene6, Scene7, Scene8, Scene9, Scene10,
    Rigid1, Rigid2, Rigid3, Rigid4, Rigid5, Rigid6, Rigid7, Rigid8, Rigid9, Rigid10,
    Rigid11, Rigid12, Rigid13, Rigid14, Rigid15, Rigid16, Rigid17,
    Rigid18, Rigid19, Rigid20, Rigid21, Rigid22, Rigid23, Rigid24, Rigid25,
    Gas1, Gas2, Gas3, Gas4, Gas5, Gas6, Gas7, Gas8,
    AdvectNone, AdvectFou, AdvectNsl, AdvectSl, AdvectMacc, AdvectBfecc,
    VorticityMinus, VorticityPlus,
    CatchUpMinus, CatchUpPlus,
    ThermalIntervalMinus, ThermalIntervalPlus,
    RigidGravity,
    NoiseMinus, NoisePlus,
    CustomFluidP, CustomFluidB, CustomFluidA,
    CustomGasP, CustomGasB, CustomGasA,
    CustomRigidP, CustomRigidB, CustomRigidA,
    CustomThermalP, CustomThermalB, CustomThermalA,
    CustomChemP, CustomChemB, CustomChemA,
    CustomPhaseP, CustomPhaseB, CustomPhaseA
};

enum class SettingsMouseResult : uint8_t { Miss = 0, Consume, Command, Drag };

enum class SettingsFlyout : uint8_t {
    None = 0,
    Advection,
    FluidScenes,
    RigidScenes,
    GasScenes
};

enum class Category : int {
    Tools = 0,
    Fluids,
    Solids,
    Gases,
    Plasma,
    Energy,
    Misc,
    Count
};

enum class PaletteId : int {
    None = -1,
    Water = 0,
    Honey,
    Wood,
    Stone,
    Glass,
    Metal,
    Erase,
    Grab,
    Brush,
    Touch,
    Wall,
    Heat,
    Cool,
    Pressurize,
    Depressurize,
    Hydrogen,
    Oxygen,
    Carbon,
    CarbonDioxide
};

// Floating-window / canvas editing mode. Separate from selected content
// (Water, Carbon, Oxygen, …). PaletteId::Brush is a Tools item, not a material.
enum class EditTool : uint8_t {
    Brush = 0,
    Erase,
    Grab,
    Touch,
    Heat
};

inline PaletteId paletteForEditTool(EditTool tool) {
    switch (tool) {
        case EditTool::Erase: return PaletteId::Erase;
        case EditTool::Grab: return PaletteId::Grab;
        case EditTool::Touch: return PaletteId::Touch;
        case EditTool::Heat: return PaletteId::Heat;
        case EditTool::Brush:
        default: return PaletteId::Brush;
    }
}

inline bool isPlaceablePalette(PaletteId id) {
    switch (id) {
        case PaletteId::Water:
        case PaletteId::Honey:
        case PaletteId::Wood:
        case PaletteId::Stone:
        case PaletteId::Glass:
        case PaletteId::Metal:
        case PaletteId::Carbon:
        case PaletteId::Wall:
        case PaletteId::Hydrogen:
        case PaletteId::Oxygen:
        case PaletteId::CarbonDioxide:
            return true;
        default:
            return false;
    }
}

inline bool isEnergyPalette(PaletteId id) {
    return id == PaletteId::Heat || id == PaletteId::Cool
        || id == PaletteId::Pressurize || id == PaletteId::Depressurize;
}

enum class HitId : int {
    None = 0,
    SpeedMinus,
    SpeedPlus,
    Settings,
    Pause,
    Clear,
    Pal0,
    Pal1,
    Pal2,
    Pal3,
    Pal4,
    Pal5,
    Pal6,
    Pal7,
    Cat0,
    Cat1,
    Cat2,
    Cat3,
    Cat4,
    Cat5,
    Cat6,
    BrushMinus,
    BrushPlus,
    Anchored,
    Sleeping,
    Powder,
    GrainMinus,
    GrainPlus,
    PowerMinus,
    PowerPlus,
    Search,
    Console,
    SettingsHit,
    ViewNorm,
    ViewChnk,
    ViewFill,
    ViewLiqp,
    ViewLvel,
    ViewLdiv,
    ViewRgdn,
    ViewRgdo,
    ViewMois,
    ViewGasp,
    ViewGasa,
    ViewGasv,
    ViewTemp,
    Dye0,
    Dye1,
    Dye2,
    Dye3,
    Dye4,
    Dye5,
    Dye6,
    Dye7,
    DyeOnly,
    DyeStrMinus,
    DyeStrPlus
};

struct MenuHit {
    RECT rc{};
    MenuCmd cmd = MenuCmd::None;
    char const *tip = nullptr;
};

} // namespace ui
