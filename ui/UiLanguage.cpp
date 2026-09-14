#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "ui/UiLanguage.h"

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace ui {
namespace {

std::map<std::string, std::wstring> gStrings;

std::wstring utf8ToWide(std::string const &utf8) {
    if (utf8.empty()) return L"";
    int n = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), static_cast<int>(utf8.size()), nullptr, 0);
    if (n <= 0) return L"";
    std::wstring out(static_cast<size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), static_cast<int>(utf8.size()), out.data(), n);
    return out;
}

void set(char const *key, wchar_t const *value) {
    gStrings[key] = value;
}

void seedEnglish() {
    set("window_title", L"PIPACE - Pixel Physics Sandbox");
    set("settings", L"SETTINGS");
    set("pause", L"PAUSE");
    set("resume", L"RESUME");
    set("clear", L"CLEAR");
    set("pause_word", L"Pause");
    set("resume_word", L"Resume");
    set("clear_word", L"Clear");
    set("step", L"Step");
    set("reset", L"Reset");
    set("slosh", L"Slosh");
    set("btn_walled", L"Walled borders");
    set("close", L"Close");
    set("settings_title", L"Settings");

    set("cat_tools", L"TOOLS");
    set("cat_fluids", L"LIQUIDS");
    set("cat_liquids", L"LIQUIDS");
    set("cat_solids", L"SOLIDS");
    set("cat_gases", L"GASES");
    set("cat_plasma", L"PLASMA");
    set("cat_energy", L"ENERGY");
    set("cat_misc", L"MISC");

    set("el_water", L"WATER");
    set("el_honey", L"HONEY");
    set("el_wood", L"WOOD");
    set("el_stone", L"STONE");
    set("el_glass", L"GLASS");
    set("el_metal", L"METAL");
    set("el_erase", L"ERASE");
    set("el_grab", L"GRAB");
    set("el_wall", L"WALL");
    set("el_heat", L"HEAT");
    set("el_cool", L"COOL");
    set("el_pressurize", L"PRESSURIZE");
    set("el_depressurize", L"DEPRESSURIZE");

    set("hint_water", L"Free-surface liquid. Paint into the world. Shift+drag previews a straight line; release to place.");
    set("hint_honey", L"Thicker, denser liquid. Mixes with water. Dye it from Properties. Shift+drag previews a straight line; release to place.");
    set("hint_wood", L"Dynamic solid. Anchored stays pinned. Sleeping starts pinned until something hits it. Absorbs water and weakens when wet. Shift+drag draws a straight line.");
    set("hint_stone", L"Heavier dynamic solid. Anchored stays pinned. Sleeping wakes when disturbed. Shift+drag draws a straight line.");
    set("hint_glass", L"Hard, brittle solid. Shatters from hard impacts. Does not absorb water. Shift+drag draws a straight line.");
    set("hint_metal", L"Dense, tough solid. Survives drops that break glass. Does not absorb water. Shift+drag draws a straight line.");
    set("hint_erase", L"Removes liquid, walls, and rigid pixels. Shift+drag previews a straight erase; release to apply.");
    set("hint_grab", L"Click a loose body and drag. Force is applied at the grab point, so off-center grabs rotate it. Shift for a stronger hold.");
    set("hint_wall", L"Paints static world walls. Anchored wood/stone also stay put. Shift+drag previews a straight line; release to place.");
    set("hint_heat", L"Adds thermal energy under the brush. Works on liquid, solids, and air. Not a substance.");
    set("hint_cool", L"Removes thermal energy under the brush. Works on liquid, solids, and air. Not a substance.");
    set("hint_pressurize", L"Adds air under the brush, raising gas pressure. Does not add liquid. Not a substance.");
    set("hint_depressurize", L"Removes air under the brush, lowering gas pressure. Does not remove liquid. Not a substance.");
    set("hint_none", L"Nothing in this tab yet.");

    set("view_normal", L"NORM");
    set("view_fill", L"FILL");
    set("view_pressure", L"LIQP");
    set("view_velocity", L"LVEL");
    set("view_divergence", L"LDIV");
    set("view_chunks", L"CHNK");
    set("view_rigid", L"RGDN");
    set("view_gas_pressure", L"GASP");
    set("view_gas_amount", L"GASA");
    set("view_gas_velocity", L"GASV");

    set("bar_norm", L"NORM");
    set("bar_chnk", L"CHNK");
    set("bar_fill", L"FILL");
    set("bar_liqp", L"LIQP");
    set("bar_lvel", L"LVEL");
    set("bar_ldiv", L"LDIV");
    set("bar_temp", L"TEMP");
    set("bar_rgdn", L"RGDN");
    set("bar_rgdo", L"RGDO");
    set("bar_mois", L"MOIS");
    set("bar_gasp", L"GASP");
    set("bar_gasa", L"GASA");
    set("bar_gasv", L"GASV");

    set("view_title_norm", L"NORM - Normal view");
    set("view_help_norm", L"Usual world render. Air is not drawn. Liquids, walls, and bodies look as they do in play.");
    set("view_title_chnk", L"CHNK - Chunk activity");
    set("view_help_chnk", L"Shows which 16x16 chunks the liquid solver is currently updating. Dark cells are asleep.");
    set("view_title_fill", L"FILL - Liquid fill");
    set("view_help_fill", L"How full each cell is of liquid, from empty to packed. Useful for spotting thin films and leftover crumbs.");
    set("view_title_liqp", L"LIQP - Liquid pressure");
    set("view_help_liqp", L"Pressure inside liquid cells from the fluid solver. Not the same as air pressure.");
    set("view_title_lvel", L"LVEL - Liquid velocity");
    set("view_help_lvel", L"How fast the liquid is moving. Brighter means faster flow.");
    set("view_title_ldiv", L"LDIV - Liquid divergence");
    set("view_help_ldiv", L"Compression / expansion leftover in the liquid velocity field. Helps debug the pressure solve.");
    set("view_title_temp", L"TEMP - Temperature");
    set("view_help_temp", L"Ambient-centered temperature overlay. Quiet near 293 K. Blue is colder, yellow/red is hotter. Ambient air stays faint. Does not change the simulation.");
    set("view_title_rgdn", L"RGDN - Rigid occupancy");
    set("view_help_rgdn", L"Highlights rigid-body pixels and pending paint. Use this to see exactly which cells a body occupies.");
    set("view_title_rgdo", L"RGDO - Rigid overlay");
    set("view_help_rgdo", L"Toggle. Draws debug marks on rigid bodies on top of whatever view is already showing. Can stay on while you use other views.");
    set("view_title_mois", L"MOIS - Rigid moisture");
    set("view_help_mois", L"Diagnostic wetness field. Dry absorbent material is dark. Wetter pixels shift toward cyan. Non-absorbent solids stay gray. Does not change the simulation.");
    set("view_title_gasp", L"GASP - Gas pressure");
    set("view_help_gasp", L"Air pressure in atmospheres. Dark is vacuum. Mid gray-blue is about 1 atm. Bright is over-pressure.");
    set("view_title_gasa", L"GASA - Gas amount");
    set("view_help_gasa", L"How much Air is in each cell. Vacuum is dark. More Air reads brighter, independent of the cell's free volume.");
    set("view_title_gasv", L"GASV - Gas velocity");
    set("view_help_gasv", L"Which way the air is moving. Quiet uniform Air stays dark.");

    set("info_brush", L"Brush Size: {0} px");
    set("info_line", L"Shift+LMB: straight line");
    set("info_fill", L"Fill: {0}%");
    set("info_cell", L"Cell [{0},{1}]");
    set("info_view", L"View: {0}");

    set("console_hint", L"see console or enter command...");
    set("search_hint", L"search for an element");
    set("none_yet", L"none yet");

    set("properties", L"PROPERTIES");
    set("inspector", L"INSPECTOR");
    set("prop_intro", L"How this placement behaves");
    set("prop_brush", L"Brush Size: {0} px");
    set("prop_anchored", L"Anchored");
    set("prop_sleeping", L"Sleeping");
    set("prop_powder", L"Powder");
    set("prop_grain", L"Particle Size: {0} px");
    set("prop_on_pinned", L"on (stays pinned)");
    set("prop_on_wakes", L"on (wakes on hit)");
    set("prop_on", L"on");
    set("prop_off", L"off");
    set("prop_solids_note", L"Anchored never moves. Sleeping wakes when disturbed. Powder paints loose grains instead of one body.");
    set("prop_material_wood", L"Material: wood");
    set("prop_material_stone", L"Material: stone");
    set("prop_material_glass", L"Material: glass");
    set("prop_material_metal", L"Material: metal");
    set("prop_substance_water", L"Substance: water");
    set("prop_substance_honey", L"Substance: honey");
    set("prop_fill_paint", L"Fill: paint to 100%");
    set("prop_dye_mode", L"Paint: dye only (no new fill)");
    set("prop_dye", L"Dye");
    set("prop_dye_only", L"Dye only");
    set("prop_dye_strength", L"Dye strength: {0}");
    set("prop_dye_note", L"Dye tints existing liquid. Painting dye into empty air or solid does nothing.");
    set("prop_tool_eraser", L"Tool: eraser");
    set("prop_tool_grab", L"Tool: grab (spring hold)");
    set("prop_tool_wall", L"Tool: static wall");
    set("prop_tool_heat", L"Tool: heat (energy)");
    set("prop_tool_cool", L"Tool: cool (energy)");
    set("prop_tool_pressurize", L"Tool: pressurize (air)");
    set("prop_tool_depressurize", L"Tool: depressurize (air)");
    set("prop_power", L"Power: {0}");
    set("prop_empty_tab", L"No elements in this tab.");

    set("ins_nothing", L"Nothing selected.");
    set("ins_hover_1", L"Hover the canvas to");
    set("ins_hover_2", L"observe a cell.");
    set("ins_kind_wall", L"Kind: static wall");
    set("ins_kind_rigid", L"Kind: rigid body");
    set("ins_kind_liquid", L"Kind: liquid / void");
    set("ins_material", L"Material: {0}");
    set("ins_id", L"Id: {0}");
    set("ins_mass", L"Mass: {0}");
    set("ins_pos", L"Pos: {0}, {1}");
    set("ins_vel", L"Vel: {0}, {1}");
    set("ins_omega", L"Omega: {0}");
    set("ins_supported", L"Supported: {0}");
    set("ins_anchored", L"Anchored: {0}");
    set("ins_sleeping", L"Sleeping: {0}");
    set("ins_at_rest", L"At rest: {0}");
    set("ins_pen", L"Pen: {0}");
    set("ins_jnjt", L"Jn/Jt: {0} / {1}");
    set("ins_poscorr", L"PosCorr: {0}, {1}");
    set("ins_damage", L"Damage: {0}");
    set("ins_component", L"Component: {0}");
    set("ins_bond_body", L"Bonds: dmg {0}, broken {1}");
    set("ins_mat_damage", L"Material dmg: {0}");
    set("ins_bond_damage", L"Bond dmg: {0}");
    set("ins_broken", L"Broken neighbors: {0}");
    set("ins_strength", L"Effective strength: {0}");
    set("ins_crack", L"Crack: {0}");
    set("ins_saturation", L"Saturation: {0}");
    set("ins_bond_broken", L"Bond broken at hover");
    set("ins_moisture", L"Absorbed (body): {0} ({1})");
    set("ins_local_moisture", L"Moisture: {0}");
    set("ins_moisture_cap", L"Capacity: {0}");
    set("ins_liquid_total", L"Liquid free/abs/splash/total {0} / {1} / {2} / {3}");
    set("ins_pending_drip", L"Pending drip: {0}");
    set("ins_fill", L"Fill: {0}");
    set("ins_comp_water", L"{0}% water");
    set("ins_comp_honey", L"{0}% honey");
    set("ins_dye", L"Dye: {0}%");
    set("ins_temp", L"Temperature: {0} K");
    set("ins_temp_header", L"TEMP");
    set("ins_dT", L"dT: {0} K");
    set("ins_thermal_energy", L"Heat: {0} J");
    set("ins_thermal_active", L"Thermally active: {0}");
    set("ins_mat_water", L"Water");
    set("ins_mat_honey", L"Honey");
    set("ins_mat_air", L"Air");
    set("ins_mat_none", L"(none)");
    set("ins_substance", L"Substance: {0}");
    set("ins_substance_id", L"Substance ID: {0}");
    set("ins_density", L"Density: {0}");
    set("ins_viscosity", L"Viscosity: {0}");
    set("ins_specific_heat", L"Specific heat: {0}");
    set("ins_melting_point", L"Melting point: {0} K");
    set("ins_boiling_point", L"Boiling point: {0} K");
    set("ins_pressure", L"Pressure: {0}");
    set("ins_surface", L"Surface: {0}");
    set("ins_volume", L"Volume {0} / {1}");
    set("ins_gas_material", L"Material: Air");
    set("ins_gas_amount", L"Gas Amount: {0}");
    set("ins_gas_vol", L"Gas Vol: {0}");
    set("ins_gas_pa", L"Pressure: {0} Pa");
    set("ins_gas_atm", L"Pressure: {0} atm");
    set("ins_gas_vel", L"Gas Vel: {0}, {1}");
    set("ins_gas_comp", L"Composition: Air 100%");
    set("ins_gas_world", L"World gas {0} / {1}");
    set("yes", L"yes");
    set("no", L"no");
    set("wood", L"wood");
    set("stone", L"stone");
    set("glass", L"glass");
    set("metal", L"metal");
    set("wood_btn", L"Wood");
    set("stone_btn", L"Stone");
    set("glass_btn", L"Glass");
    set("metal_btn", L"Metal");

    set("sec_world", L"World");
    set("sec_view", L"View");
    set("sec_look", L"Appearance");
    set("lbl_style", L"Style");
    set("lbl_effects", L"Effects");
    set("style_flat", L"Flat");
    set("style_noisy", L"Noisy");
    set("style_detailed", L"Detail");
    set("style_realistic", L"Real");
    set("style_alpha", L"Alpha");
    set("style_legacy", L"Legacy");
    set("noise_str", L"Noise {0}%");
    set("sec_quality", L"Quality");
    set("sec_solver", L"Solver");
    set("sec_gas_sim", L"Air solver");
    set("sec_sim_rate", L"Sim rate");
    set("sec_threads", L"Simulation threads");
    set("sec_speed", L"Speed");
    set("sec_fluid_scenes", L"Fluid scenes");
    set("sec_rigid_scenes", L"Rigid scenes");
    set("sec_gas_scenes", L"Gas scenes");
    set("sec_scenes", L"Test scenes");
    set("menu_fluid_scenes", L"Fluid scenes");
    set("menu_rigid_scenes", L"Rigid scenes");
    set("menu_gas_scenes", L"Gas scenes");
    set("btn_adv_menu", L"{0}");
    set("adv_trigger", L"Adv: {0}");
    set("adv_menu_none", L"None");
    set("adv_menu_fou", L"First-order upwind");
    set("adv_menu_nsl", L"Nearest SL");
    set("adv_menu_sl", L"Semi-Lagrangian");
    set("adv_menu_macc", L"MacCormack");
    set("adv_menu_bfecc", L"BFECC");
    set("btn_gravity", L"Gravity");
    set("heat_every", L"Heat / {0}");
    set("catch_up", L"Catch-up {0}");
    set("vort_str", L"Vort {0}");

    set("btn_normal", L"Normal");
    set("btn_fill", L"Fill");
    set("btn_pressure", L"Liq. P");
    set("btn_velocity", L"Velocity");
    set("btn_divergence", L"Divergence");
    set("btn_chunks", L"Chunks");
    set("btn_rigid", L"Rigid");
    set("btn_gas_pressure", L"Gas P");
    set("btn_gas_amount", L"Gas Amt");
    set("btn_gas_velocity", L"Gas Vel");
    set("btn_glow_liquids", L"Glowing Liquids");
    set("btn_outlines", L"Outlines");
    set("btn_rigid_overlay", L"Rigid overlay");
    set("btn_vorticity", L"Vorticity");
    set("btn_adv_none", L"Adv None");
    set("btn_adv_fou", L"Adv FOU");
    set("btn_adv_nsl", L"Adv NSL");
    set("btn_adv_sl", L"Adv SL");
    set("btn_adv_macc", L"Adv MacC");
    set("btn_adv_bfecc", L"Adv BFECC");
    set("btn_res_merge", L"Res. merge");
    set("btn_tension", L"Tension");
    set("btn_spray", L"Spray");
    set("quality_low", L"Low");
    set("quality_med", L"Medium");
    set("quality_high", L"High");
    set("quality_auto", L"Auto");
    set("thermal_on", L"Thermal ON");
    set("thermal_off", L"Thermal OFF");
    set("quality_using", L"{0}");
    set("quality_auto_using", L"Auto (using {0})");
    set("p_iter", L"P iter {0}");
    set("substeps", L"Substeps {0}");
    set("limiter", L"Limiter {0}");
    set("brush", L"Brush {0}");
    set("gas_off", L"Off");
    set("gas_half", L"Half");
    set("gas_full", L"Full");
    set("hz_20", L"20 Hz");
    set("hz_30", L"30 Hz");
    set("threads_auto", L"Auto");
    set("threads_auto_using", L"Auto (using {0})");
    set("threads_using", L"Using {0}");
    set("threads_parallel", L"  pressure parallel");
    set("threads_serial", L"  pressure serial");

    set("fluid_scene_1", L"1 Still pool");
    set("fluid_scene_2", L"2 Vessels");
    set("fluid_scene_3", L"3 Dam break");
    set("fluid_scene_4", L"4 Waterfall");
    set("fluid_scene_5", L"5 Obstacle");
    set("fluid_scene_6", L"6 Droplet");
    set("fluid_scene_7", L"7 Puddle");
    set("fluid_scene_8", L"8 Incline");
    set("fluid_scene_9", L"9 Nozzle");
    set("fluid_scene_10", L"10 Channel");

    set("rigid_scene_1", L"1 Fall block");
    set("rigid_scene_2", L"2 Irregular");
    set("rigid_scene_3", L"3 Hollow");
    set("rigid_scene_4", L"4 Rotation");
    set("rigid_scene_5", L"5 Into pool");
    set("rigid_scene_6", L"6 Float");
    set("rigid_scene_7", L"7 Sink");
    set("rigid_scene_8", L"8 Cup");
    set("rigid_scene_9", L"9 Thin geo");
    set("rigid_scene_10", L"10 Conserve");
    set("rigid_scene_11", L"11 Flat round");
    set("rigid_scene_12", L"12 Slope");
    set("rigid_scene_13", L"13 Rev slope");
    set("rigid_scene_14", L"14 Peanut");
    set("rigid_scene_15", L"15 No grav");
    set("rigid_scene_16", L"16 Slight pen");
    set("rigid_scene_17", L"17 Off-center");
    set("rigid_scene_18", L"18 Glass drop");
    set("rigid_scene_19", L"19 Metal drop");
    set("rigid_scene_20", L"20 Glass/metal");
    set("rigid_scene_21", L"21 Glass beam");
    set("rigid_scene_22", L"22 Wet wood");
    set("rigid_scene_23", L"23 Anchored");
    set("rigid_scene_24", L"24 Spin glass");
    set("rigid_scene_25", L"25 Moisture soak");

    set("gas_scene_1", L"A Uniform air");
    set("gas_scene_2", L"B Sealed box");
    set("gas_scene_3", L"C Two rooms");
    set("gas_scene_4", L"D Solid block");
    set("gas_scene_5", L"E Compress");
    set("gas_scene_6", L"F Expand");
    set("gas_scene_7", L"G Fall block");
    set("gas_scene_8", L"H Sealed air");

    set("log_ready", L"[PHYS] ready");
    set("log_cleared", L"[PHYS] world cleared");
    set("log_reset", L"[PHYS] world reset");
    set("log_fluid_scene", L"[PHYS] fluid scene loaded");
    set("log_rigid_scene", L"[PHYS] rigid scene loaded");
    set("log_gas_scene", L"[PHYS] gas scene loaded");
    set("log_rigid_deleted", L"[PHYS] rigid body deleted");
    set("log_grab", L"[PHYS] grabbing body");
    set("log_no_command", L"[WARN] no command interpreter yet");

    set("tip_close", L"Close the settings overlay. Esc does the same.");
    set("tip_pause", L"Pause or resume the world. Space also toggles this.");
    set("tip_step", L"While paused, advance exactly one physics tick.");
    set("tip_clear", L"Wipe liquid, air, walls, and bodies from the map.");
    set("tip_reset", L"Reload the default walls and an empty basin.");
    set("tip_slosh", L"Give the liquid a sideways shove to test sloshing.");
    set("tip_walled", L"Close the grid edges. Off: liquid, spray, air, and solids can leave into the void. On: the borders act as walls.");
    set("tip_glow", L"Adds a soft glow around liquids.");
    set("tip_outlines", L"Draw outlines around liquids and solids.");
    set("tip_style_flat", L"One flat color per material.");
    set("tip_style_noisy", L"Flat colors with subtle texture.");
    set("tip_style_detailed", L"Adds depth and shape shading.");
    set("tip_style_realistic", L"Shows richer depth, flow, and material detail.");
    set("tip_style_alpha", L"Flat colors with opacity based on amount.");
    set("tip_style_legacy", L"Classic PIPACE shading with depth, bright surfaces, motion, and noise.");
    set("tip_noise", L"Adjusts texture variation in noise-based appearance styles.");
    set("tip_overlay", L"Draw extra rigid-body debug marks on top of the current view.");
    set("tip_mat_wood", L"Set the next drawn body to wood: light, absorbent, weakens when wet.");
    set("tip_mat_stone", L"Set the next drawn body to stone: heavier than wood.");
    set("tip_mat_glass", L"Set the next drawn body to glass: hard and brittle, shatters from hard hits.");
    set("tip_mat_metal", L"Set the next drawn body to metal: dense and tough.");
    set("tip_vorticity", L"Add a little swirl back into the liquid. Off by default. Extra cost when on.");
    set("tip_advection", L"How velocity is carried with the flow. Semi-Lagrangian is the usual choice. BFECC/MacCormack cost more.");
    set("tip_residual", L"Merge leftover thin liquid into neighbors. Keeps conservation. Turn off only to inspect crumbs.");
    set("tip_quality_low", L"Cheaper sim for slow machines: fewer pressure/substep/limiter passes, no air, 20 Hz. Appearance stays as you set it.");
    set("tip_quality_med", L"Default fidelity. Same world laws as Low/High, with the usual solver budget.");
    set("tip_quality_high", L"More pressure iterations. Same gravity and conservation, slightly stiffer water, more CPU.");
    set("tip_quality_auto", L"Starts at Medium. If a tick stays overloaded, drops to Low knobs until the machine catches up.");
    set("tip_thermal_on", L"Run heat conduction and interface exchange. Stored temperatures are kept if you turn this off.");
    set("tip_thermal_off", L"Skip thermal solver work. The rest of the simulation continues. Temperatures stay as they are.");
    set("tip_tension", L"Surface tension at free surfaces. Makes beads and menisci. Skip it to save the surface pass.");
    set("tip_spray", L"Spawn splash droplets from energetic free surfaces. Off saves particle work.");
    set("tip_substeps", L"Max CFL substeps per tick. Lower is cheaper; fast jets may clip more.");
    set("tip_limiter", L"How hard the donor/receiver flux limiter works. Lower is cheaper on big pools.");
    set("tip_gas_off", L"Skip the air solver and air forces on bodies. Fastest. Gas scenes will turn air back on.");
    set("tip_gas_half", L"Update air every other tick. Cheaper than Full, still applies the last pressure field.");
    set("tip_gas_full", L"Run air every physics tick.");
    set("tip_hz_20", L"20 physics ticks per second and no catch-up. Smoother on weak CPUs; the world runs a bit slower.");
    set("tip_hz_30", L"30 physics ticks per second (default). Can catch up one extra tick if a frame runs long.");
    set("tip_threads_auto", L"Pick a worker count for this grid. On 200x120 this stays at 1; extra cores usually slow pressure down.");
    set("tip_threads", L"Force this many simulation workers. On the default grid, more than 1 is often slower.");
    set("tip_pressure", L"Cap on liquid pressure-solver iterations. Lower is cheaper and water compresses a little more.");
    set("tip_brush", L"Paint radius in pixels.");
    set("tip_speed", L"How fast wall-clock time feeds the simulator. 1x is realtime at the chosen sim rate.");
    set("tip_settings", L"Open or close the settings window. Drag the title bar to move it. Esc closes it.");
    set("tip_gravity", L"Turn gravity on rigid bodies on or off. Liquid gravity is unchanged.");
    set("tip_vort_str", L"Strength of vorticity confinement when that toggle is on. Does not create new swirl by itself.");
    set("tip_catch_up", L"How many extra physics ticks may run if a frame is late. 1 is gentlest on slow machines.");
    set("tip_heat_every", L"Thermal conduction runs every N physics ticks. Higher skips more heat work.");
    set("tip_anchored", L"Newly drawn solids stay pinned in the world and never fall. Sleeping cannot be on at the same time.");
    set("tip_sleeping", L"Newly drawn solids start asleep (pinned) until a hit, water, or the eraser wakes them. Anchored cannot be on at the same time.");
    set("tip_powder", L"Paint loose grains instead of one connected solid. Particle size is the grain width in pixels.");
    set("tip_grain", L"Width of each powder grain in pixels. 1 is dust. Larger grains are chunkier and cheaper.");
    set("tip_fluid_scene", L"Load this prepared fluid test world. Replaces the current liquid and walls.");
    set("tip_rigid_scene", L"Load this prepared rigid-body test. Replaces current bodies.");
    set("tip_gas_scene", L"Load this prepared air test. Turns the air solver on if it was off.");
    set("tip_cat_tools", L"Eraser and grab. Tools change how you edit the world, not what material you spawn.");
    set("tip_cat_fluids", L"Paintable liquids. Water is the current free-surface engine.");
    set("tip_cat_solids", L"Drawable rigid bodies. Wood, stone, glass, and metal share the same solver with different density and strength.");
    set("tip_cat_gases", L"Gas elements are not paintable yet. Air still fills the map and has its own views and scenes.");
    set("tip_cat_plasma", L"Nothing in this tab yet.");
    set("tip_cat_energy", L"Heat, cool, pressurize, and depressurize. These add or remove energy or air; they are not substances.");
    set("tip_cat_misc", L"Static world walls. These do not move.");
}

std::string unescape(std::string const &s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size()) {
            char n = s[++i];
            if (n == 'n') out.push_back('\n');
            else if (n == 't') out.push_back('\t');
            else if (n == '"') out.push_back('"');
            else if (n == '\\') out.push_back('\\');
            else out.push_back(n);
        } else out.push_back(s[i]);
    }
    return out;
}

bool parseFlatJson(std::string const &text) {
    size_t i = 0;
    auto skip = [&]() {
        while (i < text.size() && (text[i] == ' ' || text[i] == '\t' || text[i] == '\r' || text[i] == '\n')) ++i;
    };
    auto parseString = [&](std::string &out) -> bool {
        skip();
        if (i >= text.size() || text[i] != '"') return false;
        ++i;
        std::string raw;
        while (i < text.size() && text[i] != '"') {
            if (text[i] == '\\' && i + 1 < text.size()) {
                raw.push_back(text[i]);
                raw.push_back(text[i + 1]);
                i += 2;
            } else raw.push_back(text[i++]);
        }
        if (i >= text.size()) return false;
        ++i;
        out = unescape(raw);
        return true;
    };
    skip();
    if (i >= text.size() || text[i] != '{') return false;
    ++i;
    for (;;) {
        skip();
        if (i < text.size() && text[i] == '}') break;
        std::string key, value;
        if (!parseString(key)) return false;
        skip();
        if (i >= text.size() || text[i] != ':') return false;
        ++i;
        if (!parseString(value)) return false;
        gStrings[key] = utf8ToWide(value);
        skip();
        if (i < text.size() && text[i] == ',') { ++i; continue; }
        if (i < text.size() && text[i] == '}') break;
        return false;
    }
    return true;
}

std::wstring exeDir() {
    wchar_t buf[MAX_PATH]{};
    GetModuleFileNameW(nullptr, buf, MAX_PATH);
    std::wstring path = buf;
    size_t slash = path.find_last_of(L"\\/");
    if (slash == std::wstring::npos) return L".";
    return path.substr(0, slash);
}

bool tryLoadFile(std::wstring const &path) {
    FILE *f = _wfopen(path.c_str(), L"rb");
    if (!f) return false;
    std::string text;
    char buf[4096];
    size_t n = 0;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) text.append(buf, n);
    fclose(f);
    if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF) text.erase(0, 3);
    return parseFlatJson(text);
}

} // namespace

bool loadLanguage(char const *name) {
    seedEnglish();
    if (!name || !name[0]) return false;
    std::wstring file = utf8ToWide(std::string(name) + ".json");
    std::wstring root = exeDir();
    std::wstring candidates[] = {
        L"languages\\" + file,
        root + L"\\..\\languages\\" + file,
        root + L"\\languages\\" + file,
        root + L"\\..\\..\\languages\\" + file
    };
    for (std::wstring const &p : candidates)
        if (tryLoadFile(p)) return true;
    return false;
}

wchar_t const *tr(char const *key) {
    if (gStrings.empty()) seedEnglish();
    auto it = gStrings.find(key);
    if (it == gStrings.end()) return L"";
    return it->second.c_str();
}

std::wstring trf(char const *key) {
    return tr(key);
}

std::wstring replaceAll(std::wstring s, std::wstring const &from, std::wstring const &to) {
    size_t at = 0;
    while ((at = s.find(from, at)) != std::wstring::npos) {
        s.replace(at, from.size(), to);
        at += to.size();
    }
    return s;
}

std::wstring trf(char const *key, std::wstring const &a0) {
    return replaceAll(tr(key), L"{0}", a0);
}

std::wstring trf(char const *key, std::wstring const &a0, std::wstring const &a1) {
    return replaceAll(replaceAll(tr(key), L"{0}", a0), L"{1}", a1);
}

std::wstring trf(char const *key, std::wstring const &a0, std::wstring const &a1, std::wstring const &a2) {
    return replaceAll(trf(key, a0, a1), L"{2}", a2);
}

std::wstring trf(char const *key, std::wstring const &a0, std::wstring const &a1, std::wstring const &a2, std::wstring const &a3) {
    return replaceAll(trf(key, a0, a1, a2), L"{3}", a3);
}

} // namespace ui
