#pragma once

#include <cstdint>
#include <vector>

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

// Local, conservative reaction execution. Scans liquid and gas-occupied cells.
// Physical reactions convert storage -> mass -> moles using chemical.molarMass.
// Extent 1 is one mole of the written reaction. Synthetic inventories keep a
// diagnostic unit path. HomogeneousCell stays liquid/gas in one cell.
// SolidGasSurface uses one rigid source pixel plus an adjacent gas cell.
//
// activity[] is rendering feedback only (recent local reaction intensity).
// It is not chemistry state and is not advected with gas.

struct ReactionEngine {
    std::vector<float> activity;
    int reactedCellsLastTick = 0;
    float extentLastTick = 0.0f;
    float heatReleasedLastTick = 0.0f;

    ReactionEngine();

    void simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
        ThermalEngine &thermal, float dt);
    void clearActivity();
    bool hasVisibleActivity() const { return !activityCells.empty(); }

private:
    std::vector<int> activityCells;
    void decayActivity(float dt);
    void addActivity(int index, float heatJ);
};

// Diagnostic-only. Uses synthetic inventories, not SubstanceRegistry / UI.
// Writes misc/reaction_engine_sanity.tsv. Not a gameplay reaction list.
void runReactionEngineSanityCheck();
