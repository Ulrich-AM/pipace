#pragma once

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

// Local, conservative reaction execution. Scans liquid and gas-occupied cells.
// Physical reactions convert storage -> mass -> moles using chemical.molarMass.
// Extent 1 is one mole of the written reaction. Synthetic inventories keep a
// diagnostic unit path. Liquid and gas stay separate; solid/plasma still skip.

struct ReactionEngine {
    void simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
        ThermalEngine &thermal, float dt);
};

// Diagnostic-only. Uses synthetic inventories, not SubstanceRegistry / UI.
// Writes misc/reaction_engine_sanity.tsv. Not a gameplay reaction list.
void runReactionEngineSanityCheck();
