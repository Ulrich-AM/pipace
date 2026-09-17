#pragma once

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

// Local, conservative reaction execution. Empty ReactionRegistry is a no-op.
// Liquid-only for this milestone. Gas/solid participants skip the reaction
// without mutating the cell (GasEngine is Air + water vapor, not generic gas).

struct ReactionEngine {
    void simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
        ThermalEngine &thermal, float dt);
};

// Diagnostic-only. Uses synthetic inventories, not SubstanceRegistry / UI.
// Writes misc/reaction_engine_sanity.tsv. Not a gameplay reaction list.
void runReactionEngineSanityCheck();
