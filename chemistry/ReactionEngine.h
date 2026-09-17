#pragma once

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;
struct ThermalEngine;

// Local, conservative reaction execution. Empty ReactionRegistry is a no-op.
// Liquid and gas participants use separate inventories (FluidEngine fill vs
// GasEngine amount). Solid/plasma still skip. No player-facing gas reactions yet.

struct ReactionEngine {
    void simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
        ThermalEngine &thermal, float dt);
};

// Diagnostic-only. Uses synthetic inventories, not SubstanceRegistry / UI.
// Writes misc/reaction_engine_sanity.tsv. Not a gameplay reaction list.
void runReactionEngineSanityCheck();
