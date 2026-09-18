#pragma once

#include "chemistry/ReactionTypes.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "substance/SubstanceTypes.h"

// Chemistry talks to the world as SubstanceId + MatterPhase.
//
// Liquid uses FluidEngine composition slots. Gas uses GasEngine generic
// composition. Solid reactants for SolidGasSurface live on RigidBodyEngine
// source pixels, not in the gas/liquid cell. Plasma execution is unsupported.

bool reactionPhaseStorageSupported(SubstanceId id, MatterPhase requiredPhase);
float reactionQueryMatter(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId id, MatterPhase phase);
bool reactionCatalystPresent(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId catalyst);
// Cheap occupancy gate. True only if every used reactant has storage in the
// required phase (Gas / Liquid / either if None). No temperature, moles, or
// inventory work. Does not inspect products or catalysts. SolidGasSurface
// checks gas reactants only.
bool reactionCellHasRequiredReactants(FluidEngine const &fluid, GasEngine const &gas, int index,
    ReactionDefinition const &def);
// Gas-only reactant presence. Ignores solid participants. No T / moles / heat.
bool reactionCellHasRequiredGasReactants(FluidEngine const &fluid, GasEngine const &gas, int index,
    ReactionDefinition const &def);
bool reactionReadLiquidOccupancy(FluidEngine const &fluid, int index, LiquidComponentView &out);
bool reactionCommitLiquidOccupancy(FluidEngine &fluid, int index, LiquidComponentView const &view);
bool reactionReadGasOccupancy(GasEngine const &gas, int index, GasComponentView &out);
bool reactionCommitGasOccupancy(GasEngine &gas, int index, GasComponentView const &view);
