#pragma once

#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "substance/SubstanceTypes.h"

// Chemistry talks to the world as SubstanceId + MatterPhase.
//
// Current GasEngine stores Air and water vapor only. It does NOT provide
// generic SubstanceId gas composition. Arbitrary gas-phase consume/produce is
// therefore unsupported here (skip, do not fabricate). Solid execution is
// also unsupported. Liquid uses FluidEngine's generic composition slots.

bool reactionPhaseStorageSupported(SubstanceId id, MatterPhase requiredPhase);
float reactionQueryMatter(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId id, MatterPhase phase);
bool reactionCatalystPresent(FluidEngine const &fluid, int index, SubstanceId catalyst);
bool reactionReadLiquidOccupancy(FluidEngine const &fluid, int index, LiquidComponentView &out);
bool reactionCommitLiquidOccupancy(FluidEngine &fluid, int index, LiquidComponentView const &view);
