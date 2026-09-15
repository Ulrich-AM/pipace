#pragma once

#include "substance/SubstanceTypes.h"

struct FluidEngine;
struct RigidBodyEngine;
struct GasEngine;

// Read-only occupancy query. Occupancy order matches thermal sampling:
// rigid pixel, static wall, substantial liquid, then ambient gas.
// Out of bounds and empty cells return hasMatter=false (identity NONE/None).
// Mixtures report the dominant SubstanceId plus water/honey fractions; there
// is no mixture SubstanceId. Not for hot simulation loops.
MatterSample sampleMatterAt(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas, int x, int y);

void runSubstancePhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas);
