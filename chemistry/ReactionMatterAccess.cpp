#include "chemistry/ReactionMatterAccess.h"

#include "fluid/FluidTypes.h"

bool reactionPhaseStorageSupported(SubstanceId id, MatterPhase requiredPhase) {
    if (requiredPhase == MatterPhase::Gas || requiredPhase == MatterPhase::Solid
        || requiredPhase == MatterPhase::Plasma)
        return false;
    if (requiredPhase != MatterPhase::Liquid && requiredPhase != MatterPhase::None)
        return false;
    return validLiquidComponentId(id);
}

float reactionQueryMatter(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId id, MatterPhase phase)
{
    (void)gas;
    if (index < 0 || index >= GW * GH) return 0.0f;
    if (!reactionPhaseStorageSupported(id, phase)) return 0.0f;
    return fluid.liquidComponentAmount(index, id);
}

bool reactionCatalystPresent(FluidEngine const &fluid, int index, SubstanceId catalyst) {
    if (catalyst == SUBSTANCE_NONE) return true;
    if (index < 0 || index >= GW * GH) return false;
    if (fluid.liquidComponentAmount(index, catalyst) > kMinLiquidComponent) return true;
    return false;
}

bool reactionReadLiquidOccupancy(FluidEngine const &fluid, int index, LiquidComponentView &out) {
    out = {};
    if (index < 0 || index >= GW * GH) return false;
    out = fluid.liquidComponents(index);
    return true;
}

bool reactionCommitLiquidOccupancy(FluidEngine &fluid, int index, LiquidComponentView const &view) {
    return fluid.tryCommitLiquidOccupancy(index, view);
}
