#include "chemistry/ReactionMatterAccess.h"

#include "fluid/FluidTypes.h"

bool reactionPhaseStorageSupported(SubstanceId id, MatterPhase requiredPhase) {
    if (requiredPhase == MatterPhase::Solid || requiredPhase == MatterPhase::Plasma)
        return false;
    if (requiredPhase == MatterPhase::Gas)
        return validGasComponentId(id);
    if (requiredPhase == MatterPhase::Liquid)
        return validLiquidComponentId(id);
    if (requiredPhase == MatterPhase::None)
        return validLiquidComponentId(id) || validGasComponentId(id);
    return false;
}

float reactionQueryMatter(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId id, MatterPhase phase)
{
    if (index < 0 || index >= GW * GH) return 0.0f;
    if (!reactionPhaseStorageSupported(id, phase)) return 0.0f;
    if (phase == MatterPhase::Liquid) return fluid.liquidComponentAmount(index, id);
    if (phase == MatterPhase::Gas) return gas.gasComponentAmount(index, id);
    // None: do not sum liquid and gas into one fake occupancy.
    float liquid = fluid.liquidComponentAmount(index, id);
    if (liquid > kMinLiquidComponent) return liquid;
    return gas.gasComponentAmount(index, id);
}

bool reactionCatalystPresent(FluidEngine const &fluid, GasEngine const &gas, int index,
    SubstanceId catalyst)
{
    if (catalyst == SUBSTANCE_NONE) return true;
    if (index < 0 || index >= GW * GH) return false;
    if (fluid.liquidComponentAmount(index, catalyst) > kMinLiquidComponent) return true;
    if (gas.gasComponentAmount(index, catalyst) > kMinGasComponent) return true;
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

bool reactionReadGasOccupancy(GasEngine const &gas, int index, GasComponentView &out) {
    out = {};
    if (index < 0 || index >= GW * GH) return false;
    out = gas.gasComponents(index);
    return true;
}

bool reactionCommitGasOccupancy(GasEngine &gas, int index, GasComponentView const &view) {
    return gas.tryCommitGasOccupancy(index, view);
}
