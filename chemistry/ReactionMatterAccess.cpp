#include "chemistry/ReactionMatterAccess.h"
#include "chemistry/ReactionRegistry.h"

#include "fluid/FluidTypes.h"
#include "substance/SubstanceRegistry.h"

bool reactionPhaseStorageSupported(SubstanceId id, MatterPhase requiredPhase) {
    if (requiredPhase == MatterPhase::Plasma)
        return false;
    if (requiredPhase == MatterPhase::Solid) {
        if (!validSubstance(id) || !supportsPhase(id, MatterPhase::Solid)) return false;
        if (!hasMechanicalProperties(id)) return false;
        return rigidMaterialForSubstance(id) != MATERIAL_EMPTY;
    }
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

bool reactionCellHasRequiredReactants(FluidEngine const &fluid, GasEngine const &gas, int index,
    ReactionDefinition const &def)
{
    if (!validReaction(def.id)) return false;
    if (index < 0 || index >= GW * GH) return false;
    if (def.topology == ReactionTopology::SolidGasSurface)
        return reactionCellHasRequiredGasReactants(fluid, gas, index, def);
    bool any = false;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        any = true;
        bool ok = false;
        if (p.requiredPhase == MatterPhase::Solid || p.requiredPhase == MatterPhase::Plasma)
            return false;
        if (p.requiredPhase == MatterPhase::Gas || p.requiredPhase == MatterPhase::None) {
            if (gas.gasComponentAmount(index, p.substance) > kMinGasComponent) ok = true;
        }
        if (!ok && (p.requiredPhase == MatterPhase::Liquid || p.requiredPhase == MatterPhase::None)) {
            if (fluid.liquidComponentAmount(index, p.substance) > kMinLiquidComponent) ok = true;
        }
        if (!ok) return false;
    }
    return any;
}

bool reactionCellHasRequiredGasReactants(FluidEngine const &fluid, GasEngine const &gas, int index,
    ReactionDefinition const &def)
{
    (void)fluid;
    if (index < 0 || index >= GW * GH) return false;
    bool anyGas = false;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        if (p.requiredPhase == MatterPhase::Solid || p.requiredPhase == MatterPhase::Plasma)
            continue;
        if (p.requiredPhase != MatterPhase::Gas && p.requiredPhase != MatterPhase::None)
            return false;
        anyGas = true;
        if (!(gas.gasComponentAmount(index, p.substance) > kMinGasComponent)) return false;
    }
    return anyGas;
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
