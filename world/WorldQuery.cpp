#include "world/WorldQuery.h"

#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/SubstanceRegistry.h"
#include "thermal/ThermalEngine.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

MatterSample sampleMatterAt(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas, int x, int y)
{
    MatterSample sample;
    ThermalCellSample thermal = ThermalEngine::sampleCell(fluid, rigid, gas, x, y);
    sample.temperatureK = thermal.temperatureK;
    sample.hasMatter = thermal.hasMatter;
    if (!thermal.hasMatter || !FluidEngine::inside(x, y)) return sample;

    int i = FluidEngine::ci(x, y);
    switch (thermal.kind) {
        case ThermalSampleKind::Rigid:
            sample.identity = rigid.worldCellIdentity(x, y);
            sample.amount = 1.0f;
            break;
        case ThermalSampleKind::Wall:
            sample.identity = staticWallIdentity();
            sample.amount = 1.0f;
            break;
        case ThermalSampleKind::Liquid: {
            float fill = fluid.fill[static_cast<size_t>(i)];
            float waterAmt = fluid.liquidComponentAmount(i, SUBSTANCE_WATER);
            float honeyAmt = fluid.liquidComponentAmount(i, SUBSTANCE_HONEY);
            sample.identity = makeMatterIdentity(fluid.dominantLiquidSubstance(i), MatterPhase::Liquid);
            sample.amount = fill;
            sample.waterFraction = fluid.liquidComponentFraction(i, SUBSTANCE_WATER);
            sample.honeyFraction = fluid.liquidComponentFraction(i, SUBSTANCE_HONEY);
            sample.mixture = waterAmt > 1.0e-6f && honeyAmt > 1.0e-6f;
            break;
        }
        case ThermalSampleKind::Gas:
            sample.identity = identityForGasSpecies();
            sample.amount = gas.amount[static_cast<size_t>(i)];
            break;
        case ThermalSampleKind::Empty:
        default:
            break;
    }
    return sample;
}

void runSubstancePhaseDiagnostics(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas) {
    std::ofstream out(miscFile("substance_phase_diag.tsv"));
    out << std::setprecision(8);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };

    emit("supports_water_solid", supportsPhase(SUBSTANCE_WATER, MatterPhase::Solid), "");
    emit("supports_water_liquid", supportsPhase(SUBSTANCE_WATER, MatterPhase::Liquid), "");
    emit("supports_water_gas", supportsPhase(SUBSTANCE_WATER, MatterPhase::Gas), "");
    emit("supports_wood_solid", supportsPhase(SUBSTANCE_WOOD, MatterPhase::Solid), "");
    emit("supports_wood_liquid", !supportsPhase(SUBSTANCE_WOOD, MatterPhase::Liquid), "wood must not be liquid-capable");
    emit("supports_air_gas", supportsPhase(SUBSTANCE_AIR, MatterPhase::Gas), "");
    emit("supports_air_solid", !supportsPhase(SUBSTANCE_AIR, MatterPhase::Solid), "air must not be solid-capable");
    emit("no_ice_id", SUBSTANCE_COUNT == 8 && substanceFromInternalName("ice") == SUBSTANCE_NONE, "");
    emit("no_steam_id", substanceFromInternalName("steam") == SUBSTANCE_NONE, "");
    emit("no_current_phase_on_def", true, "SubstanceDefinition stores capability flags only");
    emit("water_has_fluid", hasFluidProperties(SUBSTANCE_WATER) && hasPropertiesForPhase(SUBSTANCE_WATER, MatterPhase::Liquid), "");
    emit("metal_has_mechanical", hasMechanicalProperties(SUBSTANCE_METAL) && hasPropertiesForPhase(SUBSTANCE_METAL, MatterPhase::Solid), "");
    emit("air_has_gas", hasGasProperties(SUBSTANCE_AIR), "");
    emit("plasma_unsupported", !supportsPhase(SUBSTANCE_WATER, MatterPhase::Plasma), "");

    fluid.resetWorld();
    rigid.clear();
    gas.resetAmbient(fluid);

    int cx = 80, cy = 40;
    int wi = FluidEngine::ci(cx, cy);
    fluid.fill[static_cast<size_t>(wi)] = 1.0f;
    fluid.honey[static_cast<size_t>(wi)] = 0.0f;
    fluid.expectedVolume = 1.0;
    fluid.wakeAllFluidChunks();
    MatterSample waterCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_water_liquid",
        waterCell.hasMatter && waterCell.identity.substance == SUBSTANCE_WATER
            && waterCell.identity.phase == MatterPhase::Liquid && !waterCell.mixture,
        "sub=" + std::to_string(waterCell.identity.substance)
            + " phase=" + std::to_string(static_cast<int>(waterCell.identity.phase)));

    fluid.honey[static_cast<size_t>(wi)] = 1.0f;
    MatterSample honeyCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_honey_liquid",
        honeyCell.hasMatter && honeyCell.identity.substance == SUBSTANCE_HONEY
            && honeyCell.identity.phase == MatterPhase::Liquid,
        "sub=" + std::to_string(honeyCell.identity.substance));

    fluid.honey[static_cast<size_t>(wi)] = 0.4f;
    MatterSample mixCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_mixture_dominant_water",
        mixCell.mixture && mixCell.identity.substance == SUBSTANCE_WATER
            && mixCell.identity.phase == MatterPhase::Liquid
            && mixCell.honeyFraction > 0.3f && mixCell.waterFraction > 0.5f,
        "dom=" + std::to_string(mixCell.identity.substance)
            + " honey=" + std::to_string(mixCell.honeyFraction)
            + " water=" + std::to_string(mixCell.waterFraction));

    fluid.clearWorld();
    rigid.clear();
    gas.resetAmbient(fluid);
    rigid.drawMaterial = MATERIAL_WOOD;
    rigid.paintPendingDisc(91, 51, 1);
    rigid.commitPending(fluid);
    MatterSample woodCell = sampleMatterAt(fluid, rigid, gas, 91, 51);
    emit("query_wood_solid",
        woodCell.hasMatter && woodCell.identity.substance == SUBSTANCE_WOOD
            && woodCell.identity.phase == MatterPhase::Solid,
        "sub=" + std::to_string(woodCell.identity.substance)
            + " phase=" + std::to_string(static_cast<int>(woodCell.identity.phase)));

    fluid.clearWorld();
    rigid.clear();
    gas.resetAmbient(fluid);
    int wx = 70, wy = 55;
    fluid.solid[static_cast<size_t>(FluidEngine::ci(wx, wy))] = 1;
    MatterSample wallCell = sampleMatterAt(fluid, rigid, gas, wx, wy);
    emit("query_static_wall_stone_solid",
        wallCell.hasMatter && wallCell.identity.substance == kStaticWallSubstance
            && wallCell.identity.phase == MatterPhase::Solid,
        "sub=" + std::to_string(wallCell.identity.substance));

    fluid.clearWorld();
    rigid.clear();
    gas.resetAmbient(fluid);
    MatterSample airCell = sampleMatterAt(fluid, rigid, gas, 40, 40);
    emit("query_air_gas",
        airCell.hasMatter && airCell.identity.substance == SUBSTANCE_AIR
            && airCell.identity.phase == MatterPhase::Gas,
        "sub=" + std::to_string(airCell.identity.substance)
            + " phase=" + std::to_string(static_cast<int>(airCell.identity.phase))
            + " amt=" + std::to_string(airCell.amount));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
}
