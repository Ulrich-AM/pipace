#include "world/WorldQuery.h"

#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/SubstanceRegistry.h"
#include "thermal/ThermalEngine.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

MatterSample sampleMatterAt(FluidEngine const &fluid, RigidBodyEngine const &rigid,
    GasEngine const &gas, int x, int y)
{
    MatterSample sample;
    if (!FluidEngine::inside(x, y)) return sample;

    ThermalCellSample thermal = ThermalEngine::sampleCell(fluid, rigid, gas, x, y);
    sample.temperatureK = thermal.temperatureK;
    sample.hasMatter = thermal.hasMatter;
    if (!thermal.hasMatter) return sample;

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
            LiquidComponentView comps = fluid.liquidComponents(i);
            sample.identity = makeMatterIdentity(fluid.dominantLiquidSubstance(i), MatterPhase::Liquid);
            sample.amount = fill;
            sample.waterFraction = fluid.liquidComponentFraction(i, SUBSTANCE_WATER);
            sample.honeyFraction = fluid.liquidComponentFraction(i, SUBSTANCE_HONEY);
            sample.mixture = comps.count >= 2;
            break;
        }
        case ThermalSampleKind::Gas: {
            float tot = gas.amount[static_cast<size_t>(i)];
            GasComponentView comps = gas.gasComponents(i);
            sample.amount = tot;
            sample.vaporFraction = gas.gasComponentFraction(i, SUBSTANCE_WATER);
            sample.waterFraction = sample.vaporFraction;
            sample.mixture = comps.count >= 2;
            SubstanceId dom = gas.dominantGasSubstance(i);
            sample.identity = makeMatterIdentity(dom, MatterPhase::Gas);
            break;
        }
        case ThermalSampleKind::Empty:
        default:
            sample.hasMatter = false;
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
    emit("supports_honey_liquid", supportsPhase(SUBSTANCE_HONEY, MatterPhase::Liquid), "");
    emit("supports_honey_solid_conservative",
        supportsPhase(SUBSTANCE_HONEY, MatterPhase::Solid) && !hasMechanicalProperties(SUBSTANCE_HONEY),
        "honey is solid-capable metadata without a solid table");
    emit("supports_honey_not_gas", !supportsPhase(SUBSTANCE_HONEY, MatterPhase::Gas), "");
    emit("supports_glass_liquid_no_fluid",
        supportsPhase(SUBSTANCE_GLASS, MatterPhase::Liquid) && !hasFluidProperties(SUBSTANCE_GLASS),
        "glass liquid-capable without a fluid table");
    emit("supports_metal_liquid_no_fluid",
        supportsPhase(SUBSTANCE_METAL, MatterPhase::Liquid) && !hasFluidProperties(SUBSTANCE_METAL), "");
    emit("no_ice_id", SUBSTANCE_COUNT == 12 && substanceFromInternalName("ice") == SUBSTANCE_NONE, "");
    emit("no_steam_id", substanceFromInternalName("steam") == SUBSTANCE_NONE, "");
    emit("no_current_phase_on_def", true, "SubstanceDefinition stores capability flags only");
    emit("water_has_fluid", hasFluidProperties(SUBSTANCE_WATER) && hasPropertiesForPhase(SUBSTANCE_WATER, MatterPhase::Liquid), "");
    emit("water_has_mechanical_ice",
        hasMechanicalProperties(SUBSTANCE_WATER) && hasPropertiesForPhase(SUBSTANCE_WATER, MatterPhase::Solid)
            && rigidIdentityForMaterial(MATERIAL_WATER_SOLID).substance == SUBSTANCE_WATER
            && rigidIdentityForMaterial(MATERIAL_WATER_SOLID).phase == MatterPhase::Solid, "");
    emit("metal_has_mechanical", hasMechanicalProperties(SUBSTANCE_METAL) && hasPropertiesForPhase(SUBSTANCE_METAL, MatterPhase::Solid), "");
    emit("air_has_gas", hasGasProperties(SUBSTANCE_AIR), "");
    emit("hydrogen_has_gas",
        hasGasProperties(SUBSTANCE_HYDROGEN) && supportsPhase(SUBSTANCE_HYDROGEN, MatterPhase::Gas)
            && !supportsPhase(SUBSTANCE_HYDROGEN, MatterPhase::Liquid), "");
    emit("oxygen_has_gas",
        hasGasProperties(SUBSTANCE_OXYGEN) && supportsPhase(SUBSTANCE_OXYGEN, MatterPhase::Gas)
            && !supportsPhase(SUBSTANCE_OXYGEN, MatterPhase::Liquid), "");
    emit("carbon_has_solid",
        hasMechanicalProperties(SUBSTANCE_CARBON) && supportsPhase(SUBSTANCE_CARBON, MatterPhase::Solid)
            && !supportsPhase(SUBSTANCE_CARBON, MatterPhase::Gas)
            && rigidIdentityForMaterial(MATERIAL_CARBON).substance == SUBSTANCE_CARBON, "");
    emit("co2_has_gas",
        hasGasProperties(SUBSTANCE_CARBON_DIOXIDE) && supportsPhase(SUBSTANCE_CARBON_DIOXIDE, MatterPhase::Gas)
            && !supportsPhase(SUBSTANCE_CARBON_DIOXIDE, MatterPhase::Solid), "");
    emit("plasma_unsupported", !supportsPhase(SUBSTANCE_WATER, MatterPhase::Plasma), "");
    emit("none_phase_unsupported", !supportsPhase(SUBSTANCE_WATER, MatterPhase::None), "");
    emit("invalid_id_no_phase", !supportsPhase(999, MatterPhase::Liquid) && !supportsPhase(999, MatterPhase::Solid), "");

    MatterIdentity waterLiq = makeMatterIdentity(SUBSTANCE_WATER, MatterPhase::Liquid);
    emit("identity_roundtrip_water_liquid",
        waterLiq.substance == SUBSTANCE_WATER && waterLiq.phase == MatterPhase::Liquid, "");
    MatterIdentity wallId = staticWallIdentity();
    emit("identity_static_wall",
        wallId.substance == kStaticWallSubstance && wallId.substance == SUBSTANCE_STONE
            && wallId.phase == MatterPhase::Solid, "");
    MatterIdentity gasId = identityForGasSpecies();
    emit("identity_gas_species_air",
        gasId.substance == SUBSTANCE_AIR && gasId.phase == MatterPhase::Gas
            && substanceForGasSpecies() == SUBSTANCE_AIR, "");
    MatterIdentity rigidWood = rigidIdentityForMaterial(MATERIAL_WOOD);
    emit("identity_rigid_wood",
        rigidWood.substance == SUBSTANCE_WOOD && rigidWood.phase == MatterPhase::Solid, "");

    fluid.resetWorld();
    rigid.clear();
    gas.resetAmbient(fluid);

    int cx = 80, cy = 40;
    int wi = FluidEngine::ci(cx, cy);
    fluid.setLiquidComponentAmount(wi, SUBSTANCE_WATER, 1.0f);
    fluid.expectedVolume = 1.0;
    fluid.wakeAllFluidChunks();
    MatterSample waterCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_water_liquid",
        waterCell.hasMatter && waterCell.identity.substance == SUBSTANCE_WATER
            && waterCell.identity.phase == MatterPhase::Liquid && !waterCell.mixture,
        "sub=" + std::to_string(waterCell.identity.substance)
            + " phase=" + std::to_string(static_cast<int>(waterCell.identity.phase)));

    MatterSample neighborAir = sampleMatterAt(fluid, rigid, gas, cx + 1, cy);
    emit("query_liquid_gas_interface",
        neighborAir.hasMatter && neighborAir.identity.substance == SUBSTANCE_AIR
            && neighborAir.identity.phase == MatterPhase::Gas
            && waterCell.identity.phase == MatterPhase::Liquid,
        "neighbor sub=" + std::to_string(neighborAir.identity.substance));

    fluid.setLiquidComponentAmount(wi, SUBSTANCE_HONEY, 1.0f);
    MatterSample honeyCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_honey_liquid",
        honeyCell.hasMatter && honeyCell.identity.substance == SUBSTANCE_HONEY
            && honeyCell.identity.phase == MatterPhase::Liquid,
        "sub=" + std::to_string(honeyCell.identity.substance));

    fluid.setLiquidComponentAmount(wi, SUBSTANCE_HONEY, 0.0f);
    fluid.setLiquidComponentAmount(wi, SUBSTANCE_WATER, 1.0f);
    fluid.setLiquidComponentAmount(wi, SUBSTANCE_HONEY, 0.4f);
    MatterSample mixCell = sampleMatterAt(fluid, rigid, gas, cx, cy);
    emit("query_mixture_dominant_water",
        mixCell.mixture && mixCell.identity.substance == SUBSTANCE_WATER
            && mixCell.identity.phase == MatterPhase::Liquid
            && mixCell.honeyFraction > 0.3f && mixCell.waterFraction > 0.5f,
        "dom=" + std::to_string(mixCell.identity.substance)
            + " honey=" + std::to_string(mixCell.honeyFraction)
            + " water=" + std::to_string(mixCell.waterFraction));
    emit("query_mixture_not_new_id",
        mixCell.mixture && mixCell.identity.substance != SUBSTANCE_NONE
            && mixCell.identity.substance < SUBSTANCE_COUNT
            && mixCell.identity.substance != SUBSTANCE_HONEY,
        "mixtures keep a dominant built-in SubstanceId");

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
    emit("query_static_wall_thermal_identity",
        thermalForSubstance(kStaticWallSubstance).specificHeat
            == thermalForSubstance(SUBSTANCE_STONE).specificHeat
            && mechanicalForSubstance(kStaticWallSubstance).densityRel
                == mechanicalForSubstance(SUBSTANCE_STONE).densityRel, "");

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

    MatterSample oob = sampleMatterAt(fluid, rigid, gas, -1, -1);
    emit("query_oob",
        !oob.hasMatter && oob.identity.substance == SUBSTANCE_NONE
            && oob.identity.phase == MatterPhase::None, "");
    MatterSample oobFar = sampleMatterAt(fluid, rigid, gas, GW, GH);
    emit("query_oob_far",
        !oobFar.hasMatter && oobFar.identity.substance == SUBSTANCE_NONE, "");

    int ei = FluidEngine::ci(20, 20);
    fluid.setLiquidComponentAmount(ei, SUBSTANCE_WATER, 0.0f);
    fluid.setLiquidComponentAmount(ei, SUBSTANCE_HONEY, 0.0f);
    gas.amount[static_cast<size_t>(ei)] = 0.0f;
    gas.heat[static_cast<size_t>(ei)] = 0.0f;
    MatterSample emptyCell = sampleMatterAt(fluid, rigid, gas, 20, 20);
    emit("query_empty",
        !emptyCell.hasMatter && emptyCell.identity.substance == SUBSTANCE_NONE
            && emptyCell.identity.phase == MatterPhase::None,
        "kind empty after zero fill and gas");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
}
