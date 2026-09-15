#include "substance/SubstanceRegistry.h"

#include "fluid/DiagOutput.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>

namespace {

// Table authoring still writes transition metadata on ThermalProperties.
// PhaseProperties is the query path (phaseForSubstance). Copied once at init.
void copyPhaseTransitionFromThermal(SubstanceDefinition &s) {
    s.phase.meltingPointK = s.thermal.meltingPointK;
    s.phase.boilingPointK = s.thermal.boilingPointK;
    s.phase.latentHeatFusion = s.thermal.latentFusion;
    s.phase.latentHeatVaporization = s.thermal.latentVapor;
    s.phase.referencePressurePa = 101325.0f;
}

SubstanceDefinition makeNone() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_NONE;

    // Thermal
    s.thermal.valid = false;
    s.thermal.specificHeat = 1.0f;
    s.thermal.conductivity = 0.0f;
    s.thermal.meltingPointK = 0.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 0.0f;
    s.thermal.softeningTempK = 0.0f;

    // Visual
    s.visual.colorR = 0;
    s.visual.colorG = 0;
    s.visual.colorB = 0;

    return s;
}

SubstanceDefinition makeWater() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_WATER;
    s.internalName = "water";
    s.displayName = "Water";
    s.displayNameKey = "ins_mat_water";

    // Classification / composition
    s.classification = SubstanceClass::PureSubstance;
    s.compositionKind = CompositionKind::PureChemical;
    s.formulaHint = "H2O";

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 4184.0f;
    s.thermal.conductivity = 0.598f;
    s.thermal.meltingPointK = 273.15f;
    s.thermal.boilingPointK = 373.15f;
    s.thermal.latentFusion = 3.34e5f;
    s.thermal.latentVapor = 2.26e6f;
    s.thermal.expansionCoeff = 2.07e-4f;
    s.thermal.softeningTempK = 0.0f;
    s.thermal.gasSpecificHeat = 2080.0f; // water vapor; liquid Cp stays on specificHeat

    // Phase capability. Solid/gas flags are metadata only — no ice body and no
    // steam SubstanceId. Water has no MechanicalProperties; do not borrow stone.
    // Gas-phase mass uses chemical.molarMass (see substance/PhaseTransfer.h).
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = true;
    s.phase.gasCapable = true;
    copyPhaseTransitionFromThermal(s);

    // Fluid
    s.fluid.valid = true;
    // Base liquid properties
    s.fluid.density = 1.0f;
    s.fluid.viscosity = 0.006f;
    s.fluid.surfaceTension = 0.055f;
    // Temperature response
    s.fluid.viscRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.viscArrheniusK = 1800.0f;
    s.fluid.densityRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.densityExpansivity = 2.07e-4f;

    // Chemical
    s.chemical.valid = true;
    s.chemical.molarMass = 18.015f;
    s.chemical.flammable = false;
    s.chemical.oxidizer = false;
    s.chemical.polarity = 1.0f;
    s.chemical.corrosiveness = 0.0f;

    return s;
}

SubstanceDefinition makeHoney() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_HONEY;
    s.internalName = "honey";
    s.displayName = "Honey";
    s.displayNameKey = "ins_mat_honey";

    // Classification / composition
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    s.formulaHint = nullptr;

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 2200.0f;
    s.thermal.conductivity = 0.50f;
    s.thermal.meltingPointK = 0.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 2.07e-4f;
    s.thermal.softeningTempK = 0.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = true;
    s.phase.gasCapable = false;
    copyPhaseTransitionFromThermal(s);

    // Fluid
    s.fluid.valid = true;
    // Base liquid properties
    s.fluid.density = 1.42f;
    s.fluid.viscosity = 0.35f;
    s.fluid.surfaceTension = 0.080f;
    // Temperature response
    s.fluid.viscRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.viscArrheniusK = 2500.0f;
    s.fluid.densityRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.densityExpansivity = 2.07e-4f;

    // Chemical — mixture, no single molecular formula
    s.chemical.valid = false;

    return s;
}

SubstanceDefinition makeWood() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_WOOD;
    s.internalName = "wood";
    s.displayName = "wood";
    s.displayNameKey = "wood";

    // Classification / composition
    s.classification = SubstanceClass::Biological;
    s.compositionKind = CompositionKind::Composite;

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 1700.0f;
    s.thermal.conductivity = 0.12f;
    s.thermal.meltingPointK = 0.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 5.0e-6f;
    s.thermal.softeningTempK = 450.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = false;
    s.phase.gasCapable = false;
    copyPhaseTransitionFromThermal(s);

    // Mechanical
    s.mechanical.valid = true;
    // Mass / contact
    s.mechanical.densityRel = 0.45f;
    s.mechanical.friction = 0.48f;
    s.mechanical.restitution = 0.08f;
    // Material character
    s.mechanical.hardness = 0.42f;
    s.mechanical.toughness = 0.52f;
    s.mechanical.brittleness = 0.34f;
    // Structural strength
    s.mechanical.tensileStrength = 0.55f;
    s.mechanical.compressiveStrength = 0.70f;
    s.mechanical.shearStrength = 0.48f;
    s.mechanical.fractureToughness = 0.50f;

    // Porous / moisture
    s.porous.valid = true;
    // Capacity / structure
    s.porous.porosity = 0.55f;
    s.porous.moistureCapacity = 0.62f;
    // Transport
    s.porous.absorptionRate = 0.26f;
    s.porous.permeability = 0.68f;
    s.porous.dryingRate = 0.0010f;
    // Wet modifiers
    s.porous.wetStrengthMultiplier = 0.62f;
    s.porous.wetFractureToughnessMultiplier = 0.70f;
    s.porous.wetFrictionMultiplier = 1.30f;
    s.porous.wetDensityContribution = 0.00f;

    // Chemical
    s.chemical.valid = true;
    s.chemical.molarMass = 0.0f;
    s.chemical.flammable = true;

    // Visual
    s.visual.valid = true;
    s.visual.colorR = 158;
    s.visual.colorG = 112;
    s.visual.colorB = 62;

    return s;
}

SubstanceDefinition makeStone() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_STONE;
    s.internalName = "stone";
    s.displayName = "stone";
    s.displayNameKey = "stone";

    // Classification / composition
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 880.0f;
    s.thermal.conductivity = 1.70f;
    s.thermal.meltingPointK = 1473.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 8.0e-6f;
    s.thermal.softeningTempK = 0.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = false;
    s.phase.gasCapable = false;
    copyPhaseTransitionFromThermal(s);

    // Mechanical
    s.mechanical.valid = true;
    // Mass / contact
    s.mechanical.densityRel = 2.20f;
    s.mechanical.friction = 0.55f;
    s.mechanical.restitution = 0.06f;
    // Material character
    s.mechanical.hardness = 1.35f;
    s.mechanical.toughness = 1.05f;
    s.mechanical.brittleness = 0.28f;
    // Structural strength
    s.mechanical.tensileStrength = 1.35f;
    s.mechanical.compressiveStrength = 1.80f;
    s.mechanical.shearStrength = 1.10f;
    s.mechanical.fractureToughness = 0.88f;

    // Porous / moisture
    s.porous.valid = true;
    // Capacity / structure
    s.porous.porosity = 0.08f;
    s.porous.moistureCapacity = 0.12f;
    // Transport
    s.porous.absorptionRate = 0.040f;
    s.porous.permeability = 0.08f;
    s.porous.dryingRate = 0.0006f;
    // Wet modifiers
    s.porous.wetStrengthMultiplier = 0.92f;
    s.porous.wetFractureToughnessMultiplier = 0.95f;
    s.porous.wetFrictionMultiplier = 1.05f;
    s.porous.wetDensityContribution = 0.00f;

    // Chemical
    s.chemical.valid = false;

    // Visual
    s.visual.valid = true;
    s.visual.colorR = 118;
    s.visual.colorG = 122;
    s.visual.colorB = 128;

    return s;
}

SubstanceDefinition makeGlass() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_GLASS;
    s.internalName = "glass";
    s.displayName = "glass";
    s.displayNameKey = "glass";

    // Classification / composition
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 840.0f;
    s.thermal.conductivity = 1.00f;
    s.thermal.meltingPointK = 1700.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 9.0e-6f;
    s.thermal.softeningTempK = 800.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = true;
    s.phase.gasCapable = false;
    copyPhaseTransitionFromThermal(s);

    // Mechanical
    s.mechanical.valid = true;
    // Mass / contact
    s.mechanical.densityRel = 2.50f;
    s.mechanical.friction = 0.22f;
    s.mechanical.restitution = 0.05f;
    // Material character
    s.mechanical.hardness = 1.75f;
    s.mechanical.toughness = 0.16f;
    s.mechanical.brittleness = 0.92f;
    // Structural strength
    s.mechanical.tensileStrength = 0.32f;
    s.mechanical.compressiveStrength = 1.40f;
    s.mechanical.shearStrength = 0.28f;
    s.mechanical.fractureToughness = 0.12f;

    // Chemical
    s.chemical.valid = false;

    // Visual
    s.visual.valid = true;
    s.visual.colorR = 168;
    s.visual.colorG = 204;
    s.visual.colorB = 214;

    return s;
}

SubstanceDefinition makeMetal() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_METAL;
    s.internalName = "metal";
    s.displayName = "metal";
    s.displayNameKey = "metal";

    // Classification / composition
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;

    // Thermal — generic metal, steel-ish table, not a specific element
    s.thermal.valid = true;
    s.thermal.specificHeat = 490.0f;
    s.thermal.conductivity = 50.0f;
    s.thermal.meltingPointK = 1811.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 1.2e-5f;
    s.thermal.softeningTempK = 1000.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = true;
    s.phase.liquidCapable = true;
    s.phase.gasCapable = false;
    copyPhaseTransitionFromThermal(s);

    // Mechanical
    s.mechanical.valid = true;
    // Mass / contact
    s.mechanical.densityRel = 7.80f;
    s.mechanical.friction = 0.38f;
    s.mechanical.restitution = 0.12f;
    // Material character
    s.mechanical.hardness = 2.10f;
    s.mechanical.toughness = 1.85f;
    s.mechanical.brittleness = 0.12f;
    // Structural strength
    s.mechanical.tensileStrength = 2.40f;
    s.mechanical.compressiveStrength = 2.80f;
    s.mechanical.shearStrength = 1.90f;
    s.mechanical.fractureToughness = 1.65f;

    // Chemical
    s.chemical.valid = false;

    // Visual
    s.visual.valid = true;
    s.visual.colorR = 148;
    s.visual.colorG = 152;
    s.visual.colorB = 158;

    return s;
}

SubstanceDefinition makeAir() {
    SubstanceDefinition s;

    // Identity
    s.id = SUBSTANCE_AIR;
    s.internalName = "air";
    s.displayName = "Air";
    s.displayNameKey = "ins_mat_air";

    // Classification / composition
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;

    // Thermal
    s.thermal.valid = true;
    s.thermal.specificHeat = 1005.0f;
    s.thermal.conductivity = 0.026f;
    s.thermal.meltingPointK = 0.0f;
    s.thermal.boilingPointK = 0.0f;
    s.thermal.latentFusion = 0.0f;
    s.thermal.latentVapor = 0.0f;
    s.thermal.expansionCoeff = 1.0f / AMBIENT_TEMPERATURE_K;
    s.thermal.softeningTempK = 0.0f;

    // Phase capability
    s.phase.valid = true;
    s.phase.solidCapable = false;
    s.phase.liquidCapable = false;
    s.phase.gasCapable = true;
    copyPhaseTransitionFromThermal(s);

    // Chemical
    s.chemical.valid = false;

    return s;
}

} // namespace

SubstanceDefinition const *builtinSubstanceTable() {
    static SubstanceDefinition const table[SUBSTANCE_COUNT] = {
        makeNone(),
        makeWater(),
        makeHoney(),
        makeWood(),
        makeStone(),
        makeGlass(),
        makeMetal(),
        makeAir(),
    };
    return table;
}

SubstanceId substanceFromInternalName(char const *name) {
    if (!name || !name[0]) return SUBSTANCE_NONE;
    SubstanceDefinition const *table = builtinSubstanceTable();
    for (SubstanceId id = 0; id < SUBSTANCE_COUNT; ++id) {
        char const *internal = table[id].internalName;
        if (internal && std::strcmp(internal, name) == 0) return id;
    }
    return SUBSTANCE_NONE;
}

namespace {

char const *classificationName(SubstanceClass c) {
    switch (c) {
        case SubstanceClass::PureSubstance: return "pure";
        case SubstanceClass::Mixture: return "mixture";
        case SubstanceClass::Composite: return "composite";
        case SubstanceClass::Biological: return "biological";
        case SubstanceClass::Unknown: return "unknown";
    }
    return "unknown";
}

char const *compositionName(CompositionKind k) {
    switch (k) {
        case CompositionKind::PureChemical: return "pure_chemical";
        case CompositionKind::Mixture: return "mixture";
        case CompositionKind::Composite: return "composite";
        case CompositionKind::Unknown: return "unknown";
    }
    return "unknown";
}

} // namespace

void runSubstanceRegistryDiagnostics() {
    std::ofstream out(miscFile("substance_registry_diag.tsv"));
    out << std::setprecision(8);
    SubstanceDefinition const *table = builtinSubstanceTable();

    out << "id\tinternal_name\tclassification\tcomposition_kind"
        << "\tsolid_capable\tliquid_capable\tgas_capable"
        << "\tmechanical_valid\tfluid_valid\tthermal_valid"
        << "\tporous_valid\tchemical_valid\telectrical_valid\n";
    for (SubstanceId id = 0; id < SUBSTANCE_COUNT; ++id) {
        SubstanceDefinition const &s = table[id];
        out << s.id << '\t' << s.internalName
            << '\t' << classificationName(s.classification)
            << '\t' << compositionName(s.compositionKind)
            << '\t' << (s.phase.solidCapable ? 1 : 0)
            << '\t' << (s.phase.liquidCapable ? 1 : 0)
            << '\t' << (s.phase.gasCapable ? 1 : 0)
            << '\t' << (s.mechanical.valid ? 1 : 0)
            << '\t' << (s.fluid.valid ? 1 : 0)
            << '\t' << (s.thermal.valid ? 1 : 0)
            << '\t' << (s.porous.valid ? 1 : 0)
            << '\t' << (s.chemical.valid ? 1 : 0)
            << '\t' << (s.electrical.valid ? 1 : 0) << '\n';
    }

    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed; else ++failed;
    };
    auto near = [](float a, float b, float tol = 1.0e-4f) {
        return std::fabs(a - b) <= tol;
    };

    bool idsUnique = true;
    bool namesUnique = true;
    bool indexMatches = true;
    for (SubstanceId id = 0; id < SUBSTANCE_COUNT; ++id) {
        if (table[id].id != id) indexMatches = false;
        for (SubstanceId other = static_cast<SubstanceId>(id + 1); other < SUBSTANCE_COUNT; ++other) {
            if (table[id].id == table[other].id) idsUnique = false;
            if (table[id].internalName && table[other].internalName
                && std::strcmp(table[id].internalName, table[other].internalName) == 0)
                namesUnique = false;
        }
    }
    emit("ids_unique", idsUnique, "");
    emit("names_unique", namesUnique, "");
    emit("table_index_matches_id", indexMatches, "");
    emit("builtin_count", SUBSTANCE_COUNT == 8, std::to_string(SUBSTANCE_COUNT));
    emit("none_is_slot_zero", table[SUBSTANCE_NONE].id == SUBSTANCE_NONE
        && std::strcmp(table[SUBSTANCE_NONE].internalName, "none") == 0, "");
    emit("none_no_fake_physics",
        !table[SUBSTANCE_NONE].mechanical.valid
            && !table[SUBSTANCE_NONE].fluid.valid
            && !table[SUBSTANCE_NONE].thermal.valid
            && !table[SUBSTANCE_NONE].porous.valid
            && !table[SUBSTANCE_NONE].phase.valid, "");

    SubstanceDefinition const &bad = substanceDef(999);
    emit("invalid_lookup_fallback_none", bad.id == SUBSTANCE_NONE, "");
    emit("invalid_lookup_no_crash", true, "substanceDef(999) returned");
    emit("invalid_lookup_no_fake_valid",
        !bad.mechanical.valid && !bad.fluid.valid && !bad.thermal.valid
            && !bad.porous.valid && !bad.chemical.valid && !bad.electrical.valid, "");
    SubstanceDefinition const &badWide = substanceDef(65535);
    emit("invalid_lookup_wide", badWide.id == SUBSTANCE_NONE && !badWide.fluid.valid, "");

    SubstanceDefinition const &water = substanceDef(SUBSTANCE_WATER);
    emit("water_fluid_valid", water.fluid.valid, "");
    emit("water_thermal_valid", water.thermal.valid, "");
    emit("water_mechanical_invalid", !water.mechanical.valid, "");
    emit("water_porous_invalid", !water.porous.valid, "");
    emit("water_electrical_inactive", !water.electrical.valid, "");
    emit("water_density", near(water.fluid.density, 1.0f), std::to_string(water.fluid.density));
    emit("water_viscosity", near(water.fluid.viscosity, 0.006f), std::to_string(water.fluid.viscosity));
    emit("water_surface_tension", near(water.fluid.surfaceTension, 0.055f), std::to_string(water.fluid.surfaceTension));
    emit("water_specific_heat", near(water.thermal.specificHeat, 4184.0f, 0.01f), std::to_string(water.thermal.specificHeat));
    emit("water_conductivity", near(water.thermal.conductivity, 0.598f), std::to_string(water.thermal.conductivity));
    emit("water_melting_phase", near(phaseForSubstance(SUBSTANCE_WATER).meltingPointK, 273.15f), "");
    emit("water_boiling_phase", near(phaseForSubstance(SUBSTANCE_WATER).boilingPointK, 373.15f), "");
    emit("sandbox_reference_is_water_fluid",
        sandboxReferenceLiquid().density == water.fluid.density
            && sandboxReferenceLiquid().viscosity == water.fluid.viscosity, "");

    SubstanceDefinition const &wood = substanceDef(SUBSTANCE_WOOD);
    emit("wood_mechanical_valid", wood.mechanical.valid, "");
    emit("wood_porous_valid", wood.porous.valid, "");
    emit("wood_fluid_invalid", !wood.fluid.valid, "");
    emit("wood_density_rel", near(wood.mechanical.densityRel, 0.45f), std::to_string(wood.mechanical.densityRel));
    emit("wood_porosity", near(wood.porous.porosity, 0.55f), std::to_string(wood.porous.porosity));
    emit("wood_adapter_density",
        near(materialDef(MATERIAL_WOOD).density, wood.mechanical.densityRel), "");

    SubstanceDefinition const &metal = substanceDef(SUBSTANCE_METAL);
    emit("metal_mechanical_valid", metal.mechanical.valid, "");
    emit("metal_porous_invalid", !metal.porous.valid, "");
    emit("metal_fluid_invalid", !metal.fluid.valid, "");
    emit("metal_density_rel", near(metal.mechanical.densityRel, 7.80f), std::to_string(metal.mechanical.densityRel));
    emit("metal_adapter_density",
        near(materialDef(MATERIAL_METAL).density, metal.mechanical.densityRel), "");

    SubstanceDefinition const &air = substanceDef(SUBSTANCE_AIR);
    emit("air_gas_capable", air.phase.gasCapable && !air.phase.solidCapable && !air.phase.liquidCapable, "");
    emit("air_fluid_invalid", !air.fluid.valid, "");
    emit("air_mechanical_invalid", !air.mechanical.valid, "");
    emit("air_thermal_valid", air.thermal.valid, "");
    emit("air_specific_heat", near(air.thermal.specificHeat, 1005.0f, 0.01f), std::to_string(air.thermal.specificHeat));

    SubstanceDefinition const &honey = substanceDef(SUBSTANCE_HONEY);
    emit("honey_fluid_valid", honey.fluid.valid, "");
    emit("honey_mechanical_invalid", !honey.mechanical.valid, "");
    emit("honey_density", near(honey.fluid.density, 1.42f), std::to_string(honey.fluid.density));
    emit("honey_chemical_unknown", !honey.chemical.valid, "mixture; no fake formula");

    SubstanceDefinition const &glass = substanceDef(SUBSTANCE_GLASS);
    emit("glass_porous_invalid", !glass.porous.valid, "");
    emit("glass_mechanical_valid", glass.mechanical.valid, "");

    emit("static_wall_is_stone", kStaticWallSubstance == SUBSTANCE_STONE, "");
    emit("static_wall_thermal_matches_stone",
        thermalForSubstance(kStaticWallSubstance).specificHeat
            == thermalForSubstance(SUBSTANCE_STONE).specificHeat, "");
    emit("all_electrical_inactive",
        !water.electrical.valid && !wood.electrical.valid && !metal.electrical.valid
            && !air.electrical.valid && !honey.electrical.valid, "");
    emit("no_ice_or_steam_slots",
        substanceFromInternalName("ice") == SUBSTANCE_NONE
            && substanceFromInternalName("steam") == SUBSTANCE_NONE, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, " << failed << " failed\n";
}
