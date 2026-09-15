#include "substance/SubstanceRegistry.h"

#include <cstring>

namespace {

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

    // Phase capability
    s.canExistAsSolid = true;
    s.canExistAsLiquid = true;
    s.canExistAsGas = true;
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
    s.canExistAsSolid = true;
    s.canExistAsLiquid = true;
    s.canExistAsGas = false;
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
    s.canExistAsSolid = true;
    s.canExistAsLiquid = false;
    s.canExistAsGas = false;
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
    s.canExistAsSolid = true;
    s.canExistAsLiquid = false;
    s.canExistAsGas = false;
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
    s.canExistAsSolid = true;
    s.canExistAsLiquid = true;
    s.canExistAsGas = false;
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
    s.canExistAsSolid = true;
    s.canExistAsLiquid = true;
    s.canExistAsGas = false;
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
    s.canExistAsSolid = false;
    s.canExistAsLiquid = false;
    s.canExistAsGas = true;
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
