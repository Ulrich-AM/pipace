#include "substance/SubstanceRegistry.h"

#include <cstring>

namespace {

ThermalProperties makeThermal(float cp, float k, float melt, float boil,
    float Lf, float Lv, float alpha, float soft, bool valid)
{
    ThermalProperties t;
    t.specificHeat = cp;
    t.conductivity = k;
    t.meltingPointK = melt;
    t.boilingPointK = boil;
    t.latentFusion = Lf;
    t.latentVapor = Lv;
    t.expansionCoeff = alpha;
    t.softeningTempK = soft;
    t.valid = valid;
    return t;
}

void setPhaseFromThermal(SubstanceDefinition &s, bool solidC, bool liquidC, bool gasC) {
    s.canExistAsSolid = solidC;
    s.canExistAsLiquid = liquidC;
    s.canExistAsGas = gasC;
    s.phase.valid = true;
    s.phase.solidCapable = solidC;
    s.phase.liquidCapable = liquidC;
    s.phase.gasCapable = gasC;
    s.phase.meltingPointK = s.thermal.meltingPointK;
    s.phase.boilingPointK = s.thermal.boilingPointK;
    s.phase.latentHeatFusion = s.thermal.latentFusion;
    s.phase.latentHeatVaporization = s.thermal.latentVapor;
    s.phase.referencePressurePa = 101325.0f;
}

void setVisual(SubstanceDefinition &s, int r, int g, int b) {
    s.visual.valid = true;
    s.visual.colorR = r;
    s.visual.colorG = g;
    s.visual.colorB = b;
}

void setMechanical(SubstanceDefinition &s,
    float densityRel, float friction, float restitution,
    float hardness, float toughness, float brittleness,
    float tensile, float compressive, float shear, float fracture)
{
    s.mechanical.valid = true;
    s.mechanical.densityRel = densityRel;
    s.mechanical.friction = friction;
    s.mechanical.restitution = restitution;
    s.mechanical.hardness = hardness;
    s.mechanical.toughness = toughness;
    s.mechanical.brittleness = brittleness;
    s.mechanical.tensileStrength = tensile;
    s.mechanical.compressiveStrength = compressive;
    s.mechanical.shearStrength = shear;
    s.mechanical.fractureToughness = fracture;
}

void setPorous(SubstanceDefinition &s,
    bool valid, float porosity, float permeability, float capacity,
    float absorption, float drying,
    float wetStr, float wetFric, float wetFrac, float wetDens)
{
    s.porous.valid = valid;
    s.porous.porosity = porosity;
    s.porous.permeability = permeability;
    s.porous.moistureCapacity = capacity;
    s.porous.absorptionRate = absorption;
    s.porous.dryingRate = drying;
    s.porous.wetStrengthMultiplier = wetStr;
    s.porous.wetFrictionMultiplier = wetFric;
    s.porous.wetFractureToughnessMultiplier = wetFrac;
    s.porous.wetDensityContribution = wetDens;
}

SubstanceDefinition makeNone() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_NONE;
    s.thermal = makeThermal(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, false);
    s.visual.colorR = 0;
    s.visual.colorG = 0;
    s.visual.colorB = 0;
    return s;
}

SubstanceDefinition makeWater() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_WATER;
    s.internalName = "water";
    s.displayName = "Water";
    s.displayNameKey = "ins_mat_water";
    s.classification = SubstanceClass::PureSubstance;
    s.compositionKind = CompositionKind::PureChemical;
    s.formulaHint = "H2O";
    s.thermal = makeThermal(4184.0f, 0.598f, 273.15f, 373.15f, 3.34e5f, 2.26e6f, 2.07e-4f, 0.0f, true);
    setPhaseFromThermal(s, true, true, true);
    s.fluid.valid = true;
    s.fluid.density = 1.0f;
    s.fluid.viscosity = 0.006f;
    s.fluid.surfaceTension = 0.055f;
    s.fluid.viscRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.viscArrheniusK = 1800.0f;
    s.fluid.densityRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.densityExpansivity = 2.07e-4f;
    s.fluid.thermal = s.thermal;
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
    s.id = SUBSTANCE_HONEY;
    s.internalName = "honey";
    s.displayName = "Honey";
    s.displayNameKey = "ins_mat_honey";
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    s.formulaHint = nullptr;
    s.thermal = makeThermal(2200.0f, 0.50f, 0.0f, 0.0f, 0.0f, 0.0f, 2.07e-4f, 0.0f, true);
    setPhaseFromThermal(s, true, true, false);
    s.fluid.valid = true;
    s.fluid.density = 1.42f;
    s.fluid.viscosity = 0.35f;
    s.fluid.surfaceTension = 0.080f;
    s.fluid.viscRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.viscArrheniusK = 2500.0f;
    s.fluid.densityRefTempK = AMBIENT_TEMPERATURE_K;
    s.fluid.densityExpansivity = 2.07e-4f;
    s.fluid.thermal = s.thermal;
    // Mixture: no fake molar mass / formula.
    s.chemical.valid = false;
    return s;
}

SubstanceDefinition makeWood() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_WOOD;
    s.internalName = "wood";
    s.displayName = "wood";
    s.displayNameKey = "wood";
    s.classification = SubstanceClass::Biological;
    s.compositionKind = CompositionKind::Composite;
    s.thermal = makeThermal(1700.0f, 0.12f, 0.0f, 0.0f, 0.0f, 0.0f, 5.0e-6f, 450.0f, true);
    setPhaseFromThermal(s, true, false, false);
    setMechanical(s, 0.45f, 0.48f, 0.08f, 0.42f, 0.52f, 0.34f, 0.55f, 0.70f, 0.48f, 0.50f);
    setPorous(s, true, 0.55f, 0.68f, 0.62f, 0.26f, 0.0010f, 0.62f, 1.30f, 0.70f, 0.00f);
    setVisual(s, 158, 112, 62);
    s.chemical.valid = true;
    s.chemical.molarMass = 0.0f;
    s.chemical.flammable = true;
    return s;
}

SubstanceDefinition makeStone() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_STONE;
    s.internalName = "stone";
    s.displayName = "stone";
    s.displayNameKey = "stone";
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    s.thermal = makeThermal(880.0f, 1.70f, 1473.0f, 0.0f, 0.0f, 0.0f, 8.0e-6f, 0.0f, true);
    setPhaseFromThermal(s, true, false, false);
    setMechanical(s, 2.20f, 0.55f, 0.06f, 1.35f, 1.05f, 0.28f, 1.35f, 1.80f, 1.10f, 0.88f);
    setPorous(s, true, 0.08f, 0.08f, 0.12f, 0.040f, 0.0006f, 0.92f, 1.05f, 0.95f, 0.00f);
    setVisual(s, 118, 122, 128);
    s.chemical.valid = false;
    return s;
}

SubstanceDefinition makeGlass() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_GLASS;
    s.internalName = "glass";
    s.displayName = "glass";
    s.displayNameKey = "glass";
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    s.thermal = makeThermal(840.0f, 1.00f, 1700.0f, 0.0f, 0.0f, 0.0f, 9.0e-6f, 800.0f, true);
    setPhaseFromThermal(s, true, true, false);
    setMechanical(s, 2.50f, 0.22f, 0.05f, 1.75f, 0.16f, 0.92f, 0.32f, 1.40f, 0.28f, 0.12f);
    setPorous(s, false, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f, 1.00f, 1.00f, 0.00f);
    setVisual(s, 168, 204, 214);
    s.chemical.valid = false;
    return s;
}

SubstanceDefinition makeMetal() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_METAL;
    s.internalName = "metal";
    s.displayName = "metal";
    s.displayNameKey = "metal";
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    // Generic metal, not a specific element. Thermal numbers keep the existing steel-ish table.
    s.thermal = makeThermal(490.0f, 50.0f, 1811.0f, 0.0f, 0.0f, 0.0f, 1.2e-5f, 1000.0f, true);
    setPhaseFromThermal(s, true, true, false);
    setMechanical(s, 7.80f, 0.38f, 0.12f, 2.10f, 1.85f, 0.12f, 2.40f, 2.80f, 1.90f, 1.65f);
    setPorous(s, false, 0.00f, 0.00f, 0.00f, 0.00f, 0.00f, 1.00f, 1.00f, 1.00f, 0.00f);
    setVisual(s, 148, 152, 158);
    s.chemical.valid = false;
    return s;
}

SubstanceDefinition makeAir() {
    SubstanceDefinition s;
    s.id = SUBSTANCE_AIR;
    s.internalName = "air";
    s.displayName = "Air";
    s.displayNameKey = "ins_mat_air";
    s.classification = SubstanceClass::Mixture;
    s.compositionKind = CompositionKind::Mixture;
    s.thermal = makeThermal(1005.0f, 0.026f, 0.0f, 0.0f, 0.0f, 0.0f,
        1.0f / AMBIENT_TEMPERATURE_K, 0.0f, true);
    setPhaseFromThermal(s, false, false, true);
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
