#pragma once

#include <algorithm>
#include <cmath>

// Grouped intrinsic / reference properties for SubstanceDefinition.
// World state (temperature, pressure, moisture amount, damage, velocity,
// liquid fraction, composition ratio, current phase) must NOT live here.

// Thermal units (internally consistent SI on a 2D grid with implied depth):
//   Temperature: Kelvin
//   Specific heat: J/(kg·K)
//   Conductivity: W/(m·K)  — physical; sandbox scale is ThermalConfig::conductivityScale
//   Energy: Joules
constexpr float AMBIENT_TEMPERATURE_K = 293.15f;
constexpr float WATER_DENSITY_KG_M3 = 1000.0f;
constexpr float AIR_DENSITY_KG_M3 = 1.204f;
// Dry-air molar mass that reproduces AIR_DENSITY_KG_M3 at 1 atm and ambient T
// via the ideal-gas law (used to scale other vapors, not a second air table).
constexpr float AIR_MOLAR_MASS_G_MOL = 28.97f;
constexpr float UNIVERSAL_GAS_R_J_MOL_K = 8.314462618f;
constexpr float MIN_THERMAL_CAPACITY = 1.0e-6f; // J/K
constexpr float MIN_THERMAL_MASS_KG = 1.0e-9f;
constexpr float MIN_SAFE_TEMPERATURE_K = 0.05f;
constexpr float MAX_SAFE_TEMPERATURE_K = 1.0e7f;

// Heat storage / conduction. Melting/boiling/latent here are table-authoring
// copies; PhaseProperties is the query path. No phase-change solver yet.
struct ThermalProperties {
    float specificHeat = 1000.0f;   // J/(kg·K)
    float conductivity = 0.026f;    // W/(m·K) physical
    float meltingPointK = 0.0f;
    float boilingPointK = 0.0f;
    float latentFusion = 0.0f;      // J/kg
    float latentVapor = 0.0f;       // J/kg
    float expansionCoeff = 0.0f;    // 1/K (linear solids, volumetric fluids)
    float softeningTempK = 0.0f;
    bool valid = false;
};

// Rigid / solid sandbox mechanics. Values are NOT SI; they match the existing
// MaterialDefinition table (densityRel 1 = water). This is SOLID-PHASE density
// for the rigid solver. Do not reuse it as liquid fill density or gas amount.
struct MechanicalProperties {
    bool valid = false;
    float densityRel = 0.0f;        // relative to water = 1.0
    float friction = 0.0f;
    float restitution = 0.0f;
    float hardness = 0.0f;
    float toughness = 0.0f;
    float brittleness = 0.0f;
    float tensileStrength = 0.0f;
    float compressiveStrength = 0.0f;
    float shearStrength = 0.0f;
    float fractureToughness = 0.0f;
};

// Free-surface liquid sandbox properties. Viscosity/tension are sandbox units,
// not SI Pa·s or N/m. density is relative to water at densityRefTempK.
// This is LIQUID-PHASE density for FluidEngine. Do not reuse it as rigid mass
// densityRel or as gas amount/EoS state.
struct FluidProperties {
    bool valid = false;
    float density = 1.0f;
    float viscosity = 0.0f;
    float surfaceTension = 0.0f;
    float viscRefTempK = AMBIENT_TEMPERATURE_K;
    // mu(T) = viscosity * exp(A * (1/T - 1/Tref)). A ≈ 1800 K roughly doubles
    // water viscosity from 20 °C to 0 °C. Not a lab curve; no negative viscosity.
    float viscArrheniusK = 1800.0f;
    float densityRefTempK = AMBIENT_TEMPERATURE_K;
    // Linear expansivity around lab T. Water's density maximum near 4 °C is NOT modeled.
    float densityExpansivity = 2.07e-4f;

    float viscosityAtTemperature(float temperatureK) const {
        float T = temperatureK;
        if (!(T > 1.0f) || !std::isfinite(T)) T = viscRefTempK > 1.0f ? viscRefTempK : AMBIENT_TEMPERATURE_K;
        float Tref = viscRefTempK > 1.0f ? viscRefTempK : AMBIENT_TEMPERATURE_K;
        float mu = viscosity * std::exp(viscArrheniusK * (1.0f / T - 1.0f / Tref));
        if (!std::isfinite(mu) || mu < 0.0f) return 0.0f;
        return mu;
    }

    // Data-layer only until the pressure solver can take variable density safely.
    // |ΔT| is clamped so a plasma-range T cannot collapse mass if this is called early.
    float densityAtTemperature(float temperatureK) const {
        float T = std::isfinite(temperatureK) ? temperatureK : densityRefTempK;
        float Tref = densityRefTempK > 1.0f ? densityRefTempK : AMBIENT_TEMPERATURE_K;
        float dT = std::clamp(T - Tref, -80.0f, 200.0f);
        float rho = density * (1.0f - densityExpansivity * dT);
        if (!std::isfinite(rho)) return density;
        return std::max(0.05f * std::max(0.0f, density), rho);
    }
};

// Intrinsic moisture / pore data. Local moisture amount is world state on the body.
struct PorousProperties {
    bool valid = false;
    float porosity = 0.0f;
    float moistureCapacity = 0.0f;
    float absorptionRate = 0.0f;
    float permeability = 0.0f;
    float dryingRate = 0.0f;
    float wetStrengthMultiplier = 1.0f;
    float wetFractureToughnessMultiplier = 1.0f;
    float wetFrictionMultiplier = 1.0f;
    float wetDensityContribution = 0.0f;
};

// Capability + transition metadata. Does not trigger phase changes.
// Current world phase is NOT stored here — see MatterPhase / MatterIdentity.
struct PhaseProperties {
    bool valid = false;
    bool solidCapable = false;
    bool liquidCapable = false;
    bool gasCapable = false;
    float meltingPointK = 0.0f;
    float boilingPointK = 0.0f;
    float latentHeatFusion = 0.0f;       // J/kg
    float latentHeatVaporization = 0.0f; // J/kg
    float referencePressurePa = 101325.0f;
};

// Future-use only. Not simulated. Prefer valid=false over invented precision.
struct ElectricalProperties {
    bool valid = false;
    float conductivity = 0.0f;          // S/m-ish placeholder, unused
    float relativePermittivity = 1.0f;
};

// Future-use only. Not simulated. Unknown is better than a fake formula.
struct ChemicalProperties {
    bool valid = false;
    float molarMass = 0.0f;             // g/mol; 0 = unknown
    bool flammable = false;
    bool oxidizer = false;
    float polarity = 0.0f;              // 0 = unknown / unset
    float corrosiveness = 0.0f;
};

// RGB seed for solids. Not a renderer type; physics must not include WorldVisual.
struct SubstanceVisualMetadata {
    bool valid = false;
    int colorR = 128;
    int colorG = 128;
    int colorB = 128;
};
