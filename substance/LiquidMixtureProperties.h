#pragma once

#include "fluid/FluidTypes.h"
#include "substance/SubstanceTypes.h"

// Effective liquid properties from a cell's composition.
//
// Mixing laws (approximations, not chemistry / EOS / effective-medium theory):
//   density, conductivity, surface tension — volume-fraction (phi) linear mix
//   specific heat — mass-weighted, m_i ∝ amount_i * density_i
//   viscosity — Arrhenius-style log mix: ln(mu) = Σ phi_i ln(mu_i(T))
//
// Fractions use component amounts, not cell fill. A 0.2-full cell with the same
// composition ratios has the same intrinsic properties as a full cell.
//
// Empty or invalid composition does not invent Water identity. The returned
// numbers fall back to sandboxReferenceLiquid() / Water thermal tables so
// solvers stay finite. Callers must not write that fallback into composition.

struct LiquidMixtureProperties {
    float density = 1.0f;
    float specificHeat = 4184.0f;
    float conductivity = 0.60f;
    float viscosity = 0.0f;
    float surfaceTension = 0.0f;
};

// Numerical safety floor for log-viscosity. Not a physical cutoff.
constexpr float kMinMixtureViscosity = 1.0e-8f;

LiquidMixtureProperties referenceLiquidMixture(float temperatureK);
LiquidMixtureProperties evaluateLiquidMixture(
    LiquidComponentView const &composition,
    float temperatureK);
