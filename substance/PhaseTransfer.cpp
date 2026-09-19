#include "substance/PhaseTransfer.h"

#include "fluid/DiagOutput.h"
#include "substance/SubstanceRegistry.h"
#include "world/PhaseChangeEngine.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace {

constexpr double kCpm = 4.0;
constexpr double kAbsTol = 1.0e-9;
constexpr double kRelTol = 1.0e-9;
constexpr double kFloatRelTol = 1.0e-5;

bool nearAbs(double a, double b, double tol = kAbsTol) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    return std::abs(a - b) <= tol;
}

bool nearRel(double a, double b, double rel = kRelTol, double absFloor = kAbsTol) {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    double scale = std::max({std::abs(a), std::abs(b), 1.0});
    return std::abs(a - b) <= std::max(absFloor, rel * scale);
}

std::string f8(double v) {
    std::ostringstream o;
    o << std::setprecision(12) << v;
    return o.str();
}

} // namespace

void runPhaseTransferDiagnostics() {
    std::ofstream out(miscFile("phase_transfer_diag.tsv"));
    out << std::setprecision(12);
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };

    double dx = sandboxCellLengthM(kCpm);
    double vol = sandboxCellVolumeM3(kCpm);
    double waterCellMass = liquidFillToMassKg(SUBSTANCE_WATER, 1.0, kCpm);
    double airCellMass = gasAmountToMassKg(SUBSTANCE_AIR, 1.0, kCpm);
    double vaporRho = referenceGasDensityKgM3(SUBSTANCE_WATER);
    double vaporCellMass = gasAmountToMassKg(SUBSTANCE_WATER, 1.0, kCpm);

    emit("unit_cell_edge_m", nearAbs(dx, 0.25), "dx=" + f8(dx));
    emit("unit_cell_volume_m3", nearAbs(vol, 0.015625), "V=" + f8(vol));
    emit("unit_full_water_cell_kg", nearAbs(waterCellMass, 15.625),
        "mass=" + f8(waterCellMass) + " fill=1 densityRel=1");
    emit("unit_air_amount_matches_thermal",
        nearRel(airCellMass, static_cast<double>(gasMassKg(1.0f, 4.0f)), 1.0e-6),
        "phase_transfer=" + f8(airCellMass)
            + " gasMassKg=" + f8(static_cast<double>(gasMassKg(1.0f, 4.0f))));
    emit("unit_water_vapor_density", vaporRho > 0.5 && vaporRho < 0.9,
        "rho=" + f8(vaporRho) + " kg/m3 at ambient T, reference P, molarMass");
    emit("unit_water_vapor_not_air_density",
        std::abs(vaporRho - static_cast<double>(AIR_DENSITY_KG_M3)) > 0.2,
        "vapor is not treated as air");

    double fill0 = 1.0;
    double massFromFill = liquidFillToMassKg(SUBSTANCE_WATER, fill0, kCpm);
    double fillBack = massKgToLiquidFill(SUBSTANCE_WATER, massFromFill, kCpm);
    emit("water_liquid_fill_roundtrip", nearRel(fillBack, fill0),
        "fill0=" + f8(fill0) + " fill1=" + f8(fillBack) + " mass=" + f8(massFromFill));

    double gasAmt0 = 1.0;
    double massFromGas = gasAmountToMassKg(SUBSTANCE_WATER, gasAmt0, kCpm);
    double gasBack = massKgToGasAmount(SUBSTANCE_WATER, massFromGas, kCpm);
    emit("water_gas_amount_roundtrip", nearRel(gasBack, gasAmt0),
        "amt0=" + f8(gasAmt0) + " amt1=" + f8(gasBack) + " mass=" + f8(massFromGas));

    PhaseTransferResult liqToGas = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Gas, fill0, kCpm);
    emit("liquid_mass_to_vapor_amount",
        liqToGas.success && nearRel(liqToGas.massTransferred, waterCellMass)
            && liqToGas.destinationAmountAdded > 1.0,
        "mass=" + f8(liqToGas.massTransferred)
            + " vapor_amount=" + f8(liqToGas.destinationAmountAdded)
            + " vapor_cell_mass=" + f8(vaporCellMass));

    PhaseTransferResult gasToLiq = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Gas, MatterPhase::Liquid,
        liqToGas.destinationAmountAdded, kCpm);
    emit("vapor_amount_to_liquid_mass",
        gasToLiq.success && nearRel(gasToLiq.massTransferred, waterCellMass)
            && nearRel(gasToLiq.destinationAmountAdded, fill0),
        "mass=" + f8(gasToLiq.massTransferred)
            + " fill=" + f8(gasToLiq.destinationAmountAdded));

    emit("liquid_gas_mass_error",
        nearRel(liqToGas.massTransferred, gasToLiq.massTransferred),
        "dmass=" + f8(liqToGas.massTransferred - gasToLiq.massTransferred));

    PhaseProperties const &waterPhase = phaseForSubstance(SUBSTANCE_WATER);
    double fusion = fusionEnergyJ(SUBSTANCE_WATER, waterCellMass);
    double vapor = vaporizationEnergyJ(SUBSTANCE_WATER, waterCellMass);
    emit("latent_fusion_matches_registry",
        nearRel(fusion, waterCellMass * static_cast<double>(waterPhase.latentHeatFusion)),
        "J=" + f8(fusion) + " L_f=" + f8(waterPhase.latentHeatFusion));
    emit("latent_vaporization_matches_registry",
        nearRel(vapor, waterCellMass * static_cast<double>(waterPhase.latentHeatVaporization)),
        "J=" + f8(vapor) + " L_v=" + f8(waterPhase.latentHeatVaporization));
    emit("latent_not_hardcoded_in_helper",
        nearRel(fusion, waterCellMass * static_cast<double>(waterPhase.latentHeatFusion))
            && std::abs(waterPhase.latentHeatFusion - 3.34e5f) < 1.0f
            && std::abs(waterPhase.latentHeatVaporization - 2.26e6f) < 1.0f,
        "helpers read PhaseProperties");
    emit("convert_liquid_gas_latent_sign",
        liqToGas.success && nearRel(liqToGas.energyTransferred, vapor) && liqToGas.energyTransferred > 0.0,
        "E=" + f8(liqToGas.energyTransferred));
    emit("convert_gas_liquid_latent_sign",
        gasToLiq.success && nearRel(gasToLiq.energyTransferred, -vapor) && gasToLiq.energyTransferred < 0.0,
        "E=" + f8(gasToLiq.energyTransferred));

    emit("water_liquid_to_gas_supported",
        canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Gas), "");
    emit("water_gas_to_liquid_supported",
        canTransition(SUBSTANCE_WATER, MatterPhase::Gas, MatterPhase::Liquid), "");
    emit("water_liquid_to_solid_supported",
        canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Solid),
        "live freeze uses this metadata");
    emit("wood_solid_to_liquid_unsupported",
        !canTransition(SUBSTANCE_WOOD, MatterPhase::Solid, MatterPhase::Liquid), "");
    emit("invalid_id_fails_safe",
        !canTransition(999, MatterPhase::Liquid, MatterPhase::Gas)
            && !canTransition(SUBSTANCE_NONE, MatterPhase::Liquid, MatterPhase::Gas), "");
    emit("same_phase_not_a_transition",
        !canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Liquid), "");
    emit("plasma_unsupported",
        !canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Plasma), "");
    emit("honey_liquid_to_gas_unsupported",
        !canTransition(SUBSTANCE_HONEY, MatterPhase::Liquid, MatterPhase::Gas), "");

    PhaseTransferResult iceAttempt = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Solid, 1.0, kCpm);
    double icePixelMass = solidFractionToMassKg(SUBSTANCE_WATER, 1.0, kCpm);
    emit("water_solid_convert_uses_ice_density",
        iceAttempt.success && hasMechanicalProperties(SUBSTANCE_WATER)
            && nearRel(iceAttempt.massTransferred, waterCellMass)
            && nearRel(iceAttempt.destinationAmountAdded, waterCellMass / icePixelMass)
            && icePixelMass < waterCellMass,
        "mass=" + f8(iceAttempt.massTransferred)
            + " solid_frac=" + f8(iceAttempt.destinationAmountAdded)
            + " ice_pixel_kg=" + f8(icePixelMass));
    emit("convert_liquid_solid_latent_sign",
        iceAttempt.success && iceAttempt.energyTransferred < 0.0
            && nearRel(iceAttempt.energyTransferred, -fusion),
        "E=" + f8(iceAttempt.energyTransferred));
    PhaseTransferResult iceBack = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Solid, MatterPhase::Liquid,
        iceAttempt.destinationAmountAdded, kCpm);
    emit("water_solid_liquid_roundtrip",
        iceBack.success && nearRel(iceBack.destinationAmountAdded, fill0)
            && nearRel(iceBack.massTransferred, waterCellMass),
        "fill=" + f8(iceBack.destinationAmountAdded));

    emit("no_ice_substance_id",
        SUBSTANCE_COUNT == 12 && substanceFromInternalName("ice") == SUBSTANCE_NONE, "");
    emit("no_steam_substance_id",
        substanceFromInternalName("steam") == SUBSTANCE_NONE, "");

    double startFill = 10.0;
    double startMass = liquidFillToMassKg(SUBSTANCE_WATER, startFill, kCpm);
    PhaseTransferResult tenCells = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Gas, startFill, kCpm);
    PhaseTransferResult tenBack = convertPhaseAmount(
        SUBSTANCE_WATER, MatterPhase::Gas, MatterPhase::Liquid,
        tenCells.destinationAmountAdded, kCpm);
    emit("closed_10_cells_mass",
        tenCells.success && tenBack.success
            && nearRel(startMass, tenBack.massTransferred)
            && nearRel(startFill, tenBack.destinationAmountAdded),
        "start_mass=" + f8(startMass)
            + " end_mass=" + f8(tenBack.massTransferred)
            + " end_fill=" + f8(tenBack.destinationAmountAdded));

    double cycleMass = startMass;
    double cycleAmount = 0.0;
    int const kCycles = 1000;
    for (int i = 0; i < kCycles; ++i) {
        cycleAmount = massKgToGasAmount(SUBSTANCE_WATER, cycleMass, kCpm);
        cycleMass = gasAmountToMassKg(SUBSTANCE_WATER, cycleAmount, kCpm);
    }
    double drift = cycleMass - startMass;
    emit("repeat_1000_liquid_gas_mass_drift",
        nearRel(cycleMass, startMass, 1.0e-12, 1.0e-12),
        "cycles=" + std::to_string(kCycles)
            + " start=" + f8(startMass)
            + " end=" + f8(cycleMass)
            + " drift=" + f8(drift));

    double cycleFill = startFill;
    for (int i = 0; i < kCycles; ++i) {
        double m = liquidFillToMassKg(SUBSTANCE_WATER, cycleFill, kCpm);
        cycleFill = massKgToLiquidFill(SUBSTANCE_WATER, m, kCpm);
    }
    emit("repeat_1000_liquid_fill_drift",
        nearRel(cycleFill, startFill, 1.0e-12, 1.0e-12),
        "start=" + f8(startFill) + " end=" + f8(cycleFill)
            + " drift=" + f8(cycleFill - startFill));

    float fFill = 1.0f;
    float fMass = static_cast<float>(liquidFillToMassKg(SUBSTANCE_WATER, fFill, kCpm));
    float fFill2 = static_cast<float>(massKgToLiquidFill(SUBSTANCE_WATER, fMass, kCpm));
    emit("float_liquid_fill_roundtrip",
        nearRel(static_cast<double>(fFill2), static_cast<double>(fFill), kFloatRelTol),
        "fill=" + f8(fFill2));

    float fAmt = static_cast<float>(liqToGas.destinationAmountAdded);
    float fMassG = static_cast<float>(gasAmountToMassKg(SUBSTANCE_WATER, fAmt, kCpm));
    float fAmt2 = static_cast<float>(massKgToGasAmount(SUBSTANCE_WATER, fMassG, kCpm));
    emit("float_vapor_amount_roundtrip",
        nearRel(static_cast<double>(fAmt2), static_cast<double>(fAmt), kFloatRelTol),
        "amt0=" + f8(fAmt) + " amt1=" + f8(fAmt2));

    emit("air_roundtrip_still_air_density",
        nearRel(gasAmountToMassKg(SUBSTANCE_AIR, 1.0, kCpm),
            static_cast<double>(AIR_DENSITY_KG_M3) * vol),
        "air helpers unchanged");
    emit("specific_gas_constant_water",
        specificGasConstantJPerKgK(SUBSTANCE_WATER) > 400.0
            && specificGasConstantJPerKgK(SUBSTANCE_WATER) < 500.0,
        "R_spec=" + f8(specificGasConstantJPerKgK(SUBSTANCE_WATER)));

    double P0 = static_cast<double>(waterPhase.referencePressurePa);
    double Tb = static_cast<double>(waterPhase.boilingPointK);
    double PsatTb = saturationVaporPressurePa(SUBSTANCE_WATER, Tb);
    double TsatP0 = saturationTemperatureK(SUBSTANCE_WATER, P0);
    emit("sat_pressure_at_tb_is_reference",
        nearRel(PsatTb, P0, 1.0e-5, 1.0),
        "Psat(Tb)=" + f8(PsatTb) + " P0=" + f8(P0));
    emit("sat_temperature_at_pref_is_tb",
        nearAbs(TsatP0, Tb, 0.05),
        "Tsat(P0)=" + f8(TsatP0) + " Tb=" + f8(Tb));
    emit("sat_pressure_rises_with_temperature",
        saturationVaporPressurePa(SUBSTANCE_WATER, Tb + 25.0)
            > saturationVaporPressurePa(SUBSTANCE_WATER, Tb),
        "Psat(Tb+25)=" + f8(saturationVaporPressurePa(SUBSTANCE_WATER, Tb + 25.0)));
    emit("sat_temperature_rises_with_pressure",
        saturationTemperatureK(SUBSTANCE_WATER, P0 * 2.0) > Tb + 10.0,
        "Tsat(2atm)=" + f8(saturationTemperatureK(SUBSTANCE_WATER, P0 * 2.0)));
    emit("low_pressure_lowers_boiling_temperature",
        saturationTemperatureK(SUBSTANCE_WATER, P0 * 0.1) < Tb - 20.0,
        "Tsat(0.1atm)=" + f8(saturationTemperatureK(SUBSTANCE_WATER, P0 * 0.1)));
    emit("sat_helpers_roundtrip_near_tb",
        nearRel(saturationVaporPressurePa(SUBSTANCE_WATER,
                saturationTemperatureK(SUBSTANCE_WATER, P0 * 1.5)), P0 * 1.5, 1.0e-3, 50.0),
        "P=" + f8(P0 * 1.5)
            + " Tsat=" + f8(saturationTemperatureK(SUBSTANCE_WATER, P0 * 1.5))
            + " Psat(Tsat)=" + f8(saturationVaporPressurePa(SUBSTANCE_WATER,
                saturationTemperatureK(SUBSTANCE_WATER, P0 * 1.5))));

    emit("live_lg_water_liquid_gas_eligible",
        supportsLiveLiquidGasTransition(SUBSTANCE_WATER)
            && canTransition(SUBSTANCE_WATER, MatterPhase::Liquid, MatterPhase::Gas)
            && canTransition(SUBSTANCE_WATER, MatterPhase::Gas, MatterPhase::Liquid), "");
    emit("live_lg_honey_not_eligible",
        !supportsLiveLiquidGasTransition(SUBSTANCE_HONEY), "");
    emit("live_lg_co2_not_eligible",
        !supportsLiveLiquidGasTransition(SUBSTANCE_CARBON_DIOXIDE), "");
    emit("live_lg_carbon_not_eligible",
        !supportsLiveLiquidGasTransition(SUBSTANCE_CARBON), "");
    emit("live_lg_none_not_eligible",
        !supportsLiveLiquidGasTransition(SUBSTANCE_NONE)
            && !supportsLiveLiquidGasTransition(static_cast<SubstanceId>(0xFFFF)), "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t'
        << passed << " passed, " << failed << " failed\n";
}
