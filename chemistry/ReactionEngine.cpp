#include "chemistry/ReactionEngine.h"

#include "chemistry/ReactionMatterAccess.h"
#include "chemistry/ReactionRegistry.h"
#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "fluid/FluidTypes.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "thermal/ThermalEngine.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

namespace {

constexpr float kExtentEps = 1.0e-8f;
constexpr float kDefaultMaxExtentPerSecond = 1.0f;

struct Inventory {
    LiquidComponent items[kMaxLiquidComponents]{};
    int count = 0;
    float fill = 0.0f;
    float heat = 0.0f;
    float maxFill = 1.0f;
    bool allowUnregisteredIds = false;
    bool gas = false;
};

bool idAllowed(Inventory const &inv, SubstanceId id) {
    if (inv.allowUnregisteredIds) return id != SUBSTANCE_NONE;
    if (inv.gas) return validGasComponentId(id);
    return validLiquidComponentId(id);
}

float slotMin(Inventory const &inv) {
    return inv.gas ? kMinGasComponent : kMinLiquidComponent;
}

int findSlot(Inventory const &inv, SubstanceId id) {
    for (int n = 0; n < inv.count; ++n)
        if (inv.items[n].id == id) return n;
    return -1;
}

float amountOf(Inventory const &inv, SubstanceId id) {
    int n = findSlot(inv, id);
    return n >= 0 ? inv.items[n].amount : 0.0f;
}

Inventory *pickInv(ReactionParticipant const &p, Inventory &liquid, Inventory &gas) {
    if (p.requiredPhase == MatterPhase::Gas) return &gas;
    if (p.requiredPhase == MatterPhase::Liquid) return &liquid;
    if (amountOf(liquid, p.substance) > kExtentEps) return &liquid;
    if (amountOf(gas, p.substance) > kExtentEps) return &gas;
    if (liquid.allowUnregisteredIds || validLiquidComponentId(p.substance)) return &liquid;
    return &gas;
}

Inventory const *pickInvConst(ReactionParticipant const &p, Inventory const &liquid,
    Inventory const &gas)
{
    return pickInv(p, const_cast<Inventory &>(liquid), const_cast<Inventory &>(gas));
}

void recountFill(Inventory &inv) {
    float s = 0.0f;
    int w = 0;
    for (int n = 0; n < inv.count; ++n) {
        if (inv.items[n].amount > slotMin(inv) && inv.items[n].id != SUBSTANCE_NONE) {
            if (w != n) inv.items[w] = inv.items[n];
            s += inv.items[w].amount;
            ++w;
        }
    }
    for (int n = w; n < inv.count; ++n) inv.items[n] = {};
    inv.count = w;
    inv.fill = s;
}

bool consume(Inventory &inv, SubstanceId id, float amount) {
    if (!(amount > 0.0f)) return true;
    int n = findSlot(inv, id);
    if (n < 0 || inv.items[n].amount + 1.0e-7f < amount) return false;
    inv.items[n].amount -= amount;
    recountFill(inv);
    return true;
}

bool produce(Inventory &inv, SubstanceId id, float amount) {
    if (!(amount > 0.0f)) return true;
    if (!idAllowed(inv, id)) return false;
    if (inv.fill + amount > inv.maxFill + 1.0e-5f) return false;
    int n = findSlot(inv, id);
    if (n >= 0) {
        inv.items[n].amount += amount;
        recountFill(inv);
        return true;
    }
    if (inv.count >= kMaxLiquidComponents) return false;
    inv.items[inv.count++] = {id, amount};
    recountFill(inv);
    return true;
}

bool participantStorageOk(ReactionParticipant const &p, bool synthetic) {
    if (!reactionParticipantUsed(p)) return true;
    if (synthetic) {
        if (p.requiredPhase == MatterPhase::Solid || p.requiredPhase == MatterPhase::Plasma)
            return false;
        return p.substance != SUBSTANCE_NONE;
    }
    return reactionPhaseStorageSupported(p.substance, p.requiredPhase);
}

bool reactionStorageOk(ReactionDefinition const &def, bool synthetic) {
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n)
        if (reactionParticipantUsed(def.reactants[n]) && !participantStorageOk(def.reactants[n], synthetic))
            return false;
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n)
        if (reactionParticipantUsed(def.products[n]) && !participantStorageOk(def.products[n], synthetic))
            return false;
    return def.reactantCount > 0 && def.productCount > 0;
}

float limitingExtent(ReactionDefinition const &def, Inventory const &liquid, Inventory const &gas) {
    float extent = 1.0e30f;
    bool any = false;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        any = true;
        Inventory const *inv = pickInvConst(p, liquid, gas);
        float avail = amountOf(*inv, p.substance);
        if (!(p.coefficient > 0.0f) || !(avail > kExtentEps)) return 0.0f;
        extent = std::min(extent, avail / p.coefficient);
    }
    if (!any || !std::isfinite(extent) || extent <= kExtentEps) return 0.0f;
    return extent;
}

float limitingExtent(ReactionDefinition const &def, Inventory const &inv) {
    Inventory other;
    other.gas = !inv.gas;
    other.allowUnregisteredIds = inv.allowUnregisteredIds;
    other.maxFill = inv.gas ? 1.0f : 1.0e6f;
    return inv.gas ? limitingExtent(def, other, inv) : limitingExtent(def, inv, other);
}

float rateLimit(ReactionDefinition const &def, float dt, float maxExtent) {
    float rate = def.maxExtentPerSecond;
    if (!(rate > 0.0f) || !std::isfinite(rate)) rate = kDefaultMaxExtentPerSecond;
    if (!(dt > 0.0f) || !std::isfinite(dt)) return 0.0f;
    return std::min(maxExtent, rate * dt);
}

float heatLimit(ReactionDefinition const &def, float heat, float extent) {
    if (!(extent > kExtentEps)) return 0.0f;
    float dH = def.energyChangeJPerExtent;
    if (!(dH > 0.0f)) return extent; // exo or neutral: no heat budget
    if (!(heat > 0.0f) || !std::isfinite(heat)) return 0.0f;
    return std::min(extent, heat / dH);
}

bool applyExtent(ReactionDefinition const &def, Inventory &liquid, Inventory &gas, float extent) {
    if (!(extent > kExtentEps) || !std::isfinite(extent)) return false;
    Inventory trialL = liquid;
    Inventory trialG = gas;
    bool usesLiquid = false;
    bool usesGas = false;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        Inventory *inv = pickInv(p, trialL, trialG);
        if (inv->gas) usesGas = true;
        else usesLiquid = true;
        if (!consume(*inv, p.substance, p.coefficient * extent)) return false;
    }
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.products[n];
        if (!reactionParticipantUsed(p)) continue;
        Inventory *inv = pickInv(p, trialL, trialG);
        if (inv->gas) usesGas = true;
        else usesLiquid = true;
        if (!produce(*inv, p.substance, p.coefficient * extent)) return false;
    }
    float q = reactionHeatReleasedJ(def, extent);
    Inventory *heatInv = (usesLiquid && trialL.fill > slotMin(trialL)) ? &trialL : &trialG;
    if (!usesLiquid && usesGas) heatInv = &trialG;
    if (usesLiquid && !usesGas) heatInv = &trialL;
    float heat1 = heatInv->heat + q;
    if (!std::isfinite(heat1) || heat1 < 0.0f) return false;
    heatInv->heat = heat1;
    liquid = trialL;
    gas = trialG;
    return true;
}

bool applyExtent(ReactionDefinition const &def, Inventory &inv, float extent) {
    Inventory other;
    other.gas = !inv.gas;
    other.allowUnregisteredIds = inv.allowUnregisteredIds;
    other.maxFill = inv.gas ? 1.0f : 1.0e6f;
    bool ok = inv.gas ? applyExtent(def, other, inv, extent) : applyExtent(def, inv, other, extent);
    return ok;
}

float planExtent(ReactionDefinition const &def, Inventory const &liquid, Inventory const &gas, float dt) {
    if (!reactionStorageOk(def, liquid.allowUnregisteredIds || gas.allowUnregisteredIds)) return 0.0f;
    float extent = limitingExtent(def, liquid, gas);
    extent = rateLimit(def, dt, extent);
    float heat = std::max(liquid.heat, 0.0f) + std::max(gas.heat, 0.0f);
    extent = heatLimit(def, heat, extent);
    if (!(extent > kExtentEps)) return 0.0f;
    Inventory trialL = liquid;
    Inventory trialG = gas;
    for (int attempt = 0; attempt < 8; ++attempt) {
        trialL = liquid;
        trialG = gas;
        if (applyExtent(def, trialL, trialG, extent)) return extent;
        extent *= 0.5f;
        if (!(extent > kExtentEps)) return 0.0f;
    }
    return 0.0f;
}

Inventory fromView(LiquidComponentView const &view, float fill, float heat, bool synthetic) {
    Inventory inv;
    inv.allowUnregisteredIds = synthetic;
    inv.maxFill = 1.0f;
    inv.heat = heat;
    inv.fill = fill;
    inv.count = 0;
    inv.gas = false;
    for (int n = 0; n < view.count && n < kMaxLiquidComponents; ++n) {
        if (view.items[n].amount > kMinLiquidComponent)
            inv.items[inv.count++] = view.items[n];
    }
    recountFill(inv);
    return inv;
}

Inventory fromGasView(GasComponentView const &view, float amt, float heat, bool synthetic) {
    Inventory inv;
    inv.allowUnregisteredIds = synthetic;
    inv.maxFill = 1.0e6f;
    inv.heat = heat;
    inv.fill = amt;
    inv.count = 0;
    inv.gas = true;
    for (int n = 0; n < view.count && n < kMaxGasComponents; ++n) {
        if (view.items[n].amount > kMinGasComponent)
            inv.items[inv.count++] = {view.items[n].id, view.items[n].amount};
    }
    recountFill(inv);
    return inv;
}

LiquidComponentView toView(Inventory const &inv) {
    LiquidComponentView view;
    view.count = inv.count;
    for (int n = 0; n < inv.count; ++n) view.items[n] = inv.items[n];
    return view;
}

GasComponentView toGasView(Inventory const &inv) {
    GasComponentView view;
    view.count = std::min(inv.count, kMaxGasComponents);
    for (int n = 0; n < view.count; ++n)
        view.items[n] = {inv.items[n].id, inv.items[n].amount};
    return view;
}

float planExtent(ReactionDefinition const &def, Inventory const &inv, float dt) {
    Inventory other;
    other.gas = !inv.gas;
    other.allowUnregisteredIds = inv.allowUnregisteredIds;
    other.maxFill = inv.gas ? 1.0f : 1.0e6f;
    return inv.gas ? planExtent(def, other, inv, dt) : planExtent(def, inv, other, dt);
}

bool tryReactCell(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int index, ReactionDefinition const &def, float dt)
{
    if (!validReaction(def.id)) return false;
    if (!reactionStorageOk(def, false)) return false;
    if (index < 0 || index >= GW * GH) return false;
    size_t i = static_cast<size_t>(index);
    if (fluid.solid[i] || fluid.dynamicSolid[i]) return false;
    bool hasLiquid = fluid.fill[i] > kMinLiquidComponent;
    bool hasGas = gas.amount[i] > GAS_MIN_AMOUNT;
    if (!hasLiquid && !hasGas) return false;
    if (!reactionCatalystPresent(fluid, gas, index, def.conditions.catalyst)) return false;

    int x = index % GW, y = index / GW;
    float T = hasLiquid ? ThermalEngine::liquidTempK(fluid, index) : ThermalEngine::gasTempK(gas, index);
    if (!(T > 0.0f) || !std::isfinite(T))
        T = ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
    float P = gas.pressurePa(index);
    if (!std::isfinite(T) || !std::isfinite(P)) return false;
    if (!reactionConditionsMatch(def.conditions, T, P)) return false;

    LiquidComponentView liquidView;
    GasComponentView gasView;
    if (!reactionReadLiquidOccupancy(fluid, index, liquidView)) return false;
    if (!reactionReadGasOccupancy(gas, index, gasView)) return false;
    Inventory liquidInv = fromView(liquidView, fluid.fill[i], fluid.liquidHeat[i], false);
    Inventory gasInv = fromGasView(gasView, gas.amount[i], gas.heat[i], false);
    float extent = planExtent(def, liquidInv, gasInv, dt);
    if (!(extent > kExtentEps)) return false;
    Inventory committedL = liquidInv;
    Inventory committedG = gasInv;
    if (!applyExtent(def, committedL, committedG, extent)) return false;

    LiquidComponentView nextL = toView(committedL);
    GasComponentView nextG = toGasView(committedG);
    float oldFill = fluid.fill[i];
    float oldGas = gas.amount[i];
    if (!reactionCommitLiquidOccupancy(fluid, index, nextL)) return false;
    if (!reactionCommitGasOccupancy(gas, index, nextG)) {
        (void)reactionCommitLiquidOccupancy(fluid, index, liquidView);
        fluid.fill[i] = oldFill;
        return false;
    }
    fluid.expectedVolume += static_cast<double>(fluid.fill[i] - oldFill);
    gas.expectedAmount += static_cast<double>(gas.amount[i] - oldGas);
    if (fluid.fill[i] > kMinLiquidComponent)
        fluid.liquidHeat[i] = committedL.heat;
    else
        fluid.liquidHeat[i] = 0.0f;
    if (!(fluid.liquidHeat[i] >= 0.0f) || !std::isfinite(fluid.liquidHeat[i]))
        fluid.liquidHeat[i] = 0.0f;
    if (gas.amount[i] > GAS_MIN_AMOUNT)
        gas.heat[i] = committedG.heat;
    else
        gas.heat[i] = 0.0f;
    if (!(gas.heat[i] >= 0.0f) || !std::isfinite(gas.heat[i]))
        gas.heat[i] = 0.0f;
    fluid.wakeChunkAtCell(x, y);
    gas.wakeAt(x, y);
    thermal.wakeCell(x, y);
    return true;
}

} // namespace

void ReactionEngine::simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, float dt)
{
    if (reactionCount() <= 0) return;
    if (!(dt > 0.0f) || !std::isfinite(dt)) return;
    int nTable = reactionTableSize();
    ReactionDefinition const *table = builtinReactionTable();
    int x0 = std::max(0, fluid.solveX0);
    int y0 = std::max(0, fluid.solveY0);
    int x1 = std::min(GW - 1, fluid.solveX1);
    int y1 = std::min(GH - 1, fluid.solveY1);
    for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
        int index = FluidEngine::ci(x, y);
        if (fluid.fill[static_cast<size_t>(index)] < kMinLiquidComponent
            && gas.amount[static_cast<size_t>(index)] < GAS_MIN_AMOUNT)
            continue;
        for (int r = 0; r < nTable; ++r) {
            if (table[r].id == REACTION_NONE) continue;
            (void)tryReactCell(fluid, rigid, gas, thermal, index, table[r], dt);
        }
    }
}

void runReactionEngineSanityCheck() {
    std::ofstream out(miscFile("reaction_engine_sanity.tsv"));
    int passed = 0, failed = 0;
    auto emit = [&](char const *name, bool ok, std::string const &detail) {
        out << "test\t" << name << '\t' << (ok ? "PASS" : "FAIL") << '\t' << detail << '\n';
        if (ok) ++passed;
        else ++failed;
    };
    auto near = [](float a, float b) {
        return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 1.0e-4f;
    };

    constexpr SubstanceId kA = 100;
    constexpr SubstanceId kB = 101;
    constexpr SubstanceId kC = 102;

    ReactionDefinition rx{};
    rx.id = 1;
    rx.internalName = "synthetic_2A_B_C";
    rx.reactants[0] = {kA, MatterPhase::Liquid, 2.0f};
    rx.reactants[1] = {kB, MatterPhase::Liquid, 1.0f};
    rx.products[0] = {kC, MatterPhase::Liquid, 1.0f};
    rx.reactantCount = 2;
    rx.productCount = 1;
    rx.energyChangeJPerExtent = -10.0f;
    rx.maxExtentPerSecond = 1000.0f;

    Inventory inv{};
    inv.allowUnregisteredIds = true;
    inv.maxFill = 1.0e9f;
    inv.items[0] = {kA, 10.0f};
    inv.items[1] = {kB, 3.0f};
    inv.count = 2;
    inv.heat = 100.0f;
    recountFill(inv);

    float lim = limitingExtent(rx, inv);
    emit("limiting_extent", near(lim, 3.0f), "extent=" + std::to_string(lim));

    Inventory after = inv;
    bool ok = applyExtent(rx, after, 3.0f);
    emit("extent_3_commit",
        ok && near(amountOf(after, kA), 4.0f) && near(amountOf(after, kB), 0.0f)
            && near(amountOf(after, kC), 3.0f),
        "A=" + std::to_string(amountOf(after, kA))
            + " B=" + std::to_string(amountOf(after, kB))
            + " C=" + std::to_string(amountOf(after, kC)));
    emit("exothermic_heat", ok && near(after.heat, 130.0f), "heat=" + std::to_string(after.heat));

    rx.maxExtentPerSecond = 1.0f;
    Inventory rateInv = inv;
    float eRate = planExtent(rx, rateInv, 1.0f);
    Inventory rateAfter = rateInv;
    bool rateOk = applyExtent(rx, rateAfter, eRate);
    emit("rate_limit_extent_1",
        rateOk && near(eRate, 1.0f) && near(amountOf(rateAfter, kA), 8.0f)
            && near(amountOf(rateAfter, kB), 2.0f) && near(amountOf(rateAfter, kC), 1.0f),
        "e=" + std::to_string(eRate));

    Inventory shortB = inv;
    shortB.maxFill = 1.0e9f;
    shortB.items[1].amount = 0.0f;
    recountFill(shortB);
    Inventory shortTrial = shortB;
    emit("insufficient_reactant_no_mutation",
        limitingExtent(rx, shortB) <= kExtentEps && applyExtent(rx, shortTrial, 1.0f) == false
            && near(amountOf(shortB, kA), amountOf(shortTrial, kA)),
        "");

    ReactionDefinition gasRx = rx;
    gasRx.products[0].requiredPhase = MatterPhase::Gas;
    emit("storage_ok_synthetic_gas_product", reactionStorageOk(gasRx, true), "");

    Inventory mixL{};
    mixL.allowUnregisteredIds = true;
    mixL.maxFill = 1.0e9f;
    mixL.items[0] = {kA, 4.0f};
    mixL.count = 1;
    mixL.heat = 10.0f;
    recountFill(mixL);
    Inventory mixG{};
    mixG.allowUnregisteredIds = true;
    mixG.gas = true;
    mixG.maxFill = 1.0e6f;
    mixG.items[0] = {kB, 2.0f};
    mixG.count = 1;
    mixG.heat = 5.0f;
    recountFill(mixG);
    ReactionDefinition mixRx{};
    mixRx.id = 2;
    mixRx.internalName = "synthetic_liquidA_gasB_gasC";
    mixRx.reactants[0] = {kA, MatterPhase::Liquid, 1.0f};
    mixRx.reactants[1] = {kB, MatterPhase::Gas, 1.0f};
    mixRx.products[0] = {kC, MatterPhase::Gas, 1.0f};
    mixRx.reactantCount = 2;
    mixRx.productCount = 1;
    mixRx.maxExtentPerSecond = 1000.0f;
    Inventory mixL1 = mixL, mixG1 = mixG;
    bool mixOk = applyExtent(mixRx, mixL1, mixG1, 2.0f);
    emit("mixed_phase_liquid_gas_to_gas",
        mixOk && near(amountOf(mixL1, kA), 2.0f) && near(amountOf(mixG1, kB), 0.0f)
            && near(amountOf(mixG1, kC), 2.0f)
            && near(amountOf(mixL1, kC), 0.0f),
        "liqA=" + std::to_string(amountOf(mixL1, kA))
            + " gasB=" + std::to_string(amountOf(mixG1, kB))
            + " gasC=" + std::to_string(amountOf(mixG1, kC)));

    ReactionDefinition solidRx = rx;
    solidRx.reactants[0].requiredPhase = MatterPhase::Solid;
    emit("unsupported_solid_reactant", !reactionStorageOk(solidRx, true), "");

    Inventory full{};
    full.allowUnregisteredIds = true;
    full.maxFill = 1.0e9f;
    full.heat = 50.0f;
    for (int n = 0; n < kMaxLiquidComponents; ++n)
        full.items[n] = {static_cast<SubstanceId>(200 + n), 0.2f};
    full.count = kMaxLiquidComponents;
    recountFill(full);
    Inventory fullTrial = full;
    ReactionDefinition newSlot = rx;
    newSlot.reactants[0] = {static_cast<SubstanceId>(200), MatterPhase::Liquid, 1.0f};
    newSlot.reactants[1] = {static_cast<SubstanceId>(201), MatterPhase::Liquid, 1.0f};
    newSlot.products[0] = {kC, MatterPhase::Liquid, 1.0f};
    newSlot.reactantCount = 2;
    newSlot.productCount = 1;
    emit("no_slot_no_mutation",
        !applyExtent(newSlot, fullTrial, 0.1f)
            && near(amountOf(full, static_cast<SubstanceId>(200)),
                amountOf(fullTrial, static_cast<SubstanceId>(200))),
        "count=" + std::to_string(fullTrial.count));

    ReactionDefinition endo = rx;
    endo.energyChangeJPerExtent = 80.0f;
    endo.maxExtentPerSecond = 1000.0f;
    Inventory cold{};
    cold.allowUnregisteredIds = true;
    cold.maxFill = 1.0e9f;
    cold.items[0] = {kA, 10.0f};
    cold.items[1] = {kB, 3.0f};
    cold.count = 2;
    cold.heat = 80.0f;
    recountFill(cold);
    float eEndo = planExtent(endo, cold, 10.0f);
    Inventory endoAfter = cold;
    bool endoOk = applyExtent(endo, endoAfter, eEndo);
    emit("endothermic_heat_bound",
        endoOk && near(eEndo, 1.0f) && endoAfter.heat >= -1.0e-4f,
        "e=" + std::to_string(eEndo) + " heat=" + std::to_string(endoAfter.heat));

    emit("empty_registry_fast_path", reactionCount() == 0, "count=" + std::to_string(reactionCount()));

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
