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
};

bool idAllowed(Inventory const &inv, SubstanceId id) {
    if (inv.allowUnregisteredIds) return id != SUBSTANCE_NONE;
    return validLiquidComponentId(id);
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

void recountFill(Inventory &inv) {
    float s = 0.0f;
    int w = 0;
    for (int n = 0; n < inv.count; ++n) {
        if (inv.items[n].amount > kMinLiquidComponent && inv.items[n].id != SUBSTANCE_NONE) {
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
        if (p.requiredPhase == MatterPhase::Gas || p.requiredPhase == MatterPhase::Solid
            || p.requiredPhase == MatterPhase::Plasma)
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

float limitingExtent(ReactionDefinition const &def, Inventory const &inv) {
    float extent = 1.0e30f;
    bool any = false;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        any = true;
        float avail = amountOf(inv, p.substance);
        if (!(p.coefficient > 0.0f) || !(avail > kExtentEps)) return 0.0f;
        extent = std::min(extent, avail / p.coefficient);
    }
    if (!any || !std::isfinite(extent) || extent <= kExtentEps) return 0.0f;
    return extent;
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

bool applyExtent(ReactionDefinition const &def, Inventory &inv, float extent) {
    if (!(extent > kExtentEps) || !std::isfinite(extent)) return false;
    Inventory trial = inv;
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        if (!consume(trial, p.substance, p.coefficient * extent)) return false;
    }
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.products[n];
        if (!reactionParticipantUsed(p)) continue;
        if (!produce(trial, p.substance, p.coefficient * extent)) return false;
    }
    float q = reactionHeatReleasedJ(def, extent);
    float heat1 = trial.heat + q;
    if (!std::isfinite(heat1) || heat1 < 0.0f) return false;
    trial.heat = heat1;
    inv = trial;
    return true;
}

float planExtent(ReactionDefinition const &def, Inventory const &inv, float dt) {
    if (!reactionStorageOk(def, inv.allowUnregisteredIds)) return 0.0f;
    float extent = limitingExtent(def, inv);
    extent = rateLimit(def, dt, extent);
    extent = heatLimit(def, inv.heat, extent);
    if (!(extent > kExtentEps)) return 0.0f;
    Inventory trial = inv;
    // Shrink if products cannot fit (fill cap / slot cap).
    for (int attempt = 0; attempt < 8; ++attempt) {
        trial = inv;
        if (applyExtent(def, trial, extent)) return extent;
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
    for (int n = 0; n < view.count && n < kMaxLiquidComponents; ++n) {
        if (view.items[n].amount > kMinLiquidComponent)
            inv.items[inv.count++] = view.items[n];
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

bool tryReactCell(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int index, ReactionDefinition const &def, float dt)
{
    if (!validReaction(def.id)) return false;
    if (!reactionStorageOk(def, false)) return false;
    if (index < 0 || index >= GW * GH) return false;
    size_t i = static_cast<size_t>(index);
    if (fluid.solid[i] || fluid.dynamicSolid[i]) return false;
    if (fluid.fill[i] < kMinLiquidComponent) return false;
    if (!reactionCatalystPresent(fluid, index, def.conditions.catalyst)) return false;

    int x = index % GW, y = index / GW;
    float T = ThermalEngine::liquidTempK(fluid, index);
    if (!(T > 0.0f) || !std::isfinite(T))
        T = ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
    float P = gas.pressurePa(index);
    if (!std::isfinite(T) || !std::isfinite(P)) return false;
    if (!reactionConditionsMatch(def.conditions, T, P)) return false;

    LiquidComponentView view;
    if (!reactionReadLiquidOccupancy(fluid, index, view)) return false;
    Inventory inv = fromView(view, fluid.fill[i], fluid.liquidHeat[i], false);
    float extent = planExtent(def, inv, dt);
    if (!(extent > kExtentEps)) return false;
    Inventory committed = inv;
    if (!applyExtent(def, committed, extent)) return false;

    LiquidComponentView next = toView(committed);
    float oldFill = fluid.fill[i];
    if (!reactionCommitLiquidOccupancy(fluid, index, next)) return false;
    fluid.expectedVolume += static_cast<double>(fluid.fill[i] - oldFill);
    if (fluid.fill[i] > kMinLiquidComponent)
        fluid.liquidHeat[i] = committed.heat;
    else
        fluid.liquidHeat[i] = 0.0f;
    if (!(fluid.liquidHeat[i] >= 0.0f) || !std::isfinite(fluid.liquidHeat[i]))
        fluid.liquidHeat[i] = 0.0f;
    fluid.wakeChunkAtCell(x, y);
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
        if (fluid.fill[static_cast<size_t>(index)] < kMinLiquidComponent) continue;
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
    emit("unsupported_gas_product", !reactionStorageOk(gasRx, true), "");

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
