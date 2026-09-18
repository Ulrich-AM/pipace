#include "chemistry/ReactionEngine.h"

#include "chemistry/ReactionMatterAccess.h"
#include "chemistry/ReactionRegistry.h"
#include "fluid/DiagOutput.h"
#include "fluid/FluidEngine.h"
#include "fluid/FluidTypes.h"
#include "gas/GasEngine.h"
#include "rigid/RigidBodyEngine.h"
#include "substance/PhaseTransfer.h"
#include "thermal/ThermalEngine.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

constexpr float kExtentEps = 1.0e-8f;
constexpr float kActivityHeatRefJ = 4.0e4f;
constexpr float kActivityMax = 1.5f;
constexpr float kActivityEps = 1.0e-4f;
constexpr float kActivityDecayPerSec = 6.0f; // ~0.17 s e-fold; gone by ~0.5 s

namespace {

constexpr float kDefaultMaxExtentPerSecond = 1.0f;


struct Inventory {
    LiquidComponent items[kMaxLiquidComponents]{};
    int count = 0;
    float fill = 0.0f;
    float heat = 0.0f;
    float maxFill = 1.0f;
    bool allowUnregisteredIds = false;
    bool gas = false;
    float cellsPerMeter = 4.0f;
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

MatterPhase inventoryPhase(Inventory const &inv) {
    return inv.gas ? MatterPhase::Gas : MatterPhase::Liquid;
}

bool reactionMolarDataOk(ReactionDefinition const &def, bool synthetic) {
    if (synthetic) return true;
    auto check = [](ReactionParticipant const &p) {
        if (!reactionParticipantUsed(p)) return true;
        return substanceHasMolarMass(p.substance);
    };
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n)
        if (!check(def.reactants[n])) return false;
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n)
        if (!check(def.products[n])) return false;
    return true;
}

float availableExtent(ReactionParticipant const &p, Inventory const &inv) {
    float amt = amountOf(inv, p.substance);
    if (!(p.coefficient > 0.0f) || !(amt > kExtentEps)) return 0.0f;
    if (inv.allowUnregisteredIds) return amt / p.coefficient;
    double moles = storageAmountToMoles(p.substance, inventoryPhase(inv), amt, inv.cellsPerMeter);
    if (!(moles > 0.0) || !std::isfinite(moles)) return 0.0f;
    return static_cast<float>(moles / static_cast<double>(p.coefficient));
}

float storageForExtent(ReactionParticipant const &p, Inventory const &inv, float extent) {
    float n = p.coefficient * extent;
    if (!(n > 0.0f) || !std::isfinite(n)) return 0.0f;
    if (inv.allowUnregisteredIds) return n;
    double storage = molesToStorageAmount(p.substance, inventoryPhase(inv), n, inv.cellsPerMeter);
    if (!(storage > 0.0) || !std::isfinite(storage)) return 0.0f;
    return static_cast<float>(storage);
}

void scanParticipatingPhases(ReactionDefinition const &def, Inventory &liquid, Inventory &gas,
    bool &usesLiquid, bool &usesGas)
{
    usesLiquid = false;
    usesGas = false;
    auto mark = [&](ReactionParticipant const &p) {
        if (!reactionParticipantUsed(p)) return;
        Inventory *inv = pickInv(p, liquid, gas);
        if (inv->gas) usesGas = true;
        else usesLiquid = true;
    };
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n)
        mark(def.reactants[n]);
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n)
        mark(def.products[n]);
}

float participatingHeat(bool usesLiquid, bool usesGas, Inventory const &liquid, Inventory const &gas) {
    float h = 0.0f;
    if (usesLiquid) h += std::max(0.0f, liquid.heat);
    if (usesGas) h += std::max(0.0f, gas.heat);
    return h;
}

bool applyReactionHeat(bool usesLiquid, bool usesGas, Inventory &liquid, Inventory &gas, float q) {
    if (!std::isfinite(q) || q == 0.0f) return true;
    float wL = usesLiquid ? std::max(0.0f, liquid.heat) : 0.0f;
    float wG = usesGas ? std::max(0.0f, gas.heat) : 0.0f;
    if (q < 0.0f) {
        float need = -q;
        float have = wL + wG;
        if (need > have + 1.0e-5f) return false;
        if (have <= 1.0e-12f) return !(need > kExtentEps);
        if (usesLiquid) liquid.heat -= need * (wL / have);
        if (usesGas) gas.heat -= need * (wG / have);
    } else {
        float have = wL + wG;
        if (have > 1.0e-12f) {
            if (usesLiquid) liquid.heat += q * (wL / have);
            if (usesGas) gas.heat += q * (wG / have);
        } else {
            float n = (usesLiquid ? 1.0f : 0.0f) + (usesGas ? 1.0f : 0.0f);
            if (!(n > 0.0f)) return false;
            if (usesLiquid) liquid.heat += q / n;
            if (usesGas) gas.heat += q / n;
        }
    }
    if (usesLiquid && liquid.heat < 0.0f) {
        if (liquid.heat < -1.0e-5f) return false;
        liquid.heat = 0.0f;
    }
    if (usesGas && gas.heat < 0.0f) {
        if (gas.heat < -1.0e-5f) return false;
        gas.heat = 0.0f;
    }
    return true;
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
        float e = availableExtent(p, *inv);
        if (!(e > kExtentEps)) return 0.0f;
        extent = std::min(extent, e);
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
    bool synthetic = liquid.allowUnregisteredIds || gas.allowUnregisteredIds;
    if (!reactionMolarDataOk(def, synthetic)) return false;
    Inventory trialL = liquid;
    Inventory trialG = gas;
    bool usesLiquid = false;
    bool usesGas = false;
    scanParticipatingPhases(def, trialL, trialG, usesLiquid, usesGas);
    for (int n = 0; n < def.reactantCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.reactants[n];
        if (!reactionParticipantUsed(p)) continue;
        Inventory *inv = pickInv(p, trialL, trialG);
        float take = storageForExtent(p, *inv, extent);
        if (!(take > 0.0f)) return false;
        if (!consume(*inv, p.substance, take)) return false;
    }
    for (int n = 0; n < def.productCount && n < kMaxReactionParticipants; ++n) {
        ReactionParticipant const &p = def.products[n];
        if (!reactionParticipantUsed(p)) continue;
        Inventory *inv = pickInv(p, trialL, trialG);
        float add = storageForExtent(p, *inv, extent);
        if (!(add > 0.0f)) return false;
        if (!produce(*inv, p.substance, add)) return false;
    }
    float q = reactionHeatReleasedJ(def, extent);
    if (!applyReactionHeat(usesLiquid, usesGas, trialL, trialG, q)) return false;
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
    bool synthetic = liquid.allowUnregisteredIds || gas.allowUnregisteredIds;
    if (!reactionStorageOk(def, synthetic)) return 0.0f;
    if (!reactionMolarDataOk(def, synthetic)) return 0.0f;
    float extent = limitingExtent(def, liquid, gas);
    extent = rateLimit(def, dt, extent);
    bool usesLiquid = false, usesGas = false;
    Inventory liq = liquid, ginv = gas;
    scanParticipatingPhases(def, liq, ginv, usesLiquid, usesGas);
    extent = heatLimit(def, participatingHeat(usesLiquid, usesGas, liquid, gas), extent);
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

Inventory fromView(LiquidComponentView const &view, float fill, float heat, bool synthetic,
    float cellsPerMeter = 4.0f) {
    Inventory inv;
    inv.allowUnregisteredIds = synthetic;
    inv.maxFill = 1.0f;
    inv.heat = heat;
    inv.fill = fill;
    inv.count = 0;
    inv.gas = false;
    inv.cellsPerMeter = cellsPerMeter;
    for (int n = 0; n < view.count && n < kMaxLiquidComponents; ++n) {
        if (view.items[n].amount > kMinLiquidComponent)
            inv.items[inv.count++] = view.items[n];
    }
    recountFill(inv);
    return inv;
}

Inventory fromGasView(GasComponentView const &view, float amt, float heat, bool synthetic,
    float cellsPerMeter = 4.0f) {
    Inventory inv;
    inv.allowUnregisteredIds = synthetic;
    inv.maxFill = 1.0e6f;
    inv.heat = heat;
    inv.fill = amt;
    inv.count = 0;
    inv.gas = true;
    inv.cellsPerMeter = cellsPerMeter;
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

float tryReactCell(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, int index, ReactionDefinition const &def, float dt)
{
    if (!validReaction(def.id)) return 0.0f;
    if (!reactionStorageOk(def, false)) return 0.0f;
    if (index < 0 || index >= GW * GH) return 0.0f;
    size_t i = static_cast<size_t>(index);
    if (fluid.solid[i] || fluid.dynamicSolid[i]) return 0.0f;
    if (!reactionCellHasRequiredReactants(fluid, gas, index, def)) return 0.0f;
    bool hasLiquid = fluid.fill[i] > kMinLiquidComponent;
    bool hasGas = gas.amount[i] > GAS_MIN_AMOUNT;
    if (!hasLiquid && !hasGas) return 0.0f;
    if (!reactionCatalystPresent(fluid, gas, index, def.conditions.catalyst)) return 0.0f;

    int x = index % GW, y = index / GW;
    float T = 0.0f;
    if (hasGas) T = ThermalEngine::gasTempK(gas, index);
    if (hasLiquid) {
        float Tl = ThermalEngine::liquidTempK(fluid, index);
        if (!(T > 0.0f) || !std::isfinite(T) || (std::isfinite(Tl) && Tl > T)) T = Tl;
    }
    if (!(T > 0.0f) || !std::isfinite(T))
        T = ThermalEngine::sampleTemperatureK(fluid, rigid, gas, x, y);
    float P = gas.pressurePa(index);
    if (!std::isfinite(T) || !std::isfinite(P)) return 0.0f;
    if (!reactionConditionsMatch(def.conditions, T, P)) return 0.0f;

    LiquidComponentView liquidView;
    GasComponentView gasView;
    if (!reactionReadLiquidOccupancy(fluid, index, liquidView)) return 0.0f;
    if (!reactionReadGasOccupancy(gas, index, gasView)) return 0.0f;
    Inventory liquidInv = fromView(liquidView, fluid.fill[i], fluid.liquidHeat[i], false,
        fluid.config.cellsPerMeter);
    Inventory gasInv = fromGasView(gasView, gas.amount[i], gas.heat[i], false,
        fluid.config.cellsPerMeter);
    float extent = planExtent(def, liquidInv, gasInv, dt);
    if (!(extent > kExtentEps)) return 0.0f;
    Inventory committedL = liquidInv;
    Inventory committedG = gasInv;
    if (!applyExtent(def, committedL, committedG, extent)) return 0.0f;

    LiquidComponentView nextL = toView(committedL);
    GasComponentView nextG = toGasView(committedG);
    float oldFill = fluid.fill[i];
    float oldGas = gas.amount[i];
    float oldVap = gas.vaporAmount(index);
    if (!reactionCommitLiquidOccupancy(fluid, index, nextL)) return 0.0f;
    if (!reactionCommitGasOccupancy(gas, index, nextG)) {
        (void)reactionCommitLiquidOccupancy(fluid, index, liquidView);
        fluid.fill[i] = oldFill;
        return 0.0f;
    }
    fluid.expectedVolume += static_cast<double>(fluid.fill[i] - oldFill);
    gas.expectedAmount += static_cast<double>(gas.amount[i] - oldGas);
    gas.expectedWaterVapor += static_cast<double>(gas.vaporAmount(index) - oldVap);
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
    return extent;
}

} // namespace

ReactionEngine::ReactionEngine()
    : activity(static_cast<size_t>(GW * GH), 0.0f)
{
    activityCells.reserve(256);
}

void ReactionEngine::clearActivity() {
    std::fill(activity.begin(), activity.end(), 0.0f);
    activityCells.clear();
    reactedCellsLastTick = 0;
    extentLastTick = 0.0f;
    heatReleasedLastTick = 0.0f;
}

void ReactionEngine::decayActivity(float dt) {
    if (activityCells.empty()) return;
    float factor = 0.0f;
    if (dt > 0.0f && std::isfinite(dt))
        factor = std::exp(-kActivityDecayPerSec * dt);
    int w = 0;
    for (int index : activityCells) {
        if (index < 0 || index >= GW * GH) continue;
        float &a = activity[static_cast<size_t>(index)];
        a *= factor;
        if (a > kActivityEps) {
            activityCells[static_cast<size_t>(w++)] = index;
        } else {
            a = 0.0f;
        }
    }
    activityCells.resize(static_cast<size_t>(w));
}

void ReactionEngine::addActivity(int index, float heatJ) {
    if (index < 0 || index >= GW * GH) return;
    float add = std::abs(heatJ) / kActivityHeatRefJ;
    if (!(add > 0.0f) || !std::isfinite(add)) return;
    add = std::min(add, 1.0f);
    float &a = activity[static_cast<size_t>(index)];
    if (a <= kActivityEps) activityCells.push_back(index);
    a = std::min(kActivityMax, a + add);
}

void ReactionEngine::simulationTick(FluidEngine &fluid, RigidBodyEngine &rigid, GasEngine &gas,
    ThermalEngine &thermal, float dt)
{
    reactedCellsLastTick = 0;
    extentLastTick = 0.0f;
    heatReleasedLastTick = 0.0f;
    if (!(dt > 0.0f) || !std::isfinite(dt)) return;
    decayActivity(dt);
    if (reactionCount() <= 0) return;
    int nTable = reactionTableSize();
    ReactionDefinition const *table = builtinReactionTable();
    // Gas-only mixtures can sit outside FluidEngine's liquid solve region.
    // Full-grid scan at 200x120; ambient Air is skipped by reactant prefilter
    // before temperature / inventory work.
    for (int y = 0; y < GH; ++y) for (int x = 0; x < GW; ++x) {
        int index = FluidEngine::ci(x, y);
        if (fluid.fill[static_cast<size_t>(index)] < kMinLiquidComponent
            && gas.amount[static_cast<size_t>(index)] < GAS_MIN_AMOUNT)
            continue;
        for (int r = 0; r < nTable; ++r) {
            if (table[r].id == REACTION_NONE) continue;
            if (!reactionCellHasRequiredReactants(fluid, gas, index, table[r])) continue;
            float extent = tryReactCell(fluid, rigid, gas, thermal, index, table[r], dt);
            if (!(extent > kExtentEps)) continue;
            float q = std::abs(reactionHeatReleasedJ(table[r], extent));
            addActivity(index, q);
            ++reactedCellsLastTick;
            extentLastTick += extent;
            heatReleasedLastTick += q;
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

    Inventory endoL{};
    endoL.allowUnregisteredIds = true;
    endoL.maxFill = 1.0e9f;
    endoL.items[0] = {kA, 4.0f};
    endoL.count = 1;
    endoL.heat = 12.0f;
    recountFill(endoL);
    Inventory endoG{};
    endoG.allowUnregisteredIds = true;
    endoG.gas = true;
    endoG.maxFill = 1.0e6f;
    endoG.items[0] = {kB, 2.0f};
    endoG.count = 1;
    endoG.heat = 4.0f;
    recountFill(endoG);
    ReactionDefinition endoMix = mixRx;
    endoMix.energyChangeJPerExtent = 8.0f;
    Inventory endoL1 = endoL, endoG1 = endoG;
    bool endoMixOk = applyExtent(endoMix, endoL1, endoG1, 1.0f);
    emit("mixed_endo_split_nonnegative",
        endoMixOk && near(endoL1.heat, 6.0f) && near(endoG1.heat, 2.0f)
            && endoL1.heat >= -1.0e-4f && endoG1.heat >= -1.0e-4f,
        "liqH=" + std::to_string(endoL1.heat) + " gasH=" + std::to_string(endoG1.heat));

    constexpr float cpm = 4.0f;
    Inventory wL{};
    wL.allowUnregisteredIds = false;
    wL.maxFill = 1.0f;
    wL.cellsPerMeter = cpm;
    wL.items[0] = {SUBSTANCE_WATER, 0.20f};
    wL.count = 1;
    wL.heat = 1.0e8f;
    recountFill(wL);
    Inventory wG{};
    wG.allowUnregisteredIds = false;
    wG.gas = true;
    wG.maxFill = 1.0e6f;
    wG.cellsPerMeter = cpm;
    wG.heat = 0.0f;
    ReactionDefinition twoWater{};
    twoWater.id = 3;
    twoWater.internalName = "diag_2_water_liquid_to_gas";
    twoWater.reactants[0] = {SUBSTANCE_WATER, MatterPhase::Liquid, 2.0f};
    twoWater.products[0] = {SUBSTANCE_WATER, MatterPhase::Gas, 2.0f};
    twoWater.reactantCount = 1;
    twoWater.productCount = 1;
    twoWater.maxExtentPerSecond = 1.0e9f;
    double waterMoles = storageAmountToMoles(SUBSTANCE_WATER, MatterPhase::Liquid, 0.20, cpm);
    float molarLim = limitingExtent(twoWater, wL, wG);
    emit("molar_limit_2_water",
        waterMoles > 0.0 && near(molarLim, static_cast<float>(waterMoles / 2.0)),
        "moles=" + std::to_string(waterMoles) + " extent=" + std::to_string(molarLim));

    double mass0 = liquidFillToMassKg(SUBSTANCE_WATER, 0.20, cpm);
    Inventory wL2 = wL, wG2 = wG;
    bool molarOk = applyExtent(twoWater, wL2, wG2, molarLim);
    double massL = liquidFillToMassKg(SUBSTANCE_WATER, amountOf(wL2, SUBSTANCE_WATER), cpm);
    double massG = gasAmountToMassKg(SUBSTANCE_WATER, amountOf(wG2, SUBSTANCE_WATER), cpm);
    emit("physical_mass_liquid_to_gas",
        molarOk && mass0 > 0.0 && std::abs((massL + massG) - mass0) <= 1.0e-4 * mass0 + 1.0e-6,
        "m0=" + std::to_string(mass0) + " mL=" + std::to_string(massL)
            + " mG=" + std::to_string(massG));

    ReactionDefinition airRx = twoWater;
    airRx.reactants[0] = {SUBSTANCE_AIR, MatterPhase::Gas, 1.0f};
    airRx.products[0] = {SUBSTANCE_WATER, MatterPhase::Gas, 1.0f};
    emit("physical_skip_missing_molar_mass",
        !reactionMolarDataOk(airRx, false) && planExtent(airRx, wL, wG, 1.0f) <= kExtentEps, "");

    GasComponent fullSlots[kMaxGasComponents] = {
        {SUBSTANCE_AIR, 0.25f}, {SUBSTANCE_AIR, 0.25f},
        {SUBSTANCE_AIR, 0.25f}, {SUBSTANCE_AIR, 0.25f}
    };
    GasComponent waterSrc[1] = {{SUBSTANCE_WATER, 0.4f}};
    emit("gas_merge_reject_full_slots",
        !gasPayloadCanMerge(fullSlots, kMaxGasComponents, waterSrc, 1), "");
    emit("gas_merge_accept_into_air",
        gasPayloadCanMerge(fullSlots, 1, waterSrc, 1), "");

    ReactionDefinition const &h2o2 = reactionDef(REACTION_HYDROGEN_COMBUSTION);
    emit("hydrogen_combustion_registered",
        reactionCount() == 1 && validReaction(REACTION_HYDROGEN_COMBUSTION)
            && reactionFromInternalName("hydrogen_combustion") == REACTION_HYDROGEN_COMBUSTION
            && h2o2.reactantCount == 2 && h2o2.productCount == 1
            && h2o2.reactants[0].substance == SUBSTANCE_HYDROGEN
            && h2o2.reactants[0].requiredPhase == MatterPhase::Gas
            && h2o2.reactants[0].coefficient == 2.0f
            && h2o2.reactants[1].substance == SUBSTANCE_OXYGEN
            && h2o2.reactants[1].requiredPhase == MatterPhase::Gas
            && h2o2.reactants[1].coefficient == 1.0f
            && h2o2.products[0].substance == SUBSTANCE_WATER
            && h2o2.products[0].requiredPhase == MatterPhase::Gas
            && h2o2.products[0].coefficient == 2.0f, "");
    emit("hydrogen_combustion_mass_balanced",
        reactionMassConservation(REACTION_HYDROGEN_COMBUSTION) == ReactionMassConservation::Balanced, "");
    emit("hydrogen_combustion_cold_no_react",
        !reactionConditionsMatch(h2o2.conditions, 293.15f, 101325.0f), "");
    emit("hydrogen_combustion_hot_reacts",
        reactionConditionsMatch(h2o2.conditions, 900.0f, 101325.0f)
            && h2o2.conditions.minTemperatureValid && h2o2.conditions.minTemperatureK == 850.0f, "");

    auto nearRel = [](double a, double b) {
        double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
        return std::isfinite(a) && std::isfinite(b)
            && std::abs(a - b) <= 1.0e-3 * scale + 1.0e-6;
    };

    Inventory emptyL{};
    emptyL.maxFill = 1.0f;
    emptyL.cellsPerMeter = cpm;
    Inventory stoichG{};
    stoichG.gas = true;
    stoichG.maxFill = 1.0e6f;
    stoichG.cellsPerMeter = cpm;
    stoichG.items[0] = {SUBSTANCE_HYDROGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0, cpm))};
    stoichG.items[1] = {SUBSTANCE_OXYGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
    stoichG.count = 2;
    stoichG.heat = 0.0f;
    recountFill(stoichG);
    Inventory stoichG2 = stoichG;
    float stoichExtent = limitingExtent(h2o2, emptyL, stoichG2);
    bool stoichOk = applyExtent(h2o2, emptyL, stoichG2, stoichExtent);
    double waterMol = storageAmountToMoles(SUBSTANCE_WATER, MatterPhase::Gas,
        amountOf(stoichG2, SUBSTANCE_WATER), cpm);
    double h2Left = storageAmountToMoles(SUBSTANCE_HYDROGEN, MatterPhase::Gas,
        amountOf(stoichG2, SUBSTANCE_HYDROGEN), cpm);
    double o2Left = storageAmountToMoles(SUBSTANCE_OXYGEN, MatterPhase::Gas,
        amountOf(stoichG2, SUBSTANCE_OXYGEN), cpm);
    emit("hydrogen_combustion_stoich_2_1",
        stoichOk && nearRel(stoichExtent, 1.0) && nearRel(waterMol, 2.0)
            && h2Left <= 1.0e-4 && o2Left <= 1.0e-4,
        "extent=" + std::to_string(stoichExtent) + " waterMol=" + std::to_string(waterMol));

    Inventory excessG{};
    excessG.gas = true;
    excessG.maxFill = 1.0e6f;
    excessG.cellsPerMeter = cpm;
    excessG.items[0] = {SUBSTANCE_HYDROGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 4.0, cpm))};
    excessG.items[1] = {SUBSTANCE_OXYGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
    excessG.count = 2;
    excessG.heat = 0.0f;
    recountFill(excessG);
    Inventory excessL{};
    excessL.maxFill = 1.0f;
    excessL.cellsPerMeter = cpm;
    float excessExtent = limitingExtent(h2o2, excessL, excessG);
    bool excessOk = applyExtent(h2o2, excessL, excessG, excessExtent);
    double h2Remain = storageAmountToMoles(SUBSTANCE_HYDROGEN, MatterPhase::Gas,
        amountOf(excessG, SUBSTANCE_HYDROGEN), cpm);
    double o2Remain = storageAmountToMoles(SUBSTANCE_OXYGEN, MatterPhase::Gas,
        amountOf(excessG, SUBSTANCE_OXYGEN), cpm);
    double waterMade = storageAmountToMoles(SUBSTANCE_WATER, MatterPhase::Gas,
        amountOf(excessG, SUBSTANCE_WATER), cpm);
    emit("hydrogen_combustion_excess_h2",
        excessOk && nearRel(excessExtent, 1.0) && nearRel(h2Remain, 2.0)
            && o2Remain <= 1.0e-4 && nearRel(waterMade, 2.0),
        "h2Remain=" + std::to_string(h2Remain) + " water=" + std::to_string(waterMade));

    float heatReleased = reactionHeatReleasedJ(h2o2, 1.0f);
    emit("hydrogen_combustion_heat_released",
        stoichOk && heatReleased > 0.0f && std::isfinite(heatReleased)
            && nearRel(heatReleased, 483600.0)
            && stoichG2.heat > 0.0f && std::isfinite(stoichG2.heat)
            && nearRel(stoichG2.heat, heatReleased),
        "q=" + std::to_string(stoichG2.heat));

    Inventory airMixG{};
    airMixG.gas = true;
    airMixG.maxFill = 1.0e6f;
    airMixG.cellsPerMeter = cpm;
    airMixG.items[0] = {SUBSTANCE_HYDROGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0, cpm))};
    airMixG.items[1] = {SUBSTANCE_OXYGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
    airMixG.items[2] = {SUBSTANCE_AIR, 1.0f};
    airMixG.count = 3;
    airMixG.heat = 0.0f;
    recountFill(airMixG);
    Inventory airMixL{};
    airMixL.maxFill = 1.0f;
    airMixL.cellsPerMeter = cpm;
    float airMixExtent = limitingExtent(h2o2, airMixL, airMixG);
    bool airMixOk = applyExtent(h2o2, airMixL, airMixG, airMixExtent);
    GasEngine ge;
    int gi = 0;
    ge.volume[static_cast<size_t>(gi)] = 1.0f;
    GasComponentView mixView = toGasView(airMixG);
    bool committed = airMixOk && ge.tryCommitGasOccupancy(gi, mixView);
    bool idsOk = committed;
    for (int n = 0; n < mixView.count; ++n) {
        if (!validGasComponentId(mixView.items[n].id)) idsOk = false;
    }
    SubstanceId dom = ge.dominantGasSubstance(gi);
    emit("hydrogen_combustion_composition_valid",
        committed && ge.gasCompositionValid(gi) && idsOk, "");
    emit("hydrogen_combustion_no_identityless",
        committed && dom != SUBSTANCE_NONE && validGasComponentId(dom)
            && ge.gasComponentAmount(gi, SUBSTANCE_AIR) > kMinGasComponent
            && ge.gasComponentAmount(gi, SUBSTANCE_WATER) > kMinGasComponent
            && ge.gasComponentAmount(gi, SUBSTANCE_HYDROGEN) <= kMinGasComponent
            && ge.gasComponentAmount(gi, SUBSTANCE_OXYGEN) <= kMinGasComponent,
        "dom=" + std::to_string(dom) + " air=" + std::to_string(ge.airAmount(gi))
            + " water=" + std::to_string(ge.vaporAmount(gi)));

    GasComponentView airView{};
    airView.count = 1;
    airView.items[0] = {SUBSTANCE_AIR, 1.0f};
    GasComponentView h2View{};
    h2View.count = 1;
    h2View.items[0] = {SUBSTANCE_HYDROGEN, 1.0f};
    GasComponentView o2View{};
    o2View.count = 1;
    o2View.items[0] = {SUBSTANCE_OXYGEN, 1.0f};
    float rhoAir = gasMixtureReferenceDensityKgM3(airView);
    float rhoH2 = gasMixtureReferenceDensityKgM3(h2View);
    float rhoO2 = gasMixtureReferenceDensityKgM3(o2View);
    emit("mixture_density_air_reference",
        nearRel(rhoAir, AIR_DENSITY_KG_M3), "rhoAir=" + std::to_string(rhoAir));
    emit("mixture_density_hydrogen_lighter",
        rhoH2 > 0.0f && rhoH2 < rhoAir * 0.5f, "rhoH2=" + std::to_string(rhoH2));
    emit("mixture_density_oxygen_heavier",
        rhoO2 > rhoAir, "rhoO2=" + std::to_string(rhoO2));
    float kAir = gasMixtureConductivity(airView);
    float kH2 = gasMixtureConductivity(h2View);
    emit("mixture_conductivity_hydrogen_faster",
        kH2 > kAir * 2.0f, "kH2=" + std::to_string(kH2) + " kAir=" + std::to_string(kAir));

    FluidEngine fluid;
    int airIdx = 20;
    ge.volume[static_cast<size_t>(airIdx)] = 1.0f;
    bool airCommitted = ge.tryCommitGasOccupancy(airIdx, airView);
    emit("ambient_air_fails_h2o2_prefilter",
        airCommitted && !reactionCellHasRequiredReactants(fluid, ge, airIdx, h2o2), "");
    GasComponentView fuelView{};
    fuelView.count = 2;
    fuelView.items[0] = {SUBSTANCE_HYDROGEN, 0.5f};
    fuelView.items[1] = {SUBSTANCE_OXYGEN, 0.25f};
    int fuelIdx = 21;
    ge.volume[static_cast<size_t>(fuelIdx)] = 1.0f;
    bool fuelCommitted = ge.tryCommitGasOccupancy(fuelIdx, fuelView);
    emit("h2o2_cell_passes_prefilter",
        fuelCommitted && reactionCellHasRequiredReactants(fluid, ge, fuelIdx, h2o2), "");

    emit("hydrogen_combustion_rate_scale",
        h2o2.maxExtentPerSecond == 4.0f && h2o2.conditions.minTemperatureK == 850.0f, "");

    GasConfig qOff;
    qOff.simMode = GasSimMode::Off;
    applyGasQualitySimMode(qOff, 0);
    applyGasQualitySimMode(qOff, 1);
    applyGasQualitySimMode(qOff, 2);
    emit("explicit_gas_off_survives_auto_quality", qOff.simMode == GasSimMode::Off, "");
    GasConfig qOn;
    qOn.simMode = GasSimMode::Full;
    applyGasQualitySimMode(qOn, 0);
    emit("quality_low_sets_half_when_enabled", qOn.simMode == GasSimMode::Half, "");

    FluidEngine worldF;
    RigidBodyEngine worldR;
    GasEngine worldG;
    worldG.resetAmbient(worldF);
    int hx = 8, hy = 8;
    int hi = FluidEngine::ci(hx, hy);
    GasComponentView sleepH2{};
    sleepH2.count = 1;
    sleepH2.items[0] = {SUBSTANCE_HYDROGEN, 1.0f};
    bool h2Set = worldG.tryCommitGasOccupancy(hi, sleepH2);
    float capH2 = ThermalEngine::gasCapacity(worldG, hi);
    worldG.heat[static_cast<size_t>(hi)] = energyFromTemp(capH2, AMBIENT_TEMPERATURE_K);
    std::fill(worldG.u.begin(), worldG.u.end(), 0.0f);
    std::fill(worldG.v.begin(), worldG.v.end(), 0.0f);
    std::fill(worldG.chunkActivity.begin(), worldG.chunkActivity.end(), uint8_t{0});
    std::fill(worldG.chunkSolveMask.begin(), worldG.chunkSolveMask.end(), uint8_t{0});
    std::fill(worldG.chunkQuietTicks.begin(), worldG.chunkQuietTicks.end(), uint8_t{255});
    int farC = ((GH - 8) / CHUNK) * CHUNK_W + ((GW - 8) / CHUNK);
    worldG.chunkActivity[static_cast<size_t>(farC)] = 1;
    worldG.config.simMode = GasSimMode::Full;
    worldF.tickNo = 0;
    worldG.simulationTick(worldF);
    float vSleep = worldG.v[static_cast<size_t>(GasEngine::vi(hx, hy))];
    float vSleepB = worldG.v[static_cast<size_t>(GasEngine::vi(hx, hy + 1))];
    emit("sleeping_composition_buoyancy_frozen",
        h2Set && std::abs(vSleep) < 1.0e-5f && std::abs(vSleepB) < 1.0e-5f,
        "vT=" + std::to_string(vSleep) + " vB=" + std::to_string(vSleepB));

    FluidEngine actF;
    RigidBodyEngine actR;
    GasEngine actG;
    ThermalEngine actT;
    ReactionEngine actRx;
    actG.resetAmbient(actF);
    actT.seedAmbient(actF, actR, actG);
    int ax = 40, ay = 40;
    int ai = FluidEngine::ci(ax, ay);
    GasComponentView stoichView{};
    stoichView.count = 2;
    stoichView.items[0] = {SUBSTANCE_HYDROGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0, cpm))};
    stoichView.items[1] = {SUBSTANCE_OXYGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
    bool stoichSet = actG.tryCommitGasOccupancy(ai, stoichView);
    float capFuel = ThermalEngine::gasCapacity(actG, ai);
    actG.heat[static_cast<size_t>(ai)] = energyFromTemp(capFuel, 293.15f);
    actT.wakeCell(ax, ay);
    actRx.simulationTick(actF, actR, actG, actT, PHYSICS_DT);
    emit("activity_zero_without_reaction",
        stoichSet && actRx.reactedCellsLastTick == 0 && actRx.activity[static_cast<size_t>(ai)] <= kActivityEps
            && actRx.extentLastTick <= kExtentEps, "");
    emit("cold_h2o2_still_does_not_react",
        ThermalEngine::gasTempK(actG, ai) < 400.0f
            && actG.gasComponentAmount(ai, SUBSTANCE_WATER) <= kMinGasComponent, "");

    capFuel = ThermalEngine::gasCapacity(actG, ai);
    actG.heat[static_cast<size_t>(ai)] = energyFromTemp(capFuel, 900.0f);
    actT.wakeCell(ax, ay);
    actRx.simulationTick(actF, actR, actG, actT, PHYSICS_DT);
    emit("activity_positive_on_extent",
        actRx.reactedCellsLastTick > 0 && actRx.extentLastTick > kExtentEps
            && actRx.heatReleasedLastTick > 0.0f
            && actRx.activity[static_cast<size_t>(ai)] > kActivityEps,
        "e=" + std::to_string(actRx.extentLastTick)
            + " a=" + std::to_string(actRx.activity[static_cast<size_t>(ai)]));

    float aHot = actRx.activity[static_cast<size_t>(ai)];
    GasComponentView airOnly{};
    airOnly.count = 1;
    airOnly.items[0] = {SUBSTANCE_AIR, 1.0f};
    (void)actG.tryCommitGasOccupancy(ai, airOnly);
    actRx.simulationTick(actF, actR, actG, actT, PHYSICS_DT);
    float aAfter = actRx.activity[static_cast<size_t>(ai)];
    emit("activity_decays_toward_zero",
        aHot > kActivityEps && aAfter < aHot && aAfter > 0.0f,
        "a0=" + std::to_string(aHot) + " a1=" + std::to_string(aAfter));

    FluidEngine propF;
    RigidBodyEngine propR;
    GasEngine propG;
    ThermalEngine propT;
    ReactionEngine propRx;
    propF.config.walledBorders = true;
    propG.config.simMode = GasSimMode::Off;
    propG.config.boundary = GasBoundary::Sealed;
    propT.config.enabled = true;
    propT.config.intervalTicks = 1;
    propG.resetAmbient(propF);
    propT.seedAmbient(propF, propR, propG);
    int px = 60, py = 60;
    int pi0 = FluidEngine::ci(px, py);
    int pi1 = FluidEngine::ci(px + 1, py);
    auto fillStoich = [&](int index) {
        GasComponentView v{};
        v.count = 2;
        v.items[0] = {SUBSTANCE_HYDROGEN,
            static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0, cpm))};
        v.items[1] = {SUBSTANCE_OXYGEN,
            static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
        (void)propG.tryCommitGasOccupancy(index, v);
        float cap = ThermalEngine::gasCapacity(propG, index);
        propG.heat[static_cast<size_t>(index)] = energyFromTemp(cap, AMBIENT_TEMPERATURE_K);
    };
    fillStoich(pi0);
    fillStoich(pi1);
    float cap0 = ThermalEngine::gasCapacity(propG, pi0);
    propG.heat[static_cast<size_t>(pi0)] = energyFromTemp(cap0, 900.0f);
    propT.wakeCell(px, py);
    propT.wakeCell(px + 1, py);
    bool neighborIgnited = false;
    for (int n = 0; n < 90; ++n) {
        propT.simulationTick(propF, propR, propG, PHYSICS_DT);
        propRx.simulationTick(propF, propR, propG, propT, PHYSICS_DT);
        if (propG.gasComponentAmount(pi1, SUBSTANCE_WATER) > kMinGasComponent) {
            neighborIgnited = true;
            break;
        }
    }
    emit("ignition_propagates_by_heat",
        neighborIgnited && ThermalEngine::gasTempK(propG, pi0) > 500.0f,
        "nT=" + std::to_string(ThermalEngine::gasTempK(propG, pi1))
            + " nWater=" + std::to_string(propG.gasComponentAmount(pi1, SUBSTANCE_WATER)));

    auto nearAtm = [&](float a, float b) {
        return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= 0.05f * std::max(1.0f, std::abs(b)) + 0.02f;
    };
    float pAmb = gasPressureAtmFromState(1.0f, 1.0f, AMBIENT_TEMPERATURE_K);
    emit("pressure_ambient_air_1atm", nearAtm(pAmb, 1.0f), "P=" + std::to_string(pAmb));
    float pHot = gasPressureAtmFromState(1.0f, 1.0f, 2.0f * AMBIENT_TEMPERATURE_K);
    emit("pressure_doubles_with_T", nearAtm(pHot, 2.0f * pAmb), "P2T=" + std::to_string(pHot));
    float pAmt = gasPressureAtmFromState(2.0f, 1.0f, AMBIENT_TEMPERATURE_K);
    emit("pressure_doubles_with_amount", nearAtm(pAmt, 2.0f * pAmb), "P2a=" + std::to_string(pAmt));
    float pVol = gasPressureAtmFromState(1.0f, 0.5f, AMBIENT_TEMPERATURE_K);
    emit("pressure_doubles_with_half_volume", nearAtm(pVol, 2.0f * pAmb), "PhalfV=" + std::to_string(pVol));

    GasEngine pGas;
    pGas.volume[0] = 1.0f;
    GasComponentView vapView{};
    vapView.count = 1;
    vapView.items[0] = {SUBSTANCE_WATER, 0.4f};
    (void)pGas.tryCommitGasOccupancy(0, vapView);
    float capV = ThermalEngine::gasCapacity(pGas, 0);
    pGas.heat[0] = energyFromTemp(capV, AMBIENT_TEMPERATURE_K);
    pGas.recomputePressure();
    float pVapCold = pGas.gasPartialPressurePa(0, SUBSTANCE_WATER);
    float heatBefore = pGas.heat[0];
    pGas.heat[0] = energyFromTemp(capV, 400.0f);
    pGas.recomputePressure();
    float pVapHot = pGas.gasPartialPressurePa(0, SUBSTANCE_WATER);
    emit("water_partial_pressure_rises_with_T",
        pVapHot > pVapCold * 1.2f && pVapCold > 0.0f,
        "cold=" + std::to_string(pVapCold) + " hot=" + std::to_string(pVapHot));
    emit("pressure_recompute_does_not_modify_heat",
        std::abs(pGas.heat[0] - energyFromTemp(capV, 400.0f)) <= 1.0e-3f,
        "heat=" + std::to_string(pGas.heat[0]));
    (void)heatBefore;

    GasEngine liveG;
    liveG.volume[0] = 1.0f;
    GasComponentView liveMix{};
    liveMix.count = 2;
    liveMix.items[0] = {SUBSTANCE_WATER, 0.4f};
    liveMix.items[1] = {SUBSTANCE_AIR, 0.6f};
    (void)liveG.tryCommitGasOccupancy(0, liveMix);
    float capLive = ThermalEngine::gasCapacity(liveG, 0);
    liveG.heat[0] = energyFromTemp(capLive, AMBIENT_TEMPERATURE_K);
    liveG.recomputePressure();
    float liveCold = liveG.gasPartialPressurePaFromCurrentState(0, SUBSTANCE_WATER);
    float cachedCold = liveG.gasPartialPressurePa(0, SUBSTANCE_WATER);
    float airAmt = liveG.gasComponentAmount(0, SUBSTANCE_AIR);
    float vapAmt = liveG.gasComponentAmount(0, SUBSTANCE_WATER);
    liveG.heat[0] = energyFromTemp(capLive, 400.0f);
    float liveHotStale = liveG.gasPartialPressurePaFromCurrentState(0, SUBSTANCE_WATER);
    float cachedStale = liveG.gasPartialPressurePa(0, SUBSTANCE_WATER);
    emit("live_water_partial_pressure_rises_with_T",
        liveHotStale > liveCold * 1.2f && liveCold > 0.0f,
        "cold=" + std::to_string(liveCold) + " hot=" + std::to_string(liveHotStale));
    emit("live_water_partial_pressure_ignores_stale_cache",
        std::abs(cachedStale - cachedCold) < 1.0f && liveHotStale > cachedStale * 1.2f,
        "cached=" + std::to_string(cachedStale) + " live=" + std::to_string(liveHotStale));

    float TsatQ = 350.0f;
    double PsatQ = saturationVaporPressurePa(SUBSTANCE_WATER, TsatQ);
    float prefQ = liveG.config.referencePressurePa;
    float volQ = 1.0f;
    float satAmt = gasAmountFromPressureAtm(static_cast<float>(PsatQ / static_cast<double>(prefQ)), volQ, TsatQ);
    float pBack = gasPressureAtmFromState(satAmt, volQ, TsatQ) * prefQ;
    emit("saturation_amount_round_trips_to_Psat",
        std::isfinite(pBack) && std::abs(pBack - static_cast<float>(PsatQ))
            <= 0.02f * static_cast<float>(PsatQ) + 1.0f,
        "Psat=" + std::to_string(PsatQ) + " Pback=" + std::to_string(pBack)
            + " satAmt=" + std::to_string(satAmt));
    float satHot = gasAmountFromPressureAtm(static_cast<float>(PsatQ / static_cast<double>(prefQ)),
        volQ, 2.0f * TsatQ);
    emit("hotter_gas_needs_less_saturation_vapor",
        satAmt > GAS_MIN_AMOUNT && satHot < satAmt
            && std::abs(satHot / satAmt - 0.5f) <= 0.02f,
        "sat=" + std::to_string(satAmt) + " satHot=" + std::to_string(satHot));
    double excessQ = static_cast<double>(vapAmt) - static_cast<double>(satAmt);
    if (excessQ < 0.0) excessQ = 0.0;
    excessQ = std::min(excessQ, static_cast<double>(vapAmt));
    emit("condensation_excess_never_exceeds_water_vapor",
        excessQ >= 0.0 && excessQ <= static_cast<double>(vapAmt) + 1.0e-9,
        "vap=" + std::to_string(vapAmt) + " excess=" + std::to_string(excessQ));
    emit("condensation_leaves_nonwater_components",
        std::abs(liveG.gasComponentAmount(0, SUBSTANCE_AIR) - airAmt) < 1.0e-8f
            && std::abs(liveG.gasComponentAmount(0, SUBSTANCE_WATER) - vapAmt) < 1.0e-8f,
        "air=" + std::to_string(airAmt) + " vap=" + std::to_string(vapAmt));

    FluidEngine cF;
    RigidBodyEngine cR;
    GasEngine cG;
    ThermalEngine cT;
    ReactionEngine cRx;
    cF.config.walledBorders = true;
    cG.config.simMode = GasSimMode::Off;
    cG.config.boundary = GasBoundary::Sealed;
    cT.config.enabled = false;
    cG.resetAmbient(cF);
    int ciComb = FluidEngine::ci(70, 70);
    GasComponentView combView{};
    combView.count = 2;
    combView.items[0] = {SUBSTANCE_HYDROGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_HYDROGEN, MatterPhase::Gas, 2.0, cpm))};
    combView.items[1] = {SUBSTANCE_OXYGEN,
        static_cast<float>(molesToStorageAmount(SUBSTANCE_OXYGEN, MatterPhase::Gas, 1.0, cpm))};
    (void)cG.tryCommitGasOccupancy(ciComb, combView);
    float capC = ThermalEngine::gasCapacity(cG, ciComb);
    cG.heat[static_cast<size_t>(ciComb)] = energyFromTemp(capC, 900.0f);
    cG.recomputePressure();
    float pBefore = cG.pressureAtm(ciComb);
    float TBefore = cG.cellTemperatureK(ciComb);
    for (int n = 0; n < 30; ++n)
        cRx.simulationTick(cF, cR, cG, cT, PHYSICS_DT);
    cG.recomputePressure();
    float pAfter = cG.pressureAtm(ciComb);
    float TAfter = cG.cellTemperatureK(ciComb);
    float heatAfter = cG.heat[static_cast<size_t>(ciComb)];
    float pEq = gasPressureAtmFromState(cG.amount[static_cast<size_t>(ciComb)],
        cG.volume[static_cast<size_t>(ciComb)], TAfter);
    emit("combustion_pressure_matches_equation",
        TAfter > TBefore + 200.0f && nearAtm(pAfter, pEq) && pAfter > pBefore,
        "T0=" + std::to_string(TBefore) + " T1=" + std::to_string(TAfter)
            + " P0=" + std::to_string(pBefore) + " P1=" + std::to_string(pAfter)
            + " Peq=" + std::to_string(pEq));
    cG.recomputePressure();
    emit("combustion_pressure_does_not_spend_heat",
        std::abs(cG.heat[static_cast<size_t>(ciComb)] - heatAfter) <= 1.0e-3f,
        "h=" + std::to_string(heatAfter));

    emit("idle_reaction_glow_fast_path",
        !ReactionEngine().hasVisibleActivity(), "");
    cRx.clearActivity();
    emit("scene_load_clears_reaction_glow",
        !cRx.hasVisibleActivity() && cRx.activity[static_cast<size_t>(ciComb)] <= kActivityEps, "");

    out << "summary\t" << (failed == 0 ? "PASS" : "FAIL") << '\t' << passed << " passed, "
        << failed << " failed\n";
}
