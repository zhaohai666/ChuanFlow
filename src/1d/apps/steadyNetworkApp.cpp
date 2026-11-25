#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/1d/Nodal.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

using namespace cf;

namespace
{
// log::info() stringifies with ostringstream defaults (6 digits); results
// are checked against closed-form values, so emit 12 significant digits.
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

const char* regimeName(int r)
{
    static const char* n[] = {"segregated", "transition", "intermittent",
                              "distributed"};
    return r >= 0 && r <= 3 ? n[r] : "none";
}

// The opt-in composition-generated PVT table echo needs round-trip precision:
// a checker rebuilds an explicit-table twin from these numbers and requires
// the solver readout to be byte-identical. 17 significant digits is a full
// double round-trip; only the 5.110 echo uses it.
std::string num17(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

// ---------------------------------------------------------------------------
// 5.160 (P21): the two opt-in steady-network optimization blocks. Both run
// AFTER the base solve and its report and print only their own lines, so a
// case that declares neither is byte-identical to every existing case. The
// library parser ignores unknown top-level keys, so the blocks are read here
// in the app layer alone: zero kernel physics moved.
// ---------------------------------------------------------------------------

std::size_t nodePosOf(const OneDNetwork& net, const std::string& name,
                      const std::string& who)
{
    for (std::size_t i = 0; i < net.nodeNames.size(); ++i)
        if (net.nodeNames[i] == name) return i;
    cfFatal("steadyNetwork: " + who + " node '" + name
            + "' is not in the network");
}

std::size_t pipePosOf(const OneDNetwork& net, const std::string& name,
                      const std::string& who)
{
    for (std::size_t i = 0; i < net.pipes.size(); ++i)
        if (net.pipes[i].name == name) return i;
    cfFatal("steadyNetwork: " + who + " pipe '" + name
            + "' is not in the network");
}

// The Dictionary parser stores every key it reads; a whitelist here is the
// only channel that turns an inert declaration into a refusal.
void rejectUnknownKeys(const Dictionary& d,
                       const std::vector<std::string>& allowed,
                       const std::string& who)
{
    for (const std::string& k : d.keys())
        if (std::find(allowed.begin(), allowed.end(), k) == allowed.end())
            cfFatal("steadyNetwork: " + who + " unknown key '" + k + "'");
}

struct DiaCombo
{
    std::vector<double> dia;
    double objective = 0; // the natural (user-facing) value
    double key = 0;       // the internal minimized value
    double cost = 0;
};

void runDiameterOptimization(const Dictionary& od, const OneDNetwork& net,
                             const SteadyControls& ctl)
{
    rejectUnknownKeys(od, {"pipes", "candidates", "objective", "objectiveNode",
                           "objectiveSense", "maxCombinations", "reportTop",
                           "costPerMeter", "costTable", "constraints"},
                      "diameterOptimization{}:");

    const std::vector<std::string> pipeNames = od.getWordList("pipes");
    if (pipeNames.empty())
        cfFatal("steadyNetwork: diameterOptimization{} declares no pipes");
    std::vector<std::size_t> pipePos;
    for (const std::string& pn : pipeNames)
        pipePos.push_back(pipePosOf(net, pn, "diameterOptimization{}:"));

    const std::vector<double> cand = od.getScalarList("candidates");
    for (std::size_t i = 0; i < cand.size(); ++i)
    {
        if (!(cand[i] > 0.0))
            cfFatal("steadyNetwork: diameterOptimization candidate diameters "
                    "must be positive");
        if (i > 0 && !(cand[i] > cand[i - 1]))
            cfFatal("steadyNetwork: diameterOptimization candidates must be "
                    "strictly ascending");
    }

    const std::string objective = od.getWord("objective");
    int objMode = -1; // 0 minTotalDp, 1 nodePressure, 2 minCost
    std::size_t objNode = 0;
    if (objective == "minTotalDp") objMode = 0;
    else if (objective == "nodePressure")
    {
        objMode = 1;
        if (!od.has("objectiveNode"))
            cfFatal("steadyNetwork: diameterOptimization objective "
                    "nodePressure needs objectiveNode");
        objNode = nodePosOf(net, od.getWord("objectiveNode"),
                            "diameterOptimization objectiveNode:");
    }
    else if (objective == "minCost") objMode = 2;
    else
        cfFatal("steadyNetwork: diameterOptimization unknown objective '"
                + objective + "' (minTotalDp, nodePressure or minCost)");

    const std::string sense = od.getWordOr("objectiveSense", "min");
    if (sense != "min" && sense != "max")
        cfFatal("steadyNetwork: diameterOptimization objectiveSense must be "
                "min or max");

    const bool hasCpm = od.has("costPerMeter");
    const bool hasTable = od.has("costTable");
    if (hasCpm && hasTable)
        cfFatal("steadyNetwork: diameterOptimization cannot declare both "
                "costPerMeter and costTable");
    double cpm = 0.0;
    std::vector<double> costTab;
    if (hasCpm)
    {
        cpm = od.getScalar("costPerMeter");
        if (!(cpm > 0.0))
            cfFatal("steadyNetwork: diameterOptimization costPerMeter must "
                    "be positive");
    }
    if (hasTable)
    {
        costTab = od.getScalarList("costTable");
        if (costTab.size() != cand.size())
            cfFatal("steadyNetwork: diameterOptimization costTable must have "
                    "one entry per candidate");
        for (double c : costTab)
            if (!(c >= 0.0))
                cfFatal("steadyNetwork: diameterOptimization costTable "
                        "entries must be non-negative");
    }
    if (objMode == 2 && !hasCpm && !hasTable)
        cfFatal("steadyNetwork: diameterOptimization objective minCost needs "
                "costPerMeter or costTable");
    const bool costOn = hasCpm || hasTable;

    const double mcd = od.getScalarOr("maxCombinations", 4096.0);
    if (!(mcd >= 1.0) || mcd != std::floor(mcd))
        cfFatal("steadyNetwork: diameterOptimization maxCombinations must be "
                "an integer >= 1");
    const std::size_t maxComb = static_cast<std::size_t>(mcd);
    const double rtd = od.getScalarOr("reportTop", 3.0);
    if (!(rtd >= 0.0) || rtd != std::floor(rtd))
        cfFatal("steadyNetwork: diameterOptimization reportTop must be an "
                "integer >= 0");
    const std::size_t reportTop = static_cast<std::size_t>(rtd);

    struct Con
    {
        std::size_t node;
        double lo, hi;
    };
    std::vector<Con> cons;
    if (od.has("constraints"))
    {
        const Dictionary& cd = od.getDict("constraints");
        for (const std::string& cn : cd.keys())
        {
            const Dictionary& x = cd.getDict(cn);
            rejectUnknownKeys(x, {"node", "minPressure", "maxPressure"},
                              "diameterOptimization constraint '" + cn + "':");
            if (!x.has("node"))
                cfFatal("steadyNetwork: diameterOptimization constraint '"
                        + cn + "' has no node");
            Con c;
            c.node = nodePosOf(net, x.getWord("node"),
                               "diameterOptimization constraint node:");
            const bool lo = x.has("minPressure");
            const bool hi = x.has("maxPressure");
            if (!lo && !hi)
                cfFatal("steadyNetwork: diameterOptimization constraint '"
                        + cn + "' declares neither minPressure nor "
                        "maxPressure");
            c.lo = x.getScalarOr("minPressure",
                                 -std::numeric_limits<double>::infinity());
            c.hi = x.getScalarOr("maxPressure",
                                 std::numeric_limits<double>::infinity());
            if (c.lo > c.hi)
                cfFatal("steadyNetwork: diameterOptimization constraint '"
                        + cn + "' has minPressure above maxPressure");
            cons.push_back(c);
        }
    }

    // overflow-safe count: pipe 0 is the slowest digit, the last the fastest
    std::size_t count = 1;
    for (std::size_t j = 0; j < pipePos.size(); ++j)
    {
        if (count > maxComb / cand.size())
            cfFatal("steadyNetwork: diameterOptimization enumeration exceeds "
                    "maxCombinations");
        count *= cand.size();
    }

    std::vector<DiaCombo> feas;
    for (std::size_t idx = 0; idx < count; ++idx)
    {
        OneDNetwork trial = net;
        std::vector<double> dia(pipePos.size());
        std::vector<std::size_t> dig(pipePos.size());
        std::size_t rest = idx;
        for (std::size_t jj = pipePos.size(); jj-- > 0;)
        {
            dig[jj] = rest % cand.size();
            rest /= cand.size();
            dia[jj] = cand[dig[jj]];
            trial.pipes[pipePos[jj]].diameter = cand[dig[jj]];
        }
        const SteadyNetworkResult r2 = solveSteadyNetwork(trial, ctl);
        if (!r2.converged) continue;
        bool ok = true;
        for (const Con& c : cons)
            if (!(r2.nodePressure[c.node] >= c.lo
                  && r2.nodePressure[c.node] <= c.hi))
            {
                ok = false;
                break;
            }
        if (!ok) continue;

        DiaCombo row;
        row.dia = dia;
        if (costOn)
        {
            for (std::size_t jj = 0; jj < pipePos.size(); ++jj)
            {
                const double L = net.pipes[pipePos[jj]].length;
                if (hasTable) row.cost += L * costTab[dig[jj]];
                else
                {
                    const double a = dia[jj] / cand[0];
                    row.cost += cpm * L * a * a;
                }
            }
        }
        if (objMode == 0)
        {
            for (const PipeFlowResult& f : r2.pipeFlows)
                row.objective += f.dpTotal;
            row.key = row.objective;
        }
        else if (objMode == 1)
        {
            row.objective = r2.nodePressure[objNode];
            row.key = sense == "min" ? row.objective : -row.objective;
        }
        else
        {
            row.objective = row.cost;
            row.key = row.cost;
        }
        feas.push_back(row);
    }

    // strict key order, ties keep the enumeration order: the first combo
    // enumerated wins a tie
    std::stable_sort(feas.begin(), feas.end(),
                     [](const DiaCombo& a, const DiaCombo& b)
                     { return a.key < b.key; });

    log::info("diameterOpt combos ", count, " feasible ", feas.size());
    auto printRow = [&](const std::string& tag, const DiaCombo& row)
    {
        std::string line = "diameterOpt " + tag;
        for (std::size_t jj = 0; jj < pipePos.size(); ++jj)
            line += " " + pipeNames[jj] + " " + num(row.dia[jj]);
        line += " objective " + num(row.objective);
        if (costOn) line += " cost " + num(row.cost);
        log::info(line);
    };
    if (feas.empty()) log::info("diameterOpt best NONE");
    else
    {
        printRow("best", feas[0]);
        const std::size_t top = std::min(reportTop, feas.size());
        for (std::size_t k = 1; k < top; ++k)
            printRow("runner " + std::to_string(k + 1), feas[k]);
    }
}

// The gas-lift rate curve of one well at injection gas g (Sm3/s): the gas
// joins the stream at bottom hole, so the effective producing GOR is
// gorStd + g/q, and q is the nodal operating point marched with that GOR -
// a fixed point in q. g = 0 is the base nodal solution itself.
double gasLiftRateAtGas(const NodalCase& nc, double g, double q0)
{
    if (g <= 0.0) return q0;
    double q = q0;
    for (int it = 0; it < 60; ++it)
    {
        pvt::BlackOilFluid f2 = nc.fluid;
        f2.gorStd = nc.fluid.gorStd + g / q;
        const NodalResult r = solveNodalAnalysis(nc.ipr, f2, nc.well,
                                                 nc.wellheadPressure,
                                                 nc.controls);
        if (!r.converged)
            cfFatal("steadyNetwork: gasLiftAllocation nodal solve did not "
                    "converge");
        if (!(r.rate > 0.0))
            cfFatal("steadyNetwork: gasLiftAllocation well went to zero rate "
                    "under injection gas");
        if (std::abs(r.rate - q) <= 1e-10 * r.rate) return r.rate;
        q = r.rate;
    }
    cfFatal("steadyNetwork: gasLiftAllocation the gas/GOR fixed point did "
            "not converge");
}

void runGasLiftAllocation(const Dictionary& gd, const Dictionary& caseDict,
                          const OneDNetwork& net,
                          const SteadyNetworkResult& res)
{
    rejectUnknownKeys(gd, {"budget", "resolution", "wells"},
                      "gasLiftAllocation{}:");
    if (!net.hasBlackOil)
        cfFatal("steadyNetwork: gasLiftAllocation{} needs a black-oil "
                "network: the lift gas rides the black-oil GOR");
    const double budget = gd.getScalar("budget");
    if (!(budget > 0.0))
        cfFatal("steadyNetwork: gasLiftAllocation budget must be positive");
    const double rd = gd.getScalarOr("resolution", 8.0);
    if (rd != std::floor(rd) || rd < 1.0 || rd > 64.0)
        cfFatal("steadyNetwork: gasLiftAllocation resolution must be an "
                "integer in [1, 64]");
    const int R = static_cast<int>(rd);
    if (!gd.has("wells"))
        cfFatal("steadyNetwork: gasLiftAllocation{} has no wells{} block");
    const Dictionary& wd = gd.getDict("wells");
    if (wd.keys().empty())
        cfFatal("steadyNetwork: gasLiftAllocation wells{} declares no wells");

    const double step = budget / R;
    struct WellCurve
    {
        std::string name;
        std::vector<double> q; // q[k] at gas k*step
    };
    std::vector<WellCurve> curves;
    for (const std::string& name : wd.keys())
    {
        const Dictionary& x = wd.getDict(name);
        rejectUnknownKeys(x, {"node", "ipr", "well"},
                          "gasLiftAllocation well '" + name + "':");
        if (!x.has("node"))
            cfFatal("steadyNetwork: gasLiftAllocation well '" + name
                    + "' has no node");
        const std::size_t ni = nodePosOf(net, x.getWord("node"),
                                         "gasLiftAllocation node:");
        // Reuse the nodal parser verbatim: hand it a sub-dictionary carrying
        // the case's own fluid{}, the well's ipr{} and well{}, and the base
        // network's node pressure as the wellhead back-pressure.
        Dictionary sub;
        sub.insert("fluid", caseDict.find("fluid"));
        sub.insert("ipr", x.find("ipr"));
        sub.insert("well", x.find("well"));
        Dictionary::Node pwh;
        pwh.kind = Dictionary::Kind::scalar;
        pwh.scalar = res.nodePressure[ni];
        sub.insert("wellheadPressure", std::move(pwh));
        const NodalCase nc = loadNodalCaseFromDict(sub);
        const NodalResult r0 = solveNodalAnalysis(nc.ipr, nc.fluid, nc.well,
                                                  nc.wellheadPressure,
                                                  nc.controls);
        if (!r0.converged)
            cfFatal("steadyNetwork: gasLiftAllocation well '" + name
                    + "' nodal solve did not converge");
        if (!(r0.rate > 0.0))
            cfFatal("steadyNetwork: gasLiftAllocation well '" + name
                    + "' has zero rate at zero gas: the allocation needs a "
                    "flowing well");
        WellCurve wc;
        wc.name = name;
        wc.q.push_back(r0.rate);
        for (int k = 1; k <= R; ++k)
            wc.q.push_back(gasLiftRateAtGas(nc, k * step, r0.rate));
        curves.push_back(std::move(wc));
    }

    // multiple-choice knapsack on the integer grid: dp[b] is the best total
    // rate the first i wells can produce from exactly b steps; the strict >
    // keeps the first argmax, so ties resolve to the earlier well.
    const int nW = static_cast<int>(curves.size());
    const double NEG = -std::numeric_limits<double>::infinity();
    std::vector<double> dp(R + 1, NEG), ndp(R + 1, NEG);
    dp[0] = 0.0;
    std::vector<std::vector<int>> parent(nW, std::vector<int>(R + 1, 0));
    for (int i = 0; i < nW; ++i)
    {
        for (int b = 0; b <= R; ++b)
        {
            double best = NEG;
            int bk = 0;
            for (int k = 0; k <= b; ++k)
            {
                if (dp[b - k] == NEG) continue;
                const double v = dp[b - k] + curves[i].q[k];
                if (v > best) { best = v; bk = k; }
            }
            ndp[b] = best;
            parent[i][b] = bk;
        }
        dp = ndp;
        std::fill(ndp.begin(), ndp.end(), NEG);
    }
    std::vector<int> alloc(nW, 0);
    int b = R;
    for (int i = nW - 1; i >= 0; --i)
    {
        alloc[i] = parent[i][b];
        b -= alloc[i];
    }

    std::vector<double> marg(nW, 0.0);
    for (int i = 0; i < nW; ++i)
    {
        const std::vector<double>& q = curves[i].q;
        if (alloc[i] == 0) marg[i] = (q[1] - q[0]) / step;
        else if (alloc[i] == R) marg[i] = (q[R] - q[R - 1]) / step;
        else marg[i] = (q[alloc[i] + 1] - q[alloc[i] - 1]) / (2.0 * step);
    }

    log::info("gasLift wells ", nW, " budget ", num(budget),
              " step ", num(step));
    double total = 0.0, baseTotal = 0.0;
    double mMax = NEG, mMin = -NEG;
    for (int i = 0; i < nW; ++i)
    {
        const double g = alloc[i] * step;
        const double q = curves[i].q[alloc[i]];
        log::info("gasLift well ", curves[i].name, " gas ", num(g),
                  " q ", num(q), " gain ", num(q - curves[i].q[0]),
                  " marginal ", num(marg[i]));
        total += q;
        baseTotal += curves[i].q[0];
        mMax = std::max(mMax, marg[i]);
        mMin = std::min(mMin, marg[i]);
    }
    log::info("gasLift total q ", num(total), " gain ", num(total - baseTotal));
    log::info("gasLift slopeSpread ", num(mMax - mMin));
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-1d-steady <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const OneDNetwork net = loadNetworkFromDict(caseDict);
        const SteadyControls ctl = loadSteadyControls(caseDict);

        // 5.110 (increment A): opt-in echo of the composition-generated PVT
        // table so a checker can verify rho(p)/mu(p) point by point against an
        // independent PR76 + LBC re-derivation. Absent the declaration, not a
        // single byte of output changes.
        if (caseDict.has("fluid"))
        {
            const Dictionary& fd = caseDict.getDict("fluid");
            if (fd.has("pvt") && fd.getDict("pvt").has("compositional")
                && !net.pvtFluid.liveFlash)
            {
                for (std::size_t i = 0; i < net.pvtFluid.pressure.size(); ++i)
                {
                    log::info("pvtTable p ", num17(net.pvtFluid.pressure[i]),
                              " rho ", num17(net.pvtFluid.rho[i]),
                              " mu ", num17(net.pvtFluid.mu[i]));
                }
            }
        }

        log::info("network: ", net.nNodes(), " nodes, ", net.pipes.size(), " pipes");

        const SteadyNetworkResult res = solveSteadyNetwork(net, ctl);

        for (std::size_t i = 0; i < net.nNodes(); ++i)
        {
            log::info("node ", net.nodeNames[i], " p ", num(res.nodePressure[i]));
        }
        std::vector<std::string> slugRiskPipes; // pipes with a SLUG-START dip
        bool anySlugDeclared = false;
        for (std::size_t pi = 0; pi < res.pipeFlows.size(); ++pi)
        {
            const PipeFlowResult& f = res.pipeFlows[pi];
            const NetworkPipe& pipe = net.pipes[pi];
            // D10/D7: per-station true-flash readout. Emitted only for a
            // live-flash pipe, so every existing case prints byte-identically.
            for (std::size_t i = 0; i < f.pvtFlashStations.size(); ++i)
            {
                const PvtFlashStation& s = f.pvtFlashStations[i];
                log::info("flashStation ", f.name, " i ", i,
                          " p ", num17(s.p), " rho ", num17(s.rho),
                          " mu ", num17(s.mu), " sigma ", num17(s.sigma));
            }
            std::string line = "pipe " + f.name + " Q " + num(f.flowRate)
                             + " v " + num(f.velocity);
            if (f.regime >= 0)
            {
                line += " HL " + num(f.holdup) + " regime "
                      + std::string(regimeName(f.regime))
                      + " m " + num(f.massRate);
                if (f.regimeOutlet >= 0)
                {
                    line += " oRegime "
                          + std::string(regimeName(f.regimeOutlet))
                          + " lambda " + num(f.lambdaOutlet)
                          + " Nfr " + num(f.nfrOutlet)
                          + " L1 " + num(f.patternL1)
                          + " L2 " + num(f.patternL2)
                          + " L3 " + num(f.patternL3)
                          + " L4 " + num(f.patternL4);
                    // 5.88 (D4 increment 2): the override's footprint. Printed
                    // only when a forced pattern actually disagreed with the
                    // map at the outlet segment, so an auto (or agreeing)
                    // declaration leaves the line byte-identical to 5.85.
                    if (f.regimeFromMapOutlet >= 0
                        && f.regimeFromMapOutlet != f.regimeOutlet)
                    {
                        line += " mapRegime "
                              + std::string(regimeName(f.regimeFromMapOutlet))
                              + " FORCED";
                    }
                }
                if (f.choked) line += " CHOKED";
            }
            else if (f.valveDisposition >= 0)
            {
                // 5.89 (D1 increment 3): a gas-lift valve has no wall law to
                // report a Reynolds number and a friction factor from. What it
                // does have is a force balance and a nozzle, so the line says
                // which side of its own thresholds the injection pressure came
                // in on, the two thresholds, the band between them and how far
                // the tubing pressure pulled the opening pressure off its
                // test-rack value, and the mass rate the port passed.
                static const char* dispositionName[] = {"CLOSED", "OPEN",
                                                        "IN-SPREAD"};
                line += " valve " + std::string(dispositionName[f.valveDisposition])
                      + " pvo " + num(f.valveOpeningPressure)
                      + " pvc " + num(f.valveClosingPressure)
                      + " spread " + num(f.valveSpread)
                      + " ppe " + num(f.valveProductionEffect)
                      + " tro " + num(f.valveTestRackOpening)
                      + " m " + num(f.massRate)
                      + (f.valveChoked ? " SONIC" : "");
            }
            else if (f.card.state >= 0)
            {
                // 5.90 (D1 increment 4): a rod pump row has no wall law either,
                // and what it does have is a card. The four loads say where the
                // polished rod and the plunger top out and bottom out, the
                // stretch and the effective stroke say what the fluid load cost
                // the string, and the rate, work and power say what came of it.
                // The parallelogram itself is not spelled out: it is fixed by
                // (stroke, stretch, wmin, wmax), and the checker draws it.
                static const char* cardStateName[] = {"FULL", "NO-LOAD",
                                                      "STALLED"};
                line += " rod " + std::string(cardStateName[f.card.state])
                      + " wmax " + num(f.card.maxLoad)
                      + " wmin " + num(f.card.minLoad)
                      + " dmax " + num(f.card.downholeMaxLoad)
                      + " dmin " + num(f.card.downholeMinLoad)
                      + " wr " + num(f.card.rodWeightInFluid)
                      + " fo " + num(f.card.fluidLoad)
                      + " fomax " + num(f.card.fluidLoadMax)
                      + " upthrust " + num(f.card.upthrust)
                      + " stretch " + num(f.card.stretch)
                      + " stroke " + num(f.card.stroke)
                      + " seff " + num(f.card.effectiveStroke)
                      + " ev " + num(f.card.efficiency)
                      + " disp " + num(f.card.displacement)
                      + " m " + num(f.card.massRate)
                      + " work " + num(f.card.workPerCycle)
                      + " power " + num(f.card.hydraulicPower)
                      + " cbe " + num(f.card.counterbalance);
            }
            else
            {
                line += " Re " + num(f.reynolds) + " f "
                      + num(f.frictionFactor);
            }
            line += " dp " + num(f.dpTotal);
            // 5.144 (G6): the passive thermal march's outlet temperature;
            // absent (0) on every isothermal path, so no existing line
            // gains a token.
            if (f.temperature > 0)
                line += " T " + num(f.temperature);
            if (!f.pvtFlashStations.empty())
                line += " sigmaOutlet " + num17(f.sigmaOutlet);
            if (f.dpValve != 0) line += " dpValve " + num(f.dpValve);
            if (f.dpPump != 0) line += " dpPump " + num(f.dpPump);
            // 5.191: the PCP's own rise, on the machine that priced it.
            if (f.dpPcp != 0) line += " dpPcp " + num(f.dpPcp);
            // 5.156 (P10): the compressor's passive readouts. Absent (ratio
            // 1) on every other path, so no existing line gains a token.
            if (f.compressorRatio != 1)
            {
                line += " cmpRatio " + num(f.compressorRatio)
                      + " cmpBoost " + num(f.compressorBoost)
                      + " cmpPower " + num(f.compressorPower)
                      + " cmpT2 " + num(f.compressorTemperatureOut);
            }
            if (f.transition.band >= 0)
            {
                line += " He " + num(f.transition.hedstrom)
                      + " ReB " + num(f.transition.rePlastic)
                      + " ReBt " + num(f.transition.reTransition)
                      + " Qt " + num(f.transition.rateAtTransition)
                      + " dpQt " + num(f.transition.dpAtTransition)
                      + " Vrelt " + num(f.transition.relTransitionSpeed)
                      + " band "
                      + num(static_cast<double>(f.transition.band))
                      + (f.transition.pastLaminarLimit ? " pastTransition"
                                                       : " belowTransition");
            }
            log::info(line);
            // 5.158 (P18): the slug-start forecast, printed only for a pipe
            // that declared slugPrediction{}; an undeclared case gains not a
            // single byte.
            if (pipe.slugDeclared)
            {
                anySlugDeclared = true;
                const SlugReadout sl = slugReadoutOf(pipe, f);
                bool anyStart = false;
                for (std::size_t k = 0; k < sl.dips.size(); ++k)
                {
                    const SlugDipReadout& d = sl.dips[k];
                    if (d.slugStart) anyStart = true;
                    log::info("slug ", f.name, " dip ", k, " s ", num(d.s),
                              " vol ", num(d.volume), " criterion ",
                              num(d.criterion), " regime ", regimeName(d.regime),
                              d.slugStart ? " SLUG-START" : " NO-SLUG");
                }
                log::info("slug ", f.name, " dips ", sl.dips.size());
                if (anyStart) slugRiskPipes.push_back(f.name);
            }
            // 5.159 (P22): the erosion / hydrate / corrosion risk lines, on a
            // pipe that declared riskReport{} only.
            if (pipe.riskDeclared)
            {
                const RiskReadout rk = riskReadoutOf(pipe, net.blackOil, f);
                if (rk.erosionOn)
                {
                    log::info("risk ", f.name, " erosion seg ", rk.erosionSeg,
                              " vm ", num(rk.erosionVm), " velim ",
                              num(rk.erosionVelim), " ratio ",
                              num(rk.erosionRatio),
                              rk.erosionRatio > 1.0 ? " OVER" : " OK",
                              " over ", rk.erosionOverCount);
                }
                if (rk.hydrateOn)
                {
                    static const char* gradeName[] = {"NONE", "LOW", "MEDIUM",
                                                      "HIGH"};
                    log::info("risk ", f.name, " hydrate seg ", rk.hydrateSeg,
                              " p ", num(rk.hydrateP), " T ", num(rk.hydrateT),
                              " Teq ", num(rk.hydrateTeq), " dTsub ",
                              num(rk.hydrateSubcool), " ",
                              gradeName[rk.hydrateGrade]);
                }
                if (rk.corrosionOn)
                {
                    static const char* co2BandName[] = {"LOW", "MEDIUM", "HIGH"};
                    log::info("risk ", f.name, " corrosion co2p ",
                              num(rk.co2Partial), " band ",
                              co2BandName[rk.co2Band], " h2sp ",
                              num(rk.h2sPartial), rk.sour ? " SOUR" : " SWEET");
                }
            }
            // 5.188 (G6b): the per-segment P/T profile of a thermal pipe, read
            // off the march that produced the pipe line above. One point per
            // segment inlet (x measured from the pipe inlet); the pipe line's
            // own `T` token is the outlet. Neither the anchors nor this
            // profile is a new physics - the hydraulics never read it.
            if (pipe.thermalOn()
                && f.marchSegX.size() == f.marchSegTemperature.size()
                && f.marchSegX.size() == f.marchSegPressure.size())
            {
                for (std::size_t k = 0; k < f.marchSegX.size(); ++k)
                    log::info("thermal ", f.name, " pt ", k, " x ",
                              num(f.marchSegX[k]), " p ",
                              num(f.marchSegPressure[k]), " T ",
                              num(f.marchSegTemperature[k]));
                // 5.188: the march's own exit point, so the outlet is read
                // off the same integration rather than off the node pressure
                // (which the outer solve may still be a tolerance away from).
                log::info("thermal ", f.name, " pt ", f.marchSegX.size(), " x ",
                          num(pipe.length), " p ", num(f.marchOutletPressure),
                          " T ", num(f.temperature));
                log::info("thermal ", f.name, " u ", num(f.marchThermalU),
                          " rout ", num(pipe.insulationOuterRadius),
                          " pts ", f.marchSegX.size(), " Tout ",
                          num(f.temperature));
                // 5.189: when the coefficient came from the EOS, say which EOS
                // and print the entry-state mu the march used, so the number
                // is re-derivable from the case without a solver rerun.
                if (f.marchJtMode != 0)
                {
                    log::info("thermal ", f.name, " jt ",
                              f.marchJtMode == 2 ? "eosSrk" : "eos",
                              " mu ", num(f.marchJouleThomson));
                    // 5.189: the coefficient at each segment's own state, the
                    // same (T, p) the march used. The accumulated drop is the
                    // trapezoid of these numbers against the p grid, so a
                    // checker can re-derive T_out and the mu curve apart.
                    for (std::size_t k = 0; k < f.marchSegX.size(); ++k)
                        log::info("thermal ", f.name, " jtpt ", k, " p ",
                                  num(f.marchSegPressure[k]), " T ",
                                  num(f.marchSegTemperature[k]), " mu ",
                                  num(net.pvtFluid.eosJouleThomson(
                                      f.marchSegTemperature[k],
                                      f.marchSegPressure[k],
                                      pipe.heatCapacity,
                                      f.marchJtMode == 2)));
                }
            }
        }
        if (anySlugDeclared)
        {
            std::string line = "slug risk pipes";
            if (slugRiskPipes.empty()) line += " none";
            else for (const std::string& s : slugRiskPipes) line += " " + s;
            log::info(line);
        }

        log::info("steady solved: iterations=", res.iterations,
                  " residual=", num(res.relativeResidual),
                  res.converged ? " (converged)" : " (NOT converged)");

        // 5.160 (P21): the opt-in optimization reports. They append their own
        // lines after the base readout and never touch it, so an undeclared
        // case stays byte-identical.
        if (caseDict.has("diameterOptimization"))
            runDiameterOptimization(caseDict.getDict("diameterOptimization"),
                                    net, ctl);
        if (caseDict.has("gasLiftAllocation"))
            runGasLiftAllocation(caseDict.getDict("gasLiftAllocation"),
                                 caseDict, net, res);

        return res.converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("1d-steady: ", e.what());
        return 1;
    }
}
