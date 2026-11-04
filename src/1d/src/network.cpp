#include "cf/1d/Network.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

#include "cf/1d/BeggsBrill.h"

#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double pipeArea(const NetworkPipe& p)
{
    return 0.25 * M_PI * p.diameter * p.diameter;
}

/// Hagen-Poiseuille slope pi D^4 / (128 mu L): the maximum conductance of a
/// pipe (laminar limit), used as the fallback when the difference quotient
/// degenerates at dp ~ 0.
double laminarConductance(const NetworkPipe& p, const FluidProperties& fluid)
{
    return M_PI * std::pow(p.diameter, 4) / (128.0 * fluid.mu * p.length);
}

struct RawFlow
{
    double velocity = 0; // m/s, signed from -> to
    double reynolds = 0;
    double friction = 0;
};

/// Pump head at rate q (m of fluid): H(q) = h0 - (h0 - h1)(q/q1)^2 clamped at
/// 0; constant h0 when q1 = 0. Even in q: the machine adds the same head
/// magnitude whichever way the flow runs (never sucks).
double pumpHeadOf(const NetworkPipe& p, double q)
{
    if (!(p.pumpHead > 0)) return 0.0;
    const double h0 = p.pumpShutoffHead > 0 ? p.pumpShutoffHead : p.pumpHead;
    if (!(p.pumpRate > 0)) return std::max(h0, 0.0);
    const double h = h0 - (h0 - p.pumpHead) * q * q / (p.pumpRate * p.pumpRate);
    return std::max(h, 0.0);
}

/// Friction + valve loss (Pa, positive) at an imposed rate magnitude q:
/// explicit Churchill evaluation, no iteration needed forward.
double pipeLossAt(const NetworkPipe& pipe, const FluidProperties& fluid,
                  double q)
{
    const double v = q / pipeArea(pipe);
    const double re = fluid.rho * v * pipe.diameter / fluid.mu;
    const double f = churchillFriction(re, pipe.roughness / pipe.diameter);
    return (f * pipe.length / pipe.diameter + pipe.valveK)
         * fluid.rho * v * v / 2.0;
}

/// Pipe + equipment law Q(dpTotal): dpTotal = rho g dz + loss(q) - rho g H(q).
/// Solved as a bracketed Illinois secant on the rate. The pump-head fixed
/// point (hydraulics need the head, the curve needs the rate) is NOT
/// contractive whenever the curve falls faster than the system rises
/// (rho g |H'| > d(loss)/dq): it oscillates and latches onto the H = 0 clamp.
/// Root-solving is monotone on the branch picked by the sign of
/// G(0) = rho g dz - rho g h0 - dpTotal: loss grows without bound in |q|
/// while -rho g H(q) stays bounded, so exactly one crossing per branch.
RawFlow equipmentFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                      double dpTotal, double& Hq)
{
    const bool forward = fluid.rho * kGravity * pipe.elevation
                       - fluid.rho * kGravity * pumpHeadOf(pipe, 0.0)
                       - dpTotal <= 0.0;
    auto g = [&](double x) // residual along the chosen branch, |q| = x >= 0
    {
        const double s = forward ? 1.0 : -1.0;
        return s * pipeLossAt(pipe, fluid, x)
             + fluid.rho * kGravity * pipe.elevation
             - fluid.rho * kGravity * pumpHeadOf(pipe, x) - dpTotal;
    };

    // bracket: forward g(0) <= 0, expand to a positive sample; reverse
    // g(0) > 0, expand to a negative sample (keep the near-zero endpoint)
    const double pScale = std::max(std::abs(dpTotal),
                                   fluid.rho * kGravity
                                     * (std::abs(pipe.elevation)
                                        + pumpHeadOf(pipe, 0.0)));
    const double xSeed = std::max(1e-3 * pipeArea(pipe), 1e-8);
    double xLo = 0.0, fLo = g(0.0);
    double xHi = xSeed, fHi = g(xHi);
    int expansions = 0;
    while (!(forward ? fHi > 0.0 : fHi < 0.0))
    {
        if (forward ? fHi <= 0.0 : fHi >= 0.0)
        {
            xLo = xHi;
            fLo = fHi;
        }
        xHi *= 4.0;
        fHi = g(xHi);
        if (++expansions > 60) break; // loss dominates long before this
    }

    double x = xLo + (xHi - xLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fx = g(x);
        if (std::abs(fx) < 1e-13 * std::max(pScale, 1.0)) break;
        if (forward ? fx > 0.0 : fx < 0.0)
        {
            xHi = x;
            fHi = fx;
            // Illinois: bias the retained endpoint against re-selection
            fHi *= 0.5;
        }
        else
        {
            xLo = x;
            fLo = fx;
            fLo *= 0.5;
        }
        const double xNew = xLo + (xHi - xLo) * fLo / (fLo - fHi);
        if (!(xNew > xLo && xNew < xHi)) x = 0.5 * (xLo + xHi);
        else x = xNew;
        if (std::abs(xHi - xLo) < 1e-9 * std::max(xHi, 1e-30)) break;
    }

    RawFlow r;
    r.velocity = (forward ? x : -x) / pipeArea(pipe);
    r.reynolds = fluid.rho * std::abs(r.velocity) * pipe.diameter / fluid.mu;
    r.friction = churchillFriction(r.reynolds, pipe.roughness / pipe.diameter);
    Hq = pumpHeadOf(pipe, std::abs(r.velocity) * pipeArea(pipe));
    return r;
}

/// Dense symmetric-ish n x n solve: see DenseSolve.h (shared with the
/// transient kernel).

} // namespace

label OneDNetwork::addNode(const std::string& name)
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    nodeNames.push_back(name);
    nodePressureFixed.push_back(std::numeric_limits<double>::quiet_NaN());
    nodeExternalRate.push_back(0.0);
    return static_cast<label>(nodeNames.size()) - 1;
}

label OneDNetwork::nodeIndex(const std::string& name) const
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    cfFatal("network: unknown node '" + name + "'");
    return -1;
}

label OneDNetwork::addPipe(const NetworkPipe& pipe)
{
    NetworkPipe p = pipe;
    p.from = nodeIndex(p.fromName);
    p.to = nodeIndex(p.toName);
    if (p.from == p.to) cfFatal("network: pipe '" + p.name + "' connects a node to itself");
    if (!(p.diameter > 0) || !(p.length > 0))
        cfFatal("network: pipe '" + p.name + "' needs positive length and diameter");
    pipes.push_back(p);
    return static_cast<label>(pipes.size()) - 1;
}

double churchillFriction(double reynolds, double relativeRoughness)
{
    const double re = std::max(reynolds, 1e-12);
    const double a1 = std::pow(7.0 / re, 0.9) + 0.27 * relativeRoughness;
    if (a1 <= 0) cfFatal("churchillFriction: invalid Re/roughness combination");
    const double a = std::pow(2.457 * std::log(1.0 / a1), 16.0);
    const double b = std::pow(37530.0 / re, 16.0);
    return 8.0 * std::pow(std::pow(8.0 / re, 12.0) + std::pow(a + b, -1.5), 1.0 / 12.0);
}

void solvePipeFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                   double dpTotal, PipeFlowResult& out)
{
    double Hq = 0;
    const RawFlow r = equipmentFlow(pipe, fluid, dpTotal, Hq);
    out.name = pipe.name;
    out.dpTotal = dpTotal;
    out.dpFriction = dpTotal - fluid.rho * kGravity * pipe.elevation;
    out.dpValve = pipe.valveK * fluid.rho * std::abs(r.velocity) * r.velocity
                / 2.0;
    out.dpPump = -fluid.rho * kGravity * Hq;
    out.velocity = r.velocity;
    out.flowRate = r.velocity * pipeArea(pipe);
    out.reynolds = r.reynolds;
    out.frictionFactor = r.friction;

    // conductance dQ/d(dpTotal) by central difference on the monotone law;
    // floor/replace with the laminar slope if the quotient degenerates
    const double scale = std::max(std::abs(dpTotal),
                                  fluid.rho * kGravity * std::abs(pipe.elevation));
    const double h = std::max(1e-6 * std::abs(out.dpFriction), 1e-7 * scale);
    const double lam = laminarConductance(pipe, fluid);
    if (h > 0)
    {
        double Hp, Hm;
        const RawFlow rp = equipmentFlow(pipe, fluid, dpTotal + h, Hp);
        const RawFlow rm = equipmentFlow(pipe, fluid, dpTotal - h, Hm);
        const double g = (rp.velocity - rm.velocity) * pipeArea(pipe) / (2.0 * h);
        out.conductance = (std::isfinite(g) && g > 0) ? std::min(g, lam) : lam;
    }
    else
    {
        out.conductance = lam;
    }
}

SteadyNetworkResult solveSteadyNetwork(const OneDNetwork& net,
                                       const SteadyControls& ctl)
{
    SteadyNetworkResult res;
    const std::size_t n = net.nNodes();
    res.nodePressure.assign(n, 0.0);
    res.pipeFlows.resize(net.pipes.size());

    // unknown node pressures
    std::vector<label> unk(n, -1);
    std::size_t nFree = 0;
    double pSum = 0;
    std::size_t nFixed = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (std::isnan(net.nodePressureFixed[i])) unk[i] = static_cast<label>(nFree++);
        else { pSum += net.nodePressureFixed[i]; ++nFixed; }
    }
    if (nFixed == 0 && nFree > 0)
        cfFatal("steadyNetwork: no fixed-pressure node (no reference)");
    const double pInit = nFixed ? pSum / static_cast<double>(nFixed) : 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        res.nodePressure[i] = std::isnan(net.nodePressureFixed[i])
                                  ? pInit : net.nodePressureFixed[i];
    }

    // total injected rate = residual normalisation scale
    double totalIn = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (net.nodeExternalRate[i] > 0) totalIn += net.nodeExternalRate[i];
    }
    const double qScale = std::max(totalIn, 1e-300);

    // Newton on free-node pressures
    std::vector<double> p(nFree, pInit);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) p[static_cast<std::size_t>(unk[i])] = res.nodePressure[i];
    }

    auto evaluate = [&](const std::vector<double>& pv,
                        std::vector<double>& r,
                        std::vector<std::vector<double>>& j,
                        bool withJacobian) -> double
    {
        std::vector<double> pAll(n, pInit);
        for (std::size_t i = 0; i < n; ++i)
        {
            pAll[i] = std::isnan(net.nodePressureFixed[i])
                        ? pv[static_cast<std::size_t>(unk[i])]
                        : net.nodePressureFixed[i];
        }
        r.assign(nFree, 0.0);
        if (withJacobian) j.assign(nFree, std::vector<double>(nFree, 0.0));

        double maxRes = 0;
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            PipeFlowResult fr;
            if (net.hasBlackOil)
            {
                solvePipeFlowTwoPhase(pipe, net.blackOil,
                                      pAll[static_cast<std::size_t>(pipe.from)],
                                      pAll[static_cast<std::size_t>(pipe.to)], fr);
            }
            else
            {
                solvePipeFlow(pipe, net.fluid,
                              pAll[static_cast<std::size_t>(pipe.from)]
                            - pAll[static_cast<std::size_t>(pipe.to)], fr);
            }
            if (withJacobian) res.pipeFlows[k] = fr;
            const label a = unk[static_cast<std::size_t>(pipe.from)];
            const label b = unk[static_cast<std::size_t>(pipe.to)];
            // R_a -= Q ; R_b += Q
            if (a >= 0) r[static_cast<std::size_t>(a)] -= fr.flowRate;
            if (b >= 0) r[static_cast<std::size_t>(b)] += fr.flowRate;
            if (withJacobian)
            {
                const double c = fr.conductance;
                if (a >= 0)
                {
                    j[a][a] -= c;
                    if (b >= 0) { j[a][b] += c; j[b][a] += c; j[b][b] -= c; }
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (unk[i] >= 0) r[static_cast<std::size_t>(unk[i])] += net.nodeExternalRate[i];
        }
        for (double v : r) maxRes = std::max(maxRes, std::abs(v));
        return maxRes;
    };

    std::vector<double> r;
    std::vector<std::vector<double>> j;
    double maxRes = evaluate(p, r, j, true);
    res.residual = maxRes;
    res.relativeResidual = maxRes / qScale;

    int iter = 0;
    while (res.relativeResidual > ctl.tol && iter < ctl.maxIter)
    {
        if (nFree > 0)
        {
            std::vector<double> negR(r.size());
            for (std::size_t i = 0; i < r.size(); ++i) negR[i] = -r[i];
            detail::solveDenseLinear(j, negR);
            const std::vector<double> pOld = p;
            // Globalization along the Newton ray. The two-phase law has a
            // static-head knee: below it the residual is a flat plateau
            // (Q = 0), so no strict decrease exists until the knee is
            // crossed. Accept equal-residual trials on the plateau - the
            // direction there consistently pushes the source node up and
            // the iterate climbs to the knee - but require strict decrease
            // once live. Trials leaving the physical domain (non-finite or
            // non-positive pressure) are skipped.
            bool moved = false;
            for (int step = 0; step < 16 && !moved; ++step)
            {
                const double t = step < 12 ? std::pow(0.5, step)
                                           : std::pow(2.0, step - 11);
                std::vector<double> pTry = pOld;
                bool feasible = true;
                for (std::size_t i = 0; i < nFree; ++i)
                {
                    pTry[i] += ctl.relax * negR[i] * t;
                    if (!(pTry[i] > 0.0) || !std::isfinite(pTry[i]))
                        feasible = false;
                }
                if (!feasible) continue;
                maxRes = evaluate(pTry, r, j, true);
                if (maxRes <= res.residual)
                {
                    p = pTry;
                    moved = true;
                }
            }
            if (!moved)
            {
                // trials overwrote r/j with the last rejected trial's
                // state; restore them for the point we keep
                p = pOld;
                maxRes = evaluate(p, r, j, true);
            }
            res.residual = maxRes;
            res.relativeResidual = maxRes / qScale;
        }
        ++iter;
        if (nFree == 0) break;
    }
    res.iterations = iter;
    res.converged = res.relativeResidual <= ctl.tol;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) res.nodePressure[i] = p[static_cast<std::size_t>(unk[i])];
    }
    // final pipe states at the converged pressures
    for (std::size_t k = 0; k < net.pipes.size(); ++k)
    {
        const NetworkPipe& pipe = net.pipes[k];
        if (net.hasBlackOil)
        {
            solvePipeFlowTwoPhase(pipe, net.blackOil,
                                  res.nodePressure[static_cast<std::size_t>(pipe.from)],
                                  res.nodePressure[static_cast<std::size_t>(pipe.to)],
                                  res.pipeFlows[k]);
        }
        else
        {
            solvePipeFlow(pipe, net.fluid,
                          res.nodePressure[static_cast<std::size_t>(pipe.from)]
                        - res.nodePressure[static_cast<std::size_t>(pipe.to)],
                          res.pipeFlows[k]);
        }
    }
    return res;
}

namespace
{

/// Largest root of g(q) = march(pUp, q).pOutlet - pLow via Illinois. The
/// BB gradient is U-shaped in q - holdup relief lowers the required drop
/// before friction raises it - so g can vanish several times and the
/// stable operating point is the largest root. Bracketing from q=0
/// instead latches onto the near-zero root whenever the static surplus is
/// smaller than the low-rate friction hump, and the pipe law then jumps
/// between branches as the pressures move. `choked` when the march hit
/// the floor at the upper bracket.
double twoPhaseRateForDrop(const NetworkPipe& pipe,
                           const pvt::BlackOilFluid& fluid, double pUp,
                           double pLow, double sinFlow, bool& choked)
{
    choked = false;
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    auto g = [&](double q)
    {
        const BeggsBrillMarchOut m = marchBeggsBrillPipe(
            fluid, q, pUp, pipe.length, pipe.diameter, pipe.roughness,
            sinFlow, 24);
        // choked marches return the floor shortfall: continuous with the
        // unchoked branch (there pOutlet -> pFloor) and correctly scaled
        // for the secant, unlike a fixed sentinel
        return m.ok ? m.pOutlet - pLow : pFloor - pLow;
    };

    // climb the rate axis until a negative sample sits above a positive
    // one; that bracket ends at the largest root
    const double qSeed = std::max(1e-3 * area, 1e-8);
    const double g0 = g(0.0);
    double qLo = 0.0, fLo = g0;
    double qHi = qSeed, fHi = g(qHi);
    int expansions = 0;
    while (!(fHi <= 0.0 && qLo > 0.0))
    {
        if (fHi > 0.0)
        {
            qLo = qHi;
            fLo = fHi;
        }
        else if (qLo <= 0.0 && expansions >= 20)
            break; // no gas-lift dip within reach
        qHi *= 4.0;
        fHi = g(qHi);
        if (++expansions > 60)
        {
            choked = true; // pathological: no rate within range
            return qHi;
        }
    }
    if (qLo <= 0.0)
    {
        // nothing positive above zero: either the static column alone
        // reaches the low side (root near zero) or the pipe is dead
        if (g0 <= 0.0) return 0.0;
        fLo = g0;
    }

    // Illinois secant between qLo (g>0) and qHi (g<0)
    double q = qLo + (qHi - qLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fq = g(q);
        if (std::abs(fq) < 1e-7 * pUp) break;
        if (fq > 0.0)
        {
            qLo = q;
            fLo = fq;
            // Illinois: bias the retained endpoint against re-selection
            fLo *= 0.5;
        }
        else
        {
            qHi = q;
            fHi = fq;
            fHi *= 0.5;
        }
        const double qNew = qLo + (qHi - qLo) * fLo / (fLo - fHi);
        if (!(qNew > qLo && qNew < qHi)) q = 0.5 * (qLo + qHi);
        else q = qNew;
        if (std::abs(qHi - qLo) < 1e-9 * std::max(qHi, 1e-30)) break;
    }
    return q;
}

} // namespace

bool solvePipeFlowTwoPhase(const NetworkPipe& pipe,
                           const pvt::BlackOilFluid& fluid, double pFrom,
                           double pTo, PipeFlowResult& out)
{
    out = PipeFlowResult{};
    out.name = pipe.name;
    out.regime = -1;
    const double dp = pFrom - pTo;
    out.dpTotal = dp;
    if (!(pipe.length > 0) || !(pipe.diameter > 0))
        return true;

    const bool forward = dp >= 0.0; // flow from -> to when pFrom is higher
    const double pUp = std::max(pFrom, pTo);
    const double pLow = std::min(pFrom, pTo);
    const double sinFlow =
        (forward ? 1.0 : -1.0) * pipe.elevation / pipe.length;

    bool choked = false;
    const double q = twoPhaseRateForDrop(pipe, fluid, pUp, pLow, sinFlow, choked);

    const BeggsBrillMarchOut m = marchBeggsBrillPipe(
        fluid, q, pUp, pipe.length, pipe.diameter, pipe.roughness, sinFlow, 24);
    out.flowRate = forward ? q : -q;
    out.velocity = m.vmOutlet * (forward ? 1.0 : -1.0);
    out.holdup = m.holdupOutlet;
    out.regime = m.regimeMax;
    out.massRate = m.massRate;
    out.frictionFactor = m.fTpOutlet;
    out.choked = choked || !m.ok;
    out.dpFriction = dp; // gravity is inside the BB gradient

    // conductance dQ/d(dpTotal). One-sided difference probing DEEPER into
    // the live branch: a central difference straddles the static-head knee
    // whenever the low side sits within h of it, and the resulting chord
    // points Newton back into the dead region. From the dead side (q == 0)
    // probe just past the knee so the step can cross it.
    const double h = std::max(1e-4 * std::abs(dp), 1.0);
    double c = 0.0;
    if (q > 0.0)
    {
        bool ch2 = false;
        const double qDeep = twoPhaseRateForDrop(
            pipe, fluid, pUp, std::max(pLow - h, 1.0), sinFlow, ch2);
        c = (qDeep - q) / h;
    }
    if (!(c > 0.0))
    {
        // below the static head the law is flat at zero flow; probe just
        // past the knee so Newton still gets a usable slope pointing the
        // residual downhill
        const BeggsBrillMarchOut dead =
            marchBeggsBrillPipe(fluid, 0.0, pUp, pipe.length, pipe.diameter,
                                pipe.roughness, sinFlow, 24);
        const double deficit = pLow - dead.pOutlet; // missing pressure
        if (deficit > 0.0)
        {
            bool ch2 = false;
            const double probeLow = std::max(dead.pOutlet - 0.01 * pUp, 1.0);
            const double qProbe = twoPhaseRateForDrop(
                pipe, fluid, pUp, probeLow, sinFlow, ch2);
            c = qProbe / std::max(pLow - probeLow, 1.0);
        }
        if (!(c > 0.0)) c = 1e-6; // keep the network Jacobian nonsingular
    }
    out.conductance = c;
    return true;
}

OneDNetwork loadNetworkFromDict(const Dictionary& caseDict)
{
    OneDNetwork net;

    if (caseDict.has("fluid"))
    {
        const Dictionary& fd = caseDict.getDict("fluid");
        if (fd.has("blackoil"))
        {
            // two-phase mode: sources are stock-tank oil rates, per-pipe
            // hydraulics come from Beggs-Brill + black-oil PVT
            const Dictionary& bd = fd.getDict("blackoil");
            net.hasBlackOil = true;
            net.blackOil.gasGravity = bd.getScalarOr("gasGravity", net.blackOil.gasGravity);
            net.blackOil.oilDensityStd = bd.getScalarOr("oilDensity", net.blackOil.oilDensityStd);
            net.blackOil.gorStd = bd.getScalarOr("gor", net.blackOil.gorStd);
            net.blackOil.waterCut = bd.getScalarOr("waterCut", net.blackOil.waterCut);
            net.blackOil.temperature = bd.getScalarOr("temperature", net.blackOil.temperature);
            net.blackOil.surfaceTension = bd.getScalarOr("surfaceTension", net.blackOil.surfaceTension);
            if (!(net.blackOil.oilDensityStd > 0))
                cfFatal("steadyNetwork: blackoil needs positive oilDensity");
        }
        else
        {
            net.fluid.rho = fd.getScalar("rho");
            net.fluid.mu = fd.getScalar("mu");
        }
    }
    if (!net.hasBlackOil && (!(net.fluid.rho > 0) || !(net.fluid.mu > 0)))
        cfFatal("steadyNetwork: fluid needs positive rho and mu");

    if (caseDict.has("nodes"))
    {
        const Dictionary& nodes = caseDict.getDict("nodes");
        for (const std::string& name : nodes.keys())
        {
            const label idx = net.addNode(name);
            const Dictionary& nd = nodes.getDict(name);
            if (nd.has("pressure")) net.nodePressureFixed[idx] = nd.getScalar("pressure");
        }
    }

    if (!caseDict.has("pipes")) cfFatal("steadyNetwork: case has no pipes{} block");
    const Dictionary& pipes = caseDict.getDict("pipes");
    for (const std::string& name : pipes.keys())
    {
        const Dictionary& pd = pipes.getDict(name);
        NetworkPipe p;
        p.name = name;
        p.fromName = pd.getWord("from");
        p.toName = pd.getWord("to");
        p.length = pd.getScalar("length");
        p.diameter = pd.getScalar("diameter");
        p.roughness = pd.getScalarOr("roughness", 0.0);
        p.elevation = pd.getScalarOr("elevation", 0.0);
        p.valveK = pd.getScalarOr("valve", 0.0);
        p.pumpHead = pd.getScalarOr("pumpHead", 0.0);
        p.pumpRate = pd.getScalarOr("pumpRate", 0.0);
        p.pumpShutoffHead = pd.getScalarOr("pumpShutoffHead", p.pumpHead);
        if (p.valveK < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative valve"
                    " coefficient");
        if (p.pumpHead < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative pump"
                    " head");
        if (p.pumpHead > 0 && p.pumpShutoffHead < p.pumpHead)
            cfFatal("steadyNetwork: pipe '" + name + "' has a pump curve that"
                    " rises with rate (pumpShutoffHead < pumpHead)");
        net.addPipe(p);
    }

    if (caseDict.has("sources"))
    {
        const Dictionary& src = caseDict.getDict("sources");
        for (const std::string& name : src.keys())
        {
            const Dictionary& sd = src.getDict(name);
            net.setNodeRate(sd.getWord("node"), sd.getScalar("rate"));
        }
    }

    return net;
}

SteadyControls loadSteadyControls(const Dictionary& caseDict)
{
    SteadyControls ctl;
    if (caseDict.has("solve"))
    {
        const Dictionary& sd = caseDict.getDict("solve");
        ctl.maxIter = static_cast<int>(sd.getScalarOr("maxIter", ctl.maxIter));
        ctl.tol = sd.getScalarOr("tol", ctl.tol);
        ctl.relax = sd.getScalarOr("relax", ctl.relax);
    }
    return ctl;
}

} // namespace cf
