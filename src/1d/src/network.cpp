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

/// Piecewise-linear table lookup clamped at both ends (NaN/p <= 0 takes the
/// first entry, so early Newton iterates can never produce a bad property).
double pvtLookup(const std::vector<double>& grid,
                 const std::vector<double>& vals, double p)
{
    const std::size_t n = grid.size();
    if (n == 0) return 0.0;
    if (!(p > grid[0])) return vals[0];
    if (p >= grid[n - 1]) return vals[n - 1];
    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        const std::size_t mid = (lo + hi) / 2;
        if (grid[mid] <= p) lo = mid;
        else hi = mid;
    }
    const double t = (p - grid[lo]) / (grid[hi] - grid[lo]);
    return vals[lo] + (vals[hi] - vals[lo]) * t;
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

/// The Reynolds number a rheology is reported with.
///
/// Newtonian (nothing declared): rho u d / mu. Yielded with an index: the
/// apparent number 8 rho u^2 / tau_w, because no material viscosity enters the
/// law (see yieldedWallStress). Otherwise the Metzner-Reed form
///     Re_MR = rho u^(2-n) d^n / (K' 8^(n-1))
/// with the consistency corrected by the wall shear rate rather than by the
/// apparent one, K' = K ((3n+1)/(4n))^n. That exponent is the whole point:
/// Rabinowitsch's wall rate is gamma_w = (3n+1)/(4n) . 8u/d, so the stress the
/// profile actually holds is tau_w = K gamma_w^n = K' (8u/d)^n, and only with
/// that K' does Churchill's laminar line f = 64/Re read back as 8 tau_w /
/// (rho u^2), i.e. as the integrated velocity profile rather than as a
/// correlation of it. The textbook constant K (3n+1)/(4n) - the exponent-1
/// version - collapses turbulent data onto the Newtonian Moody line but sits
/// ((3n+1)/(4n))^(1-1/n) off the laminar solution, which at n = 0.5 is 25% in
/// the rate; 5.47 measured that gap from the profile integral and this line is
/// the correction. At n = 1 the factor is 1/1 = 1 exactly, so a
/// declared-but-Newtonian pipe is still bit-identical to an undeclared one.
///
/// A yielded pipe with an index (Herschel-Bulkley) has no material viscosity at
/// all - mu does not enter its law - so the only Reynolds number that can be
/// reported is the apparent one: mu_app = tau_w d / (8 u), read off the wall
/// stress the fluid actually holds at this rate. With that number and the exact
/// f = 8 tau_w / (rho u^2) the product f * Re is 64 by definition, so the
/// assertion that leans on is the closed form for tau_w, not this identity.
double yieldedWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q);

double generalizedReynolds(const NetworkPipe& pipe, const FluidProperties& fluid,
                           double speed)
{
    if (pipe.flowIndex > 0 && pipe.yieldStress > 0)
    {
        const double tw = yieldedWallStress(pipe, fluid,
                                           speed * pipeArea(pipe));
        return 8.0 * fluid.rho * speed * speed / tw;
    }
    if (!(pipe.flowIndex > 0))
        return fluid.rho * speed * pipe.diameter / fluid.mu;
    const double n = pipe.flowIndex;
    const double kp = pipe.consistency
                    * std::pow((3.0 * n + 1.0) / (4.0 * n), n);
    return fluid.rho * std::pow(speed, 2.0 - n) * std::pow(pipe.diameter, n)
         / (kp * std::pow(8.0, n - 1.0));
}

/// Buckingham (1921) exact laminar solution for a Bingham plastic, inverted
/// for the wall stress. With tw = dp d / (4 L) and xi = tau0 / tw the rate is
///     Q = (pi d^3 / (32 mu)) tw [1 - (4/3) xi + (1/3) xi^4],
/// which rearranges (multiply the bracket by tw) into one monotone equation
///     F(tw) = tw - (4/3) tau0 + tau0^4 / (3 tw^3) = B,  B = 32 mu Q / (pi d^3)
/// with F(tau0) = 0 and F' = 1 - (tau0/tw)^4 > 0 for tw > tau0, so the root is
/// unique and the bracket [tau0, B + (4/3) tau0] is analytic -- the right end
/// already gives B + tau0^4/(3 tw^3) >= B. Newton with a bisection fallback
/// then lands on the root to rounding, which is what the closed-form checks
/// need: a 1e-9 slack here would be indistinguishable from a wrong law.
/// tw = tau0 (Q = 0) and tau0 -> 0 both come out right: the first is the left
/// end of the bracket, the second makes F(tw) = tw identically.
double binghamWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q)
{
    const double tau0 = pipe.yieldStress;
    const double b = 32.0 * fluid.mu * std::abs(q)
                   / (M_PI * pipe.diameter * pipe.diameter * pipe.diameter);
    if (!(b > 0.0)) return tau0; // no rate: the wall holds exactly the threshold

    double lo = tau0, hi = b + (4.0 / 3.0) * tau0;
    double tw = 0.5 * (lo + hi);
    const double tau04 = tau0 * tau0 * tau0 * tau0;
    for (int it = 0; it < 100; ++it)
    {
        const double t2 = tw * tw;
        const double t4 = t2 * t2;
        const double f = tw - (4.0 / 3.0) * tau0
                       + tau04 / (3.0 * t2 * tw) - b;
        if (f > 0.0) hi = tw;
        else lo = tw;
        // F' = 1 - tau0^4/tw^4, which vanishes at the left end: bisect there
        const double df = 1.0 - tau04 / t4;
        const double next = df > 0.0 ? tw - f / df : 0.5 * (lo + hi);
        tw = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
        if (!(hi > lo)) break;
        if (std::abs(hi - lo) <= 2.0 * std::numeric_limits<double>::epsilon() * hi)
            break;
    }
    return tw;
}

/// Herschel-Bulkley (tau = tau0 + K gdot^n) exact laminar rate in a circular
/// pipe, as a function of the wall stress. With tw = dp d / (4 L) and
/// D = tw - tau0 the profile integral is termwise integrable,
///     Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n)
///         [ tau0^2/(n+1) + 2 tau0 D/(2n+1) + D^2/(3n+1) ],
/// which is not a correlation: at n = 1, K = mu the bracket times D^2 is
/// Buckingham's [1 - (4/3) xi + (1/3) xi^4] times tw^4/4, and at tau0 = 0 only
/// the last term survives and the expression is the power-law closed form of
/// 5.45. So one formula contains both laws already shipped, and the threshold
/// stays geometric: Q -> 0 as tw -> tau0, i.e. dp_static = 4 tau0 L / d for
/// every n, while the *shape* of the landing near it depends on n (the leading
/// power of D is 1 + 1/n, which is 2 for Bingham and 3 at n = 0.5).
double herschelBulkleyRate(double tw, const NetworkPipe& pipe)
{
    const double n = pipe.flowIndex;
    const double d = tw - pipe.yieldStress;
    if (!(d > 0.0) || !(tw > 0.0)) return 0.0;
    const double r = 0.5 * pipe.diameter;
    return M_PI * r * r * r * std::pow(tw, -3.0)
         * std::pow(pipe.consistency, -1.0 / n) * n
         * std::pow(d, 1.0 + 1.0 / n)
         * (pipe.yieldStress * pipe.yieldStress / (n + 1.0)
          + 2.0 * pipe.yieldStress * d / (2.0 * n + 1.0)
          + d * d / (3.0 * n + 1.0));
}

/// The same law inverted: the wall stress a demanded rate must be held at.
/// Monotone in tw (the yield only removes fluid from the slow end of the
/// profile, so a higher wall stress always carries more), so a bracket plus
/// bisection is enough - no derivative, and no Newton stall at the left end
/// where the slope vanishes. The left end is the threshold exactly; the right
/// end starts at the threshold plus the tau0 = 0 (power-law) stress for the
/// same rate - a lower bound on the root, because yielding only costs rate -
/// and doubles until the law clears the demand, which it always does since the
/// rate grows like tw^(1/n).
double herschelBulkleyWallStress(const NetworkPipe& pipe, double q)
{
    const double tau0 = pipe.yieldStress;
    const double demand = std::abs(q);
    if (!(demand > 0.0)) return tau0; // no rate: the wall holds the threshold

    const double n = pipe.flowIndex;
    const double r = 0.5 * pipe.diameter;
    const double twPl =
        pipe.consistency
      * std::pow(demand * (3.0 * n + 1.0) / (M_PI * r * r * r * n), n);
    double lo = tau0, hi = tau0 + twPl;
    for (int expand = 0; expand < 200 && !(herschelBulkleyRate(hi, pipe) > demand);
         ++expand) hi *= 2.0;

    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (herschelBulkleyRate(mid, pipe) > demand) hi = mid;
        else lo = mid;
        if (!(hi > lo)) break;
        if (std::abs(hi - lo) <= 2.0 * std::numeric_limits<double>::epsilon() * hi)
            break;
    }
    return 0.5 * (lo + hi);
}

/// Which yielded law a pipe's wall stress comes from: Buckingham when only a
/// threshold and the fluid's own mu are declared (the plastic number keeps the
/// Reynolds convention of 5.46), Herschel-Bulkley when an index is declared
/// with the threshold.
double yieldedWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q)
{
    if (pipe.flowIndex > 0) return herschelBulkleyWallStress(pipe, q);
    return binghamWallStress(pipe, fluid, q);
}

/// Friction + valve loss (Pa, positive) at an imposed rate magnitude q:
/// explicit Churchill evaluation, no iteration needed forward. A yielded pipe
/// instead answers with the wall-stress solution, whose q = 0 limit is the
/// static threshold -- that is the whole difference, and why the rate law needs
/// an explicit dead branch (see equipmentFlow).
double pipeLossAt(const NetworkPipe& pipe, const FluidProperties& fluid,
                  double q)
{
    // A coupled element's loss is its partner's statement, not ours: the
    // detailed model already reports the whole drop of that region, so adding
    // this pipe's own friction on top would charge the network twice for the
    // same wall. Linear in q because the coupling layer refreshes the
    // coefficient every outer iteration (exact for laminar, secant-like for
    // anything else).
    if (pipe.coupledResistance > 0) return pipe.coupledResistance * q;

    const double v = q / pipeArea(pipe);

    if (pipe.yieldStress > 0)
        return 4.0 * pipe.length
             * yieldedWallStress(pipe, fluid, q) / pipe.diameter
             + pipe.valveK * fluid.rho * v * v / 2.0;

    const double re = generalizedReynolds(pipe, fluid, std::abs(v));
    const double f = churchillFriction(re, pipe.roughness / pipe.diameter);
    return (f * pipe.length / pipe.diameter + pipe.valveK)
         * fluid.rho * v * v / 2.0;
}

/// Friction factor to report for a solved rate: the Churchill number for the
/// all-regime law, and for a yielded pipe the exact factor implied by that
/// law's drop, f = 8 tau_w / (rho v^2) -- the number that reproduces the drop
/// printed next to it, not a correlation pasted onto a yielded wall. For a
/// Bingham plastic the Reynolds number printed beside it is the plastic one
/// (rho v d / mu), so the pair reads f Re g(xi) = 64; for Herschel-Bulkley
/// there is no material viscosity and the printed number is the apparent one,
/// which makes f Re = 64 an identity between the two reports.
/// At v = 0 there is no rate to normalise against, so the report is 0.
double pipeFrictionReport(const NetworkPipe& pipe, const FluidProperties& fluid,
                          double reynolds, double speed)
{
    if (pipe.yieldStress > 0)
    {
        if (!(speed > 0)) return 0.0;
        return 8.0 * yieldedWallStress(pipe, fluid, speed * pipeArea(pipe))
             / (fluid.rho * speed * speed);
    }
    return churchillFriction(reynolds, pipe.roughness / pipe.diameter);
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

    // The dead branch of a yielded pipe. loss(0) is not 0 there, it is the
    // static threshold 4 tau0 L / d, so g(0) can already sit on the far side of
    // the crossing: the asked drop never reaches the wall's threshold and the
    // exact answer is no flow at all. Falling through would expand the bracket
    // 60 times on two same-signed samples and hand the secant a garbage pair,
    // so this has to be decided before the bracket, not damped inside it.
    // Newtonian and power-law pipes keep the plateau branch untouched -- their
    // g(0) is the sign that sent them here in the first place.
    if (pipe.yieldStress > 0 && (forward ? fLo >= 0.0 : fLo <= 0.0))
    {
        Hq = pumpHeadOf(pipe, 0.0);
        return RawFlow();
    }

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
    r.reynolds = generalizedReynolds(pipe, fluid, std::abs(r.velocity));
    r.friction = pipeFrictionReport(pipe, fluid, r.reynolds,
                                    std::abs(r.velocity));
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
    // for a coupled element the slope of its own law is the exact bound; the
    // Hagen-Poiseuille limit of the pipe geometry is not its physics
    const double lam = pipe.coupledResistance > 0
                           ? 1.0 / pipe.coupledResistance
                           : laminarConductance(pipe, fluid);
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
            else if (net.hasPvtFluid)
            {
                solvePipeFlowPvt(pipe, net.pvtFluid,
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
        else if (net.hasPvtFluid)
        {
            solvePipeFlowPvt(pipe, net.pvtFluid,
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

double PvtFluid::density(double p) const { return pvtLookup(pressure, rho, p); }

double PvtFluid::viscosity(double p) const { return pvtLookup(pressure, mu, p); }

namespace
{

/// Segment-wise march of the local mechanical-energy balance
/// dp/ds = -rho(p) g sinFlow - f(rho v, mu(p)) rho v |v| / (2D)
///         + rho(p) g H(q)/L
/// from pUp toward the far end, with rho/mu evaluated at the segment-start
/// pressure and the valve K spread uniformly over the segments. The mass
/// rate is fixed by the UPSTREAM density (w = rho(pUp) q), so the local
/// velocity grows as the fill expands along the pipe; with a constant
/// table every segment reproduces the constant-property law exactly and
/// the march coincides with solvePipeFlow's implicit pipe law.
struct PvtMarchOut
{
    bool ok = false;
    double pOutlet = 0;
    double velocity = 0; // outlet velocity, flow direction
    double reynolds = 0; // outlet conditions, as in the BB march
    double friction = 0;
};

PvtMarchOut marchPvtPipe(const NetworkPipe& pipe, const PvtFluid& fluid,
                         double qIn, double pUp, double sinFlow, int nSeg)
{
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    const double massRate = fluid.density(pUp) * qIn;
    const double ds = pipe.length / nSeg;
    const double dz = sinFlow * ds;
    const double headPerM = pumpHeadOf(pipe, qIn) / pipe.length;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    PvtMarchOut m;
    double p = pUp;
    for (int i = 0; i < nSeg; ++i)
    {
        const double rho = fluid.density(p);
        const double mu = fluid.viscosity(p);
        const double v = massRate / (rho * area);
        const double re = rho * v * pipe.diameter / mu;
        const double f = churchillFriction(re, pipe.roughness / pipe.diameter);
        m.velocity = v;
        m.reynolds = re;
        m.friction = f;
        const double dp = rho * kGravity * dz
            + (f * ds / pipe.diameter + pipe.valveK / nSeg) * rho * v * v / 2.0
            - rho * kGravity * headPerM * ds;
        p -= dp;
        if (!(p > pFloor))
        {
            m.pOutlet = pFloor;
            return m; // ok stays false: the fill cannot carry this rate
        }
    }
    m.ok = true;
    m.pOutlet = p;
    return m;
}

/// Rate landing the march on pLow via Illinois. The loss grows without
/// bound in |q| (friction + the pump curve falling) while the static head
/// stays bounded, so the law is monotone with a single root: a plain
/// expanding bracket finds it, unlike the U-shaped two-phase gradient.
double pvtRateForDrop(const NetworkPipe& pipe, const PvtFluid& fluid,
                      double pUp, double pLow, double sinFlow, bool& choked)
{
    choked = false;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    auto g = [&](double q)
    {
        const PvtMarchOut m = marchPvtPipe(pipe, fluid, q, pUp, sinFlow, 24);
        // floor-crossing marches return the floor shortfall: continuous
        // with the unchoked branch and correctly scaled for the secant
        return m.ok ? m.pOutlet - pLow : pFloor - pLow;
    };

    const double qSeed = std::max(1e-3 * 0.25 * M_PI * pipe.diameter
                                      * pipe.diameter, 1e-8);
    const double g0 = g(0.0);
    if (g0 <= 0.0) return 0.0; // the static column alone reaches the low side

    double qLo = 0.0, fLo = g0;
    double qHi = qSeed, fHi = g(qHi);
    int expansions = 0;
    while (fHi > 0.0)
    {
        qLo = qHi;
        fLo = fHi;
        qHi *= 4.0;
        fHi = g(qHi);
        if (++expansions > 60)
        {
            choked = true; // no rate within reach lands on the low side
            return qHi;
        }
    }

    double q = qLo + (qHi - qLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fq = g(q);
        if (std::abs(fq) < 1e-12 * pUp) break;
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
        if (std::abs(qHi - qLo) < 1e-12 * std::max(qHi, 1e-30)) break;
    }
    return q;
}

} // namespace

bool solvePipeFlowPvt(const NetworkPipe& pipe, const PvtFluid& fluid,
                      double pFrom, double pTo, PipeFlowResult& out)
{
    out = PipeFlowResult{};
    out.name = pipe.name;
    out.regime = -1;
    const double dp = pFrom - pTo;
    out.dpTotal = dp;
    if (!(pipe.length > 0) || !(pipe.diameter > 0)
        || fluid.pressure.size() < 2)
        return true;

    const bool forward = dp >= 0.0; // flow from -> to when pFrom is higher
    const double pUp = std::max(pFrom, pTo);
    const double pLow = std::min(pFrom, pTo);
    const double sinFlow =
        (forward ? 1.0 : -1.0) * pipe.elevation / pipe.length;
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;

    bool choked = false;
    const double q = pvtRateForDrop(pipe, fluid, pUp, pLow, sinFlow, choked);

    const PvtMarchOut m = marchPvtPipe(pipe, fluid, q, pUp, sinFlow, 24);
    out.flowRate = (forward ? 1.0 : -1.0) * q;
    out.velocity = (forward ? 1.0 : -1.0) * m.velocity;
    out.reynolds = m.reynolds;
    out.frictionFactor = m.friction;
    out.choked = choked || !m.ok;
    out.dpFriction = dp; // gravity and the property variation are in the march
    out.dpValve = pipe.valveK * fluid.density(pUp) * std::abs(out.velocity)
                * out.velocity / 2.0;
    out.dpPump = -fluid.density(pUp) * kGravity * pumpHeadOf(pipe, q);

    // conductance dQ/d(dpTotal). One-sided difference probing DEEPER into
    // the live branch: below the static head the law is flat at zero flow
    // (same knee as the two-phase law), and a central difference straddling
    // the knee points Newton back into the dead region. From the dead side
    // (q == 0) probe just past the knee so the step can cross it.
    const double h = std::max(1e-4 * std::abs(dp), 1.0);
    double c = 0.0;
    if (q > 0.0)
    {
        bool ch2 = false;
        const double qDeep = pvtRateForDrop(
            pipe, fluid, pUp, std::max(pLow - h, 1.0), sinFlow, ch2);
        c = (qDeep - q) / h;
    }
    if (!(c > 0.0))
    {
        const PvtMarchOut dead =
            marchPvtPipe(pipe, fluid, 0.0, pUp, sinFlow, 24);
        const double deficit = pLow - dead.pOutlet; // missing pressure
        if (deficit > 0.0)
        {
            bool ch2 = false;
            const double probeLow = std::max(dead.pOutlet - 0.01 * pUp, 1.0);
            const double qProbe = pvtRateForDrop(
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
        if (fd.has("blackoil") && fd.has("pvt"))
            cfFatal("steadyNetwork: fluid cannot be both blackoil and pvt");
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
        else if (fd.has("pvt"))
        {
            // tabulated-property single-phase mode: rho(p)/mu(p) come from
            // the case table (piecewise-linear, clamped at the ends)
            const Dictionary& pd = fd.getDict("pvt");
            net.hasPvtFluid = true;
            net.pvtFluid.pressure = pd.getScalarList("pressure");
            net.pvtFluid.rho = pd.getScalarList("rho");
            net.pvtFluid.mu = pd.getScalarList("mu");
            const std::size_t n = net.pvtFluid.pressure.size();
            if (n < 2 || net.pvtFluid.rho.size() != n
                || net.pvtFluid.mu.size() != n)
                cfFatal("steadyNetwork: fluid pvt needs pressure/rho/mu"
                        " lists of equal length (at least 2)");
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!(net.pvtFluid.pressure[i] > 0)
                    || !(net.pvtFluid.rho[i] > 0)
                    || !(net.pvtFluid.mu[i] > 0))
                    cfFatal("steadyNetwork: fluid pvt entries must be"
                            " positive");
                if (i > 0 && !(net.pvtFluid.pressure[i]
                               > net.pvtFluid.pressure[i - 1]))
                    cfFatal("steadyNetwork: fluid pvt pressures must"
                            " strictly ascend");
            }
        }
        else
        {
            net.fluid.rho = fd.getScalar("rho");
            net.fluid.mu = fd.getScalar("mu");
        }
    }
    if (!net.hasBlackOil && !net.hasPvtFluid
        && (!(net.fluid.rho > 0) || !(net.fluid.mu > 0)))
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
        p.coupledResistance = pd.getScalarOr("coupledResistance", 0.0);
        p.flowIndex = pd.getScalarOr("flowIndex", 0.0);
        p.consistency = pd.getScalarOr("consistency", 0.0);
        p.yieldStress = pd.getScalarOr("yieldStress", 0.0);
        if (p.yieldStress < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative yield"
                    " stress (0 means no threshold)");
        if (p.yieldStress > 0 && (net.hasBlackOil || net.hasPvtFluid))
            cfFatal("steadyNetwork: pipe '" + name + "' carries a yield stress,"
                    " which only the constant-property single-phase laminar"
                    " law computes");
        if (p.yieldStress > 0 && p.coupledResistance > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " yield stress either");
        if (p.flowIndex < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative"
                    " power-law flow index (0 means Newtonian)");
        if (p.flowIndex > 0 && !(p.consistency > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' declares a flow index"
                    " without a positive consistency: at n != 1 the fluid's mu"
                    " is not the same quantity as K (Pa.s^n)");
        if (p.flowIndex > 0 && (net.hasBlackOil || net.hasPvtFluid))
            cfFatal("steadyNetwork: pipe '" + name + "' carries a rheology,"
                    " which only the constant-property single-phase law"
                    " computes");
        if (p.flowIndex > 0 && p.coupledResistance > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " rheology either");
        if (p.valveK < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative valve"
                    " coefficient");
        if (p.pumpHead < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative pump"
                    " head");
        if (p.coupledResistance < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative coupled"
                    " resistance");
        if (p.coupledResistance > 0 && (net.hasBlackOil || net.hasPvtFluid))
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element,"
                    " which the single-phase constant-property law is the only"
                    " one that carries");
        if (p.coupledResistance > 0 && (p.valveK > 0 || p.pumpHead > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " valve or a pump");
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
