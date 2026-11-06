// Unit tests for the partitioned 1D<->3D interface loop, driven against
// analytic partner models so the loop itself -- not the 3D solve behind it --
// is what fails when it is wrong.
//
// The four things these cases pin down:
//   * the loop lands on the state where the drop the network charges equals
//     the drop the partner reports at the rate the network delivers;
//   * for a linear partner one coefficient refresh is enough, and the answer
//     does not depend on the starting coefficient or on the relaxation;
//   * the numbers a caller reads out of CoupledResult describe ONE state:
//     the rate the partner was last asked for is the rate being reported;
//   * the transient loop re-solves one time step until the interface agrees
//     and then commits exactly that step, so its ramp follows the discrete
//     closed form of the inertial network it is coupled into;
//   * the region's inertance is one number with two possible sources, and the
//     network solves with whichever one the case names.
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include "cf/1d/Network.h"
#include "cf/1d/Transient.h"
#include "cf/common/Error.h"
#include "cf/coupling/CoupledNetwork.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }\
        else         { std::cout << "ok:   " << (msg) << "\n"; }             \
    } while (0)

static bool relNear(double a, double b, double tol = 1e-8)
{
    const double s = std::max(std::abs(a), std::abs(b));
    if (s == 0) return true;
    return std::abs(a - b) / s < tol;
}

/// dp = R0 q + k q^2: k = 0 is a laminar detailed region, k > 0 adds the
/// turbulent form loss a real block would carry. Records what it was last
/// asked for so the reported state can be checked against it.
class AnalyticElement : public coupling::ICoupledElement
{
public:
    AnalyticElement(double r0, double k) : r0_(r0), k_(k) {}

    double lossAt(double qMag) override
    {
        lastQ_ = qMag;
        lastDp_ = r0_ * qMag + k_ * qMag * qMag;
        ++calls_;
        return lastDp_;
    }

    double lastQ() const { return lastQ_; }
    double lastDp() const { return lastDp_; }
    long calls() const { return calls_; }

private:
    double r0_, k_;
    double lastQ_ = 0, lastDp_ = 0;
    long calls_ = 0;
};

namespace
{

NetworkPipe linearPipe(const std::string& name, const std::string& from,
                       const std::string& to, double resistance)
{
    NetworkPipe p;
    p.name = name;
    p.fromName = from;
    p.toName = to;
    p.length = 10.0;
    p.diameter = 0.05;
    // a coupled edge carries no 1D friction of its own, which makes it an
    // exact linear resistor -- the closed form the interface tests need
    p.coupledResistance = resistance;
    return p;
}

/// One fixed-pressure source, one fixed-pressure sink.
OneDNetwork twoNodeNetwork(double pFrom, double pTo)
{
    OneDNetwork net;
    net.fluid = FluidProperties{800.0, 0.1};
    net.addNode("Src");
    net.addNode("Snk");
    net.setNodePressure("Src", pFrom);
    net.setNodePressure("Snk", pTo);
    return net;
}

/// Rate-driven chain with the coupled edge in parallel with a bypass:
/// Src -> A -> {edge || By} -> B -> Snk.
OneDNetwork branchedNetwork(double qIn)
{
    OneDNetwork net;
    net.fluid = FluidProperties{800.0, 0.1};
    net.addNode("Src");
    net.addNode("J1");
    net.addNode("J2");
    net.addNode("Snk");
    net.setNodePressure("Snk", 1.0e5);
    net.setNodeRate("Src", qIn);
    net.addPipe(linearPipe("A", "Src", "J1", 1.0e6));
    net.addPipe(linearPipe("X", "J1", "J2", 1.0e6)); // the coupled edge
    net.addPipe(linearPipe("By", "J1", "J2", 4.0e6));
    net.addPipe(linearPipe("B", "J2", "Snk", 1.0e6));
    return net;
}

SteadyControls netCtl()
{
    SteadyControls c;
    c.maxIter = 200;
    c.tol = 1e-12;
    return c;
}

/// A transient pipe handed to a partner: length and diameter set its
/// inertance rho*L/A, `resistance` its loss law.
TransientPipeConfig coupledPipe(const std::string& name, const std::string& from,
                                const std::string& to, double length,
                                double resistance)
{
    TransientPipeConfig p;
    p.name = name;
    p.fromName = from;
    p.toName = to;
    p.length = length;
    p.diameter = 0.05;
    p.coupledResistance = resistance;
    return p;
}

/// Rate-driven transient chain with the coupled edge in parallel with a
/// bypass: Src -> A -> {edge || By} -> B -> Snk. Every pipe is a linear
/// resistor, so the whole system is exactly solvable: the incompressible
/// ports force Q_A = Q_B = qIn and Q_By = qIn - Q_X from the first step on, and
/// the two parallel momentum rows add up to
///   (M_X + M_By) dQ_X/dt + (R_X + R_By) Q_X = R_By qIn,
/// which implicit Euler relaxes as Q_X(n) = qInf + (Q_X(1) - qInf)(1 + dt/tau)^(1-n).
/// The very first step is different: the network starts at rest while the
/// injection is already on, so continuity cannot hold at t = 0 and the sudden
/// rate splits by inertance rather than by resistance,
///   Q_X(1) = (M_By/dt + R_By) qIn / ((M_X + M_By)/dt + R_X + R_By).
/// That is the statement the interface must reproduce -- a model that dropped
/// the region's inertance would start at qInf instead.
/// All four pipes are coupled so no Churchill friction disturbs the linearity.
TransientCase transientBranched(double qIn, double edgeSeed)
{
    TransientCase c;
    c.fluid = FluidProperties{800.0, 0.1};
    TransientNodeConfig src, j1, j2, snk;
    src.name = "Src"; j1.name = "J1"; j2.name = "J2";
    snk.name = "Snk"; snk.pressureFixed = true; snk.pressure = 1.0e5;
    // incompressible ports (no bulk modulus) that still hold liquid, because
    // the injected rate has to go somewhere in the inventory the kernel tracks
    for (TransientNodeConfig* n : {&src, &j1, &j2})
    {
        n->volume = 1.0;
        n->liquidFraction = 0.5;
    }
    c.nodes = {src, j1, j2, snk};
    c.pipes = {coupledPipe("A", "Src", "J1", 10.0, 1.0e6),
               coupledPipe("X", "J1", "J2", 10.0, edgeSeed),
               coupledPipe("By", "J1", "J2", 1.0, 1.0e6),
               coupledPipe("B", "J2", "Snk", 10.0, 1.0e6)};
    TransientSource feed;
    feed.nodeName = "Src";
    feed.rate = qIn;
    c.sources = {feed};
    return c;
}

} // namespace

int main()
{
    const SteadyControls nc = netCtl();

    // ---- linear partner: one refresh, and the answer is the closed form ----
    {
        OneDNetwork net = twoNodeNetwork(2.0e5, 1.0e5);
        const label edge = net.addPipe(linearPipe("X", "Src", "Snk", 1.0e6));
        AnalyticElement block(4.0e6, 0.0);

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {edge}, {&block}, ctl, nc);

        CHECK(res.converged, "linear partner converges");
        CHECK(res.iterations == 2, "linear partner needs one coefficient refresh");
        CHECK(relNear(res.state[0].q, 0.025, 1e-12), "rate is dp/R of the partner");
        CHECK(relNear(res.state[0].resistance, 4.0e6, 1e-12),
              "coefficient is the partner's");
        CHECK(relNear(net.pipes[edge].coupledResistance, res.state[0].resistance, 1e-15),
              "network carries the reported coefficient");
        CHECK(relNear(block.lastQ(), res.state[0].q, 1e-15),
              "the partner was measured at the reported rate");
        CHECK(relNear(block.lastDp(), res.state[0].dpPa, 1e-15),
              "the reported drop is the partner's own answer");
        CHECK(relNear(res.state[0].dpPa, res.state[0].resistance * res.state[0].q, 1e-10),
              "reported drop equals coefficient times rate");
        CHECK(res.history.size() == 2, "history holds every interface pass");
    }

    // ---- the same answer from a relaxation, and from another seed ----
    {
        OneDNetwork damped = twoNodeNetwork(2.0e5, 1.0e5);
        const label e1 = damped.addPipe(linearPipe("X", "Src", "Snk", 1.6e6));
        AnalyticElement b1(4.0e6, 0.0);
        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        ctl.relax = 0.5;
        const coupling::CoupledResult r1 =
            coupling::solveCoupled1d3d(damped, {e1}, {&b1}, ctl, nc);

        OneDNetwork quick = twoNodeNetwork(2.0e5, 1.0e5);
        const label e2 = quick.addPipe(linearPipe("X", "Src", "Snk", 4.0e6));
        AnalyticElement b2(4.0e6, 0.0);
        coupling::CoupledControls ctl2 = ctl;
        ctl2.relax = 1.0;
        const coupling::CoupledResult r2 =
            coupling::solveCoupled1d3d(quick, {e2}, {&b2}, ctl2, nc);

        CHECK(r1.converged && r2.converged, "both interface schemes converge");
        CHECK(relNear(r1.state[0].q, r2.state[0].q, 1e-10),
              "relaxation does not move the fixed point");
        CHECK(relNear(r1.state[0].resistance, 4.0e6, 1e-10), "relaxed coefficient is the same");
        CHECK(r1.iterations > r2.iterations, "relaxation costs passes to get there");
    }

    // ---- nonlinear partner: converges to the exact operating point ----
    {
        // dp = R0 q + k q^2 = 1e5 with the network charging R q
        const double R0 = 1.0e6, k = 1.0e6, dp = 1.0e5;
        const double qExact = (-R0 + std::sqrt(R0 * R0 + 4.0 * k * dp)) / (2.0 * k);

        OneDNetwork net = twoNodeNetwork(2.0e5, 1.0e5);
        const label edge = net.addPipe(linearPipe("X", "Src", "Snk", 1.0e6));
        AnalyticElement block(R0, k);

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-12;
        ctl.relax = 0.5;
        ctl.maxOuterIterations = 400;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {edge}, {&block}, ctl, nc);

        CHECK(res.converged, "quadratic partner converges");
        CHECK(relNear(res.state[0].q, qExact, 1e-9),
              "nonlinear fixed point is the exact root");
        CHECK(res.residual < 1e-10, "interface residual reaches the tolerance");
        CHECK(relNear(block.lastQ(), res.state[0].q, 1e-15),
              "nonlinear partner measured the reported rate");
    }

    // ---- a branched network: the loop splits flow with the bypass ----
    {
        const double qIn = 0.02, rBlock = 4.0e6, rBy = 4.0e6, rLead = 1.0e6;
        OneDNetwork net = branchedNetwork(qIn);
        const label edge = 1; // "X"
        CHECK(net.pipes[edge].name == "X", "the coupled edge is the one added second");
        net.pipes[edge].coupledResistance = 1.0e6;
        AnalyticElement block(rBlock, 0.0);

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {edge}, {&block}, ctl, nc);

        // current divider between the partner and the 1D bypass
        const double qBlockExact = qIn * rBy / (rBlock + rBy);
        const double dpBranchExact = qBlockExact * rBlock;
        const double pSnk = 1.0e5;
        const double pJ2Exact = pSnk + qIn * rLead;
        const double pJ1Exact = pJ2Exact + dpBranchExact;
        const double pSrcExact = pJ1Exact + qIn * rLead;

        CHECK(res.converged, "branched network converges");
        CHECK(relNear(res.state[0].q, qBlockExact, 1e-10), "coupled edge carries its branch share");
        CHECK(relNear(res.network.nodePressure[0], pSrcExact, 1e-10),
              "source node pressure matches the closed form");
        CHECK(relNear(res.network.nodePressure[1], pJ1Exact, 1e-10),
              "upstream interface node pressure matches");
        CHECK(relNear(res.network.nodePressure[2], pJ2Exact, 1e-10),
              "downstream interface node pressure matches");

        const auto& pf = res.network.pipeFlows;
        CHECK(relNear(pf[0].flowRate, qIn, 1e-10), "the lead takes the whole rate");
        CHECK(relNear(pf[0].flowRate, res.state[0].q + pf[2].flowRate, 1e-10),
              "node split: in = coupled + bypass");
        CHECK(relNear(pf[1].dpTotal, pf[2].dpTotal, 1e-12),
              "parallel branches share their drop");
        CHECK(relNear(pf[0].dpTotal + pf[1].dpTotal + pf[3].dpTotal,
                      pSrcExact - pSnk, 1e-10), "pressure chain closes");
    }

    // ---- an exhausted loop still reports a consistent state ----
    {
        OneDNetwork net = twoNodeNetwork(2.0e5, 1.0e5);
        const label edge = net.addPipe(linearPipe("X", "Src", "Snk", 1.0e6));
        AnalyticElement block(4.0e6, 0.0);

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-14;      // below what three relaxed passes can reach
        ctl.relax = 0.5;            // ... on a rate of one half per pass
        ctl.maxOuterIterations = 3;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {edge}, {&block}, ctl, nc);

        CHECK(!res.converged, "a stopped loop says it stopped");
        CHECK(res.iterations == 3, "every allowed pass ran");
        CHECK(res.history.size() == 3, "each pass is in the history");
        CHECK(relNear(block.lastQ(), res.state[0].q, 1e-15),
              "the stopped state is still one state");
        CHECK(relNear(block.lastDp(), res.state[0].dpPa, 1e-15),
              "the reported drop is the partner's own last answer");
        CHECK(relNear(res.residual,
                      std::abs(res.state[0].resistance * res.state[0].q - res.state[0].dpPa)
                          / res.state[0].dpPa, 1e-10),
              "the residual is the gap between the two solvers");
    }

    // ---- transient loop against a linear partner: the discrete closed form ----
    {
        const double qIn = 0.02, rBlock = 4.0e6, rBy = 1.0e6;
        const double dt = 0.15, endTime = 1.8;
        const int nSteps = 12;
        const double section = 0.25 * M_PI * 0.05 * 0.05;
        const double mX = 800.0 * 10.0 / section, mBy = 800.0 * 1.0 / section;
        const double tau = (mX + mBy) / (rBlock + rBy);
        const double qInf = qIn * rBy / (rBlock + rBy);
        const double decay = 1.0 + dt / tau;
        const double qFirst = (mBy / dt + rBy) * qIn / ((mX + mBy) / dt + rBlock + rBy);
        CHECK(qFirst < 0.75 * qInf,
              "the case really ramps: its inertial and resistive dividers differ");

        TransientNetwork net(transientBranched(qIn, 1.0e6));
        const label edge = 1; // "X"
        AnalyticElement block(rBlock, 0.0);

        coupling::TransientCoupledControls ctl;
        ctl.tolerance = 1e-12;
        const coupling::TransientCoupledResult res =
            coupling::solveTransientCoupled1d3d(net, {edge}, {&block}, ctl, dt, endTime);

        CHECK(res.steps.size() == static_cast<std::size_t>(nSteps),
              "transient loop runs every step of the interval");
        CHECK(relNear(block.lastQ(), res.state[0].q, 1e-15),
              "the reported state is the one the partner was measured at");
        CHECK(relNear(res.steps[0].edges[0].q, qFirst, 1e-9),
              "the first step splits the sudden rate by inertance");
        CHECK(relNear(res.state[0].resistance, rBlock, 1e-12),
              "the interface lands on the linear partner's coefficient");
        CHECK(relNear(net.coupledResistanceOf(edge), res.state[0].resistance, 1e-15),
              "the network carries the coefficient the last step was solved with");
        CHECK(relNear(net.time(), endTime, 1e-15), "the clock reaches endTime");

        // a linear partner's fixed point does not depend on the rate, so the
        // first step needs one refresh and every later step is already agreed
        CHECK(res.steps[0].passes == 2, "first step needs one coefficient refresh");
        long passes = 0;
        for (const coupling::TransientCoupledStep& s : res.steps) passes += s.passes;
        CHECK(passes == nSteps + 1, "later steps agree without another partner run");
        CHECK(block.calls() == nSteps + 1, "every partner run is one interface pass");

        bool everyStepAgreed = true, rampExact = true, identityExact = true;
        bool chainExact = true;
        double qPrev = 0.0;
        for (int n = 1; n <= nSteps; ++n)
        {
            const coupling::TransientCoupledStep& s = res.steps[static_cast<std::size_t>(n - 1)];
            const double qExact = n == 1
                ? qFirst
                : qInf + (qFirst - qInf) * std::pow(decay, 1 - n);
            everyStepAgreed = everyStepAgreed && s.edges[0].residual <= ctl.tolerance;
            rampExact = rampExact && relNear(s.edges[0].q, qExact, 1e-9);
            identityExact = identityExact
                && relNear(s.flowRate[0], qIn, 1e-10)
                && relNear(s.flowRate[3], qIn, 1e-10)
                && relNear(s.flowRate[1] + s.flowRate[2], qIn, 1e-10);
            // what the region's two statements say about the same interface:
            // the node pressure difference drives both the partner's resistance
            // and the region's own acceleration
            const double dpNodes = s.pressure[1] - s.pressure[2];
            chainExact = chainExact
                && relNear(dpNodes, s.edges[0].dpPa + mX * (s.edges[0].q - qPrev) / dt, 1e-8);
            qPrev = s.edges[0].q;
        }
        CHECK(everyStepAgreed, "every time step closes the interface gap");
        CHECK(rampExact, "the inertial ramp follows the discrete closed form");
        CHECK(identityExact, "incompressible ports hold the series rate at every step");
        CHECK(chainExact, "region drop equals partner drop plus its own inertance");

        const coupling::TransientCoupledStep& first = res.steps[0];
        CHECK(relNear(first.pressure[3], 1.0e5, 1e-15),
              "the fixed node keeps its pressure at every time level");
        CHECK(res.steps[0].history.size() == 2 && res.steps[1].history.size() == 1,
              "the history holds every pass of every step");
    }

    // ---- relaxation and seeds change the passes, not the trajectory ----
    {
        const double qIn = 0.02, dt = 0.15, endTime = 1.8;

        TransientNetwork quick(transientBranched(qIn, 4.0e6));
        AnalyticElement b1(4.0e6, 0.0);
        coupling::TransientCoupledControls exact;
        exact.tolerance = 1e-10;
        const coupling::TransientCoupledResult r1 =
            coupling::solveTransientCoupled1d3d(quick, {1}, {&b1}, exact, dt, endTime);

        TransientNetwork damped(transientBranched(qIn, 1.0e6));
        AnalyticElement b2(4.0e6, 0.0);
        coupling::TransientCoupledControls relaxed = exact;
        relaxed.relax = 0.5;
        relaxed.maxSubIterations = 60;
        const coupling::TransientCoupledResult r2 =
            coupling::solveTransientCoupled1d3d(damped, {1}, {&b2}, relaxed, dt, endTime);

        long passes1 = 0, passes2 = 0;
        bool sameRamp = true;
        for (std::size_t n = 0; n < r1.steps.size(); ++n)
        {
            passes1 += r1.steps[n].passes;
            passes2 += r2.steps[n].passes;
            sameRamp = sameRamp
                && relNear(r1.steps[n].edges[0].q, r2.steps[n].edges[0].q, 1e-10);
        }
        CHECK(sameRamp, "the relaxed interface follows the same ramp");
        CHECK(passes2 > passes1, "relaxation costs partner runs to get there");
        CHECK(relNear(r2.state[0].resistance, 4.0e6, 1e-10), "relaxed coefficient is the same");
    }

    // ---- where the detailed region's inertance comes from ----
    {
        // Same network, same linear partner, and a partner that hands back the
        // region's inertance: with the number the pipe would have used anyway
        // the two conventions must produce one trajectory, because the
        // partitioned system is one system. With a different number the ramp is
        // a different ramp - and only the inertance changed, since a linear
        // partner's coefficient is its resistance either way.
        const double section = 0.25 * M_PI * 0.05 * 0.05;
        const double mOwn = 800.0 * 10.0 / section;      // pipe X: length 10
        const double mRegion = mOwn * (M_PI / 4.0);      // a 78.5% equal-area bore
        const double mBy = 800.0 * 1.0 / section;
        const double rBlock = 4.0e6, rBy = 1.0e6, qIn = 0.02;
        const double dt = 0.15, endTime = 1.8;

        class InertialElement : public AnalyticElement
        {
        public:
            InertialElement(double r0, double m) : AnalyticElement(r0, 0.0), m_(m) {}
            double regionInertance() const override { return m_; }

        private:
            double m_;
        };

        TransientNetwork kept(transientBranched(qIn, 1.0e6));
        TransientNetwork given(transientBranched(qIn, 1.0e6));
        TransientNetwork other(transientBranched(qIn, 1.0e6));

        coupling::TransientCoupledControls base;
        base.tolerance = 1e-12;
        coupling::TransientCoupledControls fromPartner = base;
        fromPartner.inertiaFromPartner = true;

        AnalyticElement plain(rBlock, 0.0);
        InertialElement sameInertia(rBlock, mOwn), shaped(rBlock, mRegion);
        const coupling::TransientCoupledResult rKeep =
            coupling::solveTransientCoupled1d3d(kept, {1}, {&plain}, base, dt, endTime);
        const coupling::TransientCoupledResult rGive =
            coupling::solveTransientCoupled1d3d(given, {1}, {&sameInertia}, fromPartner, dt, endTime);
        const coupling::TransientCoupledResult rShape =
            coupling::solveTransientCoupled1d3d(other, {1}, {&shaped}, fromPartner, dt, endTime);

        CHECK(rGive.steps.size() == rKeep.steps.size()
                  && rShape.steps.size() == rKeep.steps.size(),
              "all three conventions run the whole series");
        bool identical = true, everyStep = true, shapeDiffers = false;
        for (std::size_t n = 0; n < rKeep.steps.size(); ++n)
        {
            identical = identical
                && relNear(rKeep.steps[n].edges[0].q, rGive.steps[n].edges[0].q, 1e-15)
                && relNear(rKeep.steps[n].pressure[1], rGive.steps[n].pressure[1], 1e-15);
            everyStep = everyStep
                && relNear(rGive.steps[n].edges[0].inertance, mOwn, 1e-15)
                && relNear(rKeep.steps[n].edges[0].inertance, mOwn, 1e-15)
                && relNear(rShape.steps[n].edges[0].inertance, mRegion, 1e-15);
            shapeDiffers = shapeDiffers
                || !relNear(rShape.steps[n].edges[0].q, rKeep.steps[n].edges[0].q, 1e-6);
        }
        CHECK(identical,
              "a partner that measures what the bore already says changes nothing:"
              " one rate history, one pressure history");
        CHECK(everyStep, "each step reports the inertance its own row was solved with");
        CHECK(given.pipeInertia(1) == mOwn && other.pipeInertia(1) == mRegion,
              "the network keeps the installed number after the run");
        CHECK(shapeDiffers,
              "and a different measurement is a different trajectory: the shape"
              " factor moves the ramp, not just the bookkeeping");

        // the reshaped arm, against the same closed forms the ramp test uses
        const double decayShape = 1.0 + dt / ((mRegion + mBy) / (rBlock + rBy));
        const double qInfShape = qIn * rBy / (rBlock + rBy);
        const double q1Shape = (mBy / dt + rBy) * qIn
                             / ((mRegion + mBy) / dt + rBlock + rBy);
        CHECK(relNear(rShape.steps[0].edges[0].q, q1Shape, 1e-9),
              "the first step divides the sudden rate by the installed inertance");
        CHECK(relNear(rShape.steps[3].edges[0].q - qInfShape,
                      (rShape.steps[2].edges[0].q - qInfShape) / decayShape, 1e-9),
              "and relaxes with the time constant that inertance implies");

        // a partner that states no region cannot leave the network half-declared
        {
            TransientNetwork net(transientBranched(qIn, 1.0e6));
            AnalyticElement silent(rBlock, 0.0);   // the interface's own default
            coupling::TransientCoupledControls asks = base;
            asks.inertiaFromPartner = true;
            bool refused = false;
            try
            {
                coupling::solveTransientCoupled1d3d(net, {1}, {&silent}, asks, dt, endTime);
            }
            catch (const Error&) { refused = true; }
            CHECK(refused,
                  "a convention that asks a partner for a region it does not state"
                  " is fatal, not silently the bore's number");
            CHECK(relNear(net.pipeInertia(1), mOwn, 1e-15) && net.time() == 0.0,
                  "and the refusal happens before anything is installed: the network"
                  " is still the case that was typed");
        }

        // the 1D side's own guards, in the order they are checked
        {
            TransientCase c = transientBranched(qIn, 1.0e6);
            c.pipes[3].coupledResistance = 0.0;      // pipe B: no partner at all
            TransientNetwork net(c);
            bool noPartner = false, zeroed = false, naned = false;
            try { net.setCoupledInertia(3, 1.0e6); }
            catch (const Error&) { noPartner = true; }
            try { net.setCoupledInertia(1, 0.0); }
            catch (const Error&) { zeroed = true; }
            try { net.setCoupledInertia(1, std::numeric_limits<double>::quiet_NaN()); }
            catch (const Error&) { naned = true; }
            CHECK(noPartner, "only a pipe with a partner can be given a measured inertance");
            CHECK(zeroed && naned,
                  "a non-positive or non-finite inertance is refused before it reaches"
                  " any row, so a degenerate region cannot quiet the ramp");
            CHECK(relNear(net.pipeInertia(1), mOwn, 1e-15),
                  "and a refused write leaves the row's number where it was");
            net.setCoupledInertia(1, mRegion);
            CHECK(relNear(net.pipeInertia(1), mRegion, 1e-15),
                  "pipeInertia reads back what the momentum row will use");
        }
    }

    // ---- a step the interface cannot agree on stops the run ----
    {
        TransientNetwork net(transientBranched(0.02, 1.0e6));
        AnalyticElement block(4.0e6, 0.0);
        coupling::TransientCoupledControls stuck;
        stuck.maxSubIterations = 1;   // one pass cannot reach a far seed
        stuck.tolerance = 1e-12;
        bool stopped = false;
        try
        {
            coupling::solveTransientCoupled1d3d(net, {1}, {&block}, stuck, 0.05, 0.5);
        }
        catch (const Error&) { stopped = true; }
        CHECK(stopped, "an unagreed step is fatal instead of committed");
        CHECK(net.time() == 0.0 && net.pipeFlows()[1] == 0.0,
              "a fatal interface leaves the network at its start");
    }

    // ---- two edges at once: each partner gets its own rate from one solve ----
    {
        // Src -> A -> { X1 || By } -> X2 -> Snk : X1 sits in a branch, so the
        // node split decides its rate, while X2 is in series, so continuity pins
        // it. A loop that solved the network per edge, or handed both partners
        // the same rate, fails here rather than in the 3D case.
        OneDNetwork net;
        net.fluid = FluidProperties{800.0, 0.1};
        net.addNode("Src"); net.addNode("J1"); net.addNode("J2"); net.addNode("Snk");
        net.setNodePressure("Snk", 1.0e5);
        net.setNodeRate("Src", 0.02);
        net.addPipe(linearPipe("A", "Src", "J1", 1.0e6));
        const label e1 = net.addPipe(linearPipe("X1", "J1", "J2", 1.0e6));
        net.addPipe(linearPipe("By", "J1", "J2", 4.0e6));
        const label e2 = net.addPipe(linearPipe("X2", "J2", "Snk", 1.0e6));

        AnalyticElement b1(2.0e6, 0.0), b2(8.0e6, 0.0);
        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-12;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {e1, e2}, {&b1, &b2}, ctl, nc);

        const double qX1 = 0.02 * 4.0e6 / (4.0e6 + 2.0e6);
        CHECK(res.converged, "two linear partners converge together");
        CHECK(res.edges.size() == 2 && res.state.size() == 2,
              "the result carries one state per edge");
        CHECK(res.edges[0] == e1 && res.edges[1] == e2, "the states are in edge order");
        CHECK(relNear(res.state[0].q, qX1, 1e-11),
              "the branch edge takes the closed-form split");
        CHECK(relNear(res.state[1].q, 0.02, 1e-11),
              "the series edge carries the whole rate");
        CHECK(relNear(res.state[0].resistance, 2.0e6, 1e-11)
              && relNear(res.state[1].resistance, 8.0e6, 1e-11),
              "each edge lands on its own partner's coefficient");
        CHECK(relNear(res.state[0].dpPa, 2.0e6 * qX1, 1e-11)
              && relNear(res.state[1].dpPa, 8.0e6 * 0.02, 1e-11),
              "each partner reports its own drop");
        CHECK(res.iterations == 2, "one shared refresh settles both edges");
        CHECK(b1.calls() == 2 && b2.calls() == 2,
              "one pass runs each partner once, never twice");
        CHECK(res.history.size() == 2 && res.history[1].edges.size() == 2,
              "each iteration records both edges");
        CHECK(relNear(res.residual, std::max(res.history[1].edges[0].residual,
                                             res.history[1].edges[1].residual), 1e-12),
              "the reported residual is the worst edge of the accepted pass");
        CHECK(relNear(net.pipes[static_cast<std::size_t>(e1)].coupledResistance,
                      res.state[0].resistance, 1e-15)
              && relNear(net.pipes[static_cast<std::size_t>(e2)].coupledResistance,
                         res.state[1].resistance, 1e-15),
              "the network carries both accepted coefficients");

        // the same two edges starting anywhere else reach the same state, and a
        // starting set that is far off on one edge does not stall the other
        AnalyticElement c1(2.0e6, 0.0), c2(8.0e6, 0.0);
        net.pipes[static_cast<std::size_t>(e1)].coupledResistance = 8.0e7;
        net.pipes[static_cast<std::size_t>(e2)].coupledResistance = 1.0e6;
        const coupling::CoupledResult far =
            coupling::solveCoupled1d3d(net, {e1, e2}, {&c1, &c2}, ctl, nc);
        CHECK(relNear(far.state[0].q, qX1, 1e-11) && relNear(far.state[1].q, 0.02, 1e-11),
              "the fixed point does not belong to the starting set");
        CHECK(relNear(far.state[0].resistance, res.state[0].resistance, 1e-9)
              && relNear(far.state[1].resistance, res.state[1].resistance, 1e-9),
              "both edges agree across starting sets to the interface tolerance");
    }

    // ---- two edges, one of them nonlinear: the split is the partners' ----
    {
        // X1 || X2 with nothing else between the nodes, so the drop equality and
        // the rate sum together decide the split: R1 (Q - q2) = (R0 + k q2) q2.
        // The quadratic partner is kept mild on purpose: the interface map
        // R2 -> R0 + k Q R1/(R1 + R2) has derivative -k Q R1/(R1+R2)^2 at the
        // fixed point. That derivative is bounded by -1, so a stiff k never
        // diverges the un-accelerated block Gauss-Seidel -- it only crawls (a
        // stiff k costs hundreds of passes). This block pins the split, so it
        // stays mild; the crawl is measured by the accelerated blocks below.
        const double Q = 0.02, r1 = 2.0e6, r0 = 1.0e6, k = 1.0e7;
        const double q2 = (-(r0 + r1) + std::sqrt((r0 + r1) * (r0 + r1)
                                                  + 4.0 * k * r1 * Q)) / (2.0 * k);
        CHECK(q2 > 0 && q2 < Q, "the quadratic partner has a real split");
        CHECK(std::abs(-k * Q * r1 / std::pow(r1 + r0 + k * q2, 2)) < 1.0,
              "this interface map is contractive");

        OneDNetwork net;
        net.fluid = FluidProperties{800.0, 0.1};
        net.addNode("Src"); net.addNode("J1"); net.addNode("J2"); net.addNode("Snk");
        net.setNodePressure("Snk", 1.0e5);
        net.setNodeRate("Src", Q);
        net.addPipe(linearPipe("A", "Src", "J1", 1.0e6));
        const label e1 = net.addPipe(linearPipe("X1", "J1", "J2", 1.0e6));
        const label e2 = net.addPipe(linearPipe("X2", "J1", "J2", 4.0e6));
        net.addPipe(linearPipe("B", "J2", "Snk", 1.0e6));

        AnalyticElement b1(r1, 0.0), b2(r0, k);
        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-8;
        ctl.maxOuterIterations = 30;
        const coupling::CoupledResult res =
            coupling::solveCoupled1d3d(net, {e1, e2}, {&b1, &b2}, ctl, nc);

        CHECK(res.converged, "a linear and a quadratic partner converge together");
        CHECK(relNear(res.state[1].q, q2, 1e-6), "the quadratic edge takes its root");
        CHECK(relNear(res.state[0].q, Q - q2, 1e-6), "the linear edge takes the rest");
        CHECK(relNear(res.state[0].dpPa, res.state[1].dpPa, 1e-7),
              "the two partners agree on one drop between the same nodes");
        CHECK(res.iterations > 2, "the nonlinear edge needs more than one refresh");
    }

    // ---- the loop refuses a malformed edge list ----
    {
        OneDNetwork net = branchedNetwork(0.02);
        const label edge = 1;
        AnalyticElement block(4.0e6, 0.0);
        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        auto fatal = [&ctl, &net, &nc](const std::vector<label>& edges,
                                       std::vector<coupling::ICoupledElement*> blocks)
        {
            try { coupling::solveCoupled1d3d(net, edges, blocks, ctl, nc); }
            catch (const Error&) { return true; }
            return false;
        };
        CHECK(fatal({edge, edge}, {&block, &block}),
              "the same edge twice is rejected, not solved twice");
        CHECK(fatal({edge}, {nullptr}), "an edge with no partner is rejected");
        CHECK(fatal({edge}, {&block, &block}),
              "a partner list that does not match the edges is rejected");
    }

    // ---- transient: two edges, one branch and one series, per step ----
    {
        const double qIn = 0.02, dt = 0.15, endTime = 0.6;
        TransientCase c;
        c.fluid = FluidProperties{800.0, 0.1};
        TransientNodeConfig src, j1, j2, j3, snk;
        src.name = "Src"; j1.name = "J1"; j2.name = "J2"; j3.name = "J3";
        snk.name = "Snk"; snk.pressureFixed = true; snk.pressure = 1.0e5;
        for (TransientNodeConfig* n : {&src, &j1, &j2, &j3})
        {
            n->volume = 1.0;
            n->liquidFraction = 0.5;
        }
        c.nodes = {src, j1, j2, j3, snk};
        c.pipes = {coupledPipe("A", "Src", "J1", 10.0, 1.0e6),
                   coupledPipe("X1", "J1", "J2", 10.0, 1.0e6),
                   coupledPipe("By", "J1", "J2", 1.0, 1.0e6),
                   coupledPipe("X2", "J2", "J3", 10.0, 1.0e6),
                   coupledPipe("B", "J3", "Snk", 10.0, 1.0e6)};
        TransientSource feed;
        feed.nodeName = "Src";
        feed.rate = qIn;
        c.sources = {feed};

        TransientNetwork net(c);
        const label e1 = 1, e2 = 3;
        AnalyticElement b1(4.0e6, 0.0), b2(2.0e6, 0.0);
        coupling::TransientCoupledControls ctl;
        ctl.tolerance = 1e-12;
        const coupling::TransientCoupledResult res =
            coupling::solveTransientCoupled1d3d(net, {e1, e2}, {&b1, &b2}, ctl, dt, endTime);

        const double section = 0.25 * M_PI * 0.05 * 0.05;
        const double mX1 = 800.0 * 10.0 / section, mBy = 800.0 * 1.0 / section;
        const double tau = (mX1 + mBy) / (4.0e6 + 1.0e6);
        const double qInf = qIn * 1.0e6 / (4.0e6 + 1.0e6);
        const double qFirst = (mBy / dt + 1.0e6) * qIn / ((mX1 + mBy) / dt + 4.0e6 + 1.0e6);
        const double mX2 = 800.0 * 10.0 / section;

        CHECK(res.edges.size() == 2 && res.steps[0].edges.size() == 2,
              "every time step reports both edges");
        CHECK(res.steps.back().history.back().edges.size() == 2,
              "the history records both edges of the accepted pass");

        bool agreed = true, seriesPinned = true, splitsExact = true, rowsExact = true;
        long totalPasses = 0;
        double qPrev1 = 0.0, qPrev2 = 0.0;
        for (std::size_t n = 0; n < res.steps.size(); ++n)
        {
            const coupling::TransientCoupledStep& s = res.steps[n];
            const double qExact = n == 0
                ? qFirst
                : qInf + (qFirst - qInf) * std::pow(1.0 + dt / tau, 1 - static_cast<int>(n) - 1);
            agreed = agreed && s.edges[0].residual <= ctl.tolerance
                && s.edges[1].residual <= ctl.tolerance;
            seriesPinned = seriesPinned && relNear(s.edges[1].q, qIn, 1e-9);
            splitsExact = splitsExact && relNear(s.edges[0].q, qExact, 1e-8);
            // each edge's own momentum row, closed by its own partner's drop
            const double dp1 = s.pressure[1] - s.pressure[2];
            const double dp2 = s.pressure[2] - s.pressure[3];
            rowsExact = rowsExact
                && relNear(dp1, s.edges[0].dpPa + mX1 * (s.edges[0].q - qPrev1) / dt, 1e-7)
                && relNear(dp2, s.edges[1].dpPa + mX2 * (s.edges[1].q - qPrev2) / dt, 1e-7);
            qPrev1 = s.edges[0].q;
            qPrev2 = s.edges[1].q;
            totalPasses += s.passes;
        }
        CHECK(agreed, "both interfaces close at every time level");
        CHECK(seriesPinned, "the series edge carries the injected rate at every step");
        CHECK(splitsExact, "the branch edge still follows its own inertial ramp");
        CHECK(rowsExact, "each edge's row closes against its own partner");
        CHECK(b1.calls() == totalPasses && b2.calls() == totalPasses,
              "one pass runs each partner exactly once");
        CHECK(relNear(net.coupledResistanceOf(e1), res.state[0].resistance, 1e-15)
              && relNear(net.coupledResistanceOf(e2), res.state[1].resistance, 1e-15),
              "the network carries both accepted coefficients");

        // an edge the interface cannot close stops the run and says which one
        TransientNetwork stuckNet(c);
        AnalyticElement s1(4.0e6, 0.0), s2(2.0e6, 0.0);
        coupling::TransientCoupledControls tight;
        tight.maxSubIterations = 1;
        tight.tolerance = 1e-12;
        bool stopped = false;
        std::string what;
        try
        {
            coupling::solveTransientCoupled1d3d(stuckNet, {e1, e2}, {&s1, &s2},
                                                tight, dt, 0.3);
        }
        catch (const Error& e) { stopped = true; what = e.what(); }
        CHECK(stopped, "an unagreed multi-edge step is fatal instead of committed");
        CHECK(what.find("X1") != std::string::npos, "the message names the worst edge");
        CHECK(stuckNet.time() == 0.0, "a fatal interface leaves the clock at its start");
    }

    // ---- Aitken on one edge: what a near-1 interface map costs the plain loop ----
    {
        // X || By with a quadratic X. The map the plain loop iterates is
        //   R -> r0 + k * Q * rBy / (rBy + R),
        // and at its fixed point its derivative is -k*Q*rBy/(rBy+R)^2, which for
        // a strongly form-lossy partner approaches -1: the iteration then walks as
        // a sawtooth and needs hundreds of passes, although it never diverges.
        // One secant pair per pass (Aitken delta-squared for a scalar interface)
        // is Newton along the measured direction, so it lands in a handful.
        const double Q = 0.02, rBy = 4.0e6, r0 = 1.0e6, k = 1.75e10;
        // the split the two laws agree on: rBy (Q - q) = (r0 + k q) q
        const double qX = (-(r0 + rBy) + std::sqrt((r0 + rBy) * (r0 + rBy)
                                                   + 4.0 * k * rBy * Q)) / (2.0 * k);
        const double rX = r0 + k * qX;
        const double slope = -k * Q * rBy / std::pow(rBy + rX, 2);
        CHECK(relNear(qX, 2.0e-3, 1e-12), "the chosen k puts the split where it was set");
        CHECK(slope < -0.8, "the plain interface map is near-1, i.e. nearly sawtooth");

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        ctl.relax = 1.0;

        ctl.maxOuterIterations = 40;
        OneDNetwork slow(branchedNetwork(Q));
        AnalyticElement slowBlock(r0, k);
        const coupling::CoupledResult plain =
            coupling::solveCoupled1d3d(slow, {1}, {&slowBlock}, ctl, netCtl());
        CHECK(!plain.converged, "40 passes are not enough for a near-1 map");
        CHECK(std::abs(plain.state[0].q - qX) > 1e-4 * qX,
              "and its 40th pass is still far from the split");

        ctl.maxOuterIterations = 400;
        OneDNetwork farNet(branchedNetwork(Q));
        AnalyticElement farBlock(r0, k);
        const coupling::CoupledResult slowButSure =
            coupling::solveCoupled1d3d(farNet, {1}, {&farBlock}, ctl, netCtl());
        CHECK(slowButSure.converged && slowButSure.iterations > 60,
              "the plain loop gets there, at a price in partner solves");

        OneDNetwork fastNet(branchedNetwork(Q));
        AnalyticElement fastBlock(r0, k);
        ctl.accelerate = true;
        const coupling::CoupledResult fast =
            coupling::solveCoupled1d3d(fastNet, {1}, {&fastBlock}, ctl, netCtl());
        CHECK(fast.converged, "the accelerated loop reaches the same fixed point");
        CHECK(fast.iterations <= 12, "in a dozen passes rather than hundreds");
        CHECK(slowButSure.iterations > 10 * fast.iterations,
              "so each partner solve the plain loop spends is ten it does not need");
        CHECK(relNear(fast.state[0].q, qX, 1e-6), "accelerated rate is the closed-form split");
        CHECK(relNear(fast.state[0].resistance, rX, 1e-6), "and its coefficient is");
        CHECK(relNear(fast.state[0].q, slowButSure.state[0].q, 1e-6),
              "both loops answer the same interface statement");
        CHECK(fastBlock.calls() == fast.iterations,
              "acceleration adds no partner solve per pass");
    }

    // ---- a linear interface is untouched: the first step is the plain step ----
    {
        OneDNetwork plainNet(branchedNetwork(0.02));
        AnalyticElement plainBlock(4.0e6, 0.0);
        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-12;
        ctl.maxOuterIterations = 8;
        const coupling::CoupledResult plain =
            coupling::solveCoupled1d3d(plainNet, {1}, {&plainBlock}, ctl, netCtl());

        ctl.accelerate = true;
        OneDNetwork qnNet(branchedNetwork(0.02));
        AnalyticElement qnBlock(4.0e6, 0.0);
        const coupling::CoupledResult qn =
            coupling::solveCoupled1d3d(qnNet, {1}, {&qnBlock}, ctl, netCtl());

        CHECK(plain.iterations == 2 && qn.iterations == 2,
              "one refresh is the fixed point with or without the quasi-Newton");
        CHECK(qn.history[0].edges[0].resistance == plain.history[0].edges[0].resistance,
              "the first pass has no secant pair yet, so it steps exactly as plain does");
        CHECK(relNear(qn.state[0].q, plain.state[0].q, 1e-12)
              && relNear(qn.state[0].dpPa, plain.state[0].dpPa, 1e-12),
              "and the linear answer is the same number");
    }

    // ---- two edges, one of them nearly non-iterable ----
    {
        // The stiff partner sits on X2, in parallel with the linear X1, so the
        // split is decided by the two partners and the network has to be solved
        // again for whichever coefficient the interface operator proposes.
        const double Q = 0.02, r1 = 2.0e6, r0 = 1.0e6, k = 1.75e10;
        const double q2 = (-(r0 + r1) + std::sqrt((r0 + r1) * (r0 + r1)
                                                  + 4.0 * k * r1 * Q)) / (2.0 * k);
        const double r2 = r0 + k * q2;

        auto makeNet = [Q]
        {
            OneDNetwork net;
            net.fluid = FluidProperties{800.0, 0.1};
            net.addNode("Src"); net.addNode("J1"); net.addNode("J2"); net.addNode("Snk");
            net.setNodePressure("Snk", 1.0e5);
            net.setNodeRate("Src", Q);
            net.addPipe(linearPipe("A", "Src", "J1", 1.0e6));
            net.addPipe(linearPipe("X1", "J1", "J2", 1.0e6));
            net.addPipe(linearPipe("X2", "J1", "J2", 4.0e6));
            net.addPipe(linearPipe("B", "J2", "Snk", 1.0e6));
            return net;
        };

        coupling::CoupledControls ctl;
        ctl.tolerance = 1e-10;
        ctl.maxOuterIterations = 400;

        OneDNetwork slowNet(makeNet());
        AnalyticElement b1(2.0e6, 0.0), b2(r0, k);
        const coupling::CoupledResult plain =
            coupling::solveCoupled1d3d(slowNet, {1, 2}, {&b1, &b2}, ctl, netCtl());
        CHECK(plain.converged, "the two-edge plain loop does reach the split");
        CHECK(plain.iterations > 2, "but the stiff edge makes it walk");

        OneDNetwork fastNet(makeNet());
        AnalyticElement qn1(2.0e6, 0.0), qn2(r0, k);
        ctl.accelerate = true;
        const coupling::CoupledResult qn =
            coupling::solveCoupled1d3d(fastNet, {1, 2}, {&qn1, &qn2}, ctl, netCtl());
        CHECK(qn.converged, "the accelerated two-edge loop converges");
        CHECK(qn.iterations <= 10 && plain.iterations > 20 * qn.iterations,
              "with the whole interface vector accelerated together");
        CHECK(relNear(qn.state[1].q, q2, 1e-6), "the stiff edge takes its closed-form rate");
        CHECK(relNear(qn.state[1].resistance, r2, 1e-6), "and its closed-form coefficient");
        CHECK(relNear(qn.state[0].q + qn.state[1].q, Q, 1e-9),
              "the two edges still carry the injected rate between them");
        CHECK(relNear(qn.state[0].dpPa, qn.state[1].dpPa, 1e-7),
              "and still agree on one drop across the same two nodes");
        CHECK(qn1.calls() == qn.iterations && qn2.calls() == qn.iterations,
              "one network solve and one solve per partner, as the plain loop promises");
    }

    // ---- the same acceleration inside every time step ----
    {
        // the bypass here is 1e6, so this k puts the interface map's derivative at
        // -0.5: the plain loop needs about twenty passes per step to reach 1e-6
        const double qIn = 0.02, dt = 0.05, endTime = 1.0;
        // k puts the same map's derivative at -0.82: the plain loop needs dozens of
        // passes per step for its 1e-6, the secant a handful
        const double r0 = 1.0e6, k = 1.12e10;

        coupling::TransientCoupledControls ctl;
        ctl.tolerance = 1e-6;
        ctl.relax = 1.0;
        ctl.maxSubIterations = 120;

        TransientNetwork slowNet(transientBranched(qIn, 1.0e6));
        AnalyticElement slowBlock(r0, k);
        const coupling::TransientCoupledResult plain =
            coupling::solveTransientCoupled1d3d(slowNet, {1}, {&slowBlock}, ctl, dt, endTime);

        long plainPasses = 0;
        for (const coupling::TransientCoupledStep& st : plain.steps)
            plainPasses += st.passes;

        TransientNetwork fastNet(transientBranched(qIn, 1.0e6));
        AnalyticElement fastBlock(r0, k);
        ctl.accelerate = true;
        ctl.maxSubIterations = 6;
        const coupling::TransientCoupledResult qn =
            coupling::solveTransientCoupled1d3d(fastNet, {1}, {&fastBlock}, ctl, dt, endTime);

        long qnPasses = 0;
        bool everyStepAgreed = qn.steps.size() == plain.steps.size();
        for (const coupling::TransientCoupledStep& st : qn.steps)
        {
            qnPasses += st.passes;
            for (const coupling::EdgeExchange& ex : st.edges)
                everyStepAgreed = everyStepAgreed && ex.residual <= ctl.tolerance;
        }
        CHECK(everyStepAgreed, "every accelerated step is one the interface agreed on");
        CHECK(qnPasses < plainPasses,
              "the secant step needs fewer passes per step than the plain one");

        // the pass budget, not just the count: with eight passes per step the
        // plain loop cannot agree somewhere, and the accelerated one does
        TransientNetwork tightNet(transientBranched(qIn, 1.0e6));
        AnalyticElement tightBlock(r0, k);
        coupling::TransientCoupledControls tight = ctl;
        tight.accelerate = false;
        tight.maxSubIterations = 8;
        bool plainStalled = false;
        try
        {
            coupling::solveTransientCoupled1d3d(tightNet, {1}, {&tightBlock},
                                               tight, dt, endTime);
        }
        catch (const Error&) { plainStalled = true; }
        CHECK(plainStalled, "eight passes per step is not enough for the plain map");
        CHECK(qn.steps.size() == 20, "the accelerated arm finishes the run on those eight");
        CHECK(fastBlock.calls() == qnPasses,
              "the transient acceleration costs no partner solve beyond one per pass");
        CHECK(relNear(qn.steps.back().edges[0].q, plain.steps.back().edges[0].q, 1e-4),
              "both arms end on the same split at the same time");
    }

    std::cout << (failures == 0 ? "coupling unit: all checks passed\n"
                                : "coupling unit: FAILURES\n");
    return failures == 0 ? 0 : 1;
}
