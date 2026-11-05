// Unit tests for the partitioned 1D<->3D interface loop, driven against
// analytic partner models so the loop itself -- not the 3D solve behind it --
// is what fails when it is wrong.
//
// The three things these cases pin down:
//   * the loop lands on the state where the drop the network charges equals
//     the drop the partner reports at the rate the network delivers;
//   * for a linear partner one coefficient refresh is enough, and the answer
//     does not depend on the starting coefficient or on the relaxation;
//   * the numbers a caller reads out of CoupledResult describe ONE state:
//     the rate the partner was last asked for is the rate being reported.
#include <cmath>
#include <iostream>
#include <string>

#include "cf/1d/Network.h"
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
            coupling::solveCoupled1d3d(net, edge, block, ctl, nc);

        CHECK(res.converged, "linear partner converges");
        CHECK(res.iterations == 2, "linear partner needs one coefficient refresh");
        CHECK(relNear(res.q, 0.025, 1e-12), "rate is dp/R of the partner");
        CHECK(relNear(res.resistance, 4.0e6, 1e-12), "coefficient is the partner's");
        CHECK(relNear(net.pipes[edge].coupledResistance, res.resistance, 1e-15),
              "network carries the reported coefficient");
        CHECK(relNear(block.lastQ(), res.q, 1e-15),
              "the partner was measured at the reported rate");
        CHECK(relNear(block.lastDp(), res.dpPa, 1e-15),
              "the reported drop is the partner's own answer");
        CHECK(relNear(res.dpPa, res.resistance * res.q, 1e-10),
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
            coupling::solveCoupled1d3d(damped, e1, b1, ctl, nc);

        OneDNetwork quick = twoNodeNetwork(2.0e5, 1.0e5);
        const label e2 = quick.addPipe(linearPipe("X", "Src", "Snk", 4.0e6));
        AnalyticElement b2(4.0e6, 0.0);
        coupling::CoupledControls ctl2 = ctl;
        ctl2.relax = 1.0;
        const coupling::CoupledResult r2 =
            coupling::solveCoupled1d3d(quick, e2, b2, ctl2, nc);

        CHECK(r1.converged && r2.converged, "both interface schemes converge");
        CHECK(relNear(r1.q, r2.q, 1e-10), "relaxation does not move the fixed point");
        CHECK(relNear(r1.resistance, 4.0e6, 1e-10), "relaxed coefficient is the same");
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
            coupling::solveCoupled1d3d(net, edge, block, ctl, nc);

        CHECK(res.converged, "quadratic partner converges");
        CHECK(relNear(res.q, qExact, 1e-9), "nonlinear fixed point is the exact root");
        CHECK(res.residual < 1e-10, "interface residual reaches the tolerance");
        CHECK(relNear(block.lastQ(), res.q, 1e-15),
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
            coupling::solveCoupled1d3d(net, edge, block, ctl, nc);

        // current divider between the partner and the 1D bypass
        const double qBlockExact = qIn * rBy / (rBlock + rBy);
        const double dpBranchExact = qBlockExact * rBlock;
        const double pSnk = 1.0e5;
        const double pJ2Exact = pSnk + qIn * rLead;
        const double pJ1Exact = pJ2Exact + dpBranchExact;
        const double pSrcExact = pJ1Exact + qIn * rLead;

        CHECK(res.converged, "branched network converges");
        CHECK(relNear(res.q, qBlockExact, 1e-10), "coupled edge carries its branch share");
        CHECK(relNear(res.network.nodePressure[0], pSrcExact, 1e-10),
              "source node pressure matches the closed form");
        CHECK(relNear(res.network.nodePressure[1], pJ1Exact, 1e-10),
              "upstream interface node pressure matches");
        CHECK(relNear(res.network.nodePressure[2], pJ2Exact, 1e-10),
              "downstream interface node pressure matches");

        const auto& pf = res.network.pipeFlows;
        CHECK(relNear(pf[0].flowRate, qIn, 1e-10), "the lead takes the whole rate");
        CHECK(relNear(pf[0].flowRate, res.q + pf[2].flowRate, 1e-10),
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
            coupling::solveCoupled1d3d(net, edge, block, ctl, nc);

        CHECK(!res.converged, "a stopped loop says it stopped");
        CHECK(res.iterations == 3, "every allowed pass ran");
        CHECK(res.history.size() == 3, "each pass is in the history");
        CHECK(relNear(block.lastQ(), res.q, 1e-15),
              "the stopped state is still one state");
        CHECK(relNear(block.lastDp(), res.dpPa, 1e-15),
              "the reported drop is the partner's own last answer");
        CHECK(relNear(res.residual,
                      std::abs(res.resistance * res.q - res.dpPa) / res.dpPa, 1e-10),
              "the residual is the gap between the two solvers");
    }

    std::cout << (failures == 0 ? "coupling unit: all checks passed\n"
                                : "coupling unit: FAILURES\n");
    return failures == 0 ? 0 : 1;
}
