// Unit tests for the 1D steady pipe-network kernel: Churchill friction,
// single-pipe Hagen-Poiseuille behaviour, gravity head, and damped-Newton
// network solutions against closed-form series/parallel results.
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

#include "cf/1d/Network.h"

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

static const double kPi = 3.14159265358979323846;

// laminar network resistance R = 128 mu L / (pi D^4), dp = R * Q
int main()
{
    // ---- Churchill friction: laminar limit is exactly 64/Re ----
    {
        const double re = 433.0;
        const double f = churchillFriction(re, 0.0);
        CHECK(relNear(f, 64.0 / re, 1e-9), "churchill laminar limit 64/Re");
        // turbulent, smooth wall: f should be in the Blasius ballpark
        const double ft = churchillFriction(1e5, 0.0);
        CHECK(ft > 0.012 && ft < 0.020, "churchill turbulent range");
        // rough wall raises friction
        CHECK(churchillFriction(1e5, 0.01) > ft, "roughness increases f");
    }

    const FluidProperties fluid{850.0, 0.1};
    const double D = 0.05;
    const double area = 0.25 * kPi * D * D;
    auto resistance = [D](double L)
    { return 128.0 * 0.1 * L / (kPi * D * D * D * D); };

    // ---- single pipe: dp -> Q reproduces Hagen-Poiseuille ----
    {
        NetworkPipe p;
        p.name = "P"; p.fromName = "A"; p.toName = "B";
        p.length = 1000; p.diameter = D; p.roughness = 0; p.elevation = 0;
        const double Q = 0.002;
        PipeFlowResult fr;
        solvePipeFlow(p, fluid, resistance(1000.0) * Q, fr);
        CHECK(relNear(fr.flowRate, Q, 1e-9), "laminar Q from dp (Hagen-Poiseuille)");
        CHECK(fr.velocity > 0 && relNear(fr.velocity * area, Q, 1e-9), "Q = v*A");
        CHECK(fr.reynolds < 2300, "test flow is laminar");
        // reverse dp -> reverse flow (odd symmetry of the friction law)
        solvePipeFlow(p, fluid, -resistance(1000.0) * Q, fr);
        CHECK(relNear(fr.flowRate, -Q, 1e-9), "reverse dp gives reverse Q");
    }

    // ---- gravity head: dp == rho g dz is pure hydrostatic, Q = 0 ----
    {
        NetworkPipe p;
        p.name = "P"; p.fromName = "A"; p.toName = "B";
        p.length = 100; p.diameter = D; p.roughness = 0; p.elevation = 5;
        const double dpHydro = fluid.rho * 9.80665 * 5.0;
        PipeFlowResult fr;
        solvePipeFlow(p, fluid, dpHydro, fr);
        CHECK(std::abs(fr.flowRate) < 1e-9, "hydrostatic balance -> no flow");
        CHECK(std::abs(fr.dpFriction) < 1e-3, "gravity head removed from friction dp");
    }

    // ---- one well + tubing + fixed terminal: Newton vs closed form ----
    {
        OneDNetwork net;
        net.fluid = fluid;
        net.addNode("W");
        net.addNode("Outlet");
        net.setNodePressure("Outlet", 1.0e6);
        net.setNodeRate("W", 0.002);
        NetworkPipe p;
        p.name = "Tubing"; p.fromName = "W"; p.toName = "Outlet";
        p.length = 1000; p.diameter = D; p.roughness = 0; p.elevation = 5;
        net.addPipe(p);

        const SteadyNetworkResult res = solveSteadyNetwork(net, SteadyControls{});
        CHECK(res.converged, "single-pipe network converged");
        const double pOutExpected = 1.0e6
            + resistance(1000.0) * 0.002 + fluid.rho * 9.80665 * 5.0;
        CHECK(relNear(res.nodePressure[0], pOutExpected, 1e-6), "well head pressure");
        CHECK(relNear(res.pipeFlows[0].flowRate, 0.002, 1e-6), "mass balance through tubing");
    }

    // ---- symmetric Y: injected flow splits equally across two identical
    //      parallel legs ----
    {
        OneDNetwork net;
        net.fluid = fluid;
        net.addNode("S");      // free junction, fed 0.002
        net.addNode("O1");     // fixed terminals
        net.addNode("O2");
        net.setNodeRate("S", 0.002);
        net.setNodePressure("O1", 1.0e6);
        net.setNodePressure("O2", 1.0e6);
        for (int k = 0; k < 2; ++k)
        {
            NetworkPipe p;
            p.name = "L" + std::to_string(k);
            p.fromName = "S"; p.toName = k == 0 ? "O1" : "O2";
            p.length = 500; p.diameter = D; p.roughness = 0; p.elevation = 0;
            net.addPipe(p);
        }
        const SteadyNetworkResult res = solveSteadyNetwork(net, SteadyControls{});
        CHECK(res.converged, "Y-split converged");
        CHECK(relNear(res.pipeFlows[0].flowRate, 0.001, 1e-6), "leg 0 carries half");
        CHECK(relNear(res.pipeFlows[1].flowRate, 0.001, 1e-6), "leg 1 carries half");
        const double pS = 1.0e6 + resistance(500.0) * 0.001;
        CHECK(relNear(res.nodePressure[0], pS, 1e-6), "junction pressure");
    }

    // ---- dictionary round trip: load the shipped case ----
    {
        const Dictionary caseDict = Dictionary::parseFile(
            CF_CASE_DIR "/steadyNetwork/case.dict");
        const OneDNetwork net = loadNetworkFromDict(caseDict);
        CHECK(net.nNodes() == 4, "case has 4 nodes");
        CHECK(net.pipes.size() == 3, "case has 3 pipes");
        CHECK(relNear(net.fluid.mu, 0.1, 1e-12) && relNear(net.fluid.rho, 850, 1e-12),
              "case fluid parsed");
        CHECK(net.pipes[0].elevation == 5, "case pipe elevation parsed");
        CHECK(relNear(net.nodeExternalRate[net.nodeIndex("W1")], 0.002, 1e-12),
              "case source parsed");
        const SteadyControls ctl = loadSteadyControls(caseDict);
        CHECK(ctl.maxIter == 100 && ctl.tol < 1e-9, "case solve controls parsed");
        const SteadyNetworkResult res = solveSteadyNetwork(net, ctl);
        CHECK(res.converged, "case network converged");
        double balance = 0;
        for (std::size_t i = 0; i < net.nNodes(); ++i)
            balance += net.nodeExternalRate[i];
        double out = 0;
        // terminal node receives the merged flow
        for (const PipeFlowResult& f : res.pipeFlows)
            if (f.name == "P3") out = f.flowRate;
        CHECK(relNear(out, 0.003, 1e-6), "case: terminal carries 0.003 m3/s");
        CHECK(relNear(balance, 0.003, 1e-12), "case: total injection is 0.003 m3/s");
    }

    if (failures)
    {
        std::cout << failures << " network test(s) FAILED\n";
        return 1;
    }
    std::cout << "all 1d network unit tests passed\n";
    return 0;
}
