// Unit tests for the M2 transient 1D kernel: exact hydrate kinetics,
// implicit-Euler momentum/compliance solve against a hand-rolled 2x2
// reference, Newton quality, and liquid-inventory conservation.
#include <cmath>
#include <iostream>

#include "cf/1d/Transient.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }\
        else         { std::cout << "ok:   " << (msg) << "\n"; }             \
    } while (0)

static bool relNear(double a, double b, double tol = 1e-9)
{
    const double s = std::max(std::abs(a), std::abs(b));
    if (s == 0) return true;
    return std::abs(a - b) / s < tol;
}

int main()
{
    // ---- hydrate kinetics: dh/dt = k(1-h) is integrated exactly ----
    {
        double vL = 1.0, vH = 0.0;
        const double k = 0.5;
        double t = 0;
        for (int i = 0; i < 7; ++i)
        {
            hydrateKineticsStep(vL, vH, k, 0.3);
            t += 0.3;
            CHECK(relNear(vL + vH, 1.0, 1e-12), "kinetics conserves liquid volume");
            CHECK(relNear(vL, std::exp(-k * t), 1e-12),
                  "kinetics matches exact exponential");
        }
        CHECK(relNear(vH, 1.0 - std::exp(-k * t), 1e-12), "hydrate fraction exact");
    }

    // ---- one implicit-Euler step vs hand-rolled 2x2 laminar reference ----
    {
        TransientCase c;
        c.fluid = FluidProperties{850.0, 0.1};
        TransientNodeConfig tank;
        tank.name = "TANK"; tank.volume = 1.0; tank.bulkModulus = 2.0e8;
        tank.pressure = 5.0e5; tank.liquidFraction = 1.0;
        TransientNodeConfig out;
        out.name = "OUT"; out.pressureFixed = true; out.pressure = 0.0;
        TransientPipeConfig pipe;
        pipe.name = "LINE"; pipe.fromName = "TANK"; pipe.toName = "OUT";
        pipe.length = 1000; pipe.diameter = 0.02;
        c.nodes = {tank, out};
        c.pipes = {pipe};

        TransientNetwork sim(c);
        const double dt = 0.05;
        const double res = sim.step(dt);
        CHECK(res < 1e-9, "Newton residual small on the linear laminar step");

        const double area = 0.25 * M_PI * 0.02 * 0.02;
        const double Rr = 128.0 * 0.1 * 1000.0 / (M_PI * 0.02 * 0.02 * 0.02 * 0.02);
        const double inertia = 850.0 * 1000.0 / area;
        // simultaneous implicit Euler: (I/dt + R + dt*K/V) q = p_prev
        const double q = 5.0e5 / (inertia / dt + Rr + dt * 2.0e8 / 1.0);
        const double p = 5.0e5 - dt * q * 2.0e8 / 1.0;
        CHECK(relNear(sim.pipeFlows()[0], q, 1e-10), "implicit-Euler Q matches 2x2 reference");
        CHECK(relNear(sim.pressure()[0], p, 1e-10), "implicit-Euler p matches 2x2 reference");
    }

    // ---- closed system: identical pressures -> nothing moves ----
    {
        TransientCase c;
        c.fluid = FluidProperties{850.0, 0.1};
        TransientNodeConfig a, b;
        a.name = "A"; b.name = "B";
        for (auto* n : {&a, &b})
        {
            n->volume = 1.0; n->bulkModulus = 1.0e8; n->pressure = 2.0e5;
            n->liquidFraction = 0.5; n->hydrateFraction = 0.2;
        }
        TransientPipeConfig p;
        p.name = "P"; p.fromName = "A"; p.toName = "B";
        p.length = 10; p.diameter = 0.02;
        c.nodes = {a, b}; c.pipes = {p};

        TransientNetwork sim(c);
        for (int i = 0; i < 20; ++i) sim.step(0.01);
        CHECK(relNear(sim.pipeFlows()[0], 0.0, 1e-9), "equilibrium stays at rest");
        const std::vector<double> vL = sim.liquidVolumes(), vH = sim.hydrateVolumes();
        CHECK(relNear(vL[0] + vH[0], 0.5, 1e-12) && relNear(vL[1] + vH[1], 0.5, 1e-12),
              "no phase change without hydrate rate");
    }

    // ---- hydrate rate shrinks the mobile fraction over a run ----
    {
        TransientCase c;
        c.fluid = FluidProperties{850.0, 0.1};
        TransientNodeConfig a;
        a.name = "A"; a.volume = 1.0; a.bulkModulus = 1.0e8;
        a.pressure = 2.0e5; a.liquidFraction = 1.0; a.hydrateRate = 0.1;
        TransientNodeConfig b = a;
        b.name = "B"; b.hydrateRate = 0.0;
        TransientPipeConfig p;
        p.name = "P"; p.fromName = "A"; p.toName = "B";
        p.length = 10; p.diameter = 0.02;
        c.nodes = {a, b}; c.pipes = {p};

        TransientNetwork sim(c);
        for (int i = 0; i < 10; ++i) sim.step(0.05);
        const double alphaA = sim.liquidFractions()[0];
        CHECK(relNear(alphaA, std::exp(-0.1 * 0.5), 1e-6),
              "mobile liquid decays exponentially under hydrate rate");
    }

    // ---- dictionary loader ----
    {
        const Dictionary caseDict = Dictionary::parseFile(
            CF_CASE_DIR "/transientSlug/case.dict");
        const TransientCase c = loadTransientFromDict(caseDict);
        const TransientTimeControls tc = loadTransientControls(caseDict);
        CHECK(c.nodes.size() == 6 && c.pipes.size() == 5, "slug case parsed");
        CHECK(tc.dt == 5000.0 && tc.reportEvery == 1, "slug time controls parsed");
        const TransientSeries s = runTransient(c, tc.dt, tc.dt * 4, 1);
        CHECK(s.samples.size() == 5, "runTransient sample count");
        CHECK(s.samples[4].liquidFraction[1] == 1.0, "slug: source node stays full");
        // The incompressible ports force one series flow; the nodes are
        // well-mixed, so the transport is exactly the discrete upwind
        // recursion alpha_i += (alpha_{i-1} - alpha_i) * Q dt / V.
        const double Rr = 128.0 * 0.1 * 100.0 / (M_PI * 1.6e-7);
        const double qAnalytic = 6.0e5 / (5.0 * Rr);
        double a2 = 0, a3 = 0, a4 = 0;
        for (std::size_t n = 1; n < s.samples.size(); ++n)
        {
            const double q = s.samples[n].flowRate[0];
            // the first step is still inside the inertial ramp
            CHECK(relNear(q, qAnalytic, n == 1 ? 1e-4 : 1e-9),
                  "slug: series flow matches the laminar resistance");
            const double d = q * tc.dt / 1.0;
            const double n2 = a2 + d * (1.0 - a2);
            const double n3 = a3 + d * (a2 - a3);
            const double n4 = a4 + d * (a3 - a4);
            a2 = n2; a3 = n3; a4 = n4;
        }
        CHECK(relNear(s.samples[4].liquidFraction[2], a2, 1e-9),
              "slug: node 2 fill matches the discrete CSTR recursion");
        CHECK(relNear(s.samples[4].liquidFraction[3], a3, 1e-9),
              "slug: node 3 fill matches the cascade");
        CHECK(relNear(s.samples[4].liquidFraction[4], a4, 1e-9),
              "slug: node 4 fill matches the cascade");
    }

    if (failures)
    {
        std::cout << failures << " transient test(s) FAILED\n";
        return 1;
    }
    std::cout << "all 1d transient unit tests passed\n";
    return 0;
}
