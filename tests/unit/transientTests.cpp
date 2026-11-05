// Unit tests for the M2 transient 1D kernel: exact hydrate kinetics,
// implicit-Euler momentum/compliance solve against a hand-rolled 2x2
// reference, Newton quality, and liquid-inventory conservation.
#include <cmath>
#include <iostream>

#include "cf/1d/Transient.h"
#include "cf/common/Error.h"

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

    // ---- a coupled pipe charges its partner's coefficient, not its own wall ----
    {
        TransientCase c;
        c.fluid = FluidProperties{800.0, 0.1};
        TransientNodeConfig a, b;
        a.name = "Up"; a.pressureFixed = true; a.pressure = 1.0e5;
        b.name = "Dn"; b.pressureFixed = true; b.pressure = 0.0;
        TransientPipeConfig p;
        p.name = "X"; p.fromName = "Up"; p.toName = "Dn";
        p.length = 10.0; p.diameter = 0.02;
        p.coupledResistance = 1.0e6;
        c.nodes = {a, b}; c.pipes = {p};

        const double area = 0.25 * M_PI * p.diameter * p.diameter;
        const double inertia = c.fluid.rho * p.length / area;   // 2.5465e7
        const double dt = 0.05;
        const double qExact = a.pressure / (inertia / dt + p.coupledResistance);

        TransientNetwork sim(c);
        CHECK(sim.pipeCoupled(0) && sim.coupledResistanceOf(0) == 1.0e6,
              "a coupled pipe reports its coefficient");
        const double newton = sim.step(dt);
        CHECK(newton < 1e-11, "coupled momentum row solves to Newton tolerance");
        // this geometry's own Hagen-Poiseuille resistance is 2.5e8 -- 250 times
        // the partner's 1e6 -- so a pipe that also charged its wall would run at
        // 1.3e-4 here, and at 3.9e-4 once its inertia is forgotten as well: both
        // are about a factor two away from the answer above. The agreement band
        // is the kernel's own: the Newton loop stops at a 1e-11 scaled residual,
        // which on a row dominated by inertia/dt is a ~1e-10 relative rate.
        CHECK(relNear(sim.pipeFlows()[0], qExact, 1e-10),
              "coupled rate follows R*Q plus inertance, not the pipe's friction");
    }

    // ---- a coupled pipe cannot also carry a valve ----
    {
        TransientCase c;
        c.fluid = FluidProperties{800.0, 0.1};
        TransientNodeConfig a, b;
        a.name = "Up"; a.pressureFixed = true; a.pressure = 1.0e5;
        b.name = "Dn"; b.pressureFixed = true; b.pressure = 0.0;
        TransientPipeConfig p;
        p.name = "X"; p.fromName = "Up"; p.toName = "Dn";
        p.length = 10.0; p.diameter = 0.02;
        p.coupledResistance = 1.0e6;
        p.valveK = 2.0;
        c.nodes = {a, b}; c.pipes = {p};
        bool rejected = false;
        try { TransientNetwork sim(c); }
        catch (const Error&) { rejected = true; }
        CHECK(rejected, "constructing a coupled pipe with a valve is refused");
    }

    // ---- solveStep leaves the state alone, commitStep adopts it ----
    {
        TransientCase c;
        c.fluid = FluidProperties{800.0, 0.1};
        TransientNodeConfig a, b;
        a.name = "Up"; a.pressureFixed = true; a.pressure = 1.0e5;
        b.name = "Dn"; b.pressureFixed = true; b.pressure = 0.0;
        TransientPipeConfig p;
        p.name = "X"; p.fromName = "Up"; p.toName = "Dn";
        p.length = 10.0; p.diameter = 0.02;
        p.coupledResistance = 1.0e6;
        c.nodes = {a, b}; c.pipes = {p};

        const double dt = 0.05;
        const double area = 0.25 * M_PI * p.diameter * p.diameter;
        const double inertia = c.fluid.rho * p.length / area;

        TransientNetwork sim(c);
        const double res1 = sim.solveStep(dt);
        const double q1 = sim.pendingPipeFlow(0);
        const double res2 = sim.solveStep(dt);
        const double q2 = sim.pendingPipeFlow(0);
        CHECK(q1 == q2 && res1 == res2 && res1 < 1e-11
              && sim.time() == 0.0 && sim.pipeFlows()[0] == 0.0,
              "re-solving a step from the same start gives the same state");
        sim.commitStep();
        CHECK(sim.pipeFlows()[0] == q2 && relNear(sim.time(), dt, 1e-15),
              "commitStep adopts the pending rate and the clock");
        bool refused = false;
        try { sim.commitStep(); }
        catch (const Error&) { refused = true; }
        CHECK(refused, "committing twice is refused");

        // a coefficient change invalidates the pending step: it was solved
        // against the partner's old statement, not this one
        sim.solveStep(dt);
        sim.setCoupledResistance(0, 2.0e6);
        bool dropped = false;
        try { sim.commitStep(); }
        catch (const Error&) { dropped = true; }
        CHECK(dropped, "refreshing the coefficient drops the pending step");

        TransientNetwork joined(c);
        joined.setCoupledResistance(0, 2.0e6);
        joined.solveStep(0.1);
        CHECK(relNear(joined.pendingPipeFlow(0), 1.0e5 / (inertia / 0.1 + 2.0e6), 1e-10),
              "a re-solved step carries the refreshed coefficient");
        CHECK(joined.pendingPressure(0) == 1.0e5 && joined.pendingPressure(1) == 0.0,
              "Dirichlet nodes read their own pressure at every time level");
    }

    if (failures)
    {
        std::cout << failures << " transient test(s) FAILED\n";
        return 1;
    }
    std::cout << "all 1d transient unit tests passed\n";
    return 0;
}
