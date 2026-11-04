// Unit tests for the 1D steady pipe-network kernel: Churchill friction,
// single-pipe Hagen-Poiseuille behaviour, gravity head, and damped-Newton
// network solutions against closed-form series/parallel results.
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>

#include "cf/1d/BeggsBrill.h"
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

    // ---- Beggs-Brill: golden anchor cross-checked against the fluids
    //      library (Beggs_Brill(0.6, 0.1, 915, 2.67, 180e-6, 14e-6, 0.0487,
    //      1e7, 0.05, 0deg, 0, 1) -> 686.9724506803469 Pa/m) ----
    {
        const double D = 0.05;
        const double m = 0.6, x = 0.1, rhol = 915.0, rhog = 2.67;
        const double area = 0.25 * kPi * D * D;
        BeggsBrillIn in;
        in.rhoL = rhol;
        in.rhoG = rhog;
        in.muL = 180e-6;
        in.muG = 14e-6;
        in.sigma = 0.0487;
        in.vsl = (1.0 - x) * m / rhol / area;
        in.vsg = x * m / rhog / area;
        in.pPa = 1.0e7;
        in.diameter = D;
        in.roughness = 0.0;
        in.angleDeg = 0.0;
        const BeggsBrillOut out = beggsBrillDpDs(in);
        CHECK(relNear(out.dpPerMeter, 686.9724506803469, 1e-6),
              "Beggs-Brill golden anchor (fluids cross-check)");
        CHECK(out.regime == 3, "anchor flow is distributed");
        CHECK(out.holdup > out.lambdaL, "slippage raises the holdup");
    }

    // ---- degenerate limits: pure liquid / pure gas collapse to Darcy-Weisbach ----
    {
        BeggsBrillIn in;
        in.rhoL = 1000.0;
        in.rhoG = 2.0;
        in.muL = 1e-3;
        in.muG = 1e-5;
        in.sigma = 0.07;
        in.vsl = 2.0;
        in.vsg = 0.0;
        in.pPa = 5.0e6;
        in.diameter = 0.1;
        in.roughness = 0.0;
        in.angleDeg = 0.0;
        // lambda = 1 clamps the holdup to 1, so the two-phase multiplier
        // S = ln(lambda/HL^2) vanishes and only DW friction remains
        const double re = 1000.0 * 2.0 * 0.1 / 1e-3;
        const double dw = pvt::colebrookDarcyFriction(re, 0.0) * 1000.0 * 4.0
                        / (2.0 * 0.1);
        const BeggsBrillOut out = beggsBrillDpDs(in);
        CHECK(relNear(out.dpPerMeter, dw, 1e-12),
              "pure liquid degenerates to Darcy-Weisbach");
        CHECK(out.holdup == 1.0 && out.regime == 3, "pure liquid fills the pipe");
        in.angleDeg = 90.0; // distributed regime takes no inclination lift
        const BeggsBrillOut vert = beggsBrillDpDs(in);
        CHECK(relNear(vert.dpPerMeter, dw + 1000.0 * 9.80665, 1e-12),
              "vertical pure liquid adds the exact gravity head");

        // vsl = 0: lambda = 0, holdup 0, no slippage correction, gas DW
        in.rhoL = 800.0;
        in.rhoG = 50.0;
        in.muL = 5e-3;
        in.muG = 1.2e-5;
        in.sigma = 0.02;
        in.vsl = 0.0;
        in.vsg = 8.0;
        in.diameter = 0.05;
        in.angleDeg = 0.0;
        const double reG = 50.0 * 8.0 * 0.05 / 1.2e-5;
        const double ek = 8.0 * 8.0 * 50.0 / 5.0e6;
        const double dwG = pvt::colebrookDarcyFriction(reG, 0.0) * 50.0 * 64.0
                         / (2.0 * 0.05) / (1.0 - ek);
        const BeggsBrillOut gas = beggsBrillDpDs(in);
        CHECK(relNear(gas.dpPerMeter, dwG, 1e-12),
              "pure gas degenerates to gas Darcy-Weisbach");
        CHECK(gas.holdup == 0.0 && gas.regime == 3, "pure gas: zero holdup");
    }

    // ---- march: dry oil (no free gas) keeps a liquid-full column ----
    {
        pvt::BlackOilFluid f;
        f.gasGravity = 0.75;
        f.oilDensityStd = 900.0;
        f.gorStd = 0.0;
        f.temperature = 320.0;
        f.surfaceTension = 0.03;
        const double L = 100.0, D = 0.1;
        const BeggsBrillMarchOut m =
            marchBeggsBrillPipe(f, 0.02, 4.0e6, L, D, 0.0, 0.0, 24);
        CHECK(m.ok, "dry-oil march clears the pressure floor");
        CHECK(m.holdupOutlet == 1.0, "dry oil is liquid full");
        CHECK(relNear(m.massRate, 0.02 * pvt::blackOilMassPerOilSm3(f), 1e-12),
              "march mass rate is stock-tank mass");
        // the undersaturated-Bo drift over 100 m stays sub-mille, so the
        // march total sits on the inlet gradient times length
        const pvt::BlackOilState st = pvt::evaluateBlackOil(f, 4.0e6);
        double rhoL = 0, muL = 0;
        pvt::blackOilLiquid(f, st, rhoL, muL);
        BeggsBrillIn in;
        in.rhoL = rhoL;
        in.rhoG = st.rhoGas;
        in.muL = muL;
        in.muG = st.muGas;
        in.sigma = f.surfaceTension;
        in.vsl = 0.02 * st.bo / (0.25 * kPi * D * D);
        in.vsg = 0.0;
        in.pPa = 4.0e6;
        in.diameter = D;
        in.roughness = 0.0;
        in.angleDeg = 0.0;
        const BeggsBrillOut g = beggsBrillDpDs(in);
        CHECK(relNear(m.dpTotal, g.dpPerMeter * L, 1e-3),
              "incompressible march tracks the inlet gradient");
    }

    // ---- uphill march: slip holds liquid above its no-slip share ----
    {
        pvt::BlackOilFluid f;
        f.gasGravity = 0.75;
        f.oilDensityStd = 850.0;
        f.gorStd = 60.0;
        f.temperature = 340.0;
        f.surfaceTension = 0.025;
        const BeggsBrillMarchOut m =
            marchBeggsBrillPipe(f, 0.005, 3.0e6, 500.0, 0.076, 4.6e-5, 1.0, 24);
        CHECK(m.ok, "uphill march clears the floor");
        CHECK(m.holdupOutlet >= m.lambdaOutlet - 1e-12,
              "uphill holdup is at least the no-slip fraction");
        CHECK(m.holdupOutlet > 0.0 && m.holdupOutlet <= 1.0,
              "holdup within (0, 1]");
        CHECK(m.regimeMin >= 0 && m.regimeMax <= 3, "regimes in range");
    }

    // ---- two-phase pipe law: bracketed-secant inversion closes the ends ----
    {
        pvt::BlackOilFluid f;
        f.gasGravity = 0.75;
        f.oilDensityStd = 850.0;
        f.gorStd = 60.0;
        f.temperature = 340.0;
        f.surfaceTension = 0.025;
        NetworkPipe p;
        p.name = "FLOWLINE";
        p.fromName = "A";
        p.toName = "B";
        p.length = 1000;
        p.diameter = 0.076;
        p.roughness = 4.6e-5;
        p.elevation = 50;
        PipeFlowResult out;
        const bool okf = solvePipeFlowTwoPhase(p, f, 5.0e6, 4.2e6, out);
        CHECK(okf && !out.choked, "two-phase pipe solves without choking");
        CHECK(out.flowRate > 0.0, "flow runs down the pressure gradient");
        CHECK(relNear(out.massRate,
                      out.flowRate * pvt::blackOilMassPerOilSm3(f), 1e-12),
              "pipe mass rate consistent with the stock-tank rate");
        // re-marching at the solved rate must land on the imposed outlet
        const BeggsBrillMarchOut m =
            marchBeggsBrillPipe(f, out.flowRate, 5.0e6, 1000.0, 0.076, 4.6e-5,
                                0.05, 24);
        CHECK(std::abs(m.pOutlet - 4.2e6) < 5.0,
              "secant inversion closes the imposed pressure drop");
        CHECK(out.conductance > 0.0, "positive conductance");
        CHECK(out.holdup > 0.0 && out.holdup <= 1.0, "holdup in range");
        CHECK(out.regime >= 0 && out.regime <= 3, "regime in range");
        PipeFlowResult softer;
        solvePipeFlowTwoPhase(p, f, 5.0e6, 4.5e6, softer);
        CHECK(softer.flowRate < out.flowRate,
              "smaller pressure drop gives less flow");
    }

    // ---- two-phase network end to end ----
    {
        OneDNetwork net;
        net.hasBlackOil = true;
        net.blackOil.gasGravity = 0.75;
        net.blackOil.oilDensityStd = 850.0;
        net.blackOil.gorStd = 60.0;
        net.blackOil.temperature = 340.0;
        net.blackOil.surfaceTension = 0.025;
        net.addNode("W");
        net.addNode("Outlet");
        net.setNodePressure("Outlet", 4.0e6);
        net.setNodeRate("W", 0.004);
        NetworkPipe p;
        p.name = "TUBING";
        p.fromName = "W";
        p.toName = "Outlet";
        p.length = 800;
        p.diameter = 0.076;
        p.roughness = 4.6e-5;
        p.elevation = 30;
        net.addPipe(p);

        SteadyControls ctl;
        ctl.tol = 1e-10;
        ctl.maxIter = 60;
        const SteadyNetworkResult res = solveSteadyNetwork(net, ctl);
        CHECK(res.converged, "two-phase network converged");
        CHECK(relNear(res.pipeFlows[0].flowRate, 0.004, 1e-8),
              "two-phase mass balance through the tubing");
        CHECK(res.pipeFlows[0].regime >= 0, "pipe reports a two-phase regime");
        // the well pressure is exactly the outlet plus the march of the rate
        const BeggsBrillMarchOut m =
            marchBeggsBrillPipe(net.blackOil, 0.004, res.nodePressure[0], 800.0,
                                0.076, 4.6e-5, 30.0 / 800.0, 24);
        CHECK(std::abs(m.pOutlet - 4.0e6) < 5.0,
              "network pressure closes the Beggs-Brill march");
    }

    // ---- tabulated properties: flat table reproduces the constant law ----
    {
        PvtFluid tf;
        tf.pressure = {1.0e5, 5.0e6};
        tf.rho = {850.0, 850.0};
        tf.mu = {0.1, 0.1};
        CHECK(relNear(tf.density(2.0e6), 850.0, 1e-12)
                  && relNear(tf.viscosity(2.0e6), 0.1, 1e-12),
              "flat table interpolates to the constant");

        PvtFluid lin;
        lin.pressure = {1.0e6, 3.0e6};
        lin.rho = {800.0, 900.0};
        lin.mu = {1e-3, 2e-3};
        CHECK(relNear(lin.density(2.0e6), 850.0, 1e-12), "midpoint rho linear");
        CHECK(relNear(lin.viscosity(2.0e6), 1.5e-3, 1e-12), "midpoint mu linear");
        CHECK(relNear(lin.density(1.0e6), 800.0, 1e-12)
                  && relNear(lin.density(9.0e6), 900.0, 1e-12)
                  && relNear(lin.density(1.0e5), 800.0, 1e-12),
              "table clamps outside the grid");

        // flat table: the 24-segment property march must reproduce the
        // implicit constant-property pipe law (they are the same function
        // of q, each solved by a secant)
        NetworkPipe p;
        p.name = "P"; p.fromName = "A"; p.toName = "B";
        p.length = 1000; p.diameter = D; p.roughness = 0; p.elevation = 0;
        for (const double dpTest : {2.0e5, 5.0e5, 1.0e6})
        {
            PipeFlowResult fc, fp;
            solvePipeFlow(p, fluid, dpTest, fc);
            solvePipeFlowPvt(p, tf, 1.0e6 + dpTest, 1.0e6, fp);
            CHECK(relNear(fp.flowRate, fc.flowRate, 1e-9),
                  "flat PVT rate matches the constant-property rate");
            CHECK(relNear(fp.velocity, fc.velocity, 1e-9),
                  "flat PVT velocity matches (outlet = inlet when rho const)");
        }
        // reverse direction
        PipeFlowResult fr1, fr2;
        solvePipeFlow(p, fluid, -5.0e5, fr1);
        solvePipeFlowPvt(p, tf, 1.0e6, 1.5e6, fr2);
        CHECK(relNear(fr2.flowRate, fr1.flowRate, 1e-9),
              "flat PVT reverse flow matches");

        // elevation + valve + pump on the flat table: equipment paths agree
        NetworkPipe e = p;
        e.name = "EQ"; e.elevation = 30; e.valveK = 2.0;
        e.pumpHead = 40; e.pumpRate = 0.05; e.pumpShutoffHead = 45;
        PipeFlowResult ec, ep;
        solvePipeFlow(e, fluid, 6.0e5, ec);
        solvePipeFlowPvt(e, tf, 1.6e6, 1.0e6, ep);
        CHECK(relNear(ep.flowRate, ec.flowRate, 1e-9),
              "flat PVT equipment pipe matches (valve + pump + lift)");
        CHECK(relNear(ep.dpValve, ec.dpValve, 1e-9)
                  && relNear(ep.dpPump, ec.dpPump, 1e-9),
              "flat PVT equipment dp reports match");
    }

    // ---- isothermal gas column: property march vs in-test restatement ----
    {
        // rho = p/c (c = R T): the static column is p0 exp(-g H / c); the
        // march discretises the same ODE, so the test restates the Euler
        // recursion and the mass-consistent velocity law
        const double c = 1.0e5;
        const double H = 100.0, Dg = 0.05, muG = 1.2e-5;
        const double x = 9.80665 * H / c;
        PvtFluid gas;
        gas.pressure = {1.0e5, 2.0e6};
        gas.rho = {1.0e5 / c, 2.0e6 / c}; // linear between = rho = p/c
        gas.mu = {muG, muG};
        NetworkPipe p;
        p.name = "RISER"; p.fromName = "B"; p.toName = "T";
        p.length = H; p.diameter = Dg; p.roughness = 0; p.elevation = H;

        auto march = [&](double q, double pUp)
        {
            const double area = 0.25 * kPi * Dg * Dg;
            const double w = pUp / c * q;
            const double ds = H / 24.0;
            double pr = pUp;
            double vOut = 0, reOut = 0, fOut = 0;
            for (int i = 0; i < 24; ++i)
            {
                const double rho = pr / c;
                vOut = w / (rho * area);
                reOut = rho * vOut * Dg / muG;
                fOut = churchillFriction(reOut, 0.0);
                pr -= rho * 9.80665 * ds + fOut * rho * vOut * vOut
                    / (2.0 * Dg) * ds;
            }
            return std::make_tuple(pr, vOut);
        };

        // static rest: the column alone overshoots the imposed outlet ->
        // the law is flat and the rate is exactly zero
        PipeFlowResult rest;
        const double pUp = 2.0e6;
        solvePipeFlowPvt(p, gas, pUp,
                         pUp * std::exp(-x) * (1.0 + 1e-4), rest);
        CHECK(rest.flowRate == 0.0, "gas column at static surplus -> no flow");
        // the Euler static column tracks the exponential closed form
        // within the first-order truncation x^2/(2*24)
        CHECK(std::abs(std::get<0>(march(0.0, pUp)) - pUp * std::exp(-x))
                  < pUp * 1.1 * x * x / 48.0,
              "static gas column tracks exp(-gH/c) within Euler truncation");

        // live branch: outlet 1e-4 below the static column -> friction
        // supplies the deficit; invert the same recursion in-test
        const double pLow = pUp * std::exp(-x) * (1.0 - 1e-4);
        PipeFlowResult live;
        solvePipeFlowPvt(p, gas, pUp, pLow, live);
        CHECK(live.flowRate > 0.0,
              "deficit below the static column gives flow");
        double qLo = 0.0, fLo = std::get<0>(march(0.0, pUp)) - pLow;
        double qHi = 1e-4, fHi = std::get<0>(march(qHi, pUp)) - pLow;
        while (fHi > 0.0)
        {
            qLo = qHi; fLo = fHi;
            qHi *= 4.0;
            fHi = std::get<0>(march(qHi, pUp)) - pLow;
        }
        double qRef = qLo + (qHi - qLo) * fLo / (fLo - fHi);
        for (int it = 0; it < 100; ++it)
        {
            const double fq = std::get<0>(march(qRef, pUp)) - pLow;
            if (std::abs(fq) < 1e-12 * pUp) break;
            if (fq > 0.0) { qLo = qRef; fLo = fq; fLo *= 0.5; }
            else { qHi = qRef; fHi = fq; fHi *= 0.5; }
            const double qn = qLo + (qHi - qLo) * fLo / (fLo - fHi);
            qRef = (qn > qLo && qn < qHi) ? qn : 0.5 * (qLo + qHi);
        }
        CHECK(relNear(live.flowRate, qRef, 1e-9),
              "inverted rate matches the in-test gas-column re-march");
        // the reported velocity is the OUTLET velocity of the expanding
        // fill: strictly above the inlet velocity q/A
        CHECK(live.velocity > live.flowRate / (0.25 * kPi * Dg * Dg),
              "outlet velocity exceeds the inlet (gas expands)");
        CHECK(relNear(live.velocity, std::get<1>(march(qRef, pUp)), 1e-9),
              "reported velocity is the march outlet velocity");
        CHECK(!live.choked, "gas column does not choke");
    }

    // ---- PVT network end to end: free junction closed by the march ----
    {
        OneDNetwork net;
        net.hasPvtFluid = true;
        net.pvtFluid.pressure = {1.0e6, 2.0e6};
        net.pvtFluid.rho = {10.0, 20.0}; // both nodes on rho = p/1e5, so the
                                         // linear table IS rho = p/1e5
        net.pvtFluid.mu = {1.2e-5, 1.2e-5};
        net.addNode("W");
        net.addNode("T");
        net.setNodePressure("T", 1.5e6);
        net.setNodeRate("W", 0.002);
        NetworkPipe p;
        p.name = "FLOWLINE"; p.fromName = "W"; p.toName = "T";
        p.length = 2000; p.diameter = 0.1; p.roughness = 1e-5; p.elevation = 40;
        net.addPipe(p);

        SteadyControls ctl;
        ctl.tol = 1e-11;
        const SteadyNetworkResult res = solveSteadyNetwork(net, ctl);
        CHECK(res.converged, "PVT network converged");
        CHECK(relNear(res.pipeFlows[0].flowRate, 0.002, 1e-9),
              "PVT mass balance through the line");
        // the solved well pressure closes the mass-consistent march
        auto march = [&](double pW)
        {
            const double area = 0.25 * kPi * 0.1 * 0.1;
            const double w = pW / 1.0e5 * 0.002;
            const double ds = 2000.0 / 24.0;
            double pr = pW;
            for (int i = 0; i < 24; ++i)
            {
                const double rho = pr / 1.0e5;
                const double v = w / (rho * area);
                const double re = rho * v * 0.1 / 1.2e-5;
                const double f = churchillFriction(re, 1e-5 / 0.1);
                pr -= rho * 9.80665 * 40.0 / 24.0
                    + f * rho * v * v / (2.0 * 0.1) * ds;
            }
            return pr;
        };
        CHECK(std::abs(march(res.nodePressure[0]) - 1.5e6) < 1.0,
              "network pressure closes the PVT march");
        CHECK(res.pipeFlows[0].regime == -1,
              "PVT pipe reports single-phase");
    }

    if (failures)
    {
        std::cout << failures << " network test(s) FAILED\n";
        return 1;
    }
    std::cout << "all 1d network unit tests passed\n";
    return 0;
}
