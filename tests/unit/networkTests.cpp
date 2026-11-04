// Unit tests for the 1D steady pipe-network kernel: Churchill friction,
// single-pipe Hagen-Poiseuille behaviour, gravity head, and damped-Newton
// network solutions against closed-form series/parallel results.
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

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

    if (failures)
    {
        std::cout << failures << " network test(s) FAILED\n";
        return 1;
    }
    std::cout << "all 1d network unit tests passed\n";
    return 0;
}
