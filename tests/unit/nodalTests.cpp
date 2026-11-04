// Unit tests for the nodal-analysis kernel: the closed-form IPR curves
// (linear / Vogel / Jones) and their inverses, the wellbore VLP inversion of
// the Beggs-Brill march (static column, friction-free limit), the operating
// point of the crossing with its stability side, and the dead-well verdict.
#include <cmath>
#include <iostream>

#include "cf/1d/Nodal.h"

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

static pvt::BlackOilFluid deadOil()
{
    pvt::BlackOilFluid f;
    f.gasGravity = 0.65;
    f.oilDensityStd = 850.0;
    f.gorStd = 0.0;
    f.temperature = 340.0;
    return f;
}

int main()
{
    // ---- the IPR curves are their own published closed forms ----
    {
        IpR ipr;
        ipr.form = IpR::Form::linear;
        ipr.reservoirPressure = 2.0e7;
        ipr.productivityIndex = 5.0e-12;
        CHECK(relNear(ipr.rateMax(), 5.0e-12 * 2.0e7, 1e-12),
              "linear capacity is J p_res");
        CHECK(relNear(ipr.pWf(0.0), 2.0e7, 1e-12), "linear pWf(0) = p_res");
        CHECK(ipr.pWf(ipr.rateMax()) < 1.0, "linear pWf(capacity) = 0");
        CHECK(relNear(ipr.rate(ipr.pWf(0.6 * ipr.rateMax())),
                      0.6 * ipr.rateMax(), 1e-12),
              "linear inverse round trip");
    }
    {
        IpR ipr;
        ipr.form = IpR::Form::vogel;
        ipr.reservoirPressure = 2.2e7;
        ipr.maxRate = 0.006;
        // at x = 0.5 the Vogel curve reads 1 - 0.1 - 0.2 = 0.7 of the AOF
        CHECK(relNear(ipr.rate(0.5 * 2.2e7), 0.7 * 0.006, 1e-12),
              "vogel rate(0.5 p_res) = 0.7 qmax");
        CHECK(relNear(ipr.pWf(0.7 * 0.006), 0.5 * 2.2e7, 1e-12),
              "vogel pWf(0.7 qmax) = 0.5 p_res");
        CHECK(relNear(ipr.pWf(0.0), 2.2e7, 1e-12), "vogel pWf(0) = p_res");
        CHECK(ipr.pWf(0.006) < 1.0, "vogel pWf(qmax) = 0");
        bool roundTrip = true;
        for (int k = 1; k < 10; ++k)
            roundTrip = roundTrip
                     && relNear(ipr.rate(ipr.pWf(0.1 * k * 0.006)),
                                0.1 * k * 0.006, 1e-12);
        CHECK(roundTrip, "vogel inverse round trip over the curve");
    }
    {
        IpR ipr;
        ipr.form = IpR::Form::jones;
        ipr.reservoirPressure = 2.0e7;
        ipr.deliverability = 4.0e-19;
        ipr.exponent = 1.0;
        // at x = 0.6: p_res^2 - p_wf^2 = 0.64 p_res^2
        CHECK(relNear(ipr.rate(0.6 * 2.0e7),
                      4.0e-19 * 0.64 * 4.0e14, 1e-12),
              "jones rate(0.6 p_res) = 0.64 C p_res^2 at n = 1");
        CHECK(relNear(ipr.rate(ipr.pWf(0.3 * ipr.rateMax())),
                      0.3 * ipr.rateMax(), 1e-12),
              "jones inverse round trip");
        ipr.exponent = 2.0;
        CHECK(relNear(ipr.rate(ipr.pWf(0.4 * ipr.rateMax())),
                      0.4 * ipr.rateMax(), 1e-12),
              "jones inverse round trip at n = 2");
    }

    // ---- the VLP inversion: static column and its horizontal limit ----
    {
        const pvt::BlackOilFluid f = deadOil();
        WellLift well;
        well.length = 1000.0;
        well.diameter = 0.1;
        well.roughness = 0.0;
        well.elevation = 1000.0;
        const double pStatic = vlpBottomHolePressure(f, well, 0.0, 2.5e6);
        CHECK(pStatic > 2.5e6 + 0.8 * 850.0 * 9.80665 * 1000.0,
              "the static column carries at least the lightest gas-free head");
        CHECK(pStatic < 2.5e6 + 860.0 * 9.80665 * 1000.0,
              "and it stays below the stock-tank density head");
        const BeggsBrillMarchOut m =
            marchBeggsBrillPipe(f, 0.0, pStatic, 1000.0, 0.1, 0.0, 1.0, 24);
        CHECK(m.ok && std::abs(m.pOutlet - 2.5e6) < 1e-4,
              "marching back down the static column lands on the wellhead");

        WellLift flat = well;
        flat.elevation = 0.0;
        CHECK(std::abs(vlpBottomHolePressure(f, flat, 0.0, 2.5e6) - 2.5e6)
                  < 1e-6,
              "a horizontal rate-0 wellbore holds the wellhead pressure");
    }

    // ---- friction-free limit: a fat pipe makes the VLP flat, so the
    //      operating rate collapses onto J (p_res - static) ----
    {
        const pvt::BlackOilFluid f = deadOil();
        IpR ipr;
        ipr.form = IpR::Form::linear;
        ipr.reservoirPressure = 1.5e7;
        ipr.productivityIndex = 5.0e-12;
        WellLift well;
        well.length = 1000.0;
        well.diameter = 0.5;
        well.roughness = 0.0;
        well.elevation = 1000.0;
        const double pStatic =
            vlpBottomHolePressure(f, well, 0.0, 2.5e6);
        const double qExpected = ipr.productivityIndex * (1.5e7 - pStatic);
        const NodalControls ctl;
        const NodalResult r = solveNodalAnalysis(ipr, f, well, 2.5e6, ctl);
        CHECK(r.flows && r.converged, "the fat-pipe well flows and converges");
        CHECK(relNear(r.rate, qExpected, 1e-6),
              "with friction switched off the rate is exactly J(p_res - p_static)");
        CHECK(relNear(r.bottomHolePressure, pStatic, 1e-6),
              "and the bottom hole sits on the static column");
    }

    // ---- the vogel gas well: closed-form IPR identity, VLP fixed point,
    //      stability side ----
    IpR vogel;
    vogel.form = IpR::Form::vogel;
    vogel.reservoirPressure = 2.2e7;
    vogel.maxRate = 0.006;
    pvt::BlackOilFluid gasOil;
    gasOil.gasGravity = 0.65;
    gasOil.oilDensityStd = 850.0;
    gasOil.gorStd = 90.0;
    gasOil.temperature = 340.0;
    gasOil.surfaceTension = 0.03;
    WellLift well;
    well.length = 2400.0;
    well.diameter = 0.062;
    well.roughness = 4.6e-5;
    well.elevation = 2200.0;
    const double pWh = 2.5e6;
    {
        const NodalControls ctl;
        const NodalResult r = solveNodalAnalysis(vogel, gasOil, well, pWh, ctl);
        CHECK(r.flows && r.converged, "the vogel gas well flows and converges");
        CHECK(r.rate > 0.0 && r.rate < vogel.rateMax(),
              "the operating rate sits inside the IPR capacity");
        CHECK(r.bottomHolePressure > pWh && r.bottomHolePressure < 2.2e7,
              "the bottom hole is between the wellhead and the reservoir");
        const double x = r.bottomHolePressure / 2.2e7;
        const double vogelQ = 0.006 * (1.0 - 0.2 * x - 0.8 * x * x);
        CHECK(std::abs(r.rate - vogelQ) < 1e-4 * r.rate,
              "the operating rate satisfies the Vogel identity to the solve "
              "residual");
        const BeggsBrillMarchOut m = marchBeggsBrillPipe(
            gasOil, r.rate, r.bottomHolePressure, 2400.0, 0.062, 4.6e-5,
            2200.0 / 2400.0, 24);
        CHECK(m.ok && std::abs(m.pOutlet - pWh) < 1.0,
              "marching the operating rate from pwf lands on the wellhead");
        // stability: the crossing sits on the rising friction branch, where
        // pushing the rate up demands more bottom-hole pressure than the
        // IPR is willing to give
        const double fLo = vlpBottomHolePressure(gasOil, well, 0.9 * r.rate, pWh)
                         - vogel.pWf(0.9 * r.rate);
        const double fHi = vlpBottomHolePressure(gasOil, well, 1.1 * r.rate, pWh)
                         - vogel.pWf(1.1 * r.rate);
        CHECK(fLo < 0.0 && fHi > 0.0,
              "f changes sign across the operating point (it is a crossing)");
        CHECK(r.holdup > 0.0 && r.holdup < 1.0,
              "the gas-bearing column reports a two-phase holdup");

        // the same well on dead oil: no gas to lighten the column, so the
        // lift is harder and the operating rate falls
        const NodalResult d = solveNodalAnalysis(vogel, deadOil(), well, pWh, ctl);
        CHECK(d.flows && d.converged, "the dead-oil twin also flows");
        CHECK(d.rate < r.rate,
              "gas lightens the column and raises the operating rate");
        CHECK(d.bottomHolePressure > r.bottomHolePressure,
              "the dead-oil twin needs more bottom-hole pressure");
    }

    // ---- dead well: the reservoir cannot lift the static column ----
    {
        IpR weak;
        weak.form = IpR::Form::vogel;
        weak.reservoirPressure = 5.0e6;
        weak.maxRate = 0.006;
        const NodalControls ctl;
        const NodalResult r =
            solveNodalAnalysis(weak, gasOil, well, pWh, ctl);
        CHECK(r.converged && !r.flows, "the dead well resolves, it just does not flow");
        CHECK(r.rate == 0.0, "the dead well produces nothing");
        CHECK(std::abs(r.bottomHolePressure - r.staticBottomPressure) < 1e-9,
              "the dead well's bottom hole is the static column");
    }

    // ---- dictionary loader against the shipped case ----
    {
        const Dictionary caseDict = Dictionary::parseFile(
            CF_CASE_DIR "/nodalVogel/case.dict");
        const NodalCase c = loadNodalCaseFromDict(caseDict);
        CHECK(c.ipr.form == IpR::Form::vogel, "the case asks for the vogel IPR");
        CHECK(relNear(c.ipr.reservoirPressure, 2.2e7, 1e-12),
              "the case's reservoir pressure is 2.2e7");
        CHECK(relNear(c.ipr.maxRate, 0.006, 1e-12), "the case's AOF is 0.006");
        CHECK(relNear(c.well.length, 2400.0, 1e-12), "the well is 2400 m");
        CHECK(relNear(c.well.elevation, 2200.0, 1e-12),
              "the well rises 2200 m");
        CHECK(relNear(c.wellheadPressure, 2.5e6, 1e-12),
              "the wellhead is held at 2.5 MPa");
        CHECK(relNear(c.fluid.gorStd, 90.0, 1e-12), "the case's GOR is 90");
    }

    if (failures)
    {
        std::cout << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "all nodal tests passed\n";
    return 0;
}
