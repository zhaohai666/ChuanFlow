// Unit tests for the unified PVT layer: Peng-Robinson invariants that follow
// from the equation itself (critical triple root, ideal-gas limit, fugacity
// equality at saturation), the TP flash (material balance, phase equilibrium,
// single-phase collapse) and the tabulated-property interpolation.
#include <cmath>
#include <iostream>
#include <vector>

#include "cf/pvt/Flash.h"
#include "cf/pvt/PengRobinson.h"
#include "cf/pvt/PvtCase.h"
#include "cf/pvt/PvtTable.h"

using namespace cf;
using namespace cf::pvt;

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

static Mixture feedOf(std::vector<std::pair<Component, double>> comps)
{
    Mixture m;
    for (const auto& c : comps) m.add(c.first, c.second);
    m.normalise();
    return m;
}

int main()
{
    // ---- the critical point is built into the EOS: PR's own Zc ----
    {
        Mixture pure;
        pure.add(methane(), 1.0);
        const Component c = methane();
        const EosState st = pengRobinsonState(pure, c.Tc, c.Pc);
        // v_c / b from the EOS constants; a numerical triple root is delicate,
        // so the tolerance is the round-off amplification of the cubic.
        CHECK(relNear(st.Z, 0.30740, 2e-2), "critical root Z = PR Zc");
        CHECK(relNear(st.molarVolumeOverB, 0.30740 / 0.0777960739, 2e-2),
              "critical molar volume ratio");
    }

    // ---- low pressure reproduces the ideal gas law ----
    {
        Mixture pure;
        pure.add(methane(), 1.0);
        const double T = 300.0, p = 100.0; // 100 Pa: essentially ideal
        const EosState st = pengRobinsonState(pure, T, p);
        const double rhoIdeal = p * pure.molarMass() / (kGasConstant * T);
        CHECK(relNear(st.Z, 1.0, 1e-4), "Z -> 1 in the ideal-gas limit");
        CHECK(relNear(st.massDensity, rhoIdeal, 1e-4), "density -> ideal gas");
    }

    // ---- pure-component saturation from fugacity equality ----
    {
        const Component c = ethane();
        const double psat = pureVaporPressure(c, 250.0);
        CHECK(psat > 0.0 && psat < c.Pc, "vapour pressure inside the range");
        Mixture pure;
        pure.add(c, 1.0);
        const EosState liq = pengRobinsonState(pure, 250.0, psat, RootPreference::Liquid);
        const EosState vap = pengRobinsonState(pure, 250.0, psat, RootPreference::Vapor);
        CHECK(std::abs(liq.lnPhi[0] - vap.lnPhi[0]) < 1e-9,
              "equal liquid/vapour fugacity at saturation");
        // monotone in temperature, and the two saturation codes agree
        const double psatHot = pureVaporPressure(c, 280.0);
        CHECK(psatHot > psat, "vapour pressure increases with temperature");
        Mixture feed;
        feed.add(c, 1.0);
        CHECK(relNear(bubblePressure(feed, 250.0), psat, 1e-5),
              "bubble pressure of a pure fluid equals its vapour pressure");
        CHECK(pureVaporPressure(c, c.Tc + 1.0) == 0.0,
              "no saturation pressure above the critical temperature");
    }

    // ---- mixture flash: equilibrium and balance ----
    {
        const Mixture feed = feedOf({{methane(), 0.34}, {ethane(), 0.18},
                                     {propane(), 0.26}, {nPentane(), 0.22}});
        const double T = 320.0, p = 4.0e6;
        const FlashResult r = tpFlash(feed, T, p);
        CHECK(r.converged, "flash converged");
        CHECK(!r.singlePhase, "flash found two phases");
        CHECK(r.vaporFraction > 0.0 && r.vaporFraction < 1.0, "beta inside (0, 1)");

        const double b = r.vaporFraction;
        for (std::size_t i = 0; i < feed.size(); ++i)
        {
            const double recon = b * r.vapourComposition[i]
                               + (1.0 - b) * r.liquidComposition[i];
            CHECK(relNear(recon, feed.fractions()[i], 1e-10),
                  "component material balance");
        }
        double sx = 0, sy = 0;
        for (std::size_t i = 0; i < feed.size(); ++i)
        {
            sx += r.liquidComposition[i];
            sy += r.vapourComposition[i];
        }
        CHECK(relNear(sx, 1.0, 1e-12) && relNear(sy, 1.0, 1e-12), "phase sums to one");

        // isothermal-isobaric equilibrium: f_i^L = f_i^V, i.e.
        // ln(x_i phi_i^L) = ln(y_i phi_i^V)
        double worst = 0.0;
        for (std::size_t i = 0; i < feed.size(); ++i)
        {
            const double fl = std::log(r.liquidComposition[i]) + r.liquid.lnPhi[i];
            const double fv = std::log(r.vapourComposition[i]) + r.vapour.lnPhi[i];
            worst = std::max(worst, std::abs(fl - fv));
        }
        CHECK(worst < 1e-8, "equal component fugacities between phases");

        CHECK(r.liquidMassDensity > r.vapourMassDensity, "liquid denser than vapour");
        const std::size_t iC1 = feed.index("methane");
        const std::size_t iC5 = feed.index("nPentane");
        const double kC1 = r.vapourComposition[iC1] / r.liquidComposition[iC1];
        const double kC5 = r.vapourComposition[iC5] / r.liquidComposition[iC5];
        CHECK(kC1 > kC5, "lighter component partitions to the vapour");
        CHECK(r.vapourComposition[iC1] > feed.fractions()[iC1],
              "vapour enriched in the light component");

        // pressure ordering: more liquid as the pressure rises
        const FlashResult rLow = tpFlash(feed, T, 1.0e6);
        const FlashResult rHigh = tpFlash(feed, T, 1.2e7);
        CHECK(rLow.vaporFraction > r.vaporFraction, "low pressure gives more vapour");
        CHECK(rHigh.singlePhase && rHigh.vaporFraction == 0.0,
              "high pressure collapses to a single liquid phase");
        CHECK(r.singlePhase == false, "two-phase flag stays set");
        const FlashResult rVap = tpFlash(feed, T, 1.0e3);
        CHECK(rVap.singlePhase && rVap.vaporFraction == 1.0,
              "vacuum collapses to a single vapour phase");
    }

    // ---- mixture bubble point brackets the envelope the sweep shows ----
    {
        const Mixture feed = feedOf({{methane(), 0.34}, {ethane(), 0.18},
                                     {propane(), 0.26}, {nPentane(), 0.22}});
        const double pb = bubblePressure(feed, 320.0);
        CHECK(pb > 8.0e6 && pb < 9.0e6, "mixture bubble pressure sits in the envelope");
        CHECK(tpFlash(feed, 320.0, 0.99 * pb).singlePhase == false,
              "just below the bubble point the feed flashes two-phase");
        CHECK(tpFlash(feed, 320.0, 1.01 * pb).singlePhase,
              "just above the bubble point the feed is single-phase");
        const FlashResult at = tpFlash(feed, 320.0, pb);
        CHECK(at.vaporFraction < 1e-3, "the bubble point holds an infinitesimal vapour");
        // raising the temperature pushes the bubble pressure up
        CHECK(bubblePressure(feed, 340.0) > pb, "bubble pressure grows with temperature");
    }

    // ---- bilinear tables ----
    {
        // a linear field must be reproduced exactly anywhere in the domain
        const std::vector<double> xs{0.0, 1.0, 2.0}, ys{10.0, 20.0, 30.0};
        std::vector<double> values;
        for (const double x : xs)
            for (const double y : ys) values.push_back(2.0 * x + 0.1 * y);
        const BilinearTable table(xs, ys, values);
        CHECK(relNear(table.sample(0.5, 15.0), 2.5, 1e-12), "linear interior exact");
        CHECK(relNear(table.sample(1.75, 27.5), 6.25, 1e-12), "linear top-right cell exact");
        CHECK(relNear(table.sample(-3.0, 5.0), table.sample(0.0, 10.0), 1e-15),
              "sample clamps below the grid");
        CHECK(relNear(table.sample(9.0, 99.0), table.sample(2.0, 30.0), 1e-15),
              "sample clamps above the grid");
    }

    // ---- tables generated by the EOS agree with a direct evaluation ----
    {
        const Mixture feed = feedOf({{methane(), 0.34}, {ethane(), 0.18},
                                     {propane(), 0.26}, {nPentane(), 0.22}});
        const std::vector<double> pGrid{2.0e6, 4.0e6, 6.0e6};
        const std::vector<double> tGrid{310.0, 320.0, 330.0};
        const BilinearTable rhoL = liquidDensityTable(feed, pGrid, tGrid);
        const BilinearTable rhoV = vapourDensityTable(feed, pGrid, tGrid);
        const BilinearTable beta = vapourFractionTable(feed, pGrid, tGrid);

        CHECK(relNear(rhoL.sample(4.0e6, 320.0),
                      tpFlash(feed, 320.0, 4.0e6).liquidMassDensity, 1e-14),
              "table node equals the direct flash");
        CHECK(relNear(beta.sample(4.0e6, 320.0),
                      tpFlash(feed, 320.0, 4.0e6).vaporFraction, 1e-14),
              "vapour fraction table node exact");

        // Bilinear interpolation of a smooth property: half a cell in
        // pressure at the fixed grid temperature.
        const double err = std::abs(rhoL.sample(3.0e6, 320.0)
                                    - tpFlash(feed, 320.0, 3.0e6).liquidMassDensity)
                         / tpFlash(feed, 320.0, 3.0e6).liquidMassDensity;
        CHECK(err < 2e-2, "half-cell interpolation error below 2 percent");

        CHECK(rhoV.sample(2.0e6, 320.0) < rhoV.sample(4.0e6, 320.0)
              && rhoV.sample(4.0e6, 320.0) < rhoV.sample(6.0e6, 320.0),
              "vapour density grows with pressure");
        CHECK(beta.sample(2.0e6, 320.0) > beta.sample(4.0e6, 320.0)
              && beta.sample(4.0e6, 320.0) > beta.sample(6.0e6, 320.0),
              "vapour fraction falls with pressure");
        // Along the two-phase envelope the retained liquid gets *lighter* as
        // pressure rises: the light ends that flashed off at low pressure
        // dissolve back into it.
        CHECK(rhoL.sample(2.0e6, 320.0) > rhoL.sample(4.0e6, 320.0)
              && rhoL.sample(4.0e6, 320.0) > rhoL.sample(6.0e6, 320.0),
              "equilibrium liquid density falls along the envelope");
        CHECK(rhoL.sample(6.0e6, 310.0) > rhoL.sample(6.0e6, 330.0),
              "liquid density falls with temperature");
    }

    // ---- a compressed single-phase fluid is monotone in both coordinates ----
    {
        const Mixture feed = feedOf({{methane(), 0.34}, {ethane(), 0.18},
                                     {propane(), 0.26}, {nPentane(), 0.22}});
        const std::vector<double> pGrid{1.0e7, 1.2e7, 1.4e7};
        const std::vector<double> tGrid{310.0, 320.0, 330.0};
        const BilinearTable rho = liquidDensityTable(feed, pGrid, tGrid);
        for (const double p : pGrid)
            CHECK(tpFlash(feed, 320.0, p).singlePhase, "grid is single-phase");
        CHECK(rho.sample(1.0e7, 320.0) < rho.sample(1.2e7, 320.0)
              && rho.sample(1.2e7, 320.0) < rho.sample(1.4e7, 320.0),
              "single-phase density grows with pressure");
        CHECK(rho.sample(1.2e7, 330.0) < rho.sample(1.2e7, 320.0)
              && rho.sample(1.2e7, 320.0) < rho.sample(1.2e7, 310.0),
              "single-phase density falls with temperature");
    }

    // ---- case dictionary ----
    {
        const Dictionary caseDict = Dictionary::parseFile(CF_CASE_DIR "/pvtFlash/case.dict");
        const FlashCase c = loadFlashCase(caseDict);
        CHECK(c.feed.size() == 4, "case composition parsed");
        CHECK(relNear(c.feed.fractions()[0] + c.feed.fractions()[1]
                    + c.feed.fractions()[2] + c.feed.fractions()[3], 1.0, 1e-12),
              "composition normalised");
        CHECK(c.pressures.size() == 3, "case pressure list parsed");
        CHECK(c.controls.maxIter == 500, "case solver controls parsed");
        CHECK(c.computeBubblePoint, "case bubble-point flag parsed");

        const FlashResult r = tpFlash(c.feed, c.temperature, c.pressures[1], c.controls);
        CHECK(r.converged, "case flash converges");

        bool threw = false;
        try { requireComponent("unobtainium"); }
        catch (const Error&) { threw = true; }
        CHECK(threw, "unknown component is a fatal error");
    }

    if (failures)
    {
        std::cout << failures << " pvt test(s) FAILED\n";
        return 1;
    }
    std::cout << "all pvt unit tests passed\n";
    return 0;
}
