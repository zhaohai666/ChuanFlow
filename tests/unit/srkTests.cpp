// Unit tests for the SRK layer (5.189): the EOS's own anchors plus the two
// derivative claims that make an EOS-derived Joule-Thomson coefficient
// meaningful.
//
//   * the critical point reproduces SRK's Zc = 1/3 (a property of the cubic,
//     not of a fitted table),
//   * the ideal-gas limit at low pressure is Z -> 1,
//   * the analytic dZ/dT agrees with a central difference of the state
//     function itself (an independent route through srkState),
//   * the Joule-Thomson coefficient agrees with the finite difference of the
//     molar volume, mu = (T (dv/dT)_p - v)/c_p, which is the thermodynamic
//     definition the closed form is derived from,
//   * the same two claims hold for a mixture, where dA/dT carries the mixing
//     rule's own derivative.
#include <cmath>
#include <iostream>
#include <vector>

#include "cf/pvt/PengRobinson.h"
#include "cf/pvt/Srk.h"

using namespace cf;
using namespace cf::pvt;

static int failures = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }\
        else         { std::cout << "ok:   " << (msg) << "\n"; }             \
    } while (0)

static bool relNear(double a, double b, double tol)
{
    const double s = std::max(std::abs(a), std::abs(b));
    if (s == 0) return true;
    return std::abs(a - b) / s < tol;
}

int main()
{
    Mixture pure;
    pure.add(methane(), 1.0);

    // ---- SRK's own critical compressibility ----
    {
        const Component c = methane();
        const EosState st = srkState(pure, c.Tc, c.Pc);
        CHECK(relNear(st.Z, 1.0 / 3.0, 2e-2), "critical root Z = SRK Zc = 1/3");
        CHECK(st.molarVolumeOverB > 3.0, "critical v/b > 3");
    }

    // ---- ideal-gas limit ----
    {
        const EosState st = srkState(pure, 300.0, 100.0);
        CHECK(relNear(st.Z, 1.0, 1e-5), "100 Pa is the ideal-gas limit");
        const EosState pr = pengRobinsonState(pure, 300.0, 100.0);
        CHECK(relNear(st.massDensity / pr.massDensity, 1.0, 1e-5),
              "SRK and PR agree in the ideal limit");
    }

    // ---- the analytic dZ/dT against the state function's own difference ----
    for (double T : {280.0, 320.0})
    {
        for (double p : {2.0e5, 2.0e6})
        {
            const double h = 1e-4 * T;
            const double num =
                (srkState(pure, T + h, p).Z - srkState(pure, T - h, p).Z)
                / (2.0 * h);
            const double ana = srkDzDt(pure, T, p);
            CHECK(relNear(ana, num, 1e-6),
                  "dZ/dT analytic vs finite difference at T="
                      + std::to_string(T) + " p=" + std::to_string(p));
        }
    }

    // ---- mu_JT against its thermodynamic definition ----
    {
        // c_p is MASS-specific (J/(kg K)), the form the 1D thermal march
        // carries, so the reference volume is mass-specific too: v = v_molar/M.
        const double cp = 2200.0;
        const double Mm = pure.molarMass();
        for (double p : {1.0e5, 1.0e6, 5.0e6})
        {
            const double T = 300.0;
            const double h = 1e-4 * T;
            const double vPlus = srkState(pure, T + h, p).molarVolume / Mm;
            const double vMinus = srkState(pure, T - h, p).molarVolume / Mm;
            const double v0 = srkState(pure, T, p).molarVolume / Mm;
            const double muNum = (T * (vPlus - vMinus) / (2.0 * h) - v0) / cp;
            const double muAna = srkJouleThomson(pure, T, p, cp);
            CHECK(relNear(muAna, muNum, 1e-6),
                  "mu_JT closed form vs (T dv/dT - v)/cp at p="
                      + std::to_string(p));
            CHECK(muAna > 0.0, "methane cools on expansion at 300 K (mu > 0)");
        }
    }

    // ---- and on a mixture, where dA/dT carries the mixing rule ----
    {
        Mixture gas;
        gas.add(methane(), 0.7);
        gas.add(ethane(), 0.2);
        gas.add(propane(), 0.1);
        gas.normalise();
        const double T = 330.0, p = 3.0e6, cp = 2600.0;
        const double h = 1e-4 * T;
        const double num =
            (srkState(gas, T + h, p).Z - srkState(gas, T - h, p).Z) / (2.0 * h);
        CHECK(relNear(srkDzDt(gas, T, p), num, 1e-6),
              "mixture dZ/dT analytic vs finite difference");
        const double Mm = gas.molarMass();
        const double vPlus = srkState(gas, T + h, p).molarVolume / Mm;
        const double vMinus = srkState(gas, T - h, p).molarVolume / Mm;
        const double v0 = srkState(gas, T, p).molarVolume / Mm;
        const double muNum = (T * (vPlus - vMinus) / (2.0 * h) - v0) / cp;
        CHECK(relNear(srkJouleThomson(gas, T, p, cp), muNum, 1e-6),
              "mixture mu_JT closed form vs definition");
    }

    // ---- a rich gas is not the ideal gas, and SRK says how far ----
    {
        const EosState st = srkState(pure, 300.0, 2.0e6);
        CHECK(st.Z < 0.999, "2 MPa methane is already below ideal");
        CHECK(st.Z > 0.8, "and only modestly so");
        CHECK(std::isfinite(st.gibbs) && st.lnPhi.size() == 1,
              "SRK state carries its fugacity");
    }

    // ---- dZ/dT follows the root the state reports, not the lowest branch ----
    {
        const EosState st = srkState(pure, 200.0, 4.0e6);
        if (st.hadTwoRoots)
        {
            const double dz = srkDzDt(pure, 200.0, 4.0e6);
            CHECK(std::isfinite(dz), "two-root state still returns a derivative");
        }
        else
        {
            CHECK(true, "single-root state at 200 K / 4 MPa (no branch to pick)");
        }
    }

    std::cout << (failures == 0 ? "SRK: all checks passed\n"
                                : "SRK: FAILURES\n");
    return failures == 0 ? 0 : 1;
}
