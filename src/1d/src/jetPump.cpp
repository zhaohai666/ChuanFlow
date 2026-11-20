#include "cf/1d/JetPump.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

namespace
{

// The matching function of the throat balance, g(r) = m^2 (p_d - p_s) /
// (rho v_n^2): the square term minus the friction line. Strictly decreasing
// on [0, m-1] -- its derivative is -2(m-1-r)/(m-1) - F(1+r), both terms of
// which are <= 0 there and never zero together -- so the operating point the
// pressures ask for is unique.
double throatMatching(double r, double m, double fLine)
{
    const double s = m - 1.0 - r;
    return s * s / (m - 1.0) - 0.5 * fLine * (1.0 + r) * (1.0 + r);
}

} // namespace

JetPumpSolution jetPumpSolve(const JetPump& pump, const JetPumpOperation& op)
{
    if (!(pump.nozzleArea > 0.0))
        cfFatal("jetPump: the nozzle area is positive (a bore cannot be zero or negative)");
    if (!(pump.throatArea > pump.nozzleArea))
        cfFatal("jetPump: the throat area must exceed the nozzle area: the "
                "annulus A_t - A_n is what the suction stream enters through, "
                "and the momentum balance divides by m - 1");
    if (!(pump.throatLength >= 0.0) || !(pump.frictionFactor >= 0.0))
        cfFatal("jetPump: the throat length and friction factor are non-negative");
    if (!(pump.fluidDensity > 0.0))
        cfFatal("jetPump: the fluid density is positive");
    if (!(op.motiveRate > 0.0))
        cfFatal("jetPump: the motive rate is positive: with no jet there is no "
                "momentum to hand to the suction stream");
    if (!(op.dischargePressure > 0.0))
        cfFatal("jetPump: the discharge pressure is absolute");

    JetPumpSolution s;
    const double m = pump.throatArea / pump.nozzleArea;
    const double vn = op.motiveRate / pump.nozzleArea;
    const double v2 = vn * vn;
    const double throatDiameter = std::sqrt(4.0 * pump.throatArea / M_PI);
    const double fLine = pump.frictionFactor * pump.throatLength / throatDiameter;

    s.areaRatio = m;
    s.nozzleVelocity = vn;
    s.motiveRate = op.motiveRate;
    s.dischargePressure = op.dischargePressure;

    double r = 0.0;
    double dp = 0.0; // p_d - p_s, Pa

    if (op.imposeSuctionRate)
    {
        if (!(op.suctionRate >= 0.0) || op.suctionRate > (m - 1.0) * op.motiveRate)
            cfFatal("jetPump: the imposed suction rate is outside [0, (m-1) q_n]: "
                    "m - 1 is the free-jet limit of the geometry, where the "
                    "throat makes no pressure rise at all");
        r = op.suctionRate / op.motiveRate;
        dp = pump.fluidDensity * v2 * throatMatching(r, m, fLine)
           / (m * m);
        if (!(dp >= 0.0))
            cfFatal("jetPump: the imposed suction rate leaves a negative throat "
                    "rise (discharge below chamber): the declared friction "
                    "cannot carry that much entrained flow");
        s.suctionPressure = op.dischargePressure - dp;
    }
    else
    {
        if (!(op.suctionPressure > 0.0))
            cfFatal("jetPump: the suction pressure is absolute");
        if (op.dischargePressure < op.suctionPressure)
            cfFatal("jetPump: the discharge pressure is below the suction "
                    "pressure: the chamber would feed the throat backwards, "
                    "which this element does not price");
        dp = op.dischargePressure - op.suctionPressure;
        s.suctionPressure = op.suctionPressure;

        const double target = dp * m * m / (pump.fluidDensity * v2);
        const double gMax = throatMatching(0.0, m, fLine);
        if (target >= gMax)
        {
            // Back-pressure stall: even with the chamber taking nothing in,
            // the jet overshoots the required rise; the suction port sees a
            // check valve and the pump passes its motive stream alone.
            r = 0.0;
            s.state = JetPumpState::noSuction;
        }
        else
        {
            // a r^2 + b r + c = 0 from the expanded matching function.
            const double a = 1.0 / (m - 1.0) - 0.5 * fLine;
            const double b = -(2.0 + fLine);
            const double c = (m - 1.0) - 0.5 * fLine - target;
            if (std::fabs(a) < 1e-12)
            {
                r = -c / b;
            }
            else
            {
                const double disc = b * b - 4.0 * a * c;
                if (!(disc >= 0.0))
                    cfFatal("jetPump: the throat balance has no real root");
                // The branch on [0, m-1]: the smaller root for a > 0 (the
                // larger one sits past the free-jet limit), the only positive
                // root for a < 0. Both are (-b - sqrt(disc)) / (2a).
                r = (-b - std::sqrt(disc)) / (2.0 * a);
            }
            if (r < 0.0) r = 0.0;
            if (r > m - 1.0) r = m - 1.0;
            if (r >= (m - 1.0) * (1.0 - 1e-12)) s.state = JetPumpState::freeJet;
        }
    }

    const double rhoV2 = pump.fluidDensity * v2;
    const double dpSolved = rhoV2 * throatMatching(r, m, fLine) / (m * m);
    const double supply = s.suctionPressure + 0.5 * rhoV2;

    s.flowRatio = r;
    s.suctionRate = r * op.motiveRate;
    s.dischargeRate = (1.0 + r) * op.motiveRate;
    s.throatVelocity = vn * s.dischargeRate / op.motiveRate / m;
    s.supplyPressure = supply;
    s.momentumResidual = dpSolved - dp;

    const double motiveRise = supply - op.dischargePressure;
    if (motiveRise > 0.0 && dp > 0.0)
    {
        s.headRatio = dp / motiveRise;
        s.efficiency = r * s.headRatio;
    }
    s.usefulPower = s.suctionRate * dp;
    s.motivePower = op.motiveRate * (motiveRise > 0.0 ? motiveRise : 0.0);
    return s;
}

} // namespace cf
