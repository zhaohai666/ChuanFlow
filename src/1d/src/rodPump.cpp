#include "cf/1d/RodPump.h"

#include <algorithm>
#include <cmath>

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double circleArea(double diameter)
{
    return 0.25 * M_PI * diameter * diameter;
}

} // namespace

RodPumpCard rodPumpCardOfLoad(const RodPump& pump, double fluidDensity,
                              double pressureRise)
{
    RodPumpCard c;
    const double ap = circleArea(pump.plungerDiameter);
    const double ar = circleArea(pump.rodDiameter);
    const double l = pump.rodLength;
    const double n = pump.strokesPerMinute / 60.0; // strokes per second

    c.plungerArea = ap;
    c.rodArea = ar;
    c.stroke = pump.strokeLength;

    // Statement 1: the string as the pumping unit feels it.
    c.rodWeightInAir = pump.rodDensity * kGravity * ar * l;
    c.buoyancy = fluidDensity * kGravity * ar * l;
    c.rodWeightInFluid = c.rodWeightInAir - c.buoyancy;

    // Statement 2: the load the plunger holds is the differential its two faces
    // actually see, and the ceiling it would hold with nothing but the column
    // between it and the surface. A differential that does not oppose the
    // upstroke is not a negative load: the valves do not let the column push
    // the plunger down, so the load clamps at zero and the card is the flat
    // rod-weight line. The machine still sweeps its barrel -- that is what a
    // positive-displacement pump does -- but it lifts nothing, and the work the
    // rods hand to the fluid goes with the load, not with the sweep.
    c.fluidLoadMax = fluidDensity * kGravity * l * ap;
    const double lift = std::max(pressureRise, 0.0);
    c.fluidLoad = lift * ap;
    c.upthrust = c.fluidLoadMax - c.fluidLoad;
    if (!(lift > 0.0)) c.state = RodPumpState::noLoad;

    // Statements 3 and 4: the stretch the load costs, and the stroke left over.
    c.stretch = c.fluidLoad * l / (ar * pump.youngModulus);
    c.effectiveStroke = std::max(pump.strokeLength - c.stretch, 0.0);
    c.efficiency = c.effectiveStroke / pump.strokeLength;
    if (c.effectiveStroke <= 0.0) c.state = RodPumpState::stalled;

    c.maxLoad = c.rodWeightInFluid + c.fluidLoad;
    c.minLoad = c.rodWeightInFluid;
    c.downholeMaxLoad = c.fluidLoad;
    c.downholeMinLoad = 0.0;
    c.counterbalance = 0.5 * (c.maxLoad + c.minLoad);

    c.displacement = ap * c.effectiveStroke * n;
    c.massRate = fluidDensity * c.displacement;
    c.workPerCycle = c.fluidLoad * c.effectiveStroke;
    c.hydraulicPower = c.workPerCycle * n;
    return c;
}

RodPumpCard rodPumpCardOf(const RodPump& pump, double fluidDensity,
                          double pressureAbove, double pressureBelow)
{
    return rodPumpCardOfLoad(pump, fluidDensity, pressureAbove - pressureBelow);
}

} // namespace cf
