#include "cf/1d/PlungerGasLift.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

namespace
{

constexpr double kGravity = 9.80665;

} // namespace

bool plungerGasLiftWindowOpenAt(const PlungerGasLift& pgl, double t)
{
    if (!(pgl.cyclePeriod > 0.0) || !(pgl.openWindow > 0.0))
        cfFatal("plungerGasLift: the cycle period and open window are positive");
    if (!(t >= 0.0))
        cfFatal("plungerGasLift: the schedule reads non-negative time");
    const double phase = std::fmod(t, pgl.cyclePeriod);
    return phase < std::min(pgl.openWindow, pgl.cyclePeriod);
}

PlungerGasLiftSchedule plungerGasLiftSolve(const PlungerGasLift& pgl)
{
    if (!(pgl.gasVolumePerCycle > 0.0))
        cfFatal("plungerGasLift: the released gas volume per cycle is positive "
                "(a batch of zero gas is no batch)");
    if (!(pgl.referencePressure > 0.0) || !(pgl.gasDensityAtReference > 0.0))
        cfFatal("plungerGasLift: the volume's reference pressure and gas "
                "density are positive: the batch is quoted at a state");
    if (!(pgl.cyclePeriod > 0.0))
        cfFatal("plungerGasLift: the cycle period is positive");
    if (!(pgl.openWindow > 0.0) || pgl.openWindow > pgl.cyclePeriod)
        cfFatal("plungerGasLift: the open window must lie in (0, period]: the "
                "plunger trips the valve inside its own cycle, not across two");
    if (!(pgl.releaseDepth >= 0.0))
        cfFatal("plungerGasLift: the release depth is a non-negative vertical "
                "distance below the surface");
    if (!(pgl.casingSurfacePressure > 0.0) || !(pgl.tubingSurfacePressure > 0.0))
        cfFatal("plungerGasLift: the two surface pressures are absolute");
    if (!(pgl.liquidDensity > 0.0))
        cfFatal("plungerGasLift: the annulus liquid density is positive");

    PlungerGasLiftSchedule s;
    const double hydro = pgl.liquidDensity * kGravity * pgl.releaseDepth;
    s.injectionPressureAtDepth = pgl.casingSurfacePressure + hydro;
    s.tubingPressureAtDepth = pgl.tubingSurfacePressure + hydro;

    const GasLiftValveThresholds th =
        gasLiftValveThresholds(pgl.valve, s.tubingPressureAtDepth);
    s.openingPressure = th.openingPressure;
    s.closingPressure = th.closingPressure;

    // The batch's demand on the window is a property of the schedule, stated
    // whether or not the valve can meet it.
    s.requiredMassRate = pgl.gasDensityAtReference
                       * pgl.gasVolumePerCycle / pgl.openWindow;

    // The 5.89 convention: open past Pvo, shut at or below Pvc, and the band
    // between them is priced shut and reported, because the per-cycle
    // statement has no memory of where the plunger came from.
    const bool open = s.injectionPressureAtDepth >= th.openingPressure
        && s.injectionPressureAtDepth > th.closingPressure;
    if (!open)
    {
        s.state = s.injectionPressureAtDepth <= th.closingPressure
                      ? PlungerGasLiftState::closed
                      : PlungerGasLiftState::inSpread;
        return s;
    }

    // Isothermal ideal gas at the valve: density rides with pressure, so the
    // batch's reference density prices the nozzle at the injection state.
    const double rhoC = pgl.gasDensityAtReference
                      * s.injectionPressureAtDepth / pgl.referencePressure;
    const GasLiftValveFlow flow = gasLiftValveFlow(
        pgl.valve, s.injectionPressureAtDepth, rhoC, s.tubingPressureAtDepth);
    s.choked = flow.choked;
    s.nozzleMassRateCapacity = flow.massRate;

    const double required = s.requiredMassRate;
    const double passed = std::min(flow.massRate, required);
    if (flow.massRate < required) s.state = PlungerGasLiftState::valveLimited;
    else s.state = PlungerGasLiftState::open;

    s.massRateInWindow = passed;
    s.windowRateReference = passed / pgl.gasDensityAtReference;
    s.deliveredVolumeReference = s.windowRateReference * pgl.openWindow;
    s.averageRateReference = s.deliveredVolumeReference / pgl.cyclePeriod;
    // Boyle at the valve: p_ref V_ref = Pc(z) V(z).
    const double boyle = pgl.referencePressure / s.injectionPressureAtDepth;
    s.windowRateDownhole = s.windowRateReference * boyle;
    s.deliveredVolumeDownhole = s.deliveredVolumeReference * boyle;
    return s;
}

} // namespace cf
