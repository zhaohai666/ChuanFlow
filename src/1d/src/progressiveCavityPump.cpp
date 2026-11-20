#include "cf/1d/ProgressiveCavityPump.h"

#include <algorithm>

#include "cf/common/Error.h"

namespace cf
{

ProgressiveCavityPumpCard progressiveCavityPumpCardOfLoad(
    const ProgressiveCavityPump& pump, double fluidDensity,
    double pressureRise)
{
    if (!(pump.cavitiesPerRev > 0.0))
        cfFatal("progressiveCavityPump: the cavity count per revolution must "
                "be positive: at least one sealed cavity is what makes the "
                "machine a pump");
    if (!(pump.cavityVolume > 0.0))
        cfFatal("progressiveCavityPump: the cavity volume is positive (a "
                "zero-volume cavity displaces nothing however fast it turns)");
    if (!(pump.rpm > 0.0))
        cfFatal("progressiveCavityPump: the rpm is positive: a non-turning "
                "PCP declared with a flow is a contradiction -- the "
                "displacement is kinematics, q_th = cavities x V x rpm / 60");
    if (!(pump.slip >= 0.0))
        cfFatal("progressiveCavityPump: the slip line is non-negative: "
                "clearance leakage only ever takes volume back off the "
                "discharge side");
    if (!(fluidDensity > 0.0))
        cfFatal("progressiveCavityPump: the fluid density is positive");

    ProgressiveCavityPumpCard c;
    c.theoreticalRate = pump.cavitiesPerRev * pump.cavityVolume * pump.rpm
                      / 60.0;

    // The differential the machine works against: a rise that does not oppose
    // the transfer (NO-LOAD, the rod pump's rule) leaks nothing back and costs
    // no work; the cavities still sweep their volume.
    const double rise = std::max(pressureRise, 0.0);
    c.pressureRise = rise;
    if (!(rise > 0.0)) c.state = ProgressiveCavityPumpState::noLoad;

    c.slipRate = pump.slip * rise;
    c.netRate = c.theoreticalRate - c.slipRate;
    if (rise > 0.0 && !(c.netRate > 0.0))
        c.state = ProgressiveCavityPumpState::stalled;

    c.massRate = fluidDensity * c.netRate;
    c.volumetricEfficiency = c.netRate / c.theoreticalRate;
    c.hydraulicPower = c.netRate * rise;
    return c;
}

ProgressiveCavityPumpCard progressiveCavityPumpCardOf(
    const ProgressiveCavityPump& pump, double fluidDensity,
    double pressureDischarge, double pressureIntake)
{
    if (!(pressureDischarge > 0.0) || !(pressureIntake > 0.0))
        cfFatal("progressiveCavityPump: the port pressures are absolute");
    return progressiveCavityPumpCardOfLoad(pump, fluidDensity,
                                           pressureDischarge - pressureIntake);
}

} // namespace cf
