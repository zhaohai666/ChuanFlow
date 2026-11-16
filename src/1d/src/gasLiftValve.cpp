#include "cf/1d/GasLiftValve.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

namespace
{

// The test rack has the tubing side open to atmosphere, so the setting a valve
// is shipped with is its opening pressure at one standard atmosphere absolute.
constexpr double kAtmosphericPressure = 101325.0;

} // namespace

GasLiftValveThresholds gasLiftValveThresholds(const GasLiftValve& valve,
                                              double tubingPressure)
{
    if (!(valve.areaRatio > 0.0) || !(valve.areaRatio < 1.0))
        cfFatal("gasLiftValve: the port-to-bellows area ratio must be in (0,1)");
    if (!(valve.domePressure > 0.0) || !(tubingPressure > 0.0))
        cfFatal("gasLiftValve: the dome and tubing pressures are absolute");

    const double r = valve.areaRatio;
    GasLiftValveThresholds t;
    t.openingPressure = valve.domePressure / (1.0 - r) + valve.springPressure
                      - (r / (1.0 - r)) * tubingPressure;
    t.closingPressure = valve.domePressure
                      + valve.springPressure * (1.0 - r);
    t.spread = t.openingPressure - t.closingPressure;
    t.productionEffect = (r / (1.0 - r)) * tubingPressure;
    t.testRackOpening = valve.domePressure / (1.0 - r) + valve.springPressure
                      - (r / (1.0 - r)) * kAtmosphericPressure;
    return t;
}

GasLiftValveFlow gasLiftValveFlow(const GasLiftValve& valve, double pUp,
                                  double rhoUp, double pDown)
{
    GasLiftValveFlow f;
    if (!(valve.portArea > 0.0) || !(valve.dischargeCoeff > 0.0))
        cfFatal("gasLiftValve: the port area and discharge coefficient are positive");
    if (!(valve.kappa > 1.0))
        cfFatal("gasLiftValve: the isentropic exponent must exceed 1");
    if (!(pUp > 0.0) || !(rhoUp > 0.0))
        cfFatal("gasLiftValve: the upstream state needs a positive pressure and density");

    // A check valve: the port cannot pull gas back out of the tubing.
    if (!(pDown > 0.0) || pDown >= pUp) return f;

    const double k = valve.kappa;
    // The ratio the flux brackets at, so the throat is sonic at or below it.
    // It is the bracket's own maximum, not a table lookup: d/dr [r^(2/k)
    // - r^((k+1)/k)] = 0 gives 2 r^((1-k)/k) = k + 1, i.e. the exponent k/(k-1)
    // - the textbook 0.528 at k = 1.4, not the 0.634 that exponent's neighbour
    // 1/(k-1), the throat density ratio's, would give. Choking a valve at the
    // wrong ratio mis-prices every passing arm by the bracket between the two,
    // 1.6% at the k = 1.3 of the lift gas on the rack.
    const double rc = std::pow(2.0 / (k + 1.0), k / (k - 1.0));

    const double r = pDown / pUp;
    f.choked = r <= rc;
    const double rs = f.choked ? rc : r;
    const double bracket = std::pow(rs, 2.0 / k) - std::pow(rs, (k + 1.0) / k);
    const double flux2 = 2.0 * rhoUp * pUp * (k / (k - 1.0)) * bracket;
    f.massRate = valve.dischargeCoeff * valve.portArea * std::sqrt(
        flux2 > 0.0 ? flux2 : 0.0);
    return f;
}

} // namespace cf
