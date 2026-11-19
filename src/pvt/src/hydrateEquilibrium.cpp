#include "cf/pvt/HydrateEquilibrium.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf::pvt
{

scalar hydrateEquilibriumPressure(scalar T, scalar Tref, scalar pRef, scalar dH)
{
    if (!(T > 0) || !(Tref > 0) || !(pRef > 0) || !(dH > 0))
        cfFatal("pvt: hydrateEquilibriumPressure needs positive T, Tref, pRef and dH");
    // Evaluate in the log first: a checker using math.exp of the same log
    // reproduces the reported value to the 12-digit floor.
    const scalar lnRatio = -(dH / kGasConstant) * (1.0 / T - 1.0 / Tref);
    return pRef * std::exp(lnRatio);
}

scalar hydrateEquilibriumTemperature(scalar p, scalar Tref, scalar pRef, scalar dH)
{
    if (!(p > 0) || !(Tref > 0) || !(pRef > 0) || !(dH > 0))
        cfFatal("pvt: hydrateEquilibriumTemperature needs positive p, Tref, pRef and dH");
    const scalar invT = 1.0 / Tref - (kGasConstant / dH) * std::log(p / pRef);
    if (!(invT > 0))
        cfFatal("pvt: hydrateEquilibriumTemperature: no physical root for this pressure");
    return 1.0 / invT;
}

} // namespace cf::pvt
