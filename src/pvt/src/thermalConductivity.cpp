#include "cf/pvt/ThermalConductivity.h"

#include "cf/common/Error.h"

namespace cf::pvt
{

scalar thermalConductivityFromPr(scalar mu, scalar cp, scalar Pr)
{
    if (!(mu > 0) || !(cp > 0) || !(Pr > 0))
        cfFatal("pvt: thermalConductivityFromPr needs positive mu, cp and Pr");
    return mu * cp / Pr;
}

} // namespace cf::pvt
