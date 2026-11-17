#include "cf/pvt/WaxSolubility.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf::pvt
{

scalar waxSaturationMoleFraction(scalar T, scalar Tf, scalar dHf)
{
    if (!(T > 0) || !(Tf > 0))
        cfFatal("pvt: waxSaturationMoleFraction needs positive T and Tf");
    // ln x_sat = -(dHf/R)(1/T - 1/Tf); evaluate in the log first so the
    // closed form is reproduced bit-for-bit before the exp (a checker using
    // math.exp of the same log agrees to the 12-digit report floor).
    const scalar lnXsat = -(dHf / kGasConstant) * (1.0 / T - 1.0 / Tf);
    return std::exp(lnXsat);
}

} // namespace cf::pvt
