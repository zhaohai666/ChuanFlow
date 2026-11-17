#include "cf/pvt/SurfaceTension.h"

#include "cf/common/Error.h"

namespace cf::pvt
{

scalar weinaugKatzSurfaceTension(const std::vector<scalar>& parachors,
                                 scalar Vml, scalar Vmg,
                                 const std::vector<scalar>& xs,
                                 const std::vector<scalar>& ys)
{
    const std::size_t n = parachors.size();
    if (n == 0 || xs.size() != n || ys.size() != n)
        cfFatal("pvt: weinaugKatz parachor/composition arrays disagree in length");
    if (!(Vml > 0.0) || !(Vmg > 0.0))
        cfFatal("pvt: weinaugKatz needs positive phase molar volumes");

    const scalar rhoL = 1.0 / Vml;
    const scalar rhoV = 1.0 / Vmg;
    scalar tot = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        tot += parachors[i] * (xs[i] * rhoL - ys[i] * rhoV);
    tot *= tot;
    tot *= tot; // fourth power
    return tot;
}

scalar weinaugKatzSurfaceTension(const Mixture& mix, scalar T, scalar p,
                                 const std::vector<scalar>& parachors)
{
    if (parachors.size() != mix.size())
        cfFatal("pvt: weinaugKatz needs one parachor per component");
    const FlashResult r = tpFlash(mix, T, p);
    if (r.singlePhase) return 0.0;
    if (!(r.liquidMassDensity > 0.0) || !(r.vapourMassDensity > 0.0)) return 0.0;
    const scalar Vml = r.liquidMolarMass / r.liquidMassDensity;
    const scalar Vmg = r.vapourMolarMass / r.vapourMassDensity;
    return weinaugKatzSurfaceTension(parachors, Vml, Vmg,
                                     r.liquidComposition, r.vapourComposition);
}

} // namespace cf::pvt