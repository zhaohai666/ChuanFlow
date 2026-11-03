#include "cf/turbulence/WallFunction.h"

#include <cmath>

namespace cf
{

scalar spaldingYPlus(scalar uPlus)
{
    if (uPlus <= 0) return 0;
    const scalar k = turb::kappa;
    const scalar ku = k * uPlus;
    // e^(-k C) = 1/E with C the log-law intercept
    const scalar damp = std::exp(-k * turb::spaldingC);
    return uPlus + damp * (std::exp(ku) - 1.0 - ku - 0.5 * ku * ku - ku * ku * ku / 6.0);
}

namespace
{

/// Bisection of a strictly increasing f on [0, hi] with a doubling bracket.
/// Shared by the two inversions of the law of the wall, which differ only in f.
template<typename F>
scalar bisect(F f, scalar target)
{
    scalar lo = 0.0;
    scalar hi = 1.0;
    while (f(hi) < target)
    {
        lo = hi;
        hi *= 2.0;
        if (hi > 1e6) break; // f overflows to +inf long before this
    }
    for (int i = 0; i < 100 && hi - lo > 1e-13 * hi; ++i)
    {
        const scalar mid = 0.5 * (lo + hi);
        if (f(mid) > target) hi = mid;
        else lo = mid;
    }
    return 0.5 * (lo + hi);
}

} // namespace

scalar lawOfWall(scalar yPlus)
{
    if (yPlus <= 0) return 0;
    return bisect(spaldingYPlus, yPlus);
}

WallState wallState(scalar uParallel, scalar y, scalar nu)
{
    WallState s;
    if (!(uParallel > 0) || !(y > 0) || !(nu > 0)) return s;

    // u_par = u_tau u+ and y+ = u_tau y/nu give u+ y+(u+) = u_par y/nu
    const scalar a = uParallel * y / nu;
    const scalar uPlus = bisect([&](scalar up) { return up * spaldingYPlus(up); }, a);

    s.uTau = uParallel / uPlus;
    s.yPlus = a / uPlus;
    // nu + nut = u_tau^2 y / uParallel by construction of the wall function
    s.nut = s.uTau * s.uTau * y / uParallel - nu;
    if (s.nut < 0) s.nut = 0;
    return s;
}

} // namespace cf
