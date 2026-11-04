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

scalar jayatillekeP(scalar Pr, scalar Prt)
{
    const scalar r = Pr / Prt;
    return 9.24 * (std::pow(r, 0.75) - 1.0)
         * (1.0 + 0.28 * std::exp(-0.007 * r));
}

scalar thermalSublayerYPlus(scalar Pr, scalar Prt)
{
    const scalar P = jayatillekeP(Pr, Prt);
    // f(y) = Pr y - Prt/k ln(E y) - Prt P: falls to its minimum at
    // y = Prt/(k Pr) (where the log still wins) and rises from there, so any
    // crossing below the minimum is spurious and the physical sublayer
    // thickness is bracketed upward from the minimum
    auto f = [&](scalar y)
    {
        return Pr * y - Prt / turb::kappa * std::log(turb::logE * y) - Prt * P;
    };
    scalar lo = Prt / (turb::kappa * Pr);
    scalar hi = 2.0 * lo;
    while (f(hi) < 0.0)
    {
        lo = hi;
        hi *= 2.0;
        if (hi > 1e9) break; // the linear term wins long before this
    }
    for (int i = 0; i < 100 && hi - lo > 1e-13 * hi; ++i)
    {
        const scalar mid = 0.5 * (lo + hi);
        if (f(mid) < 0.0) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

scalar thermalTPlus(scalar yPlus, scalar Pr, scalar Prt, scalar yPlusSublayer)
{
    if (yPlus <= yPlusSublayer) return Pr * yPlus;
    return Prt / turb::kappa * std::log(turb::logE * yPlus)
         + Prt * jayatillekeP(Pr, Prt);
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
