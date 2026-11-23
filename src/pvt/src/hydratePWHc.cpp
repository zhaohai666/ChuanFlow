#include <cmath>

#include "cf/common/Error.h"
#include "cf/pvt/Flash.h"
#include "cf/pvt/HydratePWHc.h"

namespace cf::pvt
{

namespace
{

// Search domain of the equilibrium pressure. The lower bracket is effectively
// p -> 0 (theta -> 0 there, so the residual is exactly Delta_mu_w(T)); the
// ceiling is a declared-domain guard: a fully-occupied lattice saturates the
// stabilisation term, so an unreachably large Langmuir-to-pressure ratio means
// no physical root exists below the ceiling and the caller is refused loudly.
constexpr scalar kPressureFloor = 1.0e-3; // Pa
constexpr scalar kPressureCeiling = 1.0e10; // Pa

/// Residual of the vW-P equilibrium equation at (T, p), negative when the
/// hydrate stabilisation RT sum nu_j ln(1 - theta_j) overcompensates the
/// driving force. theta comes from the flashed gas fugacity at (T, p).
scalar hydratePWHcResidual(const Mixture& feed, scalar T, scalar p,
                           const HydratePWHcSpec& spec,
                           const std::vector<scalar>& nu,
                           scalar* thetaSmallOut, scalar* thetaLargeOut)
{
    const FlashResult r = tpFlash(feed, T, p);
    // Gas-side fugacity: vapour phase of a two-phase flash, the single-phase
    // state otherwise (its composition is the feed and its lnPhi is whatever
    // phase PR's density rule picked - declared simplification).
    const std::vector<scalar>& y = r.vapourComposition;
    const std::vector<scalar>& lnPhi = r.vapour.lnPhi;

    scalar theta[2] = {0.0, 0.0};
    for (std::size_t cage = 0; cage < 2; ++cage)
    {
        scalar sumCf = 0.0;
        for (std::size_t g = 0; g < spec.guests.size(); ++g)
        {
            const std::size_t i = feed.index(spec.guests[g]);
            const scalar f = std::exp(lnPhi[i]) * y[i] * p;
            sumCf += spec.langmuir[g][cage] * f;
        }
        theta[cage] = sumCf / (1.0 + sumCf);
    }
    if (thetaSmallOut) *thetaSmallOut = theta[0];
    if (thetaLargeOut) *thetaLargeOut = theta[1];

    const scalar stabilisation = kGasConstant * T
        * (nu[0] * std::log(1.0 - theta[0]) + nu[1] * std::log(1.0 - theta[1]));
    const scalar deltaMu = spec.muRef
        + spec.muSlopeA * (T - spec.referenceTemperature);
    return deltaMu + stabilisation;
}

} // namespace

scalar hydratePWHcCagePerWater(const std::string& structure, int cage)
{
    if (structure == "sI")
        return cage == 0 ? 2.0 / 46.0 : 6.0 / 46.0;
    if (structure == "sII")
        return cage == 0 ? 16.0 / 136.0 : 8.0 / 136.0;
    cfFatal("pvt: hydratePWHc unknown structure '" + structure
            + "' (expected sI or sII)");
    return 0.0;
}

HydratePWHcPoint hydratePWHcEquilibrium(const Mixture& feed, scalar T,
                                        const HydratePWHcSpec& spec)
{
    if (!(T > 0.0)) cfFatal("pvt: hydratePWHc needs a positive temperature");
    if (spec.guests.empty()) cfFatal("pvt: hydratePWHc needs at least one guest");
    if (spec.langmuir.size() != spec.guests.size())
        cfFatal("pvt: hydratePWHc needs Langmuir constants for every guest");
    for (const auto& row : spec.langmuir)
    {
        if (row.size() != 2)
            cfFatal("pvt: hydratePWHc needs two Langmuir constants per guest (small, large)");
        for (const scalar c : row)
            if (!(c > 0.0))
                cfFatal("pvt: hydratePWHc Langmuir constants must be positive");
    }
    if (!(spec.referencePressure > 0.0))
        cfFatal("pvt: hydratePWHc needs a positive referencePressure");
    if (!(spec.muRef > 0.0))
        cfFatal("pvt: hydratePWHc needs a positive muRef driving force");
    if (!(spec.muSlopeA > 0.0))
        cfFatal("pvt: hydratePWHc needs a positive muSlopeA");

    const scalar nuSmall = hydratePWHcCagePerWater(spec.structure, 0);
    const scalar nuLarge = hydratePWHcCagePerWater(spec.structure, 1);
    const std::vector<scalar> nu = {nuSmall, nuLarge};

    // Bracket: the residual falls monotonically as theta grows with p, so
    // double the upper bound from the declared anchor until it turns negative.
    scalar lo = kPressureFloor, hi = spec.referencePressure;
    scalar thetaS = 0.0, thetaL = 0.0;
    scalar gLo = hydratePWHcResidual(feed, T, lo, spec, nu, &thetaS, &thetaL);
    scalar gHi = hydratePWHcResidual(feed, T, hi, spec, nu, &thetaS, &thetaL);
    if (!(gLo > 0.0))
        cfFatal("pvt: hydratePWHc at T = " + std::to_string(T)
                + ": hydrate is stable at vanishing pressure (Delta_mu_w <= 0), "
                  "no equilibrium point");
    for (int step = 0; step < 60 && gHi > 0.0; ++step)
    {
        hi *= 2.0;
        if (hi > kPressureCeiling)
        {
            hi = kPressureCeiling;
            gHi = hydratePWHcResidual(feed, T, hi, spec, nu, &thetaS, &thetaL);
            if (gHi > 0.0)
                cfFatal("pvt: hydratePWHc at T = " + std::to_string(T)
                        + ": no equilibrium pressure at or below 1e10 Pa "
                          "(Langmuir driving ratio too weak)");
            break;
        }
        gHi = hydratePWHcResidual(feed, T, hi, spec, nu, &thetaS, &thetaL);
    }

    // Bisection on the residual: ~50 halvings reach the 12-digit report floor.
    for (int iter = 0; iter < 60 && hi - lo > 1e-12 * hi; ++iter)
    {
        const scalar mid = 0.5 * (lo + hi);
        const scalar gMid = hydratePWHcResidual(feed, T, mid, spec, nu, &thetaS, &thetaL);
        if (gMid > 0.0) lo = mid;
        else hi = mid;
    }

    HydratePWHcPoint out;
    out.pEq = 0.5 * (lo + hi);
    hydratePWHcResidual(feed, T, out.pEq, spec, nu, &thetaS, &thetaL);
    out.thetaSmall = thetaS;
    out.thetaLarge = thetaL;
    return out;
}

} // namespace cf::pvt
