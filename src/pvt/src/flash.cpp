#include <algorithm>
#include <cmath>
#include <sstream>

#include "cf/common/Error.h"
#include "cf/pvt/Flash.h"

namespace cf::pvt
{

namespace
{

/// PR critical molar volume in units of b: v_c / b = Zc / Bc.
constexpr scalar kCriticalVolumeOverB = 0.30740 / 0.07779607390388846;

void normaliseInPlace(std::vector<scalar>& v)
{
    scalar sum = 0.0;
    for (const scalar x : v) sum += x;
    if (sum <= 0.0) cfFatal("pvt: flash produced a zero composition");
    for (scalar& x : v) x /= sum;
}

} // namespace

scalar rachfordRiceResidue(const std::vector<scalar>& z,
                           const std::vector<scalar>& K, scalar beta)
{
    scalar f = 0.0;
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        const scalar d = K[i] - 1.0;
        f += z[i] * d / (1.0 + beta * d);
    }
    return f;
}

scalar solveRachfordRice(const std::vector<scalar>& z, const std::vector<scalar>& K)
{
    // f is strictly decreasing on (0, 1) whenever f(0) > 0 > f(1), so
    // Newton guarded by a bisection bracket is unconditionally convergent.
    scalar lo = 0.0, hi = 1.0;
    scalar beta = 0.5;
    for (int iter = 0; iter < 100; ++iter)
    {
        const scalar f = rachfordRiceResidue(z, K, beta);
        if (f > 0.0) lo = beta;
        else hi = beta;

        scalar df = 0.0;
        for (std::size_t i = 0; i < z.size(); ++i)
        {
            const scalar d = K[i] - 1.0;
            const scalar w = 1.0 + beta * d;
            df -= z[i] * d * d / (w * w);
        }
        scalar next = (df != 0.0) ? beta - f / df : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        if (std::abs(next - beta) < 1e-15) { beta = next; break; }
        beta = next;
    }
    return beta;
}

bool isLiquidLike(const EosState& state)
{
    return state.molarVolumeOverB < kCriticalVolumeOverB;
}

FlashResult tpFlash(const Mixture& feed, scalar T, scalar p,
                    const FlashControls& controls)
{
    if (feed.size() == 0) cfFatal("pvt: flash of an empty mixture");
    const std::size_t n = feed.size();
    const std::vector<scalar> z = feed.fractions();

    FlashResult res;
    res.liquidComposition = z;
    res.vapourComposition = z;

    Mixture liq = feed, vap = feed;
    std::vector<scalar> K(n);
    for (std::size_t i = 0; i < n; ++i) K[i] = wilsonK(feed.components()[i], T, p);

    scalar beta = 0.5;
    scalar maxChange = 0.0;

    for (int iter = 1; iter <= controls.maxIter; ++iter)
    {
        const scalar f0 = rachfordRiceResidue(z, K, 0.0);
        const scalar f1 = rachfordRiceResidue(z, K, 1.0);

        std::vector<scalar> x(n), y(n);
        if (f0 <= 0.0)
        {
            // no vapour can form at this pressure: all liquid, vapour is incipient
            beta = 0.0;
            x = z;
            scalar sumKy = 0.0;
            for (std::size_t i = 0; i < n; ++i) { y[i] = K[i] * z[i]; sumKy += y[i]; }
            if (sumKy <= 0.0) cfFatal("pvt: flash produced a zero vapour composition");
            for (std::size_t i = 0; i < n; ++i) y[i] /= sumKy;
        }
        else if (f1 >= 0.0)
        {
            // past the dew point: all vapour, liquid is incipient
            beta = 1.0;
            y = z;
            scalar sumZk = 0.0;
            for (std::size_t i = 0; i < n; ++i) { x[i] = z[i] / K[i]; sumZk += x[i]; }
            if (sumZk <= 0.0) cfFatal("pvt: flash produced a zero liquid composition");
            for (std::size_t i = 0; i < n; ++i) x[i] /= sumZk;
        }
        else
        {
            beta = solveRachfordRice(z, K);
            for (std::size_t i = 0; i < n; ++i)
            {
                x[i] = z[i] / (1.0 + beta * (K[i] - 1.0));
                y[i] = K[i] * x[i];
            }
            normaliseInPlace(x);
            normaliseInPlace(y);
        }

        liq.setFractions(x);
        vap.setFractions(y);
        const EosState sl = pengRobinsonState(liq, T, p, RootPreference::Liquid);
        const EosState sv = pengRobinsonState(vap, T, p, RootPreference::Vapor);

        maxChange = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const scalar newK = std::exp(sl.lnPhi[i] - sv.lnPhi[i]);
            maxChange = std::max(maxChange, std::abs(std::log(newK / K[i])));
            K[i] = newK;
        }

        res.iterations = iter;
        res.residual = maxChange;
        if (maxChange < controls.tol)
        {
            res.converged = true;
            break;
        }
    }

    if (!res.converged)
    {
        std::ostringstream os;
        os << "flash did not reach |delta ln K| < " << controls.tol
           << " in " << controls.maxIter << " iterations (residual "
           << maxChange << ")";
        res.message = os.str();
    }

    // Rebuild the final phases from the converged K-values so the reported
    // states are exactly consistent with x, y and beta. The phase count is
    // decided here, not by how the last iteration was clamped: successive
    // substitution often starts clamped and leaves the clamp well before it
    // converges, and inside a one-phase region it collapses onto the trivial
    // K = 1 solution where the Rachford-Rice residue is only numerical noise.
    const scalar f0 = rachfordRiceResidue(z, K, 0.0);
    const scalar f1 = rachfordRiceResidue(z, K, 1.0);
    scalar spread = 0.0;
    for (const scalar k : K) spread = std::max(spread, std::abs(k - 1.0));

    bool twoPhase = false;
    if (f0 > 0.0 && f1 < 0.0 && spread > 1e-6)
    {
        beta = solveRachfordRice(z, K);
        if (beta > 0.0 && beta < 1.0)
        {
            std::vector<scalar> xt(n), yt(n);
            for (std::size_t i = 0; i < n; ++i)
            {
                xt[i] = z[i] / (1.0 + beta * (K[i] - 1.0));
                yt[i] = K[i] * xt[i];
            }
            normaliseInPlace(xt);
            normaliseInPlace(yt);
            scalar distinct = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                distinct = std::max(distinct, std::abs(yt[i] - xt[i]));
            twoPhase = distinct > 1e-8;
        }
    }

    std::vector<scalar> x(n), y(n);
    if (!twoPhase)
    {
        // Single phase: the feed stays as it is. Which label it carries is
        // decided by PR's own density, not by the (collapsed) K-values.
        x = z;
        y = z;
        res.singlePhase = true;
        const EosState one = pengRobinsonState(feed, T, p, RootPreference::LowestGibbs);
        beta = isLiquidLike(one) ? 0.0 : 1.0;
        res.liquid = one;
        res.vapour = one;
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            x[i] = z[i] / (1.0 + beta * (K[i] - 1.0));
            y[i] = K[i] * x[i];
        }
        normaliseInPlace(x);
        normaliseInPlace(y);
        liq.setFractions(x);
        vap.setFractions(y);
        res.liquid = pengRobinsonState(liq, T, p, RootPreference::Liquid);
        res.vapour = pengRobinsonState(vap, T, p, RootPreference::Vapor);
    }

    res.vaporFraction = beta;
    res.liquidComposition = x;
    res.vapourComposition = y;
    liq.setFractions(x);
    vap.setFractions(y);
    res.liquidMassDensity = res.liquid.massDensity;
    res.vapourMassDensity = res.vapour.massDensity;
    res.liquidMolarMass = liq.molarMass();
    res.vapourMolarMass = vap.molarMass();
    return res;
}

scalar bubblePressure(const Mixture& feed, scalar T)
{
    // A pure fluid's bubble point is its saturation pressure: Rachford-Rice
    // has no interior root for one component, so use fugacity equality.
    if (feed.size() == 1) return pureVaporPressure(feed.components()[0], T);

    // For a mixture the boundary is located on the flash verdict itself.
    // Successive substitution collapses onto the trivial K = 1 solution once
    // the feed is single phase, so the sum(K z) - 1 residue is unreliable at
    // pressures where no liquid root of the feed exists; "the flash finds two
    // distinct phases" is a monotone indicator of being below the bubble point.
    const auto twoPhaseAt = [&](scalar p) { return !tpFlash(feed, T, p).singlePhase; };

    scalar pEst = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        const Component& c = feed.components()[i];
        pEst += feed.fractions()[i] * c.Pc
              * std::exp(5.373 * (1.0 + c.acentric) * (1.0 - c.Tc / T));
    }
    if (!(pEst > 0.0) || !std::isfinite(pEst)) pEst = 1.0e6;

    scalar pTwo = pEst;
    bool found = twoPhaseAt(pTwo);
    for (int step = 0; step < 80 && !found; ++step)
    {
        pTwo *= 0.6;
        found = twoPhaseAt(pTwo);
    }
    if (!found) cfFatal("pvt: no two-phase region at or below the bubble estimate");

    scalar pOne = pTwo;
    found = false;
    for (int step = 0; step < 100 && !found; ++step)
    {
        pOne *= 1.25;
        found = !twoPhaseAt(pOne);
    }
    if (!found) cfFatal("pvt: could not bracket the bubble pressure from above");

    for (int iter = 0; iter < 80; ++iter)
    {
        const scalar pMid = 0.5 * (pTwo + pOne);
        if (twoPhaseAt(pMid)) pTwo = pMid;
        else pOne = pMid;
        if (pOne - pTwo < 1e-10 * pOne) break;
    }
    return 0.5 * (pOne + pTwo);
}

} // namespace cf::pvt
