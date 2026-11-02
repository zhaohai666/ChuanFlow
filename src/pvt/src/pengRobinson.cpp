#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "cf/common/Error.h"
#include "cf/pvt/PengRobinson.h"

namespace cf::pvt
{

namespace
{

constexpr scalar kSqrt2 = 1.4142135623730951;
constexpr scalar kOmegaA = 0.45723552892138218; // PR (1976) omega_a
constexpr scalar kOmegaB = 0.07779607390388846; // PR (1976) omega_b

/// Real roots of z^3 + c2 z^2 + c1 z + c0 = 0 (Cardano / trigonometric).
/// The three-real-root branch is guarded for the near-triple root at a pure
/// fluid's critical point, where the discriminant is ~0 to round-off.
std::vector<double> cubicRealRoots(double c2, double c1, double c0)
{
    const double p = c1 - c2 * c2 / 3.0;
    const double q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0;
    const double disc = q * q / 4.0 + p * p * p / 27.0;
    const double shift = -c2 / 3.0;

    std::vector<double> roots;
    if (disc > 0.0)
    {
        const double s = std::sqrt(disc);
        roots.push_back(std::cbrt(-q / 2.0 + s) + std::cbrt(-q / 2.0 - s) + shift);
        return roots;
    }
    if (p >= 0.0)
    {
        roots.push_back(shift); // triple root
        return roots;
    }
    const double r = std::sqrt(-p * p * p / 27.0);
    const double m = 2.0 * std::sqrt(-p / 3.0);
    const double theta = std::acos(std::clamp(-q / (2.0 * r), -1.0, 1.0));
    for (int k = 0; k < 3; ++k)
        roots.push_back(m * std::cos((theta + 2.0 * k * std::numbers::pi) / 3.0) + shift);
    std::sort(roots.begin(), roots.end());
    return roots;
}

} // namespace

scalar wilsonK(const Component& c, scalar T, scalar p)
{
    return c.Pc / p * std::exp(5.373 * (1.0 + c.acentric) * (1.0 - c.Tc / T));
}

PengRobinsonParams mixtureParameters(const Mixture& mix, scalar T, scalar p)
{
    if (mix.size() == 0) cfFatal("pvt: empty mixture");
    if (!(T > 0.0)) cfFatal("pvt: temperature must be positive");
    if (!(p > 0.0)) cfFatal("pvt: pressure must be positive");

    const auto& comps = mix.components();
    const auto& x = mix.fractions();
    const std::size_t n = comps.size();
    if (x.size() != n) cfFatal("pvt: mixture fractions/components mismatch");

    PengRobinsonParams prm;
    prm.aI.resize(n);
    prm.bI.resize(n);
    prm.aIJ.assign(n, std::vector<scalar>(n, 0.0));

    const scalar RT = kGasConstant * T;
    for (std::size_t i = 0; i < n; ++i)
    {
        const scalar w = comps[i].acentric;
        const scalar kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w;
        const scalar alpha = (1.0 + kappa * (1.0 - std::sqrt(T / comps[i].Tc)))
                           * (1.0 + kappa * (1.0 - std::sqrt(T / comps[i].Tc)));
        // a_i = 0.45724 (R Tc)^2 / Pc * alpha(T),  b_i = 0.07780 R Tc / Pc
        const scalar Rtc = kGasConstant * comps[i].Tc;
        prm.aI[i] = kOmegaA * Rtc * Rtc / comps[i].Pc * alpha;
        prm.bI[i] = kOmegaB * Rtc / comps[i].Pc;
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        prm.bMix += x[i] * prm.bI[i];
        for (std::size_t j = 0; j < n; ++j)
        {
            // van der Waals one-fluid rule, binary interaction kij = 0
            const scalar aij = std::sqrt(prm.aI[i] * prm.aI[j]);
            prm.aIJ[i][j] = aij;
            prm.aMix += x[i] * x[j] * aij;
        }
    }

    prm.A = prm.aMix * p / (RT * RT);
    prm.B = prm.bMix * p / (RT);
    return prm;
}

std::vector<scalar> compressibilityRoots(const PengRobinsonParams& prm)
{
    const scalar A = prm.A, B = prm.B;
    const std::vector<double> raw = cubicRealRoots(
        -(1.0 - B), A - 3.0 * B * B - 2.0 * B, -(A * B - B * B - B * B * B));

    std::vector<scalar> roots;
    for (const double z : raw)
        if (z > B * (1.0 + 1e-12)) roots.push_back(z);
    return roots;
}

EosState pengRobinsonState(const Mixture& mix, scalar T, scalar p, RootPreference pref)
{
    const PengRobinsonParams prm = mixtureParameters(mix, T, p);
    const std::vector<scalar> roots = compressibilityRoots(prm);
    if (roots.empty())
        cfFatal("pvt: Peng-Robinson cubic has no physical root");

    const auto& x = mix.fractions();
    const std::size_t n = x.size();
    const scalar RT = kGasConstant * T;
    const scalar logTermBase = 2.0 * kSqrt2 * prm.B;

    auto buildRoot = [&](scalar Z)
    {
        EosState st;
        st.Z = Z;
        st.molarVolume = Z * RT / p;
        st.numberDensityMolar = 1.0 / st.molarVolume;
        st.massDensity = mix.molarMass() / st.molarVolume;
        st.molarVolumeOverB = st.molarVolume / prm.bMix;
        st.lnPhi.resize(n);
        const scalar logTerm =
            std::log((Z + (1.0 + kSqrt2) * prm.B) / (Z + (1.0 - kSqrt2) * prm.B));
        st.gibbs = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            scalar sumA = 0.0;
            for (std::size_t j = 0; j < n; ++j) sumA += x[j] * prm.aIJ[i][j];
            const scalar biB = prm.bI[i] / prm.bMix;
            st.lnPhi[i] = biB * (Z - 1.0) - std::log(Z - prm.B)
                        - prm.A / logTermBase * (2.0 * sumA / prm.aMix - biB) * logTerm;
            st.gibbs += x[i] * st.lnPhi[i];
        }
        return st;
    };

    // Root selection: at fixed T, p and composition the stable phase is the
    // root that minimises sum(x_i ln phi_i).
    std::vector<scalar> distinct;
    for (const scalar Z : roots)
        if (std::find(distinct.begin(), distinct.end(), Z) == distinct.end())
            distinct.push_back(Z);

    scalar pick = distinct.front();
    if (pref == RootPreference::Liquid)
        pick = *std::min_element(distinct.begin(), distinct.end());
    else if (pref == RootPreference::Vapor)
        pick = *std::max_element(distinct.begin(), distinct.end());
    else
        for (const scalar Z : distinct)
        {
            const EosState cand = buildRoot(Z);
            if (cand.gibbs < buildRoot(pick).gibbs) pick = Z;
        }

    EosState best = buildRoot(pick);
    best.hadTwoRoots = distinct.size() > 1;
    best.lowestRoot = *std::min_element(distinct.begin(), distinct.end());
    best.highestRoot = *std::max_element(distinct.begin(), distinct.end());
    return best;
}

scalar pureVaporPressure(const Component& c, scalar T)
{
    if (T >= c.Tc) return 0.0;

    Mixture pure;
    pure.add(c, 1.0);

    auto g = [&](scalar p)
    {
        return pengRobinsonState(pure, T, p, RootPreference::Liquid).lnPhi[0]
             - pengRobinsonState(pure, T, p, RootPreference::Vapor).lnPhi[0];
    };

    // Wilson estimate as the starting point, then a trusted-region Newton on
    // g(p) = ln phi_liq - ln phi_vap.
    scalar p = c.Pc * std::exp(5.373 * (1.0 + c.acentric) * (1.0 - c.Tc / T));
    if (!(p > 0.0) || !std::isfinite(p)) p = 1.0e4;

    scalar pg = g(p);
    for (int iter = 0; iter < 200 && std::abs(pg) > 1e-13; ++iter)
    {
        const scalar dp = std::max(1.0, 1e-4 * p);
        const scalar slope = (g(p + dp) - pg) / dp;
        if (!std::isfinite(slope) || slope == 0.0) break;
        // Newton with a trust region of a factor of 2 in ln p
        const scalar step = std::clamp(-pg / slope, -0.5 * p, 0.5 * p);
        p = std::max(p + step, 1.0);
        pg = g(p);
        if (std::abs(step) < 1e-13 * p) break;
    }
    return p;
}

} // namespace cf::pvt
