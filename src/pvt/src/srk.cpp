#include "cf/pvt/Srk.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "cf/common/Error.h"

namespace cf::pvt
{

namespace
{

constexpr scalar kOmegaA = 0.42748023354034; // SRK (1972) omega_a
constexpr scalar kOmegaB = 0.08664007007534; // SRK (1972) omega_b

/// Real roots of z^3 + c2 z^2 + c1 z + c0 = 0, the same Cardano/trigonometric
/// construction the PR path uses (guarded near the triple root at a pure
/// fluid's critical point).
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
        roots.push_back(shift);
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

SrkParams mixtureParametersSrk(const Mixture& mix, scalar T, scalar p)
{
    if (mix.size() == 0) cfFatal("pvt: empty mixture");
    if (!(T > 0.0)) cfFatal("pvt: temperature must be positive");
    if (!(p > 0.0)) cfFatal("pvt: pressure must be positive");

    const auto& comps = mix.components();
    const auto& x = mix.fractions();
    const std::size_t n = comps.size();
    if (x.size() != n) cfFatal("pvt: mixture fractions/components mismatch");

    SrkParams prm;
    prm.aI.resize(n);
    prm.bI.resize(n);
    prm.daIdT.resize(n);
    prm.aIJ.assign(n, std::vector<scalar>(n, 0.0));

    const scalar RT = kGasConstant * T;
    for (std::size_t i = 0; i < n; ++i)
    {
        const scalar w = comps[i].acentric;
        const scalar m = 0.480 + 1.574 * w - 0.176 * w * w;
        const scalar Tr = T / comps[i].Tc;
        const scalar sqTr = std::sqrt(Tr);
        const scalar alpha = (1.0 + m * (1.0 - sqTr)) * (1.0 + m * (1.0 - sqTr));
        const scalar Rtc = kGasConstant * comps[i].Tc;
        const scalar a0 = kOmegaA * Rtc * Rtc / comps[i].Pc;
        prm.aI[i] = a0 * alpha;
        prm.bI[i] = kOmegaB * Rtc / comps[i].Pc;
        // d(alpha)/dT = -m (1 + m (1 - sqrt(T/Tc))) / sqrt(T Tc)
        prm.daIdT[i] = a0 * (-m * (1.0 + m * (1.0 - sqTr))
                             / std::sqrt(T * comps[i].Tc));
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        prm.bMix += x[i] * prm.bI[i];
        for (std::size_t j = 0; j < n; ++j)
        {
            const scalar kij = mix.binaryInteraction(i, j);
            const scalar aij = std::sqrt(prm.aI[i] * prm.aI[j]) * (1.0 - kij);
            prm.aIJ[i][j] = aij;
            prm.aMix += x[i] * x[j] * aij;
            // d(a_ij)/dT on the geometric mean: 0.5 (a_i' a_j + a_i a_j') / a_ij^0
            const scalar denom = std::sqrt(prm.aI[i] * prm.aI[j]);
            const scalar daij = denom > 0.0
                ? (1.0 - kij)
                      * 0.5 * (prm.daIdT[i] * prm.aI[j] + prm.aI[i] * prm.daIdT[j])
                      / denom
                : 0.0;
            prm.daMixdT += x[i] * x[j] * daij;
        }
    }

    prm.A = prm.aMix * p / (RT * RT);
    prm.B = prm.bMix * p / (RT);
    return prm;
}

std::vector<scalar> compressibilityRootsSrk(const SrkParams& prm)
{
    const scalar A = prm.A, B = prm.B;
    // SRK: Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
    const std::vector<double> raw =
        cubicRealRoots(-1.0, A - B - B * B, -(A * B));

    std::vector<scalar> roots;
    for (const double z : raw)
        if (z > B * (1.0 + 1e-12)) roots.push_back(z);
    return roots;
}

EosState srkState(const Mixture& mix, scalar T, scalar p, RootPreference pref)
{
    const SrkParams prm = mixtureParametersSrk(mix, T, p);
    const std::vector<scalar> roots = compressibilityRootsSrk(prm);
    if (roots.empty()) cfFatal("pvt: SRK cubic has no physical root");

    const auto& x = mix.fractions();
    const std::size_t n = x.size();
    const scalar RT = kGasConstant * T;

    auto buildRoot = [&](scalar Z)
    {
        EosState st;
        st.Z = Z;
        st.molarVolume = Z * RT / p;
        st.numberDensityMolar = 1.0 / st.molarVolume;
        st.massDensity = mix.molarMass() / st.molarVolume;
        st.molarVolumeOverB = st.molarVolume / prm.bMix;
        st.lnPhi.resize(n);
        // SRK fugacity: ln phi_i = (b_i/b)(Z-1) - ln(Z-B)
        //   - (A/B)(2 sum_j x_j a_ij / a - b_i/b) ln(1 + B/Z)
        const scalar logBZ = std::log(1.0 + prm.B / Z);
        st.gibbs = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            scalar sumA = 0.0;
            for (std::size_t j = 0; j < n; ++j) sumA += x[j] * prm.aIJ[i][j];
            const scalar biB = prm.bI[i] / prm.bMix;
            st.lnPhi[i] = biB * (Z - 1.0) - std::log(Z - prm.B)
                        - prm.A / prm.B * (2.0 * sumA / prm.aMix - biB) * logBZ;
            st.gibbs += x[i] * st.lnPhi[i];
        }
        return st;
    };

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
            if (buildRoot(Z).gibbs < buildRoot(pick).gibbs) pick = Z;

    EosState best = buildRoot(pick);
    best.hadTwoRoots = distinct.size() > 1;
    best.lowestRoot = *std::min_element(distinct.begin(), distinct.end());
    best.highestRoot = *std::max_element(distinct.begin(), distinct.end());
    return best;
}

scalar srkDzDt(const Mixture& mix, scalar T, scalar p, RootPreference pref)
{
    const SrkParams prm = mixtureParametersSrk(mix, T, p);
    const std::vector<scalar> roots = compressibilityRootsSrk(prm);
    if (roots.empty()) cfFatal("pvt: SRK cubic has no physical root");

    scalar Z = roots.front();
    if (pref == RootPreference::Liquid)
        Z = *std::min_element(roots.begin(), roots.end());
    else if (pref == RootPreference::Vapor)
        Z = *std::max_element(roots.begin(), roots.end());
    else
    {
        // same Gibbs rule as srkState: dZ/dT must follow the root the state
        // reports, or the derivative is of a different branch
        const EosState st = srkState(mix, T, p, RootPreference::LowestGibbs);
        Z = st.Z;
    }

    const scalar RT = kGasConstant * T;
    const scalar A = prm.A, B = prm.B;
    const scalar dAdT = prm.daMixdT * p / (RT * RT) - 2.0 * A / T;
    const scalar dBdT = -B / T;
    const scalar fz = 3.0 * Z * Z - 2.0 * Z + (A - B - B * B);
    if (fz == 0.0) cfFatal("pvt: SRK dZ/dT hits a vanishing F_Z");
    // d/dT of c1 = A - B - B^2 is dA/dT - (1 + 2B) dB/dT: the linear B term
    // cannot be dropped just because the squared one is present.
    const scalar numerator =
        (dAdT - (1.0 + 2.0 * B) * dBdT) * Z - (A * dBdT + B * dAdT);
    return -numerator / fz;
}

scalar srkJouleThomson(const Mixture& mix, scalar T, scalar p, scalar cp,
                       RootPreference pref)
{
    if (!(cp > 0.0)) cfFatal("pvt: Joule-Thomson needs a positive c_p");
    const scalar dZ = srkDzDt(mix, T, p, pref);
    const scalar molarMass = mix.molarMass();
    if (!(molarMass > 0.0)) cfFatal("pvt: Joule-Thomson needs a molar mass");
    return kGasConstant * T * T * dZ / (p * molarMass * cp);
}

} // namespace cf::pvt
