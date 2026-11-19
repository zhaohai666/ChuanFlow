#include "cf/pvt/Enthalpy.h"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "cf/common/Error.h"

namespace cf::pvt
{

namespace
{
constexpr scalar kSqrt2 = 1.4142135623730951;
constexpr scalar kOmegaA = 0.45723552892138218;

/// daMix/dT from the van der Waals one-fluid rule with declared kij:
/// aMix = sum_ij x_i x_j sqrt(a_i a_j)(1-kij)
/// d/dT = sum_ij x_i x_j (1-kij) * 0.5 (da_i/dT sqrt(a_j/a_i) + da_j/dT sqrt(a_i/a_j))
scalar daMixDT(const Mixture& mix, const PengRobinsonParams& prm, scalar T)
{
    const auto& x = mix.fractions();
    const auto& comps = mix.components();
    const std::size_t n = comps.size();
    scalar out = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            const scalar kij = mix.binaryInteraction(i, j);
            const scalar ai = prm.aI[i], aj = prm.aI[j];
            if (!(ai > 0) || !(aj > 0)) continue;
            const scalar dai = prDaDT(comps[i], T), daj = prDaDT(comps[j], T);
            out += x[i] * x[j] * (1.0 - kij) * 0.5
                 * (dai * std::sqrt(aj / ai) + daj * std::sqrt(ai / aj));
        }
    }
    return out;
}
} // namespace

scalar idealGasEnthalpyRise(scalar T, scalar Tref, scalar cpA, scalar cpB)
{
    if (!(T > 0) || !(Tref > 0))
        cfFatal("pvt: idealGasEnthalpyRise needs positive T and Tref");
    return cpA * (T - Tref) + 0.5 * cpB * (T * T - Tref * Tref);
}

scalar prDaDT(const Component& c, scalar T)
{
    if (!(T > 0) || !(c.Tc > 0))
        cfFatal("pvt: prDaDT needs positive T and Tc");
    const scalar w = c.acentric;
    const scalar kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w;
    const scalar sqrtTr = std::sqrt(T / c.Tc);
    const scalar f = 1.0 + kappa * (1.0 - sqrtTr);
    // alpha = f^2, d alpha/dT = 2 f * (-kappa/(2 sqrt(T Tc))) = -kappa f / sqrt(T Tc)
    const scalar dadT = -kappa * f / std::sqrt(T * c.Tc);
    const scalar Rtc = kGasConstant * c.Tc;
    return kOmegaA * Rtc * Rtc / c.Pc * dadT;
}

scalar prDepartureEnthalpy(const Mixture& mix, scalar T, scalar p, RootPreference pref)
{
    const PengRobinsonParams prm = mixtureParameters(mix, T, p);
    const EosState st = pengRobinsonState(mix, T, p, pref);
    if (!(st.Z > 0) || !(prm.bMix > 0))
        cfFatal("pvt: prDepartureEnthalpy needs a physical root and positive bMix");
    const scalar dadT = daMixDT(mix, prm, T);
    const scalar B = prm.B; // bMix p /(R T)
    const scalar t1 = st.Z + (1.0 + kSqrt2) * B;
    const scalar t2 = st.Z + (1.0 - kSqrt2) * B;
    if (!(t1 > 0) || !(t2 > 0))
        cfFatal("pvt: prDepartureEnthalpy: non-positive log argument");
    const scalar logTerm = std::log(t1 / t2);
    return kGasConstant * T * (st.Z - 1.0)
         + (T * dadT - prm.aMix) / (2.0 * kSqrt2 * prm.bMix) * logTerm;
}

scalar mixtureEnthalpy(const Mixture& feed, const FlashResult& flash,
                       scalar T, scalar p, const EnthalpyModel& model)
{
    const std::size_t n = feed.size();
    if (model.cpA.size() != n || model.cpB.size() != n)
        cfFatal("pvt: EnthalpyModel cp arrays must hold one value per component");

    // Build the liquid / vapour compositions' departure enthalpies at (T,p).
    Mixture mL(feed.components()), mV(feed.components());
    mL.setFractions(flash.liquidComposition);
    mV.setFractions(flash.vapourComposition);
    // Copy declared kij onto both phase mixtures so departure uses the same
    // binary interaction table as the feed.
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
        {
            const scalar kij = feed.binaryInteraction(i, j);
            if (kij != 0.0)
            {
                mL.setBinaryInteraction(i, j, kij);
                mV.setBinaryInteraction(i, j, kij);
            }
        }
    const scalar hR_L = prDepartureEnthalpy(mL, T, p, RootPreference::Liquid);
    const scalar hR_V = prDepartureEnthalpy(mV, T, p, RootPreference::Vapor);

    scalar hL = 0.0, hV = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const scalar hg = idealGasEnthalpyRise(T, model.referenceT, model.cpA[i], model.cpB[i]);
        hL += flash.liquidComposition[i] * (hg + hR_L);
        hV += flash.vapourComposition[i] * (hg + hR_V);
    }
    const scalar beta = flash.vaporFraction;
    return (1.0 - beta) * hL + beta * hV;
}

PhResult phFlash(const Mixture& feed, scalar p, scalar hTarget,
                 const EnthalpyModel& model, scalar Tmin, scalar Tmax,
                 const FlashControls& controls, int maxOuterIter, scalar tolH)
{
    PhResult r;
    if (!(p > 0) || !(Tmin > 0) || !(Tmax > Tmin))
    {
        r.message = "phFlash needs positive p and 0 < Tmin < Tmax";
        return r;
    }

    // Enthalpy is monotone increasing with T for positive cp; bisection on
    // H(T) - hTarget. Each trial reuses the isothermal PR flash.
    scalar lo = Tmin, hi = Tmax;
    scalar hLo = 0.0, hHi = 0.0;
    {
        const FlashResult fl = tpFlash(feed, lo, p, controls);
        hLo = mixtureEnthalpy(feed, fl, lo, p, model);
    }
    {
        const FlashResult fl = tpFlash(feed, hi, p, controls);
        hHi = mixtureEnthalpy(feed, fl, hi, p, model);
    }
    if (!(hLo <= hTarget) || !(hTarget <= hHi))
    {
        r.message = "phFlash: target enthalpy outside [H(Tmin), H(Tmax)]";
        r.enthalpy = hTarget < hLo ? hLo : hHi;
        r.residual = std::abs(r.enthalpy - hTarget);
        return r;
    }

    scalar T = 0.5 * (lo + hi);
    for (int it = 0; it < maxOuterIter; ++it)
    {
        T = 0.5 * (lo + hi);
        const FlashResult fl = tpFlash(feed, T, p, controls);
        const scalar h = mixtureEnthalpy(feed, fl, T, p, model);
        const scalar res = std::abs(h - hTarget);
        if (res <= tolH)
        {
            r.converged = fl.converged;
            r.iterations = it + 1;
            r.T = T;
            r.enthalpy = h;
            r.residual = res;
            r.flash = fl;
            r.message = fl.converged ? "converged" : "outer converged but TP flash not converged";
            return r;
        }
        if (h < hTarget) lo = T; else hi = T;
    }
    r.iterations = maxOuterIter;
    r.T = T;
    r.message = "phFlash: outer bisection did not reach tolerance";
    return r;
}

} // namespace cf::pvt
