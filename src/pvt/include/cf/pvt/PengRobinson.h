#pragma once

#include <vector>

#include "cf/pvt/Components.h"

namespace cf::pvt
{

/// Peng-Robinson (1976) cubic equation of state with van der Waals one-fluid
/// mixing rules, implemented from the published defining equations:
///   p = RT/(v-b) - a(T) / (v(v+b) + b(v-b))
/// Everything downstream (density, fugacity coefficients, flash) is derived
/// from this single state function, so the 1D and 3D kernels share one set of
/// thermodynamic properties.

enum class RootPreference
{
    LowestGibbs, // stable root: minimise sum(x_i ln phi_i)
    Liquid,      // smallest real root > B
    Vapor        // largest real root
};

struct PengRobinsonParams
{
    scalar aMix = 0;              // Pa m6/mol2
    scalar bMix = 0;              // m3/mol
    std::vector<scalar> aI;       // per component
    std::vector<scalar> bI;
    // 5.230: the alpha-function ingredients, published per component so a
    // case can check the branch it got (prKappa) against the closed form
    // without re-deriving it from the roots.
    std::vector<scalar> kappaI;   // per-component kappa of the alpha function
    std::vector<scalar> alphaI;   // per-component alpha(T)
    std::vector<std::vector<scalar>> aIJ; // cross terms (symmetric)
    scalar A = 0;                 // dimensionless a p /(RT)^2
    scalar B = 0;                 // dimensionless b p /(RT)
};

struct EosState
{
    scalar Z = 0;
    scalar molarVolume = 0; // m3/mol
    scalar numberDensityMolar = 0; // mol/m3
    scalar massDensity = 0; // kg/m3
    scalar molarVolumeOverB = 0; // v/b: PR's critical value is Zc/Bc = 3.951
    std::vector<scalar> lnPhi;
    scalar gibbs = 0;       // sum(x_i ln phi_i), dimensionless residual
    bool hadTwoRoots = false;
    scalar lowestRoot = 0;  // smallest physical root at this state
    scalar highestRoot = 0; // largest physical root at this state
};

/// Wilson (1966) K-value estimate, K = (Pc/p) exp[5.373 (1+w)(1 - Tc/T)]:
/// the standard initial guess for a flash iteration.
scalar wilsonK(const Component& c, scalar T, scalar p);

/// 5.230: the alpha-function slope parameter kappa(omega), the single source
/// of the branch for BOTH a(T) (mixtureParameters) and da/dT (prDaDT), so the
/// enthalpy departure can never price a different alpha than the EOS does.
/// Peng & Robinson (1976) Eq 17 gives kappa = 0.37464 + 1.54226 w - 0.26992 w^2
/// for omega <= 0.491; Robinson & Peng (1978, GPA RR-28) extend it for the
/// heavier fractions with the cubic kappa = 0.379642 + 1.48503 w - 0.164423 w^2
/// + 0.016666 w^3, the published branch threshold being omega = 0.491. The
/// PR76 expression is evaluated with its original operation order below, so a
/// component at or under the threshold is bit-identical to the pre-5.230 code.
scalar prKappa(scalar omega);

/// a_i(T), b_i and the mixture parameters for mix.components()/fractions().
PengRobinsonParams mixtureParameters(const Mixture& mix, scalar T, scalar p);

/// Full state (Z, v, rho, fugacity coefficients) at T, p for the composition
/// stored in `mix`.
EosState pengRobinsonState(const Mixture& mix, scalar T, scalar p,
                           RootPreference pref = RootPreference::LowestGibbs);

/// Real roots > B of the compressibility cubic; `params` must already hold A,B.
std::vector<scalar> compressibilityRoots(const PengRobinsonParams& params);

/// Pure-component saturation pressure from PR fugacity equality
/// (ln phi_liquid = ln phi_vapour). Returns 0 when T >= Tc.
scalar pureVaporPressure(const Component& c, scalar T);

} // namespace cf::pvt
