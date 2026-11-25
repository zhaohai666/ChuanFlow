#pragma once

#include <vector>

#include "cf/pvt/Components.h"
#include "cf/pvt/PengRobinson.h" // EosState / RootPreference are the shared shapes

namespace cf::pvt
{

/// Soave-Redlich-Kwong (1972) cubic EOS with van der Waals one-fluid mixing,
/// implemented from the published defining equations:
///   p = RT/(v - b) - a(T) / (v (v + b))
/// with a_i = 0.42747 (R Tc)^2/Pc * alpha(T), b_i = 0.08664 R Tc/Pc and
///   alpha = (1 + m (1 - sqrt(T/Tc)))^2,  m = 0.480 + 1.574 w - 0.176 w^2.
/// SRK's own critical compressibility is Zc = 1/3 (its cubic is
/// Z^3 - Z^2 + (A - B - B^2) Z - A B = 0), which is what the unit tests anchor
/// against the EOS definition instead of against a fitted table.
///
/// Everything here is derived from the state function: the cubic roots, the
/// fugacity coefficients, the analytic compressibility derivative dZ/dT and
/// the Joule-Thomson coefficient that follows from it. The grid of the
/// implementation (root filtering, Gibbs selection) is deliberately the same
/// as PengRobinson so a caller can swap one for the other without changing
/// what "the stable root" means.

struct SrkParams
{
    scalar aMix = 0;              // Pa m6/mol2
    scalar bMix = 0;              // m3/mol
    std::vector<scalar> aI;       // per component
    std::vector<scalar> bI;
    std::vector<std::vector<scalar>> aIJ; // symmetric cross terms
    std::vector<scalar> daIdT;    // per component, Pa m6/(mol2 K)
    scalar daMixdT = 0;           // Pa m6/(mol2 K)
    scalar A = 0;                 // dimensionless a p / (RT)^2
    scalar B = 0;                 // dimensionless b p / (RT)
};

/// a_i(T), b_i(T), the mixture parameters and their T-derivatives.
SrkParams mixtureParametersSrk(const Mixture& mix, scalar T, scalar p);

/// The physical compressibility roots (Z > B) of the SRK cubic.
std::vector<scalar> compressibilityRootsSrk(const SrkParams& prm);

/// Full state at T, p for the composition stored in `mix`, with the same
/// root-selection contract as pengRobinsonState.
EosState srkState(const Mixture& mix, scalar T, scalar p,
                  RootPreference pref = RootPreference::LowestGibbs);

/// Analytic dZ/dT at constant p and composition, on the root the preference
/// picks. Differentiating the cubic implicitly:
///   F(Z, T) = Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
///   dZ/dT = -[(dA/dT - (1 + 2B) dB/dT) Z - (A dB/dT + B dA/dT)] / F_Z
/// with dB/dT = -B/T and dA/dT = p/(R T)^2 daMix/dT - 2A/T. The 1 in front of
/// 2B is the c1 = A - B - B^2 derivative, which is easy to lose.
scalar srkDzDt(const Mixture& mix, scalar T, scalar p,
               RootPreference pref = RootPreference::LowestGibbs);

/// Joule-Thomson coefficient mu = (T (dv/dT)_p - v) / c_p in K/Pa, on the
/// SRK state, for a MASS-specific c_p (J/(kg K), the form the 1D thermal
/// march carries). With v_s = Z R T / (p M) the two Z terms cancel and
///   mu = R T^2 dZ/dT / (p M c_p),
/// i.e. the EOS's own derivative, not a declared constant and not a finite
/// difference hidden inside the call. The molar mass is not optional: omit
/// it and the coefficient is off by a factor of M.
scalar srkJouleThomson(const Mixture& mix, scalar T, scalar p, scalar cp,
                       RootPreference pref = RootPreference::LowestGibbs);

} // namespace cf::pvt
