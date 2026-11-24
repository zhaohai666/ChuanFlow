#pragma once

#include <vector>

#include "cf/pvt/Components.h"
#include "cf/pvt/Flash.h"
#include "cf/pvt/PengRobinson.h"

namespace cf::pvt
{

/// Non-isothermal flash building blocks (D13 increment 5). All heat-capacity
/// data are DECLARED by the caller (a linear ideal-gas cp(T) = a + b T per
/// component); nothing is tabulated in the library. The PR departure enthalpy
/// uses only declared EOS constants (Tc, Pc, omega, kij) via an analytic da/dT,
/// so it needs no extra input.

/// Ideal-gas molar enthalpy rise of one component from Tref to T, using a
/// declared linear cp(T) = a + b T (J/mol/K):
///   h_ig(T) - h_ig(Tref) = a (T - Tref) + b/2 (T^2 - Tref^2)
scalar idealGasEnthalpyRise(scalar T, scalar Tref, scalar cpA, scalar cpB);

/// Analytic da_i/dT of the PR alpha(T) function
///   alpha = [1 + kappa (1 - sqrt(T/Tc))]^2,
///   kappa = 0.37464 + 1.54226 w - 0.26992 w^2 (PR76 standard, omega <= 0.491 branch).
scalar prDaDT(const Component& c, scalar T);

/// Mixture departure (residual) enthalpy from the PR EOS at T,p for the
/// composition in `mix` and a chosen phase root `pref`:
///   H^R = R T (Z - 1) + (T daMix/dT - aMix)/(2 sqrt2 bMix) * ln[(Z+(1+sqrt2)B)/(Z+(1-sqrt2)B)]
/// Returns J/mol. Positive-T/p + a physical root required.
scalar prDepartureEnthalpy(const Mixture& mix, scalar T, scalar p,
                           RootPreference pref = RootPreference::LowestGibbs);

struct EnthalpyModel
{
    scalar referenceT = 298.15; // K, enthalpy zero point for the ideal-gas part
    std::vector<scalar> cpA;    // per component, J/(mol K)  (cp = a + b T)
    std::vector<scalar> cpB;    // per component, J/(mol K^2)
};

/// Molar enthalpy of a two-phase TP-flash result at T,p:
///   H = (1-beta) [sum x_i h_ig,i + H^R_liq] + beta [sum y_i h_ig,i + H^R_vap]
/// J/mol of feed. Requires cpA/cpB to hold one value per feed component.
scalar mixtureEnthalpy(const Mixture& feed, const FlashResult& flash,
                       scalar T, scalar p, const EnthalpyModel& model);

struct PhResult
{
    bool converged = false;
    int iterations = 0;
    scalar T = 0;          // solved temperature, K
    scalar enthalpy = 0;   // final molar enthalpy at (T, p), J/mol
    scalar residual = 0;   // |H - Htarget|
    FlashResult flash;     // TP flash at the solved T
    std::string message;
};

/// Isobaric isenthalpic (PH) flash: find T such that mixtureEnthalpy(T,p)=hTarget
/// at the given pressure. Brackets the temperature in [Tmin,Tmax] (bisection,
/// monotone H(T) in the single-phase regimes) and reuses tpFlash at each trial.
/// All heat-capacity data come from `model` (declared); no built-in constants.
PhResult phFlash(const Mixture& feed, scalar p, scalar hTarget,
                 const EnthalpyModel& model, scalar Tmin, scalar Tmax,
                 const FlashControls& controls = FlashControls{},
                 int maxOuterIter = 200, scalar tolH = 1.0e-6);

} // namespace cf::pvt
