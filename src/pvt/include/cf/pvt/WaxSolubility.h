#pragma once

#include "cf/common/Labels.h"
#include "cf/pvt/Components.h" // kGasConstant

namespace cf::pvt
{

/// Ideal-solution solid-liquid equilibrium (van't Hoff / Schröder–van Laar)
/// solubility of a single precipitating wax component:
///
///   ln x_sat(T) = -(dHf / R) (1/T - 1/Tf)
///   x_sat(T)    = exp[-(dHf / R) (1/T - 1/Tf)]
///
/// with Tf the pure-wax melting temperature (K), dHf its molar enthalpy of
/// fusion (J/mol) and R = kGasConstant. This is the ideal solution, dCp = 0
/// limit of the solid-liquid equilibrium line: the activity coefficient is
/// unity and the pressure dependence is neglected. x_sat = 1 at T = Tf,
/// x_sat < 1 below the melting point (supersaturation possible), and x_sat > 1
/// above it (the component is fully miscible).
///
/// Reference constants for n-tetracosane (n-C24H50, CAS 646-31-1), the wax
/// probe used by the verification case:
///   Tf   = 323.75 K   (NIST Chemistry WebBook, Domanska 1999, TRC +-0.3 K;
///                      NIST also reports the 35-value average Tfus = 324+-3 K)
///   dHf  = 54890 J/mol @ 324.1 K (NIST, Domalski & Hearing 1996; other
///                      tabulations: 53800 Mondieig 2004, 59310 Hafsaoui 2007)
///
/// The function is a pure, stateless transcription of the closed form; it is
/// deliberately independent of the PR flash so a checker can reproduce it
/// with nothing but math.exp and the two declared constants.
///
/// @param T        absolute temperature, K (must be > 0)
/// @param Tf       pure-wax melting temperature, K (must be > 0)
/// @param dHf      molar fusion enthalpy, J/mol
/// @return         equilibrium liquid mole fraction of the wax component,
///                 in (0, inf); clamp to 1 for a miscibility interpretation.
scalar waxSaturationMoleFraction(scalar T, scalar Tf, scalar dHf);

} // namespace cf::pvt
