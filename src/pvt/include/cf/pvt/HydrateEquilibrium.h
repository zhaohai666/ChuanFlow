#pragma once

#include "cf/common/Labels.h"
#include "cf/pvt/Components.h" // kGasConstant

namespace cf::pvt
{

/// Gas-hydrate dissociation equilibrium line in the Clausius–Clapeyron
/// (constant-enthalpy) limit, written about a declared reference point:
///
///   ln(p / pRef) = -(dH / R) (1/T - 1/Tref)
///   p_eq(T)      = pRef * exp[-(dH / R) (1/T - 1/Tref)]
///
/// with pRef the dissociation pressure (Pa) at Tref (K), dH the molar
/// dissociation enthalpy (J/mol) and R = kGasConstant. This is the same
/// integrated van't Hoff / Clausius–Clapeyron form used by WaxSolubility,
/// applied to the solid hydrate <-> water + gas equilibrium; it assumes ideal
/// gas and temperature-independent dH. p_eq = pRef at T = Tref and rises
/// steeply with T for positive dH.
///
/// Nothing is tabulated: Tref, pRef and dH are declared by the caller, so a
/// checker can reproduce the curve with nothing but math.exp and the three
/// declared constants. Distinct from the 1D kinetic hydrate model
/// (`TwoFluidHydrate`, which stores a curveA/curveB ln(p) line): this is a
/// standalone pvt-layer equilibrium function usable by any caller.
///
/// @param T     absolute temperature, K (must be > 0)
/// @param Tref  reference temperature, K (must be > 0)
/// @param pRef  reference dissociation pressure, Pa (must be > 0)
/// @param dH    molar dissociation enthalpy, J/mol
/// @return      equilibrium dissociation pressure, Pa
scalar hydrateEquilibriumPressure(scalar T, scalar Tref, scalar pRef, scalar dH);

/// Inverse: equilibrium temperature at a given pressure (same closed form),
/// useful for a temperature-driven kinetic model.
///   T_eq(p) = 1 / (1/Tref - (R/dH) ln(p/pRef))
scalar hydrateEquilibriumTemperature(scalar p, scalar Tref, scalar pRef, scalar dH);

} // namespace cf::pvt
