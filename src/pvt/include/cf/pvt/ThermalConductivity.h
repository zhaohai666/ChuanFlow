#pragma once

#include "cf/common/Labels.h"

namespace cf::pvt
{

/// Thermal conductivity from the Prandtl-number closure:
///
///   k = mu * cp / Pr
///
/// with mu the dynamic viscosity (Pa.s), cp the specific heat at constant
/// pressure (J/(kg.K)) and Pr the (molecular) Prandtl number (dimensionless).
/// The result k has SI units W/(m.K): (Pa.s)(J/kg/K) = (kg/m/s)(m2/s2/kg/K)
/// = m kg s^-3 K^-1 = W/(m.K). This is the standard turbulent/transport
/// relation k = rho*nu*cp/Pr written with mu = rho*nu.
///
/// All three inputs are declared by the caller (nothing is tabulated inside
/// the library, matching the project's "declare, never invent constants"
/// policy). The function is pure and stateless, so a checker can reproduce it
/// with a single multiply and divide.
///
/// @param mu   dynamic viscosity, Pa.s (must be > 0)
/// @param cp   specific heat, J/(kg.K) (must be > 0)
/// @param Pr   molecular Prandtl number (must be > 0)
/// @return     thermal conductivity k, W/(m.K)
scalar thermalConductivityFromPr(scalar mu, scalar cp, scalar Pr);

} // namespace cf::pvt
