#pragma once

#include <algorithm>
#include <cmath>

namespace cf
{

/// A generalized-Newtonian fluid, shared by the 1D case parser, the 3D FVM
/// kernel and the coupling interface. The wall constitutive law is
///
///     tau = tau_0 + K * gdot^n                       (Herschel-Bulkley)
///
/// with `flowIndex == 0` meaning "no rheology was declared": the fluid then
/// falls back to its own Newtonian `mu`, so n = 1, K = mu, tau_0 = 0 and the
/// law is exactly the one every case before this file ran. `consistency` is
/// in Pa.s^n -- at n != 1 it is NOT a viscosity, which is why the 1D parsers
/// refuse a `flowIndex` without it.
///
/// The same struct carries a *coupled edge's* rheology. There the numbers
/// describe the fluid inside the partner region, not a wall law the 1D row
/// evaluates: a coupled pipe's loss is its partner's statement, so the 1D
/// kernel must never charge it a second time (§5.43). `NetworkPipe` keeps the
/// declaration in `partnerRheology` and clears its own wall-law fields.
struct Rheology
{
    double flowIndex = 0;     // n; 0 = Newtonian (the fluid's own mu)
    double consistency = 0;   // K, Pa.s^n; required when flowIndex > 0
    double yieldStress = 0;   // tau_0, Pa; 0 = no threshold

    bool isNewtonian() const { return flowIndex == 0 && yieldStress == 0; }
    double index() const { return flowIndex > 0 ? flowIndex : 1.0; }
    double consistencyOr(double muNewtonian) const
    {
        return consistency > 0 ? consistency : muNewtonian;
    }
};

/// Apparent (dynamic) viscosity of the law at a shear rate, Pa.s. At low
/// shear rate a yielded fluid's tau_0/gdot diverges; Papanastasiou's
/// exponential regularisation keeps the algebraic law finite without
/// changing it wherever gdot >> 1/regularization:
///
///     mu_app = K gdot^(n-1) + tau_0 (1 - exp(-m gdot)) / gdot
///
/// `shearFloor` only bounds the arithmetic at gdot -> 0. For an undeclared
/// rheology the expression is K * gdot^0 / ... == muNewtonian exactly, so a
/// declared-but-Newtonian fluid is bit-identical to an undeclared one.
inline double apparentViscosityPaS(const Rheology& rheology, double gdot,
                                   double muNewtonian, double shearFloor,
                                   double regularization)
{
    const double g = std::max(gdot, shearFloor);
    const double n = rheology.index();
    double mu = rheology.consistencyOr(muNewtonian) * std::pow(g, n - 1.0);
    if (rheology.yieldStress > 0)
        mu += rheology.yieldStress
              * (-std::expm1(-regularization * g)) / g;
    return mu;
}

} // namespace cf
