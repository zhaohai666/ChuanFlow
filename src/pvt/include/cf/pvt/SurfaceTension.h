#pragma once

#include <vector>

#include "cf/pvt/Components.h"
#include "cf/pvt/Flash.h"

namespace cf::pvt
{

/// Weinaug & Katz (1943, Ind. Eng. Chem. 35:239) parachor surface tension of a
/// liquid-vapour mixture, in the SI form reproduced by the `chemicals` library
/// (this implementation is validated against that library's worked example to
/// 1e-12 in tests/check_pvt_surface_tension.py):
///   sigma[N/m] = [ sum_i P_i (x_i / Vml - y_i / Vmg) ]^4
/// P_i are parachors [N^0.25 m^2.75 / mol], Vml / Vmg the liquid / vapour
/// mixture molar volumes [m3/mol], x / y the phase mole fractions. The classic
/// cgs parachor [cm3/mol (dyn/cm)^0.25] converts as P_SI = P_cgs * 1e-6.75
/// (sigma: dyn/cm -> N/m is 1e-3, molar density: mol/cm3 -> mol/m3 is 1e6).
/// Vmg may be set very large to neglect the vapour density (low pressure).
/// Parachors are supplied per component; this unit transcribes only the
/// verified equation, not a tabulated parachor dataset.
scalar weinaugKatzSurfaceTension(const std::vector<scalar>& parachors,
                                 scalar Vml, scalar Vmg,
                                 const std::vector<scalar>& xs,
                                 const std::vector<scalar>& ys);

/// Flash-consistent convenience: tpFlash the feed at (T, p), take Vml / Vmg
/// from the phase molar mass / density (= PR molar volume) and x / y from the
/// flash, then apply Weinaug-Katz with the supplied per-component parachors
/// (feed order). A single-phase state has no liquid-vapour interface, so its
/// surface tension is 0.
scalar weinaugKatzSurfaceTension(const Mixture& mix, scalar T, scalar p,
                                 const std::vector<scalar>& parachors);

} // namespace cf::pvt