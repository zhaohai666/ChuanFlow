#pragma once

#include <vector>

#include "cf/pvt/Components.h"
#include "cf/pvt/PengRobinson.h"

namespace cf::pvt
{

/// Stiel-Thodos (1961) dilute-gas viscosity of a pure nonpolar component,
/// returned in Pa s. T, Tc [K], Pc [Pa], MW [g/mol]:
///   Tr > 1.5 : mu[cP] = 17.78e-5 (4.58 Tr - 1.67)^0.625 / xi
///   else     : mu[cP] = 34e-5 Tr^0.94 / xi
///   xi = Tc^(1/6) / (sqrt(MW) Pc[atm]^(2/3)).
scalar stielThodosViscosity(scalar T, scalar Tc, scalar Pc, scalar MW);

/// Lohrenz-Bray-Clark (1964, SPE 915) residual viscosity, in the form given
/// by Whitson & Brule (2000) and reproduced by the `chemicals` library
/// (this implementation is validated against that library's worked example
/// to 1e-12 in tests/check_pvt_lbc.py):
///   mu[cP] = mu0[cP] + (poly^4 - 1e-4) / xi
///   poly = 0.1023 + 0.023364 rr + 0.058533 rr^2 - 0.040758 rr^3
///          + 0.0093724 rr^4
///   xi = Tc^(1/6) / (sqrt(MW) Pc[atm]^(2/3)),   rr = Vc / Vm
/// mu0 is the Herning-Zipperer sqrt(MW)-weighted mix of the per-component
/// Stiel-Thodos dilute-gas viscosities; Tc, Pc, MW, Vc are Kay (mole-fraction)
/// averages. Inputs: T [K], Vm [m3/mol], zs [-], MWs [g/mol], Tcs [K],
/// Pcs [Pa], Vcs [m3/mol]. Returns mu [Pa s].
scalar lbcViscosity(scalar T, scalar Vm, const std::vector<scalar>& zs,
                    const std::vector<scalar>& MWs,
                    const std::vector<scalar>& Tcs,
                    const std::vector<scalar>& Pcs,
                    const std::vector<scalar>& Vcs);

/// EOS-consistent convenience: the PR76 molar volume supplies Vm and PR's
/// universal critical compressibility (Zc = 0.3074) supplies the critical
/// volume Vc_i = Zc R Tc_i / Pc_i, so the reduced density needs no tabulated
/// critical volumes. `pref` selects the phase root (Liquid / Vapor).
scalar lbcViscosity(const Mixture& mix, scalar T, scalar p,
                    RootPreference pref = RootPreference::LowestGibbs);

} // namespace cf::pvt