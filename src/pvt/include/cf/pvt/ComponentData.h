#pragma once

#include <string>

#include "cf/common/Error.h"
#include "cf/pvt/Components.h"

namespace cf::pvt
{

/// Opt-in built-in physical-property tables for components.
///
/// POLICY: nothing here is ever used implicitly. `Component` default
/// construction and `Mixture` never read these tables; a case must ask for
/// them explicitly (e.g. `surfaceTension { parachors fromTable; }`). This keeps
/// the standing invariant that an undeclared parachor / kij path is
/// bit-identical to the legacy path.
///
/// Provenance of the numbers
/// -------------------------
/// Parachors are the published Weinaug & Katz (1943) pure-substance values as
/// reproduced in Table 18.1 of Adewumi, "Phase Relations in Reservoir
/// Engineering" (Penn State / LibreTexts, §18.7 Surface Tension), which cites
/// the original: H. T. Weinaug and D. L. Katz, "Surface Tension of Methane-
/// Propane Mixtures", Ind. Eng. Chem. 35 (2): 239-246 (1943), DOI
/// 10.1021/ie50407a023. Units: cgs parachor [cm3/mol (dyn/cm)^0.25]. The SI
/// conversion P_SI = P_cgs * 1e-6.75 is applied by the caller exactly as it
/// is for a declared parachor (see SurfaceTension.h), so a table value and an
/// equal declared value produce the same surface tension.
///
/// Table:
///   CO2 78.0 | N2 41.0 | C1 77.0 | C2 108.0 | C3 150.3 | iC4 181.5
///   nC4 189.9 | iC5 225.0 | nC5 231.5 | nC6 271.0 | nC7 312.5 | nC8 351.5
///
/// Binary-interaction parameters (kij) are deliberately NOT tabulated here:
/// there is no single published kij table that is independent of the mixing
/// rule and temperature, so kij stays a declared, per-case input. Inventing a
/// "standard" kij table would violate the verifiable-source requirement.

/// Look up the Weinaug-Katz parachor (cgs) by dictionary name. Returns found
/// = false for a component not in the published table rather than guessing.
scalar parachorFromTableCgs(const std::string& name, bool& found);

/// Throwing variant for user-facing case paths: a typo or a component absent
/// from the published table is a loud error, never a silent default.
scalar requireParachorCgs(const std::string& name);

/// Build the SI parachor vector in feed order from the published table.
/// `P_SI = P_cgs * 1e-6.75` (see SurfaceTension.h). cfFatal if any component
/// is missing from the table.
std::vector<scalar> parachorsFromTableSi(const std::vector<Component>& comps);

} // namespace cf::pvt
