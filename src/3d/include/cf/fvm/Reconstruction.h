#pragma once

#include <string>

#include "cf/common/Field.h"
#include "cf/common/Labels.h"
#include "cf/common/Vector.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{
namespace reconstruction
{

/// Slope limiter family for MUSCL face reconstruction. Every limiter reduces
/// to the exact slope on a linear field (both one-sided slopes equal) and is
/// TVD: the limited slope never exceeds the magnitude of the downwind slope.
enum class Limiter
{
    Minmod,   ///< min(|s1|,|s2|) -- the most diffusive, always TVD
    VanLeer,  ///< harmonic mean of the two slopes
    MC,       ///< monotonized central: minmod(2s1, 2s2, (s1+s2)/2)
    Superbee  ///< max(min(2s1,s2), min(s1,2s2))
};

/// Parse a limiter name; cfFatals on anything unknown.
Limiter parseLimiter(const std::string& name);
const char* limiterName(Limiter limiter);

/// Limiter function phi(r) in the audit's ratio form:
///   vanLeer(r)  = r <= 0 ? 0 : 2r/(1+r)
///   MC(r)       = r <= 0 ? 0 : min(2r, 0.5*(1+r), 2)
///   minmod(r)   = r <= 0 ? 0 : min(r, 1)
///   superbee(r) = r <= 0 ? 0 : max(min(2r,1), min(r,2))
/// phi(1) = 1 for every limiter, which is what makes linear fields exact.
scalar limiterPhi(scalar r, Limiter limiter);

/// Limited slope from three consecutive cell values along the flow direction:
///   a = far-upwind cell value, b = upwind cell value, c = downwind cell value.
/// The upwind-biased slope is s1 = b - a (the projected cell gradient gP*dx)
/// and the downwind slope is s2 = c - b (the face slope gC*L). The ratio
/// r = s2/s1 is exactly the audit's rP = gC/gP, so the returned value
/// phi(r)*s1 is the audit's limited slope phi(rP)*gP in discretised form.
/// Returns 0 when the values are not monotone (a sign change or a flat).
scalar limitedSlope(scalar a, scalar b, scalar c, Limiter limiter);

/// First-order upwind face values (reconstruction `none`): each face takes the
/// value of the cell it is upwind of; boundary faces use the owner value.
Field<scalar> upwindFaceValues(const polyMesh& mesh, const Field<scalar>& cellData,
                               const Field<scalar>& phi);

/// MUSCL second-order face values. For every internal face the upwind cell's
/// slope is limited with `limiter` and extrapolated to the face centre along
/// the owner-neighbour centre line; the result is clamped to the values of the
/// two cells sharing the face, so it is bounded on discontinuities. Boundary
/// faces use the owner (zeroGradient). Pure function: it never touches
/// fvMatrix or any equation assembly.
Field<scalar> musclFaceValues(const polyMesh& mesh, const Field<scalar>& cellData,
                              const Field<scalar>& phi, Limiter limiter);

} // namespace reconstruction
} // namespace cf
