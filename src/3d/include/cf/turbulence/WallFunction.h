#pragma once

#include "cf/common/Labels.h"

namespace cf
{

/// Wall-layer constants shared by the turbulence models and their tests.
namespace turb
{
inline constexpr scalar kappa = 0.41;   // von Karman constant
inline constexpr scalar logE = 9.8;     // log-law roughness constant, u+ = ln(E y+)/kappa
/// Spalding's intercept: the additive constant of u+ = ln(y+)/kappa + C.
inline constexpr scalar spaldingC = 5.5676;
} // namespace turb

/// Spalding's full-range law of the wall, written the way Spalding wrote it -
/// y+ as an explicit function of u+:
///   y+ = u+ + e^(-k C) [e^(k u+) - 1 - k u+ - (k u+)^2/2 - (k u+)^3/6]
/// It tends to y+ = u+ in the viscous sublayer and to the log law
/// u+ = ln(y+)/kappa + C for large u+, so one formula serves both resolved
/// (y+ < 1) and wall-function (y+ > 30) near-wall meshes.
scalar spaldingYPlus(scalar uPlus);

/// u+ at the given y+: the inverse of Spalding's law (bisection, since the
/// relation is strictly increasing).
scalar lawOfWall(scalar yPlus);

/// State of the wall layer at the first cell centre of a wall face.
struct WallState
{
    scalar yPlus = 0;  // u_tau y / nu
    scalar uTau = 0;   // friction velocity, sqrt(tau_wall/rho)
    scalar nut = 0;    // eddy viscosity at the wall face, u_tau y / u_parallel - nu
};

/// Invert the law of the wall: find u_tau > 0 with
///   uParallel = u_tau u+(u_tau y / nu)
/// Rather than invert the law and then invert again, the two relations are
/// combined into one monotone equation for u+:
///   u+ y+(u+) = uParallel y / nu
/// whose root gives u_tau = uParallel/u+ and y+ = uParallel y/(nu u+). In the
/// viscous limit y+ = u+ makes it u_parallel = u_tau^2 y/nu, so the wall-face
/// viscosity degenerates to the molecular one and a mesh that resolves the
/// sublayer is treated as such. nut = u_tau^2 y/u_parallel - nu is then the
/// eddy viscosity that makes the momentum diffusion through the wall face
/// equal to the prescribed shear.
/// Zero tangential speed, zero distance or zero viscosity give the trivial
/// state (no shear to report).
WallState wallState(scalar uParallel, scalar y, scalar nu);

} // namespace cf
