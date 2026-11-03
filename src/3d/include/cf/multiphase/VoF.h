#pragma once

#include <string>

#include "cf/fvm/FvMatrix.h"

namespace cf::vof
{

/// Volume-of-Fluid first increment: the phase-fraction transport, and the
/// interface geometry a continuum surface force will need. Both are built from
/// the existing FVM operators and inherit their properties from them --
/// conservation because the transport uses the same face fluxes as the
/// momentum equation, boundedness because upwind convection of a
/// divergence-free flux is an M-matrix (row sum V/dt, non-positive
/// off-diagonals), so alpha stays in [0,1] with no clipping and no
/// redistribution.

struct Controls
{
    scalar tol = 1e-10;
    int maxIter = 2000;
};

/// `ddt(alpha) + div(phi, alpha) = 0` advanced by one backward-Euler step.
/// `phi` is the owner-to-neighbour face flux of a divergence-free velocity;
/// boundary faces take their alpha from the BCs already set on the field, so an
/// inflow patch has to carry the phase fraction it feeds.
SolverReport advectAlpha(volScalarField& alpha, const Field<scalar>& phi,
                         scalar dt, const Schemes& schemes, const Controls& ctl);

/// Volume occupied by phase 1, sum(alpha V).
scalar phaseVolume(const polyMesh& mesh, const Field<scalar>& alpha);

/// Alpha-weighted centre of the phase, i.e. the centroid of phase 1's volume.
Vector phaseCentroid(const polyMesh& mesh, const Field<scalar>& alpha);

/// Unit interface normal `grad(alpha)/|grad(alpha)|`. It points towards
/// increasing alpha, which for a droplet of phase 1 is inwards. A cell carries
/// an interface only if its |grad alpha| reaches `thresholdFraction` of the
/// field maximum; the rest get an exactly zero vector, which is what tells
/// `interfaceCurvature` which cells to read. The threshold is relative because
/// the magnitude of grad(alpha) is set by how the interface was initialised
/// (a cell-filled step and a tanh profile of the same interface differ by an
/// order of magnitude in it), while the *location* of the interface is not.
Field<Vector> interfaceNormal(const volScalarField& alpha, const std::string& scheme,
                              scalar thresholdFraction);

/// `kappa = div(n)` at cell centres, over the faces whose two cells both carry
/// an interface: the face normal is the linear interpolation of the two cell
/// normals, and each cell's own normal is subtracted before the dot product, so
/// a face between two cells that agree contributes nothing. Subtracting the cell
/// normal is exact in the continuum -- `n_cell * sum(Sf) = 0` for a closed cell
/// -- and it is what makes the zeroed normals of the bulk cells harmless: an
/// interface cell's face towards a non-interface neighbour simply drops out of
/// the stencil instead of injecting the +-|Sf|/(2 V) jump a plain
/// average-with-zero would put on either side of the band. For the same reason a
/// plane interface reads kappa = 0 to round-off, whatever the mesh.
///
/// Only faces *between* interface cells carry curvature, so the interface has to
/// span at least two cells in the direction along itself -- which a cell-filled
/// step does not (it is one cell thick across the band, and the normal has
/// nowhere to vary). An initialisation smoothed over a couple of cells does, and
/// the dynamics do not care how the step is smeared: the same interface advects.
///
/// A circular droplet of radius R then reads kappa = -1/R at cells a distance r
/// from the centre it really is -1/r, the difference being the band's radial
/// thickness; the sign follows the normal, which is what keeps a later CSF term
/// consistent.
Field<scalar> interfaceCurvature(const polyMesh& mesh, const Field<Vector>& normal);

} // namespace cf::vof
