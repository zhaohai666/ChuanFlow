#pragma once

#include <string>
#include <vector>

#include "cf/common/Vector.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PolyMesh.h"

namespace cf::vof
{

/// Piecewise-linear interface calculation (PLIC): the interface inside a cell is
/// not a smeared ramp of alpha but a *plane*, positioned so that the plane cuts
/// the cell's volume into the fraction the field says it holds. Everything here
/// is therefore geometry, not transport: given alpha it returns the normal, the
/// plane and the sub-cell position of the interface, and it is exact for a
/// cell-filled (one-cell-thick, unsmoothed) field, which is the state a VOF
/// field is actually in after advection.
///
/// The clipping is done by cutting the cell's own boundary loops with the plane
/// and integrating over the closed result, so it holds for any convex cell -
/// hexahedra, tetrahedra and prisms - and needs no assumption about the mesh
/// being structured. Non-convex cells are refused rather than silently wrong.

/// A face of a cell: the vertex loop wound counter-clockwise as seen from
/// outside the cell, and the outward area vector (its magnitude is the area).
struct CellFace
{
    std::vector<Vector> loop;
    Vector areaVector = Vector::zero();
};

/// The oriented boundary of one cell. CellFace::loop is rewound so that the
/// fan triangulation of the loop has its normal along areaVector, which is the
/// outward normal of the face as seen from this cell (a face stored the other
/// way round - the cell being the neighbour - is reversed here).
std::vector<CellFace> cellBoundary(const polyMesh& mesh, label cell);

/// The solid left when a convex cell is cut by the plane `n . x = offset`,
/// keeping the side `n . x >= offset`.
struct Cut
{
    scalar volume = 0;
    Vector centroid = Vector::zero();
    /// Area of the planar cut face; 0 when the cell is whole or gone.
    scalar area = 0;
    /// Centre of that face, i.e. the sub-cell position of the interface.
    Vector areaCentroid = Vector::zero();
    bool whole = false; ///< the kept side is all of the cell
    bool void_ = false; ///< the kept side is empty
};

/// Cut `cell` by the plane and report the kept side. `unitNormal` must be
/// normalised; the cell is the oriented boundary from `cellBoundary`.
Cut clipConvexCell(const std::vector<CellFace>& cell, const Vector& unitNormal,
                   scalar offset);

/// The plane constant whose `n . x >= offset` side holds `fraction` of the cell.
/// Monotone in offset (the kept volume falls as the plane sweeps through), so a
/// bracketed Newton on the cut - whose derivative is exactly the cut area -
/// reaches machine precision in a handful of steps. `fraction` outside [0,1] is
/// refused.
scalar offsetForFraction(const std::vector<CellFace>& cell, const Vector& unitNormal,
                         scalar fraction);

/// Which normal the reconstruction is built with. `youngs` is the divergence of
/// alpha at cell centres (the estimate the smeared-band path uses), `centroid`
/// is the direction from the centroid of the phase-0 part of the cell to the
/// centroid of its phase-1 part, and `mixed` is the normalised mean of the two
/// - the combination that keeps working when alpha is close to 0 or 1, where the
/// centroid difference degenerates.
enum class NormalChoice
{
    Youngs,
    Centroid,
    Mixed,
};

/// The reconstruction of every band cell of an alpha field.
struct Reconstruction
{
    Field<Vector> normal;       ///< unit, towards increasing alpha, 0 outside the band
    Field<scalar> offset;       ///< plane constant of the reconstructed interface
    Field<scalar> area;         ///< interface area inside the cell
    Field<Vector> areaCentroid; ///< its centre: the sub-cell interface position
    Field<scalar> fraction;     ///< alpha as the reconstruction reads it back

    /// Cells the interface is in: the band, by construction the same set the
    /// fields above are non-zero on.
    label bandCells = 0;
};

/// Reconstruct the interface of `alpha`. A cell carries an interface when
/// `bandTolerance < alpha < 1 - bandTolerance`, which needs no threshold tuned
/// against the magnitude of grad(alpha) - the property that lets a cell-filled
/// interface, whose gradient is one cell thick, be reconstructed at all. The
/// normal is `choice`, refined once against the cut it produces, and `fraction`
/// is the volume fraction the final plane cuts, so `fraction - alpha` measures
/// how well the reconstruction reproduces the field.
Reconstruction reconstruct(const polyMesh& mesh, const Field<scalar>& alpha,
                           const Field<Vector>& youngsGradient, NormalChoice choice,
                           scalar bandTolerance = 1e-6);

/// The same, taking the starting normal from the field's own gradient operator.
Reconstruction reconstruct(const volScalarField& alpha, const std::string& gradScheme,
                           NormalChoice choice, scalar bandTolerance = 1e-6);

/// Alpha of a disk extruded along z - the phase region is
/// `(x - c.x)^2 + (y - c.y)^2 <= R^2`, unconstrained in z - filled into every
/// cell by its EXACT overlap area, not by whether the cell centre falls inside.
/// This is what makes the field a VOF field of that geometry: alpha is then the
/// volume fraction a reconstruction has to recover, so the interface position,
/// normal and curvature measured from it are the scheme's, not the
/// initialisation's. Needs axis-aligned cells (the cell's extent in x and y has
/// to be its bounding box).
Field<scalar> exactDiskAlpha(const polyMesh& mesh, Vector centre, scalar radius);

/// The exact area of the part of a disk of radius `radius` centred at the origin
/// that falls inside the axis-aligned rectangle [x1,x2] x [y1,y2]. Analytic: the
/// integral of the chord is split at the x where the circle crosses the two
/// horizontal sides, and each piece is taken with its closed-form primitive.
scalar diskRectangleArea(scalar x1, scalar x2, scalar y1, scalar y2, scalar radius);

} // namespace cf::vof
