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
///
/// Those three read the cell on its own, so on a cell-filled interface they read
/// the staircase: their direction does not converge as the mesh refines. `fit`
/// reads the neighbourhood instead, and reads it as volume-fraction data rather
/// than as positions. Every neighbour that holds the interface has one plane that
/// cuts its own fraction exactly; if the interface really is one plane, those
/// planes coincide, and the tilt that makes them coincide is the interface's. Fitting
/// that tilt - an area-weighted least square, solved by a few damped Gauss-Newton
/// rounds - uses what alpha says about the whole neighbourhood and none of what the
/// staircase says about one cell, so unlike a plane put through the reconstructed
/// interface *points* - whose sub-cell positions step by about a cell along a
/// cell-filled front, which is more than the arc they are meant to trace - it
/// converges. On a disk at R/dx = 8 the fitted normal is 1.2 deg off the analytic
/// radial where `mixed` is 8.3, and the two refinements after that take it to 0.49
/// and 0.30 while `mixed` sits at 7.6.
///
/// The damping is not there for a neighbourhood with nothing in it. A front crossing
/// the mesh diagonally puts its cut-face centroids nearly along one line, which leaves
/// one of the two tilt directions holding only the other's slope; an undamped fit walks
/// that way and stops where its neighbourhood disagrees with it by a tenth of a cell,
/// and the same disk's worst cell is 23 deg off the radial instead of 3.1. The rounds
/// stop where a curved front asks them to: past four the fit begins explaining the
/// staircase rather than the interface it is under, and at sixteen the coarse disk's
/// mean angle is back up at 2.9 deg.
///
/// The cell's own fraction then places the reported plane, which is what keeps the
/// read-back exact. On a front that really is one plane the four rounds leave nothing
/// to place: the normal recovered is the analytic one and its constant sits 5e-8 of a
/// cell off it, where the cell-local estimators are 8 deg away.
enum class NormalChoice
{
    Youngs,
    Centroid,
    Mixed,
    Fit,
};

/// The reconstruction of every band cell of an alpha field.
struct Reconstruction
{
    Field<Vector> normal;       ///< unit, towards increasing alpha, 0 outside the band
    Field<scalar> offset;       ///< plane constant of the reconstructed interface
    Field<scalar> area;         ///< interface area inside the cell
    Field<Vector> areaCentroid; ///< its centre: the sub-cell interface position
    Field<scalar> fraction;     ///< alpha as the reconstruction reads it back

    /// What the neighbourhood's planes still disagree by after the fit: the
    /// area-weighted root-mean-square distance between them, in length units. Zero is
    /// a neighbourhood a single plane explains exactly; on a curved interface it is
    /// the sagitta the plane model has to give up, which is what the fit cannot know.
    /// Only `fit` fills it.
    Field<scalar> fitScatter;

    /// Cells the interface is in: the band, by construction the same set the
    /// fields above are non-zero on.
    label bandCells = 0;

    /// The tolerance the band was taken at, so that what lies outside it - a cell
    /// that is full or one that is empty - is known without a second threshold that
    /// can drift from this one.
    scalar bandTolerance = 1e-6;

    /// What the same neighbourhood disagreed by before any tilt, averaged over the
    /// band: the starting normals from the gradient, the neighbours each cut to their
    /// own fraction, and the planes they put those fractions at. The fit drives this
    /// to `fitScatter`, so the pair says what the neighbourhood was worth. Only
    /// `fit` fills them.
    scalar fitScatterBefore = 0;

    /// The most rounds any cell's fit took before it stopped tilting, so a case can
    /// tell a fit that converged early from one that ran out of the budget.
    label fitRounds = 0;
};

/// Reconstruct the interface of `alpha`. A cell carries an interface when
/// `bandTolerance < alpha < 1 - bandTolerance`, which needs no threshold tuned
/// against the magnitude of grad(alpha) - the property that lets a cell-filled
/// interface, whose gradient is one cell thick, be reconstructed at all. The
/// normal is `choice`, refined once against the cut it produces, and `fraction`
/// is the volume fraction the final plane cuts, so `fraction - alpha` measures
/// how well the reconstruction reproduces the field.
///
/// For `fit` the normal is not the cell's own business: the gradient starts it, then
/// every neighbour that holds the interface is cut by the candidate plane to its own
/// fraction, and the candidate is tilted by the least-squares amount that makes those
/// cuts coincide. That needs the neighbours' bounds, so it happens after the band is
/// known, and it reads alpha - not the first plane's positions - for its information.
Reconstruction reconstruct(const polyMesh& mesh, const Field<scalar>& alpha,
                           const Field<Vector>& youngsGradient, NormalChoice choice,
                           scalar bandTolerance = 1e-6);

/// The same, taking the starting normal from the field's own gradient operator.
Reconstruction reconstruct(const volScalarField& alpha, const std::string& gradScheme,
                           NormalChoice choice, scalar bandTolerance = 1e-6);

/// A convex region of space, carried as closed outward-oriented planar face
/// loops. It starts as a cell, or as a cell cut by the reconstructed plane, and
/// the point of the type is that it can be cut again - by a list of half-spaces,
/// which is how much of a phase region lies inside a given box or swept prism is
/// measured. Every cut keeps the surface closed, so the volume stays an integral
/// over its own boundary rather than an assumption about the shape.
struct ConvexRegion
{
    std::vector<std::vector<Vector>> loops;
    bool empty = false; ///< the region is nothing: a cut removed all of it

    scalar volume() const;
    Vector centroid() const;
};

/// The whole of one cell as a region.
ConvexRegion cellRegion(const polyMesh& mesh, label cell);

/// The part of `region` on the `unitNormal . x >= offset` side of the plane.
ConvexRegion clipRegion(const ConvexRegion& region, const Vector& unitNormal,
                        scalar offset);

/// One half-space: the side `unitNormal . x >= offset`.
struct HalfSpace
{
    Vector unitNormal = Vector::zero();
    scalar offset = 0;
};

/// The intersection with every half-space in the list, applied in order.
ConvexRegion clipRegion(const ConvexRegion& region,
                        const std::vector<HalfSpace>& halfSpaces);

/// Piecewise-linear interface *advection*: alpha is carried by the planes, not by
/// a flux of alpha. The phase leaving a cell through a face in a step is the
/// volume that cell's reconstructed phase region shares with the prism the face
/// sweeps through `dt`, so a cell-filled interface stays cell-filled - what the
/// finite-volume transport of alpha cannot do, because a ramp of alpha moving
/// through an upwind flux spreads at the truncation error of the scheme. Each
/// face's flux is computed once, from the upwind cell's geometry, and is the same
/// number subtracted there and added to the neighbour, so on a flux field that
/// closes per cell the phase volume is conserved by construction rather than by
/// the solver's tolerance. A cell whose several outflow faces together would take
/// more than it holds gives them all the same fraction of what they asked for, and a
/// cell that would receive more than its own void plus what it is passing on gives
/// back the same fraction to every inflowing face - which is what lets a full cell in
/// the middle of the phase keep trading, instead of being told it has no room and
/// starving the step. The factor is a property of the face, so both cells see the same
/// number and the limiter costs boundedness without costing conservation, except where
/// two limiters meet along one face and leave a cell a hair outside [0,1]: that is
/// clipped, and the clip is counted rather than hidden.
///
/// Needs convex cells, like the reconstruction. The boundary is a reservoir whose
/// phase fraction is the patch's own alpha value: its interface is not
/// reconstructed, so a case with geometric transport has to keep the interface
/// off the patches that carry flux.
struct GeometricStep
{
    /// faces that moved phase
    label fluxFaces = 0;
    /// faces whose geometric flux the limiter reduced
    label limitedFaces = 0;
    /// rounds of the two limits it took to reach a fixed point: 1 when one pass
    /// already satisfies both bounds
    label limiterRounds = 0;
    /// the smallest factor the limiter applied; 1 when nothing was limited
    scalar minFactor = 1;
    /// largest |phi| dt / V over the faces that carried phase: the sweep as a
    /// fraction of the donor cell, the same measure the finite-volume step reports
    scalar maxCourant = 0;
    /// phase volume the [0, 1] clip took back, as a fraction of the whole
    scalar clampVolume = 0;
    /// sum(alpha V) after the step
    scalar phaseVolume = 0;
    /// its change over the step, which is round-off and nothing else
    scalar volumeError = 0;
    scalar minAlpha = 0, maxAlpha = 0;
};

/// One step of geometric transport of `alpha` by the flux `phi` (owner to
/// neighbour) over `dt`, using the planes of `r` as reconstructed from the field
/// as it is now.
GeometricStep advectGeometric(volScalarField& alpha, const Field<scalar>& phi,
                              scalar dt, const Reconstruction& r);

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
