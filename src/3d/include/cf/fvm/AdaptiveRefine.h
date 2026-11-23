#pragma once

#include <optional>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Field.h"
#include "cf/common/Vector.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Analytic scalar field that drives adaptive refinement.
///
/// Two profiles are supported (the case key is `field linear|step`):
///   linear : phi(x) = dot(gradient, x) + offset
///   step   : phi(x) = low if x[axis] < interface else high
///
/// Both profiles have closed-form cell averages and box integrals, so the
/// integral of phi over the assembled mesh is independent of how it is
/// refined. That is the quantity the verifier checks for conservation.
class AdaptField
{
public:
    enum class Kind { Linear, Step };

    Kind kind = Kind::Linear;
    Vector gradient = Vector::zero(); // linear
    scalar offset = 0;                // linear
    int axis = 0;                     // step
    scalar face = 0;                  // step: coordinate of the interface
    scalar low = 0;                   // step
    scalar high = 1;                  // step

    /// Point value of the analytic field.
    scalar at(const Vector& x) const;

    /// Exact average of the field over the axis-aligned box [lo, hi].
    scalar cellAverage(const Vector& lo, const Vector& hi) const;

    /// Exact integral of the field over the axis-aligned box [lo, hi].
    scalar boxIntegral(const Vector& lo, const Vector& hi) const;

    /// Directional gradient magnitude along axis a: the change of phi between
    /// the two faces normal to a, divided by that cell edge. For a linear
    /// field this is |gradient[a]| for every cell; for a step field it is
    /// (high-low)/edge on straddling cells and zero elsewhere.
    scalar directionalGradient(const Vector& lo, const Vector& hi, int a) const;

    static AdaptField fromDict(const Dictionary& adapt, const Vector& boxMin,
                               const Vector& boxMax);
};

/// Refinement controls and the field, parsed from `adapt{}` plus the `block{}`
/// that supplies the initial Cartesian hex mesh.
struct AdaptControls
{
    scalar threshold = 0;
    label maxDepth = 0;
    AdaptField field;
    bool dumpCells = false;

    static AdaptControls fromDict(const Dictionary& caseDict, const Vector& boxMin,
                                  const Vector& boxMax);
};

/// One leaf cell of the refined Cartesian mesh. Cells stay axis-aligned boxes
/// so a conforming polyMesh (identical faces on both sides of an internal
/// face) can be rebuilt without hanging nodes.
struct AdaptCell
{
    Vector lo = Vector::zero();
    Vector hi = Vector::zero();
    label level = 0;
    scalar phi = 0;      // cell-average of the analytic field
    scalar volume = 0;   // box volume
};

/// Summary of one adaptive-refinement run.
struct AdaptResult
{
    std::optional<polyMesh> mesh;
    std::vector<AdaptCell> cells;
    label level = 0;          // deepest level actually reached
    label nCells = 0;
    label nSplit = 0;         // number of cells that were split at least once
    scalar volumeSum = 0;     // sum of cell volumes
    scalar integral = 0;      // sum of phi_c * V_c
    scalar integralError = 0; // |integral - exact analytic box integral|
    scalar minVolume = 0;
    scalar maxVolume = 0;
};

/// Adaptive Cartesian refinement operator.
///
/// The base block mesh is a structured hex grid. A cell is marked when the
/// Euclidean norm of its directional gradient exceeds the declared
/// threshold; only the axes along which the field actually varies are split.
/// To keep the rebuilt polyMesh strictly conforming, refinement is applied
/// per base bin along an axis (a slab), so all cells sharing that bin are
/// subdivided together. This is a structured analogue of octree refinement
/// with full face conformity; it is localised around detected gradients for
/// a step field and saturates at `maxDepth`.
class AdaptiveRefine
{
public:
    static AdaptResult refine(const Dictionary& caseDict, const std::string& casePath);
};

} // namespace cf
