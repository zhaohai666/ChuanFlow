#pragma once

#include <string>
#include <vector>

#include "cf/common/Field.h"
#include "cf/common/Labels.h"
#include "cf/common/Vector.h"

namespace cf
{

/// Polyhedral unstructured mesh: the geometric foundation of the FVM kernel.
/// Follows the owner/neighbour face-orientation idea from OpenFOAM's polyMesh
/// (faces carry vertex lists; internal faces see owner < neighbour) but is an
/// independent implementation.
class polyMesh
{
public:
    struct Patch
    {
        std::string name;
        std::string type;   // wall | symmetry | inlet | outlet | ...
        label start = 0;    // first face index
        label size = 0;     // number of faces
    };

    polyMesh(Field<Vector> points,
             std::vector<std::vector<label>> faces,
             Field<label> owner,
             Field<label> neighbour,
             std::vector<Patch> patches);

    label nPoints() const { return static_cast<label>(points_.size()); }
    label nCells() const { return nCells_; }
    label nFaces() const { return static_cast<label>(faces_.size()); }
    label nInternalFaces() const { return static_cast<label>(neighbour_.size()); }
    label nBoundaryFaces() const { return nFaces() - nInternalFaces(); }
    label nPatches() const { return static_cast<label>(patches_.size()); }

    const Field<Vector>& points() const { return points_; }
    const std::vector<std::vector<label>>& faces() const { return faces_; }
    const Field<label>& owner() const { return owner_; }
    const Field<label>& neighbour() const { return neighbour_; }
    const std::vector<Patch>& patches() const { return patches_; }

    label patchIdByName(const std::string& name) const;
    std::vector<std::string> patchNames() const;

    // ------------------------------------------------------------------
    // Cached geometry
    // ------------------------------------------------------------------
    const Field<scalar>& cellVolumes() const;
    const Field<Vector>& cellCentres() const;
    const Field<Vector>& faceCentres() const;
    /// Signed area vector, owner -> neighbour (outward for boundary faces).
    const Field<Vector>& faceSf() const;
    const Field<scalar>& faceAreas() const;

    /// Owner-to-neighbour centre vector (boundary faces: owner centre -> face centre).
    const Field<Vector>& deltaCoeffVec() const { return distVec_; }
    /// 1/|distVec| for orthogonal diffusion coefficients.
    const Field<scalar>& deltaCoeff() const { return invDist_; }
    /// Per-cell number of faces / cell-face addressing.
    const std::vector<std::vector<label>>& cells() const { return cells_; }

    /// Face -> patch id for boundary faces, -1 for internal faces.
    const Field<label>& facePatchMap() const { return facePatchMap_; }

    /// Bounding box corner (for reporting).
    Vector minX() const;
    Vector maxX() const;

private:
    void calcGeometry() const;

    Field<Vector> points_;
    std::vector<std::vector<label>> faces_;
    Field<label> owner_;
    Field<label> neighbour_;
    std::vector<Patch> patches_;
    label nCells_ = 0;
    std::vector<std::vector<label>> cells_;
    Field<label> facePatchMap_;

    mutable bool geomCalculated_ = false;
    mutable Field<scalar> cellVolumes_;
    mutable Field<Vector> cellCentres_;
    mutable Field<Vector> faceCentres_;
    mutable Field<Vector> sf_;
    mutable Field<scalar> magSf_;
    mutable Field<Vector> distVec_;
    mutable Field<scalar> invDist_;
};

} // namespace cf
