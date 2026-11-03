#pragma once

#include <string>
#include <vector>

#include "cf/common/Field.h"
#include "cf/common/Labels.h"
#include "cf/common/Vector.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Boundary condition kinds carried by a patch (roles of OpenFOAM's
/// fvPatchField types, implemented per-field with a compact table):
///   fixedValue     - face value prescribed
///   zeroGradient   - face value equals owner cell value
///   fixedGradient  - prescribed gradient along the face normal:
///                    Tf = Tp + grad * normalDistance
///   symmetry       - scalar: zeroGradient; vector: zero normal velocity,
///                    zeroGradient tangential (axis-aligned normals)
///   empty          - face excluded from surface sums (2D planes)
inline bool bcIsSurfaceEmpty(const std::string& type) { return type == "empty"; }

// ------------------------------------------------------------------
// vol fields (cell centred) with per-patch boundary conditions
// ------------------------------------------------------------------

class volScalarField
{
public:
    volScalarField(const polyMesh& mesh, const std::string& name, scalar init = 0);

    const std::string& name() const { return name_; }
    const polyMesh& mesh() const { return mesh_; }

    Field<scalar>& internalField() { return internal_; }
    const Field<scalar>& internalField() const { return internal_; }

    label patchIdOfFace(label face) const; // -1 for internal faces

    void setPatchValue(label patchId, const std::string& bcType, scalar value);
    std::string patchBcType(label patchId) const;

    /// Per-face value for a boundary face, overriding the patch-wide value of
    /// setPatchValue for that face (the type still comes from the patch). Used
    /// where the boundary value varies along a patch, e.g. the turbulence
    /// dissipation implied by the law of the wall. Ignored for faces that are
    /// internal or whose patch type does not carry a value.
    void setFaceValue(label face, scalar value);
    void clearFaceValues();

    /// face value from the current field (explicit use)
    scalar faceValue(label face) const;
    /// implicit linearisation: value = coeff * T[owner] + konst
    void linearise(label face, scalar& coeff, scalar& konst) const;

    scalar min() const;
    scalar max() const;

private:
    struct BC
    {
        std::string type = "zeroGradient";
        scalar value = 0;
    };

    const polyMesh& mesh_;
    std::string name_;
    Field<scalar> internal_;
    std::vector<BC> bc_;
    Field<scalar> faceValueOverride_;
    Field<char> faceValueSet_;
};

class volVectorField
{
public:
    volVectorField(const polyMesh& mesh, const std::string& name, const Vector& init = Vector::zero());

    const std::string& name() const { return name_; }
    const polyMesh& mesh() const { return mesh_; }

    Field<Vector>& internalField() { return internal_; }
    const Field<Vector>& internalField() const { return internal_; }

    label patchIdOfFace(label face) const;

    void setPatchValue(label patchId, const std::string& bcType, const Vector& value);
    std::string patchBcType(label patchId) const;
    /// prescribed value of the patch condition (zero for gradient types)
    Vector patchValue(label patchId) const;

    Vector faceValue(label face) const;
    /// component-wise linearisation against owner cell value of that component
    void linearise(label face, int comp, scalar& coeff, scalar& konst) const;

private:
    struct BC
    {
        std::string type = "zeroGradient";
        Vector value = Vector::zero();
    };

    const polyMesh& mesh_;
    std::string name_;
    Field<Vector> internal_;
    std::vector<BC> bc_;
};

/// Face-centred flux field: phi_f = U & Sf, owner -> neighbour orientation.
class surfaceScalarField
{
public:
    surfaceScalarField(const polyMesh& mesh, const std::string& name, scalar init = 0);

    const std::string& name() const { return name_; }
    const polyMesh& mesh() const { return mesh_; }

    Field<scalar>& values() { return values_; }
    const Field<scalar>& values() const { return values_; }

    /// Recompute boundary phi from a velocity field: phi_b = Uf & Sf.
    /// Patches with zero-normal (wall/symmetry) get phi = 0 automatically.
    void correctBoundaryFlux(const volVectorField& U);
    /// Set boundary phi from a uniform velocity vector (inlet-type patches);
    /// walls/symmetry stay zero.
    void fixBoundaryVelocity(const Vector& U0);

private:
    const polyMesh& mesh_;
    std::string name_;
    Field<scalar> values_;
};

} // namespace cf
