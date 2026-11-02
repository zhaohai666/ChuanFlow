#pragma once

#include <functional>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Matrix.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Read interface to the unknown field of an equation, so schemes can
/// treat scalar and vector matrices with one code path.
struct FieldView
{
    int nComp = 1;
    std::function<scalar(label cell, int comp)> cell;
    /// boundary linearisation at face (must be a boundary face):
    /// faceValue = coeff * cellValue(owner, comp) + konst
    std::function<void(label face, int comp, scalar& coeff, scalar& konst)> lin;
};

FieldView viewOf(volScalarField& f);
FieldView viewOf(volVectorField& f);

/// Schemes configuration, mirrors the role of OpenFOAM's fvSchemes file.
struct Schemes
{
    std::string ddt = "steady";       // steady | Euler
    std::string div = "upwind";       // upwind | linear | blend
    std::string laplacian = "corrected"; // none | corrected
    std::string grad = "corrected";   // linear | corrected
    scalar blendCoeff = 0.9;          // for div=blend (fraction of central)
    scalar relaxation = 1.0;          // under-relaxation applied after solve
};

/// Multi-component implicit matrix: A x = b assembled face by face.
/// The role of OpenFOAM's fvMatrix/lduMatrix, on a COO->CSR backend.
class fvMatrix
{
public:
    fvMatrix(const polyMesh& mesh, FieldView unknown, std::string name);

    const polyMesh& mesh() const { return mesh_; }
    const FieldView& unknown() const { return unknown_; }
    const std::string& name() const { return name_; }
    const Schemes& schemes() const { return schemes_; }
    Schemes& schemes() { return schemes_; }

    void addDiag(label cell, int comp, scalar v);
    void addRhs(label cell, int comp, scalar v);
    void addOffDiag(label row, label col, int comp, scalar v);

    /// Terms that break symmetry (e.g. convection) switch to BiCGStab.
    void markNonSymmetric() { preferCG_ = false; }

    /// Explicit contributions that are recomputed from the current solution
    /// between non-orthogonal correction iterations (e.g. snGrad correction).
    /// Correctors add into the rhs only; the static rhs is restored first.
    void registerCorrector(std::function<void(fvMatrix& m, const std::vector<Field<scalar>>& solution)> fn);

    /// Add to rhs without recording in the static base (used by correctors).
    void addCorrection(label cell, int comp, scalar v);

    struct Controls
    {
        scalar tol = 1e-7;
        int maxIter = 2000;
        int nonOrth = 2;
    };

    /// Assemble and solve; writes the result into `solution` (per comp) and
    /// returns the report of the last linear solve.
    SolverReport solve(const Controls&);

    std::vector<Field<scalar>>& solution() { return x_; }
    const std::vector<Field<scalar>>& solution() const { return x_; }

    /// Matrix diagonal (needed for SIMPLE: rAU = 1/A()).
    const Field<scalar>& aDiag(int comp) const { return diag_[static_cast<std::size_t>(comp)]; }

    /// A x_off + rhs excluding the diagonal contribution: H(U) for SIMPLE.
    void explicitOperator(const Field<scalar>& xIn, int comp, Field<scalar>& h) const;

private:
    const polyMesh& mesh_;
    FieldView unknown_;
    std::string name_;
    Schemes schemes_;
    std::vector<SparseMatrix::Builder> builders_;
    std::vector<Field<scalar>> diag_;
    std::vector<Field<scalar>> rhs_;
    std::vector<SparseMatrix> matrices_;
    std::vector<Field<scalar>> x_;
    std::vector<Field<scalar>> baseRhs_;
    bool preferCG_ = true;
    std::vector<std::function<void(fvMatrix&, const std::vector<Field<scalar>>&)>> correctors_;
};

// ------------------------------------------------------------------
// fvc: explicit (evaluated on the current field) operators
// ------------------------------------------------------------------

namespace fvc
{

/// Linear interpolation weights at internal faces: wOwner such that
/// Tf = wOwner*TP + (1-wOwner)*TN. Boundary faces use the BC directly.
void interpolate(const polyMesh& mesh, const Field<scalar>& cellData,
                 const std::function<scalar(label face)>& faceValueOf,
                 Field<scalar>& faceData);

Field<scalar> interpolate(const volScalarField& T);

/// Gauss gradient; scheme "linear" or "corrected" (two correction passes).
Field<Vector> grad(const volScalarField& T, const std::string& scheme);
Field<Vector> grad(const polyMesh& mesh, const Field<scalar>& data,
                   const std::function<scalar(label face, const Field<scalar>&)>& faceValueOf,
                   const std::string& scheme);

/// Explicit divergence of a flux field: sum over faces of phi * Tf.
Field<scalar> div(const volScalarField& T, const Field<scalar>& phi);

/// Net face flux per cell (continuity residual): sum over faces of +/- phi.
Field<scalar> divPhi(const polyMesh& mesh, const Field<scalar>& phi);

} // namespace fvc

/// Per-face bounded linear weight of the owner cell (internal faces only).
Field<scalar> makeInternalOwnerWeights(const polyMesh& mesh);

/// Face-value evaluator for an equation unknown: internal faces interpolate
/// the passed cell data, boundary faces apply the BC linearisation.
std::function<scalar(label face, const Field<scalar>&)> makeFaceValueFn(
    const polyMesh& mesh, const FieldView& view, int comp);

// ------------------------------------------------------------------
// fvm: implicit terms added into an fvMatrix. `sign` is the sign with
// which the term appears in the equation (all terms live on the lhs).
// ------------------------------------------------------------------

namespace fvm
{

void ddt(fvMatrix& eqn, const Field<scalar>& oldField, scalar dt, int sign = 1);
void ddt(fvMatrix& eqn, const Field<Vector>& oldField, scalar dt, int sign = 1);

/// scheme: upwind | linear | blend (uses eqn.schemes().blendCoeff)
void div(fvMatrix& eqn, const Field<scalar>& phi, const std::string& scheme, int sign = 1);

/// Diffusion term grad·(gamma grad). gamma is a face-independent scalar
/// per cell; harmonic face averaging is used when gammaCell varies.
void laplacian(fvMatrix& eqn, scalar gammaUniform, const std::string& scheme, int sign = 1);
void laplacian(fvMatrix& eqn, const Field<scalar>& gammaCell, const std::string& scheme, int sign = 1);

/// Face-centred gamma (e.g. SIMPLE's interpolated rAUf): used verbatim,
/// no cell-to-face interpolation. gammaFace must have mesh.nFaces() entries.
void laplacianFace(fvMatrix& eqn, const Field<scalar>& gammaFace, const std::string& scheme, int sign = 1);

/// Implicit source linearisation: sign * Sp * T.
void sp(fvMatrix& eqn, const Field<scalar>& Sp, int sign = 1);

} // namespace fvm

/// Load scheme settings from a case dictionary block:
///   schemes { ddt steady; div upwind; laplacian corrected; grad corrected; }
Schemes loadSchemes(const Dictionary& caseDict);

} // namespace cf
