#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{
namespace
{

/// Face-transmissibility coefficient: (Sf . d)/|d|^2, d = owner->neighbour
/// (boundary: owner->face centre). Equals |Sf|/|d| on orthogonal meshes.
scalar orthoCoeff(const polyMesh& mesh, label face)
{
    const std::size_t u = static_cast<std::size_t>(face);
    const Vector& d = mesh.deltaCoeffVec()[u];
    const scalar d2 = dot(d, d);
    if (d2 <= 0) cfFatal("laplacian: zero face distance");
    return dot(mesh.faceSf()[u], d) / d2;
}

std::vector<scalar> faceOrthoCoeffs(const polyMesh& mesh)
{
    std::vector<scalar> c(static_cast<std::size_t>(mesh.nFaces()));
    for (label f = 0; f < mesh.nFaces(); ++f) c[static_cast<std::size_t>(f)] = orthoCoeff(mesh, f);
    return c;
}

/// Cell-face gamma for a possibly-uniform diffusion coefficient.
scalar cellGammaAt(const Field<scalar>& gammaCell, label cell)
{
    return gammaCell.size() == 1 ? gammaCell[0] : gammaCell[static_cast<std::size_t>(cell)];
}

} // namespace

namespace
{

void assembleLaplacian(fvMatrix& eqn, std::vector<scalar> gammaFace, const std::string& scheme, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    const FieldView& T = eqn.unknown();
    const bool corrected = scheme == "corrected";
    const std::vector<scalar> oc = faceOrthoCoeffs(mesh);
    Field<scalar> ow = makeInternalOwnerWeights(mesh);

    // static implicit part
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar coef = gammaFace[u] * oc[u];
        const label o = mesh.owner()[u];

        if (f < mesh.nInternalFaces())
        {
            const label n = mesh.neighbour()[u];
            for (int comp = 0; comp < T.nComp; ++comp)
            {
                eqn.addDiag(o, comp, -sign * coef);
                eqn.addOffDiag(o, n, comp, sign * coef);
                eqn.addDiag(n, comp, -sign * coef);
                eqn.addOffDiag(n, o, comp, sign * coef);
            }
        }
        else
        {
            for (int comp = 0; comp < T.nComp; ++comp)
            {
                scalar c = 0, k = 0;
                T.lin(f, comp, c, k);
                // flux = coef * (Tf - Tp) = coef*((c-1)Tp + k)
                eqn.addDiag(o, comp, sign * coef * (c - 1.0));
                eqn.addRhs(o, comp, -sign * coef * k);
            }
        }
    }

    if (!corrected) return;

    // non-orthogonal correction: corr_f = gamma * (Sf - (Sf.dhat)dhat) . gradTf
    const polyMesh* msh = &mesh;
    const FieldView view = T;
    eqn.registerCorrector(
        [msh, view, gf = std::move(gammaFace), ow = std::move(ow), sign](
            fvMatrix& eq, const std::vector<Field<scalar>>& x)
        {
            for (int comp = 0; comp < view.nComp; ++comp)
            {
                auto faceVal = makeFaceValueFn(*msh, view, comp);
                Field<Vector> gCell = fvc::grad(
                    *msh, x[static_cast<std::size_t>(comp)], faceVal, "linear");

                for (label f = 0; f < msh->nFaces(); ++f)
                {
                    const std::size_t u = static_cast<std::size_t>(f);
                    Vector d;
                    Vector gradFace;
                    if (f < msh->nInternalFaces())
                    {
                        const Vector gp = gCell[static_cast<std::size_t>(msh->owner()[u])];
                        const Vector gn = gCell[static_cast<std::size_t>(msh->neighbour()[u])];
                        gradFace = ow[u] * gp + (1.0 - ow[u]) * gn;
                        d = msh->deltaCoeffVec()[u];
                    }
                    else
                    {
                        gradFace = gCell[static_cast<std::size_t>(msh->owner()[u])];
                        d = msh->deltaCoeffVec()[u];
                    }
                    const scalar dm = mag(d);
                    if (dm <= 0) continue;
                    const Vector dh = d / dm;
                    const Vector sfv = msh->faceSf()[u];
                    const Vector cfVec = sfv - dh * dot(sfv, dh);
                    const scalar corr = gf[u] * dot(cfVec, gradFace);
                    eq.addCorrection(msh->owner()[u], comp, -sign * corr);
                    if (f < msh->nInternalFaces())
                    {
                        eq.addCorrection(msh->neighbour()[u], comp, sign * corr);
                    }
                }
            }
        });
}

} // namespace

namespace fvm
{

void laplacian(fvMatrix& eqn, scalar gammaUniform, const std::string& scheme, int sign)
{
    const Field<scalar> g(1, gammaUniform);
    laplacian(eqn, g, scheme, sign);
}

void laplacian(fvMatrix& eqn, const Field<scalar>& gammaCell, const std::string& scheme, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    std::vector<scalar> gammaFace(static_cast<std::size_t>(mesh.nFaces()));
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        if (f < mesh.nInternalFaces())
        {
            const scalar gp = cellGammaAt(gammaCell, mesh.owner()[u]);
            const scalar gn = cellGammaAt(gammaCell, mesh.neighbour()[u]);
            const Vector& cP = mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[u])];
            const Vector& cN = mesh.cellCentres()[static_cast<std::size_t>(mesh.neighbour()[u])];
            const Vector& cf = mesh.faceCentres()[u];
            const scalar dP = mag(cf - cP), dN = mag(cN - cf);
            const scalar tot = dP + dN;
            gammaFace[u] = (tot > 0 && gp > 0 && gn > 0)
                ? tot * gp * gn / (dN * gp + dP * gn)
                : 0.5 * (gp + gn);
        }
        else
        {
            gammaFace[u] = cellGammaAt(gammaCell, mesh.owner()[u]);
        }
    }
    assembleLaplacian(eqn, std::move(gammaFace), scheme, sign);
}

void laplacianFace(fvMatrix& eqn, const Field<scalar>& gammaFace, const std::string& scheme, int sign)
{
    if (static_cast<label>(gammaFace.size()) != eqn.mesh().nFaces())
    {
        cfFatal("fvm::laplacianFace: gamma size != nFaces");
    }
    assembleLaplacian(eqn, std::vector<scalar>(gammaFace.begin(), gammaFace.end()), scheme, sign);
}

} // namespace fvm
} // namespace cf
