#include <cmath>
#include <functional>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/parallel/Backend.h"

namespace cf
{
namespace
{

void gaussGrad(const polyMesh& mesh, const Field<scalar>& data,
               const std::function<scalar(label face, const Field<scalar>&)>& faceValueOf,
               const Field<scalar>& ow, Field<Vector>& gradT, bool corrected)
{
    const label nC = mesh.nCells();
    gradT.assign(static_cast<std::size_t>(nC), Vector::zero());
    const auto& vol = mesh.cellVolumes();
    const auto& sf = mesh.faceSf();

    // pass 0: V * grad = sum Tf * Sf
    Field<Vector> sum(static_cast<std::size_t>(nC), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const Vector flux = sf[u] * faceValueOf(f, data);
        sum[static_cast<std::size_t>(mesh.owner()[u])] += flux;
        if (f < mesh.nInternalFaces())
        {
            sum[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
        }
    }
    // Pure per-cell map (V*grad / volume). Disjoint writes only.
    backend::parallelFor(nC, [&](label c) {
        gradT[static_cast<std::size_t>(c)] =
            sum[static_cast<std::size_t>(c)] / vol[static_cast<std::size_t>(c)];
    });

    if (!corrected) return;

    // correction passes: V*grad = sum [Tf*Sf - Cf*(grad_f . dhat)]
    // with d = faceCentre - ownerCellCentre, Cf = Sf - (Sf.dhat)dhat
    for (int pass = 0; pass < 2; ++pass)
    {
        std::vector<Vector> gradFace(static_cast<std::size_t>(mesh.nFaces()));
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            if (f < mesh.nInternalFaces())
            {
                const Vector gp = gradT[static_cast<std::size_t>(mesh.owner()[u])];
                const Vector gn = gradT[static_cast<std::size_t>(mesh.neighbour()[u])];
                gradFace[u] = ow[u] * gp + (1.0 - ow[u]) * gn;
            }
            else
            {
                gradFace[u] = gradT[static_cast<std::size_t>(mesh.owner()[u])];
            }
        }

        sum.assign(static_cast<std::size_t>(nC), Vector::zero());
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const Vector sfv = sf[u];
            const Vector d = mesh.faceCentres()[u] - mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[u])];
            const scalar dm = mag(d);
            if (dm <= 0) cfFatal("grad: degenerate face-cell distance");
            const Vector dh = d / dm;
            const Vector cfVec = sfv - dh * dot(sfv, dh);
            const Vector flux = sfv * faceValueOf(f, data) - cfVec * dot(gradFace[u], dh);
            sum[static_cast<std::size_t>(mesh.owner()[u])] += flux;
            if (f < mesh.nInternalFaces())
            {
                sum[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
            }
        }
        backend::parallelFor(nC, [&](label c) {
            gradT[static_cast<std::size_t>(c)] =
                sum[static_cast<std::size_t>(c)] / vol[static_cast<std::size_t>(c)];
        });
    }
}

} // namespace

namespace
{

/// Symmetric 3x3 stored as its six upper entries; solved by cofactors with a
/// small ridge, so a sliver stencil without independent directions degrades to
/// a projected gradient instead of a NaN.
struct Sym3
{
    scalar xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;

    void add(const Vector& d, scalar w)
    {
        xx += w * d.x * d.x;
        xy += w * d.x * d.y;
        xz += w * d.x * d.z;
        yy += w * d.y * d.y;
        yz += w * d.y * d.z;
        zz += w * d.z * d.z;
    }

    Vector solve(const Vector& b) const
    {
        const scalar t = std::max((xx + yy + zz) / 3.0, 1e-300);
        const scalar e = 1e-12 * t;
        const scalar XX = xx + e, YY = yy + e, ZZ = zz + e;
        const scalar c11 = YY * ZZ - yz * yz;
        const scalar c12 = xz * yz - xy * ZZ;
        const scalar c13 = xy * yz - YY * xz;
        const scalar c22 = XX * ZZ - xz * xz;
        const scalar c23 = xy * xz - XX * yz;
        const scalar c33 = XX * YY - xy * xy;
        const scalar det = XX * c11 + xy * c12 + xz * c13;
        if (det <= 0 || !std::isfinite(det)) return Vector::zero();
        return {(c11 * b.x + c12 * b.y + c13 * b.z) / det,
                (c12 * b.x + c22 * b.y + c23 * b.z) / det,
                (c13 * b.x + c23 * b.y + c33 * b.z) / det};
    }
};

/// The `leastSquares` gradient (Mimura's cell-face least-squares fit, the same
/// form as OpenFOAM's leastSquares scheme): per cell P the gradient solves
///   min_g  sum_f w_f (delta_f . g - (phi_f - phi_P))^2,   delta_f = Cf - CP,
/// whose normal equations are P g = sum_f w_f delta_f (phi_f - phi_P) with
/// P = sum_f w_f delta_f delta_f^T and the OpenFOAM weight w_f = 1/|delta_f|.
/// On an orthogonal mesh the face values of a linear field are exact and the
/// solution coincides with the Green-Gauss gradient; on a skewed mesh the
/// clamped-linear face values are not, so the scheme builds them from
/// gradients instead of clamped weights:
///   pass 1: bootstrap gradient from NEIGHBOUR CELL values
///     (delta = CN - CP, target phi_N - phi_P -- no interpolation, so a linear
///     field is recovered exactly whatever the skewness);
///   pass 2: face values phi_f = 0.5[(phi_P + g1_P . (Cf-CP))
///                                   + (phi_N + g1_N . (Cf-CN))], exact for a
///     linear field by the pass-1 exactness, then the face-based LS above.
void leastSquaresGrad(const polyMesh& mesh, const Field<scalar>& data,
                      const std::function<scalar(label, const Field<scalar>&)>& faceValueOf,
                      Field<Vector>& gradT)
{
    const label nC = mesh.nCells();
    const auto& cc = mesh.cellCentres();
    const auto& fc = mesh.faceCentres();
    const auto& owner = mesh.owner();
    const auto& nbr = mesh.neighbour();
    const label nInt = mesh.nInternalFaces();
    const std::size_t nU = static_cast<std::size_t>(nC);

    std::vector<Sym3> cellMat(nU), faceMat(nU);
    Field<Vector> cellRhs(nU, Vector::zero());

    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = owner[u];
        const std::size_t ou = static_cast<std::size_t>(o);
        const Vector df = fc[u] - cc[ou]; // owner -> face centre
        const scalar wf = 1.0 / std::max(mag(df), 1e-300);
        faceMat[ou].add(df, wf);
        if (f < nInt)
        {
            const std::size_t nu = static_cast<std::size_t>(nbr[u]);
            const Vector dn = cc[nu] - cc[ou]; // owner -> neighbour centre
            const scalar wn = 1.0 / std::max(mag(dn), 1e-300);
            cellMat[ou].add(dn, wn);
            cellMat[nu].add(dn, wn);
            // both sides of the face see delta = C_other - C_self, so the
            // neighbour's (delta, target) = (-dn, -dPhi/wn) multiply to the
            // SAME signed term w dn (phi_N - phi_P) the owner gets
            const scalar dPhi = (data[nu] - data[ou]) * wn;
            cellRhs[ou] += dn * dPhi;
            cellRhs[nu] += dn * dPhi;
            faceMat[nu].add(fc[u] - cc[nu], 1.0 / std::max(mag(fc[u] - cc[nu]), 1e-300));
        }
        else
        {
            // boundary: the bootstrap stencil reuses the face centre and the
            // BC-linearised face value
            cellMat[ou].add(df, wf);
            const scalar dPhi = (faceValueOf(f, data) - data[ou]) * wf;
            cellRhs[ou] += df * dPhi;
        }
    }

    Field<Vector> g0(nU, Vector::zero());
    backend::parallelFor(nC, [&](label c) {
        const std::size_t u = static_cast<std::size_t>(c);
        g0[u] = cellMat[u].solve(cellRhs[u]);
    });

    Field<Vector> faceRhs(nU, Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = owner[u];
        const std::size_t ou = static_cast<std::size_t>(o);
        scalar vf;
        if (f < nInt)
        {
            const std::size_t nu = static_cast<std::size_t>(nbr[u]);
            vf = 0.5 * ((data[ou] + dot(g0[ou], fc[u] - cc[ou]))
                      + (data[nu] + dot(g0[nu], fc[u] - cc[nu])));
            // the neighbour runs its own face stencil over the same face
            const Vector dN = fc[u] - cc[nu];
            const scalar wN = 1.0 / std::max(mag(dN), 1e-300);
            faceRhs[nu] += dN * ((vf - data[nu]) * wN);
        }
        else
        {
            vf = faceValueOf(f, data);
        }
        const Vector d = fc[u] - cc[ou];
        const scalar w = 1.0 / std::max(mag(d), 1e-300);
        faceRhs[ou] += d * ((vf - data[ou]) * w);
    }

    gradT.assign(nU, Vector::zero());
    backend::parallelFor(nC, [&](label c) {
        const std::size_t u = static_cast<std::size_t>(c);
        gradT[u] = faceMat[u].solve(faceRhs[u]);
    });
}

} // namespace

namespace fvc
{

Field<Vector> grad(const volScalarField& T, const std::string& scheme)
{
    const polyMesh& mesh = T.mesh();
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    Field<Vector> out;
    if (scheme == "leastSquares")
    {
        leastSquaresGrad(mesh, T.internalField(),
                         [&T](label f, const Field<scalar>&) { return T.faceValue(f); }, out);
        return out;
    }
    gaussGrad(mesh, T.internalField(),
              [&T](label f, const Field<scalar>&) { return T.faceValue(f); },
              ow, out, scheme == "corrected");
    return out;
}

Field<Vector> grad(const polyMesh& mesh, const Field<scalar>& data,
                   const std::function<scalar(label face, const Field<scalar>&)>& faceValueOf,
                   const std::string& scheme)
{
    if (scheme == "leastSquares")
    {
        Field<Vector> out;
        leastSquaresGrad(mesh, data, faceValueOf, out);
        return out;
    }
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    Field<Vector> out;
    gaussGrad(mesh, data, faceValueOf, ow, out, scheme == "corrected");
    return out;
}

} // namespace fvc

} // namespace cf
