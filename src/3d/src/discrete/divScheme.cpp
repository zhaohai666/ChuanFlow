#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{
namespace fvm
{

void div(fvMatrix& eqn, const Field<scalar>& phi, const std::string& scheme, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    const FieldView& T = eqn.unknown();
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    eqn.markNonSymmetric();

    const scalar alpha = scheme == "blend" ? eqn.schemes().blendCoeff : (scheme == "linear" ? 1.0 : 0.0);

    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label o = mesh.owner()[u];
        const scalar p = phi[u];

        if (f < mesh.nInternalFaces())
        {
            const label n = mesh.neighbour()[u];
            // upwind weights
            scalar wPo = 0.0, wNo = 0.0; // weights on T_owner / T_neighbour
            if (p >= 0) { wPo = 1.0; } else { wNo = 1.0; }
            if (alpha > 0)
            {
                // blend central + upwind (alpha = 1 pure central)
                const scalar c = ow[u];
                wPo = alpha * c + (1.0 - alpha) * wPo;
                wNo = alpha * (1.0 - c) + (1.0 - alpha) * wNo;
            }
            // owner equation: +flux; neighbour equation: -flux
            //   flux = p*(wPo*TP + wNo*TN)
            const scalar fp = sign * p * wPo; // coefficient on owner value
            const scalar fn = sign * p * wNo; // coefficient on neighbour value
            for (int comp = 0; comp < T.nComp; ++comp)
            {
                eqn.addDiag(o, comp, fp);
                eqn.addOffDiag(o, n, comp, fn);
                eqn.addDiag(n, comp, -fn);
                eqn.addOffDiag(n, o, comp, -fp);
            }
        }
        else
        {
            for (int comp = 0; comp < T.nComp; ++comp)
            {
                scalar c = 0, k = 0;
                T.lin(f, comp, c, k);
                eqn.addDiag(o, comp, sign * p * c);
                eqn.addRhs(o, comp, -sign * p * k);
            }
        }
    }
}

void ddt(fvMatrix& eqn, const Field<scalar>& oldField, scalar dt, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    if (dt <= 0) cfFatal("fvm::ddt: non-positive time step");
    if (eqn.unknown().nComp != 1) cfFatal("fvm::ddt: scalar old field on vector equation");
    const auto& vol = mesh.cellVolumes();
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar a = sign * vol[u] / dt;
        eqn.addDiag(c, 0, a);
        eqn.addRhs(c, 0, a * oldField[u]);
    }
}

void ddt(fvMatrix& eqn, const Field<Vector>& oldField, scalar dt, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    if (dt <= 0) cfFatal("fvm::ddt: non-positive time step");
    const auto& vol = mesh.cellVolumes();
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar a = sign * vol[u] / dt;
        for (int comp = 0; comp < eqn.unknown().nComp; ++comp)
        {
            eqn.addDiag(c, comp, a);
            eqn.addRhs(c, comp, a * oldField[u][comp]);
        }
    }
}

void sp(fvMatrix& eqn, const Field<scalar>& Sp, int sign)
{
    const polyMesh& mesh = eqn.mesh();
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        for (int comp = 0; comp < eqn.unknown().nComp; ++comp)
        {
            eqn.addDiag(c, comp, sign * Sp[static_cast<std::size_t>(c)]);
        }
    }
}

} // namespace fvm
} // namespace cf
