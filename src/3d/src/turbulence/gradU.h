#pragma once

#include <cmath>

#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{
namespace turb_detail
{

/// gradU[c][i] = gradient of the c-th velocity component at cell i, so
/// gradU[c][i][s] = d U_c / d x_s. Shared by the LES and RSM closures, which
/// both need the full tensor rather than the strain contraction the RANS base
/// evaluates.
inline void computeGradU(const polyMesh& mesh, volVectorField& U,
                         const Schemes& schemes, Field<Vector> gradU[3])
{
    const FieldView view = viewOf(U);
    for (int comp = 0; comp < 3; ++comp)
    {
        Field<scalar> g(static_cast<std::size_t>(mesh.nCells()));
        for (std::size_t i = 0; i < g.size(); ++i) g[i] = U.internalField()[i][comp];
        gradU[comp] = fvc::grad(mesh, g, makeFaceValueFn(mesh, view, comp), schemes.grad);
    }
}

/// Strain-rate tensor S_ij = (dU_i/dx_j + dU_j/dx_i)/2 and |S| = sqrt(2 S:S)
/// at one cell, from the component gradients.
inline void strainAt(const Field<Vector> gradU[3], std::size_t c,
                     scalar S[3][3], scalar& magS)
{
    scalar s2 = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            S[i][j] = 0.5 * (gradU[i][c][j] + gradU[j][c][i]);
            s2 += S[i][j] * S[i][j];
        }
    }
    magS = std::sqrt(2.0 * s2);
}

} // namespace turb_detail
} // namespace cf
