#pragma once

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PressureCorrection.h"

namespace cf
{
namespace app
{

/// Optional `initial { U (Ux Uy Uz); }` block: the cell-centred starting
/// velocity. Without it the domain starts at rest, which is what a startup
/// case wants; an inlet-driven turbulent case needs a started field, because
/// the first pressure correction would otherwise have to build the whole flux
/// of the domain from nothing and overshoots.
inline void seedInitialVelocity(const Dictionary& caseDict, const polyMesh& mesh,
                                volVectorField& U, surfaceScalarField& phi)
{
    if (!caseDict.has("initial")) return;
    const Dictionary& d = caseDict.getDict("initial");
    if (!d.has("U")) return;

    const auto v = d.getScalarList("U");
    const Vector u0{v[0], v[1], v[2]};
    for (Vector& u : U.internalField()) u = u0;
    updatePhiInternal(mesh, U, phi);
    phi.correctBoundaryFlux(U);
    log::info("initial U = (", u0.x, ", ", u0.y, ", ", u0.z, ")");
}

} // namespace app
} // namespace cf
