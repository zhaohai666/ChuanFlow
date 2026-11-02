#pragma once

#include "cf/fvm/FvMatrix.h"

namespace cf
{

/// Explicit face-flux correction consistent with a pressure(-correction)
/// field: `phi_f -= rAUf_f * SnGrad(p)_f + non-orthogonal compensation`.
/// This is OpenFOAM's `pEqn.flux()` step, shared by SIMPLE and PISO so both
/// couplings close the same discrete continuity.
///
/// `rAUf` must be the face-interpolated operator that built the pressure
/// equation, otherwise the corrected flux is not divergence free.
void correctFluxes(const polyMesh& mesh,
                   const volScalarField& p,
                   const Field<Vector>& gradP,
                   const Field<scalar>& rAUf,
                   const Schemes& schemes,
                   surfaceScalarField& phi);

/// phi at internal faces = interp(U) & Sf; boundary faces are untouched.
void updatePhiInternal(const polyMesh& mesh, const volVectorField& U,
                       surfaceScalarField& phi);

/// Interpolate a cell-centred coefficient onto faces with the same bounded
/// owner weights used for the fluxes.
Field<scalar> interpolateToFaces(const polyMesh& mesh, const Field<scalar>& cellData);

/// max |div(phi)| normalised by the largest face area.
scalar continuityError(const polyMesh& mesh, const surfaceScalarField& phi);

} // namespace cf
