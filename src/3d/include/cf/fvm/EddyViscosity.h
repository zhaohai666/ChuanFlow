#pragma once

#include "cf/common/Field.h"
#include "cf/common/Labels.h"

namespace cf
{

/// The part of a turbulence model the pressure-velocity solvers use, declared
/// in the FVM layer so the solvers stay independent of the turbulence models:
/// the effective viscosity the momentum equation is assembled with, and a
/// hook to advance the turbulence equations once per outer iteration.
class IEddyViscosity
{
public:
    virtual ~IEddyViscosity() = default;

    /// nu + nut at every face (mesh.nFaces() entries). The boundary faces of
    /// wall patches carry the law-of-the-wall value, so the wall shear of the
    /// momentum equation is the one the wall layer prescribes instead of the
    /// one the near-wall cell viscosity would give.
    virtual const Field<scalar>& nuEffFace() const = 0;

    /// Advance the turbulence equations with the current velocity and flux.
    /// dt <= 0 solves the steady equations, under-relaxed.
    virtual void correctTurbulence(scalar dt) = 0;
};

} // namespace cf
