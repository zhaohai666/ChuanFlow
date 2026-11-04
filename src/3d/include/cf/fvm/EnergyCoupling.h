#pragma once

#include "cf/common/Field.h"
#include "cf/common/Vector.h"

namespace cf
{

/// The part of a thermal model the pressure-velocity solvers use, declared in
/// the FVM layer so the solvers stay independent of the energy equation - the
/// same role EddyViscosity.h plays for the turbulence models: the buoyancy
/// acceleration the momentum equation is forced with, and a hook to advance the
/// temperature once per outer iteration.
class IEnergyCoupling
{
public:
    virtual ~IEnergyCoupling() = default;

    /// Transport the temperature with the flux and velocity the solver just
    /// produced. `dt` <= 0 solves the steady equation, under-relaxed against
    /// the field's previous state.
    virtual void correctEnergy(scalar dt) = 0;

    /// Boussinesq acceleration (m/s2) per cell from the last solved
    /// temperature, null when the case has no buoyancy.
    virtual const Field<Vector>* buoyancyForce() const = 0;
};

} // namespace cf
