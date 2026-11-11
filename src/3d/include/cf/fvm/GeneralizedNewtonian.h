#pragma once

#include <vector>

#include "cf/common/Field.h"
#include "cf/common/Labels.h"
#include "cf/common/Rheology.h"
#include "cf/fvm/EddyViscosity.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Numerical controls of the algebraic generalized-Newtonian law. These are
/// not fitted constants: `shearFloor` only keeps gdot -> 0 finite,
/// `regularization` is Papanastasiou's m in tau_0(1 - exp(-m gdot)), and
/// `relaxation` damps the viscosity field between SIMPLE sweeps the same way
/// the turbulence models damp nut. A case with tau_0 = 0 and n = 1 has no
/// need for any of them, since its viscosity is the constant it always was.
struct GeneralizedNewtonianControls
{
    scalar shearFloor = 1e-6;       // 1/s
    scalar regularization = 100.0;  // s
    scalar relaxation = 0.7;        // 0 < r <= 1
};

/// The effective viscosity the pressure-velocity solvers assemble the momentum
/// diffusion with, for a power-law / Herschel-Bulkley fluid:
///
///     nu_eff(gdot) = [ K gdot^(n-1) + tau_0 (1-exp(-m gdot))/gdot ] / rho
///
/// It implements `IEddyViscosity` because SIMPLE and PISO already own a
/// per-face viscosity hook: `fvm::laplacianFace` takes exactly the face field
/// this class refreshes once per outer iteration. Nothing in either solver's
/// equation assembly changes, and a Newtonian declaration produces a face
/// field that is bit-identical to the constant `nu` the solvers used before.
///
/// Cell values come from the velocity gradient's symmetric part,
/// gdot = |S| = sqrt(2 S:S), the same strain measure the RANS models use.
/// Wall faces are overwritten with the shear rate the discrete wall flux is
/// built from, |U_c,parallel| / y_wall, so the wall shear stress the momentum
/// equation charges is tau_w = mu_eff(gdot_w) * gdot_w rather than a
/// cell-centre viscosity extrapolated to the wall. That distinction matters
/// for a shear-thinning fluid: the wall is where gdot is largest and nu_eff
/// smallest.
///
/// This is an algebraic law, so `correctTurbulence(dt)` ignores dt: it is
/// called once per sweep by the solver to refresh the field, not to advance a
/// transport equation.
class GeneralizedNewtonianViscosity : public IEddyViscosity
{
public:
    GeneralizedNewtonianViscosity(const polyMesh& mesh, volVectorField& U,
                                  scalar rho, scalar nuNewtonian,
                                  const Rheology& rheology,
                                  const Schemes& schemes,
                                  const GeneralizedNewtonianControls& ctl = {});

    void correctTurbulence(scalar dt) override;
    const Field<scalar>& nuEffFace() const override { return nuEffFace_; }

    const Rheology& rheology() const { return rheology_; }
    const Field<scalar>& strainRate() const { return strain_; }
    const Field<scalar>& viscosityCell() const { return nuCell_; }
    const std::vector<scalar>& wallViscosity() const { return wallNu_; }
    const std::vector<scalar>& wallShearRate() const { return wallGdot_; }
    const std::vector<label>& wallFaces() const { return wallFaces_; }

    scalar maxShearRate() const;
    scalar minViscosity() const;
    scalar maxViscosity() const;
    /// Face-average of the discrete wall shear stress rho*nu_w*gdot_w (Pa).
    /// The developed pressure gradient of a plane channel answers to exactly
    /// this number: |dp/dx| = 2 tau_w / H.
    scalar meanWallShearPa() const;
    /// Area averages of the two factors of that stress, so a regression can
    /// rebuild tau_w = rho * nu_w * gdot_w from the printed numbers and check
    /// nu_w against K*gdot_w^(n-1)/rho rather than against a fitted curve.
    scalar meanWallViscosity() const;
    scalar meanWallShearRate() const;

private:
    void computeStrainRate();
    scalar apparentNu(scalar gdot) const;

    const polyMesh& mesh_;
    volVectorField& U_;
    scalar rho_;
    scalar nuNewtonian_;
    scalar muNewtonian_;
    Rheology rheology_;
    Schemes schemes_;
    GeneralizedNewtonianControls ctl_;

    Field<scalar> strain_;
    Field<scalar> nuCell_;
    Field<scalar> nuEffFace_;
    std::vector<label> wallFaces_;
    std::vector<scalar> wallDist_;
    std::vector<scalar> wallNu_;
    std::vector<scalar> wallGdot_;
};

} // namespace cf
