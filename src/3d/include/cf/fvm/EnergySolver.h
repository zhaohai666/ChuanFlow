#pragma once

#include <string>

#include "cf/common/Dictionary.h"
#include "cf/fvm/EnergyCoupling.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{

/// Settings of the incompressible energy equation, read from the case's
/// `energy{}` block:
///
///   dT/dt + div(u, T) - div(alpha grad T) = S
///   alpha = kappa/(rho cp) + nut/turbulentPrandtl      [m2/s]
///   S     = volumetricSource/(rho cp) + nu |S_ij|^2/cp  [K/s]
///
/// Everything the equation itself needs is in these units; the solvers stay
/// kinematic (density 1), so `rho` only scales the conductivity and the source,
/// and the buoyancy enters as an acceleration.
struct EnergyProps
{
    scalar rho = 1000;             // kg/m3
    scalar cp = 4180;              // J/(kg K)
    scalar kappa = 0.6;            // W/(m K) molecular conductivity
    scalar turbulentPrandtl = 0.9; // Pr_t of the eddy diffusivity
    scalar volumetricSource = 0;   // W/m3, uniform heat generation
    bool viscousDissipation = false;
    scalar beta = 0;               // 1/K thermal expansion
    scalar tRef = 0;               // K, temperature the buoyancy force vanishes at
    Vector gravity = Vector::zero(); // m/s2, e.g. (0 0 -9.81)
    scalar relaxation = 0.8;       // steady under-relaxation of T
    scalar initialT = 0;           // K, uniform starting temperature
    int initialMode = 0;           // cosine mode of the starting perturbation
    scalar initialAmplitude = 0;   // K
    std::string initialAxis = "z"; // axis the mode is measured along

    /// Read the case's `energy{}` block, then validate(). The apps check
    /// `caseDict.has("energy")` first: no block means the case has no thermal
    /// field at all.
    static EnergyProps fromCase(const Dictionary& caseDict);

    /// Refuse settings the equation cannot integrate: a non-finite property, a
    /// non-positive rho, cp or kappa (which would freeze the diffusivity or the
    /// source), a non-positive turbulent Prandtl, a relaxation factor outside
    /// (0, 1] that either discards the equation or amplifies every iterate, and
    /// a negative starting mode.
    void validate() const;

    /// kappa/(rho cp): the molecular thermal diffusivity [m2/s].
    scalar alphaMolecular() const { return kappa / (rho * cp); }
};

/// Transport of temperature through the incompressible solvers' own flux, and
/// the Boussinesq acceleration the momentum equation is forced with.
///
/// The equation is assembled with the same fvm terms as every other transport
/// here (fvm::ddt, fvm::div over the continuity flux, fvm::laplacian over the
/// face-harmonic diffusivity), so the temperature sees exactly the flux that
/// closes continuity. Convection is upwind by that flux, which is what makes
/// the discrete total enthalpy balance hold to the linear-solver tolerance:
/// heat in through the boundaries plus heat generated equals enthalpy out.
///
/// Thermal boundary conditions are the scalar-field patch conditions:
///   fixedValue     - prescribed temperature
///   fixedGradient  - prescribed gradient along the outward face normal; the
///                    case files it as `heatFlux` (W/m2 INTO the fluid, so
///                    positive heats) and the app converts it with Fourier's
///                    law, grad = heatFlux/kappa
///   zeroGradient   - adiabatic
///
/// Buoyancy is the Boussinesq deviation only: `f = -beta (T - tRef) gravity`.
/// The reference hydrostatic (rho g) is not part of this force, so a case that
/// needs it carries it in the solver's own `force`. A uniform force is the one
/// the pressure of a projection solve can absorb - `tests/cases/energyBuoyantRest`
/// measures what "absorb" means here: in a box closed on every face the coupling
/// leaves a residual |U| of ~2.7e-4 m/s against the 2.45 m/s the same force
/// would reach unbalanced, and a case with no thermal field driven by an
/// equivalent `force` reproduces the same number, so the floor belongs to the
/// closed-box pressure-velocity coupling, not to the energy equation. A
/// stratified column whose hydrostatics vary with height needs a p_rgh
/// formulation (pending) before it can hold a conduction state without a
/// spurious circulation.
///
/// The eddy diffusivity is the plain gradient-transport analogy with Pr_t. On a
/// mesh that resolves the sublayer (y+ ~ 1) nut -> 0 at the wall and the wall
/// heat flux stays molecular, which is what this model assumes; on a
/// wall-function mesh it over-predicts the transfer, and a thermal wall function
/// is the next increment.
class EnergyModel : public IEnergyCoupling
{
public:
    /// `phi` is the solver's continuity flux and is read on every correct().
    /// `nut` and `strainRate` are optional turbulence state, read as it stands
    /// when the equation is assembled (null for laminar). The solver owns the
    /// mesh, flux and turbulence model, so this object must not outlive them.
    EnergyModel(const polyMesh& mesh, const surfaceScalarField& phi,
                const Schemes& schemes, const EnergyProps& props, scalar nu,
                const Field<scalar>* nut = nullptr,
                const Field<scalar>* strainRate = nullptr);

    const EnergyProps& props() const { return props_; }

    volScalarField& temperature() { return t_; }
    const volScalarField& temperature() const { return t_; }

    /// Starting field: `initialT` everywhere, plus
    /// `initialAmplitude * cos(mode * pi * xi)` along `initialAxis`, xi the
    /// normalised coordinate over the mesh bounding box. With zeroGradient walls
    /// such a cosine is an exact eigenvector of the discrete diffusion operator,
    /// which is what makes a decay case analytically checkable.
    void initialize();

    void correctEnergy(scalar dt) override;
    const Field<Vector>* buoyancyForce() const override;

    /// Thermal diffusivity [m2/s] per cell: molecular plus nut/Pr_t.
    Field<scalar> alphaCell() const;

    /// Heat rate (W) entering the fluid through a patch, from the same face
    /// diffusion flux the equation is assembled with. Positive heats the
    /// fluid, so over all patches it balances the volumetric source and the
    /// net enthalpy flux.
    scalar patchHeatRate(label patch) const;
    /// integral of the volumetric source over the domain (W).
    scalar volumetricHeatRate() const;
    /// rho cp times the net volumetric flux leaving through the boundary faces,
    /// each weighted by its own face temperature (W): the enthalpy the fluid
    /// carries out. With the patch heat rates it closes the discrete balance.
    scalar netEnthalpyRate() const;

  private:
    /// Sources [K/s] per cell: volumetric generation and viscous dissipation.
    Field<scalar> sourceRate() const;
    /// f = -beta (T - tRef) gravity over the current temperature.
    void refreshBuoyancy();

    const polyMesh& mesh_;
    const surfaceScalarField& phi_;
    Schemes schemes_;
    EnergyProps props_;
    scalar nu_;
    const Field<scalar>* nut_;
    const Field<scalar>* strainRate_;
    volScalarField t_;
    Field<Vector> buoyancy_;
    bool hasBuoyancy_ = false;

    static constexpr scalar linearTol = 1e-13;
    static constexpr int linearMaxIter = 20000;
};

} // namespace cf
