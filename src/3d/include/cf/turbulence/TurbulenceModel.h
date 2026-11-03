#pragma once

#include <memory>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Field.h"
#include "cf/fvm/EddyViscosity.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/turbulence/WallFunction.h"

namespace cf
{

/// Numerical floors of the transport equations. k and the dissipation variable
/// are only ever produced in pairs, and both linearisations divide one by the
/// other, so both need a floor well below any physical value of the cases here
/// (1e-12 m2/s3 for epsilon corresponds to a friction velocity of 1e-4 m/s).
namespace turb
{
inline constexpr scalar kFloor = 1e-12;
inline constexpr scalar dissFloor = 1e-12;
/// The turbulence equations are inner iterations of the pressure-velocity
/// loop, so they carry their own linear-solver settings rather than the
/// case's. 1e-8 is far below the accuracy the models themselves have.
inline constexpr scalar linearTol = 1e-8;
inline constexpr int linearMaxIter = 2000;
} // namespace turb

/// Turbulence settings of a case:
///   turbulence
///   {
///       model        kOmegaSST;   // laminar | kEpsilon | kOmegaSST
///       intensity    0.05;        // inlet turbulence intensity, u'/U
///       lengthScale  0.02;        // inlet turbulence length scale
///       relaxation   0.7;         // under-relaxation of k and of the
///                                 // dissipation variable
///   }
/// `lengthScale` <= 0 falls back to 7% of the mesh bounding box, the usual
/// hydraulic-diameter rule for a case that states no inlet condition.
struct TurbulenceConfig
{
    std::string model = "laminar";
    scalar intensity = 0.05;
    scalar lengthScale = -1;
    scalar relaxation = 0.7;

    /// A missing `turbulence` block means laminar.
    static TurbulenceConfig fromCase(const Dictionary& caseDict);

    /// Refuse a configuration the models cannot integrate with: a negative
    /// intensity seeds a negative k, and a relaxation factor outside (0, 1]
    /// either discards the equation (0) or amplifies every iterate (> 1).
    /// Called by fromCase() and by the factory, so a configuration built in
    /// code is checked the same way as one read from a case file.
    void validate() const;
};

/// Base of the two-equation RANS models. It owns what k-epsilon and k-omega
/// share - the turbulence-kinetic-energy transport, the wall layer, the eddy
/// viscosity and the face-wise effective viscosity the momentum equation is
/// assembled with - and leaves the dissipation-variable equation and the nut
/// formula to the subclasses.
///
/// Everything is kinematic (density 1), matching the incompressible solvers:
/// k and the dissipation variable are per unit mass and nut has the units of nu.
class RansModel : public IEddyViscosity
{
public:
    /// Factory over the case's turbulence{} block. The model keeps references
    /// to U and phi and reads them on every correctTurbulence(), so the caller
    /// only has to hand this object to the solver as the eddy-viscosity source.
    static std::unique_ptr<RansModel> create(const polyMesh& mesh,
                                             volVectorField& U,
                                             const surfaceScalarField& phi,
                                             scalar nu,
                                             const Schemes& schemes,
                                             const TurbulenceConfig& cfg);

    virtual ~RansModel() = default;

    /// Set the initial k and dissipation from the case's intensity and length
    /// scale. create() calls it, so an app that takes the factory does not have
    /// to; a test that wants a specific state calls it before it overwrites the
    /// fields.
    void init() { initializeFields(); }

    /// k0 and diss0 through the whole domain, prescribed at the inlet patches;
    /// walls keep k zeroGradient (its value follows from the wall layer) and get
    /// the dissipation Dirichlet per face, every other patch zeroGradient.
    /// initializeFields() uses it, and a test can seed an exact turbulence
    /// state with it.
    void initializeUniform(scalar k0, scalar diss0);

    void correctTurbulence(scalar dt) final;
    const Field<scalar>& nuEffFace() const final { return nuEffFace_; }

    virtual std::string name() const = 0;
    /// name of the dissipation variable, for reports and the vtu field
    virtual std::string dissipationName() const = 0;

    const polyMesh& mesh() const { return mesh_; }
    const Field<scalar>& nut() const { return nut_; }
    const Field<scalar>& k() const { return k_.internalField(); }
    const Field<scalar>& dissipation() const { return diss_.internalField(); }
    const volScalarField& kField() const { return k_; }
    const volScalarField& dissipationField() const { return diss_; }
    /// |S| = sqrt(2 S:S) per cell, from the current velocity gradient
    const Field<scalar>& strainRate() const { return strain_; }
    /// P_k = nut |S|^2 per cell, before any model limiter
    const Field<scalar>& production() const { return prod_; }

    /// Wall faces in a fixed order: every boundary face of the patches whose
    /// mesh type is "wall". The arrays below share that order.
    const std::vector<label>& wallFaces() const { return wallFaces_; }
    const Field<scalar>& wallUTau() const { return uTau_; }
    const Field<scalar>& wallYPlus() const { return yPlus_; }
    const Field<scalar>& wallNut() const { return nutWall_; }
    const Field<scalar>& wallDistance() const { return wallDist_; }
    /// y+ over the wall faces, excluding faces with no shear
    scalar maxWallYPlus() const;
    scalar meanWallYPlus() const;
    /// bulk friction velocity: the area-weighted mean of u_tau over the walls
    scalar meanFrictionVelocity() const;

    /// Distance from every cell centre to the nearest wall. Solved as a
    /// shortest path over the cell graph seeded at the wall-adjacent cells:
    /// exact on an orthogonal mesh, an upper bound otherwise, and cheap enough
    /// for a blend function that only needs it to within a cell size.
    const Field<scalar>& cellWallDistance() const;

protected:
    RansModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
              scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg,
              std::string dissName);

    /// k = 1.5 (uRef*intensity)^2 at the inlet patches with the matching
    /// dissipation from lengthScale, and elsewhere a value scaled by the local
    /// |U| so a cold start does not begin with a uniform k.
    virtual void initializeFields() = 0;
    /// false for the laminar model, which has neither equation
    virtual bool hasDissipation() const { return true; }
    /// transport of epsilon or omega with the current k, nut and velocity
    virtual void solveDissipation(scalar dt) = 0;
    /// nut from k and the dissipation variable
    virtual void updateNut() = 0;
    /// multiplier of nut in the k-equation diffusion coefficient
    /// nu + nut*kDiffusivity(): 1/sigma_eps for k-epsilon, sigma_k for SST
    virtual Field<scalar> kDiffusivity() const = 0;
    /// production of k as it enters the transport: the base uses the unlimited
    /// nut |S|^2, SST caps it (Menter's limiter)
    virtual Field<scalar> kProduction() const { return prod_; }
    /// volume-scaled, non-negative coefficient multiplying k on the lhs:
    /// the destruction term enters the equation as +sp(destructionSp())
    virtual Field<scalar> destructionSp() const = 0;
    /// Dirichlet value of the dissipation variable at a wall face
    virtual scalar wallDissipation(const WallState& ws, scalar kCell, scalar y) const = 0;

    /// Turbulence state the log layer is in when production and dissipation
    /// balance there: k = u_tau^2/sqrt(Cmu) with the matching dissipation.
    ///
    /// A wall-function mesh cannot reach that state by transport alone. The cell
    /// adjacent to the wall sees the boundary face value of U, so its strain
    /// rate is U_P/(dy/2) instead of the log-layer u_tau/(kappa y) -- a few
    /// times too big at y+ ~ 50 -- and its production nut |S|^2 then feeds back through
    /// nut = Cmu k^2/eps as a k^2 growth rate that the dissipation, which lags
    /// one iteration behind, cannot hold down. Models that use a wall function
    /// therefore replace the transported pair in that cell with the equilibrium
    /// one; the default replaces nothing, for a model that resolves the
    /// sublayer and needs no such closure.
    virtual void wallEquilibrium(const WallState&, scalar, scalar& kEq, scalar& dissEq) const
    {
        kEq = 0;
        dissEq = 0;
    }

    /// nu + nut*sigma at every face: harmonic mean across internal faces,
    /// owner-cell value at the other boundary faces, and at wall faces the
    /// wall-layer eddy viscosity. The momentum equation and both turbulence
    /// transports use this, so the diffusion out of the first cell always sees
    /// the shear the wall model prescribes instead of the cell's own nut.
    Field<scalar> diffusionFace(const Field<scalar>& gammaCell,
                                const Field<scalar>& sigmaCell) const;

    const polyMesh& mesh_;
    volVectorField& U_;
    const surfaceScalarField& phi_;
    scalar nu_;
    Schemes schemes_;
    TurbulenceConfig cfg_;

    volScalarField k_;
    volScalarField diss_;

    Field<scalar> nut_;
    Field<scalar> nutFace_;
    Field<scalar> nuEffFace_;
    Field<scalar> strain_;
    Field<scalar> prod_;

    std::vector<label> wallFaces_;
    Field<scalar> uTau_, yPlus_, nutWall_, wallDist_;

    /// reference velocity for the intensity/length-scale initialisation: the
    /// prescribed |U| of the inlet patches, or the current maximum |U|
    scalar uReference() const;
    /// turbulence length scale: cfg.lengthScale, or 7% of the box diagonal
    scalar lengthScale() const;

    /// nut interpolated to faces (wall faces from the wall layer) and the
    /// resulting nu + nut the momentum equation is assembled with
    void refreshFaceNut();

    /// Replace k and the dissipation variable with wallEquilibrium() in every
    /// cell that touches a wall face, after the transport has run and before
    /// nut is rebuilt from it. The transport still assembles over those cells --
    /// it is what carries the shear into the outer layer -- but its value in the
    /// wall cell is discarded, which is what makes the cell algebraic. Where a
    /// cell touches two walls the stronger shear sets both values.
    void imposeWallEquilibrium();

    mutable bool cellDistCalculated_ = false;
    mutable Field<scalar> cellWallDist_;

private:
    void computeStrainRate();
    void computeWallLayer();
    void solveK(scalar dt);
};

/// No modelling: nut stays zero, so the momentum equation keeps the molecular
/// viscosity, but the wall layer is still evaluated and a case can report the
/// y+ its mesh produces.
class LaminarModel : public RansModel
{
public:
    LaminarModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                 scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg);

    std::string name() const override { return "laminar"; }
    std::string dissipationName() const override { return "none"; }

protected:
    void initializeFields() override;
    bool hasDissipation() const override { return false; }
    void solveDissipation(scalar) override {}
    void updateNut() override;
    Field<scalar> kDiffusivity() const override;
    Field<scalar> destructionSp() const override;
    scalar wallDissipation(const WallState&, scalar, scalar) const override;
};

} // namespace cf
