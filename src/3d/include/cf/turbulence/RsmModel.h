#pragma once

#include <string>

#include "cf/common/Field.h"
#include "cf/fvm/EddyViscosity.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Second-moment closure: transport of the six independent Reynolds stresses
/// R_ij = <u_i u_j> (kinematic, per unit mass) plus the dissipation rate.
///
///   dR_ij/dt + div(phi, R_ij) - laplacian(nu + nut/sigmaR, R_ij)
///        = P_ij - C1 (eps/k)(R_ij - 2/3 k delta_ij) - 2/3 eps delta_ij
///   deps/dt + div(phi, eps) - laplacian(nu + nut/sigmaE, eps)
///        = C1e (eps/k) P_k - C2e eps^2/k
///
/// The pressure-strain correlation is the LRR return-to-isotropy term: the
/// slow linear regression -C1 (eps/k) a_ij on the anisotropy a_ij, with the
/// wall-reflection correction omitted (declared: the closure is used on
/// meshes and cases without a wall layer, see the case notes). Production is
/// the exact stress contraction P_ij = -(R_jk dU_i/dx_k + R_ik dU_j/dx_k), so
/// the model carries the anisotropy the eddy-viscosity models cannot.
///
/// The momentum equation still sees a scalar eddy viscosity (the pressure-
/// velocity solvers take IEddyViscosity): nut = -(R:S)/|S|^2, the value that
/// reproduces the modelled shear stress of the current mean strain exactly;
/// where the strain vanishes it falls back to 0 (no stress to represent).
/// Realizability is kept by clamping each off-diagonal stress to the Cauchy-
/// Schwarz bound sqrt(R_ii R_jj) and the diagonal stresses to non-negativity.
class RsmModel : public IEddyViscosity
{
public:
    /// The `rsm{}` block of a case:
    ///   rsm { C1 1.8; C2 0.55; C1eps 1.44; C2eps 1.92;
    ///         sigmaR 1.0; sigmaEps 1.3; relaxation 0.7;
    ///         k0 0.01; eps0 0.001; }
    /// C2 (the rapid pressure-strain regression) is accepted but held at 0 in
    /// the closure below: only the slow LRR term is claimed.
    struct Config
    {
        scalar C1 = 1.8;
        scalar C2 = 0.0;
        scalar C1eps = 1.44;
        scalar C2eps = 1.92;
        scalar sigmaR = 1.0;
        scalar sigmaEps = 1.3;
        scalar relaxation = 0.7;
        scalar k0 = 1e-2;   // initial turbulent kinetic energy
        scalar eps0 = 1e-3; // initial dissipation rate

        static Config fromCase(const Dictionary& caseDict);
        void validate() const;
    };

    RsmModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
             scalar nu, const Schemes& schemes, const Config& cfg);

    void correctTurbulence(scalar dt) override;
    const Field<scalar>& nuEffFace() const override { return nuEffFace_; }

    std::string name() const { return "rsm"; }

    /// Stress order: 0=11 1=22 2=33 3=12 4=13 5=23.
    const volScalarField& stress(int comp) const { return R_[comp]; }
    const Field<scalar>& epsilon() const { return eps_.internalField(); }
    const Field<scalar>& k() const { return k_; }
    const Field<scalar>& nut() const { return nut_; }
    const Field<scalar>& production() const { return pk_; }
    const Config& config() const { return cfg_; }

private:
    void solveStress(int comp, scalar dt, const Field<scalar>* prod,
                     const Field<scalar>& kCell);
    void solveEpsilon(scalar dt);
    void updateNut();
    void refreshFaceNut();

    const polyMesh& mesh_;
    volVectorField& U_;
    const surfaceScalarField& phi_;
    scalar nu_;
    Schemes schemes_;
    Config cfg_;

    volScalarField R_[6];
    volScalarField eps_;

    Field<scalar> k_;
    Field<scalar> nut_;
    Field<scalar> nuEffFace_;
    Field<scalar> pk_;

    static constexpr scalar kFloor = 1e-12;
    static constexpr scalar dissFloor = 1e-12;
    static constexpr scalar linearTol = 1e-9;
    static constexpr int linearMaxIter = 2000;
};

} // namespace cf
