#pragma once

#include <string>

#include "cf/common/Field.h"
#include "cf/fvm/EddyViscosity.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Large-eddy simulation sub-grid viscosity: the filtered momentum equation
/// keeps the molecular viscosity and models the unresolved transport with an
/// eddy viscosity built from the *resolved* strain rate and the local filter
/// width Delta = V^(1/3) (the cube root of the cell volume, the standard
/// box filter of a finite-volume mesh).
///
///   Smagorinsky: nut = (Cs Delta)^2 |S|
///   WALE:        nut = (Cw Delta)^2 (Sd:Sd)^(3/2) /
///                       ( (S:S)^(5/2) + (Sd:Sd)^(5/4) )
///                with Sd = sym(grad U)^2_traceless, the tensor the wall-
///                adapted local eddy viscosity is built from. Sd vanishes
///                identically for a unidirectional shear (grad U squared is
///                zero there), which is WALE's quadratic near-wall damping:
///                the model returns nut = 0 for pure shear instead of the
///                Smagorinsky (Cs Delta)^2 |S| that over-drags a wall layer.
///
/// The model carries no transported equation, so correctTurbulence() only
/// re-evaluates the algebraic viscosity from the current velocity; it plugs
/// into the pressure-velocity solvers through IEddyViscosity exactly like the
/// RANS models. There is no wall function: a wall face keeps the owner-cell
/// viscosity, so the mesh has to resolve the layer the LES is asked to carry.
class LesModel : public IEddyViscosity
{
public:
    /// The `les{}` block of a case:
    ///   les { model smagorinsky | wale; Cs 0.17; Cw 0.325; }
    struct Config
    {
        std::string model = "smagorinsky";
        scalar Cs = 0.17;   // Smagorinsky constant
        scalar Cw = 0.325;  // WALE constant

        static Config fromCase(const Dictionary& caseDict);
        void validate() const;
    };

    LesModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
             scalar nu, const Schemes& schemes, const Config& cfg);

    void correctTurbulence(scalar dt) override;
    const Field<scalar>& nuEffFace() const override { return nuEffFace_; }

    std::string name() const { return "les-" + cfg_.model; }
    const Config& config() const { return cfg_; }
    const Field<scalar>& nut() const { return nut_; }
    const Field<scalar>& strainRate() const { return strain_; }
    /// filter width Delta = V^(1/3) per cell
    const Field<scalar>& filterWidth() const { return delta_; }

private:
    void refreshFaceNut();

    const polyMesh& mesh_;
    volVectorField& U_;
    scalar nu_;
    Schemes schemes_;
    Config cfg_;

    Field<scalar> nut_;
    Field<scalar> nuEffFace_;
    Field<scalar> strain_;
    Field<scalar> delta_;
};

} // namespace cf
