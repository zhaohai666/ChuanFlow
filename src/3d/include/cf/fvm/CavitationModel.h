#pragma once

#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Field.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Homogeneous transport cavitation (F12): one volume fraction alpha_v of
/// vapour is advected with the solver's own flux and driven by a mass
/// transfer rate that reads the local pressure against the declared
/// saturation pressure,
///
///   dalpha/dt + div(phi, alpha) - laplacian(D, alpha)
///       = c(p) alpha (1 - alpha),
///   c(p) = (3/Rb) [ Ce sqrt(2 max(pSat - p, 0)/rhoL)
///                 - Cc sqrt(2 max(p - pSat, 0)/rhoL) ],
///
/// the Zwart-Gerber-Bernard form with both branches sharing the alpha(1-alpha)
/// saturation factor: vapour grows where the liquid is stretched below the
/// saturation pressure and collapses again where the diffuser recovers it.
/// The rate is linearised implicitly on alpha (the (1-alpha) factor frozen at
/// the old field), so the fraction stays in [0, 1] for any dt with c*dt < 1
/// and the solver clamps the remainder.
///
/// The mixture density rho = alpha rhoV + (1-alpha) rhoL feeds the pressure
/// correction through PisoControls::cellInvRho (the Poisson equation's
/// response 1/rho), while the momentum equation itself stays at unit density
/// - the declared simplification of the transport model at the small void
/// fractions these regressions reach.
class CavitationModel
{
public:
    /// The `cavitation{}` block of a case:
    ///   cavitation { pSat -80; rhoL 998; rhoV 1; Ce 2897; Cc 80;
    ///                Rb 1e-5; alpha0 1e-6; diffusivity 0; }
    /// pSat is the saturation pressure in the solver's kinematic units
    /// (m2/s2 = Pa/rhoL).
    struct Config
    {
        scalar pSat = 0;
        scalar rhoL = 998;
        scalar rhoV = 1;
        scalar Ce = 2897;
        scalar Cc = 80;
        scalar Rb = 1e-5;
        scalar alpha0 = 1e-6;
        scalar diffusivity = 0;

        static Config fromCase(const Dictionary& caseDict);
        void validate() const;
    };

    CavitationModel(const polyMesh& mesh, const surfaceScalarField& phi,
                    const Schemes& schemes, const Config& cfg);

    /// Transport alpha_v one step with the flux the coupling closed and the
    /// pressure the step ended on.
    void correct(scalar dt, const Field<scalar>& pCell);

    const volScalarField& alpha() const { return alpha_; }
    const Field<scalar>& alphaCell() const { return alpha_.internalField(); }
    /// alpha rhoV + (1-alpha) rhoL per cell
    Field<scalar> mixtureDensity() const;
    scalar alphaMax() const;
    scalar alphaMin() const;
    const Config& config() const { return cfg_; }

private:
    const polyMesh& mesh_;
    const surfaceScalarField& phi_;
    Schemes schemes_;
    Config cfg_;
    volScalarField alpha_;

    static constexpr scalar linearTol = 1e-10;
    static constexpr int linearMaxIter = 5000;
};

} // namespace cf
