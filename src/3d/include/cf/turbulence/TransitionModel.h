#pragma once

#include "cf/common/Field.h"
#include "cf/fvm/Fields.h"
#include "cf/turbulence/KEpsilon.h"

namespace cf
{

/// Intermittency transport on the k-epsilon base: a single scalar gamma in
/// [0, 1] marks the fraction of time a control volume is turbulent, and the
/// eddy viscosity is the RANS one scaled by it,
///
///   nut_model = gamma * nut_kEpsilon
///
/// so the layer stays laminar (nut = 0) until gamma grows and then relaxes to
/// the full model. The growth source is the declared onset correlation: the
/// transition location follows from the momentum-thickness Reynolds number
/// threshold reOnset through the flat-plate reading Re_x = U_ref x / nu,
///
///   x_t = nu * reOnset / U_ref,
///   S_gamma = gammaRate (U_ref / x_t) tanh((x - x_t)/width) (1 - gamma),
///
/// a smooth ramp centred on x_t. The tanh is odd about x_t and the source
/// saturates at (1 - gamma), so the steady profile passes gamma = 0.5 at x_t
/// whatever the ramp rate; width defaults to 5% of x_t. gamma is transported
/// with the same upwind flux as k and diffused with nu + nut.
class TransitionModel : public KEpsilon
{
public:
    TransitionModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                    scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg);

    std::string name() const override { return "gammaTransition"; }
    const Field<scalar>& gamma() const { return gamma_.internalField(); }
    /// declared onset location x_t (0 when the case declares no threshold)
    scalar onsetX() const { return onsetX_; }

protected:
    void initializeFields() override;
    void solveDissipation(scalar dt) override;
    void updateNut() override;
    const Field<scalar>* auxField() const override { return &gamma_.internalField(); }
    std::string auxFieldName() const override { return "gamma"; }

private:
    void solveGamma(scalar dt);

    volScalarField gamma_;
    scalar onsetX_ = 0;
    scalar width_ = 0;
    scalar rate_ = 4.0;
};

} // namespace cf
