#pragma once

#include "cf/common/Field.h"
#include "cf/turbulence/KEpsilon.h"

namespace cf
{

/// Detached-eddy simulation on the k-epsilon base: the RANS model transports
/// k and epsilon as usual, and the eddy viscosity is then rebuilt from the
/// DES length-scale shield
///
///   L = min(kappa d_wall, C_DES Delta),   Delta = V^(1/3)
///   eps = k^(3/2)/L,                      nut = Cmu k^2/eps = Cmu k^(1/2) L
///
/// Near a wall (small d_wall) the shield selects the RANS length kappa d and
/// the model behaves as the wall-modeled k-epsilon it is; where the mesh is
/// finer than the local turbulent length (small Delta, away from the wall) it
/// selects the LES length C_DES Delta and the eddy viscosity becomes the
/// Smagorinsky-like Cmu k^(1/2) C_DES Delta of a one-equation LES. The shield
/// is monotone in both arguments by construction: refining the mesh can only
/// reduce nut, never increase it past the RANS value.
///
/// The transported epsilon keeps its own equation (the shield is applied to
/// the dissipation used by nut and reported through turbDissipation()), which
/// is the standard delayed-detached-eddy reading of the hybrid as a length
/// limiter rather than a replacement of the RANS destruction.
class DesModel : public KEpsilon
{
public:
    DesModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
             scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg);

    std::string name() const override { return "des"; }
    /// the shielded length L per cell, written to the vtu as `desLength`
    const Field<scalar>& desLength() const { return length_; }

protected:
    void updateNut() override;
    const Field<scalar>* auxField() const override { return &length_; }
    std::string auxFieldName() const override { return "desLength"; }

private:
    scalar cDES_ = 0.61;
    Field<scalar> length_;
};

} // namespace cf
