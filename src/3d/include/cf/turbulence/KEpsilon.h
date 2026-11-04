#pragma once

#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{

/// Standard k-epsilon (Launder & Sharma) in kinematic form with algebraic wall
/// functions:
///   dk/dt + div(phi,k) - laplacian(nu + nut/sigmaK, k) = P - eps
///   deps/dt + div(phi,eps) - laplacian(nu + nut/sigmaE, eps)
///        = C1 eps P/k - C2 eps^2/k
///   nut = Cmu k^2/eps
/// The wall faces carry the law-of-the-wall values: k keeps a zero gradient
/// (its cell value follows from the wall shear) while eps is prescribed as
/// u_tau^3/(kappa y), the homogeneous-decay value of the log layer, so the mesh
/// needs no sublayer resolution.
class KEpsilon : public RansModel
{
public:
    KEpsilon(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
             scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg);

    std::string name() const override { return "kEpsilon"; }
    std::string dissipationName() const override { return "epsilon"; }

    static constexpr scalar Cmu = 0.09;
    static constexpr scalar C1 = 1.44;
    static constexpr scalar C2 = 1.92;
    static constexpr scalar sigmaK = 1.0;
    static constexpr scalar sigmaE = 1.3;

protected:
    void initializeFields() override;
    void solveDissipation(scalar dt) override;
    void updateNut() override;
    void updateTurbDissipation() override;
    Field<scalar> kDiffusivity() const override;
    Field<scalar> destructionSp() const override;
    scalar wallDissipation(const WallState& ws, scalar kCell, scalar y) const override;
    /// log-layer equilibrium of the standard model: k = u_tau^2/sqrt(Cmu) with
    /// eps = u_tau^3/(kappa y), the same shear the wall function prescribes
    void wallEquilibrium(const WallState& ws, scalar y,
                         scalar& kEq, scalar& dissEq) const override;
};

} // namespace cf
