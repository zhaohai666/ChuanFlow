#pragma once

#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{

/// Menter's shear-stress-transport k-omega model in kinematic form:
///   dk/dt + div(phi,k) - laplacian(nu + nut*sigmaK, k) = Pk - betaStar k omega
///   dw/dt + div(phi,w) - laplacian(nu + nut*sigmaW, w)
///        = gamma |S|^2 - beta w^2 + (1-F1) 2 sigmaW2/(betaStar w) grad(w).grad(k)
///   nut = a1 k / max(a1 w, |S| F2)
/// The blend F1 selects the k-omega formulation near walls (F1 = 1) and the
/// transformed k-epsilon one in the free stream (F1 = 0), which is what makes
/// the model usable without a wall function; the wall faces still get the
/// omega value implied by the law of the wall, blended with the viscous-
/// sublayer limit, so a mesh that resolves the sublayer and one that does not
/// both work.
class KOmegaSST : public RansModel
{
public:
    KOmegaSST(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
              scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg);

    std::string name() const override { return "kOmegaSST"; }
    std::string dissipationName() const override { return "omega"; }

    static constexpr scalar betaStar = 0.09;
    static constexpr scalar beta1 = 0.075;
    static constexpr scalar beta2 = 0.0828;
    static constexpr scalar sigmaK1 = 0.85;
    static constexpr scalar sigmaK2 = 1.0;
    static constexpr scalar sigmaW1 = 0.5;
    static constexpr scalar sigmaW2 = 0.856;
    static constexpr scalar a1 = 0.31;
    /// zeta* of the blended production/destruction relation gamma = beta/betaStar
    /// - zeta*/sigmaW, which keeps the free-stream eddy viscosity at Cmu = 0.09
    static constexpr scalar zetaStar = 0.3;
    static constexpr scalar gamma1 = beta1 / betaStar - zetaStar / sigmaW1;
    static constexpr scalar gamma2 = beta2 / betaStar - zetaStar / sigmaW2;
    /// Menter's production limiter: Pk <= 10 betaStar k w
    static constexpr scalar productionLimit = 10.0;

    /// blend functions, exposed for the unit tests. `cellSize` is the largest
    /// cell edge of the mesh, Menter's free-stream guard.
    static scalar f1(scalar k, scalar omega, scalar y, scalar nu, scalar cellSize);
    static scalar f2(scalar k, scalar omega, scalar y, scalar nu);

    /// largest cell volume to the 1/3, used as the cell size of F1's third
    /// argument (a cube root is a mesh-size estimate that needs no face data)
    scalar cellSize() const { return cellSize_; }

protected:
    void initializeFields() override;
    void solveDissipation(scalar dt) override;
    void updateNut() override;
    void updateTurbDissipation() override;
    Field<scalar> kDiffusivity() const override;
    Field<scalar> kProduction() const override;
    Field<scalar> destructionSp() const override;
    scalar wallDissipation(const WallState& ws, scalar kCell, scalar y) const override;
    /// log-layer equilibrium of the wall-function form: k = u_tau^2/sqrt(beta*)
    /// with omega = u_tau/(sqrt(beta*) kappa y), the state at which the
    /// production and destruction of both equations balance in the log layer
    void wallEquilibrium(const WallState& ws, scalar y,
                         scalar& kEq, scalar& dissEq) const override;

private:
    scalar cellSize_ = 1;
};

} // namespace cf
