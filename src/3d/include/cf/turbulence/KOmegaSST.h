#pragma once

#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{

/// Menter's shear-stress-transport k-omega model in kinematic form:
///   dk/dt + div(phi,k) - laplacian(nu + nut*sigmaK, k) = Pk - betaStar k omega
///   dw/dt + div(phi,w) - laplacian(nu + nut*sigmaW, w)
///        = gamma |S|^2 - beta w^2 + (1-F1) 2 sigmaW2/w grad(w).grad(k)
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
    /// Menter's gamma relation (Menter 2003, AIAA J 41(2)): gamma_i = beta_i/betaStar
    /// - sigmaWi kappa^2/sqrt(betaStar), which ties the production/destruction
    /// balance to the von Karman constant so the free-stream eddy viscosity
    /// settles at Cmu = 0.09
    static constexpr scalar sqrtBetaStar = 0.3; // sqrt(0.09), exact in decimal
    static constexpr scalar gamma1 = beta1 / betaStar - sigmaW1 * turb::kappa * turb::kappa / sqrtBetaStar;
    static constexpr scalar gamma2 = beta2 / betaStar - sigmaW2 * turb::kappa * turb::kappa / sqrtBetaStar;
    /// Menter's production limiter: Pk <= 10 betaStar k w
    static constexpr scalar productionLimit = 10.0;

    /// blend functions, exposed for the unit tests. `cdKw` is Menter's
    /// cross-diffusion guard CD_komega = max(2 sigmaW2/omega grad(w).grad(k),
    /// 1e-10) of the cell (see cdKOmega()); F1's third argument caps arg1 in
    /// the free stream, where the guard is active, and opens when the
    /// gradients vanish (CD at its floor).
    static scalar f1(scalar k, scalar omega, scalar y, scalar nu, scalar cdKw);
    static scalar f2(scalar k, scalar omega, scalar y, scalar nu);

    /// Menter's CD_komega field from the current k/omega gradients, floored
    /// at 1e-10 so a zero-gradient region leaves the free-stream guard open
    Field<scalar> cdKOmega() const;

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
};

} // namespace cf
