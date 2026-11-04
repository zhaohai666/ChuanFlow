#include "cf/turbulence/KOmegaSST.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Log.h"

namespace cf
{
namespace
{

/// Cmu of the k-epsilon side of the blend, used to turn a length scale into an
/// omega scale: omega = k^(1/2)/(Cmu^(1/4) L).
constexpr scalar cmuBlend = 0.09;

} // namespace

KOmegaSST::KOmegaSST(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                     scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg)
    : RansModel(mesh, U, phi, nu, schemes, cfg, "omega")
{
    scalar v = 0;
    for (const scalar c : mesh.cellVolumes()) v = std::max(v, c);
    cellSize_ = std::cbrt(std::max(v, 1e-300));
}

scalar KOmegaSST::f1(scalar k, scalar omega, scalar y, scalar nu, scalar cellSize)
{
    const scalar kk = std::max(k, turb::kFloor);
    const scalar w = std::max(omega, turb::dissFloor);
    const scalar yy = std::max(y, 1e-12);
    // Menter's three-argument blend: the eddy-viscosity and viscous-sublayer
    // guards, plus the free-stream guard that needs the near-wall cell size so
    // F1 does not switch to the epsilon formulation inside a coarse shear layer
    const scalar arg = std::max(
        std::max(std::sqrt(kk) / (betaStar * w * yy), 500.0 * nu / (yy * yy * w)),
        4.0 * sigmaW2 * kk / (w * yy * yy * cellSize));
    return std::tanh(std::pow(arg, 4.0));
}

scalar KOmegaSST::f2(scalar k, scalar omega, scalar y, scalar nu)
{
    (void)nu;
    const scalar kk = std::max(k, turb::kFloor);
    const scalar w = std::max(omega, turb::dissFloor);
    const scalar yy = std::max(y, 1e-12);
    const scalar arg = 2.0 * std::sqrt(kk) / (betaStar * w * yy);
    return std::tanh(arg * arg);
}

void KOmegaSST::initializeFields()
{
    const scalar u = uReference();
    const scalar k0 = 1.5 * std::pow(u * cfg_.intensity, 2);
    const scalar w0 = std::sqrt(std::max(k0, turb::kFloor))
                    / (std::pow(cmuBlend, 0.25) * lengthScale());
    initializeUniform(k0, w0);
    log::info("turbulence kOmegaSST: uRef=", u, " k0=", k0, " omega0=", w0,
              " lengthScale=", lengthScale());
}

Field<scalar> KOmegaSST::kProduction() const
{
    // Menter's limiter: the eddy viscosity may not extract more than ten times
    // the equilibrium destruction rate, which caps nut in strong shear
    Field<scalar> p(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < p.size(); ++c)
    {
        const scalar k = std::max(k_.internalField()[c], turb::kFloor);
        const scalar w = std::max(diss_.internalField()[c], turb::dissFloor);
        p[c] = std::min(prod_[c], productionLimit * betaStar * k * w);
    }
    return p;
}

Field<scalar> KOmegaSST::kDiffusivity() const
{
    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    const auto& y = cellWallDistance();
    Field<scalar> s(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar f1c = f1(k_.internalField()[c], diss_.internalField()[c], y[c], nu_, cellSize_);
        s[c] = f1c * sigmaK1 + (1.0 - f1c) * sigmaK2;
    }
    return s;
}

Field<scalar> KOmegaSST::destructionSp() const
{
    const auto& vol = mesh_.cellVolumes();
    Field<scalar> sp(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < sp.size(); ++c)
    {
        sp[c] = vol[c] * betaStar * std::max(diss_.internalField()[c], turb::dissFloor);
    }
    return sp;
}

void KOmegaSST::solveDissipation(scalar dt)
{
    const Field<scalar> wOld = diss_.internalField();
    const auto& vol = mesh_.cellVolumes();
    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    const auto& y = cellWallDistance();

    // blend functions of the current field
    Field<scalar> f1c(nC), sigmaW(nC), gammaW(nC), betaW(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        f1c[c] = f1(k_.internalField()[c], wOld[c], y[c], nu_, cellSize_);
        sigmaW[c] = f1c[c] * sigmaW1 + (1.0 - f1c[c]) * sigmaW2;
        gammaW[c] = f1c[c] * gamma1 + (1.0 - f1c[c]) * gamma2;
        betaW[c] = f1c[c] * beta1 + (1.0 - f1c[c]) * beta2;
    }

    fvMatrix eqn(mesh_, viewOf(diss_), "omegaEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, wOld, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);

    Field<scalar> gamma(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        gamma[c] = nu_ + nut_[c] * sigmaW[c];
    }
    fvm::laplacianFace(eqn, diffusionFace(gamma, sigmaW), schemes_.laplacian, -1);

    // cross diffusion (1-F1) 2 sigmaW2 /(betaStar w) grad(w).grad(k): explicit
    // and floored at zero, so an anti-parallel gradient pair cannot turn it
    // into a destruction term the linear solver would have to fight
    const Field<Vector> gradW = fvc::grad(diss_, schemes_.grad);
    const Field<Vector> gradK = fvc::grad(k_, schemes_.grad);

    Field<scalar> destr(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar w = std::max(wOld[c], turb::dissFloor);
        eqn.addRhs(static_cast<label>(c), 0, vol[c] * gammaW[c] * strain_[c] * strain_[c]);
        const scalar cd = (1.0 - f1c[c]) * 2.0 * sigmaW2 / (betaStar * w)
                        * dot(gradW[c], gradK[c]);
        if (cd > 0) eqn.addRhs(static_cast<label>(c), 0, vol[c] * cd);
        destr[c] = vol[c] * betaW[c] * w;
    }
    fvm::sp(eqn, destr, 1);

    eqn.solve(fvMatrix::Controls{turb::linearTol, turb::linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    for (std::size_t c = 0; c < diss_.internalField().size(); ++c)
    {
        const scalar wNew = wOld[c] + rl * (x[c] - wOld[c]);
        diss_.internalField()[c] = std::max(wNew, turb::dissFloor);
    }
}

void KOmegaSST::updateNut()
{
    const auto& y = cellWallDistance();
    for (std::size_t c = 0; c < nut_.size(); ++c)
    {
        const scalar k = std::max(k_.internalField()[c], turb::kFloor);
        const scalar w = std::max(diss_.internalField()[c], turb::dissFloor);
        const scalar f2c = f2(k, w, y[c], nu_);
        const scalar denom = std::max(a1 * w, strain_[c] * f2c);
        nut_[c] = std::max(a1 * k / denom, 0.0);
    }
}

void KOmegaSST::updateTurbDissipation()
{
    // the model transports omega; epsilon is the k-equation's destruction term
    const Field<scalar>& k = k_.internalField();
    const Field<scalar>& w = diss_.internalField();
    for (std::size_t c = 0; c < eps_.size(); ++c) eps_[c] = betaStar * k[c] * w[c];
}

scalar KOmegaSST::wallDissipation(const WallState& ws, scalar kCell, scalar y) const
{
    // viscous-sublayer limit omega -> 6 nu/(beta1 y^2) and log-layer value
    // omega -> u_tau/(sqrt(betaStar) kappa y); combining them in quadrature
    // gives the right asymptote in both limits on any near-wall spacing
    const scalar wVis = 6.0 * nu_ / (beta1 * y * y);
    const scalar wLog = ws.uTau / (std::sqrt(betaStar) * turb::kappa * y);
    const scalar wBlend = std::sqrt(wVis * wVis + wLog * wLog);
    (void)kCell;
    return std::max(wBlend, turb::dissFloor);
}

void KOmegaSST::wallEquilibrium(const WallState& ws, scalar y,
                                scalar& kEq, scalar& dissEq) const
{
    kEq = ws.uTau * ws.uTau / std::sqrt(betaStar);
    dissEq = ws.uTau / (std::sqrt(betaStar) * turb::kappa * y);
}

} // namespace cf
