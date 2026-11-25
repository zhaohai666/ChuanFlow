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
}

scalar KOmegaSST::f1(scalar k, scalar omega, scalar y, scalar nu, scalar cdKw)
{
    const scalar kk = std::max(k, turb::kFloor);
    const scalar w = std::max(omega, turb::dissFloor);
    const scalar yy = std::max(y, 1e-12);
    // Menter 1994 (Eq. 10): the outer operator is a MIN, not a max -- the
    // third argument caps arg1 in the free stream, where the cross-diffusion
    // guard CD_komega is active, and opens (CD at its floor, the argument
    // huge) where the gradients vanish so the two wall arguments decide.
    const scalar cd = std::max(cdKw, 1e-10);
    const scalar arg = std::min(
        std::max(std::sqrt(kk) / (betaStar * w * yy), 500.0 * nu / (yy * yy * w)),
        4.0 * sigmaW2 * kk / (cd * yy * yy));
    return std::tanh(std::pow(arg, 4.0));
}

scalar KOmegaSST::f2(scalar k, scalar omega, scalar y, scalar nu)
{
    const scalar kk = std::max(k, turb::kFloor);
    const scalar w = std::max(omega, turb::dissFloor);
    const scalar yy = std::max(y, 1e-12);
    // Menter 1994 (Eq. 8): arg2 = max(2 sqrt(k)/(beta* w y), 500 nu/(y^2 w));
    // the viscous term keeps F2 = 1 in the sublayer on any near-wall spacing
    const scalar arg = std::max(2.0 * std::sqrt(kk) / (betaStar * w * yy),
                                500.0 * nu / (yy * yy * w));
    return std::tanh(arg * arg);
}

Field<scalar> KOmegaSST::cdKOmega() const
{
    const Field<Vector> gw = fvc::grad(diss_, schemes_.grad);
    const Field<Vector> gk = fvc::grad(k_, schemes_.grad);
    Field<scalar> cd(gw.size());
    for (std::size_t c = 0; c < cd.size(); ++c)
    {
        const scalar w = std::max(diss_.internalField()[c], turb::dissFloor);
        cd[c] = std::max(2.0 * sigmaW2 / w * dot(gw[c], gk[c]), 1e-10);
    }
    return cd;
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

scalar KOmegaSST::omegaProduction(scalar gamma, scalar strain, scalar omega,
                                  scalar f2)
{
    // OpenFOAM kOmegaSSTBase.C / Menter 2003 consistent production:
    //   P_omega = gamma min(S^2, (c1/a1) betaStar omega
    //                      max(a1 omega, F2 S)).
    // The F2 S branch is required when the Bradshaw limiter makes
    // nut = a1 k/max(a1 omega, F2 S) smaller than k/omega.
    const scalar s = std::max(strain, 0.0);
    const scalar w = std::max(omega, turb::dissFloor);
    const scalar bradshaw = std::max(a1 * w, f2 * s);
    const scalar cap = productionLimit / a1 * betaStar * w * bradshaw;
    return gamma * std::min(s * s, cap);
}

Field<scalar> KOmegaSST::kDiffusivity() const
{
    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    const auto& y = cellWallDistance();
    const Field<scalar> cd = cdKOmega();
    Field<scalar> s(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar f1c = f1(k_.internalField()[c], diss_.internalField()[c], y[c], nu_, cd[c]);
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

    // blend functions of the current field; the CD_komega guard of F1's
    // third argument comes from the same gradients as the cross-diffusion
    // source below, so they are computed once here
    const Field<Vector> gradW = fvc::grad(diss_, schemes_.grad);
    const Field<Vector> gradK = fvc::grad(k_, schemes_.grad);
    Field<scalar> f1c(nC), f2c(nC), cdKw(nC), sigmaW(nC), gammaW(nC), betaW(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar wGuard = std::max(wOld[c], turb::dissFloor);
        cdKw[c] = std::max(2.0 * sigmaW2 / wGuard * dot(gradW[c], gradK[c]), 1e-10);
        f1c[c] = f1(k_.internalField()[c], wOld[c], y[c], nu_, cdKw[c]);
        f2c[c] = f2(k_.internalField()[c], wOld[c], y[c], nu_);
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

    // cross diffusion (1-F1) 2 sigmaW2/w grad(w).grad(k) (Menter 2003):
    // explicit, with the full sign kept -- a negative CD is a legitimate
    // additional sink of omega in the k-epsilon branch, matching the
    // reference implementation's signed source. The gradients are the ones
    // taken above for the F1 guard.

    Field<scalar> destr(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar w = std::max(wOld[c], turb::dissFloor);
        eqn.addRhs(static_cast<label>(c), 0,
                   vol[c] * omegaProduction(gammaW[c], strain_[c], w, f2c[c]));
        const scalar cd = (1.0 - f1c[c]) * 2.0 * sigmaW2 / w
                        * dot(gradW[c], gradK[c]);
        eqn.addRhs(static_cast<label>(c), 0, vol[c] * cd);
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
