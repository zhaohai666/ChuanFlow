#include "cf/turbulence/TransitionModel.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

TransitionModel::TransitionModel(const polyMesh& mesh, volVectorField& U,
                                 const surfaceScalarField& phi, scalar nu,
                                 const Schemes& schemes, const TurbulenceConfig& cfg)
    : KEpsilon(mesh, U, phi, nu, schemes, cfg), gamma_(mesh, "gamma", 0.0)
{
    if (!(cfg.reOnset > 0))
        cfFatal("gammaTransition: turbulence{ reOnset } must be a positive Re threshold");
    if (!(cfg.gammaRate > 0)) cfFatal("gammaTransition: gammaRate must be positive");
    rate_ = cfg.gammaRate;

    for (const polyMesh::Patch& p : mesh_.patches())
    {
        const label patch = mesh_.patchIdByName(p.name);
        gamma_.setPatchValue(patch, p.type == "inlet" ? "fixedValue" : "zeroGradient", 0.0);
    }
}

void TransitionModel::initializeFields()
{
    KEpsilon::initializeFields();
    for (scalar& v : gamma_.internalField()) v = 0.0;

    // x_t = nu Re_crit / U_ref: the declared Re_theta threshold read on the
    // flat plate through Re_x. Recomputed per run from the inlet velocity.
    const scalar u = uReference();
    if (!(u > 0)) cfFatal("gammaTransition: no velocity scale to place the onset");
    onsetX_ = nu_ * cfg_.reOnset / u;
    width_ = cfg_.onsetWidth > 0 ? cfg_.onsetWidth : 0.05 * onsetX_;
}

void TransitionModel::solveDissipation(scalar dt)
{
    KEpsilon::solveDissipation(dt);
    solveGamma(dt);
}

void TransitionModel::solveGamma(scalar dt)
{
    const Field<scalar> old = gamma_.internalField();
    const auto& vol = mesh_.cellVolumes();
    const auto& cc = mesh_.cellCentres();
    const std::size_t nC = static_cast<std::size_t>(mesh_.nCells());
    const scalar u = uReference();
    const scalar amp = rate_ * u / onsetX_;

    fvMatrix eqn(mesh_, viewOf(gamma_), "gammaEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, old, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);
    Field<scalar> gamma(nC);
    for (std::size_t c = 0; c < nC; ++c)
    {
        gamma[c] = nu_ + nut_[c]; // sigma_gamma = 1
    }
    fvm::laplacianFace(eqn, harmonicFaceGamma(mesh_, gamma), schemes_.laplacian, -1);

    // S = a(x) (1 - gamma): the growing branch (a > 0) is linearised
    // implicitly on gamma, the decaying branch upstream stays explicit.
    Field<scalar> sp(nC, 0.0);
    Field<scalar> rhs(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
    {
        const scalar a = amp * std::tanh((cc[c].x - onsetX_) / width_);
        if (a > 0)
        {
            sp[c] = vol[c] * a;
            rhs[c] = vol[c] * a;
        }
        else
        {
            rhs[c] = vol[c] * a * (1.0 - old[c]);
        }
    }
    fvm::sp(eqn, sp, 1);
    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        eqn.addRhs(c, 0, rhs[static_cast<std::size_t>(c)]);
    }

    eqn.solve(fvMatrix::Controls{turb::linearTol, turb::linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    auto& g = gamma_.internalField();
    for (std::size_t c = 0; c < nC; ++c)
    {
        g[c] = std::min(std::max(old[c] + rl * (x[c] - old[c]), 0.0), 1.0);
    }
}

void TransitionModel::updateNut()
{
    KEpsilon::updateNut();
    for (std::size_t c = 0; c < nut_.size(); ++c)
    {
        nut_[c] *= gamma_.internalField()[c];
    }
}

} // namespace cf
