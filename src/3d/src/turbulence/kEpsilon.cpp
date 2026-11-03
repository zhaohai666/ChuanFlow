#include "cf/turbulence/KEpsilon.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Log.h"

namespace cf
{

KEpsilon::KEpsilon(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                   scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg)
    : RansModel(mesh, U, phi, nu, schemes, cfg, "epsilon")
{
}

void KEpsilon::initializeFields()
{
    // k from the turbulence intensity, eps from the length scale through the
    // equilibrium of the log layer: k ~ u_tau^2, eps ~ u_tau^3/(kappa y) and
    // L ~ k^1.5/eps give the standard Cmu^0.25 k^1.5/L.
    const scalar u = uReference();
    const scalar k0 = 1.5 * std::pow(u * cfg_.intensity, 2);
    const scalar e0 = std::pow(Cmu, 0.25) * std::pow(std::max(k0, turb::kFloor), 1.5)
                    / lengthScale();
    initializeUniform(k0, e0);
    log::info("turbulence kEpsilon: uRef=", u, " k0=", k0, " eps0=", e0,
              " lengthScale=", lengthScale());
}

Field<scalar> KEpsilon::kDiffusivity() const
{
    return Field<scalar>(static_cast<std::size_t>(mesh_.nCells()), 1.0 / sigmaK);
}

Field<scalar> KEpsilon::destructionSp() const
{
    const auto& vol = mesh_.cellVolumes();
    Field<scalar> sp(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < sp.size(); ++c)
    {
        // eps = (eps/k) k on the lhs: the ratio is bounded because a cell with
        // no turbulence has both k and eps at their floors
        sp[c] = vol[c] * diss_.internalField()[c] / std::max(k_.internalField()[c], turb::kFloor);
    }
    return sp;
}

void KEpsilon::solveDissipation(scalar dt)
{
    const Field<scalar> eOld = diss_.internalField();
    const auto& vol = mesh_.cellVolumes();

    fvMatrix eqn(mesh_, viewOf(diss_), "epsilonEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, eOld, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);

    Field<scalar> sigma(static_cast<std::size_t>(mesh_.nCells()), 1.0 / sigmaE);
    Field<scalar> gamma(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < gamma.size(); ++c)
    {
        gamma[c] = nu_ + nut_[c] * sigma[c];
    }
    fvm::laplacianFace(eqn, diffusionFace(gamma, sigma), schemes_.laplacian, -1);

    const Field<scalar> pk = kProduction();
    Field<scalar> destr(static_cast<std::size_t>(mesh_.nCells()));
    for (std::size_t c = 0; c < gamma.size(); ++c)
    {
        const scalar k = std::max(k_.internalField()[c], turb::kFloor);
        const scalar e = std::max(eOld[c], turb::dissFloor);
        eqn.addRhs(static_cast<label>(c), 0, vol[c] * C1 * e * pk[c] / k);
        destr[c] = vol[c] * C2 * e / k;
    }
    fvm::sp(eqn, destr, 1);

    eqn.solve(fvMatrix::Controls{turb::linearTol, turb::linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    const scalar rl = cfg_.relaxation;
    for (std::size_t c = 0; c < diss_.internalField().size(); ++c)
    {
        const scalar eNew = eOld[c] + rl * (x[c] - eOld[c]);
        diss_.internalField()[c] = std::max(eNew, turb::dissFloor);
    }
}

void KEpsilon::updateNut()
{
    for (std::size_t c = 0; c < nut_.size(); ++c)
    {
        const scalar k = std::max(k_.internalField()[c], turb::kFloor);
        const scalar e = std::max(diss_.internalField()[c], turb::dissFloor);
        nut_[c] = std::max(Cmu * k * k / e, 0.0);
    }
}

scalar KEpsilon::wallDissipation(const WallState& ws, scalar /*kCell*/, scalar y) const
{
    // homogeneous decay of the log layer: eps = u_tau^3/(kappa y)
    const scalar e = std::pow(ws.uTau, 3) / (turb::kappa * y);
    return std::max(e, turb::dissFloor);
}

void KEpsilon::wallEquilibrium(const WallState& ws, scalar y,
                               scalar& kEq, scalar& dissEq) const
{
    kEq = ws.uTau * ws.uTau / std::sqrt(Cmu);
    dissEq = std::pow(ws.uTau, 3) / (turb::kappa * y);
}

} // namespace cf
