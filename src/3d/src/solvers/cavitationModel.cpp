#include "cf/fvm/CavitationModel.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

CavitationModel::Config CavitationModel::Config::fromCase(const Dictionary& caseDict)
{
    Config c;
    if (!caseDict.has("cavitation")) cfFatal("case has no cavitation{} block");
    const Dictionary& d = caseDict.getDict("cavitation");
    c.pSat = d.getScalar("pSat");
    c.rhoL = d.getScalarOr("rhoL", c.rhoL);
    c.rhoV = d.getScalarOr("rhoV", c.rhoV);
    c.Ce = d.getScalarOr("Ce", c.Ce);
    c.Cc = d.getScalarOr("Cc", c.Cc);
    c.Rb = d.getScalarOr("Rb", c.Rb);
    c.alpha0 = d.getScalarOr("alpha0", c.alpha0);
    c.diffusivity = d.getScalarOr("diffusivity", c.diffusivity);
    c.validate();
    return c;
}

void CavitationModel::Config::validate() const
{
    if (!(rhoL > rhoV) || !(rhoV > 0)) cfFatal("cavitation: need rhoL > rhoV > 0");
    if (!(Ce > 0) || !(Cc > 0)) cfFatal("cavitation: Ce and Cc must be positive");
    if (!(Rb > 0)) cfFatal("cavitation: Rb must be positive");
    if (!(alpha0 >= 0 && alpha0 <= 1)) cfFatal("cavitation: alpha0 must be in [0, 1]");
    if (!(diffusivity >= 0)) cfFatal("cavitation: diffusivity must be non-negative");
    if (!std::isfinite(pSat)) cfFatal("cavitation: pSat must be finite");
}

CavitationModel::CavitationModel(const polyMesh& mesh, const surfaceScalarField& phi,
                                 const Schemes& schemes, const Config& cfg)
    : mesh_(mesh), phi_(phi), schemes_(schemes), cfg_(cfg),
      alpha_(mesh, "alpha_v", cfg.alpha0)
{
    for (const polyMesh::Patch& p : mesh_.patches())
    {
        const label patch = mesh_.patchIdByName(p.name);
        alpha_.setPatchValue(patch, p.type == "inlet" ? "fixedValue" : "zeroGradient",
                             cfg_.alpha0);
    }
}

void CavitationModel::correct(scalar dt, const Field<scalar>& pCell)
{
    const label nC = mesh_.nCells();
    const auto& vol = mesh_.cellVolumes();
    const Field<scalar> old = alpha_.internalField();

    fvMatrix eqn(mesh_, viewOf(alpha_), "alphaEqn");
    eqn.schemes() = schemes_;
    eqn.schemes().relaxation = 1.0;
    if (dt > 0) fvm::ddt(eqn, old, dt);
    fvm::div(eqn, phi_.values(), schemes_.div);
    if (cfg_.diffusivity > 0)
    {
        Field<scalar> d(static_cast<std::size_t>(nC), cfg_.diffusivity);
        fvm::laplacianFace(eqn, harmonicFaceGamma(mesh_, d), schemes_.laplacian, -1);
    }

    // c(p) alpha (1 - alpha_old) moved to the lhs: the diagonal gains
    // -vol c (1 - alpha_old), positive wherever the vapour collapses
    Field<scalar> sp(static_cast<std::size_t>(nC), 0.0);
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar dp = cfg_.pSat - pCell[u];
        const scalar ce = cfg_.Ce * std::sqrt(2.0 * std::max(dp, 0.0) / cfg_.rhoL);
        const scalar cc = cfg_.Cc * std::sqrt(2.0 * std::max(-dp, 0.0) / cfg_.rhoL);
        const scalar rate = 3.0 * (ce - cc) / cfg_.Rb;
        sp[u] = -vol[u] * rate * (1.0 - old[u]);
    }
    fvm::sp(eqn, sp, 1);

    eqn.solve(fvMatrix::Controls{linearTol, linearMaxIter, 2});

    const Field<scalar>& x = eqn.solution()[0];
    auto& a = alpha_.internalField();
    for (std::size_t u = 0; u < a.size(); ++u)
    {
        a[u] = std::min(std::max(x[u], 0.0), 1.0);
    }
}

Field<scalar> CavitationModel::mixtureDensity() const
{
    Field<scalar> rho(alpha_.internalField().size());
    for (std::size_t u = 0; u < rho.size(); ++u)
    {
        const scalar a = alpha_.internalField()[u];
        rho[u] = a * cfg_.rhoV + (1.0 - a) * cfg_.rhoL;
    }
    return rho;
}

scalar CavitationModel::alphaMax() const
{
    scalar m = 0;
    for (const scalar a : alpha_.internalField()) m = std::max(m, a);
    return m;
}

scalar CavitationModel::alphaMin() const
{
    scalar m = 1;
    for (const scalar a : alpha_.internalField()) m = std::min(m, a);
    return m;
}

} // namespace cf
