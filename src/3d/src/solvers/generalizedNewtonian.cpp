#include "cf/fvm/GeneralizedNewtonian.h"

#include <algorithm>
#include <cmath>
#include <functional>

#include "cf/common/Error.h"

namespace cf
{
namespace
{

Field<scalar> component(const Field<Vector>& u, int comp)
{
    Field<scalar> out(static_cast<std::size_t>(u.size()));
    for (std::size_t i = 0; i < u.size(); ++i) out[i] = u[i][comp];
    return out;
}

bool isWallPatch(const polyMesh::Patch& p) { return p.type == "wall"; }

} // namespace

GeneralizedNewtonianViscosity::GeneralizedNewtonianViscosity(
    const polyMesh& mesh, volVectorField& U, scalar rho, scalar nuNewtonian,
    const Rheology& rheology, const Schemes& schemes,
    const GeneralizedNewtonianControls& ctl)
    : mesh_(mesh),
      U_(U),
      rho_(rho),
      nuNewtonian_(nuNewtonian),
      muNewtonian_(nuNewtonian * rho),
      rheology_(rheology),
      schemes_(schemes),
      ctl_(ctl),
      strain_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      nuCell_(static_cast<std::size_t>(mesh.nCells()), nuNewtonian),
      nuEffFace_(static_cast<std::size_t>(mesh.nFaces()), nuNewtonian)
{
    if (!(rho_ > 0) || !(nuNewtonian_ > 0))
        cfFatal("generalizedNewtonian: rho and the Newtonian nu must be positive");
    if (!(ctl_.shearFloor > 0) || !(ctl_.regularization > 0))
        cfFatal("generalizedNewtonian: shearFloor and regularization must be positive");
    if (!(ctl_.relaxation > 0) || ctl_.relaxation > 1.0)
        cfFatal("generalizedNewtonian: relaxation must lie in (0, 1]");
    if (rheology_.flowIndex < 0 || rheology_.yieldStress < 0)
        cfFatal("generalizedNewtonian: a negative flow index or yield stress"
                " has no constitutive reading");
    if (rheology_.flowIndex > 0 && !(rheology_.consistency > 0))
        cfFatal("generalizedNewtonian: a flow index needs a positive"
                " consistency K (Pa.s^n)");

    for (const polyMesh::Patch& p : mesh_.patches())
    {
        if (!isWallPatch(p)) continue;
        for (label f = p.start; f < p.start + p.size; ++f)
            wallFaces_.push_back(f);
    }
    wallDist_.assign(wallFaces_.size(), 0.0);
    wallNu_.assign(wallFaces_.size(), nuNewtonian_);
    wallGdot_.assign(wallFaces_.size(), 0.0);
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const std::size_t f = static_cast<std::size_t>(wallFaces_[i]);
        const Vector& sf = mesh_.faceSf()[f];
        const scalar area = mag(sf);
        if (!(area > 0))
            cfFatal("generalizedNewtonian: zero-area wall face");
        // Normal distance from the owner cell centre to the wall face: exactly
        // the denominator the discrete diffusion term divides U_c by.
        wallDist_[i] = dot(mesh_.deltaCoeffVec()[f], sf / area);
        if (!(wallDist_[i] > 0))
            cfFatal("generalizedNewtonian: wall face not outside its owner cell");
    }
}

scalar GeneralizedNewtonianViscosity::apparentNu(scalar gdot) const
{
    return apparentViscosityPaS(rheology_, gdot, muNewtonian_,
                                ctl_.shearFloor, ctl_.regularization) / rho_;
}

void GeneralizedNewtonianViscosity::computeStrainRate()
{
    const FieldView view = viewOf(U_);
    Field<scalar> u[3];
    for (int comp = 0; comp < 3; ++comp)
        u[comp] = component(U_.internalField(), comp);

    // gradU[c][s] = d U_c / d x_s; the symmetric part's magnitude
    // sqrt(2 S:S) is the shear rate a scalar constitutive law is written in.
    Field<Vector> gradU[3];
    for (int comp = 0; comp < 3; ++comp)
        gradU[comp] = fvc::grad(mesh_, u[comp],
                               makeFaceValueFn(mesh_, view, comp), schemes_.grad);

    for (std::size_t c = 0; c < strain_.size(); ++c)
    {
        scalar s2 = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                const scalar sij = 0.5 * (gradU[i][c][j] + gradU[j][c][i]);
                s2 += sij * sij;
            }
        strain_[c] = std::sqrt(2.0 * s2);
    }
}

void GeneralizedNewtonianViscosity::correctTurbulence(scalar dt)
{
    (void)dt;   // an algebraic law, not a transport equation
    computeStrainRate();

    const scalar r = ctl_.relaxation;
    for (std::size_t c = 0; c < nuCell_.size(); ++c)
    {
        const scalar fresh = apparentNu(strain_[c]);
        nuCell_[c] = r * fresh + (1.0 - r) * nuCell_[c];
        if (!(nuCell_[c] > 0) || !std::isfinite(nuCell_[c]))
            cfFatal("generalizedNewtonian: non-finite effective viscosity");
    }
    nuEffFace_ = harmonicFaceGamma(mesh_, nuCell_);

    // The wall face is where the discrete momentum equation charges shear:
    // tau_w = rho * nu_w * |U_c,parallel| / y_wall. Overwriting the
    // extrapolated cell viscosity with the value at THAT shear rate is what
    // makes a shear-thinning channel resist less than a plug would suggest.
    const auto& own = mesh_.owner();
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const std::size_t f = static_cast<std::size_t>(wallFaces_[i]);
        const Vector& sf = mesh_.faceSf()[f];
        const Vector n = sf / mag(sf);
        const Vector& up = U_.internalField()[static_cast<std::size_t>(own[f])];
        const scalar parallel = mag(up - n * dot(up, n));
        const scalar gdot = parallel / wallDist_[i];
        const scalar fresh = apparentNu(gdot);
        wallNu_[i] = r * fresh + (1.0 - r) * wallNu_[i];
        wallGdot_[i] = gdot;
        nuEffFace_[f] = wallNu_[i];
    }
}

scalar GeneralizedNewtonianViscosity::maxShearRate() const
{
    scalar v = 0;
    for (const scalar s : strain_) v = std::max(v, s);
    for (const scalar s : wallGdot_) v = std::max(v, s);
    return v;
}

scalar GeneralizedNewtonianViscosity::minViscosity() const
{
    scalar v = nuCell_.empty() ? nuNewtonian_ : nuCell_[0];
    for (const scalar n : nuCell_) v = std::min(v, n);
    return v;
}

scalar GeneralizedNewtonianViscosity::maxViscosity() const
{
    scalar v = nuCell_.empty() ? nuNewtonian_ : nuCell_[0];
    for (const scalar n : nuCell_) v = std::max(v, n);
    return v;
}

scalar GeneralizedNewtonianViscosity::meanWallShearPa() const
{
    scalar weighted = 0, area = 0;
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const std::size_t f = static_cast<std::size_t>(wallFaces_[i]);
        const scalar a = mag(mesh_.faceSf()[f]);
        weighted += a * rho_ * wallNu_[i] * wallGdot_[i];
        area += a;
    }
    return area > 0 ? weighted / area : 0.0;
}

scalar GeneralizedNewtonianViscosity::meanWallViscosity() const
{
    scalar weighted = 0, area = 0;
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const scalar a = mag(mesh_.faceSf()[static_cast<std::size_t>(wallFaces_[i])]);
        weighted += a * wallNu_[i];
        area += a;
    }
    return area > 0 ? weighted / area : 0.0;
}

scalar GeneralizedNewtonianViscosity::meanWallShearRate() const
{
    scalar weighted = 0, area = 0;
    for (std::size_t i = 0; i < wallFaces_.size(); ++i)
    {
        const scalar a = mag(mesh_.faceSf()[static_cast<std::size_t>(wallFaces_[i])]);
        weighted += a * wallGdot_[i];
        area += a;
    }
    return area > 0 ? weighted / area : 0.0;
}

} // namespace cf
