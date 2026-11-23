#include "cf/turbulence/DesModel.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{

DesModel::DesModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                   scalar nu, const Schemes& schemes, const TurbulenceConfig& cfg)
    : KEpsilon(mesh, U, phi, nu, schemes, cfg),
      length_(static_cast<std::size_t>(mesh.nCells()), 0.0)
{
    cDES_ = cfg.cDES;
    if (!(cDES_ > 0)) cfFatal("des: C_DES must be positive");
}

void DesModel::updateNut()
{
    // The RANS transport has run; replace the dissipation with the shielded
    // one and rebuild nut from it. cellWallDistance() is the same shortest-
    // path field the SST blend uses.
    const Field<scalar>& d = cellWallDistance();
    const auto& vol = mesh_.cellVolumes();
    const std::size_t n = nut_.size();
    for (std::size_t c = 0; c < n; ++c)
    {
        const scalar delta = std::cbrt(vol[c]);
        const scalar ell = std::min(turb::kappa * d[c], cDES_ * delta);
        length_[c] = ell;
        const scalar k = std::max(k_.internalField()[c], turb::kFloor);
        const scalar epsShielded = k * std::sqrt(k) / ell;
        const scalar eps = std::max(epsShielded, turb::dissFloor);
        diss_.internalField()[c] = eps;
        nut_[c] = KEpsilon::Cmu * k * k / eps;
    }
}

} // namespace cf
