#include "cf/turbulence/LesModel.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "gradU.h"

namespace cf
{

LesModel::Config LesModel::Config::fromCase(const Dictionary& caseDict)
{
    Config c;
    if (!caseDict.has("les")) return c;
    const Dictionary& d = caseDict.getDict("les");
    c.model = d.getWordOr("model", c.model);
    c.Cs = d.getScalarOr("Cs", c.Cs);
    c.Cw = d.getScalarOr("Cw", c.Cw);
    c.validate();
    return c;
}

void LesModel::Config::validate() const
{
    if (model != "smagorinsky" && model != "wale")
    {
        cfFatal("les: unknown model '" + model + "' (smagorinsky | wale)");
    }
    if (!(Cs >= 0) || !std::isfinite(Cs)) cfFatal("les: Cs must be finite and non-negative");
    if (!(Cw >= 0) || !std::isfinite(Cw)) cfFatal("les: Cw must be finite and non-negative");
}

LesModel::LesModel(const polyMesh& mesh, volVectorField& U, const surfaceScalarField& phi,
                   scalar nu, const Schemes& schemes, const Config& cfg)
    : mesh_(mesh), U_(U), nu_(nu), schemes_(schemes), cfg_(cfg),
      nut_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      nuEffFace_(static_cast<std::size_t>(mesh.nFaces()), nu),
      strain_(static_cast<std::size_t>(mesh.nCells()), 0.0),
      delta_(static_cast<std::size_t>(mesh.nCells()), 0.0)
{
    if (!(nu_ > 0)) cfFatal("les: molecular viscosity must be positive");
    cfg_.validate();
    (void)phi; // algebraic model: no transport, the flux is not read
    const auto& vol = mesh_.cellVolumes();
    for (std::size_t c = 0; c < delta_.size(); ++c)
    {
        delta_[c] = std::cbrt(vol[c]);
    }
    refreshFaceNut();
}

void LesModel::correctTurbulence(scalar /*dt*/)
{
    Field<Vector> gradU[3];
    turb_detail::computeGradU(mesh_, U_, schemes_, gradU);

    const bool wale = cfg_.model == "wale";
    for (std::size_t c = 0; c < nut_.size(); ++c)
    {
        scalar S[3][3];
        scalar magS = 0;
        turb_detail::strainAt(gradU, c, S, magS);
        strain_[c] = magS;

        if (!wale)
        {
            nut_[c] = cfg_.Cs * cfg_.Cs * delta_[c] * delta_[c] * magS;
            continue;
        }

        // Sd_ij = (g^2 + (g^2)^T)/2 - tr/3, g_ij = dU_i/dx_j
        scalar g2[3][3] = {};
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    g2[i][j] += gradU[i][c][k] * gradU[k][c][j];
                }
            }
        }
        scalar sd[3][3];
        const scalar tr = (g2[0][0] + g2[1][1] + g2[2][2]) / 3.0;
        scalar sd2 = 0, s2 = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                sd[i][j] = 0.5 * (g2[i][j] + g2[j][i]);
                sd[i][j] -= (i == j ? tr : 0.0);
                sd2 += sd[i][j] * sd[i][j];
                s2 += S[i][j] * S[i][j];
            }
        }
        // (S:S) = |S|^2/2; the WALE denominator is written in S:S and Sd:Sd
        const scalar num = std::pow(sd2, 1.5);
        const scalar den = std::pow(s2, 2.5) + std::pow(sd2, 1.25);
        nut_[c] = den > 0 ? cfg_.Cw * cfg_.Cw * delta_[c] * delta_[c] * num / den : 0.0;
    }

    refreshFaceNut();
}

void LesModel::refreshFaceNut()
{
    const Field<scalar> nutFace = harmonicFaceGamma(mesh_, nut_);
    for (std::size_t f = 0; f < nuEffFace_.size(); ++f)
    {
        nuEffFace_[f] = nu_ + nutFace[f];
    }
}

} // namespace cf
