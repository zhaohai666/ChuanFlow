#pragma once

#include <memory>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/VtkWriter.h"
#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{
namespace app
{

/// The turbulence{} block of a case, or nullptr when the case does not ask for
/// a model. The model holds references to U and phi, so the caller must keep it
/// alive for as long as the solver runs.
inline std::unique_ptr<RansModel> makeTurbulence(const Dictionary& caseDict,
                                                 const polyMesh& mesh,
                                                 volVectorField& U,
                                                 const surfaceScalarField& phi,
                                                 scalar nu,
                                                 const Schemes& schemes)
{
    if (!caseDict.has("turbulence")) return nullptr;
    const TurbulenceConfig cfg = TurbulenceConfig::fromCase(caseDict);
    std::unique_ptr<RansModel> model =
        RansModel::create(mesh, U, phi, nu, schemes, cfg);
    if (cfg.lengthScale > 0)
    {
        log::info("turbulence: ", model->name(), " intensity=", cfg.intensity,
                  " lengthScale=", cfg.lengthScale);
    }
    else
    {
        log::info("turbulence: ", model->name(), " intensity=", cfg.intensity,
                  " lengthScale=auto (7% of the box diagonal)");
    }
    return model;
}

/// u_tau and y+ of the wall layer the model last evaluated: the numbers that say
/// whether the mesh is in the wall-function band (30 < y+ < 300) or resolves the
/// sublayer (y+ ~ 1), and therefore whether the model is applicable.
inline void logWallReport(const RansModel& model)
{
    if (model.wallFaces().empty())
    {
        log::info("turbulence: no wall patch, y+ not reported");
        return;
    }
    log::info("turbulence: u_tau=", model.meanFrictionVelocity(),
              " y+ mean=", model.meanWallYPlus(), " max=", model.maxWallYPlus());
}

/// k, the dissipation variable and nut, so a case can be inspected in ParaView
/// and a checker can read the near-wall profile out of the vtu.
inline void writeTurbulence(VtkWriter& writer, const RansModel& model)
{
    writer.addCellScalar("nut", model.nut());
    if (model.dissipationName() == "none") return;
    writer.addCellScalar("k", model.k());
    writer.addCellScalar(model.dissipationName(), model.dissipation());
}

} // namespace app
} // namespace cf
