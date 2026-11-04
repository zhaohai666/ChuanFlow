#pragma once

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/EnergySolver.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/VtkWriter.h"
#include "cf/turbulence/TurbulenceModel.h"

namespace cf
{
namespace app
{

namespace energy_detail
{
/// Reports carry the heat balance, and a checker asserts it to well below the
/// six digits a plain stream write gives.
inline std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
} // namespace energy_detail

/// The thermal field of a case: the model plus the patches that carry a heat
/// condition, so the balance report knows which ones to add up.
struct CaseEnergy
{
    const polyMesh& mesh;
    std::unique_ptr<EnergyModel> model;
    std::vector<std::string> thermalPatches;

    CaseEnergy(const polyMesh& m, std::unique_ptr<EnergyModel> e,
               std::vector<std::string> names)
        : mesh(m), model(std::move(e)), thermalPatches(std::move(names))
    {
    }

    const volScalarField& temperature() const { return model->temperature(); }

    /// Temperature range, the heat rate through every patch that carries a
    /// condition, and the residual of the discrete balance
    /// heat in + generated = enthalpy out. For a steady case that residual is
    /// zero to the linear-solver tolerance; for a transient one it is the rate
    /// the fluid stores.
    void logBalance() const
    {
        log::info("energy: Tmin ", energy_detail::num(model->temperature().min()),
                  " Tmax ", energy_detail::num(model->temperature().max()));
        scalar in = 0;
        for (const std::string& name : thermalPatches)
        {
            const label pid = mesh.patchIdByName(name);
            const scalar q = model->patchHeatRate(pid);
            log::info("heatPatch ", name, " Q ", energy_detail::num(q));
            const auto wr = model->wallLayerReport(pid);
            if (wr.active)
            {
                log::info("thermalWall ", name, " uTau ", energy_detail::num(wr.uTau),
                          " yPlus ", energy_detail::num(wr.yPlus), " tPlus ",
                          energy_detail::num(wr.tPlus), " q ", energy_detail::num(wr.heatRate));
            }
            in += q;
        }
        const scalar src = model->volumetricHeatRate();
        const scalar out = model->netEnthalpyRate();
        log::info("energy: in ", energy_detail::num(in), " source ",
                  energy_detail::num(src), " enthalpyOut ", energy_detail::num(out),
                  " balance ", energy_detail::num(in + src - out));
    }

    void write(VtkWriter& writer) const
    {
        writer.addCellScalar("T", model->temperature().internalField());
    }
};

/// The energy{} block of a case, or a null model when the case does not ask for
/// a thermal field. The model owns T and reads the flux - and, through `turb`,
/// the eddy viscosity, strain rate and turbulent dissipation - on every
/// advance, so build it once the fields it reads exist and keep it alive for
/// as long as the solver runs.
///
/// Thermal patch conditions come from the same patches{} block as the flow
/// ones:
///   patches
///   {
///       hotWall  { T 340; }                 // fixedValue temperature
///       fluxWall { heatFlux 1.0e3; }         // W/m2 into the fluid
///       lid      { TType zeroGradient; }     // adiabatic, or any scalar BC
///   }
/// A patch with none of those keys keeps the scalar default, which is
/// adiabatic, so a case only lists the patches that carry heat.
inline CaseEnergy makeEnergy(const Dictionary& caseDict, const polyMesh& mesh,
                             const surfaceScalarField& phi, scalar nu,
                             const Schemes& schemes, RansModel* turb)
{
    std::vector<std::string> names;
    if (!caseDict.has("energy")) return CaseEnergy{mesh, nullptr, std::move(names)};

    const EnergyProps props = EnergyProps::fromCase(caseDict);
    EnergyWallLayer wl;
    if (turb)
    {
        wl.faces = &turb->wallFaces();
        wl.uTau = &turb->wallUTau();
        wl.yPlus = &turb->wallYPlus();
        wl.wallDist = &turb->wallDistance();
    }
    std::unique_ptr<EnergyModel> model = std::make_unique<EnergyModel>(
        mesh, phi, schemes, props, nu, turb ? &turb->nut() : nullptr,
        turb ? &turb->strainRate() : nullptr,
        turb ? &turb->turbDissipation() : nullptr, turb ? &wl : nullptr);

    if (caseDict.has("patches"))
    {
        const Dictionary& ps = caseDict.getDict("patches");
        volScalarField& t = model->temperature();
        for (const auto& name : ps.keys())
        {
            const Dictionary& pd = ps.getDict(name);
            const label patch = mesh.patchIdByName(name);
            if (pd.has("T"))
            {
                t.setPatchValue(patch, pd.getWordOr("TType", "fixedValue"), pd.getScalar("T"));
            }
            else if (pd.has("heatFlux"))
            {
                // Fourier's law: the gradient along the outward normal that
                // sends heatFlux W/m2 into the fluid
                t.setPatchValue(patch, "fixedGradient", pd.getScalar("heatFlux") / props.kappa);
            }
            else if (pd.has("TType"))
            {
                t.setPatchValue(patch, pd.getWord("TType"), pd.getScalarOr("TValue", 0.0));
            }
            else
            {
                continue;
            }
            names.push_back(name);
        }
    }

    // after the boundary conditions: the starting mode has to feel like the
    // field the case will relax from
    model->initialize();

    log::info("energy: rho=", props.rho, " cp=", props.cp, " kappa=", props.kappa,
              " alpha=", props.alphaMolecular(), " expansion=", props.beta,
              " dissipation=", props.viscousDissipation ? "on" : "off",
              " turbulentDissipation=", props.turbulentDissipation ? "on" : "off",
              " thermalPatches=", names.size());
    return CaseEnergy{mesh, std::move(model), std::move(names)};
}

} // namespace app
} // namespace cf
