// RSM family binding: second-moment (Reynolds stress) closure advanced by the
// transient incompressible PISO coupling.
//
// Replicates src/3d/apps/rsmApp.cpp's driver: mesh load, patch BCs, analytic
// profile seed, phi closure, RsmModel::Config, probe-cell selection, then the
// PISO transient loop with the model as eddy-viscosity hook. The six stress
// transports, epsilon equation and realizability clamping live in the kernel
// RsmModel; this file drives it and copies the state out to NumPy / plain views.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PisoSolver.h"
#include "cf/turbulence/RsmModel.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

using Vec3 = std::array<double, 3>;

struct RsmRow
{
    double t = 0;
    std::array<double, 6> stress{0, 0, 0, 0, 0, 0};
    double epsilon = 0;
    double k = 0;
    double nut = 0;
    double production = 0;
    double energy = 0;
};

struct RsmCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    std::unique_ptr<cf::RsmModel> rsm;
    cf::RsmModel::Config cfg;
    std::vector<RsmRow> history;
    std::size_t probe = 0;
    double endTime = 1.0;
    double dt = 0.01;
};

void applyPatches(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                  cf::volVectorField& U, cf::volScalarField& p)
{
    if (!caseDict.has("patches")) return;
    const cf::Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const cf::Dictionary& pd = ps.getDict(name);
        const cf::label patch = mesh.patchIdByName(name);
        if (patch < 0)
            throw std::runtime_error("rsm: patch '" + name + "' is not in the mesh");
        if (pd.has("U"))
        {
            const auto ul = pd.getScalarList("U");
            U.setPatchValue(patch, pd.getWordOr("UType", "fixedValue"),
                            cf::Vector{ul[0], ul[1], ul[2]});
        }
        else
        {
            U.setPatchValue(patch, pd.getWordOr("UType", "zeroGradient"),
                            cf::Vector::zero());
        }
        if (pd.has("p"))
            p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"),
                            pd.getScalar("p"));
        else
            p.setPatchValue(patch, "zeroGradient", 0.0);
    }
}

void seedProfile(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                 cf::volVectorField& U)
{
    if (!caseDict.has("profile")) return;
    const cf::Dictionary& pd = caseDict.getDict("profile");
    const std::string type = pd.getWord("type");
    const cf::scalar rate = pd.getScalarOr("rate", 1.0);
    const cf::Vector lo = mesh.minX();
    for (std::size_t c = 0; c < U.internalField().size(); ++c)
    {
        const cf::Vector& x = mesh.cellCentres()[c];
        if (type == "linearShear")
            U.internalField()[c] = cf::Vector{rate * (x.y - lo.y), 0, 0};
        else if (type == "zero")
            U.internalField()[c] = cf::Vector::zero();
        else
            throw std::runtime_error("profile: unknown type '" + type + "'");
    }
}

RsmRow reportRow(const cf::polyMesh& mesh, const cf::volVectorField& U,
                 const cf::RsmModel& rsm, std::size_t probe, cf::scalar t)
{
    RsmRow row;
    row.t = t;
    for (int c = 0; c < 6; ++c)
        row.stress[c] = rsm.stress(c).internalField()[probe];
    row.epsilon = rsm.epsilon()[probe];
    row.k = rsm.k()[probe];
    row.nut = rsm.nut()[probe];
    row.production = rsm.production()[probe];
    cf::scalar energy = 0, volSum = 0;
    for (std::size_t c = 0; c < U.internalField().size(); ++c)
    {
        energy += 0.5 * cf::mag2(U.internalField()[c]) * mesh.cellVolumes()[c];
        volSum += mesh.cellVolumes()[c];
    }
    row.energy = energy / volSum;
    return row;
}

std::unique_ptr<RsmCase> solveRsmCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<RsmCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    const cf::scalar nu = caseDict.getScalarOr("nu", -1.0);
    if (!(nu > 0))
        throw std::runtime_error("rsm: case needs a positive top-level nu");

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);

    applyPatches(caseDict, *out->mesh, *out->U, *out->p);
    seedProfile(caseDict, *out->mesh, *out->U);
    phi.correctBoundaryFlux(*out->U);

    cf::Schemes schemes = cf::loadSchemes(caseDict);
    schemes.ddt = "Euler";

    out->cfg = cf::RsmModel::Config::fromCase(caseDict);
    out->rsm = std::make_unique<cf::RsmModel>(*out->mesh, *out->U, phi, nu,
                                             schemes, out->cfg);

    cf::Vector target = 0.5 * (*out->mesh).minX() + 0.5 * (*out->mesh).maxX();
    if (caseDict.has("probe"))
    {
        const auto pl = caseDict.getDict("probe").getScalarList("xyz");
        target = cf::Vector{pl[0], pl[1], pl[2]};
    }
    cf::scalar best = 1e300;
    for (std::size_t c = 0; c < out->mesh->cellCentres().size(); ++c)
    {
        const cf::scalar d2 = cf::mag2(out->mesh->cellCentres()[c] - target);
        if (d2 < best) { best = d2; out->probe = c; }
    }

    cf::PisoControls ctl;
    cf::scalar endTime = 1.0;
    int reportEvery = 10;
    if (caseDict.has("piso"))
    {
        const cf::Dictionary& sd = caseDict.getDict("piso");
        ctl.dt = sd.getScalar("dt");
        endTime = sd.getScalar("endTime");
        ctl.correctors = static_cast<int>(sd.getScalarOr("correctors", ctl.correctors));
        ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
        ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
        ctl.linearMaxIter =
            static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
        ctl.nonOrthCorrectors =
            static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
        reportEvery = static_cast<int>(sd.getScalarOr("reportEvery", reportEvery));
        if (sd.has("force"))
        {
            const auto f = sd.getScalarList("force");
            ctl.bodyForce = cf::Vector{f[0], f[1], f[2]};
        }
    }
    out->endTime = endTime;
    out->dt = ctl.dt;

    out->history.push_back(reportRow(*out->mesh, *out->U, *out->rsm, out->probe, 0.0));

    cf::scalar t = 0;
    int step = 0;
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const cf::scalar dt = std::min(ctl.dt, endTime - t);
        cf::PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        cf::solveTransientIncompressiblePiso(*out->mesh, *out->U, *out->p, phi,
                                            nu, schemes, stepCtl, out->rsm.get());
        t += dt;
        ++step;
        if (step % reportEvery == 0)
            out->history.push_back(reportRow(*out->mesh, *out->U, *out->rsm,
                                             out->probe, t));
    }
    out->history.push_back(reportRow(*out->mesh, *out->U, *out->rsm, out->probe, t));
    return out;
}

} // namespace

void registerRsmBindings(py::module_& m)
{
    py::class_<RsmRow>(m, "RsmRow",
        "One reported RSM state at the probe cell: six Reynolds stresses "
        "(11,22,33,12,13,23), epsilon, k, nut, production and mean kinetic energy.")
        .def_readonly("t", &RsmRow::t)
        .def_readonly("stress", &RsmRow::stress)
        .def_readonly("epsilon", &RsmRow::epsilon)
        .def_readonly("k", &RsmRow::k)
        .def_readonly("nut", &RsmRow::nut)
        .def_readonly("production", &RsmRow::production)
        .def_readonly("energy", &RsmRow::energy);

    py::class_<RsmCase, std::unique_ptr<RsmCase>>(m, "RsmCase",
        "Result of one RSM run: solved velocity/pressure, six stress fields, "
        "epsilon, k, nut, production, the probe-cell history and the model config.")
        .def_property_readonly("U", [](const RsmCase& r)
            { return cfbind::vectorFieldToNumpy(r.U->internalField()); })
        .def_property_readonly("p", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.p->internalField()); })
        .def_property_readonly("stress_11", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(0).internalField()); })
        .def_property_readonly("stress_22", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(1).internalField()); })
        .def_property_readonly("stress_33", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(2).internalField()); })
        .def_property_readonly("stress_12", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(3).internalField()); })
        .def_property_readonly("stress_13", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(4).internalField()); })
        .def_property_readonly("stress_23", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->stress(5).internalField()); })
        .def_property_readonly("epsilon", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->epsilon()); })
        .def_property_readonly("k", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->k()); })
        .def_property_readonly("nut", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->nut()); })
        .def_property_readonly("production", [](const RsmCase& r)
            { return cfbind::scalarFieldToNumpy(r.rsm->production()); })
        .def_readonly("history", &RsmCase::history)
        .def_readonly("probe", &RsmCase::probe)
        .def_property_readonly("C1", [](const RsmCase& r) { return r.cfg.C1; })
        .def_property_readonly("C1eps", [](const RsmCase& r) { return r.cfg.C1eps; })
        .def_property_readonly("C2eps", [](const RsmCase& r) { return r.cfg.C2eps; })
        .def_property_readonly("k0", [](const RsmCase& r) { return r.cfg.k0; })
        .def_property_readonly("eps0", [](const RsmCase& r) { return r.cfg.eps0; })
        .def_property_readonly("end_time", [](const RsmCase& r) { return r.endTime; })
        .def_property_readonly("dt", [](const RsmCase& r) { return r.dt; })
        .def("n_cells", [](const RsmCase& r) { return r.mesh->nCells(); })
        .def("cell_centres", [](const RsmCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("solve_rsm_case", &solveRsmCase, py::arg("case_dict_path"),
          "Read an RSM case.dict and advance the transient incompressible PISO "
          "coupling with the second-moment Reynolds-stress closure.");
}
