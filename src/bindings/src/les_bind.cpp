// LES family binding: algebraic Smagorinsky / WALE sub-grid viscosity advanced
// by the transient incompressible PISO coupling.
//
// Replicates src/3d/apps/lesApp.cpp's driver: mesh load, patch BCs, analytic
// profile seed, phi closure, LesModel::Config, one algebraic evaluation on the
// seeded field, then the PISO transient loop with the model as eddy-viscosity
// hook. The SGS closures and the stress/filter-width math live in the kernel
// LesModel; this file drives it and copies the state out to NumPy / plain views.
#include <algorithm>
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
#include "cf/turbulence/LesModel.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct LesRow
{
    double t = 0;
    double nutMin = 0;
    double nutMax = 0;
    double nutMean = 0;
    double strainMax = 0;
    double energy = 0;
};

struct LesCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    std::unique_ptr<cf::LesModel> les;
    cf::LesModel::Config cfg;
    std::vector<LesRow> history;
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
            throw std::runtime_error("les: patch '" + name + "' is not in the mesh");
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
    const cf::scalar amp = pd.getScalarOr("amplitude", 1.0);
    const cf::Vector lo = mesh.minX();
    const cf::Vector hi = mesh.maxX();
    const cf::scalar H = hi.y - lo.y;
    for (std::size_t c = 0; c < U.internalField().size(); ++c)
    {
        const cf::Vector& x = mesh.cellCentres()[c];
        if (type == "linearShear")
        {
            U.internalField()[c] = cf::Vector{rate * (x.y - lo.y), 0, 0};
        }
        else if (type == "cosineShear")
        {
            U.internalField()[c] =
                cf::Vector{amp * std::cos(M_PI * (x.y - lo.y) / H), 0, 0};
        }
        else if (type == "uniaxialStrain")
        {
            U.internalField()[c] = cf::Vector{rate * (x.x - lo.x),
                                              -0.5 * rate * (x.y - lo.y),
                                              -0.5 * rate * (x.z - lo.z)};
        }
        else
        {
            throw std::runtime_error("profile: unknown type '" + type + "'");
        }
    }
}

LesRow reportRow(const cf::polyMesh& mesh, const cf::volVectorField& U,
                 const cf::LesModel& les, cf::scalar t)
{
    const auto& nut = les.nut();
    const auto& s = les.strainRate();
    cf::scalar nutMin = 1e300, nutMax = 0, nutSum = 0, sMax = 0;
    cf::scalar energy = 0, volSum = 0;
    for (std::size_t c = 0; c < nut.size(); ++c)
    {
        nutMin = std::min(nutMin, nut[c]);
        nutMax = std::max(nutMax, nut[c]);
        nutSum += nut[c];
        sMax = std::max(sMax, s[c]);
        energy += 0.5 * cf::mag2(U.internalField()[c]) * mesh.cellVolumes()[c];
        volSum += mesh.cellVolumes()[c];
    }
    LesRow row;
    row.t = t;
    row.nutMin = nutMin;
    row.nutMax = nutMax;
    row.nutMean = nutSum / static_cast<cf::scalar>(nut.size());
    row.strainMax = sMax;
    row.energy = energy / volSum;
    return row;
}

std::unique_ptr<LesCase> solveLesCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<LesCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    const cf::scalar nu = caseDict.getScalarOr("nu", -1.0);
    if (!(nu > 0))
        throw std::runtime_error("les: case needs a positive top-level nu");

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);

    applyPatches(caseDict, *out->mesh, *out->U, *out->p);
    seedProfile(caseDict, *out->mesh, *out->U);
    phi.correctBoundaryFlux(*out->U);

    cf::Schemes schemes = cf::loadSchemes(caseDict);
    schemes.ddt = "Euler";

    out->cfg = cf::LesModel::Config::fromCase(caseDict);
    out->les = std::make_unique<cf::LesModel>(*out->mesh, *out->U, phi, nu,
                                              schemes, out->cfg);

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

    out->les->correctTurbulence(0.0);
    out->history.push_back(reportRow(*out->mesh, *out->U, *out->les, 0.0));

    cf::scalar t = 0;
    int step = 0;
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const cf::scalar dt = std::min(ctl.dt, endTime - t);
        cf::PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        cf::solveTransientIncompressiblePiso(*out->mesh, *out->U, *out->p, phi,
                                            nu, schemes, stepCtl, out->les.get());
        t += dt;
        ++step;
        if (step % reportEvery == 0)
            out->history.push_back(reportRow(*out->mesh, *out->U, *out->les, t));
    }
    out->history.push_back(reportRow(*out->mesh, *out->U, *out->les, t));
    return out;
}

} // namespace

void registerLesBindings(py::module_& m)
{
    py::class_<LesRow>(m, "LesRow",
        "One reported LES state: time, sub-grid viscosity min/max/mean, "
        "resolved strain-rate max and volume-averaged kinetic energy.")
        .def_readonly("t", &LesRow::t)
        .def_readonly("nut_min", &LesRow::nutMin)
        .def_readonly("nut_max", &LesRow::nutMax)
        .def_readonly("nut_mean", &LesRow::nutMean)
        .def_readonly("strain_max", &LesRow::strainMax)
        .def_readonly("energy", &LesRow::energy);

    py::class_<LesCase, std::unique_ptr<LesCase>>(m, "LesCase",
        "Result of one LES SGS run: solved velocity/pressure, per-cell nut, "
        "strain rate and filter width, plus the reported history.")
        .def_property_readonly("U", [](const LesCase& r)
            { return cfbind::vectorFieldToNumpy(r.U->internalField()); })
        .def_property_readonly("p", [](const LesCase& r)
            { return cfbind::scalarFieldToNumpy(r.p->internalField()); })
        .def_property_readonly("nut", [](const LesCase& r)
            { return cfbind::scalarFieldToNumpy(r.les->nut()); })
        .def_property_readonly("strain_rate", [](const LesCase& r)
            { return cfbind::scalarFieldToNumpy(r.les->strainRate()); })
        .def_property_readonly("filter_width", [](const LesCase& r)
            { return cfbind::scalarFieldToNumpy(r.les->filterWidth()); })
        .def_readonly("history", &LesCase::history)
        .def_property_readonly("model", [](const LesCase& r)
            { return r.les->name(); })
        .def_property_readonly("Cs", [](const LesCase& r) { return r.cfg.Cs; })
        .def_property_readonly("Cw", [](const LesCase& r) { return r.cfg.Cw; })
        .def_property_readonly("end_time", [](const LesCase& r) { return r.endTime; })
        .def_property_readonly("dt", [](const LesCase& r) { return r.dt; })
        .def("n_cells", [](const LesCase& r) { return r.mesh->nCells(); })
        .def("cell_centres", [](const LesCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("solve_les_case", &solveLesCase, py::arg("case_dict_path"),
          "Read a LES SGS case.dict (Smagorinsky or WALE) and advance the "
          "transient incompressible PISO coupling with the algebraic SGS model.");
}
