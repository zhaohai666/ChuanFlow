// Phase-change melting / solidification family (F12, enthalpy-porosity).
//
// Replicates src/3d/apps/meltingApp.cpp's driver on the kernel side: PISO
// momentum steps with the Darcy porous damping read from the current liquid
// fraction, an EnergyModel that carries the latent heat as an enthalpy-based
// apparent capacity, and the same report the CLI app prints (time, f_l=0.5
// front, enthalpy offset, cumulative boundary heat, liquid/solid speed bounds,
// f_l min/max).
//
// The per-cell latent-heat bookkeeping lives in the kernel EnergyModel; this
// file only drives it and copies the state out to NumPy.
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/EnergySolver.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PisoSolver.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct MeltingReport
{
    double t = 0;
    double yF = 0;
    double H = 0;
    double Qin = 0;
    double uLiquid = 0;
    double uSolid = 0;
    double flMin = 1;
    double flMax = 0;
};

struct MeltingResult
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    std::unique_ptr<cf::EnergyModel> energy;
    std::vector<MeltingReport> reports;
    int steps = 0;
};

cf::PisoControls readPisoControls(const cf::Dictionary& caseDict,
                                  double& endTime, int& reportEvery)
{
    cf::PisoControls ctl;
    endTime = 1.0;
    reportEvery = 10;
    if (!caseDict.has("piso"))
        throw std::runtime_error("melting: case needs a piso{ dt; endTime; } block");
    const cf::Dictionary& sd = caseDict.getDict("piso");
    ctl.dt = sd.getScalar("dt");
    endTime = sd.getScalar("endTime");
    ctl.correctors = static_cast<int>(sd.getScalarOr("correctors", ctl.correctors));
    ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
    ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
    ctl.linearMaxIter = static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
    ctl.nonOrthCorrectors = static_cast<int>(
        sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
    reportEvery = static_cast<int>(sd.getScalarOr("reportEvery", reportEvery));
    if (sd.has("force"))
    {
        const auto f = sd.getScalarList("force");
        ctl.bodyForce = cf::Vector{f[0], f[1], f[2]};
    }
    return ctl;
}

void applyFlowPatches(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                       cf::volVectorField& U, cf::volScalarField& p)
{
    if (!caseDict.has("patches")) return;
    const cf::Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const cf::label patch = mesh.patchIdByName(name);
        const cf::Dictionary& pd = ps.getDict(name);
        if (pd.has("U"))
        {
            const auto ul = pd.getScalarList("U");
            const std::string type = pd.getWordOr("UType", "fixedValue");
            U.setPatchValue(patch, type, cf::Vector{ul[0], ul[1], ul[2]});
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

std::unique_ptr<MeltingResult> solveMeltingCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);

    auto out = std::make_unique<MeltingResult>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    const double nu = caseDict.getScalarOr("nu", -1.0);
    if (!(nu > 0))
        throw std::runtime_error("melting: case needs a positive top-level nu");

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);
    applyFlowPatches(caseDict, *out->mesh, *out->U, *out->p);

    cf::Schemes schemes = cf::loadSchemes(caseDict);
    schemes.ddt = "Euler";

    if (!caseDict.has("energy"))
        throw std::runtime_error("melting: case needs an energy{} block");
    const cf::EnergyProps props = cf::EnergyProps::fromCase(caseDict);
    if (!props.phase.enabled)
        throw std::runtime_error("melting: energy{} needs a phaseChange sub-block");

    phi.correctBoundaryFlux(*out->U);
    out->energy = std::make_unique<cf::EnergyModel>(
        *out->mesh, phi, schemes, props, nu);

    // Thermal BCs share the patches{} block with the flow ones.
    if (caseDict.has("patches"))
    {
        const cf::Dictionary& ps = caseDict.getDict("patches");
        cf::volScalarField& t = out->energy->temperature();
        for (const std::string& name : ps.keys())
        {
            const cf::label patch = out->mesh->patchIdByName(name);
            const cf::Dictionary& pd = ps.getDict(name);
            if (pd.has("T"))
                t.setPatchValue(patch, pd.getWordOr("TType", "fixedValue"),
                                 pd.getScalar("T"));
            else if (pd.has("heatFlux"))
                t.setPatchValue(patch, "fixedGradient",
                                pd.getScalar("heatFlux") / props.kappa);
            else if (pd.has("TType"))
                t.setPatchValue(patch, pd.getWord("TType"),
                                pd.getScalarOr("TValue", 0.0));
        }
    }
    out->energy->initialize();

    double endTime = 1.0;
    int reportEvery = 10;
    cf::PisoControls ctl = readPisoControls(caseDict, endTime, reportEvery);

    const std::size_t nC = static_cast<std::size_t>(out->mesh->nCells());
    cf::Field<cf::Vector> darcy(nC, cf::Vector::zero());
    ctl.cellForce = &darcy;

    std::vector<std::size_t> order(nC);
    for (std::size_t c = 0; c < nC; ++c) order[c] = c;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
              { return out->mesh->cellCentres()[a].y <
                       out->mesh->cellCentres()[b].y; });

    const double A = props.phase.mushyConstant;
    const double eps = props.phase.mushyEpsilon;

    double t = 0;
    double qCum = 0;
    double qPrev = 0;
    const double h0 = out->energy->totalEnthalpy();

    auto patchHeat = [&]()
    {
        double q = 0;
        for (cf::label pt = 0; pt < out->mesh->nPatches(); ++pt)
            q += out->energy->patchHeatRate(pt);
        return q;
    };

    auto report = [&]()
    {
        const cf::Field<cf::scalar>& fl = out->energy->liquidFraction();
        MeltingReport row;
        row.t = t;
        // front: f_l crosses 0.5 upward along y
        row.yF = out->mesh->minX().y;
        for (std::size_t i = 1; i < nC; ++i)
        {
            const std::size_t a = order[i - 1], b = order[i];
            if (fl[a] < 0.5 && fl[b] >= 0.5)
            {
                const double ya = out->mesh->cellCentres()[a].y;
                const double yb = out->mesh->cellCentres()[b].y;
                row.yF = ya + 0.5 * (yb - ya);
                break;
            }
        }
        double uLiquid = 0, uSolid = 0, flMin = 1, flMax = 0;
        for (std::size_t c = 0; c < nC; ++c)
        {
            flMin = std::min(flMin, fl[c]);
            flMax = std::max(flMax, fl[c]);
            const double um = cf::mag(out->U->internalField()[c]);
            if (fl[c] > 0.9) uLiquid = std::max(uLiquid, um);
            if (fl[c] < 0.1) uSolid = std::max(uSolid, um);
        }
        row.H = out->energy->totalEnthalpy() - h0;
        row.Qin = qCum;
        row.uLiquid = uLiquid;
        row.uSolid = uSolid;
        row.flMin = flMin;
        row.flMax = flMax;
        out->reports.push_back(row);
    };

    qPrev = patchHeat();
    report();
    int step = 0;
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const double dt = std::min(ctl.dt, endTime - t);
        const cf::Field<cf::scalar>& fl = out->energy->liquidFraction();
        const double cap = 0.5 / dt;
        for (std::size_t c = 0; c < nC; ++c)
        {
            const double damp = std::min(
                A * (1.0 - fl[c]) * (1.0 - fl[c]) / (fl[c] * fl[c] * fl[c] + eps),
                cap);
            darcy[c] = out->U->internalField()[c] * (-damp);
        }
        cf::PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        solveTransientIncompressiblePiso(*out->mesh, *out->U, *out->p, phi, nu,
                                         schemes, stepCtl, nullptr,
                                         out->energy.get());
        t += dt;
        ++step;

        const double qNow = patchHeat();
        qCum += 0.5 * dt * (qPrev + qNow);
        qPrev = qNow;
        if (step % reportEvery == 0) report();
    }
    report();
    out->steps = step;
    return out;
}

} // namespace

void registerMeltingBindings(py::module_& m)
{
    py::class_<MeltingReport>(m, "MeltingReport",
        "One reported state of the melting run: time, f_l=0.5 front position,"
        " enthalpy offset, cumulative boundary heat, liquid/solid max speed and"
        " liquid-fraction bounds.")
        .def_readonly("t", &MeltingReport::t)
        .def_readonly("yF", &MeltingReport::yF)
        .def_readonly("H", &MeltingReport::H)
        .def_readonly("Qin", &MeltingReport::Qin)
        .def_readonly("uLiquid", &MeltingReport::uLiquid)
        .def_readonly("uSolid", &MeltingReport::uSolid)
        .def_readonly("flMin", &MeltingReport::flMin)
        .def_readonly("flMax", &MeltingReport::flMax);

    py::class_<MeltingResult, std::unique_ptr<MeltingResult>>(
        m, "MeltingResult",
        "Result of one enthalpy-porosity melting run. `reports` holds the"
        " series the CLI app prints; `T`/`liquid_fraction`/`U`/`p` are the"
        " final cell-centred fields.")
        .def_property_readonly("reports",
            [](const MeltingResult& r) { return r.reports; })
        .def_property_readonly("n_steps", [](const MeltingResult& r) { return r.steps; })
        .def("n_cells", [](const MeltingResult& r) { return r.mesh->nCells(); })
        .def("T", [](const MeltingResult& r)
             { return cfbind::scalarFieldToNumpy(r.energy->temperature().internalField()); },
             "Final temperature, shape (n,), float64 C-order copy.")
        .def("liquid_fraction", [](const MeltingResult& r)
             { return cfbind::scalarFieldToNumpy(r.energy->liquidFraction()); },
             "Final liquid fraction f_l, shape (n,), float64 C-order copy.")
        .def("U", [](const MeltingResult& r)
             { return cfbind::vectorFieldToNumpy(r.U->internalField()); },
             "Final velocity, shape (n, 3), float64 C-order copy.")
        .def("p", [](const MeltingResult& r)
             { return cfbind::scalarFieldToNumpy(r.p->internalField()); },
             "Final pressure, shape (n,), float64 C-order copy.")
        .def("cell_centres",
            [](const MeltingResult& r)
             { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("solve_melting_case", &solveMeltingCase,
          py::arg("case_dict_path"),
          "Read a phase-change melting case.dict and march the PISO + enthalpy"
          " energy coupling to endTime; returns the report series and final"
          " fields.");
}
