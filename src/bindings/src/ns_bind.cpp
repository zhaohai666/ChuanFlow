// Incompressible Navier-Stokes family binding (grid 1i): laminar transient
// PISO and steady SIMPLE pressure-velocity coupling, driven off the case dict
// the 3d apps read. Replicates the patch BCs and initial seed of
// pisoFoamApp / simpleFoamApp, but with the turbulence and energy hooks left
// null (laminar only; see src/bindings/README_ns.md for the contract and the
// deferred scope). The solver operators live in the kernel; this file drives
// them and copies the state out to NumPy / plain views.
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
#include "cf/fvm/PressureCorrection.h"
#include "cf/fvm/PisoSolver.h"
#include "cf/fvm/SimpleSolver.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct NsPisoStep
{
    double t = 0;
    double u_max = 0;
    double continuity = 0;
    int correctors = 0;
};

struct NsPisoCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    std::vector<NsPisoStep> history;
    int steps = 0;
    double worstContinuity = 0;
    double nu = 0;
    double endTime = 0;
    double dt = 0;
};

struct NsSimpleCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    cf::SimpleResult flow;
    double nu = 0;
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
            throw std::runtime_error("ns: patch '" + name + "' is not in the mesh");
        if (pd.has("U"))
        {
            const auto ul = pd.getScalarList("U");
            U.setPatchValue(patch, pd.getWordOr("UType", "fixedValue"),
                            cf::Vector{ul[0], ul[1], ul[2]});
        }
        else
        {
            // wall patches are no-slip unless the case says otherwise; the old
            // zeroGradient default silently turned them into free-slip walls.
            const bool isWall =
                mesh.patches()[static_cast<std::size_t>(patch)].type == "wall";
            U.setPatchValue(patch,
                            isWall ? "fixedValue"
                                   : pd.getWordOr("UType", "zeroGradient"),
                            cf::Vector::zero());
        }
        if (pd.has("p"))
            p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"),
                            pd.getScalar("p"));
        else
            p.setPatchValue(patch, "zeroGradient", 0.0);
    }
}

void seedInitialVelocity(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                         cf::volVectorField& U, cf::surfaceScalarField& phi)
{
    if (!caseDict.has("initial")) return;
    const cf::Dictionary& d = caseDict.getDict("initial");
    if (!d.has("U")) return;
    const auto v = d.getScalarList("U");
    const cf::Vector u0{v[0], v[1], v[2]};
    for (cf::Vector& u : U.internalField()) u = u0;
    cf::updatePhiInternal(mesh, U, phi);
    phi.correctBoundaryFlux(U);
}

cf::scalar readNu(const cf::Dictionary& caseDict)
{
    const cf::scalar nu = caseDict.getScalar("nu");
    if (!(nu > 0))
        throw std::runtime_error("ns: nu must be finite and positive");
    return nu;
}

std::unique_ptr<NsPisoCase> solvePisoNsCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<NsPisoCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));
    out->nu = readNu(caseDict);

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);

    applyPatches(caseDict, *out->mesh, *out->U, *out->p);
    seedInitialVelocity(caseDict, *out->mesh, *out->U, phi);

    cf::Schemes schemes = cf::loadSchemes(caseDict);
    schemes.ddt = "Euler";

    cf::PisoControls ctl;
    if (!caseDict.has("piso"))
        throw std::runtime_error("ns: case needs a piso { dt; endTime; ... } block");
    const cf::Dictionary& sd = caseDict.getDict("piso");
    ctl.dt = sd.getScalar("dt");
    out->endTime = sd.getScalar("endTime");
    ctl.correctors = static_cast<int>(sd.getScalarOr("correctors", ctl.correctors));
    ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
    ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
    ctl.linearMaxIter = static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
    ctl.nonOrthCorrectors =
        static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
    out->dt = ctl.dt;
    if (sd.has("force"))
    {
        const auto f = sd.getScalarList("force");
        ctl.bodyForce = cf::Vector{f[0], f[1], f[2]};
    }

    if (!(ctl.dt > 0))
        throw std::runtime_error("ns: piso dt must be finite and positive");
    if (!(out->endTime > 0))
        throw std::runtime_error("ns: piso endTime must be finite and positive");

    phi.correctBoundaryFlux(*out->U);

    cf::scalar t = 0;
    int reportEvery = 1;
    if (sd.has("reportEvery"))
        reportEvery = static_cast<int>(sd.getScalar("reportEvery"));

    while (t < out->endTime - 1e-12 * std::max(1.0, out->endTime))
    {
        cf::PisoControls stepCtl = ctl;
        stepCtl.dt = std::min(ctl.dt, out->endTime - t);
        const cf::PisoStepResult res =
            cf::solveTransientIncompressiblePiso(*out->mesh, *out->U, *out->p, phi,
                                                out->nu, schemes, stepCtl);
        t += stepCtl.dt;
        ++out->steps;

        cf::scalar uMax = 0;
        for (const cf::Vector& u : out->U->internalField())
            uMax = std::max(uMax, u.x);
        out->worstContinuity = std::max(out->worstContinuity, res.continuityError);

        if (reportEvery <= 0 || out->steps % reportEvery == 0
            || t >= out->endTime - 1e-12 * out->endTime)
        {
            NsPisoStep row;
            row.t = t;
            row.u_max = uMax;
            row.continuity = res.continuityError;
            row.correctors = res.correctorsUsed;
            out->history.push_back(row);
        }
    }

    if (!std::isfinite(out->worstContinuity))
        throw std::runtime_error("ns: non-finite continuity error during PISO solve");
    return out;
}

std::unique_ptr<NsSimpleCase> solveSimpleNsCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<NsSimpleCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));
    out->nu = readNu(caseDict);

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);

    applyPatches(caseDict, *out->mesh, *out->U, *out->p);
    seedInitialVelocity(caseDict, *out->mesh, *out->U, phi);

    const cf::Schemes schemes = cf::loadSchemes(caseDict);

    cf::SimpleControls ctl;
    if (!caseDict.has("simple"))
        throw std::runtime_error("ns: case needs a simple { ... } block");
    const cf::Dictionary& sd = caseDict.getDict("simple");
    ctl.maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", ctl.maxIterations));
    ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
    ctl.pcorrTol = sd.getScalarOr("pcorrTol", ctl.pcorrTol);
    ctl.momentumTol = sd.getScalarOr("momentumTol", ctl.momentumTol);
    ctl.relaxU = sd.getScalarOr("relaxU", ctl.relaxU);
    ctl.relaxP = sd.getScalarOr("relaxP", ctl.relaxP);
    ctl.pRelaxStart = sd.getScalarOr("pRelaxStart", ctl.pRelaxStart);
    ctl.nonOrthCorrectors =
        static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
    if (sd.has("force"))
    {
        const auto f = sd.getScalarList("force");
        ctl.bodyForce = cf::Vector{f[0], f[1], f[2]};
    }

    phi.correctBoundaryFlux(*out->U);
    out->flow = cf::solveIncompressibleSimple(*out->mesh, *out->U, *out->p, phi,
                                              out->nu, schemes, ctl);
    return out;
}

} // namespace

void registerNsBindings(py::module_& m)
{
    py::class_<NsPisoStep>(m, "NsPisoStep",
        "One reported PISO step: time, centreline velocity probe, continuity "
        "error after the last corrector and correctors used.")
        .def_readonly("t", &NsPisoStep::t)
        .def_readonly("u_max", &NsPisoStep::u_max)
        .def_readonly("continuity", &NsPisoStep::continuity)
        .def_readonly("correctors", &NsPisoStep::correctors);

    py::class_<NsPisoCase, std::unique_ptr<NsPisoCase>>(m, "NsPisoCase",
        "Result of one laminar transient PISO incompressible run: solved U/p, "
        "step history and worst continuity error.")
        .def_property_readonly("U", [](const NsPisoCase& r)
            { return cfbind::vectorFieldToNumpy(r.U->internalField()); })
        .def_property_readonly("p", [](const NsPisoCase& r)
            { return cfbind::scalarFieldToNumpy(r.p->internalField()); })
        .def_readonly("history", &NsPisoCase::history)
        .def_readonly("steps", &NsPisoCase::steps)
        .def_readonly("worst_continuity", &NsPisoCase::worstContinuity)
        .def_readonly("nu", &NsPisoCase::nu)
        .def_readonly("end_time", &NsPisoCase::endTime)
        .def_readonly("dt", &NsPisoCase::dt)
        .def("n_cells", [](const NsPisoCase& r) { return r.mesh->nCells(); })
        .def("cell_centres", [](const NsPisoCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    py::class_<NsSimpleCase, std::unique_ptr<NsSimpleCase>>(m, "NsSimpleCase",
        "Result of one laminar steady SIMPLE incompressible solve: solved U/p "
        "and the outer-loop stop signals.")
        .def_property_readonly("U", [](const NsSimpleCase& r)
            { return cfbind::vectorFieldToNumpy(r.U->internalField()); })
        .def_property_readonly("p", [](const NsSimpleCase& r)
            { return cfbind::scalarFieldToNumpy(r.p->internalField()); })
        .def_property_readonly("iterations",
            [](const NsSimpleCase& r) { return r.flow.iterations; })
        .def_property_readonly("continuity",
            [](const NsSimpleCase& r) { return r.flow.finalContinuityError; })
        .def_property_readonly("momentum_residual",
            [](const NsSimpleCase& r) { return r.flow.finalMomentumResidual; })
        .def_property_readonly("converged",
            [](const NsSimpleCase& r) { return r.flow.converged; })
        .def_readonly("nu", &NsSimpleCase::nu)
        .def("n_cells", [](const NsSimpleCase& r) { return r.mesh->nCells(); })
        .def("cell_centres", [](const NsSimpleCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("solve_piso_ns_case", &solvePisoNsCase, py::arg("case_dict_path"),
          "Read a laminar PISO incompressible case.dict, advance the transient "
          "solve and return the final U/p and step history.");
    m.def("solve_simple_ns_case", &solveSimpleNsCase, py::arg("case_dict_path"),
          "Read a laminar steady SIMPLE case.dict and return the solved U/p "
          "with the outer-loop stop signals.");
}
