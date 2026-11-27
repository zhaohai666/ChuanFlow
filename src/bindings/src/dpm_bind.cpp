// Discrete-particle (DPM) family binding: a steady carrier solved once, then
// Lagrangian spheres tracked through the frozen field (one-way), with the
// two-way Picard sweep the app runs when the case declares `twoWay`.
//
// Replicates src/3d/apps/dpmApp.cpp's driver: carrier SIMPLE solve, patch BCs,
// optional initial velocity seed, then runDpm. The drag law, particle
// integration and impulse accumulation live in the kernel; this file drives
// them and copies the state out to NumPy / plain views.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Dpm.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PressureCorrection.h"
#include "cf/fvm/SimpleSolver.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

using Vec3 = std::array<double, 3>;

Vec3 arr(const cf::Vector& v) { return {v.x, v.y, v.z}; }

struct DpmParticleView
{
    int id = 0;
    int stream = 0;
    std::string status;
    Vec3 x0{};
    Vec3 u0{};
    Vec3 x{};
    Vec3 u{};
    double age = 0;
    int cell = -1;
    double diameter = 0;
    double rhoParticle = 0;
};

struct DpmSampleView
{
    int id = 0;
    int step = 0;
    double t = 0;
    Vec3 x{};
    Vec3 u{};
    Vec3 uFluid{};
    double reynolds = 0;
    double cd = 0;
    Vec3 dragAcceleration{};
};

struct DpmCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> U;
    std::unique_ptr<cf::volScalarField> p;
    cf::SimpleResult flow;
    cf::DpmConfig cfg;
    cf::DpmResult result;
    bool twoWay = false;
    int outerUsed = 0;
    double twoWayResidual = 0;
    std::vector<DpmParticleView> particles;
    std::vector<DpmSampleView> samples;
};

void applyCarrierPatches(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                         cf::volVectorField& U, cf::volScalarField& p)
{
    if (!caseDict.has("patches")) return;
    const cf::Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const cf::Dictionary& pd = ps.getDict(name);
        const cf::label patch = mesh.patchIdByName(name);
        if (patch < 0)
            throw std::runtime_error("dpm: patch '" + name + "' is not in the mesh");
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

cf::Field<cf::Vector> forceFromImpulse(const cf::polyMesh& mesh,
                                       const cf::DpmResult& result,
                                       cf::scalar rhoF)
{
    cf::Field<cf::Vector> f(static_cast<std::size_t>(mesh.nCells()), cf::Vector::zero());
    if (!(result.sourceWindow > 0)) return f;
    const auto& vol = mesh.cellVolumes();
    for (cf::label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        f[u] = -1.0 * result.cellDragImpulse[u]
               / (rhoF * vol[u] * result.sourceWindow);
    }
    return f;
}

std::unique_ptr<DpmCase> solveDpmCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<DpmCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    const cf::scalar nu = caseDict.getScalar("nu");
    if (!(nu > 0))
        throw std::runtime_error("dpm: nu must be positive");

    out->U = std::make_unique<cf::volVectorField>(*out->mesh, "U");
    out->p = std::make_unique<cf::volScalarField>(*out->mesh, "p");
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);

    applyCarrierPatches(caseDict, *out->mesh, *out->U, *out->p);
    seedInitialVelocity(caseDict, *out->mesh, *out->U, phi);

    const cf::Schemes schemes = cf::loadSchemes(caseDict);
    cf::SimpleControls ctl;
    if (caseDict.has("simple"))
    {
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
    }
    out->flow = cf::solveIncompressibleSimple(*out->mesh, *out->U, *out->p, phi,
                                              nu, schemes, ctl);
    if (!out->flow.converged)
        throw std::runtime_error("dpm: carrier flow did not converge");

    out->cfg = cf::readDpmConfig(caseDict);
    const cf::scalar mu = out->cfg.carrierRho * nu;
    out->twoWay = out->cfg.twoWay;

    if (out->twoWay)
    {
        cf::Field<cf::Vector> force(static_cast<std::size_t>(out->mesh->nCells()),
                                    cf::Vector::zero());
        cf::scalar lastResidual = std::numeric_limits<cf::scalar>::infinity();
        int usedOuter = 0;
        for (int outer = 1; outer <= out->cfg.outerIterations; ++outer)
        {
            cf::SimpleControls twoCtl = ctl;
            twoCtl.cellForce = &force;
            const cf::SimpleResult coupled =
                cf::solveIncompressibleSimple(*out->mesh, *out->U, *out->p, phi,
                                              nu, schemes, twoCtl);
            if (!coupled.converged)
                throw std::runtime_error("dpm: two-way carrier flow did not converge");
            const cf::DpmResult trial = cf::runDpm(*out->mesh, *out->U, out->cfg, mu);
            const cf::Field<cf::Vector> source =
                forceFromImpulse(*out->mesh, trial, out->cfg.carrierRho);
            cf::scalar delta = 0, scale = 0;
            for (cf::label c = 0; c < out->mesh->nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                delta = std::max(delta, static_cast<cf::scalar>(cf::mag(source[u] - force[u])));
                scale = std::max(scale, static_cast<cf::scalar>(cf::mag(source[u])));
            }
            lastResidual = delta / std::max(scale, 1.0e-30);
            for (cf::label c = 0; c < out->mesh->nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                force[u] = (1.0 - out->cfg.twoWayRelax) * force[u]
                           + out->cfg.twoWayRelax * source[u];
            }
            usedOuter = outer;
            if (lastResidual <= out->cfg.twoWayTol) break;
        }
        out->outerUsed = usedOuter;
        out->twoWayResidual = lastResidual;
        out->result = cf::runDpm(*out->mesh, *out->U, out->cfg, mu);
    }
    else
    {
        out->result = cf::runDpm(*out->mesh, *out->U, out->cfg, mu);
    }

    for (const cf::DpmParticle& p : out->result.particles)
    {
        DpmParticleView v;
        v.id = p.id;
        v.stream = p.stream;
        v.status = cf::dpmStatusName(p.status);
        v.x0 = arr(p.x0);
        v.u0 = arr(p.u0);
        v.x = arr(p.x);
        v.u = arr(p.u);
        v.age = p.age;
        v.cell = p.cell;
        v.diameter = p.diameter;
        v.rhoParticle = p.rhoParticle;
        out->particles.push_back(v);
    }
    for (const cf::DpmSample& s : out->result.samples)
    {
        DpmSampleView v;
        v.id = s.id;
        v.step = s.step;
        v.t = s.t;
        v.x = arr(s.x);
        v.u = arr(s.u);
        v.uFluid = arr(s.uFluid);
        v.reynolds = s.reynolds;
        v.cd = s.cd;
        v.dragAcceleration = arr(s.dragAcceleration);
        out->samples.push_back(v);
    }
    return out;
}

} // namespace

void registerDpmBindings(py::module_& m)
{
    py::class_<DpmParticleView>(m, "DpmParticle",
        "One launched particle: launch state, terminal state and status.")
        .def_readonly("id", &DpmParticleView::id)
        .def_readonly("stream", &DpmParticleView::stream)
        .def_readonly("status", &DpmParticleView::status)
        .def_readonly("x0", &DpmParticleView::x0)
        .def_readonly("u0", &DpmParticleView::u0)
        .def_readonly("x", &DpmParticleView::x)
        .def_readonly("u", &DpmParticleView::u)
        .def_readonly("age", &DpmParticleView::age)
        .def_readonly("cell", &DpmParticleView::cell)
        .def_readonly("diameter", &DpmParticleView::diameter)
        .def_readonly("rho_particle", &DpmParticleView::rhoParticle);

    py::class_<DpmSampleView>(m, "DpmSample",
        "One frozen-flow trajectory sample, before the velocity update.")
        .def_readonly("id", &DpmSampleView::id)
        .def_readonly("step", &DpmSampleView::step)
        .def_readonly("t", &DpmSampleView::t)
        .def_readonly("x", &DpmSampleView::x)
        .def_readonly("u", &DpmSampleView::u)
        .def_readonly("u_fluid", &DpmSampleView::uFluid)
        .def_readonly("reynolds", &DpmSampleView::reynolds)
        .def_readonly("cd", &DpmSampleView::cd)
        .def_readonly("drag_acceleration", &DpmSampleView::dragAcceleration);

    py::class_<DpmCase, std::unique_ptr<DpmCase>>(m, "DpmCase",
        "Result of one DPM run: solved carrier field, particles, trajectory "
        "samples and outlet/wall/expired counts.")
        .def_property_readonly("U", [](const DpmCase& r)
            { return cfbind::vectorFieldToNumpy(r.U->internalField()); })
        .def_property_readonly("p", [](const DpmCase& r)
            { return cfbind::scalarFieldToNumpy(r.p->internalField()); })
        .def_readonly("particles", &DpmCase::particles)
        .def_readonly("samples", &DpmCase::samples)
        .def_property_readonly("outlet_count",
            [](const DpmCase& r) { return r.result.outletCount; })
        .def_property_readonly("wall_count",
            [](const DpmCase& r) { return r.result.wallCount; })
        .def_property_readonly("expired_count",
            [](const DpmCase& r) { return r.result.expiredCount; })
        .def_readonly("two_way", &DpmCase::twoWay)
        .def_readonly("outer_used", &DpmCase::outerUsed)
        .def_readonly("two_way_residual", &DpmCase::twoWayResidual)
        .def_property_readonly("carrier_converged",
            [](const DpmCase& r) { return r.flow.converged; })
        .def_property_readonly("carrier_iterations",
            [](const DpmCase& r) { return r.flow.iterations; })
        .def_property_readonly("carrier_continuity_error",
            [](const DpmCase& r) { return r.flow.finalContinuityError; })
        .def_property_readonly("dt",
            [](const DpmCase& r) { return r.cfg.dt; })
        .def_property_readonly("carrier_rho",
            [](const DpmCase& r) { return r.cfg.carrierRho; })
        .def_property_readonly("drag_law",
            [](const DpmCase& r) { return r.cfg.drag; })
        .def_property_readonly("cell_drag_impulse", [](const DpmCase& r)
            {
                cf::Field<cf::Vector> imp(r.result.cellDragImpulse.begin(),
                                          r.result.cellDragImpulse.end());
                return cfbind::vectorFieldToNumpy(imp);
            })
        .def_property_readonly("source_window",
            [](const DpmCase& r) { return r.result.sourceWindow; })
        .def("n_cells", [](const DpmCase& r) { return r.mesh->nCells(); })
        .def("cell_centres", [](const DpmCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("dpm_sphere_cd", &cf::sphereDragCoefficient, py::arg("reynolds"),
          "Schiller-Naumann sphere drag coefficient Cd(Re).");
    m.def("dpm_drag_acceleration",
          [](const Vec3& uf, const Vec3& up, double rhoF, double rhoP,
             double diameter, double mu)
          {
              const cf::Vector a = cf::particleDragAcceleration(
                  cf::Vector{uf[0], uf[1], uf[2]},
                  cf::Vector{up[0], up[1], up[2]}, rhoF, rhoP, diameter, mu);
              return arr(a);
          },
          py::arg("u_fluid"), py::arg("u_particle"), py::arg("rho_f"),
          py::arg("rho_p"), py::arg("diameter"), py::arg("mu"),
          "Drag acceleration a = 3 rho_f Cd/(4 rho_p d) (U_f - U_p).");

    m.def("solve_dpm_case", &solveDpmCase, py::arg("case_dict_path"),
          "Read a DPM case.dict, solve the steady carrier and track the "
          "Lagrangian particles; returns the carrier field and trajectories.");
}
