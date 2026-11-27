// Dynamic-mesh 6DoF family binding: GCL-corrected relative-flux scalar transport
// on a rigidly translating mesh, integrated alongside a rigid body with one
// translational and one rotational DOF.
//
// Replicates src/3d/apps/dynamicMesh6DofApp.cpp's driver: readDynamicMeshConfig,
// mesh load, uniform scalar seed, prescribed rigid translation, relative face
// flux, GCL sweep/residual, transportStep, and the RK4 body update. The
// transport correction, body ODE integration and GCL identity live in the kernel
// DynamicMesh6Dof module; this file drives it and copies state out.
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/DynamicMesh6Dof.h"
#include "cf/fvm/MeshLoader.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct DynRow
{
    double t = 0;
    double x = 0;
    double v = 0;
    double theta = 0;
    double omega = 0;
    double scalarMean = 0;
    double massBalance = 0;
};

struct DynCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    cf::DynamicMeshConfig cfg;
    cf::BodyState body;
    std::vector<DynRow> history;
    cf::Field<cf::scalar> phi;
};

std::unique_ptr<DynCase> solveDynMeshCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<DynCase>();
    out->cfg = cf::readDynamicMeshConfig(caseDict);
    if (!out->cfg.enabled)
        throw std::runtime_error("dynamicMesh: needs a dynamicMesh{} declaration");
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    out->phi = cf::Field<cf::scalar>(static_cast<std::size_t>(out->mesh->nCells()),
                                     out->cfg.initialValue);
    cf::Field<cf::scalar> oldVol = out->mesh->cellVolumes();

    out->body.x = 0;
    out->body.v = 0;
    out->body.theta = out->cfg.theta0;
    out->body.omega = out->cfg.omega0;

    // Initial state (t = 0): no motion has happened, so the GCL residual is 0.
    {
        DynRow row;
        row.t = 0.0;
        row.x = out->body.x;
        row.v = out->body.v;
        row.theta = out->body.theta;
        row.omega = out->body.omega;
        row.scalarMean = out->cfg.initialValue;
        row.massBalance = 0.0;
        out->history.push_back(row);
    }

    for (int step = 1; step <= out->cfg.steps; ++step)
    {
        out->body = cf::advanceBody(out->body, out->cfg, out->cfg.dt);

        const cf::Field<cf::Vector>& pts = out->mesh->points();
        cf::Field<cf::Vector> newPoints(static_cast<std::size_t>(pts.size()));
        for (std::size_t i = 0; i < pts.size(); ++i)
            newPoints[i] = pts[i] + out->cfg.meshVelocity * out->cfg.dt;
        out->mesh->movePoints(newPoints);

        const cf::Field<cf::scalar> phiRel =
            cf::relativeFaceFlux(*out->mesh, out->cfg.fluidVelocity, out->cfg.meshVelocity);
        const cf::Field<cf::scalar> sweep =
            cf::gclGeometricSweep(*out->mesh, out->cfg.meshVelocity);
        const cf::Field<cf::scalar> gclRes =
            cf::gclCellResidual(*out->mesh, out->cfg.meshVelocity, oldVol, out->cfg.dt);
        cf::scalar massBalance = 0;
        for (cf::scalar r : gclRes) massBalance = std::max(massBalance, std::abs(r));

        out->phi = cf::transportStep(*out->mesh, out->phi, phiRel, sweep, oldVol,
                                    out->cfg.dt);
        oldVol = out->mesh->cellVolumes();

        if (step % out->cfg.reportEvery == 0 || step == out->cfg.steps)
        {
            // Recompute the mean on the current geometry after the transport.
            cf::scalar totalMass = 0, totalVol = 0;
            const auto& vol = out->mesh->cellVolumes();
            for (cf::label c = 0; c < out->mesh->nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                totalMass += out->phi[u] * vol[u];
                totalVol += vol[u];
            }
            DynRow row;
            row.t = step * out->cfg.dt;
            row.x = out->body.x;
            row.v = out->body.v;
            row.theta = out->body.theta;
            row.omega = out->body.omega;
            row.scalarMean = totalMass / totalVol;
            row.massBalance = massBalance;
            out->history.push_back(row);
        }
    }
    return out;
}

} // namespace

void registerDynMeshBindings(py::module_& m)
{
    py::class_<DynRow>(m, "DynMeshRow",
        "One reported 6DoF state: body translation (x, v), rotation "
        "(theta, omega), uniform scalar mean and the max per-cell GCL residual.")
        .def_readonly("t", &DynRow::t)
        .def_readonly("x", &DynRow::x)
        .def_readonly("v", &DynRow::v)
        .def_readonly("theta", &DynRow::theta)
        .def_readonly("omega", &DynRow::omega)
        .def_readonly("scalar_mean", &DynRow::scalarMean)
        .def_readonly("mass_balance", &DynRow::massBalance);

    py::class_<DynCase, std::unique_ptr<DynCase>>(m, "DynMeshCase",
        "Result of one moving-mesh 6DoF run: final transported scalar, body "
        "history and the model configuration.")
        .def_property_readonly("phi", [](const DynCase& r)
            { return cfbind::scalarFieldToNumpy(r.phi); })
        .def_readonly("history", &DynCase::history)
        .def_property_readonly("mesh_velocity", [](const DynCase& r)
            { return std::array<double, 3>{r.cfg.meshVelocity.x, r.cfg.meshVelocity.y,
                                          r.cfg.meshVelocity.z}; })
        .def_property_readonly("body_force", [](const DynCase& r)
            { return std::array<double, 3>{r.cfg.bodyForce.x, r.cfg.bodyForce.y,
                                          r.cfg.bodyForce.z}; })
        .def_property_readonly("mass", [](const DynCase& r) { return r.cfg.mass; })
        .def_property_readonly("drag", [](const DynCase& r) { return r.cfg.drag; })
        .def_property_readonly("inertia", [](const DynCase& r) { return r.cfg.inertia; })
        .def_property_readonly("rot_drag", [](const DynCase& r) { return r.cfg.rotDrag; })
        .def_property_readonly("rot_stiffness",
            [](const DynCase& r) { return r.cfg.rotStiffness; })
        .def_property_readonly("theta0", [](const DynCase& r) { return r.cfg.theta0; })
        .def_property_readonly("omega0", [](const DynCase& r) { return r.cfg.omega0; })
        .def_property_readonly("initial_value",
            [](const DynCase& r) { return r.cfg.initialValue; })
        .def_property_readonly("dt", [](const DynCase& r) { return r.cfg.dt; })
        .def_property_readonly("steps", [](const DynCase& r) { return r.cfg.steps; })
        .def("n_cells", [](const DynCase& r) { return r.mesh->nCells(); });

    m.def("solve_dynmesh_case", &solveDynMeshCase, py::arg("case_dict_path"),
          "Read a dynamic-mesh 6DoF case.dict, transport the uniform scalar on the "
          "rigidly moving mesh and advance the translational/rotational body ODE.");
}
