// ALE moving-mesh family binding: rigid translation / uniform expansion vertex
// motion and the Geometric Conservation Law residual per cell.
//
// Replicates src/3d/apps/aleApp.cpp's driver: readAleConfig, mesh load, then
// per step vertex velocity at the current positions, advanceVertices,
// movePoints, gclResidual on the new geometry. The motion and the GCL identity
// live in the kernel Ale module; this file drives it and copies the report rows
// and final geometry out.
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Ale.h"
#include "cf/fvm/MeshLoader.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct AleRow
{
    int step = 0;
    double t = 0;
    double totalVolume = 0;
    double maxGcl = 0;
    double sumGcl = 0;
    double xMin = 0;
    double xMax = 0;
};

struct AleCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    cf::AleConfig cfg;
    std::vector<AleRow> history;
    double initialVolume = 0;
};

std::unique_ptr<AleCase> solveAleCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<AleCase>();
    out->cfg = cf::readAleConfig(caseDict);
    if (!out->cfg.enabled)
        throw std::runtime_error("ale: needs an ale{} declaration");
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    cf::Field<cf::scalar> oldVol = out->mesh->cellVolumes();
    cf::scalar totalV0 = 0;
    for (cf::scalar v : oldVol) totalV0 += v;
    out->initialVolume = totalV0;

    for (int step = 1; step <= out->cfg.steps; ++step)
    {
        const cf::Field<cf::Vector> vVert = cf::meshVertexVelocity(*out->mesh, out->cfg);
        const cf::Field<cf::Vector> newPoints = cf::advanceVertices(*out->mesh, out->cfg);
        out->mesh->movePoints(newPoints);
        const cf::Field<cf::scalar> res =
            cf::gclResidual(*out->mesh, oldVol, vVert, out->cfg.dt);

        cf::scalar maxRes = 0, sumRes = 0, totalV = 0;
        for (cf::scalar r : res)
        {
            maxRes = std::max(maxRes, std::abs(r));
            sumRes += r;
        }
        for (cf::scalar v : out->mesh->cellVolumes()) totalV += v;

        if (step % out->cfg.reportEvery == 0 || step == out->cfg.steps)
        {
            AleRow row;
            row.step = step;
            row.t = step * out->cfg.dt;
            row.totalVolume = totalV;
            row.maxGcl = maxRes;
            row.sumGcl = sumRes;
            row.xMin = out->mesh->minX().x;
            row.xMax = out->mesh->maxX().x;
            out->history.push_back(row);
        }
        oldVol = out->mesh->cellVolumes();
    }
    return out;
}

} // namespace

void registerAleBindings(py::module_& m)
{
    py::class_<AleRow>(m, "AleRow",
        "One reported ALE step: time, total mesh volume, per-cell GCL residual "
        "max/sum and the current x bounding box.")
        .def_readonly("step", &AleRow::step)
        .def_readonly("t", &AleRow::t)
        .def_readonly("total_volume", &AleRow::totalVolume)
        .def_readonly("max_gcl", &AleRow::maxGcl)
        .def_readonly("sum_gcl", &AleRow::sumGcl)
        .def_readonly("x_min", &AleRow::xMin)
        .def_readonly("x_max", &AleRow::xMax);

    py::class_<AleCase, std::unique_ptr<AleCase>>(m, "AleCase",
        "Result of one ALE moving-mesh run: motion history, final cell volumes "
        "and vertex/centre coordinates, and the motion configuration.")
        .def_readonly("history", &AleCase::history)
        .def_property_readonly("initial_volume",
            [](const AleCase& r) { return r.initialVolume; })
        .def_property_readonly("volumes", [](const AleCase& r)
            { return cfbind::scalarFieldToNumpy(r.mesh->cellVolumes()); })
        .def_property_readonly("points", [](const AleCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->points()); })
        .def_property_readonly("cell_centres", [](const AleCase& r)
            { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); })
        .def_property_readonly("motion", [](const AleCase& r) { return r.cfg.motion; })
        .def_property_readonly("translation", [](const AleCase& r)
            { return std::array<double, 3>{r.cfg.translation.x, r.cfg.translation.y,
                                          r.cfg.translation.z}; })
        .def_property_readonly("expansion_rate",
            [](const AleCase& r) { return r.cfg.expansionRate; })
        .def_property_readonly("dt", [](const AleCase& r) { return r.cfg.dt; })
        .def_property_readonly("steps", [](const AleCase& r) { return r.cfg.steps; })
        .def("n_cells", [](const AleCase& r) { return r.mesh->nCells(); })
        .def("n_points", [](const AleCase& r) { return r.mesh->nPoints(); });

    m.def("solve_ale_case", &solveAleCase, py::arg("case_dict_path"),
          "Read an ALE moving-mesh case.dict (translate / expand), advance the "
          "vertices and report the per-step GCL residual and mesh geometry.");
}
