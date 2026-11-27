// Scalar transport family: read a steady (or Euler-transient) advection-diffusion
// case.dict, assemble the scalar equation exactly like scalarTransportApp, and
// hand the solved cell-centred scalar field back to Python as NumPy arrays.
//
// Conversion contract lives in field_numpy.h: float64, C-order copies, scalar
// (n,) arrays; the mesh points / cell centres are read-only copies as well.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/common/Matrix.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct ScalarTransportResult
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volScalarField> T;
    cf::SolverReport report;
};

std::unique_ptr<ScalarTransportResult> solveScalarTransport(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);

    auto out = std::make_unique<ScalarTransportResult>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    // ---- velocity-driven face flux phi = U0 & Sf ----
    const std::vector<double> vlist = caseDict.getScalarList("velocity");
    if (vlist.size() != 3)
        throw std::runtime_error("scalar transport: velocity must be a 3-vector");
    const cf::Vector U0{vlist[0], vlist[1], vlist[2]};

    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);
    for (cf::label f = 0; f < out->mesh->nInternalFaces(); ++f)
    {
        phi.values()[static_cast<std::size_t>(f)] =
            cf::dot(U0, out->mesh->faceSf()[static_cast<std::size_t>(f)]);
    }
    phi.fixBoundaryVelocity(U0);

    // ---- scalar field and its boundary conditions ----
    out->T = std::make_unique<cf::volScalarField>(
        *out->mesh, "T", caseDict.getScalarOr("initialValue", 0.0));
    if (caseDict.has("patches"))
    {
        const cf::Dictionary& ps = caseDict.getDict("patches");
        for (const auto& name : ps.keys())
        {
            const cf::Dictionary& pd = ps.getDict(name);
            const cf::label patch = out->mesh->patchIdByName(name);
            out->T->setPatchValue(patch, pd.getWord("type"),
                                   pd.getScalarOr("value", 0.0));
        }
    }

    // ---- equation ----
    const cf::Schemes schemes = cf::loadSchemes(caseDict);
    const double D = caseDict.getScalar("diffusivity");

    cf::fvMatrix eqn(*out->mesh, cf::viewOf(*out->T), "TEqn");
    eqn.schemes() = schemes;
    if (schemes.ddt == "Euler")
    {
        const double dt = caseDict.getScalar("deltaT");
        cf::fvm::ddt(eqn, out->T->internalField(), dt);
    }
    cf::fvm::div(eqn, phi.values(), schemes.div);
    cf::fvm::laplacian(eqn, D, schemes.laplacian, -1);

    cf::fvMatrix::Controls ctl;
    ctl.tol = caseDict.getScalarOr("solverTol", 1e-8);
    ctl.maxIter = static_cast<int>(caseDict.getScalarOr("solverMaxIter", 5000));
    ctl.nonOrth = 2;
    out->report = eqn.solve(ctl);

    out->T->internalField() = eqn.solution()[0];
    if (!out->report.converged)
        throw std::runtime_error("solve_scalar_transport_case: scalar solve did "
                                 "not converge");
    return out;
}

} // namespace

void registerScalarTransportBindings(py::module_& m)
{
    py::class_<ScalarTransportResult, std::unique_ptr<ScalarTransportResult>>(
        m, "ScalarTransportResult",
        "Result of one scalar advection-diffusion solve. `T` is the cell-centred"
        " solved scalar; `cell_centres`/`points` give the geometry the analytic"
        " profile is compared against.")
        .def_property_readonly("converged",
            [](const ScalarTransportResult& r) { return r.report.converged; })
        .def_property_readonly("iterations",
            [](const ScalarTransportResult& r) { return r.report.iterations; })
        .def_property_readonly("residual",
            [](const ScalarTransportResult& r) { return r.report.finalResidual; })
        .def("n_cells", [](const ScalarTransportResult& r)
             { return r.mesh->nCells(); })
        .def("T", [](const ScalarTransportResult& r)
             { return cfbind::scalarFieldToNumpy(r.T->internalField()); },
             "Solved cell-centred scalar, shape (n,), float64 C-order copy.")
        .def("cell_centres",
            [](const ScalarTransportResult& r)
             { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); },
             "Cell-centre coordinates, shape (n, 3), float64 C-order copy.")
        .def("points",
            [](const ScalarTransportResult& r)
             { return cfbind::vectorFieldToNumpy(r.mesh->points()); },
             "Mesh vertex coordinates, shape (nP, 3), float64 C-order copy.")
        .def("min", [](const ScalarTransportResult& r) { return r.T->min(); })
        .def("max", [](const ScalarTransportResult& r) { return r.T->max(); });

    m.def("solve_scalar_transport_case", &solveScalarTransport,
          py::arg("case_dict_path"),
          "Read a scalar transport case.dict and solve the advection-diffusion"
          " equation; returns the converged cell-centred field.");
}
