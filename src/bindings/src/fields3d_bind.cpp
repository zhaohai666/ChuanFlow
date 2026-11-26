// 3D family: load a mesh, build cell-centred fields, run a built-in solve
// (the steady P1 radiation operator, the smallest scalar FVM solve in the
// tree) and read the field back as NumPy arrays.
//
// Conversion contract lives in field_numpy.h: float64, C-order copies,
// scalar (n,) / vector (n, 3), strict shape and dtype checks on write-back.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/fvm/Radiation.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct Mesh3d
{
    std::shared_ptr<cf::polyMesh> mesh;
};

struct ScalarField3d
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volScalarField> field;
};

struct VectorField3d
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volVectorField> field;
};

struct P1Result
{
    Mesh3d mesh;
    cf::RadiationResult result;
};

} // namespace

void registerFields3dBindings(py::module_& m)
{
    py::class_<Mesh3d>(m, "Mesh3d",
        "Polyhedral mesh loaded from a case.dict `block{}` or `mesh{}` "
        "declaration. Geometry is read-only through the binding.")
        .def(py::init([](const std::string& casePath)
             {
                 const cf::Dictionary dict = cf::Dictionary::parseFile(casePath);
                 Mesh3d out;
                 out.mesh = std::make_shared<cf::polyMesh>(
                     cf::MeshLoader::fromCase(dict, casePath));
                 return out;
             }), py::arg("case_dict_path"))
        .def("n_cells", [](const Mesh3d& self) { return self.mesh->nCells(); })
        .def("n_faces", [](const Mesh3d& self) { return self.mesh->nFaces(); })
        .def("n_patches", [](const Mesh3d& self) { return self.mesh->nPatches(); })
        .def("patch_names", [](const Mesh3d& self) { return self.mesh->patchNames(); })
        .def("cell_centres",
             [](const Mesh3d& self)
             { return cfbind::vectorFieldToNumpy(self.mesh->cellCentres()); },
             "Cell-centre coordinates, shape (n, 3), float64 C-order copy.")
        .def("cell_volumes",
             [](const Mesh3d& self)
             { return cfbind::scalarFieldToNumpy(self.mesh->cellVolumes()); },
             "Cell volumes, shape (n,), float64 C-order copy.");

    py::class_<ScalarField3d>(m, "ScalarField3d",
        "Cell-centred scalar field over a Mesh3d. `to_numpy()` returns a "
        "float64 C-order copy; `from_numpy()` writes a (n,) float64 array back "
        "into kernel storage.")
        .def(py::init([](const Mesh3d& mesh, const std::string& name, double init)
             {
                 auto f = std::make_unique<ScalarField3d>();
                 f->mesh = mesh.mesh;
                 f->field = std::make_unique<cf::volScalarField>(*mesh.mesh, name, init);
                 return f;
             }), py::arg("mesh"), py::arg("name"), py::arg("init") = 0.0)
        .def("name", [](const ScalarField3d& f) { return f.field->name(); })
        .def("n_cells", [](const ScalarField3d& f) { return f.field->internalField().size(); })
        .def("to_numpy",
             [](const ScalarField3d& f)
             { return cfbind::scalarFieldToNumpy(f.field->internalField()); },
             "Field values, shape (n,), float64 C-order copy.")
        .def("from_numpy",
             [](ScalarField3d& f, const py::array& arr)
             { cfbind::scalarFieldFromNumpy(f.field->internalField(), arr); },
             py::arg("values"),
             "Write a (n,) float64 C-order array back into the field.")
        .def("min", [](const ScalarField3d& f) { return f.field->min(); })
        .def("max", [](const ScalarField3d& f) { return f.field->max(); });

    py::class_<VectorField3d>(m, "VectorField3d",
        "Cell-centred vector field over a Mesh3d. `to_numpy()` returns a "
        "float64 C-order copy of shape (n, 3); `from_numpy()` writes a "
        "(n, 3) float64 array back.")
        .def(py::init([](const Mesh3d& mesh, const std::string& name)
             {
                 auto f = std::make_unique<VectorField3d>();
                 f->mesh = mesh.mesh;
                 f->field = std::make_unique<cf::volVectorField>(*mesh.mesh, name);
                 return f;
             }), py::arg("mesh"), py::arg("name"))
        .def("name", [](const VectorField3d& f) { return f.field->name(); })
        .def("n_cells", [](const VectorField3d& f) { return f.field->internalField().size(); })
        .def("to_numpy",
             [](const VectorField3d& f)
             { return cfbind::vectorFieldToNumpy(f.field->internalField()); },
             "Field vectors, shape (n, 3), float64 C-order copy.")
        .def("from_numpy",
             [](VectorField3d& f, const py::array& arr)
             { cfbind::vectorFieldFromNumpy(f.field->internalField(), arr); },
             py::arg("values"),
             "Write a (n, 3) float64 C-order array back into the field.");

    py::class_<P1Result>(m, "P1Result",
        "Result of one steady P1 radiation solve. `G` is the incident "
        "radiation field, `q_rad` the radiative flux vector, `net_source` "
        "the per-cell radiative heating; `balance` closes the discrete "
        "energy balance (volume source + boundary flux ~ 0).")
        .def_property_readonly("solved", [](const P1Result& r) { return r.result.solved; })
        .def_property_readonly("iterations", [](const P1Result& r) { return r.result.iterations; })
        .def_property_readonly("residual", [](const P1Result& r) { return r.result.residual; })
        .def_property_readonly("balance", [](const P1Result& r) { return r.result.balance; })
        .def_property_readonly("volume_source", [](const P1Result& r) { return r.result.volumeSource; })
        .def_property_readonly("boundary_flux", [](const P1Result& r) { return r.result.boundaryFlux; })
        .def_property_readonly("mesh", [](const P1Result& r) { return r.mesh; })
        .def("n_cells", [](const P1Result& r) { return r.result.G.size(); })
        .def("G", [](const P1Result& r)
             { return cfbind::scalarFieldToNumpy(r.result.G); },
             "Incident radiation [W/m2], shape (n,), float64 C-order copy.")
        .def("q_rad", [](const P1Result& r)
             { return cfbind::vectorFieldToNumpy(r.result.qRad); },
             "Radiative flux [W/m2], shape (n, 3), float64 C-order copy.")
        .def("net_source", [](const P1Result& r)
             { return cfbind::scalarFieldToNumpy(r.result.netSource); },
             "Per-cell net radiative heating [W/m3], shape (n,).");

    m.def("load_mesh", [](const std::string& casePath)
          {
              const cf::Dictionary dict = cf::Dictionary::parseFile(casePath);
              Mesh3d out;
              out.mesh = std::make_shared<cf::polyMesh>(
                  cf::MeshLoader::fromCase(dict, casePath));
              return out;
          },
          py::arg("case_dict_path"),
          "Load the mesh declared in a case dictionary (block{} or mesh{}).");

    m.def("solve_p1_case", [](const std::string& casePath)
          {
              const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
              const cf::RadiationConfig rad = cf::readRadiationConfig(caseDict);
              if (!rad.enabled)
                  throw std::runtime_error("solve_p1_case: needs a radiation{} declaration");

              P1Result out;
              out.mesh.mesh = std::make_shared<cf::polyMesh>(
                  cf::MeshLoader::fromCase(caseDict, casePath));

              cf::volScalarField G(*out.mesh.mesh, "G", 0.0);
              if (caseDict.has("patches"))
              {
                  const cf::Dictionary& ps = caseDict.getDict("patches");
                  for (const auto& name : ps.keys())
                  {
                      const cf::Dictionary& pd = ps.getDict(name);
                      const cf::label patch = out.mesh.mesh->patchIdByName(name);
                      if (pd.has("wallTemperature"))
                      {
                          const double tw = pd.getScalar("wallTemperature");
                          G.setPatchValue(patch, "fixedValue",
                                          cf::blackbodyIncidentRadiation(tw, rad.sbConstant));
                      }
                      else
                      {
                          G.setPatchValue(patch, "zeroGradient", 0.0);
                      }
                  }
              }

              const std::size_t n = static_cast<std::size_t>(out.mesh.mesh->nCells());
              cf::Field<cf::scalar> mediumT(n, rad.mediumTemperature);
              const cf::Schemes schemes = cf::loadSchemes(caseDict);

              out.result = cf::solveP1(*out.mesh.mesh, G, rad, schemes, mediumT);
              if (!out.result.solved)
                  throw std::runtime_error("solve_p1_case: P1 solve did not converge");
              return out;
          },
          py::arg("case_dict_path"),
          "Run the steady P1 radiation solve declared in a case dictionary and "
          "return the converged fields as NumPy arrays.");
}
