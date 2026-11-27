// Steady one-step combustion family (F + s O -> products, Arrhenius rate).
//
// Replicates src/3d/apps/combustionApp.cpp's driver: constant-velocity face
// flux, species patch BCs, solveCombustion (linear fuel equation when the
// oxidiser is pinned, coupled Picard otherwise), and the boundary-flux /
// sink diagnostics the app reports. The chemistry and the rate law live in the
// kernel CombustionModel; this file drives it and copies the state out.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Combustion.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct CombustionCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    cf::CombustionConfig cfg;
    cf::CombustionResult res;
    double fuelIn = 0, fuelOut = 0, oxIn = 0, oxOut = 0;
    bool fixedOx = false;
};

void applyFluxBC(cf::surfaceScalarField& phi, const cf::polyMesh& mesh,
                 const cf::Vector& U0)
{
    for (cf::label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        phi.values()[static_cast<std::size_t>(f)] =
            cf::dot(U0, mesh.faceSf()[static_cast<std::size_t>(f)]);
    }
    for (cf::label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        std::string type;
        for (const auto& pa : mesh.patches())
            if (f >= pa.start && f < pa.start + pa.size) { type = pa.type; break; }
        phi.values()[u] = (type == "wall" || type == "symmetry" || type == "empty")
                              ? 0.0 : cf::dot(U0, mesh.faceSf()[u]);
    }
}

std::unique_ptr<CombustionCase> solveCombustionCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    const cf::CombustionConfig cfg = cf::readCombustionConfig(caseDict);
    if (!cfg.enabled)
        throw std::runtime_error("combustion: needs a combustion{} declaration");

    auto out = std::make_unique<CombustionCase>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));
    out->cfg = cfg;

    const auto vl = caseDict.getScalarList("velocity");
    const cf::Vector U0{vl[0], vl[1], vl[2]};
    cf::surfaceScalarField phi(*out->mesh, "phi", 0.0);
    applyFluxBC(phi, *out->mesh, U0);

    cf::volScalarField Yf(*out->mesh, "Yf", 0.0);
    cf::volScalarField Yo(*out->mesh, "Yo", 0.0);
    if (caseDict.has("patches"))
    {
        const cf::Dictionary& ps = caseDict.getDict("patches");
        for (const std::string& name : ps.keys())
        {
            const cf::Dictionary& pd = ps.getDict(name);
            const cf::label patch = out->mesh->patchIdByName(name);
            Yf.setPatchValue(patch, pd.has("fuelValue") ? "fixedValue" : "zeroGradient",
                             pd.getScalarOr("fuelValue", 0.0));
            Yo.setPatchValue(patch, pd.has("oxValue") ? "fixedValue" : "zeroGradient",
                             pd.getScalarOr("oxValue", 0.0));
        }
    }

    out->fixedOx = cfg.fixedOxidiser >= 0;
    out->res = cf::solveCombustion(*out->mesh, phi, Yf, Yo, cfg,
                                   cf::loadSchemes(caseDict), out->fixedOx);

    // Boundary fluxes from the face BC linearisation, mirroring the app.
    const auto& sf = out->mesh->faceSf();
    const auto& dc = out->mesh->deltaCoeffVec();
    for (cf::label f = out->mesh->nInternalFaces(); f < out->mesh->nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const cf::label owner = out->mesh->owner()[u];
        const cf::Vector& d = dc[u];
        const double d2 = cf::dot(d, d);
        if (!(d2 > 0)) continue;
        const double oc = cf::dot(sf[u], d) / d2;
        double cf, kf, co, ko;
        Yf.linearise(f, cf, kf);
        Yo.linearise(f, co, ko);
        const double yfOwner = out->res.Yf[static_cast<std::size_t>(owner)];
        const double yoOwner = out->res.Yo[static_cast<std::size_t>(owner)];
        const double phiF = phi.values()[u];
        const double fFace = cf * yfOwner + kf;
        const double oFace = co * yoOwner + ko;
        const double fuelFlux = phiF * fFace - cfg.diffusivity * oc * (fFace - yfOwner);
        const double oxFlux = phiF * oFace - cfg.diffusivity * oc * (oFace - yoOwner);
        if (fuelFlux < 0) out->fuelIn += -fuelFlux; else out->fuelOut += fuelFlux;
        if (oxFlux < 0) out->oxIn += -oxFlux; else out->oxOut += oxFlux;
    }
    return out;
}

} // namespace

void registerCombustionBindings(py::module_& m)
{
    py::class_<CombustionCase, std::unique_ptr<CombustionCase>>(
        m, "CombustionCase",
        "Result of one steady one-step combustion solve. `Yf`/`Yo` are the"
        " solved species mass fractions; `reaction_rate` is k Yf Yo per cell;"
        " the flux/sink diagnostics close the discrete species balance.")
        .def_property_readonly("Yf", [](const CombustionCase& r)
            { return cfbind::scalarFieldToNumpy(r.res.Yf); })
        .def_property_readonly("Yo", [](const CombustionCase& r)
            { return cfbind::scalarFieldToNumpy(r.res.Yo); })
        .def_property_readonly("reaction_rate", [](const CombustionCase& r)
            { return cfbind::scalarFieldToNumpy(r.res.reactionRate); })
        .def_property_readonly("fuel_in", [](const CombustionCase& r) { return r.fuelIn; })
        .def_property_readonly("fuel_out", [](const CombustionCase& r) { return r.fuelOut; })
        .def_property_readonly("ox_in", [](const CombustionCase& r) { return r.oxIn; })
        .def_property_readonly("ox_out", [](const CombustionCase& r) { return r.oxOut; })
        .def_property_readonly("fuel_consumed",
            [](const CombustionCase& r) { return r.res.fuelConsumed; })
        .def_property_readonly("ox_consumed",
            [](const CombustionCase& r) { return r.res.oxConsumed; })
        .def_property_readonly("fuel_iters",
            [](const CombustionCase& r) { return r.res.fuelIters; })
        .def_property_readonly("ox_iters",
            [](const CombustionCase& r) { return r.res.oxIters; })
        .def_property_readonly("solved", [](const CombustionCase& r) { return r.res.solved; })
        .def_property_readonly("fixed_oxidiser",
            [](const CombustionCase& r) { return r.fixedOx; })
        .def_property_readonly("diffusivity",
            [](const CombustionCase& r) { return r.cfg.diffusivity; })
        .def_property_readonly("stoich",
            [](const CombustionCase& r) { return r.cfg.stoich; })
        .def_property_readonly("temperature",
            [](const CombustionCase& r) { return r.cfg.temperature; })
        .def_property_readonly("fixed_oxidiser_value",
            [](const CombustionCase& r) { return r.cfg.fixedOxidiser; })
        .def_property_readonly("fuel_inlet",
            [](const CombustionCase& r) { return r.cfg.fuelInlet; })
        .def_property_readonly("rate_constant",
            [](const CombustionCase& r) { return r.cfg.rateConstant(); })
        .def("n_cells", [](const CombustionCase& r) { return r.mesh->nCells(); })
        .def("cell_centres",
            [](const CombustionCase& r)
             { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); })
        .def("arrhenius_rate",
            [](const CombustionCase& r, double T)
            { return cf::arrheniusRate(r.cfg.preExp, r.cfg.activationE,
                                        r.cfg.gasConstant, T); },
            py::arg("temperature"),
            "k(T) = A exp(-Ea/(R T)) for an arbitrary temperature.");

    m.def("solve_combustion_case", &solveCombustionCase,
          py::arg("case_dict_path"),
          "Read a steady one-step combustion case.dict and solve the species"
          " transport with the Arrhenius reaction rate.");
}
