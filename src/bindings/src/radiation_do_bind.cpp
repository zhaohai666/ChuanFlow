// Grey radiation family binding beyond P1: discrete ordinates (S_N), the
// surface-to-surface (S2S) enclosure limit, and the delegated P1 trend branch.
//
// Replicates src/3d/apps/radiationDoS2sApp.cpp's driver: readDoS2sConfig,
// slab geometry/wall temperatures, then solveDoSlab / solveS2sSlab / the P1
// delegate. The RTE sweep, view-factor algebra and P1 diffusion live in the
// kernel; this file drives it and copies the fields and balance out to NumPy.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/Radiation.h"
#include "cf/fvm/RadiationDoS2s.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct RadCase
{
    std::shared_ptr<cf::polyMesh> mesh;
    cf::DoS2sConfig cfg;
    cf::RadSlabResult res;
    std::string model;
    int nCells = 0;
    double length = 0;
    double TL = 0;
    double TR = 0;
};

struct Slab
{
    int nCells = 0;
    cf::scalar length = 0;
    cf::scalar TL = 0;
    cf::scalar TR = 0;
};

Slab readSlab(const cf::Dictionary& caseDict, const cf::polyMesh& mesh)
{
    Slab s;
    const cf::Dictionary& blk = caseDict.getDict("block");
    const auto mn = blk.getScalarList("min");
    const auto mx = blk.getScalarList("max");
    const auto nc = blk.getScalarList("cells");
    if (mn.size() < 3 || mx.size() < 3 || nc.size() < 3)
        throw std::runtime_error("radiation DO: block needs min/max/cells vectors");
    s.length = mx[0] - mn[0];
    s.nCells = static_cast<int>(nc[0]);
    if (!(s.length > 0) || !(s.nCells > 0))
        throw std::runtime_error("radiation DO: slab needs positive length and cell count");
    if (static_cast<int>(nc[1]) != 1 || static_cast<int>(nc[2]) != 1 ||
        mesh.nCells() != s.nCells)
        throw std::runtime_error("radiation DO: only a 1D slab (cells y=z=1) is supported");

    if (!caseDict.has("patches"))
        throw std::runtime_error("radiation DO: needs patches{} with wallTemperature");
    const cf::Dictionary& ps = caseDict.getDict("patches");
    if (!ps.has("left") || !ps.has("right"))
        throw std::runtime_error("radiation DO: slab needs patches 'left' and 'right'");
    s.TL = ps.getDict("left").getScalar("wallTemperature");
    s.TR = ps.getDict("right").getScalar("wallTemperature");
    if (!(s.TL > 0) || !(s.TR > 0))
        throw std::runtime_error("radiation DO: wallTemperature must be finite and positive");
    return s;
}

const char* modelName(cf::RadiationModel m)
{
    switch (m)
    {
    case cf::RadiationModel::DO: return "do";
    case cf::RadiationModel::S2S: return "s2s";
    case cf::RadiationModel::P1: return "p1";
    }
    return "unknown";
}

cf::RadSlabResult solveP1Branch(const cf::Dictionary& caseDict,
                                const cf::polyMesh& mesh,
                                const cf::DoS2sConfig& cfg, const Slab& s)
{
    cf::RadiationConfig rp1 = cf::readRadiationConfig(caseDict);
    if (!rp1.enabled)
        throw std::runtime_error("radiation DO: p1 needs a radiation{} declaration");
    cf::volScalarField G(mesh, "G", 0.0);
    if (caseDict.has("patches"))
    {
        const cf::Dictionary& ps = caseDict.getDict("patches");
        for (const auto& name : ps.keys())
        {
            const cf::Dictionary& pd = ps.getDict(name);
            const cf::label patch = mesh.patchIdByName(name);
            if (pd.has("wallTemperature"))
            {
                const cf::scalar tw = pd.getScalar("wallTemperature");
                if (!(tw > 0))
                    throw std::runtime_error("radiation DO: wallTemperature must be positive");
                G.setPatchValue(patch, "fixedValue",
                                cf::blackbodyIncidentRadiation(tw, cfg.sbConstant));
            }
            else
            {
                G.setPatchValue(patch, "zeroGradient", 0.0);
            }
        }
    }
    const cf::Field<cf::scalar> mediumT(
        static_cast<std::size_t>(mesh.nCells()), rp1.mediumTemperature);
    const cf::Schemes schemes = cf::loadSchemes(caseDict);
    const cf::RadiationResult rr = cf::solveP1(mesh, G, rp1, schemes, mediumT);
    if (!rr.solved)
        throw std::runtime_error("radiation DO: delegated P1 solve did not converge");

    cf::RadSlabResult out;
    const int N = s.nCells;
    out.G.assign(N, 0.0);
    out.q.assign(N, cf::Vector::zero());
    out.S.assign(N, 0.0);
    for (int i = 0; i < N; ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        out.G[u] = rr.G[u];
        out.q[u] = rr.qRad[u];
        out.S[u] = rr.netSource[u];
    }
    out.qLeft = rr.qRad.empty() ? 0 : rr.qRad[0].x;
    out.qRight = rr.qRad.empty() ? 0 : -rr.qRad[static_cast<std::size_t>(N - 1)].x;
    out.volSrc = rr.volumeSource;
    out.balance = rr.balance;
    out.iterations = rr.iterations;
    out.residual = rr.residual;
    out.converged = rr.solved;
    return out;
}

std::unique_ptr<RadCase> solveRadCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<RadCase>();
    out->cfg = cf::readDoS2sConfig(caseDict);
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));
    const Slab slab = readSlab(caseDict, *out->mesh);
    out->nCells = slab.nCells;
    out->length = slab.length;
    out->TL = slab.TL;
    out->TR = slab.TR;
    out->model = modelName(out->cfg.model);

    switch (out->cfg.model)
    {
    case cf::RadiationModel::DO:
        out->res = cf::solveDoSlab(slab.nCells, slab.length, slab.TL, slab.TR, out->cfg);
        break;
    case cf::RadiationModel::S2S:
        out->res = cf::solveS2sSlab(slab.nCells, slab.length, slab.TL, slab.TR, out->cfg);
        break;
    case cf::RadiationModel::P1:
        out->res = solveP1Branch(caseDict, *out->mesh, out->cfg, slab);
        break;
    }
    return out;
}

} // namespace

void registerRadiationDoBindings(py::module_& m)
{
    py::class_<RadCase, std::unique_ptr<RadCase>>(m, "RadCase",
        "Result of one grey-radiation slab solve: incident radiation G, net flux "
        "q and volumetric source S per cell, wall fluxes, discrete balance and "
        "the model/quadrature constants.")
        .def_property_readonly("G", [](const RadCase& r)
            { return cfbind::scalarFieldToNumpy(r.res.G); })
        .def_property_readonly("q", [](const RadCase& r)
            { return cfbind::vectorFieldToNumpy(r.res.q); })
        .def_property_readonly("S", [](const RadCase& r)
            { return cfbind::scalarFieldToNumpy(r.res.S); })
        .def_property_readonly("q_left", [](const RadCase& r) { return r.res.qLeft; })
        .def_property_readonly("q_right", [](const RadCase& r) { return r.res.qRight; })
        .def_property_readonly("vol_source", [](const RadCase& r) { return r.res.volSrc; })
        .def_property_readonly("balance", [](const RadCase& r) { return r.res.balance; })
        .def_property_readonly("iterations", [](const RadCase& r) { return r.res.iterations; })
        .def_property_readonly("residual", [](const RadCase& r) { return r.res.residual; })
        .def_property_readonly("converged", [](const RadCase& r) { return r.res.converged; })
        .def_property_readonly("view_factor", [](const RadCase& r) { return r.res.viewFactor; })
        .def_readonly("model", &RadCase::model)
        .def_property_readonly("absorption", [](const RadCase& r) { return r.cfg.absorption; })
        .def_property_readonly("scattering", [](const RadCase& r) { return r.cfg.scattering; })
        .def_property_readonly("medium_temperature",
            [](const RadCase& r) { return r.cfg.mediumTemperature; })
        .def_property_readonly("quadrature_order",
            [](const RadCase& r) { return r.cfg.quadratureOrder; })
        .def_property_readonly("n_cells", [](const RadCase& r) { return r.nCells; })
        .def_property_readonly("length", [](const RadCase& r) { return r.length; })
        .def_property_readonly("T_left", [](const RadCase& r) { return r.TL; })
        .def_property_readonly("T_right", [](const RadCase& r) { return r.TR; });

    m.def("dos2s_blackbody_g", &cf::dos2sBlackbodyG,
          py::arg("temperature"), py::arg("sigma_sb"),
          "Blackbody incident radiation G_w = 4 sigma T^4.");

    m.def("solve_radiation_case", &solveRadCase, py::arg("case_dict_path"),
          "Read a grey-radiation slab case.dict (do / s2s / p1) and solve the "
          "radiative transfer on the 1D black-wall slab.");
}
