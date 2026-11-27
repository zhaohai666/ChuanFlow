// Coupling family: read a 1D<->3D coupled case.dict, run the partitioned
// steady interface solve, and hand the structured result back to Python.
//
// The construction path mirrors src/coupling/apps/coupledApp.cpp and
// multiCoupledApp.cpp: one DetailBlock per coupled edge, one fresh network
// per run. Both case layouts are accepted:
//   * single edge:  couple { edge X; seeds ( ... ); }  detail { ... }
//   * multi edge:   couple { edges ( X1 X2 ); starts { ... }; }  details { ... }
//
// The exported per-edge shape follows the kernel `EdgeExchange` / `BlockMeasure`
// one-for-one (snake_case spelling), so multi-edge cases index their edges the
// same way the GUI side does.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/coupling/CoupledNetwork.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/SimpleSolver.h"
#include "cf/turbulence/TurbulenceModel.h"

#include "coupling_views.h"
#include "field_numpy.h"

namespace py = pybind11;

namespace
{

using cfbind::BlockMeasureView;
using cfbind::EdgeExchangeView;
using cfbind::toView;

cf::SimpleControls loadSimple(const cf::Dictionary& d)
{
    cf::SimpleControls ctl;
    if (!d.has("simple")) return ctl;
    const cf::Dictionary& sd = d.getDict("simple");
    ctl.maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", ctl.maxIterations));
    ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
    ctl.pcorrTol = sd.getScalarOr("pcorrTol", ctl.pcorrTol);
    ctl.momentumTol = sd.getScalarOr("momentumTol", ctl.momentumTol);
    ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
    ctl.linearMaxIter = static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
    ctl.nonOrthCorrectors = static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
    ctl.relaxU = sd.getScalarOr("relaxU", ctl.relaxU);
    ctl.relaxP = sd.getScalarOr("relaxP", ctl.relaxP);
    ctl.pRelaxStart = sd.getScalarOr("pRelaxStart", ctl.pRelaxStart);
    return ctl;
}


/// Owns the network, the partner blocks and the result of one coupled run.
struct CoupledCase
{
    cf::OneDNetwork net;
    std::vector<cf::label> edges;
    std::vector<std::string> edgeNames;
    std::vector<std::unique_ptr<cf::coupling::DetailBlock>> blocks;
    cf::coupling::CoupledResult result;
};

std::unique_ptr<CoupledCase> runCoupled(const std::string& casePath,
                                        std::size_t seedIndex)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);

    const cf::Dictionary& fd = caseDict.getDict("fluid");
    const double rho = fd.getScalar("rho");
    const double mu = fd.getScalar("mu");
    if (rho <= 0 || mu <= 0)
        throw std::runtime_error("coupled: fluid rho and mu must both be positive");
    const double nu = mu / rho;

    const cf::SteadyControls netCtl = cf::loadSteadyControls(caseDict);
    const cf::Dictionary& cd = caseDict.getDict("couple");

    // Edge list: multi layout `edges ( ... )`, single layout `edge X`.
    std::vector<std::string> edgeNames;
    if (cd.has("edges"))
        edgeNames = cd.getWordList("edges");
    else
        edgeNames.push_back(cd.getWord("edge"));
    if (edgeNames.empty())
        throw std::runtime_error("couple{} names no coupled edge");

    // Detail dict per edge: top-level `detail{}` for the single layout,
    // `details{ <edge> {} }` for the multi one.
    const bool multiDetails = caseDict.has("details");

    // Starting coefficients: single `seeds ( R ... )`, multi
    // `starts { <set> { <edge> R; ... } }`. One run per seed/start set; the
    // binding takes the requested one (default 0).
    std::vector<double> seedRs(edgeNames.size(), 0.0);
    if (cd.has("starts"))
    {
        const cf::Dictionary& starts = cd.getDict("starts");
        const std::vector<std::string> names = starts.keys();
        if (names.empty())
            throw std::runtime_error("couple.starts{} declares no starting set");
        if (seedIndex >= names.size())
            throw std::out_of_range("couple: seed index " + std::to_string(seedIndex)
                                    + " out of range (" + std::to_string(names.size())
                                    + " starting sets)");
        const cf::Dictionary& sd = starts.getDict(names[seedIndex]);
        for (std::size_t i = 0; i < edgeNames.size(); ++i)
            seedRs[i] = sd.getScalar(edgeNames[i]);
    }
    else
    {
        const std::vector<double> seeds = cd.getScalarList("seeds");
        if (seeds.empty())
            throw std::runtime_error("couple{} needs seeds ( ) or starts{}");
        if (seedIndex >= seeds.size())
            throw std::out_of_range("couple: seed index " + std::to_string(seedIndex)
                                    + " out of range (" + std::to_string(seeds.size())
                                    + " seeds)");
        for (std::size_t i = 0; i < edgeNames.size(); ++i)
            seedRs[i] = seeds[seedIndex];
    }

    cf::coupling::CoupledControls ctl;
    ctl.maxOuterIterations = static_cast<cf::label>(
        cd.getScalarOr("maxOuterIterations", ctl.maxOuterIterations));
    ctl.relax = cd.getScalarOr("relax", ctl.relax);
    ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
    ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");

    const cf::GeneralizedNewtonianControls gnCtl =
        cf::coupling::readGeneralizedNewtonianControls(caseDict);

    auto run = std::make_unique<CoupledCase>();
    run->edgeNames = edgeNames;

    // One block per coupled edge.
    for (const std::string& edgeName : edgeNames)
    {
        const cf::Dictionary& dd = multiDetails
            ? caseDict.getDict("details").getDict(edgeName)
            : caseDict.getDict("detail");

        cf::polyMesh mesh = cf::MeshLoader::fromCase(dd, casePath);
        const cf::Schemes schemes = cf::loadSchemes(dd);
        const cf::SimpleControls simple = loadSimple(dd);

        if (!dd.has("probe"))
            throw std::runtime_error("detail{} needs probe{ min max; }");
        const cf::Dictionary& pd = dd.getDict("probe");
        const double probeMin = pd.getScalar("min");
        const double probeMax = pd.getScalar("max");
        if (!(probeMin < probeMax))
            throw std::runtime_error("detail probe min must be below max");

        const cf::TurbulenceConfig turb = cf::TurbulenceConfig::fromCase(dd);
        const cf::Rheology rheology =
            cf::coupling::readPartnerRheology(caseDict, edgeName);
        const cf::coupling::PhaseConfig phases = cf::coupling::readPhases(dd);
        const cf::coupling::CompressibleConfig compressible =
            cf::coupling::readCompressible(dd);

        run->blocks.push_back(std::make_unique<cf::coupling::DetailBlock>(
            std::move(mesh), nu, rho, schemes, simple, probeMin, probeMax, turb,
            "flat", rheology, gnCtl, phases, compressible));
    }

    run->net = cf::loadNetworkFromDict(caseDict);
    for (std::size_t i = 0; i < edgeNames.size(); ++i)
    {
        cf::label edge = -1;
        for (std::size_t k = 0; k < run->net.pipes.size(); ++k)
        {
            if (run->net.pipes[k].name == edgeNames[i])
            {
                edge = static_cast<cf::label>(k);
                break;
            }
        }
        if (edge < 0)
            throw std::runtime_error("coupled: no pipe named '" + edgeNames[i]
                                     + "' in the network");
        run->edges.push_back(edge);
        run->net.pipes[static_cast<std::size_t>(edge)].coupledResistance = seedRs[i];
    }

    std::vector<cf::coupling::ICoupledElement*> partners;
    for (auto& b : run->blocks) partners.push_back(b.get());

    run->result = cf::coupling::solveCoupled1d3d(
        run->net, run->edges, partners, ctl, netCtl);
    return run;
}

py::list exchangesOf(const std::vector<cf::coupling::EdgeExchange>& src)
{
    py::list out;
    for (const auto& ex : src) out.append(py::cast(toView(ex)));
    return out;
}

} // namespace

void registerCouplingBindings(py::module_& m)
{
    py::class_<EdgeExchangeView>(m, "EdgeExchange",
        "One edge's side of one interface exchange. Field names track the "
        "kernel `EdgeExchange`; `q` is the network's delivered rate, `dpPa` "
        "the partner's reported drop, `residual` the relative cross-solver gap.")
        .def_readonly("edge", &EdgeExchangeView::edge)
        .def_readonly("name", &EdgeExchangeView::name)
        .def_readonly("resistance", &EdgeExchangeView::resistance)
        .def_readonly("q", &EdgeExchangeView::q)
        .def_readonly("dp_pa", &EdgeExchangeView::dpPa)
        .def_readonly("residual", &EdgeExchangeView::residual)
        .def_readonly("sweeps", &EdgeExchangeView::sweeps)
        .def_readonly("block_converged", &EdgeExchangeView::blockConverged)
        .def_readonly("inertance", &EdgeExchangeView::inertance)
        .def_readonly("unsteady", &EdgeExchangeView::unsteady)
        .def_readonly("dp_friction_pa", &EdgeExchangeView::dpFrictionPa)
        .def_readonly("dp_inertial_pa", &EdgeExchangeView::dpInertialPa)
        .def_readonly("q_prev", &EdgeExchangeView::qPrev)
        .def_readonly("dq", &EdgeExchangeView::dq)
        .def_readonly("row_drop_pa", &EdgeExchangeView::rowDropPa)
        .def_readonly("row_error", &EdgeExchangeView::rowError)
        .def_readonly("chord_pa", &EdgeExchangeView::chordPa)
        .def_readonly("inlet_shape", &EdgeExchangeView::inletShape)
        .def_readonly("flux_imposed", &EdgeExchangeView::fluxImposed)
        .def_readonly("beta_imposed", &EdgeExchangeView::betaImposed)
        .def_readonly("beta_in", &EdgeExchangeView::betaIn)
        .def_readonly("beta_out", &EdgeExchangeView::betaOut)
        .def_readonly("inlet_area", &EdgeExchangeView::inletArea)
        .def_readonly("p_rise_pa", &EdgeExchangeView::pRisePa)
        .def_readonly("p_rise_x", &EdgeExchangeView::pRiseX)
        .def_readonly("dp_convective_pa", &EdgeExchangeView::dpConvectivePa)
        .def_readonly("convective_row_pa", &EdgeExchangeView::convectiveRowPa)
        .def_readonly("thermal_on", &EdgeExchangeView::thermalOn)
        .def_readonly("inlet_temperature", &EdgeExchangeView::inletTemperature)
        .def_readonly("wall_heat_rate", &EdgeExchangeView::wallHeatRate)
        .def_readonly("inlet_heat_rate", &EdgeExchangeView::inletHeatRate)
        .def_readonly("enthalpy_out", &EdgeExchangeView::enthalpyOut)
        .def_readonly("line_heat_flux", &EdgeExchangeView::lineHeatFlux)
        .def_readonly("t1d", &EdgeExchangeView::t1d)
        .def_readonly("thermal_residual", &EdgeExchangeView::thermalResidual)
        .def_readonly("thermal_balance", &EdgeExchangeView::thermalBalance);

    py::class_<BlockMeasureView>(m, "BlockMeasure",
        "What one 3D partner measured on its own last solve (kernel "
        "`BlockMeasure`). `fluxIn`/`fluxOut` are the section integrals at the "
        "two ends; `probeSlope` the developed-region dp/dx.")
        .def_readonly("cells", &BlockMeasureView::cells)
        .def_readonly("sweeps", &BlockMeasureView::sweeps)
        .def_readonly("converged", &BlockMeasureView::converged)
        .def_readonly("dp_pa", &BlockMeasureView::dpPa)
        .def_readonly("imposed_q", &BlockMeasureView::imposedQ)
        .def_readonly("flux_in", &BlockMeasureView::fluxIn)
        .def_readonly("flux_out", &BlockMeasureView::fluxOut)
        .def_readonly("plane_flux_min", &BlockMeasureView::planeFluxMin)
        .def_readonly("plane_flux_max", &BlockMeasureView::planeFluxMax)
        .def_readonly("probe_slope", &BlockMeasureView::probeSlope)
        .def_readonly("probe_spread", &BlockMeasureView::probeSpread)
        .def_readonly("probe_sections", &BlockMeasureView::probeSections)
        .def_readonly("inlet_area", &BlockMeasureView::inletArea)
        .def_readonly("plug_speed", &BlockMeasureView::plugSpeed)
        .def_readonly("ux_max", &BlockMeasureView::uxMax)
        .def_readonly("ux_centre", &BlockMeasureView::uxCentre)
        .def_readonly("turb_model", &BlockMeasureView::turbModel)
        .def_readonly("nut_max", &BlockMeasureView::nutMax)
        .def_readonly("wall_y_plus_max", &BlockMeasureView::wallYPlusMax)
        .def_readonly("wall_u_tau", &BlockMeasureView::wallUTau)
        .def_readonly("inlet_shape", &BlockMeasureView::inletShape)
        .def_readonly("flux_imposed", &BlockMeasureView::fluxImposed)
        .def_readonly("beta_imposed", &BlockMeasureView::betaImposed)
        .def_readonly("beta_in", &BlockMeasureView::betaIn)
        .def_readonly("beta_out", &BlockMeasureView::betaOut)
        .def_readonly("beta_stations", &BlockMeasureView::betaStations)
        .def_readonly("pressure_stations", &BlockMeasureView::pressureStations);

    py::class_<CoupledCase, std::unique_ptr<CoupledCase>>(m, "CoupledCase",
        "Result of one steady 1D<->3D coupled run. `exchanges` is one entry "
        "per coupled edge in interface order; `pipe_flows`/`node_pressure` "
        "describe the converged 1D network; `block_measures` and the "
        "`block_velocity`/`block_pressure` arrays describe the 3D partners.")
        .def_property_readonly("converged", [](const CoupledCase& c) { return c.result.converged; })
        .def_property_readonly("iterations", [](const CoupledCase& c) { return c.result.iterations; })
        .def_property_readonly("residual", [](const CoupledCase& c) { return c.result.residual; })
        .def_property_readonly("edge_names", [](const CoupledCase& c) { return c.edgeNames; })
        .def_property_readonly("exchanges",
            [](const CoupledCase& c) { return exchangesOf(c.result.state); })
        .def_property_readonly("history",
            [](const CoupledCase& c)
            {
                py::list hist;
                for (const auto& it : c.result.history)
                {
                    py::dict d;
                    d["iteration"] = it.iteration;
                    d["edges"] = exchangesOf(it.edges);
                    hist.append(std::move(d));
                }
                return hist;
            })
        .def_property_readonly("node_pressure",
            [](const CoupledCase& c) { return c.result.network.nodePressure; })
        .def_property_readonly("pipe_flows",
            [](const CoupledCase& c) { return c.result.network.pipeFlows; })
        .def_property_readonly("block_measures",
            [](const CoupledCase& c)
            {
                py::list out;
                for (const auto& b : c.blocks)
                    out.append(py::cast(toView(b->measure())));
                return out;
            })
        .def("block_velocity",
            [](const CoupledCase& c, std::size_t i)
            {
                if (i >= c.blocks.size())
                    throw std::out_of_range("block index out of range");
                return cfbind::vectorFieldToNumpy(c.blocks[i]->U().internalField());
            },
            py::arg("edge_index"),
            "Cell-centred partner velocity, shape (n, 3), float64 C-order copy.")
        .def("block_pressure",
            [](const CoupledCase& c, std::size_t i)
            {
                if (i >= c.blocks.size())
                    throw std::out_of_range("block index out of range");
                return cfbind::scalarFieldToNumpy(c.blocks[i]->p().internalField());
            },
            py::arg("edge_index"),
            "Cell-centred partner pressure, shape (n,), float64 C-order copy.");

    m.def("solve_coupled_case", &runCoupled,
          py::arg("case_dict_path"), py::arg("seed_index") = 0,
          "Read a coupled case.dict and run the steady 1D<->3D interface "
          "solve. `seed_index` selects which starting coefficient set runs "
          "(single `seeds` list index, or multi `starts` set order).");
}
