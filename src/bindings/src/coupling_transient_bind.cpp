// Transient 1D<->3D coupling family: read a time-stepping coupled case.dict,
// march the partitioned interface solve, and hand the per-step structured
// result back to Python.
//
// The construction path mirrors src/coupling/apps/transientCoupledApp.cpp
// (single edge, `detail{}`) and transientCoupledMultiApp.cpp (one `details{}`
// block per edge, `starts{}` sets). Both layouts are accepted:
//   * single edge:  couple { edge X; seeds ( ... ); }  detail { ... }
//   * multi edge:   couple { edges ( X1 X2 ); starts { ... }; }  details { ... }
//
// Every committed time level is copied out of the kernel
// `TransientCoupledStep` (edge exchanges, node pressures, pipe flows) so Python
// owns the numbers; the last accepted partner blocks stay alive for the U/p
// field read-back.
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/1d/Transient.h"
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

/// One block's own ingredients, built once per edge before the selected arm
/// runs. The app copies the mesh into a fresh block per arm; the binding keeps
/// the spec and copies the mesh the same way when it constructs its block.
struct BlockSpec
{
    cf::polyMesh mesh;
    cf::Schemes schemes;
    cf::SimpleControls simple;
    double probeMin = 0;
    double probeMax = 0;
    cf::TurbulenceConfig turb;
    cf::Rheology rheology;
    cf::coupling::PhaseConfig phases;
    cf::coupling::CompressibleConfig compressible;
    cf::coupling::ThermalCouplingConfig thermal;
};

/// Per-arm (per starting set) choices parsed straight off the couple{} block,
/// mirroring the two apps' own defaults for every missing key.
struct ArmSpec
{
    std::vector<double> startRs;       // one starting resistance per edge
    std::vector<std::string> shapes;   // inlet shape per edge
    bool accelerate = false;
    std::string inertia = "network";  // network | partner | slab
    std::string ddt;                   // "" -> the detail's own schemes
    std::vector<std::string> edgeInertias;   // per-edge source words, may be empty
    std::vector<std::string> edgeDdts;       // per-edge ddt words, may be empty
    std::vector<int> edgeConvective;          // 0 report / 1 charge, may be empty
};

int inertiaSourceIndex(const std::string& w)
{
    if (w == "network") return 0;
    if (w == "partner") return 1;
    if (w == "slab") return 2;
    throw std::runtime_error("coupled transient: inertia source '" + w
                              + "' is none of network, partner, slab");
}

/// Owns the case, the partner blocks and the result of one transient run.
struct TransientCoupledCase
{
    cf::TransientCase baseCase;
    std::vector<std::string> edgeNames;
    std::vector<std::unique_ptr<cf::coupling::DetailBlock>> blocks;
    cf::coupling::TransientCoupledResult result;
};

struct StepView
{
    double t = 0;
    cf::label passes = 0;
    std::vector<cfbind::EdgeExchangeView> edges;
    std::vector<double> pressure;
    std::vector<double> flowRate;
    double thermalResidual = 0;
};

std::unique_ptr<TransientCoupledCase> runTransientCoupled(const std::string& casePath,
                                                          std::size_t seedIndex)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);

    // ---- the 1D side ----
    const cf::TransientCase c = cf::loadTransientFromDict(caseDict);
    if (c.fluid.rho <= 0 || c.fluid.mu <= 0)
        throw std::runtime_error("coupled transient: fluid rho and mu must both be positive");
    const cf::TransientTimeControls tc = cf::loadTransientControls(caseDict);
    const double nu = c.fluid.mu / c.fluid.rho;

    const cf::Dictionary& cd = caseDict.getDict("couple");
    const bool multiLayout = cd.has("edges");

    std::vector<std::string> edgeNames;
    if (multiLayout)
        edgeNames = cd.getWordList("edges");
    else
        edgeNames.push_back(cd.getWord("edge"));
    if (edgeNames.empty())
        throw std::runtime_error("couple{} names no coupled edge");

    cf::coupling::TransientCoupledControls ctl;
    ctl.maxSubIterations = static_cast<cf::label>(
        cd.getScalarOr("maxSubIterations", ctl.maxSubIterations));
    ctl.relax = cd.getScalarOr("relax", ctl.relax);
    ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
    ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");

    // ---- starting sets: `starts{}` (multi) or `seeds ( )` (single) ----
    std::vector<ArmSpec> arms;
    if (cd.has("starts"))
    {
        const cf::Dictionary& startsDict = cd.getDict("starts");
        const std::vector<std::string> startNames = startsDict.keys();
        if (startNames.empty())
            throw std::runtime_error("couple.starts{} declares no starting set");

        std::vector<bool> accelerates(startNames.size(), ctl.accelerate);
        if (cd.has("accelerates"))
        {
            const std::vector<double> a = cd.getScalarList("accelerates");
            if (a.size() == 1)
                std::fill(accelerates.begin(), accelerates.end(), a[0] != 0.0);
            else if (a.size() != accelerates.size())
                throw std::runtime_error("couple.accelerates needs one entry per"
                                         " starting set, or exactly one for all");
            else
                for (std::size_t i = 0; i < a.size(); ++i)
                    accelerates[i] = a[i] != 0.0;
        }

        std::vector<std::string> inertias(startNames.size(), "network");
        if (cd.has("inertias"))
        {
            const std::vector<std::string> w = cd.getWordList("inertias");
            if (w.size() == 1 && startNames.size() > 1)
                inertias.assign(startNames.size(), w[0]);
            else if (w.size() != inertias.size())
                throw std::runtime_error("couple.inertias needs one entry per"
                                         " starting set, or exactly one for all");
            else
                inertias = w;
        }
        for (const std::string& w : inertias) inertiaSourceIndex(w);

        std::vector<std::string> ddts(startNames.size(), "");
        if (cd.has("ddts"))
        {
            const std::vector<std::string> w = cd.getWordList("ddts");
            if (w.size() == 1 && startNames.size() > 1)
                ddts.assign(startNames.size(), w[0]);
            else if (w.size() != ddts.size())
                throw std::runtime_error("couple.ddts needs one entry per starting"
                                         " set, or exactly one for all");
            else
                ddts = w;
            for (const std::string& w : ddts)
                if (!w.empty() && w != "steady" && w != "Euler")
                    throw std::runtime_error("couple.ddts entry '" + w
                                              + "' is neither steady nor Euler");
        }

        const cf::Dictionary* ddtFrom = cd.has("ddtFrom") ? &cd.getDict("ddtFrom") : nullptr;
        const cf::Dictionary* inertiaFrom = cd.has("inertiaFrom") ? &cd.getDict("inertiaFrom") : nullptr;
        const cf::Dictionary* convectiveToRow = cd.has("convectiveToRow") ? &cd.getDict("convectiveToRow") : nullptr;

        for (std::size_t s = 0; s < startNames.size(); ++s)
        {
            const std::string& sname = startNames[s];
            const cf::Dictionary& sd = startsDict.getDict(sname);
            ArmSpec arm;
            arm.accelerate = accelerates[s];
            arm.inertia = inertias[s];
            arm.ddt = ddts[s];
            for (const std::string& edgeName : edgeNames)
            {
                if (!sd.has(edgeName))
                    throw std::runtime_error("couple.starts." + sname + " gives no"
                                             " starting resistance for edge '"
                                             + edgeName + "'");
                arm.startRs.push_back(sd.getScalar(edgeName));
            }
            arm.shapes.assign(edgeNames.size(), "flat");
            if (sd.has("profile"))
            {
                const cf::Dictionary::Kind kind = sd.find("profile").kind;
                if (kind == cf::Dictionary::Kind::word)
                    std::fill(arm.shapes.begin(), arm.shapes.end(), sd.getWord("profile"));
                else if (kind == cf::Dictionary::Kind::list)
                {
                    arm.shapes = sd.getWordList("profile");
                    if (arm.shapes.size() != edgeNames.size())
                        throw std::runtime_error("couple.starts." + sname
                            + " profile names " + std::to_string(arm.shapes.size())
                            + " shapes for " + std::to_string(edgeNames.size())
                            + " coupled edges");
                }
                else
                    throw std::runtime_error("couple.starts." + sname
                        + " profile must be one shape word or a list of them");
            }
            if (ddtFrom && ddtFrom->has(sname))
            {
                const cf::Dictionary& ed = ddtFrom->getDict(sname);
                for (const std::string& edgeName : edgeNames)
                {
                    const std::string e = ed.getWord(edgeName);
                    if (e != "steady" && e != "Euler")
                        throw std::runtime_error("ddtFrom." + sname + "." + edgeName
                            + " is '" + e + "', neither steady nor Euler");
                    arm.edgeDdts.push_back(e);
                }
            }
            if (inertiaFrom && inertiaFrom->has(sname))
            {
                const cf::Dictionary& ed = inertiaFrom->getDict(sname);
                for (const std::string& edgeName : edgeNames)
                {
                    const std::string e = ed.getWord(edgeName);
                    inertiaSourceIndex(e);
                    arm.edgeInertias.push_back(e);
                }
            }
            if (convectiveToRow && convectiveToRow->has(sname))
            {
                const cf::Dictionary& ed = convectiveToRow->getDict(sname);
                for (const std::string& edgeName : edgeNames)
                {
                    const std::string e = ed.getWord(edgeName);
                    if (e == "report") arm.edgeConvective.push_back(0);
                    else if (e == "charge") arm.edgeConvective.push_back(1);
                    else
                        throw std::runtime_error("convectiveToRow." + sname + "."
                            + edgeName + " is '" + e + "', neither report nor charge");
                }
            }
            arms.push_back(std::move(arm));
        }
    }
    else
    {
        const std::vector<double> seeds = cd.getScalarList("seeds");
        if (seeds.empty())
            throw std::runtime_error("couple{} needs seeds ( ) or starts{}");
        std::vector<std::string> inertias(seeds.size(), "network");
        if (cd.has("inertias"))
        {
            const std::vector<std::string> w = cd.getWordList("inertias");
            if (w.size() == 1 && seeds.size() > 1)
                inertias.assign(seeds.size(), w[0]);
            else if (w.size() != inertias.size())
                throw std::runtime_error("couple.inertias needs one entry per seed,"
                                         " or one entry for all of them");
            else
                inertias = w;
        }
        for (const std::string& w : inertias) inertiaSourceIndex(w);

        for (std::size_t s = 0; s < seeds.size(); ++s)
        {
            ArmSpec arm;
            arm.accelerate = ctl.accelerate;
            arm.inertia = inertias[s];
            arm.startRs.assign(edgeNames.size(), seeds[s]);
            arm.shapes.assign(edgeNames.size(), "flat");
            arms.push_back(std::move(arm));
        }
    }

    if (seedIndex >= arms.size())
        throw std::out_of_range("coupled transient: seed index "
            + std::to_string(seedIndex) + " out of range ("
            + std::to_string(arms.size()) + " starting sets)");
    const ArmSpec& arm = arms[seedIndex];

    // ---- the 3D side: one spec per coupled edge ----
    const bool multiDetails = caseDict.has("details");
    const cf::GeneralizedNewtonianControls gnCtl =
        cf::coupling::readGeneralizedNewtonianControls(caseDict);

    std::vector<BlockSpec> specs;
    for (const std::string& edgeName : edgeNames)
    {
        const cf::Dictionary& dd = multiDetails
            ? caseDict.getDict("details").getDict(edgeName)
            : caseDict.getDict("detail");

        BlockSpec spec{cf::MeshLoader::fromCase(dd, casePath)};
        spec.schemes = cf::loadSchemes(dd);
        spec.simple = loadSimple(dd);
        if (!dd.has("probe"))
            throw std::runtime_error("detail{} needs probe{ min max; }");
        const cf::Dictionary& pd = dd.getDict("probe");
        spec.probeMin = pd.getScalar("min");
        spec.probeMax = pd.getScalar("max");
        if (!(spec.probeMin < spec.probeMax))
            throw std::runtime_error("detail probe min must be below max");
        spec.turb = cf::TurbulenceConfig::fromCase(dd);
        spec.rheology = cf::coupling::readPartnerRheology(caseDict, edgeName);
        spec.phases = cf::coupling::readPhases(dd);
        spec.compressible = cf::coupling::readCompressible(dd);
        spec.thermal = cf::coupling::readThermalCoupling(dd);
        specs.push_back(std::move(spec));
    }

    auto run = std::make_unique<TransientCoupledCase>();
    run->baseCase = c;
    run->edgeNames = edgeNames;

    cf::coupling::TransientCoupledControls runCtl = ctl;
    runCtl.accelerate = arm.accelerate;
    runCtl.inertiaFromPartner = (arm.inertia == "partner");
    if (!arm.edgeInertias.empty())
    {
        std::vector<int> src;
        for (const std::string& w : arm.edgeInertias) src.push_back(inertiaSourceIndex(w));
        runCtl.inertiaFromPartnerPerEdge = std::move(src);
    }
    else if (arm.inertia == "slab")
    {
        runCtl.inertiaFromPartnerPerEdge = std::vector<int>(edgeNames.size(), 2);
    }
    if (!arm.edgeConvective.empty())
        runCtl.convectiveIntoRowPerEdge = arm.edgeConvective;
    // The single-edge app installs the conjugate-heat pass on the interface;
    // the multi-edge app builds hydraulic-only blocks.
    if (!multiLayout)
        runCtl.thermal = specs[0].thermal;

    std::vector<cf::label> edges;
    cf::TransientNetwork net(c);
    for (std::size_t i = 0; i < edgeNames.size(); ++i)
    {
        cf::label edge = -1;
        const std::vector<std::string> names = net.pipeNames();
        for (std::size_t k = 0; k < names.size(); ++k)
            if (names[k] == edgeNames[i]) { edge = static_cast<cf::label>(k); break; }
        if (edge < 0)
            throw std::runtime_error("coupled transient: no pipe named '" + edgeNames[i]
                                      + "' in the network");
        edges.push_back(edge);
        net.setCoupledResistance(edge, arm.startRs[i]);
    }

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        cf::Schemes sch = specs[i].schemes;
        if (!arm.ddt.empty()) sch.ddt = arm.ddt;
        if (!arm.edgeDdts.empty()) sch.ddt = arm.edgeDdts[i];
        const cf::coupling::ThermalCouplingConfig thermal =
            multiLayout ? cf::coupling::ThermalCouplingConfig{} : specs[i].thermal;
        run->blocks.push_back(std::make_unique<cf::coupling::DetailBlock>(
            cf::polyMesh(specs[i].mesh), nu, c.fluid.rho, sch, specs[i].simple,
            specs[i].probeMin, specs[i].probeMax, specs[i].turb, arm.shapes[i],
            specs[i].rheology, gnCtl, specs[i].phases, specs[i].compressible, thermal));
    }

    std::vector<cf::coupling::ICoupledElement*> partners;
    for (auto& b : run->blocks) partners.push_back(b.get());

    run->result = cf::coupling::solveTransientCoupled1d3d(
        net, edges, partners, runCtl, tc.dt, tc.endTime);
    return run;
}

py::list exchangesOf(const std::vector<cf::coupling::EdgeExchange>& src)
{
    py::list out;
    for (const auto& ex : src) out.append(py::cast(cfbind::toView(ex)));
    return out;
}

} // namespace

void registerCouplingTransientBindings(py::module_& m)
{
    py::class_<StepView>(m, "TransientCoupledStep",
        "One committed time level of a transient coupled run. `edges` is one"
        " exchange per coupled edge in interface order; `pressure`/`flow_rate`"
        " are the whole 1D network state at that level.")
        .def_readonly("t", &StepView::t)
        .def_readonly("passes", &StepView::passes)
        .def_readonly("edges", &StepView::edges)
        .def_readonly("pressure", &StepView::pressure)
        .def_readonly("flow_rate", &StepView::flowRate)
        .def_readonly("thermal_residual", &StepView::thermalResidual);

    py::class_<TransientCoupledCase, std::unique_ptr<TransientCoupledCase>>(
        m, "TransientCoupledCase",
        "Result of one transient 1D<->3D coupled run. `steps` holds every"
        " committed time level; `state` the last level's per-edge exchanges;"
        " `block_measures` and the `block_velocity`/`block_pressure` arrays"
        " describe the partner blocks at the last accepted rate.")
        .def_property_readonly("edge_names",
            [](const TransientCoupledCase& c) { return c.edgeNames; })
        .def_property_readonly("node_names",
            [](const TransientCoupledCase& c) { return c.result.nodeNames; })
        .def_property_readonly("pipe_names",
            [](const TransientCoupledCase& c) { return c.result.pipeNames; })
        .def_property_readonly("steps",
            [](const TransientCoupledCase& c)
            {
                py::list out;
                for (const auto& st : c.result.steps)
                {
                    StepView v;
                    v.t = st.t;
                    v.passes = st.passes;
                    v.pressure = st.pressure;
                    v.flowRate = st.flowRate;
                    v.thermalResidual = st.thermalResidual;
                    for (const auto& ex : st.edges)
                        v.edges.push_back(cfbind::toView(ex));
                    out.append(py::cast(std::move(v)));
                }
                return out;
            })
        .def_property_readonly("state",
            [](const TransientCoupledCase& c) { return exchangesOf(c.result.state); })
        .def_property_readonly("residual",
            [](const TransientCoupledCase& c) { return c.result.residual; })
        .def_property_readonly("block_measures",
            [](const TransientCoupledCase& c)
            {
                py::list out;
                for (const auto& b : c.blocks)
                    out.append(py::cast(cfbind::toView(b->measure())));
                return out;
            })
        .def("block_velocity",
            [](const TransientCoupledCase& c, std::size_t i)
            {
                if (i >= c.blocks.size())
                    throw std::out_of_range("block index out of range");
                return cfbind::vectorFieldToNumpy(c.blocks[i]->U().internalField());
            },
            py::arg("edge_index"),
            "Cell-centred partner velocity, shape (n, 3), float64 C-order copy.")
        .def("block_pressure",
            [](const TransientCoupledCase& c, std::size_t i)
            {
                if (i >= c.blocks.size())
                    throw std::out_of_range("block index out of range");
                return cfbind::scalarFieldToNumpy(c.blocks[i]->p().internalField());
            },
            py::arg("edge_index"),
            "Cell-centred partner pressure, shape (n,), float64 C-order copy.");

    m.def("solve_transient_coupled_case", &runTransientCoupled,
          py::arg("case_dict_path"), py::arg("seed_index") = 0,
          "Read a transient coupled case.dict and march the 1D<->3D interface"
          " solve. `seed_index` selects the starting coefficient set (single"
          " `seeds` list index, or multi `starts` set order).");
}
