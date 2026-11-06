// chflow-coupled-1d3d-multi: a steady 1D pipe network with several 3D detailed
// regions embedded in it. Each region is a `coupledResistance` pipe on the 1D
// side and a SIMPLE block on the 3D side; one interface pass solves the network
// once and asks every block for the drop at the rate that pass delivers through
// its own edge, so the edges are refreshed together and the pass stops on the
// worst of them.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/coupling/CoupledNetwork.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;
using namespace cf::coupling;

namespace
{
// log::info() stringifies with ostringstream defaults (6 digits); the checker
// asserts cross-solver identities, so emit every reported number in full.
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

// simpleFoam reads its own simple{} block and no other app needs these keys,
// so the coupling app parses its own copy rather than a shared loader.
SimpleControls loadSimple(const Dictionary& d)
{
    SimpleControls ctl;
    if (!d.has("simple")) return ctl;
    const Dictionary& sd = d.getDict("simple");
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

/// One block's own dictionary: the mesh, the discretisation, the SIMPLE
/// controls, the probe window the developed gradient is measured over, and the
/// file its field is written to.
struct BlockSpec
{
    polyMesh mesh;
    Schemes schemes;
    SimpleControls simple;
    scalar probeMin = 0;
    scalar probeMax = 0;
    std::string vtk;
    TurbulenceConfig turb;
};

BlockSpec loadBlockSpec(const Dictionary& dd, const char* argv0,
                        const std::string& edgeName)
{
    if (!dd.has("probe"))
        cfFatal("coupled1d3d: details." + edgeName + "{} needs probe{ min max; },"
                " the cross-section window the developed gradient is measured over");
    const Dictionary& pd = dd.getDict("probe");
    const scalar probeMin = pd.getScalar("min");
    const scalar probeMax = pd.getScalar("max");
    if (!(probeMin < probeMax))
        cfFatal("coupled1d3d: details." + edgeName + " probe min must be below max");
    return BlockSpec{MeshLoader::fromCase(dd, argv0), loadSchemes(dd), loadSimple(dd),
                     probeMin, probeMax, dd.getWordOr("vtk", ""),
                     TurbulenceConfig::fromCase(dd)};
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-coupled-1d3d-multi <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        // ---- the 1D side ----
        const Dictionary& fd = caseDict.getDict("fluid");
        const double rho = fd.getScalar("rho");
        const double mu = fd.getScalar("mu");
        if (rho <= 0 || mu <= 0)
            cfFatal("coupled1d3d: fluid rho and mu must both be positive");
        const double nu = mu / rho;

        const SteadyControls netCtl = loadSteadyControls(caseDict);

        const Dictionary& cd = caseDict.getDict("couple");
        const std::vector<std::string> edgeNames = cd.getWordList("edges");
        if (edgeNames.empty())
            cfFatal("coupled1d3d: couple{} needs the coupled pipe names in"
                    " edges ( )");
        CoupledControls ctl;
        ctl.maxOuterIterations = static_cast<label>(
            cd.getScalarOr("maxOuterIterations", ctl.maxOuterIterations));
        ctl.relax = cd.getScalarOr("relax", ctl.relax);
        ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
        ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");

        // Each named set is one whole-network run, and it quotes a starting
        // coefficient for every coupled edge at once: a single edge starting far
        // off must not be able to pull its neighbour's answer along with it.
        if (!cd.has("starts"))
            cfFatal("coupled1d3d: couple{} needs starts{ name { <edge> <R>; ... } },"
                    " one whole-network starting coefficient set per run");
        const Dictionary& startsDict = cd.getDict("starts");
        std::vector<std::string> startNames;
        std::vector<std::vector<double>> startRs;
        for (const std::string& sname : startsDict.keys())
        {
            const Dictionary& sd = startsDict.getDict(sname);
            std::vector<double> rs;
            for (const std::string& edgeName : edgeNames)
            {
                if (!sd.has(edgeName))
                    cfFatal("coupled1d3d: starts." + sname + " gives no starting"
                            " resistance for coupled edge '" + edgeName + "'");
                rs.push_back(sd.getScalar(edgeName));
            }
            startNames.push_back(sname);
            startRs.push_back(std::move(rs));
        }
        if (startNames.empty())
            cfFatal("coupled1d3d: couple.starts{} needs at least one set");

        // Which interface operator each starting set is run through is a per-run
        // choice, so one case can send the same network and the same blocks
        // through plain Gauss-Seidel and through the accelerated one and the two
        // arms differ in the loop and in nothing else.
        std::vector<bool> accelerates(startNames.size(), ctl.accelerate);
        if (cd.has("accelerates"))
        {
            const std::vector<double> a = cd.getScalarList("accelerates");
            if (a.size() == 1)
                std::fill(accelerates.begin(), accelerates.end(), a[0] != 0.0);
            else if (a.size() != accelerates.size())
                cfFatal("coupled1d3d: couple.accelerates needs one entry per"
                        " starting set, or exactly one for all of them (got "
                        + std::to_string(a.size()) + " for "
                        + std::to_string(accelerates.size()) + " sets)");
            else
                for (std::size_t i = 0; i < a.size(); ++i)
                    accelerates[i] = a[i] != 0.0;
        }

        // ---- the 3D side: one block per coupled edge ----
        const Dictionary& allDetails = caseDict.getDict("details");
        std::vector<BlockSpec> specs;
        for (const std::string& edgeName : edgeNames)
        {
            if (!allDetails.has(edgeName))
                cfFatal("coupled1d3d: edge '" + edgeName + "' has no details."
                        + edgeName + "{} block");
            specs.push_back(loadBlockSpec(allDetails.getDict(edgeName), argv[1],
                                          edgeName));
        }
        for (std::size_t i = 0; i < specs.size(); ++i)
        {
            log::info("block ", edgeNames[i], ": ", specs[i].mesh.nCells(), " cells,",
                      " ", specs[i].mesh.nPatches(), " patches, probe ",
                      num(specs[i].probeMin), " to ", num(specs[i].probeMax),
                      ", turbulence ", specs[i].turb.model);
        }
        log::info("network: ", edgeNames.size(), " coupled edges, ",
                  num(startRs.size()), " starting sets; nu ", num(nu));

        bool allConverged = true;
        for (std::size_t s = 0; s < startNames.size(); ++s)
        {
            const std::string& startName = startNames[s];

            // one fresh set of blocks per starting set: each block warm-starts
            // from its own previous solution, so sharing them between starting
            // sets would let run 2 inherit run 1's interface history
            std::vector<std::unique_ptr<DetailBlock>> owners;
            std::vector<ICoupledElement*> partners;
            for (std::size_t i = 0; i < specs.size(); ++i)
            {
                owners.push_back(std::make_unique<DetailBlock>(
                    specs[i].mesh, nu, rho, specs[i].schemes, specs[i].simple,
                    specs[i].probeMin, specs[i].probeMax, specs[i].turb));
                partners.push_back(owners.back().get());
            }

            // and one fresh network: the loop writes the converged coefficients
            // into the coupled pipes, so a reused instance would carry them as
            // starting values and the runs would no longer be independent
            OneDNetwork net = loadNetworkFromDict(caseDict);
            std::vector<label> edges;
            for (const std::string& edgeName : edgeNames)
            {
                label edge = -1;
                for (std::size_t i = 0; i < net.pipes.size(); ++i)
                {
                    if (net.pipes[i].name == edgeName) edge = static_cast<label>(i);
                }
                if (edge < 0)
                    cfFatal("coupled1d3d: no pipe named '" + edgeName
                            + "' in the network");
                edges.push_back(edge);
            }
            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                net.pipes[static_cast<std::size_t>(edges[i])].coupledResistance =
                    startRs[s][i];
                log::info("start ", startName, " edge ", edgeNames[i],
                          " resistance ", num(startRs[s][i]));
            }

            CoupledControls runCtl = ctl;
            runCtl.accelerate = accelerates[s];
            log::info("start ", startName, " interface ",
                      runCtl.accelerate ? "quasi-Newton" : "plain Gauss-Seidel");

            const CoupledResult res =
                solveCoupled1d3d(net, edges, partners, runCtl, netCtl);

            for (const OuterIteration& it : res.history)
            {
                for (std::size_t i = 0; i < it.edges.size(); ++i)
                {
                    const EdgeExchange& ex = it.edges[i];
                    log::info("interface ", startName, " iter ", it.iteration,
                              " edge ", edgeNames[i], " R ", num(ex.resistance),
                              " q ", num(ex.q), " dpBlock ", num(ex.dpPa),
                              " sweeps ", ex.sweeps, " blockConverged ",
                              (ex.blockConverged ? 1 : 0),
                              " residual ", num(ex.residual));
                }
            }

            for (std::size_t i = 0; i < res.state.size(); ++i)
            {
                const EdgeExchange& st = res.state[i];
                log::info("settled ", startName, " edge ", edgeNames[i],
                          " q ", num(st.q), " dp ", num(st.dpPa),
                          " R ", num(st.resistance), " residual ", num(st.residual),
                          " sweeps ", st.sweeps);
            }
            log::info("coupled ", startName, ": iterations ", res.iterations,
                      " worstResidual ", num(res.residual),
                      res.converged ? " (converged)" : " (NOT converged)",
                      " accelerate ", (runCtl.accelerate ? 1 : 0),
                      " blockSolves ", res.iterations * res.edges.size());

            // name-keyed so a checker reads the network it got rather than the
            // one it expects: the coupled edges' own drops come from the solve
            std::ostringstream dpNet;
            dpNet << std::setprecision(17);
            for (std::size_t i = 0; i < res.state.size(); ++i)
            {
                const PipeFlowResult& pf =
                    res.network.pipeFlows[static_cast<std::size_t>(edges[i])];
                dpNet << " " << edgeNames[i] << " " << pf.dpTotal
                      << " " << pf.flowRate;
            }
            log::info("coupledEdges ", startName, " dp q:", dpNet.str());

            std::ostringstream nodes;
            nodes << std::setprecision(17);
            for (std::size_t i = 0; i < res.network.nodePressure.size(); ++i)
                nodes << " " << net.nodeNames[i] << " "
                      << res.network.nodePressure[i];
            log::info("nodes ", startName, ":", nodes.str());

            std::ostringstream pipes;
            pipes << std::setprecision(17);
            for (std::size_t i = 0; i < res.network.pipeFlows.size(); ++i)
                pipes << " " << net.pipes[i].name << " "
                      << res.network.pipeFlows[i].flowRate;
            log::info("pipes ", startName, ":", pipes.str());

            // each block's measure describes its own last solve, so it is
            // reported per edge and per starting set
            for (std::size_t i = 0; i < owners.size(); ++i)
            {
                const BlockMeasure& m = owners[i]->measure();
                log::info("block ", startName, " edge ", edgeNames[i],
                          ": cells ", m.cells, " sweeps ", m.sweeps,
                          " blockConverged ", (m.converged ? 1 : 0),
                          " fluxIn ", num(m.fluxIn), " fluxOut ", num(m.fluxOut),
                          " planeFluxMin ", num(m.planeFluxMin),
                          " planeFluxMax ", num(m.planeFluxMax),
                          " inletArea ", num(m.inletArea), " plugSpeed ",
                          num(m.plugSpeed), " uxMax ", num(m.uxMax), " uxCentre ",
                          num(m.uxCentre), " probeSections ", m.probeSections,
                          " probeSlope ", num(m.probeSlope), " probeSpread ",
                          num(m.probeSpread), " dpPa ", num(m.dpPa), " turb ",
                          (m.turbModel.empty() ? std::string("laminar") : m.turbModel),
                          " nutMax ", num(m.nutMax), " wallYPlusMax ",
                          num(m.wallYPlusMax), " wallUTau ", num(m.wallUTau));
                if (s + 1 == startNames.size() && !specs[i].vtk.empty())
                {
                    VtkWriter w(owners[i]->mesh());
                    w.addCellVector("U", owners[i]->U().internalField());
                    w.addCellScalar("p", owners[i]->p().internalField());
                    w.write(specs[i].vtk);
                }
            }
            allConverged = allConverged && res.converged;
        }

        return allConverged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("coupled1d3d: ", e.what());
        return 1;
    }
}
