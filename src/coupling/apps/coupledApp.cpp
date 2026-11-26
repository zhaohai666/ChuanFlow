// chflow-coupled-1d3d: a steady 1D pipe network with one or more 3D detailed
// regions embedded in it. Each detailed region is a `coupledResistance` pipe
// on the 1D side and a SIMPLE block on the 3D side; the two exchange
// rate-for-drop until the network charges the drop the block reports.
//
// The legacy single-edge declaration is `couple { edge X; ... }` with one
// `detail { ... }` block. A second, multi-edge declaration `couple { edges
// ( X1 X2 ); ... }` pairs each named pipe with its own `details { X1 { ... }
// }` block; every edge is solved on the same network pass and reported on its
// own `coupled solved: edge ...` line.
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
/// controls, the probe window the developed gradient is measured over, the
/// fluid laws, and the file its field is written to.
struct BlockSpec
{
    polyMesh mesh;
    Schemes schemes;
    SimpleControls simple;
    scalar probeMin = -std::numeric_limits<scalar>::infinity();
    scalar probeMax = probeMin;
    TurbulenceConfig turb;
    Rheology rheology;
    PhaseConfig phases;
    CompressibleConfig compressible;
    std::string vtk;
};

BlockSpec loadBlockSpec(const Dictionary& dd, const char* argv0,
                        const Dictionary& caseDict, const std::string& edgeName)
{
    const TurbulenceConfig turb = TurbulenceConfig::fromCase(dd);
    const Rheology rheology = readPartnerRheology(caseDict, edgeName);
    const PhaseConfig phases = readPhases(dd);
    const CompressibleConfig compressible = readCompressible(dd);
    scalar probeMin = -std::numeric_limits<scalar>::infinity();
    scalar probeMax = probeMin;
    if (dd.has("probe"))
    {
        const Dictionary& pd = dd.getDict("probe");
        probeMin = pd.getScalar("min");
        probeMax = pd.getScalar("max");
        if (!(probeMin < probeMax))
            cfFatal("coupled1d3d: detail probe min must be below max");
    }
    return BlockSpec{
        MeshLoader::fromCase(dd, argv0),
        loadSchemes(dd),
        loadSimple(dd),
        probeMin,
        probeMax,
        turb,
        rheology,
        phases,
        compressible,
        dd.getWordOr("vtk", ""),
    };
}

/// Write the block's own solved fields besides U/p. Only fields the solver
/// actually built are added: T from the energy model, alpha from the frozen
/// two-phase field, nut from the RANS model. A case that declared none of
/// them keeps its old U/p-only VTU byte-for-byte.
void writeBlockFields(VtkWriter& w, const DetailBlock& block)
{
    if (const EnergyModel* energy = block.energy())
        w.addCellScalar("T", energy->temperature().internalField());
    if (const PhaseField* ph = block.phases())
        w.addCellScalar("alpha", ph->alphaCell());
    if (const RansModel* turb = block.turbulence())
        w.addCellScalar("nut", turb->nut());
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-coupled-1d3d <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        // ---- the 1D side ----
        const FluidProperties fluid = [&caseDict]
        {
            const Dictionary& fd = caseDict.getDict("fluid");
            return FluidProperties{fd.getScalar("rho"), fd.getScalar("mu")};
        }();
        if (fluid.rho <= 0 || fluid.mu <= 0)
            cfFatal("coupled1d3d: fluid rho and mu must both be positive");

        const SteadyControls netCtl = loadSteadyControls(caseDict);

        const Dictionary& cd = caseDict.getDict("couple");
        // Backward compatible: the legacy `edge <name>;` word, plus the new
        // `edges ( N1 N2 );` list. One edge keeps the old single-block path
        // byte-for-byte; several edges pair one block per edge.
        std::vector<std::string> edgeNames;
        if (cd.has("edges"))
            edgeNames = cd.getWordList("edges");
        else
            edgeNames.push_back(cd.getWord("edge"));
        if (edgeNames.empty())
            cfFatal("coupled1d3d: couple{} needs `edge <name>;` or"
                    " `edges ( <name> ... );`");
        const bool multiEdge = edgeNames.size() > 1;

        CoupledControls ctl;
        ctl.maxOuterIterations = static_cast<label>(cd.getScalarOr("maxOuterIterations", ctl.maxOuterIterations));
        ctl.relax = cd.getScalarOr("relax", ctl.relax);
        ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");
        ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
        // the starting interface coefficients: each one runs the whole
        // partitioned solve on its own, so a seed that lands anywhere else is
        // not a coupling but a fixed point the initial value still owns
        std::vector<double> seeds;
        if (cd.has("seeds")) seeds = cd.getScalarList("seeds");
        if (seeds.empty())
            cfFatal("coupled1d3d: couple{} needs at least one starting"
                    " resistance in seeds ( )");

        // Per-seed interface operator: the same seed run once plainly and once
        // accelerated differs only in how the next coefficient is chosen, so
        // the two pass counts say what the quasi-Newton did to this block.
        std::vector<bool> accelerates(seeds.size(), ctl.accelerate);
        if (cd.has("accelerates"))
        {
            const std::vector<double> a = cd.getScalarList("accelerates");
            if (a.size() == 1)
                std::fill(accelerates.begin(), accelerates.end(), a[0] != 0.0);
            else if (a.size() != seeds.size())
                cfFatal("coupled1d3d: couple.accelerates needs one entry per"
                        " seed (" + std::to_string(seeds.size()) + " seeds, got "
                        + std::to_string(a.size()) + ")");
            else
                for (std::size_t i = 0; i < a.size(); ++i)
                    accelerates[i] = a[i] != 0.0;
        }

        // ---- the 3D side: one block spec per coupled edge ----
        const GeneralizedNewtonianControls gnCtl =
            readGeneralizedNewtonianControls(caseDict);
        std::vector<BlockSpec> specs;
        if (caseDict.has("details"))
        {
            const Dictionary& allDetails = caseDict.getDict("details");
            for (const std::string& edgeName : edgeNames)
            {
                if (!allDetails.has(edgeName))
                    cfFatal("coupled1d3d: edge '" + edgeName + "' has no"
                            " details." + edgeName + "{} block");
                specs.push_back(loadBlockSpec(allDetails.getDict(edgeName),
                                              argv[1], caseDict, edgeName));
            }
        }
        else
        {
            // Legacy single `detail {}` block; only valid for one edge.
            const Dictionary& dd = caseDict.getDict("detail");
            if (multiEdge)
                cfFatal("coupled1d3d: several coupled edges need a"
                        " details{} block per edge");
            specs.push_back(loadBlockSpec(dd, argv[1], caseDict, edgeNames[0]));
        }

        for (std::size_t i = 0; i < specs.size(); ++i)
        {
            log::info("network+block: ", edgeNames[i], " block has ",
                      specs[i].mesh.nCells(), " cells and ",
                      specs[i].mesh.nPatches(), " patches; nu ",
                      num(fluid.mu / fluid.rho), "; turbulence ",
                      specs[i].turb.model);
            if (!specs[i].rheology.isNewtonian())
                log::info("block fluid ", edgeNames[i], ": n ",
                          num(specs[i].rheology.index()),
                          " K ", num(specs[i].rheology.consistencyOr(fluid.mu)),
                          " tau0 ", num(specs[i].rheology.yieldStress),
                          " shearFloor ", num(gnCtl.shearFloor),
                          " regularization ", num(gnCtl.regularization),
                          " relaxation ", num(gnCtl.relaxation));
        }

        bool allConverged = true;
        CoupledResult last;
        for (std::size_t s = 0; s < seeds.size(); ++s)
        {
            // one fresh block per seed: the block warm-starts from its own
            // previous solution, so sharing it between seeds would let seed 2
            // inherit seed 1's interface history
            std::vector<std::unique_ptr<DetailBlock>> owners;
            std::vector<ICoupledElement*> partners;
            for (std::size_t i = 0; i < specs.size(); ++i)
            {
                owners.push_back(std::make_unique<DetailBlock>(
                    specs[i].mesh, fluid.mu / fluid.rho, fluid.rho,
                    specs[i].schemes, specs[i].simple,
                    specs[i].probeMin, specs[i].probeMax, specs[i].turb,
                    "flat", specs[i].rheology, gnCtl,
                    specs[i].phases, specs[i].compressible));
                partners.push_back(owners.back().get());
            }

            // and one fresh network: the loop writes the converged coefficient
            // into the coupled pipe, so a reused instance would carry it as a
            // starting value and the seed would no longer be independent
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
                net.pipes[static_cast<std::size_t>(edge)].coupledResistance = seeds[s];
            }

            log::info("seed ", s, " starting resistance ", num(seeds[s]),
                      " interface ", (accelerates[s] ? "quasi-Newton"
                                                     : "plain Gauss-Seidel"));
            CoupledControls runCtl = ctl;
            runCtl.accelerate = accelerates[s];
            const CoupledResult res =
                solveCoupled1d3d(net, edges, partners, runCtl, netCtl);

            for (std::size_t i = 0; i < res.state.size(); ++i)
            {
                const EdgeExchange& st = res.state[i];
                if (multiEdge)
                    log::info("coupled seed ", s, " edge ", edgeNames[i],
                              ": iterations ", res.iterations,
                              " q ", num(st.q), " dp ", num(st.dpPa),
                              " R ", num(st.resistance),
                              " residual ", num(st.residual),
                              res.converged ? " (converged)" : " (NOT converged)",
                              " accelerate ", (accelerates[s] ? 1 : 0));
                else
                {
                    // legacy single-edge line, unchanged
                    log::info("coupled seed ", s, ": iterations ", res.iterations,
                              " q ", num(st.q), " dp ", num(st.dpPa),
                              " R ", num(st.resistance),
                              " residual ", num(res.residual),
                              res.converged ? " (converged)" : " (NOT converged)",
                              " accelerate ", (accelerates[s] ? 1 : 0));
                }
            }

            for (std::size_t i = 0; i < owners.size(); ++i)
            {
                const DetailBlock& block = *owners[i];
                const BlockMeasure& m = block.measure();
                const std::string tag = multiEdge ? edgeNames[i] : "";
                const std::string seedPrefix = multiEdge
                    ? "block seed " + std::to_string(s) + " edge " + tag
                    : "block seed " + std::to_string(s);
                log::info(seedPrefix, ": cells ", m.cells, " sweeps ", m.sweeps,
                          " blockConverged ", (m.converged ? 1 : 0),
                          " fluxIn ", num(m.fluxIn), " fluxOut ", num(m.fluxOut),
                          " planeFluxMin ", num(m.planeFluxMin),
                          " planeFluxMax ", num(m.planeFluxMax),
                          " inletArea ", num(m.inletArea),
                          " plugSpeed ", num(m.plugSpeed),
                          " uxMax ", num(m.uxMax), " uxCentre ", num(m.uxCentre),
                          " probeSections ", m.probeSections,
                          " probeSlope ", num(m.probeSlope),
                          " probeSpread ", num(m.probeSpread),
                          " dpPa ", num(m.dpPa),
                          " turb ", (m.turbModel.empty() ? std::string("laminar")
                                                         : m.turbModel),
                          " nutMax ", num(m.nutMax),
                          " wallYPlusMax ", num(m.wallYPlusMax),
                          " wallUTau ", num(m.wallUTau));
                if (const PhaseField* ph = block.phases())
                    log::info("block phases seed ", s,
                              multiEdge ? std::string(" edge ") + tag : "",
                              ": alphaMean ", num(ph->alphaMean()),
                              " rhoMean ", num(ph->rhoMean()),
                              " muMean ", num(ph->muMean()),
                              " bodyDrop ", num(ph->bodyDropPa()));
                // 5.79: the compressible partner reports its own state. Five
                // stations are enough for a checker to rebuild p^2(x) from the
                // dict alone without trusting the solver's own integral.
                if (block.compressible())
                {
                    const CompressibleMeasure& cm = block.compressibleMeasure();
                    log::info("block compressible seed ", s,
                              multiEdge ? std::string(" edge ") + tag : "",
                              ": pIn ", num(cm.pIn), " pOut ", num(cm.pOut),
                              " dp ", num(cm.dpPa),
                              " rhoIn ", num(cm.rhoIn),
                              " rhoOut ", num(cm.rhoOut),
                              " qIn ", num(cm.qIn), " qOut ", num(cm.qOut),
                              " massRate ", num(cm.massRate),
                              " uIn ", num(cm.uIn), " uOut ", num(cm.uOut),
                              " p2Slope ", num(cm.p2Slope),
                              " sections ", cm.sections,
                              " length ", num(cm.length),
                              " height ", num(cm.height),
                              " width ", num(cm.width),
                              " area ", num(cm.area));
                    const std::size_t n = cm.xStation.size();
                    if (n >= 2)
                    {
                        const std::size_t idx[5] = {0, n / 4, n / 2,
                                                   (3 * n) / 4, n - 1};
                        for (std::size_t j = 0; j < 5; ++j)
                        {
                            const std::size_t k = idx[j];
                            log::info("block compressible station ", k,
                                      ": x ", num(cm.xStation[k]),
                                      " p ", num(cm.pStation[k]),
                                      " rho ", num(cm.rhoStation[k]),
                                      " q ", num(cm.qStation[k]),
                                      " u ", num(cm.uStation[k]));
                        }
                    }
                }
                if (const GeneralizedNewtonianViscosity* gn =
                        block.generalizedNewtonian())
                    log::info("block fluid seed ", s,
                              multiEdge ? std::string(" edge ") + tag : "",
                              ": n ", num(gn->rheology().index()),
                              " K ", num(gn->rheology().consistencyOr(fluid.mu)),
                              " tau0 ", num(gn->rheology().yieldStress),
                              " gdotMax ", num(gn->maxShearRate()),
                              " nuMin ", num(gn->minViscosity()),
                              " nuMax ", num(gn->maxViscosity()),
                              " nuWall ", num(gn->meanWallViscosity()),
                              " gdotWall ", num(gn->meanWallShearRate()),
                              " tauWall ", num(gn->meanWallShearPa()));
            }

            for (const OuterIteration& it : res.history)
            {
                for (std::size_t i = 0; i < it.edges.size(); ++i)
                {
                    const EdgeExchange& ex = it.edges[i];
                    if (multiEdge)
                        log::info("interface seed ", s, " iter ", it.iteration,
                                  " edge ", edgeNames[i],
                                  " R ", num(ex.resistance), " q ", num(ex.q),
                                  " dpBlock ", num(ex.dpPa),
                                  " sweeps ", ex.sweeps,
                                  " residual ", num(ex.residual));
                    else
                        log::info("interface seed ", s, " iter ", it.iteration,
                                  " R ", num(ex.resistance), " q ", num(ex.q),
                                  " dpBlock ", num(ex.dpPa),
                                  " sweeps ", ex.sweeps,
                                  " residual ", num(ex.residual));
                }
            }

            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                const PipeFlowResult& edgeFlow =
                    res.network.pipeFlows[static_cast<std::size_t>(edges[i])];
                if (multiEdge)
                    log::info("network seed ", s, " edge ", edgeNames[i],
                              ": dpNetwork ", num(edgeFlow.dpTotal),
                              " qNetwork ", num(edgeFlow.flowRate),
                              " v ", num(edgeFlow.velocity),
                              " Re ", num(edgeFlow.reynolds),
                              " nodeResidual ", num(res.network.relativeResidual));
                else
                    log::info("network seed ", s, ": dpNetwork ",
                              num(edgeFlow.dpTotal),
                              " qNetwork ", num(edgeFlow.flowRate),
                              " v ", num(edgeFlow.velocity),
                              " Re ", num(edgeFlow.reynolds),
                              " nodeResidual ", num(res.network.relativeResidual));
            }

            if (s + 1 == seeds.size())
            {
                last = res;
                // Legacy: the single block's VTU path comes from output{ vtk };
                // multi-edge cases carry each block's path in its details{} spec.
                if (!multiEdge && caseDict.has("output"))
                {
                    const std::string vtk =
                        caseDict.getDict("output").getWordOr("vtk", "");
                    if (!vtk.empty())
                    {
                        VtkWriter w(owners[0]->mesh());
                        w.addCellVector("U", owners[0]->U().internalField());
                        w.addCellScalar("p", owners[0]->p().internalField());
                        writeBlockFields(w, *owners[0]);
                        w.write(vtk);
                    }
                }
                else
                {
                    for (std::size_t i = 0; i < specs.size(); ++i)
                    {
                        if (specs[i].vtk.empty()) continue;
                        VtkWriter w(owners[i]->mesh());
                        w.addCellVector("U", owners[i]->U().internalField());
                        w.addCellScalar("p", owners[i]->p().internalField());
                        writeBlockFields(w, *owners[i]);
                        w.write(specs[i].vtk);
                    }
                }
            }
            allConverged = allConverged && res.converged;
        }

        // the converged state, described by the 1D network alone: the checker
        // closes the pressure chain and the node splits from these numbers
        for (std::size_t i = 0; i < last.network.nodePressure.size(); ++i)
            log::info("node final p ", num(last.network.nodePressure[i]));
        for (const PipeFlowResult& f : last.network.pipeFlows)
        {
            log::info("pipe final ", f.name, " Q ", num(f.flowRate),
                      " v ", num(f.velocity), " Re ", num(f.reynolds),
                      " dp ", num(f.dpTotal));
        }

        // One solved line per coupled edge; the single-edge legacy line keeps
        // its original key set and order. Each edge reports its own exchange
        // residual (the aggregate `last.residual` is the max over edges).
        for (std::size_t i = 0; i < last.state.size(); ++i)
        {
            const EdgeExchange& ex = last.state[i];
            log::info("coupled solved: edge ", edgeNames[i],
                      " q ", num(ex.q),
                      " dp ", num(ex.dpPa),
                      " R ", num(ex.resistance),
                      " iterations ", last.iterations,
                      " residual ", num(ex.residual));
        }

        return allConverged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("coupled1d3d: ", e.what());
        return 1;
    }
}
