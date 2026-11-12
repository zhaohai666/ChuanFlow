// chflow-coupled-1d3d: a steady 1D pipe network with one 3D detailed region
// embedded in it. The detailed region is a `coupledResistance` pipe on the 1D
// side and a SIMPLE block on the 3D side; the two exchange rate-for-drop until
// the network charges the drop the block reports.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

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
        const std::string edgeName = cd.getWord("edge");
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

        // ---- the 3D side ----
        const Dictionary& dd = caseDict.getDict("detail");
        polyMesh mesh = MeshLoader::fromCase(dd, argv[1]);
        const Schemes schemes = loadSchemes(dd);
        const SimpleControls blockCtl = loadSimple(dd);
        // A detail{} without a turbulence{} block is laminar, exactly as the
        // standalone 3D apps read their own case.
        const TurbulenceConfig blockTurb = TurbulenceConfig::fromCase(dd);
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
        else
        {
            cfFatal("coupled1d3d: detail{} needs probe{ min max; }, the"
                    " cross-section window the developed gradient is measured"
                    " over");
        }

        const Rheology blockRheology = readPartnerRheology(caseDict, edgeName);
        const GeneralizedNewtonianControls gnCtl =
            readGeneralizedNewtonianControls(caseDict);
        // 5.78: a two-phase partner block. Absent in every older case.
        const PhaseConfig blockPhases = readPhases(dd);
        log::info("network+block: block has ", mesh.nCells(), " cells and ",
                  mesh.nPatches(), " patches; nu ", num(fluid.mu / fluid.rho),
                  "; turbulence ", blockTurb.model);
        if (!blockRheology.isNewtonian())
            log::info("block fluid: n ", num(blockRheology.index()),
                      " K ", num(blockRheology.consistencyOr(fluid.mu)),
                      " tau0 ", num(blockRheology.yieldStress),
                      " shearFloor ", num(gnCtl.shearFloor),
                      " regularization ", num(gnCtl.regularization),
                      " relaxation ", num(gnCtl.relaxation));

        bool allConverged = true;
        CoupledResult last;
        for (std::size_t s = 0; s < seeds.size(); ++s)
        {
            // one fresh block per seed: the block warm-starts from its own
            // previous solution, so sharing it between seeds would let seed 2
            // inherit seed 1's interface history
            DetailBlock block(mesh, fluid.mu / fluid.rho, fluid.rho, schemes,
                               blockCtl, probeMin, probeMax, blockTurb, "flat",
                               blockRheology, gnCtl, blockPhases);

            // and one fresh network: the loop writes the converged coefficient
            // into the coupled pipe, so a reused instance would carry it as a
            // starting value and the seed would no longer be independent
            OneDNetwork net = loadNetworkFromDict(caseDict);
            label edge = -1;
            for (std::size_t i = 0; i < net.pipes.size(); ++i)
            {
                if (net.pipes[i].name == edgeName) edge = static_cast<label>(i);
            }
            if (edge < 0)
                cfFatal("coupled1d3d: no pipe named '" + edgeName
                        + "' in the network");
            net.pipes[static_cast<std::size_t>(edge)].coupledResistance = seeds[s];

            log::info("seed ", s, " starting resistance ", num(seeds[s]),
                      " interface ", (accelerates[s] ? "quasi-Newton"
                                                     : "plain Gauss-Seidel"));
            CoupledControls runCtl = ctl;
            runCtl.accelerate = accelerates[s];
            const CoupledResult res =
                solveCoupled1d3d(net, {edge}, {&block}, runCtl, netCtl);
            const EdgeExchange& st = res.state[0];

            log::info("coupled seed ", s, ": iterations ", res.iterations,
                      " q ", num(st.q), " dp ", num(st.dpPa),
                      " R ", num(st.resistance), " residual ", num(res.residual),
                      res.converged ? " (converged)" : " (NOT converged)",
                      " accelerate ", (accelerates[s] ? 1 : 0));

            const BlockMeasure& m = block.measure();
            log::info("block seed ", s, ": cells ", m.cells, " sweeps ", m.sweeps,
                      " blockConverged ", (m.converged ? 1 : 0),
                      " fluxIn ", num(m.fluxIn), " fluxOut ", num(m.fluxOut),
                      " planeFluxMin ", num(m.planeFluxMin),
                      " planeFluxMax ", num(m.planeFluxMax),
                      " inletArea ", num(m.inletArea), " plugSpeed ", num(m.plugSpeed),
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
                          ": alphaMean ", num(ph->alphaMean()),
                          " rhoMean ", num(ph->rhoMean()),
                          " muMean ", num(ph->muMean()),
                          " bodyDrop ", num(ph->bodyDropPa()));
            if (const GeneralizedNewtonianViscosity* gn =
                    block.generalizedNewtonian())
                log::info("block fluid seed ", s, ": n ",
                          num(gn->rheology().index()),
                          " K ", num(gn->rheology().consistencyOr(fluid.mu)),
                          " tau0 ", num(gn->rheology().yieldStress),
                          " gdotMax ", num(gn->maxShearRate()),
                          " nuMin ", num(gn->minViscosity()),
                          " nuMax ", num(gn->maxViscosity()),
                          " nuWall ", num(gn->meanWallViscosity()),
                          " gdotWall ", num(gn->meanWallShearRate()),
                          " tauWall ", num(gn->meanWallShearPa()));

            for (const OuterIteration& it : res.history)
            {
                const EdgeExchange& ex = it.edges[0];
                log::info("interface seed ", s, " iter ", it.iteration,
                          " R ", num(ex.resistance), " q ", num(ex.q),
                          " dpBlock ", num(ex.dpPa),
                          " sweeps ", ex.sweeps,
                          " residual ", num(ex.residual));
            }

            const PipeFlowResult& edgeFlow =
                res.network.pipeFlows[static_cast<std::size_t>(edge)];
            log::info("network seed ", s, ": dpNetwork ", num(edgeFlow.dpTotal),
                      " qNetwork ", num(edgeFlow.flowRate), " v ", num(edgeFlow.velocity),
                      " Re ", num(edgeFlow.reynolds),
                      " nodeResidual ", num(res.network.relativeResidual));

            if (s + 1 == seeds.size())
            {
                last = res;
                if (caseDict.has("output"))
                {
                    const std::string vtk =
                        caseDict.getDict("output").getWordOr("vtk", "");
                    if (!vtk.empty())
                    {
                        VtkWriter w(block.mesh());
                        w.addCellVector("U", block.U().internalField());
                        w.addCellScalar("p", block.p().internalField());
                        w.write(vtk);
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

        const EdgeExchange& ex = last.state[0];
        log::info("coupled solved: edge ", edgeName, " q ", num(ex.q),
                  " dp ", num(ex.dpPa), " R ", num(ex.resistance),
                  " iterations ", last.iterations,
                  " residual ", num(last.residual));

        return allConverged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("coupled1d3d: ", e.what());
        return 1;
    }
}
