// chflow-coupled-1d3d-transient: a transient 1D pipe network with one 3D
// detailed region embedded in it. The region is a `coupledResistance` pipe on
// the 1D side and a SIMPLE block on the 3D side; inside every time step the two
// exchange rate-for-drop until the drop the network charges is the drop the
// block reports at the rate that step delivers.
//
// Every step is logged with its own state, so a regression can check the
// interface at each time level rather than only at the end: the ramp the two
// solvers produce together is the result, not a number the case was given.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include "cf/1d/Transient.h"
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
        std::cerr << "usage: chflow-coupled-1d3d-transient <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        // ---- the 1D side ----
        const TransientCase c = loadTransientFromDict(caseDict);
        if (c.fluid.rho <= 0 || c.fluid.mu <= 0)
            cfFatal("coupled1d3dTransient: fluid rho and mu must both be positive");
        const TransientTimeControls tc = loadTransientControls(caseDict);

        const Dictionary& cd = caseDict.getDict("couple");
        const std::string edgeName = cd.getWord("edge");
        TransientCoupledControls ctl;
        ctl.maxSubIterations = static_cast<label>(cd.getScalarOr("maxSubIterations", ctl.maxSubIterations));
        ctl.relax = cd.getScalarOr("relax", ctl.relax);
        ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");
        ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
        // the starting interface coefficients: each one runs the whole transient
        // series on its own, so a seed that lands anywhere else is not a
        // coupling but a trajectory the initial value still owns
        std::vector<double> seeds;
        if (cd.has("seeds")) seeds = cd.getScalarList("seeds");
        if (seeds.empty())
            cfFatal("coupled1d3dTransient: couple{} needs at least one starting"
                    " resistance in seeds ( )");

        // ---- the 3D side ----
        const Dictionary& dd = caseDict.getDict("detail");
        polyMesh mesh = MeshLoader::fromCase(dd, argv[1]);
        const Schemes schemes = loadSchemes(dd);
        const SimpleControls blockCtl = loadSimple(dd);
        scalar probeMin = -std::numeric_limits<scalar>::infinity();
        scalar probeMax = probeMin;
        if (dd.has("probe"))
        {
            const Dictionary& pd = dd.getDict("probe");
            probeMin = pd.getScalar("min");
            probeMax = pd.getScalar("max");
            if (!(probeMin < probeMax))
                cfFatal("coupled1d3dTransient: detail probe min must be below max");
        }
        else
        {
            cfFatal("coupled1d3dTransient: detail{} needs probe{ min max; }, the"
                    " cross-section window the developed gradient is measured"
                    " over");
        }

        log::info("transient network+block: block has ", mesh.nCells(), " cells and ",
                  mesh.nPatches(), " patches; nu ", num(c.fluid.mu / c.fluid.rho),
                  "; dt ", num(tc.dt), " endTime ", num(tc.endTime),
                  "; interface ", (ctl.accelerate ? "quasi-Newton" : "plain Gauss-Seidel"));

        for (std::size_t s = 0; s < seeds.size(); ++s)
        {
            // one fresh block per seed: the block warm-starts from its own
            // previous solution, so sharing it between seeds would let seed 2
            // inherit seed 1's interface history
            DetailBlock block(mesh, c.fluid.mu / c.fluid.rho, c.fluid.rho, schemes,
                               blockCtl, probeMin, probeMax);

            // and one fresh network: the agreed coefficient carries from step to
            // step, so a reused instance would start the series where the
            // previous seed left it
            TransientNetwork net(c);
            label edge = -1;
            for (std::size_t i = 0; i < net.pipeNames().size(); ++i)
            {
                if (net.pipeNames()[i] == edgeName) edge = static_cast<label>(i);
            }
            if (edge < 0)
                cfFatal("coupled1d3dTransient: no pipe named '" + edgeName
                        + "' in the network");
            net.setCoupledResistance(edge, seeds[s]);

            log::info("transient seed ", s, " starting resistance ", num(seeds[s]));
            const TransientCoupledResult res =
                solveTransientCoupled1d3d(net, {edge}, {&block}, ctl, tc.dt, tc.endTime);

            for (std::size_t i = 0; i < res.steps.size(); ++i)
            {
                const TransientCoupledStep& st = res.steps[i];
                const TransientInterfacePass& p = st.history.back();
                const EdgeExchange& ex = st.edges[0];
                log::info("transient step ", i + 1, ": t ", num(st.t),
                          " q ", num(ex.q), " dpNetwork ",
                          num(ex.resistance * ex.q),
                          " dpBlock ", num(ex.dpPa), " R ", num(ex.resistance),
                          " residual ", num(ex.residual),
                          " passes ", st.passes, " subPasses ", st.history.size(),
                          " sweeps ", ex.sweeps, " newton ", num(p.newtonResidual),
                          " blockConverged ", (ex.blockConverged ? 1 : 0));

                std::ostringstream nodes;
                nodes << std::setprecision(17);
                for (std::size_t i2 = 0; i2 < res.nodeNames.size(); ++i2)
                    nodes << " " << res.nodeNames[i2] << " " << st.pressure[i2];
                log::info("nodes step ", i + 1, ":", nodes.str());

                std::ostringstream pipes;
                pipes << std::setprecision(17);
                for (std::size_t i2 = 0; i2 < res.pipeNames.size(); ++i2)
                    pipes << " " << res.pipeNames[i2] << " " << st.flowRate[i2];
                log::info("pipes step ", i + 1, ":", pipes.str());
            }

            // the block's own field describes the last rate it was asked for,
            // which is the rate the series ended at
            const BlockMeasure& m = block.measure();
            log::info("transient block ", s, ": cells ", m.cells,
                      " sweeps ", m.sweeps,
                      " blockConverged ", (m.converged ? 1 : 0),
                      " fluxIn ", num(m.fluxIn), " fluxOut ", num(m.fluxOut),
                      " planeFluxMin ", num(m.planeFluxMin),
                      " planeFluxMax ", num(m.planeFluxMax),
                      " inletArea ", num(m.inletArea), " plugSpeed ", num(m.plugSpeed),
                      " uxMax ", num(m.uxMax), " uxCentre ", num(m.uxCentre),
                      " probeSections ", m.probeSections,
                      " probeSlope ", num(m.probeSlope),
                      " probeSpread ", num(m.probeSpread),
                      " dpPa ", num(m.dpPa));

            log::info("transient coupled seed ", s, ": edge ", edgeName,
                      " steps ", res.steps.size(), " q ", num(res.state[0].q),
                      " dp ", num(res.state[0].dpPa),
                      " R ", num(res.state[0].resistance),
                      " residual ", num(res.residual),
                      " t ", num(net.time()));

            if (s + 1 == seeds.size() && caseDict.has("output"))
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

        // a step the interface cannot agree on is fatal inside the loop, so
        // reaching here means every time level was solved consistently
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("coupled1d3dTransient: ", e.what());
        return 1;
    }
}
