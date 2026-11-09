// chflow-coupled-1d3d-transient-multi: a transient 1D pipe network with several
// 3D detailed regions embedded in it. One interface pass solves the network once
// and asks every block for the drop at the rate that pass delivers through its
// own edge, so all edges are refreshed against one system state and the pass
// stops on the worst of them -- inside every time step.
//
// The steady multi-edge app (chflow-coupled-1d3d-multi) answers the same network
// without time. Sharing the topology and the two block descriptions is what
// makes the two comparable: both reduce their partners to a per-edge resistance,
// and each checker asks for its own duct's planar Poiseuille value, so the two
// published deficits are the same physical statement about the same two blocks.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

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
        cfFatal("coupled1d3dTransient: details." + edgeName + "{} needs probe"
                " { min max; }, the cross-section window the developed gradient"
                " is measured over");
    const Dictionary& pd = dd.getDict("probe");
    const scalar probeMin = pd.getScalar("min");
    const scalar probeMax = pd.getScalar("max");
    if (!(probeMin < probeMax))
        cfFatal("coupled1d3dTransient: details." + edgeName + " probe min must"
                " be below max");
    return BlockSpec{MeshLoader::fromCase(dd, argv0), loadSchemes(dd), loadSimple(dd),
                     probeMin, probeMax, dd.getWordOr("vtk", ""),
                     TurbulenceConfig::fromCase(dd)};
}

/// Where one coupled edge's inertance comes from, in the order the interface
/// takes its decision: the drawing's bore, the partner's whole region, or the
/// slab of that region whose pressure difference the partner actually reports.
/// The last is only a different number from the second for a block that steps
/// in time -- and only a different *claim* once one does, because a
/// quasi-steady partner can measure both and cannot say which its drop is worth.
int inertiaSourceIndex(const std::string& w)
{
    if (w == "network") return 0;
    if (w == "partner") return 1;
    if (w == "slab") return 2;
    cfFatal("coupled1d3dTransient: inertia source '" + w
            + "' is none of network, partner, slab");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-coupled-1d3d-transient-multi <case.dict>\n";
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
        const double nu = c.fluid.mu / c.fluid.rho;

        const Dictionary& cd = caseDict.getDict("couple");
        const std::vector<std::string> edgeNames = cd.getWordList("edges");
        if (edgeNames.empty())
            cfFatal("coupled1d3dTransient: couple{} needs the coupled pipe names"
                    " in edges ( )");

        TransientCoupledControls ctl;
        ctl.maxSubIterations = static_cast<label>(
            cd.getScalarOr("maxSubIterations", ctl.maxSubIterations));
        ctl.relax = cd.getScalarOr("relax", ctl.relax);
        ctl.tolerance = cd.getScalarOr("tolerance", ctl.tolerance);
        ctl.accelerate = cd.has("accelerate") && cd.getBool("accelerate");

        // Each named set is one whole-series run and quotes a starting
        // coefficient for every coupled edge at once: an edge that starts far off
        // must not be able to move its neighbour's trajectory.
        if (!cd.has("starts"))
            cfFatal("coupled1d3dTransient: couple{} needs starts{ name"
                    " { <edge> <R>; ... } }, one whole-network starting set per run");
        const Dictionary& startsDict = cd.getDict("starts");
        std::vector<std::string> startNames;
        std::vector<std::vector<double>> startRs;
        std::vector<std::vector<std::string>> startShapes;
        for (const std::string& sname : startsDict.keys())
        {
            const Dictionary& sd = startsDict.getDict(sname);
            std::vector<double> rs;
            for (const std::string& edgeName : edgeNames)
            {
                if (!sd.has(edgeName))
                    cfFatal("coupled1d3dTransient: starts." + sname + " gives no"
                            " starting resistance for coupled edge '" + edgeName + "'");
                rs.push_back(sd.getScalar(edgeName));
            }
            startNames.push_back(sname);
            startRs.push_back(std::move(rs));

            // The same arm also says what shape its rate arrives in, exactly as
            // on the steady side. An edge without a `profile` keeps the flat plug
            // the interface has always written, so a transient case that names
            // none is the case it was -- and inside a time step the shape is not
            // a one-off: every pass of every level re-writes the inlet faces.
            std::vector<std::string> shapes(edgeNames.size(), "flat");
            if (sd.has("profile"))
            {
                const Dictionary::Kind kind = sd.find("profile").kind;
                if (kind == Dictionary::Kind::word)
                {
                    const std::string one = sd.getWord("profile");
                    std::fill(shapes.begin(), shapes.end(), one);
                }
                else if (kind == Dictionary::Kind::list)
                {
                    shapes = sd.getWordList("profile");
                    if (shapes.size() != edgeNames.size())
                        cfFatal("coupled1d3dTransient: starts." + sname
                                + " profile names " + std::to_string(shapes.size())
                                + " shapes for "
                                + std::to_string(edgeNames.size()) + " coupled edges");
                }
                else
                    cfFatal("coupled1d3dTransient: starts." + sname + " profile must"
                            " be one shape word or a list of them, one per coupled"
                            " edge");
            }
            startShapes.push_back(std::move(shapes));
        }
        if (startNames.empty())
            cfFatal("coupled1d3dTransient: couple.starts{} needs at least one set");

        // Which interface operator, and where each region's inertance comes from,
        // are per-run choices: one case can then send the same network and the
        // same blocks through Gauss-Seidel and through the quasi-Newton, or
        // through the bore's inertance and the block's own measurement, and the
        // arms differ in that one thing and nothing else.
        std::vector<bool> accelerates(startNames.size(), ctl.accelerate);
        if (cd.has("accelerates"))
        {
            const std::vector<double> a = cd.getScalarList("accelerates");
            if (a.size() == 1)
                std::fill(accelerates.begin(), accelerates.end(), a[0] != 0.0);
            else if (a.size() != accelerates.size())
                cfFatal("coupled1d3dTransient: couple.accelerates needs one entry"
                        " per starting set, or exactly one for all of them (got "
                        + std::to_string(a.size()) + " for "
                        + std::to_string(accelerates.size()) + " sets)");
            else
                for (std::size_t i = 0; i < a.size(); ++i) accelerates[i] = a[i] != 0.0;
        }

        std::vector<std::string> inertias(startNames.size(), "network");
        if (cd.has("inertias"))
        {
            const std::vector<std::string> w = cd.getWordList("inertias");
            if (w.size() == 1 && startNames.size() > 1)
                inertias.assign(startNames.size(), w[0]);
            else if (w.size() != inertias.size())
                cfFatal("coupled1d3dTransient: couple.inertias needs one entry"
                        " per starting set, or exactly one for all of them (got "
                        + std::to_string(w.size()) + " for "
                        + std::to_string(inertias.size()) + " sets)");
            else
                inertias = w;
        }
        for (const std::string& w : inertias) inertiaSourceIndex(w);

        // Whether the partner carries its own acceleration, arm by arm. The
        // blocks are built fresh for every arm, so one run can put the same
        // mesh, the same seed and the same tolerance under a quasi-steady
        // partner and under a time-accurate one: the only difference between
        // those two arms is whether the region's `du/dt` was in the equation
        // that produced its drop. Omitting the entry leaves every arm at what
        // its own `details.<edge>.schemes{ ddt ... }` says.
        std::vector<std::string> ddts(startNames.size(), "");
        if (cd.has("ddts"))
        {
            const std::vector<std::string> w = cd.getWordList("ddts");
            if (w.size() == 1 && startNames.size() > 1)
                ddts.assign(startNames.size(), w[0]);
            else if (w.size() != ddts.size())
                cfFatal("coupled1d3dTransient: couple.ddts needs one entry per"
                        " starting set, or exactly one for all of them (got "
                        + std::to_string(w.size()) + " for "
                        + std::to_string(ddts.size()) + " sets)");
            else
                ddts = w;
        }
        for (const std::string& w : ddts)
            if (!w.empty() && w != "steady" && w != "Euler")
                cfFatal("coupled1d3dTransient: ddts entry '" + w
                        + "' is neither steady nor Euler");

        // The same decision made per edge, so one run can trust the block on the
        // edge whose region it can measure and the drawing on the edge whose
        // bore was typed in. Only the runs that need mixing declare anything;
        // the rest keep answering to the per-arm switch above, unchanged.
        std::vector<std::vector<std::string>> armEdgeInertias(startNames.size());
        bool anyPerEdge = false;
        if (cd.has("inertiaFrom"))
        {
            const Dictionary& from = cd.getDict("inertiaFrom");
            for (std::size_t s = 0; s < startNames.size(); ++s)
            {
                if (!from.has(startNames[s])) continue;
                const Dictionary& sd = from.getDict(startNames[s]);
                std::vector<std::string> w;
                for (const std::string& edgeName : edgeNames)
                {
                    if (!sd.has(edgeName))
                        cfFatal("coupled1d3dTransient: inertiaFrom." + startNames[s]
                                + " gives no source for coupled edge '" + edgeName + "'");
                    const std::string e = sd.getWord(edgeName);
                    inertiaSourceIndex(e);
                    w.push_back(e);
                }
                armEdgeInertias[s] = std::move(w);
                anyPerEdge = true;
            }
        }

        // ---- the 3D side: one block per coupled edge ----
        const Dictionary& allDetails = caseDict.getDict("details");
        std::vector<BlockSpec> specs;
        for (const std::string& edgeName : edgeNames)
        {
            if (!allDetails.has(edgeName))
                cfFatal("coupled1d3dTransient: edge '" + edgeName + "' has no"
                        " details." + edgeName + "{} block");
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
        log::info("transient network+blocks: ", edgeNames.size(), " coupled edges,",
                  " ", startRs.size(), " starting sets; nu ", num(nu),
                  "; dt ", num(tc.dt), " endTime ", num(tc.endTime),
                  "; interface ", (ctl.accelerate ? "quasi-Newton" : "plain Gauss-Seidel"),
                  "; region inertia ", inertias[0],
                  (inertias.size() > 1 ? " (arms)" : ""),
                  (anyPerEdge ? " (per-edge on some arms)" : ""));

        for (std::size_t s = 0; s < startNames.size(); ++s)
        {
            // one fresh set of blocks per starting set: each warm-starts from its
            // own previous solution, so sharing them between sets would let the
            // second run inherit the first one's interface history
            std::vector<std::unique_ptr<DetailBlock>> owners;
            std::vector<ICoupledElement*> partners;
            std::vector<std::string> blockDdt;
            for (std::size_t i = 0; i < specs.size(); ++i)
            {
                Schemes sch = specs[i].schemes;
                if (!ddts[s].empty()) sch.ddt = ddts[s];
                blockDdt.push_back(sch.ddt);
                owners.push_back(std::make_unique<DetailBlock>(
                    specs[i].mesh, nu, c.fluid.rho, sch, specs[i].simple,
                    specs[i].probeMin, specs[i].probeMax, specs[i].turb,
                    startShapes[s][i]));
                partners.push_back(owners.back().get());
            }

            // and one fresh network: the loop writes the agreed coefficients and
            // the agreed inertances into the coupled pipes, and the step-to-step
            // ramp lives in the pipes' own qPrev, so a reused instance would not
            // start where the case says it starts
            TransientNetwork net(c);
            std::vector<label> edges;
            for (const std::string& edgeName : edgeNames)
            {
                label edge = -1;
                for (std::size_t i = 0; i < net.pipeNames().size(); ++i)
                {
                    if (net.pipeNames()[i] == edgeName) edge = static_cast<label>(i);
                }
                if (edge < 0)
                    cfFatal("coupled1d3dTransient: no pipe named '" + edgeName
                            + "' in the network");
                edges.push_back(edge);
            }
            for (std::size_t i = 0; i < edges.size(); ++i)
                net.setCoupledResistance(edges[i], startRs[s][i]);

            TransientCoupledControls run = ctl;
            run.accelerate = accelerates[s];
            run.inertiaFromPartner = (inertias[s] == "partner");
            const std::vector<std::string>& edgeWords = armEdgeInertias[s];
            if (!edgeWords.empty())
            {
                std::vector<int> src;
                for (const std::string& w : edgeWords)
                    src.push_back(inertiaSourceIndex(w));
                run.inertiaFromPartnerPerEdge = std::move(src);
            }
            else if (inertias[s] == "slab")
            {
                // one answer for the whole interface still travels as the
                // per-edge decision, because the switch underneath it is a
                // yes-or-no and the slab is a third answer
                run.inertiaFromPartnerPerEdge = std::vector<int>(edges.size(), 2);
            }

            log::info("transient start ", startNames[s], ": inertia ", inertias[s],
                      " accelerate ", (run.accelerate ? "yes" : "no"));
            // Which source each edge of this arm actually ends up using - the
            // arm's own switch when the case names no table for it, the table's
            // entry when it does. `transient arm` rather than `transient start`
            // because this line is per edge, and the runs' readers split the log
            // into arms on the `transient start <name>: inertia ...` header.
            for (std::size_t i = 0; i < edges.size(); ++i)
                log::info("transient arm ", startNames[s], " edge ", edgeNames[i],
                          " inertia ",
                          (edgeWords.empty() ? inertias[s] : edgeWords[i]),
                          " index ", edges[i]);
            for (std::size_t i = 0; i < edges.size(); ++i)
                log::info("transient start ", startNames[s], " edge ", edgeNames[i],
                          " index ", edges[i], " resistance ", num(startRs[s][i]));
            // and what shape this arm's rate arrives in on each edge. Per edge,
            // because the interface asks each region separately, and per arm,
            // because one case can hand the same mesh a plug in one of its
            // starting sets and a Poiseuille profile in another.
            for (std::size_t i = 0; i < edges.size(); ++i)
                log::info("transient arm ", startNames[s], " edge ", edgeNames[i],
                          " profile ", startShapes[s][i], " index ", edges[i]);
            // What this arm's partner is going to be asked to carry: the case's
            // own ddt word, and the two inertances its mesh measures -- the
            // region's and that of only the slab whose pressure difference it
            // reports. `slab` installs the second, every arm is graded on it.
            for (std::size_t i = 0; i < edges.size(); ++i)
                log::info("transient partner ", startNames[s], " edge ", edgeNames[i],
                          " ddt ", blockDdt[i],
                          " regionInertance ", num(owners[i]->regionInertance()),
                          " probeInertance ", num(owners[i]->probeInertance()),
                          " index ", edges[i]);

            const TransientCoupledResult res =
                solveTransientCoupled1d3d(net, edges, partners, run, tc.dt, tc.endTime);

            for (std::size_t n = 0; n < res.steps.size(); ++n)
            {
                const TransientCoupledStep& st = res.steps[n];

                // the schedule the step's passes ran through, edge by edge: the
                // stop is taken on the worst of them, so which edge a pass was
                // spent for is a fact about this loop rather than an assumption
                for (const TransientInterfacePass& pass : st.history)
                {
                    label worst = 0;
                    for (std::size_t i = 1; i < pass.edges.size(); ++i)
                    {
                        if (pass.edges[i].residual > pass.edges[worst].residual)
                            worst = static_cast<label>(i);
                    }
                    std::string rel, coef, rate;
                    for (const EdgeExchange& ex : pass.edges)
                    {
                        rel += " " + num(ex.residual);
                        coef += " " + num(ex.resistance);
                        rate += " " + num(ex.q);
                    }
                    log::info("transient step ", n + 1, " pass ", pass.pass,
                              " worst ", worst, " residuals", rel, " coeff", coef,
                              " rates", rate);
                }

                log::info("transient step ", n + 1, ": t ", num(st.t),
                          " passes ", st.passes, " subPasses ", st.history.size(),
                          " newton ", num(st.history.back().newtonResidual));
                for (std::size_t i = 0; i < st.edges.size(); ++i)
                {
                    const EdgeExchange& ex = st.edges[i];
                    log::info("transient step ", n + 1, " edge ", i, " name ", ex.name,
                              " q ", num(ex.q), " dpNetwork ", num(ex.resistance * ex.q),
                              " dpBlock ", num(ex.dpPa), " R ", num(ex.resistance),
                              " residual ", num(ex.residual), " sweeps ", ex.sweeps,
                              " M ", num(ex.inertance),
                              " blockConverged ", (ex.blockConverged ? 1 : 0),
                              // an unsteady partner splits its drop: the
                              // resistance is charged with the friction half and
                              // the row's inertance term is graded against the
                              // other half, which is what `rowDrop`/`rowError`
                              // report. A quasi-steady one leaves them at its
                              // whole drop and an error of 0: there is nothing to
                              // grade it on.
                              " dpFriction ", num(ex.dpFrictionPa),
                              " dpInertia ", num(ex.dpInertialPa),
                              " qPrev ", num(ex.qPrev), " dq ", num(ex.dq),
                              " rowDrop ", num(ex.rowDropPa),
                              " rowError ", num(ex.rowError),
                              // the shape travels with the rate at this level,
                              // not with the run: what was asked for, what the
                              // inlet plane of the solved field reads, and what
                              // its outlet reads. `fluxImposed` is the section
                              // integral of the condition the inlet really got,
                              // so it equals this level's `q` however the rate
                              // was spread. A plug that has not finished
                              // developing shows up here as betaOut above
                              // betaIn, level after level.
                              " inletShape ", ex.inletShape,
                              " fluxImposed ", num(ex.fluxImposed),
                              " betaImposed ", num(ex.betaImposed),
                              " betaIn ", num(ex.betaIn), " betaOut ", num(ex.betaOut),
                              " inletArea ", num(ex.inletArea),
                              // How much of this level's drop its own solve handed
                              // back: the highest station pressure of that same
                              // field above its first station. A constant-area duct
                              // can only pull pressure down, so a positive rise is
                              // the inlet's excess momentum flux relaxing toward the
                              // region's own shape -- nil on a plug, and here on the
                              // level that asked for it rather than once per arm.
                              " pRisePa ", num(ex.pRisePa), " pRiseX ", num(ex.pRiseX),
                              " unsteady ", (ex.unsteady ? 1 : 0));
                }

                std::ostringstream nodes;
                nodes << std::setprecision(17);
                for (std::size_t i = 0; i < res.nodeNames.size(); ++i)
                    nodes << " " << res.nodeNames[i] << " " << st.pressure[i];
                log::info("nodes step ", n + 1, ":", nodes.str());

                std::ostringstream pipes;
                pipes << std::setprecision(17);
                for (std::size_t i = 0; i < res.pipeNames.size(); ++i)
                    pipes << " " << res.pipeNames[i] << " " << st.flowRate[i];
                log::info("pipes step ", n + 1, ":", pipes.str());
            }

            // each block's field describes the last rate it was asked for, which
            // is the rate the series ended with on that edge
            for (std::size_t i = 0; i < owners.size(); ++i)
            {
                const BlockMeasure& m = owners[i]->measure();
                log::info("transient block ", i, " name ", edgeNames[i],
                          ": cells ", m.cells, " sweeps ", m.sweeps,
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
                          " inertance ", num(m.inertance),
                          " turb ", (m.turbModel.empty() ? std::string("laminar")
                                                         : m.turbModel),
                          " nutMax ", num(m.nutMax),
                          " wallYPlusMax ", num(m.wallYPlusMax),
                          " wallUTau ", num(m.wallUTau),
                          // what the last solve made of its own acceleration:
                          // the force the slab's history asks for, what that is
                          // worth per unit rate-of-change, the inertance of the
                          // slab it was measured over, and the profile's own
                          // momentum-flux correction
                          " unsteady ", (m.unsteady ? 1 : 0),
                          " dpInertialPa ", num(m.dpInertialPa),
                          " dpFrictionPa ", num(m.dpFrictionPa),
                          " inertanceDynamic ", num(m.inertanceDynamic),
                          " inertanceProbe ", num(m.inertanceProbe),
                          " betaProfile ", num(m.betaProfile),
                          // the same shape ledger the steady multi-edge app
                          // prints, now for the last level the run accepted: what
                          // was handed over, as a rate and as a distribution, and
                          // what the region's own field ended up carrying
                          " inletShape ", m.inletShape, " fluxImposed ",
                          num(m.fluxImposed), " betaImposed ", num(m.betaImposed),
                          " momentumFluxPa ", num(m.momentumFluxImposedPa),
                          " betaIn ", num(m.betaIn), " betaOut ", num(m.betaOut),
                          " dtLevel ", num(m.dtLevel));
                std::ostringstream betas;
                betas << std::setprecision(17);
                for (const std::pair<scalar, scalar>& sb : m.betaStations)
                    betas << " " << sb.first << " " << sb.second;
                if (!m.betaStations.empty())
                    log::info("transient betaProfile ", startNames[s], " edge ",
                              edgeNames[i], ":", betas.str());
                // the same stations, now carrying pressure: the entry region's
                // bill for the shape it was handed, as a curve instead of the
                // single `pRisePa` number above. Reported for the last level the
                // run accepted, exactly as the beta sequence is.
                std::ostringstream pStations;
                pStations << std::setprecision(17);
                for (const std::pair<scalar, scalar>& sp : m.pressureStations)
                    pStations << " " << sp.first << " " << sp.second;
                if (!m.pressureStations.empty())
                    log::info("transient pressureProfile ", startNames[s], " edge ",
                              edgeNames[i], ":", pStations.str());
            }

            for (std::size_t i = 0; i < res.state.size(); ++i)
            {
                const EdgeExchange& ex = res.state[i];
                log::info("transient coupled start ", startNames[s], " edge ",
                          edgeNames[i], " name ", ex.name, " q ", num(ex.q),
                          " dp ", num(ex.dpPa), " R ", num(ex.resistance),
                          " residual ", num(ex.residual), " M ", num(ex.inertance));
            }
            log::info("transient coupled start ", startNames[s], ": steps ",
                      res.steps.size(), " worst residual ", num(res.residual),
                      " t ", num(net.time()));

            // the fields are written once per case, for the last run: every
            // starting set solves the same blocks and would overwrite the files
            if (s + 1 == startNames.size())
            {
                for (std::size_t i = 0; i < specs.size(); ++i)
                {
                    if (specs[i].vtk.empty()) continue;
                    VtkWriter w(owners[i]->mesh());
                    w.addCellVector("U", owners[i]->U().internalField());
                    w.addCellScalar("p", owners[i]->p().internalField());
                    w.write(specs[i].vtk);
                }
            }
        }

        // a step any edge cannot agree on is fatal inside the loop, so reaching
        // here means every time level was solved consistently on every edge
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("coupled1d3dTransient: ", e.what());
        return 1;
    }
}
