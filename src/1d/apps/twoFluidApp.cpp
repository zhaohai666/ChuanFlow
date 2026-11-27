#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/TwoFluid.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

using namespace cf;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

// G1: the per-cell composition columns - " comp <beta>" (PR flash vapour
// fraction diagnostic) plus the transported mole fractions " z0 ... zN-1"
// (empty string when no composition{} feed was declared)
std::string compositionReport(const TwoFluidPipe& s, int i)
{
    if (!s.compositionOn()) return std::string();
    std::string out = " comp " + num(s.cellVaporFraction()[i]);
    const int nc = s.compositionCount();
    for (int k = 0; k < nc; ++k)
        out += " z" + std::to_string(k) + " " + num(s.cellComposition(k)[i]);
    return out;
}

// G1: the declaration line naming the feed components and their z0
std::string compositionDecl(const TwoFluidCase& c)
{
    std::string line = "twoFluid composition: "
                     + std::to_string(c.compositionNames.size())
                     + " components";
    for (std::size_t k = 0; k < c.compositionNames.size(); ++k)
        line += " " + c.compositionNames[k] + " " + num(c.compositionZ0[k]);
    return line;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-1d-twofluid <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        const double dt = caseDict.getDict("time").getScalar("dt");
        const double endTime = caseDict.getDict("time").getScalar("endTime");
        const int reportEvery = static_cast<int>(
            caseDict.getDict("time").getScalarOr("reportEvery", 1.0));
        const int report = reportEvery > 0 ? reportEvery : 1;
        const int nSteps = static_cast<int>(std::round(endTime / dt));

        if (caseDict.has("network"))
        {
            const TwoFluidNetworkCase nc = loadTwoFluidNetworkCase(caseDict);
            int sections = 0;
            for (const TwoFluidGeo& g : nc.pipes) sections += g.sections;
            log::info("twoFluid network: ", nc.pipes.size(), " pipes ",
                      sections, " sections, ", nc.junctions.size(),
                      " junctions, dt=", num(dt), " endTime=", num(endTime));
            // 5.75 (A5): same self-naming as the serial branch
            if (nc.perPipe.fluid.flowIndex > 0)
                log::info("twoFluid liquid rheology: flowIndex ",
                          num(nc.perPipe.fluid.flowIndex), " consistency ",
                          num(nc.perPipe.fluid.consistency));
            if (nc.perPipe.compositionOn)
                log::info(compositionDecl(nc.perPipe)); // G1

            TwoFluidNetwork net(nc);
            auto emit = [&](const TwoFluidNetwork& n)
            {
                log::info("t ", num(n.time()), " net");
                for (const std::string& j : n.junctionNames())
                    log::info("junction ", j, " p ",
                              num(n.junctionPressure(j)), " net ",
                              num(n.junctionNetRate(j)));
                for (const std::string& pn : n.pipeNames())
                {
                    const TwoFluidPipe& s = n.pipe(pn);
                    log::info("pipe ", pn, " vm ", num(s.mixtureVelocity()),
                              " dp ", num(s.pressureDrop()));
                    // 5.77: the per-pipe latch flag, same wording as the
                    // serial report but named (a network has many rows)
                    if (nc.perPipe.fluid.yieldStress > 0)
                        log::info("static hold ", pn, " ",
                                  s.staticHold() ? 1 : 0);
                    const std::vector<double> p = s.nodePressure();
                    for (std::size_t q = 0; q < p.size(); ++q)
                        log::info("pnode ", pn, " ", q, " p ", num(p[q]));
                    const std::vector<double> a = s.cellAlpha();
                    const std::vector<double> w = s.cellWater();
                    const std::vector<double> vsl = s.cellVsl();
                    const std::vector<double> vsg = s.cellVsg();
                    const std::vector<double> T = s.cellTemperature();
                    const std::vector<double> ep
                        = nc.perPipe.fluid.gasR > 0 ? s.cellPressure()
                                                    : std::vector<double>();
                    const std::vector<double> hy
                        = nc.perPipe.hydrateOn ? s.cellHydrate()
                                               : std::vector<double>();
                    const std::vector<double> wx
                        = nc.perPipe.waxOn ? s.cellWax()
                                           : std::vector<double>();
                    const std::vector<double> de
                        = nc.perPipe.waxOn && nc.perPipe.wax.feedback
                              ? s.cellEffectiveDiameter()
                              : std::vector<double>();
                    const std::vector<double> az
                        = nc.perPipe.asphalteneOn ? s.cellAsphaltene()
                                                  : std::vector<double>();
                    const std::vector<double> sc
                        = nc.perPipe.scaleOn ? s.cellScale()
                                             : std::vector<double>();
                    const std::vector<double> ent
                        = nc.perPipe.entrainmentOn ? s.cellEntrained()
                                                   : std::vector<double>();
                    const std::vector<double> inh
                        = nc.perPipe.inhibitorOn ? s.cellInhibitor()
                                                 : std::vector<double>();
                    // O4 (5.163): same gated rigorous-equilibrium columns
                    // as the serial report
                    const std::vector<double> teq = nc.perPipe
                        .inhibitorEquilibriumOn
                            ? s.cellEquilibriumTemperature()
                            : std::vector<double>();
                    const std::vector<double> dts = nc.perPipe
                        .inhibitorEquilibriumOn
                            ? s.cellSubcooling()
                            : std::vector<double>();
                    for (int i = 0; i < s.sections(); ++i)
                        log::info("pcell ", pn, " ", i, " alpha ",
                                  num(a[i]), " water ", num(w[i]), " vsl ",
                                  num(vsl[i]), " vsg ", num(vsg[i]),
                                  T.empty() ? "" : (" T " + num(T[i])),
                                  ep.empty() ? "" : (" p " + num(ep[i])),
                                  hy.empty() ? "" : (" h " + num(hy[i])),
                                  wx.empty() ? "" : (" wax " + num(wx[i])),
                                  de.empty() ? "" : (" dEff " + num(de[i])),
                                  az.empty() ? "" : (" asph " + num(az[i])),
                                  sc.empty() ? "" : (" scale " + num(sc[i])),
                                  ent.empty() ? ""
                                              : (" entr " + num(ent[i])),
                                  inh.empty() ? ""
                                              : (" inh " + num(inh[i])),
                                  teq.empty() ? ""
                                              : (" Teq " + num(teq[i])
                                                 + " dTsub " + num(dts[i])),
                                  compositionReport(s, i));
                }
            };

            emit(net);
            for (int i = 1; i <= nSteps; ++i)
            {
                net.step(dt);
                if (i % report == 0 || i == nSteps) emit(net);
            }
            log::info("twoFluid network solved: steps=", nSteps);
            return 0;
        }

        const TwoFluidCase c = loadTwoFluidCase(caseDict);

        int sections = 0;
        for (const TwoFluidGeo& g : c.pipes) sections += g.sections;
        log::info("twoFluid: ", sections, " sections, dt=", num(dt),
                  " endTime=", num(endTime));
        // 5.75 (A5): a declared liquid rheology names itself and its
        // Rabinowitsch-corrected consistency, so the report says which law
        // the wall stress came from
        if (c.fluid.flowIndex > 0)
            log::info("twoFluid liquid rheology: flowIndex ",
                      num(c.fluid.flowIndex), " consistency ",
                      num(c.fluid.consistency), " (K' ",
                      num(c.fluid.consistency * std::pow(
                          (3.0 * c.fluid.flowIndex + 1.0)
                              / (4.0 * c.fluid.flowIndex),
                          c.fluid.flowIndex)), ")");
        // emulsion/inhibitor declarations name themselves, so the report
        // says which mixing law and which depression line the run used
        if (c.emulsionOn)
            log::info("twoFluid emulsion: model ",
                      c.emulsion.model == 1 ? "brinkman" : "richardson",
                      " muWater ", num(c.emulsion.muWater), " muOil ",
                      num(c.fluid.muL));
        if (c.inhibitorOn)
            log::info("twoFluid inhibitor: kDepression ",
                      num(c.inhibitor.kDepression), " initial ",
                      num(c.inhibitor.initial));
        // O4 (5.163): the rigorous activity root names its own constant
        if (c.inhibitorEquilibriumOn)
            log::info("twoFluid inhibitorEquilibrium: activityK ",
                      num(c.inhibitorEquilibrium.activityK), " basis ",
                      c.inhibitorEquilibrium.basisWater ? "water"
                                                        : "liquid");
        if (c.sandOn)
            log::info("twoFluid sand: vt ", num(c.sand.terminalVelocity),
                      " uCrit ", num(c.sand.criticalVelocity), " erodeRate ",
                      num(c.sand.erodeRate), " feedback ",
                      c.sand.feedback ? "true" : "false");
        // 5.229 (H1b): the distributed phase momentum names its own mode
        if (c.sixEquationOn)
            log::info("twoFluid sixEquation: writeBack ",
                      c.sixEquation.writeBack ? "true" : "false",
                      c.sixEquation.hasU0
                          ? (" uG0 " + num(c.sixEquation.uG0)
                             + " uL0 " + num(c.sixEquation.uL0))
                          : std::string());
        if (c.compositionOn) log::info(compositionDecl(c)); // G1

        TwoFluidPipe sim(c);

        // O2 (5.161): cursor over the slug tracking event journal, so every
        // GEN/MERGE/OUT is reported even between two sampled frames
        std::size_t slugTrackCursor = 0;
        // O3 (5.162): cursor over the purge wave TRIGGER/COMPLETE journal
        std::size_t purgeCursor = 0;
        // E6 (5.176): cursor over the pig ARRIVE journal
        std::size_t pigCursor = 0;

        auto emit = [&](const TwoFluidPipe& s)
        {
            log::info("t ", num(s.time()), " vm ", num(s.mixtureVelocity()),
                      " dp ", num(s.pressureDrop()));
            // 5.77: the latch the step that produced this sample started
            // with -- 1 means this velocity is the row's own u = 0, not a
            // value the moving-wall law was solved for. At t = 0 no step
            // has been taken, so the flag reads 0 even for a held case.
            if (c.fluid.yieldStress > 0)
                log::info("static hold ", s.staticHold() ? 1 : 0);
            const std::vector<double> p = s.nodePressure();
            for (std::size_t j = 0; j < p.size(); ++j)
                log::info("node ", j, " p ", num(p[j]));
            const std::vector<double> a = s.cellAlpha();
            const std::vector<double> w = s.cellWater();
            const std::vector<double> vsl = s.cellVsl();
            const std::vector<double> vsg = s.cellVsg();
            const std::vector<double> sl
                = c.closure.model == 2 || c.closure.model == 3
                      ? s.cellSlip()
                      : std::vector<double>();
            const std::vector<double> T = s.cellTemperature();
            const std::vector<double> ep
                = c.fluid.gasR > 0 ? s.cellPressure() : std::vector<double>();
            const std::vector<double> hy
                = c.hydrateOn ? s.cellHydrate() : std::vector<double>();
            const std::vector<double> wx
                = c.waxOn ? s.cellWax() : std::vector<double>();
            const std::vector<double> de
                = (c.waxOn && c.wax.feedback)
                      || (c.sandOn && c.sand.feedback)
                      || c.blockageOn || c.hydraulicFeedbackOn
                          ? s.cellEffectiveDiameter()
                          : std::vector<double>();
            // 5.227: the hydrate deposition layer thickness (m) rides the
            // cell rows as a gated `blk` column - absent unless blockage{}
            const std::vector<double> bl
                = c.blockageOn ? s.cellBlockageThickness()
                               : std::vector<double>();
            // 5.230: the scheduled wall temperature the step saw, gated on
            // ambientSchedule{}
            const std::vector<double> wt
                = c.ambientScheduleOn ? s.cellWallTemperature()
                                      : std::vector<double>();
            const std::vector<double> az
                = c.asphalteneOn ? s.cellAsphaltene()
                                 : std::vector<double>();
            const std::vector<double> sc
                = c.scaleOn ? s.cellScale() : std::vector<double>();
            const std::vector<double> ent
                = c.entrainmentOn ? s.cellEntrained()
                                  : std::vector<double>();
            const std::vector<double> inh
                = c.inhibitorOn ? s.cellInhibitor() : std::vector<double>();
            // O4 (5.163): the rigorous inhibited equilibrium line and the
            // subcooling the kinetics see - gated columns, absent when the
            // block is not declared
            const std::vector<double> teq
                = c.inhibitorEquilibriumOn ? s.cellEquilibriumTemperature()
                                           : std::vector<double>();
            const std::vector<double> dts
                = c.inhibitorEquilibriumOn ? s.cellSubcooling()
                                           : std::vector<double>();
            const std::vector<double> me
                = c.emulsionOn ? s.cellLiquidViscosity()
                               : std::vector<double>();
            const std::vector<double> sa
                = c.sandOn ? s.cellSand() : std::vector<double>();
            const std::vector<double> be
                = c.sandOn ? s.cellBed() : std::vector<double>();
            // E7b (5.182): the opt-in SLUG/FILM partition rides the cell
            // rows as a numeric `zone` key (1 = the cell centre is inside
            // a tracked slug span, 0 = film), absent unless slugTracking
            // declared zoning. The existing keys - including the closure
            // slip `s` - are untouched.
            const std::vector<int> zoned
                = c.slugTrackingOn && c.slugTracking.zoning
                      ? s.cellZones()
                      : std::vector<int>();
            // E7c (5.183): the linear acoustic layer's gated report columns
            // - the cell acoustic pressure rides the cell rows as a `pA`
            // key and the face acoustic velocities as one `aunode` line per
            // face - absent unless acoustics{} was declared and enabled.
            const std::vector<double> pAc
                = c.acousticsOn ? s.cellAcousticPressure()
                                : std::vector<double>();
            const std::vector<double> uAc
                = c.acousticsOn ? s.faceAcousticVelocity()
                                : std::vector<double>();
            // 5.229 (H1b): the distributed phase velocities ride the cell
            // rows as `uG`/`uL` keys - absent unless sixEquation{} was
            // declared and enabled
            const std::vector<double> ug
                = c.sixEquationOn ? s.cellGasVelocity()
                                  : std::vector<double>();
            const std::vector<double> ul
                = c.sixEquationOn ? s.cellLiquidVelocity()
                                  : std::vector<double>();
            // one face acoustic velocity line per face (same idiom as the
            // node pressure lines), gated on the acoustics{} block
            if (!uAc.empty())
                for (std::size_t j = 0; j < uAc.size(); ++j)
                    log::info("aunode ", j, " u ", num(uAc[j]));
            // face flux columns (EOS-gated): this step's face liquid
            // volume flux F and gas volume flux G, one line per face - an
            // independent checker asserts the model identity F + G = vM
            // from the report alone
            if (c.fluid.gasR > 0)
            {
                const std::vector<double> fl = s.faceLiquidFlux();
                const std::vector<double> fg = s.faceGasVolumeFlux();
                for (std::size_t j = 0; j < fl.size(); ++j)
                    log::info("face ", j, " F ", num(fl[j]), " G ",
                              num(fg[j]));
            }
            for (int i = 0; i < s.sections(); ++i)
            {
                log::info("cell ", i, " alpha ", num(a[i]), " water ",
                          num(w[i]), " vsl ", num(vsl[i]), " vsg ",
                          num(vsg[i]),
                          T.empty() ? "" : (" T " + num(T[i])),
                          sl.empty() ? "" : (" s " + num(sl[i])),
                          zoned.empty()
                              ? "" : (" zone " + std::to_string(zoned[i])),
                          pAc.empty() ? "" : (" pA " + num(pAc[i])),
                          ep.empty() ? "" : (" p " + num(ep[i])),
                          hy.empty() ? "" : (" h " + num(hy[i])),
                          wx.empty() ? "" : (" wax " + num(wx[i])),
                          de.empty() ? "" : (" dEff " + num(de[i])),
                          bl.empty() ? "" : (" blk " + num(bl[i])),
                          az.empty() ? "" : (" asph " + num(az[i])),
                          sc.empty() ? "" : (" scale " + num(sc[i])),
                          ent.empty() ? "" : (" entr " + num(ent[i])),
                          inh.empty() ? "" : (" inh " + num(inh[i])),
                          teq.empty() ? ""
                                      : (" Teq " + num(teq[i])
                                         + " dTsub " + num(dts[i])),
                          me.empty() ? "" : (" muEff " + num(me[i])),
                          sa.empty() ? "" : (" sand " + num(sa[i])),
                          be.empty() ? "" : (" bed " + num(be[i])),
                          wt.empty() ? "" : (" wallT " + num(wt[i])),
                          ug.empty() ? ""
                                     : (" uG " + num(ug[i])
                                        + " uL " + num(ul[i])),
                          compositionReport(s, i));
            }
            // 5.195 (G10 fix): the clamp that used to drop sand silently. It
            // prints only when it fired (the cumulative volume the transport
            // step found over-full and put into the bed), so an ordinary case
            // report is byte-identical to the pre-fix one.
            if (c.sandOn && s.sandClampToBed() > 0.0)
                log::info("sandclamp vol ", num(s.sandClampToBed()));
            // 5.187 (F-class): the transient risk report - one erosion line,
            // plus the hydrate (non-blocking) and corrosion lines when their
            // channels are declared. Pure readout: absent, nothing prints and
            // the frame is the untouched pre-5.187 report.
            if (c.riskOn)
            {
                const TwoFluidRiskReadout rk = s.riskReadout();
                log::info("risk erosion seg ", rk.erosionSeg,
                          " vm ", num(rk.erosionVm),
                          " rhom ", num(rk.erosionRhoM),
                          " velim ", num(rk.erosionVelim),
                          " ratio ", num(rk.erosionRatio),
                          " over ", rk.erosionOverCount,
                          rk.erosionRatio > 1.0 ? " OVER" : " OK");
                if (rk.hydrateOn)
                    log::info("risk hydrate seg ", rk.hydrateSeg,
                              " p ", num(rk.hydrateP),
                              " T ", num(rk.hydrateT),
                              " Teq ", num(rk.hydrateTeq),
                              " dTsub ", num(rk.hydrateSubcool),
                              " grade ", rk.hydrateGrade);
                if (rk.corrosionOn)
                    log::info("risk corrosion co2 ", num(rk.co2Partial),
                              " h2s ", num(rk.h2sPartial),
                              " band ", rk.co2Band,
                              " sour ", rk.sour ? 1 : 0);
            }
            if (c.slugOn)
            {
                const std::vector<TwoFluidSlugBody> bodies = s.slugBodies();
                log::info("slug count ", bodies.size());
                for (std::size_t k = 0; k < bodies.size(); ++k)
                    log::info("slugbody ", k, " rear ",
                              num(bodies[k].rear), " front ",
                              num(bodies[k].front), " length ",
                              num(bodies[k].length), " inv ",
                              num(bodies[k].inventory));
            }
            if (c.slugTrackingOn)
            {
                // O2 (5.161): the live Lagrangian slugs (state lines, one
                // per object) followed by the event journal lines not yet
                // reported; the event token separates the two kinds
                const std::vector<TwoFluidSlugTrack>& objs = s.slugTracks();
                log::info("slugtrack count ", objs.size());
                for (const TwoFluidSlugTrack& o : objs)
                    log::info("slugtrack ", o.id, " t ", num(s.time()),
                              " xR ", num(o.xR), " xF ", num(o.xF), " L ",
                              num(o.xF - o.xR), " inv ", num(o.inv),
                              " uF ", num(o.uF), " uR ", num(o.uR),
                              " event NONE");
                const std::vector<TwoFluidSlugTrackEvent>& evs
                    = s.slugTrackEvents();
                for (; slugTrackCursor < evs.size(); ++slugTrackCursor)
                {
                    const TwoFluidSlugTrackEvent& e = evs[slugTrackCursor];
                    const char* kind = e.kind == 1 ? "GEN"
                                     : e.kind == 2 ? "MERGE" : "OUT";
                    std::string line = "slugtrack " + std::to_string(e.id)
                        + " t " + num(e.t) + " xR " + num(e.xR)
                        + " xF " + num(e.xF) + " L " + num(e.xF - e.xR)
                        + " inv " + num(e.inv) + " uF " + num(e.uF)
                        + " uR " + num(e.uR) + " event " + kind;
                    if (e.kind == 3)
                        line += " out " + num(e.outVol) + " duration "
                              + num(e.duration);
                    log::info(line);
                }
            }
            if (c.purgeWaveOn)
            {
                // O3 (5.162): the purge wave state line every frame (the
                // phase token is NONE while the column accumulates) and the
                // TRIGGER/COMPLETE journal lines not yet reported, same
                // cursor idiom as slugtrack; COMPLETE adds the cycle
                // accounting (out pulse volume, duration, period)
                const TwoFluidPurgeState& ps = s.purgeState();
                const char* ph = ps.phase == 1 ? "TRIGGER"
                             : ps.phase == 2 ? "PROPAGATE"
                             : ps.phase == 3 ? "COMPLETE" : "NONE";
                log::info("purge t ", num(s.time()), " phase ", ph,
                          " x_wave ", num(ps.xWave), " u_wave ",
                          num(ps.uWave), " cap ", num(ps.capPressure),
                          " pbase ", num(ps.basePressure), " head ",
                          num(ps.head), " cycles ", ps.cycles, " out ",
                          num(ps.outVol));
                const std::vector<TwoFluidPurgeEvent>& evs = s.purgeEvents();
                for (; purgeCursor < evs.size(); ++purgeCursor)
                {
                    const TwoFluidPurgeEvent& e = evs[purgeCursor];
                    std::string line = "purge t " + num(e.t) + " x_wave "
                        + num(e.xWave) + " head " + num(e.head) + " cap "
                        + num(e.capPressure) + " phase "
                        + (e.kind == 1 ? "TRIGGER" : "COMPLETE");
                    if (e.kind == 2)
                        line += " cycle " + std::to_string(e.cycle)
                              + " out " + num(e.outVol) + " duration "
                              + num(e.duration) + " period " + num(e.period);
                    log::info(line);
                }
            }
            if (c.piggingOn)
            {
                // E6 (5.176): one pig state line per sampled frame (RUN while
                // travelling, NONE after arrival) and each ARRIVE journal line
                // once, even between frames -- the same cursor idiom as the
                // slug tracking / purge reports.
                const TwoFluidPigState& ps = s.pigState();
                log::info("pig t ", num(s.time()), " x ", num(ps.x),
                          " v ", num(ps.v), " holdupFront ",
                          num(ps.holdupFront), " outVol ", num(ps.outVol),
                          " event ", ps.arrived ? "NONE" : "RUN");
                const std::vector<TwoFluidPigEvent>& evs = s.pigEvents();
                for (; pigCursor < evs.size(); ++pigCursor)
                {
                    const TwoFluidPigEvent& e = evs[pigCursor];
                    log::info("pig t ", num(e.t), " x ", num(e.x),
                              " v ", num(e.v), " holdupFront ",
                              num(e.holdupFront), " outVol ", num(e.outVol),
                              " event ARRIVE");
                }
            }
            if (c.terrainOn)
            {
                const std::vector<TwoFluidTerrainDip> dips = s.terrainDips();
                log::info("terrain count ", dips.size());
                for (std::size_t k = 0; k < dips.size(); ++k)
                    log::info("terrain dip ", k, " s ",
                              num(dips[k].position), " depth ",
                              num(dips[k].depth), " len ",
                              num(dips[k].basinLength), " inv ",
                              num(dips[k].inventory), " cap ",
                              num(dips[k].capacity), " fill ",
                              num(dips[k].fill), " head ",
                              num(dips[k].head),
                              dips[k].slugReady ? " SLUG-READY" : " POOLING");
            }
            if (c.terrainOn && c.terrain.slugVelocity > 0)
            {
                const std::vector<TwoFluidTerrainSlug> slugs = s.terrainSlugs();
                log::info("terrain deposits ", s.terrainDeposits());
                for (std::size_t k = 0; k < slugs.size(); ++k)
                    log::info("terrain slug ", k, " front ",
                              num(slugs[k].front), " vol ",
                              num(slugs[k].volume), " from ",
                              num(slugs[k].capturedFrom),
                              c.terrain.dynamicUs
                                  ? (" us " + num(slugs[k].velocity)
                                     + " pcap " + num(slugs[k].capPressure)
                                     + " cfront " + num(slugs[k].calcFront))
                                  : std::string(),
                              slugs[k].active ? " ACTIVE" : " IDLE");
            }
            if (!T.empty())
                log::info("thermal: qWall ", num(s.wallHeatRate()),
                          " qInlet ", num(s.inletHeatRate()));
            if (c.formationOn)
            {
                log::info("formation: qFar ", num(s.formationHeatRate()),
                          " qIntake ", num(s.formationIntakeRate()),
                          " qStore ", num(s.formationStorageRate()),
                          " farRise ", num(s.formationFarRise()));
                if (c.formation.semiInfinite)
                {
                    // probes of the E1 implementation itself; the checker
                    // weighs these against an independent quadrature.
                    log::info("formation e1: 0.1 ",
                              num(exponentialIntegralE1(0.1)),
                              " 1 ", num(exponentialIntegralE1(1.0)),
                              " 10 ", num(exponentialIntegralE1(10.0)),
                              " 50 ", num(exponentialIntegralE1(50.0)));
                }
                for (int i = 0; i < s.sections(); ++i)
                {
                    const std::vector<double> fp = s.formationProfile(i);
                    std::string line = "form " + std::to_string(i);
                    for (std::size_t j = 0; j < fp.size(); ++j)
                        line += " " + num(fp[j]);
                    log::info(line);
                }
            }
        };

        emit(sim);
        for (int i = 1; i <= nSteps; ++i)
        {
            sim.step(dt);
            if (i % report == 0 || i == nSteps) emit(sim);
        }
        log::info("twoFluid solved: steps=", nSteps);
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("1d-twofluid: ", e.what());
        return 1;
    }
}
