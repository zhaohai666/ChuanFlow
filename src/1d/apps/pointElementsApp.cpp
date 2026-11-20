#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/Compressor.h"
#include "cf/1d/JetPump.h"
#include "cf/1d/PlungerGasLift.h"
#include "cf/1d/ProgressiveCavityPump.h"
#include "cf/1d/Separator.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

using namespace cf;

namespace
{
// Closed-form elements: results are checked against recomputed algebra to
// 1e-9, so emit 12 significant digits (the steady app's convention).
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

const Dictionary* section(const Dictionary& root, const char* name)
{
    static Dictionary empty;
    if (!root.has(name)) return &empty;
    return &root.getDict(name);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-1d-point-elements <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const Dictionary& jps = *section(caseDict, "jetPumps");
        const Dictionary& pcps = *section(caseDict, "pcpPumps");
        const Dictionary& pgls = *section(caseDict, "plungerGasLifts");
        const Dictionary& seps = *section(caseDict, "separators");
        const Dictionary& cmps = *section(caseDict, "compressors");

        log::info("point elements: ", jps.keys().size(), " jet pumps, ",
                  pcps.keys().size(), " pcp pumps, ", pgls.keys().size(),
                  " plunger gas lifts, ", seps.keys().size(), " separators, ",
                  cmps.keys().size(), " compressors");

        for (const std::string& name : jps.keys())
        {
            const Dictionary& d = jps.getDict(name);
            JetPump pump;
            pump.nozzleArea = d.getScalarOr("nozzleArea", 0.0);
            pump.throatArea = d.getScalarOr("throatArea", 0.0);
            pump.throatLength = d.getScalarOr("throatLength", 0.0);
            pump.frictionFactor = d.getScalarOr("frictionFactor", 0.0);
            pump.fluidDensity = d.getScalar("rho");

            JetPumpOperation op;
            op.motiveRate = d.getScalarOr("motiveRate", 0.0);
            op.dischargePressure = d.getScalarOr("dischargePressure", 0.0);
            op.imposeSuctionRate = d.has("suctionRate");
            if (op.imposeSuctionRate)
            {
                op.suctionRate = d.getScalar("suctionRate");
            }
            else
            {
                op.suctionPressure = d.getScalarOr("suctionPressure", 0.0);
            }

            static const char* stateName[] = {"PUMPING", "NO-SUCTION",
                                              "FREE-JET"};
            const JetPumpSolution s = jetPumpSolve(pump, op);
            log::info("jetPump ", name, " state ", stateName[static_cast<int>(s.state)],
                      " m ", num(s.areaRatio), " r ", num(s.flowRatio),
                      " vn ", num(s.nozzleVelocity), " v3 ", num(s.throatVelocity),
                      " qn ", num(s.motiveRate), " qs ", num(s.suctionRate),
                      " qd ", num(s.dischargeRate),
                      " ps ", num(s.suctionPressure),
                      " pd ", num(s.dischargePressure),
                      " psupply ", num(s.supplyPressure),
                      " N ", num(s.headRatio), " eff ", num(s.efficiency),
                      " pu ", num(s.usefulPower), " pm ", num(s.motivePower),
                      " resid ", num(s.momentumResidual));
        }

        for (const std::string& name : pcps.keys())
        {
            const Dictionary& d = pcps.getDict(name);
            ProgressiveCavityPump pump;
            pump.cavitiesPerRev = d.getScalarOr("cavitiesPerRev", 0.0);
            pump.cavityVolume = d.getScalarOr("cavityVolume", 0.0);
            pump.rpm = d.getScalarOr("rpm", 0.0);
            pump.slip = d.getScalarOr("slip", 0.0);
            const double rho = d.getScalar("rho");

            ProgressiveCavityPumpCard c;
            if (d.has("dischargePressure") && d.has("intakePressure"))
                c = progressiveCavityPumpCardOf(pump, rho,
                                                d.getScalar("dischargePressure"),
                                                d.getScalar("intakePressure"));
            else
                c = progressiveCavityPumpCardOfLoad(pump, rho,
                                                    d.getScalarOr("pressureRise", 0.0));

            static const char* stateName[] = {"PUMPING", "NO-LOAD", "STALLED"};
            log::info("pcp ", name, " state ", stateName[static_cast<int>(c.state)],
                      " qth ", num(c.theoreticalRate), " slip ", num(c.slipRate),
                      " q ", num(c.netRate), " m ", num(c.massRate),
                      " dp ", num(c.pressureRise), " ev ", num(c.volumetricEfficiency),
                      " power ", num(c.hydraulicPower));
        }

        for (const std::string& name : pgls.keys())
        {
            const Dictionary& d = pgls.getDict(name);
            PlungerGasLift pgl;
            pgl.gasVolumePerCycle = d.getScalarOr("gasVolumePerCycle", 0.0);
            pgl.referencePressure = d.getScalarOr("referencePressure", 101325.0);
            pgl.gasDensityAtReference = d.getScalarOr("gasDensityAtReference", 0.0);
            pgl.cyclePeriod = d.getScalarOr("cyclePeriod", 0.0);
            pgl.openWindow = d.getScalarOr("openWindow", 0.0);
            pgl.releaseDepth = d.getScalarOr("releaseDepth", 0.0);
            pgl.casingSurfacePressure = d.getScalarOr("casingSurfacePressure", 0.0);
            pgl.tubingSurfacePressure = d.getScalarOr("tubingSurfacePressure", 0.0);
            pgl.liquidDensity = d.getScalarOr("liquidDensity", 0.0);
            pgl.valve.domePressure = d.getScalarOr("liftDomePressure", 0.0);
            pgl.valve.areaRatio = d.getScalarOr("liftAreaRatio", 0.0);
            pgl.valve.springPressure = d.getScalarOr("liftSpringPressure", 0.0);
            pgl.valve.portArea = 0.25 * M_PI
                               * d.getScalarOr("liftPortDiameter", 0.0)
                               * d.getScalarOr("liftPortDiameter", 0.0);
            pgl.valve.dischargeCoeff = d.getScalarOr("liftDischargeCoeff", 0.0);
            pgl.valve.kappa = d.getScalarOr("liftKappa", 0.0);

            static const char* stateName[] = {"OPEN", "VALVE-LIMITED",
                                              "CLOSED", "IN-SPREAD"};
            const PlungerGasLiftSchedule s = plungerGasLiftSolve(pgl);
            log::info("pgl ", name, " state ", stateName[static_cast<int>(s.state)],
                      " pc ", num(s.injectionPressureAtDepth),
                      " pt ", num(s.tubingPressureAtDepth),
                      " pvo ", num(s.openingPressure),
                      " pvc ", num(s.closingPressure),
                      " choked ", num(s.choked ? 1.0 : 0.0),
                      " qw ", num(s.windowRateReference),
                      " qwd ", num(s.windowRateDownhole),
                      " qavg ", num(s.averageRateReference),
                      " mw ", num(s.massRateInWindow),
                      " mreq ", num(s.requiredMassRate),
                      " mcap ", num(s.nozzleMassRateCapacity),
                      " volref ", num(s.deliveredVolumeReference),
                      " voldh ", num(s.deliveredVolumeDownhole));

            // The cycle window's own footprint: probe the schedule at the
            // declared times so the checker can verify open/closed against
            // (t mod T) < w without trusting the solver's bookkeeping.
            if (d.has("probeTimes"))
            {
                for (const double t : d.getScalarList("probeTimes"))
                {
                    log::info("pglWindow ", name, " t ", num(t), " open ",
                              num(plungerGasLiftWindowOpenAt(pgl, t) ? 1.0 : 0.0));
                }
            }
        }

        for (const std::string& name : seps.keys())
        {
            const Dictionary& d = seps.getDict(name);
            const Dictionary& bd = d.getDict("fluid").getDict("blackoil");
            pvt::BlackOilFluid fluid;
            fluid.gasGravity = bd.getScalarOr("gasGravity", fluid.gasGravity);
            fluid.oilDensityStd = bd.getScalarOr("oilDensity", fluid.oilDensityStd);
            fluid.gorStd = bd.getScalarOr("gor", fluid.gorStd);
            fluid.waterCut = bd.getScalarOr("waterCut", fluid.waterCut);
            fluid.waterDensity = bd.getScalarOr("waterDensity", fluid.waterDensity);

            const SeparatorSolution s = separatorSolve(
                fluid, d.getScalar("feedRate"), d.getScalar("pressure"),
                d.getScalar("temperature"), d.getScalarOr("efficiency", 1.0),
                d.getScalarOr("dP", 0.0));
            log::info("sep ", name, " p ", num(s.pressure), " T ",
                      num(s.temperature), " eff ", num(s.efficiency),
                      " pout ", num(s.outletPressure), " rs ", num(s.rs),
                      " glib ", num(s.gasLiberated), " bo ", num(s.bo),
                      " bg ", num(s.bg), " z ", num(s.z), " rhogs ",
                      num(s.rhoGasStd), " mf ", num(s.feedMassRate), " mli ",
                      num(s.liquidMassRateIdeal), " ment ",
                      num(s.entrainedMassRate), " ml ", num(s.liquidMassRate),
                      " mg ", num(s.gasMassRate), " closure ", num(s.closure));
        }

        for (const std::string& name : cmps.keys())
        {
            const Dictionary& d = cmps.getDict(name);
            Compressor c;
            c.pressureRatio = d.getScalarOr("pressureRatio", 0.0);
            c.shutoffRatio = d.getScalarOr("shutoffRatio", c.pressureRatio);
            c.ratedRate = d.getScalarOr("ratedRate", 0.0);
            c.polytropicHead = d.getScalarOr("polytropicHead", 0.0);
            c.polytropicEfficiency = d.getScalarOr("polytropicEfficiency", 1.0);
            c.kappa = d.getScalarOr("kappa", 0.0);
            c.specificGasConstant = d.getScalarOr("specificGasConstant", 0.0);
            c.compressibility = d.getScalarOr("compressibility", 1.0);

            const CompressorSolution s = compressorSolve(
                c, d.getScalar("suctionPressure"),
                d.getScalar("suctionTemperature"), d.getScalar("rate"));
            log::info("cmp ", name, " mode ", s.mode == 0 ? "curve" : "head",
                      " q ", num(s.rate), " p1 ", num(s.pressureIn), " t1 ",
                      num(s.temperatureIn), " pr ", num(s.ratio), " p2 ",
                      num(s.pressureOut), " m ", num(s.massRate), " h ",
                      num(s.head), " w0 ", num(s.powerIdeal), " w ",
                      num(s.powerShaft), " t2 ", num(s.temperatureOut),
                      " eta ", num(c.polytropicEfficiency));
        }

        log::info("point elements solved");
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("1d-point-elements: ", e.what());
        return 1;
    }
}
