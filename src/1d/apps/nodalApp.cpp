#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/Nodal.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

using namespace cf;

namespace
{
// log::info() stringifies with ostringstream defaults (6 digits); results
// are checked against closed-form values, so emit 12 significant digits.
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

const char* regimeName(int r)
{
    static const char* names[] = {"segregated", "transition", "intermittent",
                                  "distributed"};
    return (r >= 0 && r < 4) ? names[r] : "single-phase";
}

const char* iprFormName(IpR::Form f)
{
    switch (f)
    {
        case IpR::Form::linear: return "linear";
        case IpR::Form::vogel: return "vogel";
        case IpR::Form::jones: return "jones";
    }
    return "?";
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-1d-nodal <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const NodalCase c = loadNodalCaseFromDict(caseDict);

        log::info("nodal: ipr ", iprFormName(c.ipr.form),
                  " p_res ", num(c.ipr.reservoirPressure),
                  " qmax ", num(c.ipr.rateMax()),
                  ", well ", num(c.well.length), " m x ", num(c.well.diameter),
                  " m, wellhead ", num(c.wellheadPressure), " Pa");

        if (!c.laterals.empty())
        {
            // 5.157 (P5): multi-lateral completion. The lateral lines announce
            // the switch, so the main ipr{} above is not silently unused: the
            // junction rates come from the laterals, the main bore is the lift.
            std::string names;
            for (const LateralSpec& lat : c.laterals)
            {
                if (!names.empty()) names += " ";
                names += lat.name;
            }
            log::info("nodal: laterals ", c.laterals.size(), " (", names,
                      "): the main bore lifts the summed lateral flow, the "
                      "main ipr does not feed the junction");

            const MultiLateralResult r = solveNodalLaterals(
                c.fluid, c.well, c.wellheadPressure, c.laterals, c.controls);

            log::info("column p ", num(r.main.staticBottomPressure),
                      " (the rate-0 static column)");
            for (const LateralResult& lr : r.laterals)
                log::info("lateral ", lr.name, " q ", num(lr.rate),
                          " pwf ", num(lr.bottomHolePressure), " holdup ",
                          num(lr.holdup), " regime ", regimeName(lr.regime),
                          " flows ", lr.flows ? "true" : "false");
            log::info("well q ", num(r.main.rate), " pwf ",
                      num(r.main.bottomHolePressure), " pth ",
                      num(c.wellheadPressure));
            log::info("well holdup ", num(r.main.holdup), " regime ",
                      regimeName(r.main.regime));
            log::info("nodal solved: iterations=", r.iterations,
                      " flows=", r.main.flows ? "true" : "false",
                      r.converged ? " (converged)" : " (NOT converged)");
            return r.converged ? 0 : 1;
        }

        const NodalResult r = solveNodalAnalysis(
            c.ipr, c.fluid, c.well, c.wellheadPressure, c.controls);

        log::info("column p ", num(r.staticBottomPressure),
                  " (the rate-0 static column)");
        log::info("well q ", num(r.rate), " pwf ", num(r.bottomHolePressure),
                  " pth ", num(c.wellheadPressure));
        log::info("well holdup ", num(r.holdup), " regime ",
                  regimeName(r.regime));
        log::info("nodal solved: iterations=", r.iterations,
                  " flows=", r.flows ? "true" : "false",
                  r.converged ? " (converged)" : " (NOT converged)");

        return r.converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("1d-nodal: ", e.what());
        return 1;
    }
}
