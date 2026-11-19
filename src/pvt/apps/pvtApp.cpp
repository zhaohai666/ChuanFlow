#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/pvt/ComponentData.h"
#include "cf/pvt/Enthalpy.h"
#include "cf/pvt/PvtCase.h"
#include "cf/pvt/HydrateEquilibrium.h"
#include "cf/pvt/PvtTable.h"
#include "cf/pvt/SurfaceTension.h"
#include "cf/pvt/ThermalConductivity.h"
#include "cf/pvt/Viscosity.h"

using namespace cf;
using namespace cf::pvt;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

void reportOne(const FlashCase& c, double p, bool wantVisc,
               const std::vector<double>* parachors)
{
    const FlashResult r = tpFlash(c.feed, c.temperature, p, c.controls);
    log::info("pvt flash T ", num(c.temperature), " p ", num(p),
              " beta ", num(r.vaporFraction), " single ", r.singlePhase ? 1 : 0,
              " iterations ", r.iterations,
              " converged ", r.converged ? 1 : 0, " residual ", num(r.residual));
    log::info("pvt liquid rho ", num(r.liquidMassDensity), " M ", num(r.liquidMolarMass));
    log::info("pvt vapour rho ", num(r.vapourMassDensity), " M ", num(r.vapourMolarMass));
    for (std::size_t i = 0; i < c.feed.size(); ++i)
    {
        log::info("pvt x ", c.feed.components()[i].name, " ", num(r.liquidComposition[i]));
        log::info("pvt y ", c.feed.components()[i].name, " ", num(r.vapourComposition[i]));
    }
    if (!r.converged) log::warn("flash: ", r.message);

    // D7 (5.93): LBC residual viscosity of each flashed phase, gated on a
    // `viscosity` block so untouched cases keep their exact report. Vm comes
    // from the phase molar mass / density ( = PR molar volume) and Vc from
    // PR's critical compressibility, matching the convenience overload.
    if (wantVisc)
    {
        const std::vector<Component>& comps = c.feed.components();
        std::vector<double> MWs(comps.size()), Tcs(comps.size()),
            Pcs(comps.size()), Vcs(comps.size());
        for (std::size_t i = 0; i < comps.size(); ++i)
        {
            MWs[i] = comps[i].molarMass * 1.0e3;
            Tcs[i] = comps[i].Tc;
            Pcs[i] = comps[i].Pc;
            Vcs[i] = 0.3074 * kGasConstant * comps[i].Tc / comps[i].Pc;
        }
        auto phaseMu = [&](const std::vector<scalar>& zs, scalar M, scalar rho)
        {
            if (!(rho > 0.0)) return 0.0;
            return lbcViscosity(c.temperature, M / rho, zs, MWs, Tcs, Pcs, Vcs);
        };
        log::info("pvt viscosity liquid ",
                  num(phaseMu(r.liquidComposition, r.liquidMolarMass,
                              r.liquidMassDensity)),
                  " vapour ",
                  num(phaseMu(r.vapourComposition, r.vapourMolarMass,
                              r.vapourMassDensity)));
    }

    // D7 (5.94): Weinaug-Katz parachor surface tension of the liquid-vapour
    // interface, gated on a `surfaceTension` block. Vml/Vmg come from the phase
    // molar mass / density (= PR molar volume); a single-phase flash has no
    // interface and reports 0.
    if (parachors)
    {
        scalar sigma = 0.0;
        if (!r.singlePhase && r.liquidMassDensity > 0.0 && r.vapourMassDensity > 0.0)
            sigma = weinaugKatzSurfaceTension(
                *parachors, r.liquidMolarMass / r.liquidMassDensity,
                r.vapourMolarMass / r.vapourMassDensity,
                r.liquidComposition, r.vapourComposition);
        log::info("pvt surfaceTension ", num(sigma));
    }
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-pvt <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const FlashCase c = loadFlashCase(caseDict);
        const bool wantVisc = caseDict.has("viscosity");

        const bool wantSigma = caseDict.has("surfaceTension");
        std::vector<double> parachors;
        if (wantSigma)
        {
            const Dictionary& st = caseDict.getDict("surfaceTension");
            const std::string model = st.getWordOr("model", "weinaugKatz");
            if (model != "weinaugKatz")
                cfFatal("pvt: unknown surfaceTension model '" + model + "'");
            const Dictionary::Node& pn = st.find("parachors");
            if (pn.kind == Dictionary::Kind::word
                && st.getWord("parachors") == "fromTable")
            {
                // Opt-in built-in Weinaug-Katz (1943) parachors. Explicit
                // declaration (the list branch below) is untouched and stays
                // byte-identical; a table lookup is only ever requested by the
                // `fromTable` token.
                parachors = parachorsFromTableSi(c.feed.components());
                log::info("pvt surfaceTension parachors fromTable (Weinaug-Katz 1943)");
            }
            else if (pn.kind == Dictionary::Kind::list)
            {
                parachors = st.getScalarList("parachors");
            }
            else
            {
                cfFatal("pvt: surfaceTension parachors must be a list or 'fromTable'");
            }
            if (parachors.size() != c.feed.size())
                cfFatal("pvt: surfaceTension parachors must give one value per component");
        }
        const std::vector<double>* paraPtr = wantSigma ? &parachors : nullptr;

        log::info("pvt: ", c.feed.size(), " components, feed ", c.feed.describe(),
                  " molarMass ", num(c.feed.molarMass()));

        for (const double p : c.pressures) reportOne(c, p, wantVisc, paraPtr);

        if (c.computeBubblePoint)
        {
            const double pb = bubblePressure(c.feed, c.temperature);
            log::info("pvt bubble T ", num(c.temperature), " p ", num(pb));
            reportOne(c, pb, wantVisc, paraPtr);
        }

        // D13 increment 4: Prandtl thermal-conductivity closure. Gated on a
        // `thermalConductivity{}` block; absent means this report is untouched.
        if (caseDict.has("thermalConductivity"))
        {
            const Dictionary& tc = caseDict.getDict("thermalConductivity");
            const scalar mu = tc.getScalar("mu");
            const scalar cp = tc.getScalar("cp");
            const scalar Pr = tc.getScalar("Pr");
            log::info("pvt thermalConductivity k ",
                      num(thermalConductivityFromPr(mu, cp, Pr)));
        }

        // D13 increment 4: hydrate Clausius–Clapeyron dissociation line about
        // a declared reference point; gated, all constants declared in-case.
        if (caseDict.has("hydrateEquilibrium"))
        {
            const Dictionary& he = caseDict.getDict("hydrateEquilibrium");
            const scalar Tref = he.getScalar("Tref");
            const scalar pRef = he.getScalar("pRef");
            const scalar dH = he.getScalar("dH");
            for (const scalar T : he.getScalarList("temperatures"))
            {
                log::info("pvt hydrateEquilibrium T ", num(T), " p ",
                          num(hydrateEquilibriumPressure(T, Tref, pRef, dH)));
            }
        }

        // D13 increment 5: non-isothermal (PH) flash. Gated on `nonIsothermal`;
        // all heat-capacity coefficients are declared in-case (cp = a + b T),
        // departure enthalpy uses only the declared PR constants via da/dT.
        if (caseDict.has("nonIsothermal"))
        {
            const Dictionary& ni = caseDict.getDict("nonIsothermal");
            EnthalpyModel em;
            em.cpA = ni.getScalarList("cpA");
            em.cpB = ni.getScalarList("cpB");
            if (em.cpA.size() != c.feed.size() || em.cpB.size() != c.feed.size())
                cfFatal("pvt: nonIsothermal cpA/cpB need one value per component");
            em.referenceT = ni.getScalarOr("referenceT", em.referenceT);

            const scalar pIso = ni.getScalar("pressure");
            for (const scalar T : ni.getScalarList("enthalpyT"))
            {
                const FlashResult fl = tpFlash(c.feed, T, pIso);
                log::info("pvt enthalpyAt T ", num(T), " p ", num(pIso),
                          " h ", num(mixtureEnthalpy(c.feed, fl, T, pIso, em)));
            }
            const PhResult ph = phFlash(c.feed, pIso, ni.getScalar("hTarget"), em,
                                        ni.getScalar("Tmin"), ni.getScalar("Tmax"));
            if (!ph.converged) cfFatal("pvt: nonIsothermal phFlash failed: " + ph.message);
            log::info("pvt phFlash T ", num(ph.T), " p ", num(pIso),
                      " h ", num(ph.enthalpy), " residual ", num(ph.residual),
                      " converged ", ph.converged ? 1 : 0,
                      " iterations ", ph.iterations,
                      " beta ", num(ph.flash.vaporFraction));
        }

        // A coarse table shows interpolation and EOS agree; the regression
        // script samples it back through the same keys.
        if (caseDict.has("table"))
        {
            const Dictionary& tbl = caseDict.getDict("table");
            const std::vector<double> pGrid = tbl.getScalarList("pressures");
            const std::vector<double> tGrid = tbl.getScalarList("temperatures");
            const std::string prop = tbl.getWordOr("property", "liquidDensity");
            std::optional<BilinearTable> table;
            if (prop == "liquidDensity") table = liquidDensityTable(c.feed, pGrid, tGrid);
            else if (prop == "vapourDensity") table = vapourDensityTable(c.feed, pGrid, tGrid);
            else if (prop == "vapourFraction") table = vapourFractionTable(c.feed, pGrid, tGrid);
            else cfFatal("pvt: unknown table property '" + prop + "'");

            for (const double p : pGrid)
                for (const double t : tGrid)
                    log::info("pvt table ", prop, " T ", num(t), " p ", num(p),
                              " value ", num(table->sample(p, t)));
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error(e.what());
        return 1;
    }
}
