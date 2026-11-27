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
#include "cf/pvt/HydratePWHc.h"
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

// 5.230: the prAlpha echo needs round-trip precision -- the checker compares
// kappa/alpha/a/dadT against the published PR76/PR78 closed forms at 1e-14,
// which a 12-digit print cannot carry.
std::string num17(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

void reportOne(const FlashCase& c, double p, const Dictionary* viscosityDict,
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
    // from the phase molar mass / density ( = PR molar volume). Vc comes from
    // the DECLARED `criticalVolumes` list (m3/mol, feed order) when present;
    // absent one the historical PR-Zc fallback is kept bit-for-bit. The
    // `model` word is now loud instead of a dead field.
    if (viscosityDict)
    {
        const Dictionary& vd = *viscosityDict;
        const std::string vm = vd.getWordOr("model", "lbc");
        if (vm != "lbc")
            cfFatal("pvt: unknown viscosity model '" + vm + "' (only lbc)");
        const std::vector<Component>& comps = c.feed.components();
        std::vector<double> MWs(comps.size()), Tcs(comps.size()), Pcs(comps.size());
        std::vector<double> Vcs;
        if (vd.has("criticalVolumes"))
        {
            Vcs = vd.getScalarList("criticalVolumes");
            if (Vcs.size() != comps.size())
                cfFatal("pvt: viscosity criticalVolumes must give one value"
                        " per component (" + std::to_string(Vcs.size()) + " given,"
                        " " + std::to_string(comps.size()) + " expected)");
            for (const double vc : Vcs)
                if (!(vc > 0.0))
                    cfFatal("pvt: viscosity criticalVolumes must be positive");
        }
        else
        {
            Vcs.resize(comps.size());
            for (std::size_t i = 0; i < comps.size(); ++i)
                Vcs[i] = 0.3074 * kGasConstant * comps[i].Tc / comps[i].Pc;
        }
        for (std::size_t i = 0; i < comps.size(); ++i)
        {
            MWs[i] = comps[i].molarMass * 1.0e3;
            Tcs[i] = comps[i].Tc;
            Pcs[i] = comps[i].Pc;
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
        const Dictionary* viscosityDict = caseDict.has("viscosity")
            ? &caseDict.getDict("viscosity") : nullptr;

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

        for (const double p : c.pressures) reportOne(c, p, viscosityDict, paraPtr);

        // 5.230: the alpha-function echo, gated on a `prAlpha{}` block. It
        // prints the ACTUAL per-component kappa/alpha/a the kernel built at
        // the declared temperature (mixtureParameters, not a re-evaluation),
        // plus the prDaDT the enthalpy departure consumes, so the checker can
        // pin the PR76 / PR1978 branch (threshold omega = 0.491) against the
        // published closed forms and confirm a and dadT share one kappa.
        // Absent the block, this app's report is byte-identical to before.
        if (caseDict.has("prAlpha"))
        {
            const Dictionary& pa = caseDict.getDict("prAlpha");
            const scalar tAlpha = pa.getScalar("temperature");
            if (!(tAlpha > 0.0))
                cfFatal("pvt: prAlpha temperature must be positive");
            const PengRobinsonParams prm =
                mixtureParameters(c.feed, tAlpha, c.pressures.front());
            const auto& comps = c.feed.components();
            for (std::size_t i = 0; i < comps.size(); ++i)
            {
                log::info("pvt prAlpha comp ", comps[i].name,
                          " T ", num17(tAlpha),
                          " omega ", num17(comps[i].acentric),
                          " kappa ", num17(prm.kappaI[i]),
                          " alpha ", num17(prm.alphaI[i]),
                          " a ", num17(prm.aI[i]),
                          " b ", num17(prm.bI[i]),
                          " dadT ", num17(prDaDT(comps[i], tAlpha)));
            }
        }

        if (c.computeBubblePoint)
        {
            const double pb = bubblePressure(c.feed, c.temperature);
            log::info("pvt bubble T ", num(c.temperature), " p ", num(pb));
            reportOne(c, pb, viscosityDict, paraPtr);
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

        // E8: simplified van der Waals - Platteeuw hydrate equilibrium, gated
        // on a `hydratePWHc{}` block and mutually exclusive with the
        // Clausius-Clapeyron `hydrateEquilibrium{}` block. Every constant is
        // declared in-case (Langmuir in 1/MPa, converted to 1/Pa here); the
        // kernel holds no tabulated data.
        if (caseDict.has("hydratePWHc"))
        {
            if (caseDict.has("hydrateEquilibrium"))
                cfFatal("pvt: hydratePWHc and hydrateEquilibrium are mutually "
                        "exclusive equilibrium closures; declare one only");
            const Dictionary& hp = caseDict.getDict("hydratePWHc");

            HydratePWHcSpec spec;
            spec.structure = hp.getWord("structure");
            if (spec.structure != "sI" && spec.structure != "sII")
                cfFatal("pvt: hydratePWHc structure must be 'sI' or 'sII', got '"
                        + spec.structure + "'");

            spec.guests = hp.getWordList("guests");
            if (spec.guests.empty())
                cfFatal("pvt: hydratePWHc needs a non-empty guests ( ... ) list");
            for (std::size_t g = 0; g < spec.guests.size(); ++g)
            {
                c.feed.index(spec.guests[g]); // fatal unless declared in composition
                for (std::size_t h = 0; h < g; ++h)
                    if (spec.guests[h] == spec.guests[g])
                        cfFatal("pvt: hydratePWHc guest '" + spec.guests[g]
                                + "' is declared twice");
                // Declared as langmuir<Guest> ( ... ); one value per cage type
                // (small, large), in 1/MPa.
                const std::string key = "langmuir" + spec.guests[g];
                if (!hp.has(key))
                    cfFatal("pvt: hydratePWHc is missing the '" + key
                            + "' Langmuir constants for guest '" + spec.guests[g] + "'");
                const std::vector<double> cList = hp.getScalarList(key);
                if (cList.size() != 2)
                    cfFatal("pvt: hydratePWHc '" + key
                            + "' needs exactly two values (small cage, large cage)");
                for (const double cv : cList)
                    if (!(cv > 0.0))
                        cfFatal("pvt: hydratePWHc Langmuir constants must be "
                                "positive (1/MPa), guest '" + spec.guests[g] + "'");
                // 1/MPa -> 1/Pa: a value of C [1/MPa] equals C*1e-6 [1/Pa].
                spec.langmuir.push_back({cList[0] * 1.0e-6, cList[1] * 1.0e-6});
            }

            spec.referenceTemperature = hp.getScalar("referenceTemperature");
            spec.referencePressure = hp.getScalar("referencePressure");
            spec.enthalpyOfDissociation = hp.getScalar("enthalpyOfDissociation");
            spec.muRef = hp.getScalar("muRef");
            spec.muSlopeA = hp.getScalar("muSlopeA");
            if (!(spec.referenceTemperature > 0.0))
                cfFatal("pvt: hydratePWHc needs a positive referenceTemperature");
            if (!(spec.referencePressure > 0.0))
                cfFatal("pvt: hydratePWHc needs a positive referencePressure");
            if (!(spec.enthalpyOfDissociation > 0.0))
                cfFatal("pvt: hydratePWHc needs a positive enthalpyOfDissociation");
            if (!(spec.muRef > 0.0))
                cfFatal("pvt: hydratePWHc needs a positive muRef");
            if (!(spec.muSlopeA > 0.0))
                cfFatal("pvt: hydratePWHc needs a positive muSlopeA");

            const std::vector<double> temps = hp.getScalarList("temperatures");
            if (temps.empty())
                cfFatal("pvt: hydratePWHc needs a non-empty temperatures ( ... ) list");
            for (const double T : temps)
                if (!(T > 0.0))
                    cfFatal("pvt: hydratePWHc temperatures must be positive");

            const scalar nuSmall = hydratePWHcCagePerWater(spec.structure, 0);
            const scalar nuLarge = hydratePWHcCagePerWater(spec.structure, 1);
            log::info("pvt hydratePWHc structure ", spec.structure,
                      " guests ", spec.guests.size(),
                      " cagesSmall ", num(nuSmall), " cagesLarge ", num(nuLarge));

            bool monotonic = true;
            double previous = 0.0;
            for (const double T : temps)
            {
                const HydratePWHcPoint pt = hydratePWHcEquilibrium(c.feed, T, spec);
                if (previous > 0.0 && !(pt.pEq > previous)) monotonic = false;
                previous = pt.pEq;
                log::info("pvt hydratePWHc T ", num(T), " pEq ", num(pt.pEq),
                          " thetaSmall ", num(pt.thetaSmall),
                          " thetaLarge ", num(pt.thetaLarge));
            }
            log::info("pvt hydratePWHc curve monotonic ", monotonic ? 1 : 0);
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
