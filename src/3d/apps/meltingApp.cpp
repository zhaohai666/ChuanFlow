// chflow-3d-melting: enthalpy-porosity melting/solidification (F12).
//
// The energy equation carries the latent heat as an apparent heat capacity
// (declared through energy{ phaseChange{ latentHeat; solidusT; liquidusT;
// mushyConstant; mushyEpsilon; } }), and the momentum equation gets the
// Darcy term of the porous-medium formulation,
//
//   f_Darcy = -A_mush (1 - f_l)^2 / (f_l^3 + eps) U,
//
// recomputed from the current liquid fraction every step: the solid offers
// infinite resistance, the mushy zone a resistance that diverges as f_l -> 0,
// and the melt flows freely.
//
// Reports per `piso{ reportEvery }` steps and at the end:
//   t, front position yF (f_l = 0.5 crossing along y), total enthalpy H,
//   cumulative boundary heat input Qin, max |U| in the melt and in the solid,
//   f_l min/max.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/EnergySolver.h"
#include "cf/fvm/MeshLoader.h"
#include "case3dTurb.h"

using namespace cf;

static int appRun(int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "usage: chflow-3d-melting <case.dict>\n");
        return 1;
    }
    const Dictionary caseDict = Dictionary::parseFile(argv[1]);

    polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
    const scalar nu = app3dturb::caseNu(caseDict);

    volVectorField U(mesh, "U");
    volScalarField p(mesh, "p");
    surfaceScalarField phi(mesh, "phi", 0.0);
    app3dturb::applyPatches(caseDict, mesh, U, p);

    Schemes schemes = loadSchemes(caseDict);
    schemes.ddt = "Euler";

    if (!caseDict.has("energy")) cfFatal("melting: case needs an energy{} block");
    const EnergyProps props = EnergyProps::fromCase(caseDict);
    if (!props.phase.enabled)
        cfFatal("melting: energy{} needs a phaseChange sub-block");

    phi.correctBoundaryFlux(U);

    EnergyModel energy(mesh, phi, schemes, props, nu);
    // Thermal BCs come from the same patches{} block as the flow ones; without
    // them every patch stays adiabatic and no heat can enter the slab.
    if (caseDict.has("patches"))
    {
        const Dictionary& ps = caseDict.getDict("patches");
        volScalarField& t = energy.temperature();
        for (const std::string& name : ps.keys())
        {
            const Dictionary& pd = ps.getDict(name);
            const label patch = mesh.patchIdByName(name);
            if (pd.has("T"))
            {
                t.setPatchValue(patch, pd.getWordOr("TType", "fixedValue"), pd.getScalar("T"));
            }
            else if (pd.has("heatFlux"))
            {
                t.setPatchValue(patch, "fixedGradient", pd.getScalar("heatFlux") / props.kappa);
            }
            else if (pd.has("TType"))
            {
                t.setPatchValue(patch, pd.getWord("TType"), pd.getScalarOr("TValue", 0.0));
            }
        }
    }
    energy.initialize();

    scalar endTime = 1.0;
    int reportEvery = 10;
    // 5.196: per-cell T / f_l readout, opt-in so old frames stay identical
    const bool reportCells = caseDict.has("reportCells")
                           && caseDict.getBool("reportCells");
    PisoControls ctl = app3dturb::pisoControls(caseDict, endTime, reportEvery);

    const std::size_t nC = static_cast<std::size_t>(mesh.nCells());
    Field<Vector> darcy(nC, Vector::zero());
    ctl.cellForce = &darcy;

    // cell order by y for the front tracker
    std::vector<std::size_t> order(nC);
    for (std::size_t c = 0; c < nC; ++c) order[c] = c;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
              { return mesh.cellCentres()[a].y < mesh.cellCentres()[b].y; });

    const scalar A = props.phase.mushyConstant;
    const scalar eps = props.phase.mushyEpsilon;

    scalar t = 0;
    scalar qCum = 0;
    scalar qPrev = 0;
    const scalar h0 = energy.totalEnthalpy();

    auto patchHeat = [&]()
    {
        scalar q = 0;
        for (label pt = 0; pt < mesh.nPatches(); ++pt) q += energy.patchHeatRate(pt);
        return q;
    };

    auto report = [&]()
    {
        const auto& fl = energy.liquidFraction();
        // front: f_l crosses 0.5 upward from the bottom of the column
        scalar yF = mesh.minX().y;
        for (std::size_t i = 1; i < nC; ++i)
        {
            const std::size_t a = order[i - 1], b = order[i];
            if (fl[a] < 0.5 && fl[b] >= 0.5)
            {
                const scalar ya = mesh.cellCentres()[a].y, yb = mesh.cellCentres()[b].y;
                yF = ya + 0.5 * (yb - ya); // fl is linear in the mushy band
                break;
            }
        }
        scalar uLiquid = 0, uSolid = 0, flMin = 1, flMax = 0;
        for (std::size_t c = 0; c < nC; ++c)
        {
            flMin = std::min(flMin, fl[c]);
            flMax = std::max(flMax, fl[c]);
            const scalar um = mag(U.internalField()[c]);
            if (fl[c] > 0.9) uLiquid = std::max(uLiquid, um);
            if (fl[c] < 0.1) uSolid = std::max(uSolid, um);
        }
        // 5.196: opt-in per-cell readout (`reportCells true`), off by default
        // so the pre-5.196 frames are byte-identical. It exists because the
        // latent-heat bookkeeping is a PER-CELL enthalpy statement, which the
        // global H line cannot witness.
        if (reportCells)
        {
            for (std::size_t c = 0; c < nC; ++c)
                std::printf("cell %zu T %s fl %s\n", c,
                            app3dturb::num(energy.temperature().internalField()[c]).c_str(),
                            app3dturb::num(fl[c]).c_str());
        }
        std::printf("t %s yF %s H %s Qin %s uLiquid %s uSolid %s flMin %s flMax %s\n",
                    app3dturb::num(t).c_str(), app3dturb::num(yF).c_str(),
                    app3dturb::num(energy.totalEnthalpy() - h0).c_str(),
                    app3dturb::num(qCum).c_str(), app3dturb::num(uLiquid).c_str(),
                    app3dturb::num(uSolid).c_str(), app3dturb::num(flMin).c_str(),
                    app3dturb::num(flMax).c_str());
    };

    qPrev = patchHeat();
    report();
    int step = 0;
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const scalar dt = std::min(ctl.dt, endTime - t);
        // Darcy momentum damping from the current liquid fraction. The term
        // is explicit in U, so the damping rate is capped at half the step
        // rate: the solid still loses half its velocity per step and reaches
        // rest in a handful of them, and the scheme cannot oscillate.
        const auto& fl = energy.liquidFraction();
        const scalar cap = 0.5 / dt;
        for (std::size_t c = 0; c < nC; ++c)
        {
            const scalar damp = std::min(A * (1.0 - fl[c]) * (1.0 - fl[c])
                                         / (fl[c] * fl[c] * fl[c] + eps), cap);
            darcy[c] = U.internalField()[c] * (-damp);
        }
        PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        solveTransientIncompressiblePiso(mesh, U, p, phi, nu, schemes, stepCtl,
                                         nullptr, &energy);
        t += dt;
        ++step;

        const scalar qNow = patchHeat();
        qCum += 0.5 * dt * (qPrev + qNow);
        qPrev = qNow;
        if (step % reportEvery == 0) report();
    }
    report();

    if (caseDict.has("output"))
    {
        const Dictionary& od = caseDict.getDict("output");
        const std::string vtk = od.getWordOr("vtk", "");
        if (!vtk.empty())
        {
            VtkWriter w(mesh);
            w.addCellVector("U", U.internalField());
            w.addCellScalar("p", p.internalField());
            w.addCellScalar("T", energy.temperature().internalField());
            w.addCellScalar("liquidFraction", energy.liquidFraction());
            w.write(vtk);
        }
    }
    return 0;
}

// Uncaught cfFatal dies silently on MSVC (fail-fast, nothing on the pipe the
// checkers read); catch and report, the adaptiveRefineApp convention.
int main(int argc, char** argv)
{
    try
    {
        return appRun(argc, argv);
    }
    catch (const std::exception& e)
    {
        log::error("cf::Error: ", e.what());
        return 1;
    }
}
