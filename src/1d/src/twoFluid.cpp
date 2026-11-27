#include "cf/1d/TwoFluid.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "cf/common/Error.h"

#include "cf/1d/Network.h" // churchillFriction
#include "cf/pvt/WaxSolubility.h" // van't Hoff wax saturation (5.96)
#include "cf/pvt/Components.h" // G1: named components for the feed (5.143)
#include "cf/pvt/Flash.h"      // G1: per-cell tpFlash diagnostics (5.143)
#include "cf/pvt/HydrateEquilibrium.h" // 5.186: the declared CC hydrate line
#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double mixtureMu(const TwoFluidFluid& f, double a)
{
    return a * f.muL + (1.0 - a) * f.muG;
}

/// Theis well function E1(x) = int_x^inf exp(-u)/u du. Two independent
/// branches: the convergent series for x < 1 (Euler's gamma + logarithm),
/// and the Legendre continued fraction for x >= 1. x -> inf is the s -> 0
/// limit of the line-source kernel and is exactly zero: that is what makes
/// the Duhamel superposition explicit in the flux history.
double expIntegralE1(double x)
{
    if (!(x > 0.0) || !std::isfinite(x)) return 0.0;
    if (x < 1.0)
    {
        constexpr double gamma = 0.577215664901532860606512090082402431;
        double sum = -gamma - std::log(x);
        double term = 1.0;
        for (int n = 1; n <= 80; ++n)
        {
            term *= -x / double(n);          // (-x)^n / n!
            const double add = -term / double(n);
            sum += add;
            if (std::abs(add) < 1e-18 * std::abs(sum)) break;
        }
        return sum;
    }
    // NR expint(n=1,x): b = x + 1; h = 1/b; a_i = -i^2; b_i = x + 2i + 1
    constexpr double tiny = 1e-300;
    constexpr double eps = 1e-16;
    double b = x + 1.0;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= 300; ++i)
    {
        const double a = -double(i) * double(i);
        b += 2.0;
        d = 1.0 / (a * d + b);
        if (std::abs(d) < tiny) d = tiny;
        c = b + a / c;
        if (std::abs(c) < tiny) c = tiny;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < eps) break;
    }
    return std::exp(-x) * h;
}

/// Duhamel superposition of a piecewise-constant far-field flux history.
/// q'(t) = q'_m on (t_{m-1}, t_m]; the step response is
/// F(s) = E1(r_e^2/(4 alpha s))/(4 pi k), F(0) = 0, so
///   dT(t_n) = sum_m (q'_m - q'_{m-1}) F(t_n - t_m),  q'_0 = 0.
double lineSourceRise(double tNow, const std::vector<double>& qHist,
                      const std::vector<double>& tHist, double rOuter,
                      double alpha, double k)
{
    double rise = 0.0;
    double qPrev = 0.0;
    for (std::size_t m = 0; m < qHist.size(); ++m)
    {
        const double s = tNow - tHist[m];
        // s == 0 is the current step: F(0) = 0, so its own flux cannot feed
        // back into its own boundary. No implicit loop is needed.
        const double F = (s > 0.0)
            ? expIntegralE1(rOuter * rOuter / (4.0 * alpha * s))
              / (4.0 * M_PI * k)
            : 0.0;
        rise += (qHist[m] - qPrev) * F;
        qPrev = qHist[m];
    }
    return rise;
}

/// Godunov flux of the quadratic F(a) = qa a^2 + qb a + qc between the left
/// and right states. Convex F (qa > 0) shocks on a falling jump and
/// rarefies on a rising one, concave F mirrors both, and the linear case
/// (qa = 0) is plain upwind on the constant slope. The rarefaction branch
/// evaluates F at the sonic state -qb/(2 qa) clamped into the jump.
double godunovQuadratic(double qa, double qb, double qc, double aL, double aR)
{
    const double fL = qc + aL * (qb + qa * aL);
    if (aL == aR) return fL;
    const double fR = qc + aR * (qb + qa * aR);
    const bool shock = qa > 0 ? aL > aR : aL < aR;
    if (shock)
    {
        const double s = (fL - fR) / (aL - aR);
        return s >= 0 ? fL : fR;
    }
    if (qa == 0) return qb >= 0 ? fL : fR;
    const double sonic = -qb / (2.0 * qa);
    const double lo = std::min(aL, aR), hi = std::max(aL, aR);
    const double a = std::min(std::max(sonic, lo), hi);
    return qc + a * (qb + qa * a);
}

} // namespace

double exponentialIntegralE1(double x)
{
    return expIntegralE1(x);
}

TwoFluidPipe::TwoFluidPipe(const TwoFluidCase& c)
    : cfg_(c)
{
    const TwoFluidFluid& f = cfg_.fluid;
    if (!(f.rhoL > 0) || !(f.rhoG > 0) || !(f.muL > 0) || !(f.muG > 0))
        cfFatal("twoFluid: fluid needs positive densities and viscosities");
    if (f.gasR < 0)
        cfFatal("twoFluid: negative gasR");
    if (cfg_.pipes.empty())
        cfFatal("twoFluid: case has no pipe segment");
    diameter_ = cfg_.pipes.front().diameter;
    length_ = 0;
    for (const TwoFluidGeo& g : cfg_.pipes)
    {
        if (!(g.length > 0) || !(g.diameter > 0) || g.sections < 1)
            cfFatal("twoFluid: pipe needs positive length/diameter and >= 1"
                    " sections");
        if (g.roughness < 0)
            cfFatal("twoFluid: negative roughness");
        // the mixture momentum is lumped on one velocity, which needs a
        // constant bore along the chain
        if (g.diameter != diameter_)
            cfFatal("twoFluid: a serial chain needs a uniform diameter"
                    " (the lumped mixture momentum carries one velocity)");
        length_ += g.length;
    }

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
    if (!cfg_.networkMode)
    {
        // serial-case boundary combinations; network members are checked
        // by the network loader instead (the junction end boundary is
        // owned by the TwoFluidNetwork and rewritten every step)
        const bool inSealed = in.type == TwoFluidBoundary::Type::sealed;
        const bool outSealed = out.type == TwoFluidBoundary::Type::sealed;
        if (inSealed && !outSealed)
            cfFatal("twoFluid: a sealed inlet needs a sealed outlet"
                    " (incompressible phases cannot draw through-flow)");
        if (outSealed && in.type == TwoFluidBoundary::Type::fixedPressure)
            cfFatal("twoFluid: fixed-pressure inlet into a sealed pipe has no"
                    " steady state with incompressible phases");
        if (in.type == TwoFluidBoundary::Type::phaseRates
            && in.vsl + in.vsg < 0)
            cfFatal("twoFluid: phaseRates injection with negative total rate"
                    " is not supported");
    }
    if (slipClosure() && !(cfg_.closure.interfacialDrag > 0))
        cfFatal("twoFluid: the dynamic/full closures need a positive"
                " interfacialDrag");

    n_ = 0;
    for (const TwoFluidGeo& g : cfg_.pipes) n_ += g.sections;
    dxs_.reserve(n_);
    sinT_.reserve(n_);
    relrough_.reserve(n_);
    for (const TwoFluidGeo& g : cfg_.pipes)
    {
        const double dx = g.length / g.sections;
        const double st = std::sin(g.angleDeg * M_PI / 180.0);
        for (int k = 0; k < g.sections; ++k)
        {
            dxs_.push_back(dx);
            sinT_.push_back(st);
            relrough_.push_back(g.roughness / g.diameter);
        }
    }
    alpha_.assign(n_, cfg_.alpha0);
    water_.assign(n_, cfg_.water0);
    slip_.assign(n_, 0.0);
    sDotOld_.assign(n_, 0.0);
    std::vector<double> centers(n_);
    double arc = 0;
    for (int i = 0; i < n_; ++i)
    {
        centers[i] = arc + 0.5 * dxs_[i];
        arc += dxs_[i];
    }
    for (const TwoFluidFill& fill : cfg_.fills)
    {
        if (fill.from > fill.to || fill.alpha < 0 || fill.alpha > 1)
            cfFatal("twoFluid: fill block needs from <= to and alpha in"
                    " [0,1]");
        for (int i = 0; i < n_; ++i)
        {
            if (centers[i] >= fill.from && centers[i] <= fill.to)
                alpha_[i] = fill.alpha;
        }
    }
    // incompressible phases: an imposed rate is felt along the whole pipe
    // from t = 0; pressure-driven and sealed pipes start from rest unless
    // 5.77 gave the case an explicit initial velocity (network members own
    // their velocity through the junction Newton, so vm0 is serial-only)
    if (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates)
        vM_ = vMOld_ = cfg_.inlet.vsl + cfg_.inlet.vsg;
    else if (!cfg_.networkMode && cfg_.vm0 != 0.0)
        vM_ = vMOld_ = cfg_.vm0;

    // vertical drop of each cell centre from the chain inlet (geothermal
    // depth), accumulated along the per-segment slopes
    centerDrop_.assign(n_, 0.0);
    double drop = 0;
    for (int i = 0; i < n_; ++i)
    {
        drop += 0.5 * dxs_[i] * sinT_[i];
        centerDrop_[i] = drop;
        drop += 0.5 * dxs_[i] * sinT_[i];
    }

    // gas inventory per unit area: rhoG (1 - alpha) dx; with the EOS on it
    // is the transported state, otherwise it stays at its initial value
    // and every density lookup falls back to the constant fluid rhoG
    gasMass_.resize(n_);
    for (int i = 0; i < n_; ++i)
        gasMass_[i] = f.rhoG * (1.0 - alpha_[i]) * dxs_[i];

    // O3 (5.162): the purge wave window must be a real window on this
    // chain's cell layout (the loader checked the values, only the arc
    // geometry is known here)
    if (cfg_.purgeWaveOn)
    {
        const TwoFluidPurgeWave& pw = cfg_.purgeWave;
        const double end = pw.riserEnd > 0 ? pw.riserEnd : length_;
        if (pw.riserStart >= end || end > length_ + 1e-9)
            cfFatal("twoFluid: purgeWave riser window outside the chain");
        bool any = false;
        double acc = 0;
        for (int i = 0; i < n_; ++i)
        {
            const double centre = acc + 0.5 * dxs_[i];
            if (centre >= pw.riserStart && centre <= end)
            { any = true; break; }
            acc += dxs_[i];
        }
        if (!any)
            cfFatal("twoFluid: purgeWave riser window holds no cell centre");
    }

    // E6 (5.176): the pig's launch/arrival arcs must sit on this chain. The
    // loader checked their signs and ordering; only the arc geometry is
    // known here. Launch the pig at construction so the t = 0 report already
    // shows it RUN at startX with a zero pile.
    if (cfg_.piggingOn)
    {
        const TwoFluidPigging& pg = cfg_.pigging;
        const double end = pg.endX > 0 ? pg.endX : length_;
        if (pg.startX < 0 || pg.startX > length_ + 1e-9 || end > length_ + 1e-9
            || pg.startX >= end)
            cfFatal("twoFluid: pigging launch/arrival arc outside the chain");
        pigState_.launched = true;
        pigState_.x = pg.startX;
        pigState_.v = pg.pigVelocity;
    }

    // E7c (5.183): the acoustic field starts at the uniform declared
    // overpressure with the line at rest (u' = 0 on every face); the vent
    // opens at t = 0. The state is allocated only on the declared path, so
    // the undeclared constructor keeps the old footprint verbatim.
    if (cfg_.acousticsOn)
    {
        pA_.assign(n_, cfg_.acoustics.initialGauge);
        uA_.assign(n_ + 1, 0.0);
    }

    if (cfg_.thermalOn)
    {
        const TwoFluidThermal& th = cfg_.thermal;
        if (!(f.cL > 0) || !(f.cG > 0))
            cfFatal("twoFluid: thermal on needs positive fluid cL/cG");
        if (th.wallModel != 0 && th.wallModel != 1)
            cfFatal("twoFluid: thermal wallModel is constant or geothermal");
        if (th.wallModel == 1)
        {
            if (!(th.geothermalGradient > 0))
                cfFatal("twoFluid: geothermal wall needs a positive gradient");
            double dropTotal = 0;
            for (int i = 0; i < n_; ++i) dropTotal += dxs_[i] * sinT_[i];
            const double dExit = th.inletDepth - dropTotal;
            if (dExit < 0)
                cfFatal("twoFluid: geothermal pipe leaves the ground "
                        "(depth at the outlet is " + std::to_string(dExit)
                        + " m)");
        }
        // both ghost temperatures are required: the outlet one only feeds
        // backflow, but backflow is a normal fixed-pressure outcome
        if (cfg_.inlet.type != TwoFluidBoundary::Type::sealed
            && !(cfg_.inlet.temperature > 0))
            cfFatal("twoFluid: thermal on needs an inlet boundary temperature");
        if (cfg_.outlet.type != TwoFluidBoundary::Type::sealed
            && !(cfg_.outlet.temperature > 0))
            cfFatal("twoFluid: thermal on needs an outlet boundary temperature");
        T_.assign(n_, cfg_.temperature0);
    }

    if (cfg_.formationOn)
    {
        if (!cfg_.thermalOn)
            cfFatal("twoFluid: formation{} needs a thermal{} block");
        const TwoFluidFormation& fm = cfg_.formation;
        const double rw = 0.5 * diameter_;
        if (!(fm.conductivity > 0))
            cfFatal("twoFluid: formation needs a positive conductivity");
        if (!(fm.outerRadius > rw))
            cfFatal("twoFluid: formation outerRadius must exceed the bore"
                    " radius (diameter/2)");
        if (!(fm.rhoC > 0))
            cfFatal("twoFluid: formation needs a positive rhoC");
        if (fm.shells < 1)
            cfFatal("twoFluid: formation shells must be at least one");
        const int M = fm.shells;
        // geometric node radii from the bore wall out to the far field
        formR_.resize(M + 1);
        const double ratio = fm.outerRadius / rw;
        for (int j = 0; j <= M; ++j)
            formR_[j] = rw * std::pow(ratio, double(j) / M);
        formR_[0] = rw;
        formR_[M] = fm.outerRadius;
        formG_.assign(n_, std::vector<double>(M, 0.0));
        formC_.assign(n_, std::vector<double>(M, 0.0));
        formT_.assign(n_, std::vector<double>(M + 1, 0.0));
        const double twoPiK = 2.0 * M_PI * fm.conductivity;
        for (int i = 0; i < n_; ++i)
        {
            const double dx = dxs_[i];
            for (int j = 0; j < M; ++j)
            {
                // cylindrical Fourier conductance between adjacent nodes
                formG_[i][j] = twoPiK * dx
                             / std::log(formR_[j + 1] / formR_[j]);
                // node j owns the shell between geometric midpoints; the
                // last node extends out to the far-field radius r_e
                const double lo = (j == 0)
                    ? rw : std::sqrt(formR_[j - 1] * formR_[j]);
                const double hi = (j == M - 1)
                    ? fm.outerRadius : std::sqrt(formR_[j] * formR_[j + 1]);
                formC_[i][j] = fm.rhoC * dx * M_PI * (hi * hi - lo * lo);
            }
            // the rock starts at the far-field temperature
            const double te = wallTemperature(i);
            for (int j = 0; j <= M; ++j) formT_[i][j] = te;
        }
    }

    if (cfg_.hydrateOn)
    {
        if (!eosOn())
            cfFatal("twoFluid: hydrate kinetics needs the gas EOS"
                    " (fluid gasR > 0)");
        if (!cfg_.thermalOn)
            cfFatal("twoFluid: hydrate kinetics needs a thermal{} block");
        const TwoFluidHydrate& hy = cfg_.hydrate;
        if (!(hy.growthRate >= 0) || !(hy.dissociationRate >= 0))
            cfFatal("twoFluid: negative hydrate rate");
        if (!(hy.gasContent > 0))
            cfFatal("twoFluid: hydrate needs a positive gasContent");
        if (!(hy.latentHeat > 0))
            cfFatal("twoFluid: hydrate needs a positive latentHeat");
        if (hy.initial < 0 || hy.initial > 1)
            cfFatal("twoFluid: initial hydrate fraction outside [0,1]");
        if (cfg_.water0 + hy.initial > 1.0)
            cfFatal("twoFluid: water cut plus initial hydrate fraction"
                    " exceeds 1");
        hydrate_.assign(n_, hy.initial);
        hydrateDelta_.assign(n_, 0.0);
    }

    if (cfg_.waxOn)
    {
        if (!(cfg_.thermalOn))
            cfFatal("twoFluid: wax precipitation needs a thermal{} block");
        const TwoFluidWax& wx = cfg_.wax;
        if (!(wx.meltingTemperature > 0) || !(wx.fusionEnthalpy > 0))
            cfFatal("twoFluid: wax needs positive meltingTemperature and"
                    " fusionEnthalpy");
        if (!(wx.waxFraction >= 0) || wx.waxFraction > 1)
            cfFatal("twoFluid: waxFraction outside [0,1]");
        if (!(wx.growthRate >= 0) || !(wx.dissolutionRate >= 0))
            cfFatal("twoFluid: negative wax rate");
        if (wx.initial < 0 || wx.initial > 1)
            cfFatal("twoFluid: initial wax fraction outside [0,1]");
        if (wx.latentHeat < 0)
            cfFatal("twoFluid: negative wax latentHeat");
        if (wx.depositRoughness < 0)
            cfFatal("twoFluid: negative wax depositRoughness");
        wax_.assign(n_, wx.initial);
        waxDelta_.assign(n_, 0.0);
    }

    if (cfg_.asphalteneOn)
    {
        // 5.99: the AOP and its ramp are declared experimental/model inputs;
        // the kernel only implements the onset equilibrium and the kinetics.
        const TwoFluidAsphaltene& az = cfg_.asphaltene;
        if (!(az.onsetPressure > 0))
            cfFatal("twoFluid: asphaltene needs a positive onsetPressure");
        if (az.onsetWidth < 0)
            cfFatal("twoFluid: negative asphaltene onsetWidth");
        if (!(az.asphalteneFraction >= 0) || az.asphalteneFraction > 1)
            cfFatal("twoFluid: asphalteneFraction outside [0,1]");
        if (!(az.growthRate >= 0) || !(az.dissolutionRate >= 0))
            cfFatal("twoFluid: negative asphaltene rate");
        if (az.initial < 0 || az.initial > 1)
            cfFatal("twoFluid: initial asphaltene fraction outside [0,1]");
        asph_.assign(n_, az.initial);
    }

    if (cfg_.scaleOn)
    {
        // 5.99: activities and Ksp are declared; the default Ksp cites
        // calcite at 25 C. Only the drive Omega = IAP/Ksp vs 1 is physics.
        const TwoFluidScale& sc = cfg_.scale;
        if (sc.calciumActivity < 0 || sc.carbonateActivity < 0)
            cfFatal("twoFluid: negative scale ion activity");
        if (!(sc.solubilityProduct > 0))
            cfFatal("twoFluid: scale needs a positive solubilityProduct");
        if (!(sc.scaleFraction >= 0) || sc.scaleFraction > 1)
            cfFatal("twoFluid: scaleFraction outside [0,1]");
        if (!(sc.growthRate >= 0) || !(sc.dissolutionRate >= 0))
            cfFatal("twoFluid: negative scale rate");
        if (sc.initial < 0 || sc.initial > 1)
            cfFatal("twoFluid: initial scale fraction outside [0,1]");
        scale_.assign(n_, sc.initial);
    }

    if (cfg_.sandOn)
    {
        if (!(cfg_.sand.terminalVelocity > 0))
            cfFatal("twoFluid: sand needs a positive terminalVelocity");
        if (!(cfg_.sand.criticalVelocity > 0))
            cfFatal("twoFluid: sand needs a positive criticalVelocity");
        if (cfg_.sand.erodeRate < 0)
            cfFatal("twoFluid: negative sand erodeRate");
        if (cfg_.sand.initial < 0 || cfg_.sand.initial > 1)
            cfFatal("twoFluid: initial suspended sand fraction outside [0,1]");
        if (cfg_.sand.bedInitial < 0 || cfg_.sand.bedInitial >= 1)
            cfFatal("twoFluid: initial sand bed fraction outside [0,1)");
        sandFrac_.assign(n_, cfg_.sand.initial);
        bed_.assign(n_, cfg_.sand.bedInitial);
    }

    if (cfg_.compositionOn)
    {
        // G1: the loader validated the feed (names known, 2..8 components,
        // z >= 0, sum within 1e-8 of 1, auto-normalised). The boundary
        // inflow composition defaults to the declared feed; a declared
        // zInlet is checked against the same rules.
        const int nc = static_cast<int>(cfg_.compositionZ0.size());
        if (nc < 2 || nc > 8)
            cfFatal("twoFluid: composition needs 2..8 components");
        for (const std::string& name : cfg_.compositionNames)
            pvt::requireComponent(name); // unknown names are loud
        auto prepare = [this, nc](TwoFluidBoundary& b, const char* who) {
            if (b.zInlet.empty())
            {
                b.zInlet = cfg_.compositionZ0;
                return;
            }
            if (static_cast<int>(b.zInlet.size()) != nc)
                cfFatal(std::string("twoFluid: ") + who + " zInlet needs "
                        + std::to_string(nc) + " fractions");
            double sum = 0;
            for (double z : b.zInlet)
            {
                if (!(z >= 0))
                    cfFatal(std::string("twoFluid: ") + who
                            + " zInlet has a negative fraction");
                sum += z;
            }
            if (std::fabs(sum - 1.0) > 1e-8)
                cfFatal(std::string("twoFluid: ") + who + " zInlet fractions"
                        " must sum to 1 within 1e-8");
            for (double& z : b.zInlet) z /= sum;
        };
        prepare(cfg_.inlet, "inlet");
        prepare(cfg_.outlet, "outlet");
        // flash pressure reference where the cell EOS pressure is not
        // available (gasR off, or a liquid-full cell): outlet first (the
        // gauge convention's base), then inlet
        flashPRef_ = 0.0;
        if (cfg_.outlet.type == TwoFluidBoundary::Type::fixedPressure
            && cfg_.outlet.pressure > 0)
            flashPRef_ = cfg_.outlet.pressure;
        else if (cfg_.inlet.type == TwoFluidBoundary::Type::fixedPressure
                 && cfg_.inlet.pressure > 0)
            flashPRef_ = cfg_.inlet.pressure;
        if (!eosOn() && flashPRef_ <= 0)
            cfFatal("twoFluid: composition flash needs a pressure: enable"
                    " the gas EOS (gasR) or declare a fixedPressure"
                    " boundary");
        compZ_.assign(nc - 1, std::vector<double>(n_, 0.0));
        for (int k = 0; k < nc - 1; ++k)
            std::fill(compZ_[k].begin(), compZ_[k].end(),
                      cfg_.compositionZ0[k]);
        compBeta_.assign(n_, 0.0);
        compRhoL_.assign(n_, 0.0);
        compRhoG_.assign(n_, 0.0);
        compositionStep(0.0); // the t=0 frame already reports the flash
    }

    // 5.98/5.99: start the effective diameter from the TOTAL declared
    // deposit. With feedback off dEff_ is never read by the wall law, so the
    // old expressions stay bit-for-bit; with feedback on the deposit steps
    // refresh it every step via refreshDepositDiameter(). The sand bed must
    // be sized BEFORE this call: a feedback bore starts from the initial
    // bed, not from the raw diameter.
    dEff_.assign(n_, diameter_);
    // 5.219/5.220: the blockage layer and the wear start at zero; the
    // deposit refresh below already composes them (both steps are inert
    // when the blocks are undeclared).
    blockage_.assign(n_, 0.0);
    wear_.assign(n_, 0.0);
    refreshDepositDiameter(0.0);
    // 5.222: the first frame reports the schedule's wall temperature at
    // t = 0 (the expression is time-free when the schedule is undeclared).
    if (cfg_.ambientScheduleOn)
    {
        wallTUsed_.assign(n_, 0.0);
        for (int i = 0; i < n_; ++i) wallTUsed_[i] = wallTemperatureAt(i, 0.0);
    }

    if (cfg_.entrainmentOn)
    {
        if (!(cfg_.entrainment.vRef > 0))
            cfFatal("twoFluid: entrainment needs a positive vRef ramp"
                    " width");
        if (!(cfg_.entrainment.relaxTime >= 0))
            cfFatal("twoFluid: negative entrainment relaxTime");
        if (cfg_.entrainment.initial < 0 || cfg_.entrainment.initial > 1)
            cfFatal("twoFluid: initial entrained fraction outside [0,1]");
        entrained_.assign(n_, cfg_.entrainment.initial);
    }

    if (cfg_.emulsionOn)
    {
        // the law prices the Newtonian liquid share; a power-law liquid
        // carries its own wall law and the declaration would be ignored
        if (cfg_.fluid.flowIndex > 0)
            cfFatal("twoFluid: emulsion{} declares a mixing viscosity for"
                    " the Newtonian liquid share, but the fluid declares a"
                    " power-law flowIndex: the two laws would fight over"
                    " the same wall");
        if (!(cfg_.emulsion.muWater > 0))
            cfFatal("twoFluid: emulsion needs a positive muWater");
    }

    if (cfg_.inhibitorOn)
    {
        if (!cfg_.hydrateOn)
            cfFatal("twoFluid: inhibitor{} depresses the hydrate"
                    " equilibrium line but no hydrate{} is declared: the"
                    " offset would act on nothing");
        if (!(cfg_.inhibitor.kDepression > 0))
            cfFatal("twoFluid: inhibitor needs a positive kDepression"
                    " (K per unit mass fraction)");
        if (cfg_.inhibitor.initial < 0 || cfg_.inhibitor.initial > 1)
            cfFatal("twoFluid: initial inhibitor fraction outside [0,1]");
        inhibitor_.assign(n_, cfg_.inhibitor.initial);
    }

    // O4 (5.163): the rigorous activity root rides the transported
    // inhibitor fraction, so both parent declarations must be present
    if (cfg_.inhibitorEquilibriumOn)
    {
        if (!cfg_.inhibitorOn)
            cfFatal("twoFluid: inhibitorEquilibrium{} roots the curve at"
                    " the transported inhibitor fraction but no inhibitor{}"
                    " is declared: the activity would act on nothing");
        if (!cfg_.hydrateOn)
            cfFatal("twoFluid: inhibitorEquilibrium{} shifts the hydrate"
                    " equilibrium curve but no hydrate{} is declared");
        if (!(cfg_.inhibitorEquilibrium.activityK > 0))
            cfFatal("twoFluid: inhibitorEquilibrium needs a positive"
                    " activityK (the activity constant in"
                    " aW = 1/(1 + activityK w))");
    }
}

void TwoFluidPipe::fluxCoeffs(double& qa, double& qb, double& qc) const
{
    // F(a) = a_l * vL(a): liquid volume flux per unit area (algebraic
    // closures; the dynamic closure builds face-local coefficients from
    // the slip state)
    if (cfg_.closure.model == 0) // wallis: vG = vM + alpha*us
    {
        qa = cfg_.closure.slipVelocity;
        qb = vM_ - cfg_.closure.slipVelocity;
        qc = 0.0;
    }
    else // zuber: vG = C0 vM + vD
    {
        qa = 0.0;
        qb = cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
        qc = vM_ * (1.0 - cfg_.closure.c0) - cfg_.closure.driftVelocity;
    }
}

double TwoFluidPipe::liquidFlux(double alpha, double s, int cell) const
{
    if (slipClosure()) // dynamic/full: vL = vM - alpha_g s
        return s * alpha * alpha + (vM_ - s) * alpha;
    if (cell >= 0 && cfg_.closure.gravitySlip)
    {
        // gravity-projected drift at the section (shut-in-safe variant)
        const double d = driftOf(cell);
        if (cfg_.closure.model == 0)
            return alpha * (d * alpha + vM_ - d);
        return alpha * (cfg_.closure.c0 * vM_ + d)
             + vM_ * (1.0 - cfg_.closure.c0) - d;
    }
    double qa, qb, qc;
    fluxCoeffs(qa, qb, qc);
    return qc + alpha * (qb + qa * alpha);
}

double TwoFluidPipe::gasVelocity(double alpha, double s, int cell) const
{
    if (cfg_.closure.model == 0)
    {
        if (cell >= 0 && cfg_.closure.gravitySlip)
            return vM_ + alpha * driftOf(cell);
        return vM_ + alpha * cfg_.closure.slipVelocity;
    }
    if (slipClosure()) // dynamic/full: vG = vM + alpha s
        return vM_ + alpha * s;
    if (cell >= 0 && cfg_.closure.gravitySlip)
        return cfg_.closure.c0 * vM_ + driftOf(cell);
    return cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
}

void TwoFluidPipe::gasFluxCoeffs(double& qa, double& qb, double& qc) const
{
    // G(a) = alpha_g * vG(a): gas volume flux per unit area (algebraic
    // closures; the dynamic closure builds face-local coefficients from
    // the slip state)
    if (cfg_.closure.model == 0) // wallis: vG = vM + alpha*us
    {
        qa = -cfg_.closure.slipVelocity;
        qb = cfg_.closure.slipVelocity - vM_;
        qc = vM_;
    }
    else // zuber: vG = C0 vM + vD (constant)
    {
        qa = 0.0;
        qb = -(cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity);
        qc = cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
    }
}

double TwoFluidPipe::gasVolumeFlux(double alpha, double s, int cell) const
{
    if (slipClosure()) // dynamic/full: G = (1 - alpha)(vM + alpha s)
        return vM_ + (s - vM_) * alpha - s * alpha * alpha;
    if (cell >= 0 && cfg_.closure.gravitySlip)
    {
        // gravity-projected drift at the section (shut-in-safe variant)
        const double d = driftOf(cell);
        if (cfg_.closure.model == 0)
            return vM_ + (d - vM_) * alpha - d * alpha * alpha;
        const double vg = cfg_.closure.c0 * vM_ + d;
        return (1.0 - alpha) * vg;
    }
    double qa, qb, qc;
    gasFluxCoeffs(qa, qb, qc);
    return qc + alpha * (qb + qa * alpha);
}

double TwoFluidPipe::mixtureRhoCell(int i) const
{
    double rhoG = cfg_.fluid.rhoG;
    if (eosOn())
    {
        const double ag = 1.0 - alpha_[i];
        rhoG = ag > 1e-12 ? gasMass_[i] / (ag * dxs_[i]) : cfg_.fluid.rhoG;
    }
    return alpha_[i] * cfg_.fluid.rhoL + (1.0 - alpha_[i]) * rhoG;
}

// 5.187 (F-class): the transient risk report, the transient P22 readout
// priced off this frame's own sections. PURE REPORT: nothing below writes any
// model field; every number is recomputed from the state the last step left.
TwoFluidRiskReadout TwoFluidPipe::riskReadout() const
{
    TwoFluidRiskReadout r;
    if (!cfg_.riskOn) return r;
    const TwoFluidRiskConfig& rk = cfg_.risk;
    r.declared = true;

    // ---- erosion: API RP 14E at the worst vm / Ve section ----------------
    // Ve = C / sqrt(rho_m) with C quoted in ft/s at lb/ft3; the SI density
    // converts with 1 lb/ft3 = 16.0184634 kg/m3 and the result with 0.3048,
    // the exact conversion 5.159 already uses on the steady march.
    double best = -1.0;
    int over = 0;
    for (int i = 0; i < n_; ++i)
    {
        const double rhoM = mixtureRhoCell(i);
        const double rhoLb = rhoM / 16.0184634;
        const double velim = rhoLb > 0.0
                                 ? 0.3048 * rk.erosionC / std::sqrt(rhoLb)
                                 : 0.0;
        const double ratio = velim > 0.0 ? std::fabs(vM_) / velim : 0.0;
        if (ratio > 1.0) ++over;
        if (ratio > best)
        {
            best = ratio;
            r.erosionSeg = i;
            r.erosionRhoM = rhoM;
            r.erosionVm = std::fabs(vM_);
            r.erosionVelim = velim;
            r.erosionRatio = ratio;
        }
    }
    r.erosionOverCount = over;

    // ---- hydrate: the NON-BLOCKING reading -------------------------------
    // The declared dissociation line priced at each section's own EOS
    // pressure and temperature; the worst is the LARGEST subcooling Teq - T
    // (formation needs the stream below Teq and the margin grows with it).
    // h itself never enters: the field is reported, not charged.
    if (rk.hydrateOn)
    {
        r.hydrateOn = true;
        const std::vector<double> p = cellPressure();
        double worst = 0.0;
        bool any = false;
        for (int i = 0; i < n_; ++i)
        {
            if (!(p[i] > 0.0)) continue; // sealed section: no pressure to price
            const double teq = pvt::hydrateEquilibriumTemperature(
                p[i], rk.hydrateTref, rk.hydratePRef, rk.hydrateDH);
            const double t = sectionTemperature(i);
            const double dts = teq - t;
            if (!any || dts > worst)
            {
                any = true;
                worst = dts;
                r.hydrateSeg = i;
                r.hydrateP = p[i];
                r.hydrateT = t;
                r.hydrateTeq = teq;
                r.hydrateSubcool = dts;
            }
        }
        r.hydrateGrade = r.hydrateSubcool < 0.0
                             ? 0
                             : (r.hydrateSubcool < rk.hydrateLowK
                                    ? 1
                                    : (r.hydrateSubcool < rk.hydrateMediumK
                                           ? 2
                                           : 3));
    }

    // ---- corrosion: the declared sour fractions at the highest pressure --
    if (rk.corrosionOn)
    {
        r.corrosionOn = true;
        const std::vector<double> p = cellPressure();
        double pmax = 0.0;
        for (int i = 0; i < n_; ++i) pmax = std::max(pmax, p[i]);
        // Without the gas EOS there is no section pressure at all; the loader
        // refuses that declaration, so pmax > 0 here by construction.
        r.co2Partial = rk.co2Fraction * pmax;
        r.h2sPartial = rk.h2sFraction * pmax;
        r.co2Band = r.co2Partial < rk.co2LowPa
                        ? 0
                        : (r.co2Partial < rk.co2HighPa ? 1 : 2);
        r.sour = r.h2sPartial >= 3447.0; // 0.5 psia
    }
    return r;
}

std::vector<double> TwoFluidPipe::cellPressure() const
{    std::vector<double> p(n_, 0.0);
    if (!eosOn()) return p;
    for (int i = 0; i < n_; ++i)
    {
        const double ag = 1.0 - alpha_[i];
        // a section with no gas volume has no meaningful pressure
        p[i] = ag > 1e-12 ? gasMass_[i] * cfg_.fluid.gasR
                            * sectionTemperature(i) / (ag * dxs_[i])
                          : 0.0;
    }
    return p;
}

double TwoFluidPipe::cellGasDensity(int i) const
{
    const double ag = 1.0 - alpha_[i];
    return eosOn() && ag > 1e-12 ? gasMass_[i] / (ag * dxs_[i])
                                 : cfg_.fluid.rhoG;
}

double TwoFluidPipe::ghostGasDensity(const TwoFluidBoundary& b) const
{
    if (!eosOn()) return cfg_.fluid.rhoG;
    return gasDensityOf(b.pressure,
                        b.temperature > 0 ? b.temperature : cfg_.temperature0);
}

// Emulsion mixing viscosity (G8): the transported water cut prices the
// Newtonian liquid share of the wall law. Both laws are exactly muL at
// phi = 0, so a dry line - and the whole standing suite, which never
// declares emulsion{} - reads the old number verbatim.
double TwoFluidPipe::liquidViscosityAt(int i) const
{
    const TwoFluidFluid& f = cfg_.fluid;
    if (!cfg_.emulsionOn) return f.muL;
    const double phi = std::min(std::max(water_[i], 0.0), 1.0);
    if (cfg_.emulsion.model == 1) // brinkman: mu_oil (1 + 2.5 phi)
        return f.muL * (1.0 + 2.5 * phi);
    // richardson: mu_oil^(1-phi) mu_water^phi (log-linear product rule)
    return std::pow(f.muL, 1.0 - phi)
         * std::pow(cfg_.emulsion.muWater, phi);
}

double TwoFluidPipe::cellFriction(int i, double vM) const
{
    // homogeneous Churchill law on mixture properties; |dp| per cell
    const TwoFluidFluid& f = cfg_.fluid;
    const double a = alpha_[i];
    const double rho = mixtureRhoCell(i);
    const double v = std::abs(vM);
    if (v == 0) return 0.0;
    // 5.98 (D2 increment 3): with feedback off, d == diameter_ and
    // rr == relrough_[i] EXACTLY, so every expression below is the old law
    // verbatim. With feedback on the deposit sets both the wall diameter and
    // the exposed roughness; the thin-layer approximation leaves area() and
    // the transport faces at the declared bore.
    // 5.219/5.220: a blockage layer or bore wear also moves the wall
    // diameter; the exposed roughness stays the wax deposit's only when a
    // deposit feedback owns the wall, otherwise the pipe's own roughness.
    const bool depositFb = (cfg_.waxOn && cfg_.wax.feedback)
                        || (cfg_.sandOn && cfg_.sand.feedback);
    const bool fb = depositFb || cfg_.blockageOn || cfg_.hydraulicFeedbackOn;
    const double d = fb ? dEff_[i] : diameter_;
    const double rr = depositFb ? cfg_.wax.depositRoughness / d
                                : relrough_[i];
    if (f.flowIndex <= 0 && !(f.yieldStress > 0))
    {
        // 5.75/5.77: the undeclared branch is the old law verbatim -- the
        // whole standing two-fluid suite is its bit-for-bit witness. With
        // emulsion{} the liquid share is the mixing law on the water cut
        // (exactly muL at phi = 0 and when undeclared, so this stays it).
        const double mu = a * liquidViscosityAt(i) + (1.0 - a) * f.muG;
        const double re = rho * v * d / mu;
        const double fr = churchillFriction(re, rr)
                        * (dxs_[i] / d) * 0.5 * rho * v * v;
        return fr;
    }
    // 5.75 (A5): a declared power-law liquid shares the wall with the gas.
    // The liquid's share of the wall stress is Metzner-Reed in its own wall
    // rate (K' = K ((3n+1)/4n)^n, so alpha = 1 reads back the exact
    // single-phase power-law pipe law) and the gas share stays Newtonian;
    // the apparent mixture viscosity tau_w/gamma drives the same Churchill
    // blend, whose laminar member f = 64/Re reads back exactly
    // (4 dx/d) tau_w -- the same identity the undeclared branch has always
    // had, now with a rate-dependent viscosity.
    // 5.77: the liquid share also carries the static wall stress tau_0.
    // As v -> 0+ the apparent viscosity tau_w/gamma diverges, so the
    // laminar member tends to (4 dx/d) (alpha tau_0) from ABOVE -- the
    // moving-wall law's q = 0 limit is exactly the distributed threshold
    // the held row replaces, the same shape 5.49 gave the lumped row.
    const double gamma = 8.0 * v / d;
    double tau;
    if (f.flowIndex > 0)
    {
        const double kp = f.consistency * std::pow(
            (3.0 * f.flowIndex + 1.0) / (4.0 * f.flowIndex), f.flowIndex);
        tau = a * kp * std::pow(gamma, f.flowIndex)
            + (1.0 - a) * f.muG * gamma;
    }
    else
    {
        // Newtonian liquid with a yield stress (Bingham-style): keep the
        // muL gamma share, add tau_0 below. flowIndex > 0 with tau_0 = 0
        // keeps the 5.75 expression verbatim. The muL share reads the
        // emulsion law when declared (exactly muL otherwise).
        tau = a * liquidViscosityAt(i) * gamma + (1.0 - a) * f.muG * gamma;
    }
    if (f.yieldStress > 0) tau += a * f.yieldStress;
    const double muEff = tau / gamma;
    const double re = rho * v * d / muEff;
    if (f.yieldStress > 0 && re < 1e-9)
    {
        // 5.77: churchillFriction clamps Re at 1e-12. A yielded fluid's
        // mu_eff = tau_w/gamma diverges as v -> 0, so Re ~ v^2 dives
        // below that clamp and the frozen f = 64/1e-12 makes the wall
        // loss vanish with v^2 -- losing the very threshold the static
        // wall is made of. In that neighbourhood the laminar member IS
        // the whole law, and its exact identity (4 dx/d) tau_w is finite
        // and independent of v. Above it the normal Churchill blend is
        // bit-for-bit what 5.75 already used.
        return 4.0 * (dxs_[i] / d) * tau;
    }
    const double fr = churchillFriction(re, rr)
                    * (dxs_[i] / d) * 0.5 * rho * v * v;
    return fr;
}

std::vector<double> TwoFluidPipe::cellPressureLosses(double vM) const
{
    // signed pressure loss per cell: friction (odd in vM) + gravity head
    std::vector<double> loss(n_);
    for (int i = 0; i < n_; ++i)
    {
        const double grav = mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        loss[i] = (vM >= 0 ? cellFriction(i, vM) : -cellFriction(i, vM)) + grav;
    }
    return loss;
}

// 5.77 (A5 increment 3): the distributed static wall. The single-phase row
// of 5.49 holds 4 tau_0 L / d; here only the LIQUID share of each cell's
// wall can yield, and each section contributes its own arc length, so the
// threshold is a sum over the frozen holdup field rather than a geometry
// constant. It is still geometry + holdup only: neither n nor K enters.
double TwoFluidPipe::gravityHeadPa() const
{
    double sum = 0;
    for (int i = 0; i < n_; ++i)
        sum += mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
    return sum;
}

double TwoFluidPipe::yieldThresholdPa() const
{
    if (!(cfg_.fluid.yieldStress > 0)) return 0.0;
    double sum = 0;
    const bool fb = (cfg_.waxOn && cfg_.wax.feedback)
                 || (cfg_.sandOn && cfg_.sand.feedback)
                 || cfg_.blockageOn || cfg_.hydraulicFeedbackOn;
    for (int i = 0; i < n_; ++i)
        sum += 4.0 * dxs_[i] / (fb ? dEff_[i] : diameter_) * alpha_[i]
             * cfg_.fluid.yieldStress;
    return sum;
}

void TwoFluidPipe::latchStaticHold(double dpDrive)
{
    driveLatch_ = dpDrive - gravityHeadPa();
    const double thr = yieldThresholdPa();
    if (!(thr > 0))
    {
        held_ = false;
        return;
    }
    // 5.49: decide once, from the state the step STARTS in, so the drive a
    // held pipe sees cannot flip the row while the same solve is moving the
    // node pressures. A column already held stays held until the drive
    // exceeds the wall; a moving column is NOT held here -- it parks on the
    // crossing witness and this latch takes over on the next step.
    held_ = (vMOld_ == 0.0 || held_)
         && std::abs(driveLatch_) <= thr;
}

double TwoFluidPipe::restResidualPa(double dt, double dpDrive) const
{
    // the exact u -> 0 limit of momentumResidual on the branch vMOld_
    // came from: friction does not vanish, it tends to the SIGNED
    // distributed threshold (the wall law is not a function at rest).
    double rhoBar = 0;
    double slipInertia = 0;
    for (int i = 0; i < n_; ++i)
    {
        rhoBar += mixtureRhoCell(i);
        slipInertia += alpha_[i] * (1.0 - alpha_[i])
                     * (cellGasDensity(i) - cfg_.fluid.rhoL) * sDotOld_[i]
                     * dxs_[i];
    }
    rhoBar /= n_;
    const double fr = std::copysign(yieldThresholdPa(), vMOld_);
    return rhoBar * length_ * (0.0 - vMOld_) / dt
         - slipInertia + fr + gravityHeadPa() - dpDrive;
}

bool TwoFluidPipe::yieldCrossing(double dt, double u, double dpDrive) const
{
    if (!(cfg_.fluid.yieldStress > 0) || vMOld_ == 0.0) return false;
    const double drive = dpDrive - gravityHeadPa();
    if (std::abs(drive) > yieldThresholdPa()) return false;
    // 5.49's crossing, stated on the residual: on the old branch friction
    // is bounded below by the threshold, so a sub-threshold drive leaves
    // no root once the old velocity is small enough. The plain witness
    // (an iterate that went through zero) is kept for the network Newton,
    // whose unknown really can change sign.
    if (u * vMOld_ <= 0.0) return true;
    const double r0 = restResidualPa(dt, dpDrive);
    return vMOld_ > 0 ? r0 >= 0.0 : r0 <= 0.0;
}

double TwoFluidPipe::momentumScale(double dt) const
{
    double rhoBar = 0;
    for (int i = 0; i < n_; ++i) rhoBar += mixtureRhoCell(i);
    rhoBar /= n_;
    double gravSum = 0;
    for (int i = 0; i < n_; ++i)
        gravSum += mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
    return rhoBar * length_ * (std::abs(vMOld_) + 1.0) / dt
         + std::abs(gravSum) + 1.0;
}

double TwoFluidPipe::momentumResidual(double dt, double u,
                                      double dpDrive) const
{
    // implicit Euler on the mixture momentum with the holdup frozen at
    // the old time level; the dynamic closure adds the lagged slip
    // inertia alpha alpha_g (rhoG - rhoL) ds/dt of the two phase
    // momenta (it vanishes at steady state). vM is a velocity, so the
    // inertia of the column is its mass per unit area rho_m L (not
    // divided by the area - that would be the volumetric-flow-rate form)
    double rhoBar = 0;
    for (int i = 0; i < n_; ++i) rhoBar += mixtureRhoCell(i);
    rhoBar /= n_;
    const double inertia = rhoBar * length_;
    double gravSum = 0;
    double slipInertia = 0;
    for (int i = 0; i < n_; ++i)
    {
        gravSum += mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        slipInertia += alpha_[i] * (1.0 - alpha_[i])
                     * (cellGasDensity(i) - cfg_.fluid.rhoL) * sDotOld_[i]
                     * dxs_[i];
    }
    double fr = 0;
    if (u == 0.0 && cfg_.fluid.yieldStress > 0)
    {
        // 5.77: the u -> 0 limit of the yielded wall law is the SIGNED
        // distributed threshold, not zero -- the same convention 5.49's
        // frictionDp(q = 0, drive) uses. The held row replaces this
        // equation altogether, so the value only feeds Newton's stencil
        // at the rest point. With zero declaration the branch is not
        // entered and the old fr = 0 stands bit for bit.
        fr = std::copysign(yieldThresholdPa(), dpDrive - gravSum);
    }
    else
    {
        for (int i = 0; i < n_; ++i)
            fr += (u >= 0 ? cellFriction(i, u) : -cellFriction(i, u));
    }
    return inertia * (u - vMOld_) / dt - slipInertia + fr + gravSum
         - dpDrive;
}

void TwoFluidPipe::momentumSolve(double dt)
{
    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;

    // 5.219: a standing bridging plug beats the pump - the line is held at
    // rest (vM = 0) for every step the bridge stands, whatever the boundary
    // type. Undeclared (or unbridged) this branch is not even entered.
    if (cfg_.blockageOn && blocked_)
    {
        clearStaticHold();
        vM_ = 0.0;
        return;
    }

    if (in.type == TwoFluidBoundary::Type::phaseRates)
    {
        // a rate boundary imposes the velocity: no static wall to ask
        clearStaticHold();
        vM_ = in.vsl + in.vsg;
        return;
    }
    if (in.type == TwoFluidBoundary::Type::sealed)
    {
        // sealed both ends: the mixture volume flux vanishes identically
        clearStaticHold();
        vM_ = 0.0;
        return;
    }

    const double dpDrive = in.pressure - out.pressure;
    // 5.77: latch the held set once from the step's starting state. The
    // held row is `u = 0`, NOT an inversion of the moving-wall law: at
    // u -> 0+ that law tends to the threshold from above, so a
    // sub-threshold drive would otherwise be answered by an equal reverse
    // flow (the exact error 5.49 refused in the lumped row).
    latchStaticHold(dpDrive);
    if (held_)
    {
        vM_ = 0.0;
        return;
    }
    // 5.77: if the moving branch has no root left, stop BEFORE asking a
    // Newton whose central difference would straddle the threshold jump
    // (5.49's frictionSlope note). Otherwise the iterate oscillates just
    // above rest forever instead of parking.
    if (yieldCrossing(dt, vMOld_, dpDrive))
    {
        vM_ = 0.0;
        return;
    }
    auto residual = [&](double u)
    {
        return momentumResidual(dt, u, dpDrive);
    };

    double u = vMOld_;
    const double scale = momentumScale(dt) + std::abs(dpDrive);
    for (int it = 0; it < 60; ++it)
    {
        const double r = residual(u);
        if (std::abs(r) <= 1e-12 * scale) break;
        const double h = std::max(1e-6 * std::max(std::abs(u), 1e-3), 1e-12);
        double slope;
        if (std::abs(u) <= h && cfg_.fluid.yieldStress > 0)
        {
            // 5.77: stay on ONE side of rest. The yielded wall law is
            // continuous on each side but jumps across u = 0, so a
            // centred stencil would put the threshold jump itself into
            // the quotient (5.49's frictionSlope note). With the rest
            // value above, this quotient is the true viscous slope.
            const double u2 = u + (u >= 0 ? h : -h);
            slope = (residual(u2) - r) / (u2 - u);
        }
        else
        {
            slope = (residual(u + h) - residual(u - h)) / (2.0 * h);
        }
        if (!(slope > 0))
            cfFatal("twoFluid: non-monotone momentum residual");
        u -= r / slope;
    }
    if (!std::isfinite(u)) cfFatal("twoFluid: momentum solve diverged");
    // 5.77: a moving column whose iterate crossed zero under a
    // sub-threshold drive stops this step; the next step's latch sees
    // vMOld_ == 0 and takes it onto the held branch.
    vM_ = yieldCrossing(dt, u, dpDrive) ? 0.0 : u;
}

void TwoFluidPipe::slipStep(double dt)
{
    if (cfg_.closure.model == 3)
    {
        // full gas momentum per section (the distributed pressure momentum
        // form): alpha_g rhoG dx dvG/dt = alpha_g dp (the gas share of the
        // mixture row, whose alpha_g fr part cancels the proportional gas
        // wall split) - alpha_g rhoG g sinT dx - sigma dx. Substituting
        // vG = vM + alpha s and expanding the material acceleration
        // dvG/dt = aM + alpha ds/dt + vG dvG/dx gives
        //   ds/dt = [(rho_m - rhoG)(g sinT + aM) - rhoG vG dvG/dx
        //            - sigma/alpha_g]/(alpha rhoG)
        // with sigma = C s|s| semi-implicit through the origin secant. At
        // a uniform steady state the convective term vanishes identically
        // (steady fluxes are uniform, so rhoG and vG are), leaving the
        // terminal slip s* = sign(R) sqrt(alpha_g |R| / C) with
        // R = (rho_m - rhoG)(g sinT + aM). Buoyancy, the mixture
        // acceleration and the gas convection drive the slip; wall
        // friction rides the pressure.
        const double aM = (vM_ - vMOld_) / dt;
        std::vector<double> old = slip_;
        for (int i = 0; i < n_; ++i)
        {
            const double a = alpha_[i];
            const double ag = 1.0 - a;
            if (a <= 1e-12 || ag <= 1e-12)
            {
                // single-phase cell: no interface, no drag, no slip
                slip_[i] = 0.0;
                sDotOld_[i] = -old[i] / dt;
                continue;
            }
            // per-section gas density: with the EOS the drive and the drag
            // denominator are section-local
            const double rhoG = cellGasDensity(i);
            // gas convective acceleration vG dvG/dx: forward difference
            // along the flow (vM cancels in the difference), zero gradient
            // at the outlet face - the same face convention as the outlet
            // fluxes
            double conv = 0.0;
            if (i + 1 < n_)
            {
                const double vGi = vM_ + a * old[i];
                const double vGj = vM_ + alpha_[i + 1] * old[i + 1];
                conv = vGi * (vGj - vGi) / dxs_[i];
            }
            const double r = (mixtureRhoCell(i) - rhoG)
                           * (kGravity * sinT_[i] + aM)
                           - rhoG * conv;
            const double rate = cfg_.closure.interfacialDrag
                              * std::abs(old[i]) / (ag * a * rhoG);
            slip_[i] = (old[i] / dt + r / (a * rhoG))
                     / (1.0 / dt + rate);
            sDotOld_[i] = (slip_[i] - old[i]) / dt;
        }
        return;
    }
    if (cfg_.closure.model != 2) return;

    // reduced two-momentum slip ODE per section (holdup frozen):
    //   K_s ds/dt = (1/rhoG - 1/rhoL)(rho_m du/dt + loss/dx) - D sigma,
    // sigma = C s |s| the interfacial drag. The drag is semi-implicit as
    // the secant through the origin, which is unconditionally stable and
    // lands exactly on the terminal slip s* = sign(R) sqrt(|R|/(C D)).
    const TwoFluidFluid& f = cfg_.fluid;
    const double aM = (vM_ - vMOld_) / dt;
    std::vector<double> old = slip_;
    for (int i = 0; i < n_; ++i)
    {
        const double a = alpha_[i];
        const double ag = 1.0 - a;
        if (a <= 1e-12 || ag <= 1e-12)
        {
            // single-phase cell: no interface, no drag, no slip
            slip_[i] = 0.0;
            sDotOld_[i] = -old[i] / dt;
            continue;
        }
        // per-section gas density: with the EOS the slip inertia factor,
        // drive and drag denominator are section-local
        const double rhoG = cellGasDensity(i);
        const double invRhoDiff = 1.0 / rhoG - 1.0 / f.rhoL;
        const double drho = rhoG - f.rhoL;
        const double ks = 1.0 + a * ag * drho * drho / (f.rhoL * rhoG);
        const double loss = (vM_ >= 0 ? cellFriction(i, vM_)
                                      : -cellFriction(i, vM_))
                          + mixtureRhoCell(i) * kGravity * sinT_[i] * dxs_[i];
        const double drive = invRhoDiff
                           * (mixtureRhoCell(i) * aM + loss / dxs_[i]);
        const double d = 1.0 / (f.rhoL * a) + 1.0 / (rhoG * ag);
        const double rate = cfg_.closure.interfacialDrag
                          * std::abs(old[i]) * d / ks;
        slip_[i] = (old[i] / dt + drive / ks) / (1.0 / dt + rate);
        sDotOld_[i] = (slip_[i] - old[i]) / dt;
    }
}

void TwoFluidPipe::transportStep(double dt)
{
    double qa = 0, qb = 0, qc = 0;
    if (!slipClosure()) fluxCoeffs(qa, qb, qc);

    // per-face slip of the dynamic closure (ghost slip at the open
    // boundaries, cell pair average inside); unused by the algebraic
    // closures
    auto faceSlip = [&](int j)
    {
        return j == 0 ? slip_[0]
             : (j == n_ ? slip_[n_ - 1]
                        : 0.5 * (slip_[j - 1] + slip_[j]));
    };

    // closure gravitySlip (shut-in-safe algebraic variant): the drift is
    // projected by the local inclination, so the wallis/zuber coefficients
    // become face-local - the face takes the average sinT of its two
    // sections, a boundary face the adjacent section. Without the key the
    // shared unprojected coefficients above stay bit-for-bit.
    auto faceDrift = [&](int j) -> double
    {
        const double sinFace =
            j == 0 ? sinT_[0]
                   : (j == n_ ? sinT_[n_ - 1]
                              : 0.5 * (sinT_[j - 1] + sinT_[j]));
        const double base = cfg_.closure.model == 0
                                ? cfg_.closure.slipVelocity
                                : cfg_.closure.driftVelocity;
        return base * sinFace;
    };
    auto liquidCoeffsAt = [&](int j, double& fqa, double& fqb, double& fqc)
    {
        const double d = faceDrift(j);
        if (cfg_.closure.model == 0) // wallis: vG = vM + alpha*d
        {
            fqa = d;
            fqb = vM_ - d;
            fqc = 0.0;
        }
        else // zuber: vG = C0 vM + d
        {
            fqa = 0.0;
            fqb = cfg_.closure.c0 * vM_ + d;
            fqc = vM_ * (1.0 - cfg_.closure.c0) - d;
        }
    };
    const bool faceLocal = !slipClosure() && cfg_.closure.gravitySlip;

    // CFL guard for the explicit Godunov transport (monotone under
    // dt * max|F'| / dx <= 1; F' is monotone so the extremes sit at 0 and 1;
    // the dynamic closure has a face-local flux, so take the max over faces)
    double cmax = 0;
    if (slipClosure())
    {
        for (int j = 0; j <= n_; ++j)
        {
            const double s = faceSlip(j);
            cmax = std::max(cmax, std::max(std::abs(vM_ - s),
                                           std::abs(vM_ + s)));
        }
    }
    else if (faceLocal)
    {
        for (int j = 0; j <= n_; ++j)
        {
            double fqa, fqb, fqc;
            liquidCoeffsAt(j, fqa, fqb, fqc);
            cmax = std::max(cmax, std::max(std::abs(fqb),
                                           std::abs(2.0 * fqa + fqb)));
        }
    }
    else
    {
        cmax = std::max(std::abs(qb), std::abs(2.0 * qa + qb));
    }
    // entrained droplets ride the gas itself (not alpha_g vG), a speed the
    // F' bounds already cover for the algebraic closures - keep the guard
    // explicit so a closure change cannot silently outrun the CFL
    if (cfg_.entrainmentOn)
        for (int i = 0; i < n_; ++i)
            cmax = std::max(cmax,
                            std::abs(gasVelocity(alpha_[i],
                                                 slipClosure()
                                                     ? slip_[i]
                                                     : 0.0, i)));
    const double dxMin = *std::min_element(dxs_.begin(), dxs_.end());
    if (dt * cmax / dxMin > 1.0 + 1e-12)
        cfFatal("twoFluid: CFL violated (dt * |dF/dalpha| / dx = "
                + std::to_string(dt * cmax / dxMin)
                + ") - reduce dt");

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
    std::vector<double> flux(n_ + 1), wflux(n_ + 1);
    std::vector<double> hflux(n_ + 1, 0.0);
    std::vector<double> wxf(n_ + 1, 0.0); // 5.97: deposited-wax content flux
    std::vector<double> iflux(n_ + 1, 0.0); // inhibitor content flux
    std::vector<double> sflux(n_ + 1, 0.0); // G10: suspended sand content flux
    // G1: content flux alpha*z_k per transported component (empty when no
    // composition{} feed was declared - every loop below then no-ops)
    std::vector<std::vector<double>> zflux;
    if (cfg_.compositionOn)
        zflux.assign(compZ_.size(), std::vector<double>(n_ + 1, 0.0));
    std::vector<double> gvol(n_ + 1, 0.0), gmflux(n_ + 1, 0.0);
    std::vector<double> dflux(n_ + 1, 0.0); // entrained droplet content flux

    for (int j = 0; j <= n_; ++j)
    {
        if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::sealed) flux[0] = 0;
            else if (in.type == TwoFluidBoundary::Type::phaseRates)
                flux[0] = in.vsl; // injected liquid superficial velocity
            else if (slipClosure())
            {
                const double s = faceSlip(0);
                flux[0] = godunovQuadratic(s, vM_ - s, 0.0, in.alpha,
                                           alpha_[0]);
            }
            else if (faceLocal)
            {
                double fqa, fqb, fqc;
                liquidCoeffsAt(0, fqa, fqb, fqc);
                flux[0] = godunovQuadratic(fqa, fqb, fqc, in.alpha,
                                           alpha_[0]);
            }
            else
                flux[0] = godunovQuadratic(qa, qb, qc, in.alpha, alpha_[0]);
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::sealed) flux[n_] = 0;
            else if (out.type == TwoFluidBoundary::Type::phaseRates)
                flux[n_] = out.vsl; // junction-imposed split (supplied
                                    // branch fed through its outlet)
            else
                flux[n_] = liquidFlux(alpha_[n_ - 1], faceSlip(n_), n_ - 1);
        }
        else if (slipClosure())
        {
            const double s = faceSlip(j);
            flux[j] = godunovQuadratic(s, vM_ - s, 0.0, alpha_[j - 1],
                                       alpha_[j]);
        }
        else if (faceLocal)
        {
            double fqa, fqb, fqc;
            liquidCoeffsAt(j, fqa, fqb, fqc);
            flux[j] = godunovQuadratic(fqa, fqb, fqc, alpha_[j - 1],
                                       alpha_[j]);
        }
        else
        {
            flux[j] = godunovQuadratic(qa, qb, qc, alpha_[j - 1], alpha_[j]);
        }

        // water cut rides the liquid flux, upwinded on the flux direction
        if (flux[j] == 0) wflux[j] = 0;
        else if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::phaseRates)
                wflux[0] = in.waterFraction * flux[0];
            else
                wflux[0] = (flux[0] >= 0 ? in.waterFraction : water_[0])
                         * flux[0];
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::phaseRates)
                wflux[n_] = out.waterFraction * flux[n_];
            else
                wflux[n_] = water_[n_ - 1] * flux[n_];
        }
        else
        {
            wflux[j] = (flux[j] >= 0 ? water_[j - 1] : water_[j]) * flux[j];
        }

        // 5.97 (D2): the deposited-wax fraction rides the liquid flux beside
        // the water cut and hydrate, same upwind rule. The inflow composition
        // carries waxDepositFraction (distinct from the thermodynamic z).
        if (cfg_.waxOn)
        {
            if (flux[j] == 0) wxf[j] = 0;
            else if (j == 0)
            {
                if (in.type == TwoFluidBoundary::Type::phaseRates)
                    wxf[0] = in.waxDepositFraction * flux[0];
                else
                    wxf[0] = (flux[0] >= 0 ? in.waxDepositFraction : wax_[0])
                           * flux[0];
            }
            else if (j == n_)
            {
                if (out.type == TwoFluidBoundary::Type::phaseRates)
                    wxf[n_] = out.waxDepositFraction * flux[n_];
                else
                    wxf[n_] = wax_[n_ - 1] * flux[n_];
            }
            else
            {
                wxf[j] = (flux[j] >= 0 ? wax_[j - 1] : wax_[j]) * flux[j];
            }
        }

        // the inhibitor fraction rides the liquid flux with the same
        // upwind rule as the water cut (inflow composition carries
        // inhibitorFraction); placed before the hydrate section, whose
        // `continue` would otherwise skip it on a hydrate-free line
        if (cfg_.inhibitorOn)
        {
            if (flux[j] == 0) iflux[j] = 0;
            else if (j == 0)
            {
                if (in.type == TwoFluidBoundary::Type::phaseRates)
                    iflux[0] = in.inhibitorFraction * flux[0];
                else
                    iflux[0] = (flux[0] >= 0 ? in.inhibitorFraction
                                             : inhibitor_[0]) * flux[0];
            }
            else if (j == n_)
            {
                if (out.type == TwoFluidBoundary::Type::phaseRates)
                    iflux[n_] = out.inhibitorFraction * flux[n_];
                else
                    iflux[n_] = inhibitor_[n_ - 1] * flux[n_];
            }
            else
            {
                iflux[j] = (flux[j] >= 0 ? inhibitor_[j - 1]
                                         : inhibitor_[j]) * flux[j];
            }
        }

        // G10: the suspended sand fraction rides the liquid flux with the
        // same upwind rule as the water cut (inflow composition carries
        // sandFraction); the bed is static and never enters the flux
        if (cfg_.sandOn)
        {
            if (flux[j] == 0) sflux[j] = 0;
            else if (j == 0)
            {
                if (in.type == TwoFluidBoundary::Type::phaseRates)
                    sflux[0] = in.sandFraction * flux[0];
                else
                    sflux[0] = (flux[0] >= 0 ? in.sandFraction
                                             : sandFrac_[0]) * flux[0];
            }
            else if (j == n_)
            {
                if (out.type == TwoFluidBoundary::Type::phaseRates)
                    sflux[n_] = out.sandFraction * flux[n_];
                else
                    sflux[n_] = sandFrac_[n_ - 1] * flux[n_];
            }
            else
            {
                sflux[j] = (flux[j] >= 0 ? sandFrac_[j - 1]
                                         : sandFrac_[j]) * flux[j];
            }
        }

        // G1: every transported mole fraction rides the liquid flux with
        // the same upwind rule as the sand content (the inflow composition
        // carries zInlet); the closed last fraction needs no flux
        for (size_t k = 0; k < zflux.size(); ++k)
        {
            if (flux[j] == 0) { zflux[k][j] = 0; continue; }
            if (j == 0)
            {
                if (in.type == TwoFluidBoundary::Type::phaseRates)
                    zflux[k][0] = in.zInlet[k] * flux[0];
                else
                    zflux[k][0] = (flux[0] >= 0 ? in.zInlet[k]
                                                : compZ_[k][0]) * flux[0];
            }
            else if (j == n_)
            {
                if (out.type == TwoFluidBoundary::Type::phaseRates)
                    zflux[k][n_] = out.zInlet[k] * flux[n_];
                else
                    zflux[k][n_] = compZ_[k][n_ - 1] * flux[n_];
            }
            else
            {
                zflux[k][j] = (flux[j] >= 0 ? compZ_[k][j - 1]
                                            : compZ_[k][j]) * flux[j];
            }
        }

        // hydrate fraction rides the liquid flux beside the water cut,
        // same upwind rule (the inflow composition carries hydrateFraction)
        if (!cfg_.hydrateOn) continue;
        if (flux[j] == 0) hflux[j] = 0;
        else if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::phaseRates)
                hflux[0] = in.hydrateFraction * flux[0];
            else
                hflux[0] = (flux[0] >= 0 ? in.hydrateFraction : hydrate_[0])
                         * flux[0];
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::phaseRates)
                hflux[n_] = out.hydrateFraction * flux[n_];
            else
                hflux[n_] = hydrate_[n_ - 1] * flux[n_];
        }
        else
        {
            hflux[j] = (flux[j] >= 0 ? hydrate_[j - 1] : hydrate_[j])
                     * flux[j];
        }
    }

    // gas volume flux per face. The gas inventory is the pipe volume
    // minus the liquid volume, so the gas transport equation is the
    // exact complement of the liquid one and its consistent face flux is
    // G = vM - F - NOT an independent Riemann solve of the mirrored
    // quadratic: the mirror G(a) = vM - F(a) is pointwise exact for every
    // closure, but the Godunov branch selection of the mirrored
    // quadratic does not mirror the liquid one (the qa == 0 linear leg
    // upwinds the opposite state, the rising/falling jumps swap the
    // shock and rarefaction branches), so an independent solve breaks
    // F + G = vM on every non-uniform face - the gas volume flux then
    // exceeds or falls short of the mixture flux, over/under-sweeps gas
    // at each front and drifts the EOS gas pressure. Interior faces and
    // a fixed-pressure inlet therefore take the one-line identity below.
    // The imposed faces keep their exact imposed values, which satisfy
    // the same identity by construction (verified per boundary): a
    // phaseRates inlet carries vM = vsl + vsg, so in.vsg IS vM - vsl
    // (kept verbatim rather than re-derived, which would round it off
    // the imposed rate by one ulp); a phaseRates outlet receives the
    // junction split from the pipe's own flux functions at the mixed
    // holdup, where liquidFluxAt + gasVolumeFluxAt = vM holds for every
    // closure; a sealed pipe carries vM = 0 with F = 0, so the identity
    // answers 0 there verbatim; an open outlet keeps the pointwise
    // closure, where gasVolumeFlux == vM - liquidFlux holds identically.
    // The volume flux exists regardless of the EOS (the thermal step
    // advects gas enthalpy with it); with the EOS the mass flux upwinds
    // the density from the side the volume flux comes from: the injected
    // reference density at a rate boundary, the ghost density at a
    // pressure boundary, the cell density inside.
    for (int j = 0; j <= n_; ++j)
    {
        if (j == 0 && in.type == TwoFluidBoundary::Type::phaseRates)
            gvol[0] = in.vsg; // the injected rate, the identity verbatim
        else if (j == n_
                 && out.type == TwoFluidBoundary::Type::phaseRates)
            gvol[n_] = out.vsg; // junction-imposed split, F + G = vM
        else if (j == n_
                 && out.type != TwoFluidBoundary::Type::sealed)
            gvol[n_] = gasVolumeFlux(alpha_[n_ - 1], faceSlip(n_),
                                     n_ - 1); // == vM - F pointwise
        else
            gvol[j] = vM_ - flux[j]; // interior + fixedPressure inlet
                                     // + sealed faces (vM = 0, F = 0)

        if (eosOn())
        {
            // rate-type faces inject the boundary's gas density: the
            // fluid reference rhoG, or the junction mix (rhoGInjected)
            const double rhoIn0 = cfg_.inlet.rhoGInjected > 0
                                      ? cfg_.inlet.rhoGInjected
                                      : cfg_.fluid.rhoG;
            const double rhoInN = cfg_.outlet.rhoGInjected > 0
                                      ? cfg_.outlet.rhoGInjected
                                      : cfg_.fluid.rhoG;
            const double rhoUp =
                gvol[j] >= 0
                    ? (j == 0
                           ? (in.type == TwoFluidBoundary::Type::phaseRates
                                  ? rhoIn0
                                  : ghostGasDensity(in))
                           : j == n_ ? cellGasDensity(n_ - 1)
                                     : cellGasDensity(j - 1))
                    : (j == n_ ? (out.type == TwoFluidBoundary::Type::phaseRates
                                      ? rhoInN
                                      : ghostGasDensity(out))
                               : cellGasDensity(j));
            gmflux[j] = rhoUp * gvol[j];
        }
    }

    // entrained droplets ride the gas at its velocity: the droplet content
    // flux qd vG (qd = alpha d) upwinds on the GAS direction (droplets go
    // where the gas goes), the film flux F - qd vG completes the liquid
    // flux F so film + droplet + gas still carry the mixture velocity; rate
    // faces impose the inflow composition's share of the liquid rate
    // instead (the injected liquid arrives entrainedFraction entrained)
    if (cfg_.entrainmentOn)
    {
        for (int j = 0; j <= n_; ++j)
        {
            if (j == 0)
            {
                if (in.type == TwoFluidBoundary::Type::phaseRates)
                    dflux[0] = in.entrainedFraction * flux[0];
                else if (in.type != TwoFluidBoundary::Type::sealed)
                    dflux[0] = gvol[0] > 0
                        ? in.alpha * in.entrainedFraction
                              * gasVelocity(in.alpha, faceSlip(0))
                        : (gvol[0] < 0
                               ? alpha_[0] * entrained_[0]
                                     * gasVelocity(alpha_[0], faceSlip(0))
                               : 0.0);
            }
            else if (j == n_)
            {
                if (out.type == TwoFluidBoundary::Type::phaseRates)
                    dflux[n_] = out.entrainedFraction * flux[n_];
                else if (out.type != TwoFluidBoundary::Type::sealed)
                    dflux[n_] = gvol[n_] < 0
                        ? out.alpha * out.entrainedFraction
                              * gasVelocity(out.alpha, faceSlip(n_))
                        : (gvol[n_] > 0
                               ? alpha_[n_ - 1] * entrained_[n_ - 1]
                                     * gasVelocity(alpha_[n_ - 1],
                                                   faceSlip(n_))
                               : 0.0);
            }
            else if (gvol[j] > 0)
            {
                dflux[j] = alpha_[j - 1] * entrained_[j - 1]
                         * gasVelocity(alpha_[j - 1], faceSlip(j));
            }
            else if (gvol[j] < 0)
            {
                dflux[j] = alpha_[j] * entrained_[j]
                         * gasVelocity(alpha_[j], faceSlip(j));
            }
        }
    }
    liquidFaceFlux_ = flux;
    gasFaceFlux_ = gvol;

    static constexpr double kTol = 1e-9;
    for (int i = 0; i < n_; ++i)
    {
        const double aNew = alpha_[i] + dt / dxs_[i] * (flux[i] - flux[i + 1]);
        if (aNew < -kTol || aNew > 1.0 + kTol)
            cfFatal("twoFluid: holdup left [0,1] - reduce dt");
        const double afOld = alpha_[i] * water_[i];
        const double ahOld = cfg_.hydrateOn ? alpha_[i] * hydrate_[i] : 0.0;
        const double awOld = cfg_.waxOn ? alpha_[i] * wax_[i] : 0.0;
        const double aentOld
            = cfg_.entrainmentOn ? alpha_[i] * entrained_[i] : 0.0;
        const double aiOld = cfg_.inhibitorOn
                                 ? alpha_[i] * inhibitor_[i] : 0.0;
        const double asOld = cfg_.sandOn ? alpha_[i] * sandFrac_[i] : 0.0;
        const double aOld = alpha_[i]; // G1: pre-transport holdup for the
                                       // composition contents
        double afNew = afOld + dt / dxs_[i] * (wflux[i] - wflux[i + 1]);
        const double aCap = std::max(aNew, 0.0);
        if (afNew < -kTol || afNew > aCap + kTol)
            cfFatal("twoFluid: water cut left [0, holdup] - reduce dt");
        alpha_[i] = std::min(std::max(aNew, 0.0), 1.0);
        afNew = std::min(std::max(afNew, 0.0), aCap);
        water_[i] = aCap > kTol ? afNew / aCap : 0.0;
        if (cfg_.hydrateOn)
        {
            // hydrate rides the liquid with the same upwind rule and the
            // same clamp semantics as the water cut (ahOld was taken
            // against the pre-transport holdup, beside afOld)
            double ahNew = ahOld
                         + dt / dxs_[i] * (hflux[i] - hflux[i + 1]);
            ahNew = std::min(std::max(ahNew, 0.0), aCap);
            hydrate_[i] = aCap > kTol ? ahNew / aCap : 0.0;
        }
        if (cfg_.waxOn)
        {
            // 5.97: the deposited-wax content alpha*w rides the liquid with
            // the same upwind rule and clamp semantics as water / hydrate
            // (awOld was captured against the pre-transport holdup, beside
            // afOld / ahOld).
            double awNew = awOld + dt / dxs_[i] * (wxf[i] - wxf[i + 1]);
            awNew = std::min(std::max(awNew, 0.0), aCap);
            wax_[i] = aCap > kTol ? awNew / aCap : 0.0;
        }
        if (cfg_.entrainmentOn)
        {
            // the droplet content rides its own flux with the same clamp
            // semantics (aentOld against the pre-transport holdup); the
            // film content alpha (1 - d) then follows by linearity
            double aentNew = aentOld
                           + dt / dxs_[i] * (dflux[i] - dflux[i + 1]);
            aentNew = std::min(std::max(aentNew, 0.0), aCap);
            entrained_[i] = aCap > kTol ? aentNew / aCap : 0.0;
        }
        if (cfg_.inhibitorOn)
        {
            // the inhibitor content alpha*c rides the liquid with the same
            // upwind rule and clamp semantics as the water cut
            double aiNew = aiOld + dt / dxs_[i] * (iflux[i] - iflux[i + 1]);
            aiNew = std::min(std::max(aiNew, 0.0), aCap);
            inhibitor_[i] = aCap > kTol ? aiNew / aCap : 0.0;
        }
        if (cfg_.sandOn)
        {
            // the suspended sand content alpha*S rides the liquid with the
            // upwind rule. The bath can hold at most aCap of sand, and sand
            // cannot be suspended into nothing: the previous clamp
            // min(max(asNew,0), aCap) SILENTLY DROPPED the mass that did not
            // fit, so a step that over-filled a cell lost sand out of the
            // alpha*S + B inventory with nothing said. The excess now settles
            // into that cell's bed, which is where a physically over-full
            // suspension puts it, and the deficit is refused like the water
            // cut's own negative transient (a negative content means dt is
            // too large for this dx, not a state to clamp).
            double asNew = asOld + dt / dxs_[i] * (sflux[i] - sflux[i + 1]);
            if (asNew < -kTol)
                cfFatal("twoFluid: suspended sand content went negative at cell "
                        + std::to_string(i) + " (asNew "
                        + std::to_string(asNew) + " m3/m3) - reduce dt or the "
                        "erodeRate: the clamp that used to swallow this hid the "
                        "loss, not the violation");
            if (aCap <= kTol)
            {
                // no liquid left in this cell: `asNew / aCap` was written as a
                // hard zero, so every grain the cell still carried disappeared
                // from the alpha*S + B inventory. Dry sand settles: it goes to
                // the bed and the readout says how much.
                if (asNew > 0.0) sandClampToBed_ += asNew;
                bed_[i] += std::max(asNew, 0.0);
                sandFrac_[i] = 0.0;
            }
            else
            {
                if (asNew > aCap)
                {
                    const double excess = asNew - aCap;
                    bed_[i] += excess;
                    sandClampToBed_ += excess;
                    asNew = aCap;
                }
                sandFrac_[i] = asNew / aCap;
            }
        }
        for (size_t k = 0; k < compZ_.size(); ++k)
        {
            // G1: the content alpha*z_k rides the liquid with the same
            // upwind rule and clamp semantics as the sand fraction; the
            // closed fraction 1 - sum stays in [0,1] because the upwind
            // update is a convex combination of cell and inflow states
            // (all of which carry normalised feeds)
            const double azOld = aOld * compZ_[k][i];
            double azNew = azOld + dt / dxs_[i] * (zflux[k][i]
                                                   - zflux[k][i + 1]);
            azNew = std::min(std::max(azNew, 0.0), aCap);
            compZ_[k][i] = aCap > kTol ? azNew / aCap : 0.0;
        }
        if (eosOn())
        {
            // gas mass continuity on the per-unit-area inventory
            // rhoG (1 - alpha) dx: the mass-flux divergence changes it
            // directly (no dx factor); upwinding keeps it non-negative
            const double agNew = 1.0 - alpha_[i];
            const double numer
                = gasMass_[i] - dt * (gmflux[i + 1] - gmflux[i]);
            const double trace = kTol * cfg_.fluid.rhoG * dxs_[i];
            if (agNew <= 1e-12 && numer > trace)
                cfFatal("twoFluid: gas sealed in a closing cell - reduce dt");
            if (numer < -trace)
                cfFatal("twoFluid: gas mass went negative - reduce dt");
            gasMass_[i] = std::max(numer, 0.0);
        }
    }
}

void TwoFluidPipe::entrainmentStep(double dt)
{
    if (!cfg_.entrainmentOn) return;
    static constexpr double kTol = 1e-9;
    const TwoFluidEntrainment& en = cfg_.entrainment;
    for (int i = 0; i < n_; ++i)
    {
        if (alpha_[i] <= kTol) continue; // no liquid: nothing to split
        // the equilibrium fraction reads the closure gas velocity of the
        // just-transported state; implicit Euler on the linear relaxation
        // dd/dt = (E_eq(vG) - d) / relaxTime
        const double vG = gasVelocity(alpha_[i],
                                      slipClosure() ? slip_[i] : 0.0, i);
        const double eEq = std::min(
            std::max((vG - en.vCrit) / en.vRef, 0.0), 1.0);
        entrained_[i] += (eEq - entrained_[i]) * dt / (en.relaxTime + dt);
    }
}

void TwoFluidPipe::hydrateStep(double dt)
{
    if (!cfg_.hydrateOn) return;
    static constexpr double kTol = 1e-9;
    const TwoFluidHydrate& hy = cfg_.hydrate;
    // the subcooling uses the just-transported state: the EOS pressure of
    // the transported gas inventory and the current section temperature
    const std::vector<double> p = cellPressure();
    for (int i = 0; i < n_; ++i)
    {
        const double ag = 1.0 - alpha_[i];
        if (ag <= 1e-12)
            cfFatal("twoFluid: hydrate kinetics needs gas volume in every"
                    " cell - keep a gas-bearing holdup");
        // with inhibitor{} the equilibrium line is depressed at the local
        // transported fraction - the Hammerschmidt linear term, or the
        // rigorous activity root when inhibitorEquilibrium{} is declared
        // (5.163); without, the declared curve verbatim (exactly the old
        // difference)
        const double dT = equilibriumTemperatureInhibited(p[i], i)
                        - sectionTemperature(i);
        const double hOld = hydrate_[i];
        // hydrate-able slurry material: free water plus existing hydrate
        // (oil is inert); the growth law converges to the full conversion
        const double total = water_[i] + hOld;
        double hNew;
        if (dT >= 0)
        {
            // implicit Euler on dh/dt = kg dT (total - h)
            hNew = (hOld / dt + hy.growthRate * dT * total)
                 / (1.0 / dt + hy.growthRate * dT);
            // gas availability: gasContent kg per m3 of hydrate off the EOS
            // inventory, capped so the mass stays non-negative
            const double trace = kTol * cfg_.fluid.rhoG * dxs_[i];
            const double dhMax = (gasMass_[i] - trace)
                               / (hy.gasContent * alpha_[i] * dxs_[i]);
            hNew = std::min(hNew, hOld + std::max(dhMax, 0.0));
        }
        else
        {
            // implicit Euler on dh/dt = -kd |dT| h
            hNew = hOld / (1.0 + dt * hy.dissociationRate * (-dT));
        }
        hNew = std::min(std::max(hNew, 0.0), total);
        const double dh = hNew - hOld;
        // 1:1 volume move between free water and hydrate (the
        // volume-preserving convention of the network kernel)
        water_[i] = std::min(std::max(total - hNew, 0.0), 1.0);
        hydrate_[i] = hNew;
        // stoichiometric gas mass sink on the per-area EOS inventory
        // (reappears on dissociation)
        gasMass_[i] -= hy.gasContent * alpha_[i] * dh * dxs_[i];
        // hydrate volume change per bulk m3: the latent-heat source the
        // thermal step consumes (positive = release)
        hydrateDelta_[i] = alpha_[i] * dh;
    }
}

// G10: sand / cuttings exchange between the suspended load (rides the
// liquid, transported in transportStep) and the static bed. One-sided
// switch on the mixture velocity: below the declared critical velocity the
// suspended fraction settles across half the bore (rate vt/(D/2), implicit
// Euler), above it the bed erodes at erodeRate (vM/uCrit - 1) (implicit
// Euler). The cell's total sand volume alpha*S + B changes only by
// transport; both sides are exactly reproducible closed forms.
void TwoFluidPipe::sandStep(double dt)
{
    if (!cfg_.sandOn) return;
    const TwoFluidSand& sd = cfg_.sand;
    const double vM = vM_;
    if (vM < sd.criticalVelocity)
    {
        const double k = sd.terminalVelocity / (0.5 * diameter_);
        for (int i = 0; i < n_; ++i)
        {
            const double sNew = sandFrac_[i] / (1.0 + dt * k);
            bed_[i] += alpha_[i] * (sandFrac_[i] - sNew);
            sandFrac_[i] = sNew;
        }
    }
    else if (sd.erodeRate > 0)
    {
        const double k = sd.erodeRate * (vM / sd.criticalVelocity - 1.0);
        const double kEps = 1e-12;
        for (int i = 0; i < n_; ++i)
        {
            const double bRelaxed = bed_[i] / (1.0 + dt * k);
            if (alpha_[i] <= kEps)
            {
                // no liquid to suspend into: the bed keeps the sand rather
                // than silently dropping it
                continue;
            }
            const double eroded = bed_[i] - bRelaxed;
            const double capacity = alpha_[i] * std::max(0.0, 1.0 - sandFrac_[i]);
            const double moved = std::min(eroded, capacity);
            sandFrac_[i] += moved / alpha_[i];
            bed_[i] -= moved; // whatever the suspension cannot take stays put
        }
    }
}

void TwoFluidPipe::compositionStep(double)
{
    if (!cfg_.compositionOn) return;
    // G1: per-cell PR flash of the transported feed at the cell's own
    // (T, p) - the transient counterpart of the steady liveFlash (5.121).
    // Pure diagnostic: the result lands in compBeta_/compRhoL_/compRhoG_
    // for the report and the cellVaporFraction() accessor; the momentum,
    // holdup and density fields never read it.
    const std::vector<double> p = eosOn() ? cellPressure()
                                          : std::vector<double>();
    const int nc = static_cast<int>(cfg_.compositionZ0.size());
    std::vector<double> z(nc, 0.0);
    for (int i = 0; i < n_; ++i)
    {
        double s = 0.0;
        for (int k = 0; k < nc - 1; ++k)
        {
            z[k] = compZ_[k][i];
            s += z[k];
        }
        z[nc - 1] = std::max(0.0, 1.0 - s); // closed fraction
        // flash pressure: the cell's EOS pressure where available, else
        // the declared fixedPressure reference (5.143)
        double pi = flashPRef_;
        if (eosOn() && p[i] > 0) pi = p[i];
        if (!(pi > 0))
            cfFatal("twoFluid: composition flash: no pressure in cell "
                    + std::to_string(i) + " (liquid-full cell with no"
                    " fixedPressure reference) - reduce dt or declare a"
                    " fixedPressure boundary");
        pvt::Mixture feed;
        for (int k = 0; k < nc; ++k)
            feed.add(pvt::requireComponent(cfg_.compositionNames[k]), z[k]);
        const pvt::FlashResult fr
            = pvt::tpFlash(feed, sectionTemperature(i), pi);
        // the steady liveFlash precedent (5.121) prices the returned
        // state without rejecting an unconverged RR iterate; the result
        // here is diagnostic-only, so do the same
        if (!(fr.liquidMassDensity > 0.0) || !(fr.vapourMassDensity > 0.0))
            cfFatal("twoFluid: composition flash: non-positive phase"
                    " density in cell " + std::to_string(i));
        compBeta_[i] = fr.vaporFraction;
        compRhoL_[i] = fr.liquidMassDensity;
        compRhoG_[i] = fr.vapourMassDensity;
    }
}

std::vector<double> TwoFluidPipe::cellComposition(int k) const
{
    if (!cfg_.compositionOn) return {};
    const int nc = static_cast<int>(cfg_.compositionZ0.size());
    if (k < 0 || k >= nc)
        cfFatal("twoFluid: cellComposition: component index out of range");
    if (k < nc - 1) return compZ_[k];
    std::vector<double> out(n_, 0.0); // the closed fraction
    for (int i = 0; i < n_; ++i)
    {
        double s = 0.0;
        for (int j = 0; j < nc - 1; ++j) s += compZ_[j][i];
        out[i] = 1.0 - s;
    }
    return out;
}

void TwoFluidPipe::waxStep(double dt)
{
    if (!cfg_.waxOn) return;
    const TwoFluidWax& wx = cfg_.wax;
    // The LOCAL reaction is per section on its own current temperature;
    // convection of the deposited fraction happens in transportStep (5.97).
    // The dissolved wax is a per-cell store: total wax z (the declared
    // fraction) splits into dissolved c and deposited w, c + w = z. The
    // van't Hoff equilibrium leaves c = min(z, x_sat(T)) dissolved, so the
    // equilibrium deposit is
    //   w_eq = max(0, z - x_sat(T))
    // and the cell relaxes to it at first order:
    //   dw/dt = -k (w - w_eq),  k = growthRate if w < w_eq (precipitation),
    //                                dissolutionRate otherwise.
    // Implicit Euler keeps w in [0, w_eq] <= z unconditionally and is exactly
    // reproducible by a checker using nothing but exp and the two declared
    // constants: w_new = (w/dt + k w_eq) / (1/dt + k).
    for (int i = 0; i < n_; ++i)
    {
        const double xSat = pvt::waxSaturationMoleFraction(
            sectionTemperature(i), wx.meltingTemperature, wx.fusionEnthalpy);
        const double drive = wx.waxFraction - xSat;
        const double wEq = std::max(drive, 0.0);
        const double k = drive > 0 ? wx.growthRate : wx.dissolutionRate;
        const double wNew = std::min(
            std::max((wax_[i] / dt + k * wEq) / (1.0 / dt + k), 0.0),
            wx.waxFraction);
        const double dw = wNew - wax_[i];
        wax_[i] = wNew;
        // precipitation (dw > 0) releases latent heat per m3 of liquid wax;
        // dissolution absorbs it. The thermal step consumes this source.
        waxDelta_[i] = alpha_[i] * dw;
    }
}

// 5.98/5.99: one owner for the thin-layer diameter. The deposit that shrinks
// the wall-shear bore is the TOTAL of every depositing species (wax +
// asphaltene + scale); with only wax declared this is exactly the 5.98
// behavior, bit for bit. 5.219: with blockage{} the hydrate layer runs here
// and is charged linearly on top (d_eff = d_deposit - 2 * delta).
void TwoFluidPipe::refreshDepositDiameter(double dt)
{
    const bool waxFb = cfg_.waxOn && cfg_.wax.feedback;
    const bool sandFb = cfg_.sandOn && cfg_.sand.feedback;
    if (!waxFb && !sandFb && !cfg_.blockageOn && !cfg_.hydraulicFeedbackOn)
        return;
    if (cfg_.blockageOn) hydrateBlockageStep(dt);
    for (int i = 0; i < n_; ++i)
    {
        const double w = (cfg_.waxOn ? wax_[i] : 0.0)
            + (cfg_.asphalteneOn ? asph_[i] : 0.0)
            + (cfg_.scaleOn ? scale_[i] : 0.0)
            + (sandFb ? bed_[i] : 0.0);
        dEff_[i] = diameter_ * std::sqrt(std::max(1.0 - w, 1e-12));
        if (cfg_.blockageOn)
            dEff_[i] = std::max(dEff_[i] - 2.0 * blockage_[i], 1e-12);
    }
}

// 5.219: hydrate blockage. The wall layer obeys the declared deposition /
// erosion balance, implicit Euler exactly as the wax/scale relaxations:
//   delta_new = (delta_old / dt + depositionRate * h) / (1 / dt + erosionRate)
// (erosionRate = 0 degenerates to the pure accumulation delta + dt*dep*h).
// The layer cannot eat more than the bore: it is clamped below D/2 so the
// linear charge keeps d_eff positive. The bridging arm is a latch on the
// FINAL hydrate field of the step: where any section reaches the packing
// fraction the line plugs (vM = 0 on the spot and held by momentumSolve
// for every step the bridge stands); dissociation back below the limit
// releases it.
void TwoFluidPipe::hydrateBlockageStep(double dt)
{
    if (dt <= 0) return;
    const TwoFluidBlockage& bk = cfg_.blockage;
    const double deltaMax = 0.5 * diameter_ * (1.0 - 1e-9);
    bool bridged = false;
    for (int i = 0; i < n_; ++i)
    {
        double d = (blockage_[i] / dt + bk.depositionRate * hydrate_[i])
                 / (1.0 / dt + bk.erosionRate);
        d = std::min(std::max(d, 0.0), deltaMax);
        blockage_[i] = d;
        if (hydrate_[i] >= bk.packingFraction) bridged = true;
    }
    if (bridged) vM_ = 0.0; // the step that forms the bridge ends plugged
    blocked_ = bridged;
}

// 5.220: erosion / corrosion hydraulic feedback. The bore GROWS where the
// risk is live: the section's own API RP 14E limit (the 5.187 law, priced
// on this frame's mixture density) is crossed by |vM|, or the declared H2S
// fraction at the highest section EOS pressure sits on the sour side of the
// published 0.5 psia = 3447 Pa line. The growth is the declared explicit
// accumulation (rate * dt per step), so the checker re-derives every frame
// from the reported state. Applied on top of the deposit diameter.
void TwoFluidPipe::erosionCorrosionFeedbackStep(double dt)
{
    if (!cfg_.hydraulicFeedbackOn) return;
    const TwoFluidHydraulicFeedback& fb = cfg_.hydraulicFeedback;
    const std::vector<double> p = cellPressure();
    double pmax = 0.0;
    for (int i = 0; i < n_; ++i) pmax = std::max(pmax, p[i]);
    const bool sour = fb.h2sFraction * pmax >= 3447.0; // 0.5 psia
    for (int i = 0; i < n_; ++i)
    {
        const double rhoM = mixtureRhoCell(i);
        const double rhoLb = rhoM / 16.0184634;
        const double velim = rhoLb > 0.0
            ? 0.3048 * fb.erosionC / std::sqrt(rhoLb) : 0.0;
        const double ratio = velim > 0.0 ? std::fabs(vM_) / velim : 0.0;
        wear_[i] += dt * (fb.erosionRate * (ratio > 1.0 ? 1.0 : 0.0)
                        + fb.corrosionRate * (sour ? 1.0 : 0.0));
        dEff_[i] += wear_[i];
    }
}

// 5.99: asphaltene onset. The equilibrium comes from the local cell-centre
// hydraulic pressure against the DECLARED AOP; the kinetics are the same
// implicit-Euler first-order relaxation waxStep uses. No latent heat, no
// transport, no junction mixing (wall-bound local field).
void TwoFluidPipe::asphalteneStep(double dt)
{
    if (!cfg_.asphalteneOn) return;
    const TwoFluidAsphaltene& az = cfg_.asphaltene;
    const std::vector<double> pF = nodePressure();
    for (int i = 0; i < n_; ++i)
    {
        const double p = 0.5 * (pF[i] + pF[i + 1]);
        const double drive = az.onsetWidth > 0
            ? (az.onsetPressure - p) / az.onsetWidth
            : (az.onsetPressure - p > 0 ? 1.0 : 0.0);
        const double aEq = az.asphalteneFraction
            * std::min(std::max(drive, 0.0), 1.0);
        const double k = asph_[i] < aEq ? az.growthRate
                                        : az.dissolutionRate;
        asph_[i] = std::min(
            std::max((asph_[i] / dt + k * aEq) / (1.0 / dt + k), 0.0),
            az.asphalteneFraction);
    }
}

// 5.99: CaCO3 scale. Omega = IAP/Ksp from the DECLARED brine chemistry; the
// deposited fraction relaxes toward z_s * clamp(Omega - 1, 0, 1) with the
// same implicit-Euler relaxation. No latent heat, no transport.
void TwoFluidPipe::scaleStep(double dt)
{
    if (!cfg_.scaleOn) return;
    const TwoFluidScale& sc = cfg_.scale;
    const double omega = sc.calciumActivity * sc.carbonateActivity
                       / sc.solubilityProduct;
    const double sEq = sc.scaleFraction * std::min(
        std::max(omega - 1.0, 0.0), 1.0);
    for (int i = 0; i < n_; ++i)
    {
        const double k = scale_[i] < sEq ? sc.growthRate
                                         : sc.dissolutionRate;
        scale_[i] = std::min(
            std::max((scale_[i] / dt + k * sEq) / (1.0 / dt + k), 0.0),
            sc.scaleFraction);
    }
}

void TwoFluidPipe::applySchedule(double dt)
{
    const std::vector<TwoFluidSchedulePoint>& sch = cfg_.inlet.schedule;
    if (sch.empty()) return;
    // the last point whose time the step interval [t_, t_ + dt] covers
    // governs; evaluation at the interval midpoint keeps that decision
    // exact even though the accumulated t_ drifts by ~1e-12 after a
    // thousand additions (a point at an exact step boundary would
    // otherwise fire a step late)
    const double tEval = t_ + 0.5 * dt;
    std::size_t k = 0;
    while (k + 1 < sch.size() && sch[k + 1].t <= tEval) ++k;
    if (sch[k].t <= tEval)
    {
        cfg_.inlet.vsl = sch[k].vsl;
        cfg_.inlet.vsg = sch[k].vsg;
    }
}

void TwoFluidPipe::step(double dt)
{
    if (!(dt > 0)) cfFatal("twoFluid: dt must be positive");
    applySchedule(dt);
    // junction branches have their mixture velocity imposed by the
    // network's coupled Newton (externalMomentum_ is a one-shot flag)
    if (externalMomentum_) externalMomentum_ = false;
    else momentumSolve(dt);
    slipStep(dt);
    transportStep(dt);
    terrainSlugStep(dt);
    purgeWaveStep(dt); // O3 (5.162): the blowdown front rides the final field
    entrainmentStep(dt);
    hydrateStep(dt);
    sandStep(dt);
    waxStep(dt);
    asphalteneStep(dt);
    scaleStep(dt);
    refreshDepositDiameter(dt);
    erosionCorrosionFeedbackStep(dt); // 5.220: bore growth rides the final state
    compositionStep(dt); // G1: per-cell PR flash diagnostics (reads only)
    thermalStep(dt);
    slugTrackStep(dt); // O2 (5.161): rides the step's final holdup field
    piggingStep(dt); // E6 (5.176): reduced pig rides the final holdup field
    acousticsStep(dt); // E7c (5.183): linear acoustic layer (report-only)
    vMOld_ = vM_;
    t_ += dt;
}

double TwoFluidPipe::wallTemperature(int i) const
{
    return wallTemperatureAt(i, t_);
}

// 5.222: with ambientSchedule{} the ambient level is the piecewise-linear
// interpolation of the declared (times, ambientT) at the evaluation time:
// the constant model's wall IS the ambient, the geothermal model keeps its
// per-section depth offset on top of the moving surface level. Without the
// schedule the time argument is ignored and the old expression returns
// verbatim (bit-identical history).
double TwoFluidPipe::wallTemperatureAt(int i, double t) const
{
    const TwoFluidThermal& th = cfg_.thermal;
    if (!cfg_.ambientScheduleOn)
    {
        if (th.wallModel == 0) return th.wallTemperature;
        const double depth = th.inletDepth - centerDrop_[i];
        return th.surfaceTemperature + th.geothermalGradient * depth;
    }
    const double amb = ambientAtTime(t);
    if (th.wallModel == 0) return amb;
    const double depth = th.inletDepth - centerDrop_[i];
    return amb + th.geothermalGradient * depth;
}

// 5.222: the declared schedule, held (clamped) outside the time window.
// The segment is the last point at or before t; the value is the exact
// linear blend on it.
double TwoFluidPipe::ambientAtTime(double t) const
{
    const std::vector<double>& ts = cfg_.ambientSchedule.times;
    const std::vector<double>& tv = cfg_.ambientSchedule.ambientT;
    if (t <= ts.front()) return tv.front();
    if (t >= ts.back()) return tv.back();
    std::size_t k = 0;
    while (k + 1 < ts.size() - 1 && ts[k + 1] <= t) ++k;
    const double f = (t - ts[k]) / (ts[k + 1] - ts[k]);
    return tv[k] + f * (tv[k + 1] - tv[k]);
}

void TwoFluidPipe::thermalStep(double dt)
{
    if (!cfg_.thermalOn) return;
    const TwoFluidFluid& f = cfg_.fluid;
    const TwoFluidThermal& th = cfg_.thermal;
    const double area = M_PI * diameter_ * diameter_ / 4.0;

    // implicit Euler on the section temperatures. The phase enthalpies ride
    // the phase volume fluxes of this step's transport - with slip the
    // liquid no longer moves at the mixture velocity, so a mixture-rhoc * vM
    // weighting would over-drag the liquid capacity against the slip - and
    // each phase upwinds its temperature on its own flux direction. The gas
    // density upwinds like the mass flux does (injected reference density at
    // a rate boundary, ghost density at a pressure boundary, cell density
    // inside).
    std::vector<std::vector<double>> mat(n_, std::vector<double>(n_, 0.0));
    std::vector<double> rhs(n_, 0.0);
    for (int i = 0; i < n_; ++i)
    {
        const double cap = area * dxs_[i]
                         * (alpha_[i] * f.rhoL * f.cL
                            + (1.0 - alpha_[i]) * cellGasDensity(i) * f.cG);
        const double kw = th.wallHeatTransfer * dxs_[i];
        // 5.222: the implicit source is evaluated at the step's END time;
        // without a schedule the expression ignores t (bit-identical)
        const double tw = cfg_.formationOn ? formT_[i][0]
                                           : wallTemperatureAt(i, t_ + dt);
        if (cfg_.ambientScheduleOn) wallTUsed_[i] = wallTemperatureAt(i, t_ + dt);
        mat[i][i] = cap / dt + kw;
        rhs[i] = cap / dt * T_[i] + kw * tw;
        if (cfg_.hydrateOn)
        {
            // hydrate formation releases its latent heat into the section
            // (dissociation absorbs it), J/s over the cell bulk volume
            rhs[i] += cfg_.hydrate.latentHeat * hydrateDelta_[i] * area
                    * dxs_[i] / dt;
        }
        if (cfg_.waxOn)
        {
            // wax precipitation releases its latent heat into the section
            // (dissolution absorbs it), J/s over the cell bulk volume
            rhs[i] += cfg_.wax.latentHeat * waxDelta_[i] * area
                    * dxs_[i] / dt;
        }
    }

    // rate-type faces inject the boundary's gas density: the fluid
    // reference rhoG, or the junction mix (rhoGInjected)
    auto injectedRhoG = [&](const TwoFluidBoundary& b) -> double
    {
        return b.rhoGInjected > 0 ? b.rhoGInjected : f.rhoG;
    };
    auto gasRho = [&](int j) -> double
    {
        if (!eosOn()) return f.rhoG;
        const double q = gasFaceFlux_[j];
        if (q > 0)
            return j == 0
                ? (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates
                       ? injectedRhoG(cfg_.inlet)
                       : ghostGasDensity(cfg_.inlet))
                : j == n_ ? cellGasDensity(n_ - 1) : cellGasDensity(j - 1);
        if (q < 0)
            return j == n_
                ? (cfg_.outlet.type == TwoFluidBoundary::Type::phaseRates
                       ? injectedRhoG(cfg_.outlet)
                       : ghostGasDensity(cfg_.outlet))
                : cellGasDensity(j);
        return 0.0;
    };
    // one face-phase contribution: implicit upwinding on the phase's own
    // flux direction, assembled in the conservative receive-capacity form:
    // the receiver's diagonal takes +k (the arrived mass joins the cell's
    // capacity before it exchanges heat), the donor takes nothing (its
    // capacity loss IS the outflux - it cancels against the capacity update
    // C_new = C_old + dt(in - out)), so diag - |off| = cap/dt + kw exactly
    // and the max principle holds unconditionally: every temperature stays
    // within [min(T_old, T_wall), max(T_old, T_wall, T_boundary)] - also
    // under persistent countercurrent slip at shut-in, where the historical
    // donor-drain form closed the row-sum margin and let temperatures
    // overshoot any physical bound (5.131). A boundary receive also adds +k
    // to the diagonal; a boundary drain (backflow out through an open
    // boundary) still assembles nothing.
    auto faceTerm = [&](int j, double q, double k)
    {
        if (q == 0) return;
        const int up = q > 0 ? j - 1 : j;
        const int dn = q > 0 ? j : j - 1;
        const bool upIn = up >= 0 && up < n_;
        const bool dnIn = dn >= 0 && dn < n_;
        if (dnIn)
        {
            if (upIn) mat[dn][up] -= k;
            else
                // boundary receive: arrived capacity joins the cell
                rhs[dn] += k * (j == 0 ? cfg_.inlet.temperature
                                       : cfg_.outlet.temperature);
            mat[dn][dn] += k;
        }
    };

    inletHeatRate_ = 0;
    for (int j = 0; j <= n_; ++j)
    {
        const double kl = area * f.rhoL * f.cL
                        * std::abs(liquidFaceFlux_[j]);
        const double kg = area * gasRho(j) * f.cG
                        * std::abs(gasFaceFlux_[j]);
        faceTerm(j, liquidFaceFlux_[j], kl);
        faceTerm(j, gasFaceFlux_[j], kg);
        if (j == 0)
        {
            if (liquidFaceFlux_[j] > 0)
                inletHeatRate_ += kl * cfg_.inlet.temperature;
            if (gasFaceFlux_[j] > 0)
                inletHeatRate_ += kg * cfg_.inlet.temperature;
        }
        else if (j == n_)
        {
            if (liquidFaceFlux_[j] < 0)
                inletHeatRate_ += kl * cfg_.outlet.temperature;
            if (gasFaceFlux_[j] < 0)
                inletHeatRate_ += kg * cfg_.outlet.temperature;
        }
    }
    detail::solveDenseLinear(mat, rhs);
    for (int i = 0; i < n_; ++i) T_[i] = rhs[i];

    // march the radial rock with the new fluid temperatures, then report the
    // wall exchange at the end-of-step wall temperature so the formation
    // energy balance (intake = storage + far-field rejection) closes exactly
    if (cfg_.formationOn) formationStep(dt);

    wallHeatRate_ = 0;
    for (int i = 0; i < n_; ++i)
    {
        const double tw = cfg_.formationOn ? formT_[i][0]
                                           : wallTemperatureAt(i, t_ + dt);
        wallHeatRate_ += th.wallHeatTransfer * dxs_[i] * (tw - T_[i]);
    }
}

void TwoFluidPipe::formationStep(double dt)
{
    const TwoFluidFormation& fm = cfg_.formation;
    const TwoFluidThermal& th = cfg_.thermal;
    const int M = fm.shells;
    formHeatRate_ = 0;
    formIntakeRate_ = 0;
    formStorageRate_ = 0;
    // 5.84: on the semi-infinite path the outer field answers the whole
    // well's recorded flux history through the E1/Duhamel kernel. The rise is
    // shared by every section, while each section keeps its own geothermal
    // far-field offset. The Dirichlet path never evaluates the kernel.
    const bool semi = fm.semiInfinite;
    formFarRise_ = semi
        ? lineSourceRise(t_, formQHist_, formTHist_, fm.outerRadius,
                         fm.conductivity / fm.rhoC, fm.conductivity)
        : 0.0;
    for (int i = 0; i < n_; ++i)
    {
        const double kw = th.wallHeatTransfer * dxs_[i];
        const double te = semi
            ? wallTemperatureAt(i, t_ + dt) + formFarRise_
            : wallTemperatureAt(i, t_ + dt);
        const std::vector<double>& g = formG_[i];
        const std::vector<double>& c = formC_[i];
        std::vector<double>& tf = formT_[i];
        // implicit Euler on the M unknown nodes; node M is the far field.
        // Node 0 is the wall: film exchange with the new fluid temperature
        // plus conduction to node 1. The last node conducts to the Dirichlet
        // far field. Internal faces telescope, so the section balance closes.
        std::vector<std::vector<double>> a(M, std::vector<double>(M, 0.0));
        std::vector<double> b(M, 0.0);
        for (int j = 0; j < M; ++j)
        {
            double diag = c[j] / dt;
            double rhs = c[j] / dt * tf[j];
            if (j == 0)
            {
                diag += kw + g[0];
                rhs += kw * T_[i];
            }
            else
            {
                diag += g[j - 1] + g[j];
                a[j][j - 1] -= g[j - 1];
            }
            if (j < M - 1)
                a[j][j + 1] -= g[j];
            else
                rhs += g[M - 1] * te;
            a[j][j] = diag;
            b[j] = rhs;
        }
        detail::solveDenseLinear(a, b);
        for (int j = 0; j < M; ++j)
        {
            formStorageRate_ += c[j] * (b[j] - tf[j]) / dt;
            tf[j] = b[j];
        }
        tf[M] = te;
        formIntakeRate_ += kw * (T_[i] - tf[0]);
        formHeatRate_ += g[M - 1] * (tf[M - 1] - te);
    }
    if (semi)
    {
        // Record the step-average far-field flux per metre of well. F(0) = 0
        // guarantees this entry cannot influence the rise of this same step.
        formTHist_.push_back(t_);
        formQHist_.push_back(formHeatRate_ / length_);
    }
}

std::vector<double> TwoFluidPipe::formationProfile(int i) const
{
    if (!cfg_.formationOn || i < 0 || i >= n_) return {};
    return formT_[i];
}

double TwoFluidPipe::pressureDrop() const
{
    const std::vector<double> p = nodePressure();
    return p.front() - p.back();
}

std::vector<double> TwoFluidPipe::nodePressure() const
{
    const std::vector<double> loss = cellPressureLosses(vM_);
    std::vector<double> p(n_ + 1, 0.0);
    if (cfg_.inlet.type == TwoFluidBoundary::Type::fixedPressure)
    {
        p[0] = cfg_.inlet.pressure;
        for (int i = 0; i < n_; ++i) p[i + 1] = p[i] - loss[i];
    }
    else
    {
        // rate-driven or sealed: gauge pressure relative to the outlet
        p[n_] = cfg_.outlet.type == TwoFluidBoundary::Type::fixedPressure
                    ? cfg_.outlet.pressure
                    : 0.0;
        for (int i = n_ - 1; i >= 0; --i) p[i] = p[i + 1] + loss[i];
    }
    return p;
}

std::vector<double> TwoFluidPipe::cellWater() const { return water_; }

std::vector<double> TwoFluidPipe::cellVsl() const
{
    // superficial liquid velocity = liquid volume flux per unit area
    std::vector<double> v(n_);
    for (int i = 0; i < n_; ++i)
        v[i] = liquidFlux(alpha_[i], slip_[i], i);
    return v;
}

std::vector<double> TwoFluidPipe::cellVsg() const
{
    std::vector<double> v(n_);
    for (int i = 0; i < n_; ++i)
        v[i] = (1.0 - alpha_[i]) * gasVelocity(alpha_[i], slip_[i], i);
    return v;
}

std::vector<TwoFluidSlugBody> TwoFluidPipe::slugRuns(double alphaThreshold,
                                                     double minLength) const
{
    std::vector<TwoFluidSlugBody> out;
    bool inRun = false;
    double rear = 0, len = 0, inv = 0, x = 0;
    auto closeRun = [&]()
    {
        if (inRun && len >= minLength)
            out.push_back({rear, rear + len, len, inv});
        inRun = false;
    };
    for (int i = 0; i < n_; ++i)
    {
        if (alpha_[i] >= alphaThreshold)
        {
            if (!inRun)
            {
                inRun = true;
                rear = x;
                len = 0;
                inv = 0;
            }
            len += dxs_[i];
            inv += alpha_[i] * dxs_[i];
        }
        else
        {
            closeRun();
        }
        x += dxs_[i];
    }
    closeRun();
    return out;
}

std::vector<TwoFluidSlugBody> TwoFluidPipe::slugBodies() const
{
    if (!cfg_.slugOn) return {};
    return slugRuns(cfg_.slug.alphaThreshold, cfg_.slug.minLength);
}

double TwoFluidPipe::bendiksenFrontVelocity(int frontCell) const
{
    // E7b (5.182): Bendiksen (1984) translational velocity, priced at the
    // cell holding the front face. v_front = C0 v_m + V_d with V_d either
    // the declared constant or the default drift law
    //   V_d = 0.35 sin(theta) sqrt(g D (rhoL - rhoG)/rhoL),
    // Bendiksen's turbulent-branch coefficient with the density-corrected
    // Taylor-bubble nose form (see TwoFluidSlugTracking for the full
    // citation and applicability boundary). The declared densities are
    // used (not the EOS per-cell gas densities) and the density difference
    // is clamped at zero so a mis-ordered declaration can never reach the
    // sqrt (the loader rejects it loudly anyway).
    const TwoFluidSlugTracking& st = cfg_.slugTracking;
    double vd = st.bendiksenVd;
    if (!st.bendiksenVdDeclared)
    {
        const double dr = std::max(cfg_.fluid.rhoL - cfg_.fluid.rhoG, 0.0);
        vd = 0.35 * sinT_[frontCell]
           * std::sqrt(kGravity * diameter_ * dr / cfg_.fluid.rhoL);
    }
    return st.bendiksenC0 * vM_ + vd;
}

std::vector<int> TwoFluidPipe::cellZones() const
{
    // E7b (5.182): the SLUG/FILM partition. A cell is SLUG (1) when its
    // centre lies inside [xR, xF] of a live tracked slug object, FILM (0)
    // otherwise - the zone boundaries ride the moving tracked faces frame
    // by frame. Diagnostic only: read-only over the tracking state, no
    // write-back to the holdup field.
    std::vector<int> z(n_, 0);
    if (!cfg_.slugTrackingOn || !cfg_.slugTracking.zoning) return z;
    double acc = 0;
    for (int i = 0; i < n_; ++i)
    {
        const double centre = acc + 0.5 * dxs_[i];
        for (const TwoFluidSlugTrack& o : slugTracks_)
            if (centre >= o.xR && centre <= o.xF) { z[i] = 1; break; }
        acc += dxs_[i];
    }
    return z;
}

std::vector<TwoFluidTerrainDip> TwoFluidPipe::terrainDips() const
{
    std::vector<TwoFluidTerrainDip> out;
    if (!cfg_.terrainOn || n_ < 3) return out;

    // elevation of every cell face above the chain inlet, the polyline the
    // dips are read from; face f sits at the upstream end of cell f, so face
    // n is the outlet. angleDeg is positive uphill, so elevation accumulates
    // +dx sinT along the chain (the same sign as centerDrop_)
    std::vector<double> z(n_ + 1, 0.0);
    std::vector<double> arc(n_ + 1, 0.0);
    for (int i = 0; i < n_; ++i)
    {
        z[i + 1] = z[i] + dxs_[i] * sinT_[i];
        arc[i + 1] = arc[i] + dxs_[i];
    }
    auto isCrest = [&](int f)
    {
        return f > 0 && f < n_ && z[f] > z[f - 1] && z[f] > z[f + 1];
    };

    const double rhoL = cfg_.fluid.rhoL;
    for (int fd = 1; fd < n_; ++fd)
    {
        // a dip is a strict local minimum of the elevation polyline
        if (!(z[fd] < z[fd - 1] && z[fd] < z[fd + 1])) continue;

        // basin bounds: nearest crest uphill on each side, else the chain end
        int up = 0;
        for (int f = fd - 1; f >= 1; --f)
            if (isCrest(f)) { up = f; break; }
        int dn = n_;
        for (int f = fd + 1; f < n_; ++f)
            if (isCrest(f)) { dn = f; break; }

        // the pool spills over the LOWER bounding crest
        const double zSpill = std::min(z[up], z[dn]);
        const double depth = zSpill - z[fd];
        if (depth < cfg_.terrain.minDepth) continue;

        double inv = 0, cap = 0, len = 0;
        for (int i = up; i < dn; ++i)
        {
            inv += alpha_[i] * dxs_[i];
            len += dxs_[i];
            // arc fraction of this cell submerged below the spill level: the
            // cell elevation runs linearly from z[i] to z[i+1]
            const double dz = z[i + 1] - z[i];
            double frac;
            if (std::abs(dz) < 1e-14)
                frac = (z[i] <= zSpill) ? 1.0 : 0.0;
            else
                frac = std::min(1.0, std::max(
                    0.0, (zSpill - std::min(z[i], z[i + 1])) / std::abs(dz)));
            cap += dxs_[i] * frac;
        }

        TwoFluidTerrainDip d;
        d.position = arc[fd];
        d.depth = depth;
        d.basinLength = len;
        d.inventory = inv;
        d.capacity = cap;
        d.fill = (cap > 1e-14) ? inv / cap : 0.0;
        d.head = rhoL * kGravity * depth;
        d.slugReady = d.fill >= 1.0;
        d.dipFace = fd;
        d.basinStart = up;
        d.basinEnd = dn;
        out.push_back(d);
    }
    return out;
}

std::vector<TwoFluidTerrainSlug> TwoFluidPipe::terrainSlugs() const
{
    return terrainSlugs_;
}

void TwoFluidPipe::terrainSlugStep(double dt)
{
    const double us = cfg_.terrain.slugVelocity;
    if (!cfg_.terrainOn || us <= 0 || n_ < 3) return;
    const double rhoL = cfg_.fluid.rhoL;

    // The terrain overlay moves liquid only: the captured lump and the
    // discharge redistribution are booked on the liquid side (inventory /
    // s.volume), and the gas already in the touched cells stays there.
    // gasMass_ is the continuity state the gas-flux step advances; this
    // overlay has no gas source to book a rescaling against, so the inventory
    // is left untouched - the same mass-closure idiom as purgeWaveStep.
    // Holding the gas density fixed here would scale gas mass by the ag
    // ratio without a booked transfer; the gas flux step rebalances the two
    // pockets afterwards, and the sealed-equilibrium witness still holds.
    auto setHoldup = [&](int i, double aNew)
    {
        alpha_[i] = aNew;
    };

    std::vector<TwoFluidTerrainDip> dips = terrainDips();
    if (terrainSlugs_.size() != dips.size())
        terrainSlugs_.assign(dips.size(), TwoFluidTerrainSlug{});

    double totalArc = 0;
    for (int i = 0; i < n_; ++i) totalArc += dxs_[i];

    for (std::size_t k = 0; k < dips.size(); ++k)
    {
        TwoFluidTerrainSlug& s = terrainSlugs_[k];
        const TwoFluidTerrainDip& d = dips[k];
        if (!s.active)
        {
            // initiation: the basin has reached its spill point (increment 1's
            // fill >= 1). The excess above the geometric spill capacity bridges
            // the pipe and is swept into a coherent lump; the basin is drawn
            // down to exactly capacity (proportional drawdown keeps its shape).
            if (d.fill >= 1.0 && d.inventory > d.capacity)
            {
                const double invBefore = d.inventory;
                const double excess = invBefore - d.capacity;
                const double factor = d.capacity / invBefore;
                for (int i = d.basinStart; i < d.basinEnd; ++i)
                    setHoldup(i, alpha_[i] * factor);
                s.active = true;
                s.volume = excess;
                s.front = d.position;
                s.capturedFrom = invBefore;
                s.velocity = us;
                if (cfg_.terrain.dynamicUs)
                    s.capPressure = cellPressure()[d.dipFace];
            }
        }
        else
        {
            double u = us;
            if (cfg_.terrain.dynamicUs)
            {
                // G2 gas-cap feedback: the cap pressure behind the lump
                // relaxes (implicit Euler) towards the EOS pressure of the
                // cell at the slug rear, and the travel velocity is the
                // energy-balance ejection value against the front's EOS
                // pressure plus the hydrostatic head still to be lifted
                // over the terrain ahead.
                std::vector<double> z(n_ + 1, 0.0);
                for (int i = 0; i < n_; ++i)
                    z[i + 1] = z[i] + dxs_[i] * sinT_[i];
                auto zAt = [&](double arcPos)
                {
                    double acc = 0;
                    for (int i = 0; i < n_; ++i)
                    {
                        acc += dxs_[i];
                        if (arcPos <= acc)
                        {
                            const double f = (arcPos - (acc - dxs_[i]))
                                           / dxs_[i];
                            return z[i] + f * (z[i + 1] - z[i]);
                        }
                    }
                    return z[n_];
                };
                const std::vector<double> pCell = cellPressure();
                double rear = s.front - s.volume;
                if (rear < 0) rear = 0;
                int ri = 0;
                {
                    double acc = 0;
                    for (int i = 0; i < n_; ++i)
                    {
                        acc += dxs_[i];
                        if (rear < acc) { ri = i; break; }
                        ri = i;
                    }
                }
                int fi = n_ - 1;
                {
                    double acc = 0;
                    for (int i = 0; i < n_; ++i)
                    {
                        acc += dxs_[i];
                        if (s.front < acc) { fi = i; break; }
                        fi = i;
                    }
                }
                const double pRear = pCell[ri];
                const double pFront = pCell[fi]
                    + rhoL * kGravity
                    * std::max(0.0, zAt(s.front) - z[d.dipFace]);
                const double tau = cfg_.terrain.gasCapTime;
                s.capPressure = (s.capPressure / dt + pRear / tau)
                              / (1.0 / dt + 1.0 / tau);
                u = std::sqrt(2.0 * std::max(s.capPressure - pFront, 0.0)
                              / rhoL);
                if (cfg_.terrain.maxSlugVelocity > 0)
                    u = std::min(u, cfg_.terrain.maxSlugVelocity);
                s.velocity = u;
                // the pre-travel position the reported us/pcap were computed
                // at (the report's "front" is this plus u*dt); independent
                // checkers re-derive the relation from this value to land
                // on the same cell and elevation the kernel sampled
                s.calcFront = s.front;
            }
            // kinematic travel at the (declared or gas-driven) velocity
            s.front += u * dt;
            if (s.front >= totalArc)
            {
                // discharge into the downstream reach's available holdup
                // capacity (1 - alpha) dx: alpha_i += V (1 - alpha_i) / avail,
                // which keeps alpha in [0, 1] and adds exactly V to
                // sum(alpha dx). V is capped at avail: an EOS-crushed
                // downstream gas pocket can offer less capacity than the
                // lump, and the remainder then holds at the outlet (the
                // lump stays active, conserving its volume) until the
                // transport frees space. With V <= avail (the incompressible
                // witness of increment 2) this is the old discharge verbatim.
                const int start = d.dipFace;
                double avail = 0;
                for (int i = start; i < n_; ++i)
                    avail += (1.0 - alpha_[i]) * dxs_[i];
                if (avail > 1e-14)
                {
                    const double V = std::min(s.volume, avail);
                    for (int i = start; i < n_; ++i)
                        setHoldup(i, alpha_[i]
                                     + V * (1.0 - alpha_[i]) / avail);
                    s.volume -= V;
                    if (s.volume <= 1e-14)
                    {
                        s.active = false;
                        s.volume = 0;
                        ++terrainDeposits_;
                    }
                    s.front = totalArc;
                }
                else
                {
                    // no available downstream capacity yet: hold the lump at
                    // the outlet (conserve its volume) until the transport
                    // frees space, rather than silently dropping the liquid
                    s.front = totalArc;
                }
            }
        }
    }
}

void TwoFluidPipe::slugTrackStep(double dt)
{
    if (!cfg_.slugTrackingOn) return;
    const TwoFluidSlugTracking& st = cfg_.slugTracking;
    const double totalArc = length_;
    const double evT = t_ + dt;

    // cell containing the arc position (a face position rides the cell to
    // its right; the outlet end rides the last cell)
    auto cellAt = [&](double x)
    {
        double acc = 0;
        for (int i = 0; i < n_; ++i)
        {
            if (x < acc + dxs_[i]) return i;
            acc += dxs_[i];
        }
        return n_ - 1;
    };
    // local liquid transport velocity at the section: the closure's phase
    // velocity vsl/alpha (the superficial flux divided by the holdup that
    // carries it); the mixture velocity where there is no liquid to price
    auto liquidVelocity = [&](int i)
    {
        const double a = alpha_[i];
        if (a <= 1e-9) return vM_;
        return liquidFlux(a, slip_[i], i) / a;
    };
    // field re-estimate of the inventory inside [xR, xF] clipped to the
    // domain; partial cells weigh in by their overlap fraction, so growth
    // and shrinkage emerge from the holdup the transport just solved
    auto inventory = [&](double xR, double xF)
    {
        const double lo = std::max(0.0, xR), hi = std::min(totalArc, xF);
        if (hi <= lo) return 0.0;
        double inv = 0, acc = 0;
        for (int i = 0; i < n_; ++i)
        {
            const double o = std::min(hi, acc + dxs_[i]) - std::max(lo, acc);
            if (o > 0) inv += alpha_[i] * o;
            acc += dxs_[i];
        }
        return inv;
    };
    // elevation of the chain at an arc position (the gasCap head term,
    // the 5.142 zAt verbatim on this pipe's own polyline)
    auto zAt = [&](double pos)
    {
        double acc = 0, z = 0;
        if (pos <= 0) return 0.0;
        for (int i = 0; i < n_; ++i)
        {
            const double z1 = z + dxs_[i] * sinT_[i];
            if (pos <= acc + dxs_[i])
                return z + (pos - acc) / dxs_[i] * (z1 - z);
            acc += dxs_[i];
            z = z1;
        }
        return z;
    };

    std::vector<double> pCell; // EOS profile, fetched once for gasCap mode

    // 1. advance: fronts at the migration velocity, rears at the local
    // liquid velocity; the front pins at the outlet while the rear drains
    // through it, the rear never passes the front (a drained slug
    // collapses onto its face and rides it out)
    for (auto& o : slugTracks_)
    {
        if (st.frontMode == 1)
        {
            // E7b (5.182): the Bendiksen closure prices the front at the
            // cell holding it; the reported uF is the effective value
            o.uF = bendiksenFrontVelocity(cellAt(o.xF));
        }
        else if (st.velocityMode == 0)
        {
            o.uF = st.slugVelocity;
        }
        else
        {
            if (pCell.empty()) pCell = cellPressure();
            const int ri = cellAt(o.xR), fi = cellAt(o.xF);
            const double tau = st.gasCapTime;
            o.capPressure = (o.capPressure / dt + pCell[ri] / tau)
                          / (1.0 / dt + 1.0 / tau);
            const double rhoL = cfg_.fluid.rhoL;
            const double pFront = pCell[fi]
                + rhoL * kGravity * std::max(0.0, zAt(o.xF) - zAt(o.xR));
            double u = std::sqrt(2.0 * std::max(o.capPressure - pFront, 0.0)
                                 / rhoL);
            if (st.maxSlugVelocity > 0) u = std::min(u, st.maxSlugVelocity);
            o.uF = u;
        }
        o.uR = liquidVelocity(cellAt(o.xR));
        o.xF += o.uF * dt;
        o.xR += o.uR * dt;
        if (o.xR < 0) o.xR = 0;
        if (o.xF > totalArc) o.xF = totalArc;
        if (o.xR > o.xF) o.xR = o.xF;
        if (o.xR < totalArc) o.inv = inventory(o.xR, o.xF);
        // discharge bookkeeping: while the span covers the outlet, the
        // liquid flux through the outlet face leaves with the slug
        if (o.xF >= totalArc && o.xR < totalArc)
            o.outVol += liquidFlux(alpha_[n_ - 1], slip_[n_ - 1], n_ - 1) * dt;
    }

    // 2. retire: the rear crossed the domain end, OUT
    for (std::size_t k = 0; k < slugTracks_.size();)
    {
        const TwoFluidSlugTrack& o = slugTracks_[k];
        if (o.xR >= totalArc)
        {
            TwoFluidSlugTrackEvent e;
            e.t = evT;
            e.id = o.id;
            e.kind = 3;
            e.xR = o.xR;
            e.xF = o.xF;
            e.inv = o.inv;
            e.uF = o.uF;
            e.uR = o.uR;
            e.outVol = o.outVol;
            e.duration = evT - o.tBorn;
            slugTrackJournal_.push_back(e);
            slugTracks_.erase(slugTracks_.begin() + k);
        }
        else
        {
            ++k;
        }
    }

    // 3. merge: sorted in flow order, the trailing front reaching the
    // leading rear coalesces the pair fully - volumes add, the leader's
    // front and the trailer's rear survive, the trailer keeps the id
    std::sort(slugTracks_.begin(), slugTracks_.end(),
              [](const TwoFluidSlugTrack& a, const TwoFluidSlugTrack& b)
              { return a.xR < b.xR; });
    for (std::size_t k = 0; k + 1 < slugTracks_.size();)
    {
        if (slugTracks_[k + 1].xR <= slugTracks_[k].xF)
        {
            TwoFluidSlugTrack& a = slugTracks_[k];
            const TwoFluidSlugTrack& b = slugTracks_[k + 1];
            TwoFluidSlugTrackEvent e;
            e.t = evT;
            e.id = a.id;
            e.kind = 2;
            e.xR = a.xR;
            e.xF = std::max(a.xF, b.xF);
            e.inv = a.inv + b.inv;
            e.uF = b.uF;
            e.uR = a.uR;
            slugTrackJournal_.push_back(e);
            a.xF = e.xF;
            a.inv = e.inv;
            a.uF = e.uF;
            a.outVol += b.outVol;
            slugTracks_.erase(slugTracks_.begin() + k + 1);
            // stay on k: the survivor may absorb the next object too
        }
        else
        {
            ++k;
        }
    }

    // 4. generate: a maximal run the scan finds that no object overlaps is
    // born as a new tracked slug, at the exact faces the statistics report
    const std::vector<TwoFluidSlugBody> runs =
        slugRuns(st.alphaThreshold, st.minLength);
    for (const TwoFluidSlugBody& r : runs)
    {
        bool tracked = false;
        for (const auto& o : slugTracks_)
            if (r.front > o.xR && r.rear < o.xF) { tracked = true; break; }
        if (tracked) continue;
        TwoFluidSlugTrack o;
        o.id = slugTrackNextId_++;
        o.xR = r.rear;
        o.xF = r.front;
        o.inv = r.inventory;
        o.tBorn = evT;
        o.uR = liquidVelocity(cellAt(o.xR));
        if (st.frontMode == 1)
        {
            // E7b (5.182): the Bendiksen value at birth (the front's cell)
            o.uF = bendiksenFrontVelocity(cellAt(o.xF));
        }
        else if (st.velocityMode == 0)
        {
            o.uF = st.slugVelocity;
        }
        else
        {
            if (pCell.empty()) pCell = cellPressure();
            o.capPressure = pCell[cellAt(o.xR)];
            o.uF = 0; // at birth the front sits on the field face; the
                      // ejection velocity is solved from the next step on
        }
        slugTracks_.push_back(o);
        TwoFluidSlugTrackEvent e;
        e.t = evT;
        e.id = o.id;
        e.kind = 1;
        e.xR = o.xR;
        e.xF = o.xF;
        e.inv = o.inv;
        e.uF = o.uF;
        e.uR = o.uR;
        slugTrackJournal_.push_back(e);
    }
}

void TwoFluidPipe::purgeWaveStep(double dt)
{
    if (!cfg_.purgeWaveOn) return;
    const TwoFluidPurgeWave& pw = cfg_.purgeWave;
    const double rhoL = cfg_.fluid.rhoL;
    const double evT = t_ + dt;
    const double start = pw.riserStart;
    const double end = pw.riserEnd > 0 ? pw.riserEnd : length_;

    // The purge front moves liquid only: the displaced liquid volume is
    // booked to outVol, and the gas already in the swept cell stays there.
    // gasMass_ is the continuity state the gas-flux step advances (see the
    // gas-mass continuity in the transient step); this overlay has no gas
    // source to book a rescaling against, so the inventory is left untouched.
    // Rescaling it by ag_new/ag_old would hold the gas density (and the EOS
    // pressure) fixed but manufacture unbooked gas mass as the gas volume
    // grows, and would suppress the cap blowdown the slug cycle is witnessed
    // on. The gas thins into the vacated volume; the flux steps rebalance it
    // on refill.

    // window cells: those whose centre lies inside [start, end]; the base
    // cell (the gas-cap pocket) is the first window cell
    int i0 = -1, i1 = -1;
    double acc = 0;
    for (int i = 0; i < n_; ++i)
    {
        const double centre = acc + 0.5 * dxs_[i];
        if (centre >= start && centre <= end)
        {
            if (i0 < 0) i0 = i;
            i1 = i;
        }
        acc += dxs_[i];
    }

    // the COMPLETE frame reports the completion, the next frame returns to
    // accumulation (the front position is cleared with the idle phase)
    if (purgeState_.phase == 3)
    {
        purgeState_.phase = 0;
        purgeState_.xWave = 0;
        purgeState_.uWave = 0;
    }

    // 1. cap relaxation: the 5.142 implicit-Euler idiom on an independent
    // state, targeting the EOS pressure of the riser-base cell
    const std::vector<double> pCell = cellPressure();
    const double pBase = pCell[i0];
    if (!purgeInit_)
    {
        purgeState_.capPressure = pBase;
        purgeInit_ = true;
    }
    const double tau = pw.gasCapTime;
    purgeState_.capPressure = (purgeState_.capPressure / dt + pBase / tau)
                            / (1.0 / dt + 1.0 / tau);
    purgeState_.basePressure = pBase;

    // 2. column head and continuity over the window
    double head = 0;
    bool continuous = true;
    for (int i = i0; i <= i1; ++i)
    {
        head += alpha_[i] * dxs_[i] * sinT_[i];
        if (alpha_[i] < pw.continuityAlpha) continuous = false;
    }
    purgeState_.head = head;

    // 3. trigger: the cap can no longer support the static head of a
    // continuous column -> a purge front is born at the riser base
    bool active = purgeState_.phase == 1 || purgeState_.phase == 2;
    if (!active && continuous && head > 0
        && (pw.maxCycles <= 0 || purgeState_.cycles < pw.maxCycles)
        && purgeState_.capPressure
           < pw.triggerFactor * rhoL * kGravity * head)
    {
        active = true;
        purgeState_.phase = 1;
        purgeState_.xWave = start;
        purgeHeadAtTrigger_ = head;
        purgeCycleOutStart_ = purgeState_.outVol;
        purgeCycleStartT_ = evT;
        TwoFluidPurgeEvent e;
        e.t = evT;
        e.kind = 1;
        e.xWave = start;
        e.head = head;
        e.capPressure = purgeState_.capPressure;
        purgeJournal_.push_back(e);
    }
    else if (purgeState_.phase == 1)
    {
        purgeState_.phase = 2;
    }

    if (!active) return;

    // 4. propagate: the front rides the declared constant, or the
    // energy-balance value of the imbalance priced with the column head at
    // the blowdown instant (the mass being driven down) and the LIVE cap
    // pressure resisting it
    double u = pw.uWave;
    if (pw.waveMode == 0)
        u = std::sqrt(2.0 * std::max(pw.triggerFactor * rhoL * kGravity
                                         * purgeHeadAtTrigger_
                                     - purgeState_.capPressure,
                                     0.0) / rhoL);
    purgeState_.uWave = u;
    purgeState_.xWave += u * dt;
    if (purgeState_.xWave > end) purgeState_.xWave = end;

    // 5. purge: every window cell the front has passed (cell centre below
    // the front) relaxes to the residual holdup - the wall film value the
    // blowing gas sustains, held in both directions against the transport
    // drain. The net displaced liquid leaves through the riser top into
    // the out accumulator (the flow pulse)
    acc = 0;
    for (int i = i0; i <= i1; ++i)
    {
        const double centre = acc + 0.5 * dxs_[i];
        if (centre <= purgeState_.xWave && alpha_[i] != pw.alphaResidual)
        {
            purgeState_.outVol += (alpha_[i] - pw.alphaResidual) * dxs_[i];
            alpha_[i] = pw.alphaResidual;
        }
        acc += dxs_[i];
    }

    // 6. complete: the front crossed the riser top, the column is blown
    // dry, the cycle counter grows and accumulation resumes
    if (purgeState_.xWave >= end)
    {
        purgeState_.phase = 3;
        ++purgeState_.cycles;
        TwoFluidPurgeEvent e;
        e.t = evT;
        e.kind = 2;
        e.xWave = purgeState_.xWave;
        e.head = head;
        e.capPressure = purgeState_.capPressure;
        e.cycle = purgeState_.cycles;
        e.outVol = purgeState_.outVol - purgeCycleOutStart_;
        e.duration = evT - purgeCycleStartT_;
        e.period = purgeLastCompleteT_ > 0 ? evT - purgeLastCompleteT_ : 0.0;
        purgeLastCompleteT_ = evT;
        purgeJournal_.push_back(e);
    }
}

void TwoFluidPipe::piggingStep(double dt)
{
    if (!cfg_.piggingOn || !pigState_.launched || pigState_.arrived) return;
    const TwoFluidPigging& pg = cfg_.pigging;
    const double end = pg.endX > 0 ? pg.endX : length_;
    const double evT = t_ + dt;

    // swept-liquid inventory over [x0, x1] of the step's final holdup field;
    // partial cells weigh in by their overlap fraction, exactly the idiom
    // slugTrackStep uses. H is a diagnosed volume (per m2 of bore), not a
    // write-back to alpha_.
    auto swept = [&](double x0, double x1)
    {
        const double lo = std::max(0.0, x0), hi = std::min(length_, x1);
        if (hi <= lo) return 0.0;
        double inv = 0, acc = 0;
        for (int i = 0; i < n_; ++i)
        {
            const double o = std::min(hi, acc + dxs_[i]) - std::max(lo, acc);
            if (o > 0) inv += alpha_[i] * o;
            acc += dxs_[i];
        }
        return inv;
    };

    // a growing front chokes the pig: v = u0 (1 - beta clamp(H/C, 0, 1)).
    const double choke = pg.frontCapacity > 0
        ? std::min(1.0, std::max(0.0, pigState_.holdupFront
                                 / pg.frontCapacity)) : 0.0;
    const double v = pg.pigVelocity * (1.0 - pg.blockage * choke);
    const double x0 = pigState_.x;
    double x1 = x0 + v * dt;
    if (x1 > end) x1 = end;

    pigState_.holdupFront += pg.pushEfficiency * swept(x0, x1);
    pigState_.x = x1;
    pigState_.v = v;

    if (x1 >= end)
    {
        // arrival: the whole piled front is released; freeze the pig (the
        // state line thereafter reads v = 0, event NONE) and journal the
        // ARRIVE with the velocity that carried it to the end arc.
        pigState_.outVol = pigState_.holdupFront;
        pigState_.arrived = true;
        pigState_.v = 0.0;
        TwoFluidPigEvent e;
        e.t = evT;
        e.kind = 1;
        e.x = x1;
        e.v = v;
        e.holdupFront = pigState_.holdupFront;
        e.outVol = pigState_.outVol;
        pigJournal_.push_back(e);
    }
}

// E7c (5.183): one explicit staggered step of the linear acoustic state.
// The pair (p', u') solves the linearized Euler acoustics of a quiescent
// duct (LeVeque ch. 9): dp'/dt = -rho a^2 du'/dx, du'/dt = -(1/rho) dp'/dx.
// Discretization: u lives on faces, p in cells; the u rows update FIRST on
// the current p, then the p rows update on the NEW u - the symplectic
// (leapfrog) pairing that is neutrally stable for a dt/dx <= 1, unlike
// forward Euler on both rows which amplifies every mode. The open end is
// pressure-release (ghost p' = 0, its face updates on the half control
// volume next to the boundary), the closed end pins u' = 0 (a rate or
// sealed boundary imposes its velocity, so acoustically it is a wall).
// No write-back: the mean state (momentum, holdup, gas mass) is untouched,
// the acoustic layer is superposed reporting only.
void TwoFluidPipe::acousticsStep(double dt)
{
    if (!cfg_.acousticsOn) return;
    const TwoFluidAcoustics& ac = cfg_.acoustics;

    // per-cell wave speed and inertia. Default twoPhase path: Wood's
    // homogeneous-mixture sound speed (A. B. Wood, "A Textbook of Sound",
    // 1930; the standard uniform-equilibrium compressibility rule),
    //   1 / (rho_m a^2) = alphaL / (rhoL aL^2) + alphaG / (rhoG aG^2),
    // with aG^2 = gasR T the EOS-exact gas derivative and rho_m the cell's
    // mixture density (the acoustic rows themselves carry rho_m as inertia,
    // not rhoG). aL is the DECLARED liquid sound speed; absent one the liquid
    // is taken incompressible (aL -> inf, zero compliance) - the kernel
    // invents no liquid EOS. A pure-gas cell therefore gives exactly
    // sqrt(gasR T), so an all-gas line keeps the original EOS wave speed.
    // Applicability boundary: Wood assumes uniform PRESSURE equilibrium - both
    // phases compress under the same dp with no interphase mass/heat exchange
    // and no slip, so it prices the mixture's effective bulk modulus; it is
    // the right wave speed for a bubbly/dispersed homogeneous mixture and it
    // degenerates correctly (alphaG -> 0 gives aL, alphaL -> 0 gives aG).
    // The legacy pure-gas-everywhere path is still selectable as aModel eos;
    // aModel wood is the same mixture law but requires a declared aL; a
    // declared soundSpeed replaces any model wholesale.
    std::vector<double> a2(n_), rho(n_);
    double aMax = 0.0, dxMin = 1e300;
    for (int i = 0; i < n_; ++i)
    {
        double a;
        if (ac.soundSpeed > 0) a = ac.soundSpeed;
        else if (ac.aModel == "wood")
        {
            const double gasShare = 1.0 - alpha_[i];
            const double aG2 = cfg_.fluid.gasR * sectionTemperature(i);
            const double rhoG = cellGasDensity(i);
            // the compliance sum; the degenerate single-phase cells
            // drop the vanishing phase's term by its alpha weight
            const double compliance
                = alpha_[i] / (cfg_.fluid.rhoL * ac.liquidSoundSpeed
                               * ac.liquidSoundSpeed)
                + (gasShare > 0 ? gasShare / (rhoG * aG2) : 0.0);
            a = std::sqrt(1.0 / (compliance * mixtureRhoCell(i)));
        }
        else if (ac.aModel == "twoPhase")
        {
            // Same Wood compliance, but the liquid speed is optional: absent
            // a declaration the liquid is incompressible (zero compliance),
            // which needs no property table. Pure gas then degenerates to the
            // EOS derivative exactly; two-phase cells get the mixture speed a
            // depressurization wave actually rides, not the pure-gas speed.
            const double gasShare = 1.0 - alpha_[i];
            const double aG2 = cfg_.fluid.gasR * sectionTemperature(i);
            const double rhoG = cellGasDensity(i);
            if (!(gasShare > 0))
            {
                // liquid-only cell: the incompressible limit has no finite
                // speed, so the declared aL is the only honest answer
                if (!(ac.liquidSoundSpeed > 0))
                    cfFatal("twoFluid: acoustics twoPhase reached a"
                            " liquid-only cell (alpha = 1): declare"
                            " liquidSoundSpeed for a finite liquid wave"
                            " speed");
                a = ac.liquidSoundSpeed;
            }
            else
            {
                const double liquidCompliance = ac.liquidSoundSpeed > 0
                    ? alpha_[i] / (cfg_.fluid.rhoL * ac.liquidSoundSpeed
                                   * ac.liquidSoundSpeed)
                    : 0.0;
                const double compliance = liquidCompliance + gasShare / (rhoG * aG2);
                a = std::sqrt(1.0 / (compliance * mixtureRhoCell(i)));
            }
        }
        else a = std::sqrt(cfg_.fluid.gasR * sectionTemperature(i)); // eos legacy
        a2[i] = a * a;
        rho[i] = mixtureRhoCell(i);
        if (a > aMax) aMax = a;
        if (dxs_[i] < dxMin) dxMin = dxs_[i];
    }
    // the stability bound is explicit and loud: a violation is a refusal,
    // not a silent wobble (there is no adaptive dt on this row - the app
    // drives fixed dt, so the case must declare one that fits)
    const double cflNow = aMax * dt / dxMin;
    if (cflNow > ac.cfl)
    {
        std::ostringstream os;
        os << std::setprecision(6) << cflNow;
        cfFatal("twoFluid: acoustic CFL violated - a_max dt / dx_min = "
                + os.str() + " exceeds the declared bound "
                + std::to_string(ac.cfl) + " - reduce dt");
    }

    std::vector<double> uNew(n_ + 1, 0.0);
    // interior faces: centred momentum over the two half cells
    for (int j = 1; j < n_; ++j)
    {
        const double rf = 0.5 * (rho[j - 1] + rho[j]);
        const double df = 0.5 * (dxs_[j - 1] + dxs_[j]);
        uNew[j] = uA_[j] - dt * (pA_[j] - pA_[j - 1]) / (rf * df);
    }
    if (ac.openAtInlet)
    {
        // inlet vents: ghost p' = 0 left of face 0 (half control volume)
        uNew[0] = uA_[0]
                - dt * (0.0 - pA_[0]) / (rho[0] * 0.5 * dxs_[0]);
        uNew[n_] = 0.0; // outlet end is the hard wall
    }
    else
    {
        // outlet vents: ghost p' = 0 right of face n (half control volume)
        uNew[n_] = uA_[n_]
                - dt * (0.0 - pA_[n_ - 1]) / (rho[n_ - 1] * 0.5 * dxs_[n_ - 1]);
        uNew[0] = 0.0; // inlet end is the hard wall
    }
    // cell pressures on the NEW face velocities (the symplectic pairing)
    std::vector<double> pNew(n_, 0.0);
    for (int i = 0; i < n_; ++i)
        pNew[i] = pA_[i]
                - rho[i] * a2[i] * dt / dxs_[i] * (uNew[i + 1] - uNew[i]);
    uA_.swap(uNew);
    pA_.swap(pNew);
}

namespace
{

TwoFluidBoundary loadTwoFluidBoundary(const Dictionary& d)
{
    TwoFluidBoundary b;
    const std::string type = d.getWord("type");
    if (type == "sealed") b.type = TwoFluidBoundary::Type::sealed;
    else if (type == "phaseRates")
        b.type = TwoFluidBoundary::Type::phaseRates;
    else if (type == "fixedPressure")
        b.type = TwoFluidBoundary::Type::fixedPressure;
    else cfFatal("twoFluid: unknown boundary type '" + type + "'");
    b.vsl = d.getScalarOr("vsl", 0.0);
    b.vsg = d.getScalarOr("vsg", 0.0);
    if (d.has("schedule"))
    {
        const Dictionary& sd = d.getDict("schedule");
        for (const std::string& name : sd.keys())
        {
            const Dictionary& pd = sd.getDict(name);
            TwoFluidSchedulePoint pt;
            pt.t = pd.getScalar("t");
            // missing keys hold the previous rate (the scalar rates
            // before the first point)
            pt.vsl = pd.getScalarOr("vsl", b.vsl);
            pt.vsg = pd.getScalarOr("vsg", b.vsg);
            b.schedule.push_back(pt);
        }
        std::sort(b.schedule.begin(), b.schedule.end(),
                  [](const TwoFluidSchedulePoint& x,
                     const TwoFluidSchedulePoint& y) { return x.t < y.t; });
        for (std::size_t k = 0; k < b.schedule.size(); ++k)
        {
            if (!(b.schedule[k].t > 0))
                cfFatal("twoFluid: schedule points need t > 0 (the scalar"
                        " rates govern at t = 0)");
            if (k && !(b.schedule[k - 1].t < b.schedule[k].t))
                cfFatal("twoFluid: schedule points need distinct times");
        }
    }
    b.pressure = d.getScalarOr("p", 0.0);
    b.alpha = d.getScalarOr("alpha", 0.0);
    b.waterFraction = d.getScalarOr("waterFraction", 0.0);
    b.hydrateFraction = d.getScalarOr("hydrateFraction", 0.0);
    b.waxDepositFraction = d.getScalarOr("waxDepositFraction", 0.0);
    b.entrainedFraction = d.getScalarOr("entrainedFraction", 0.0);
    b.inhibitorFraction = d.getScalarOr("inhibitorFraction", 0.0);
    b.sandFraction = d.getScalarOr("sandFraction", 0.0);
    if (d.has("zInlet")) b.zInlet = d.getScalarList("zInlet"); // G1
    b.temperature = d.getScalarOr("temperature", 0.0);
    return b;
}

// the physics blocks shared by the serial case and the branch network:
// fluid / closure / initial / thermal / hydrate
void loadTwoFluidCommon(const Dictionary& caseDict, TwoFluidCase& c)
{
    if (!caseDict.has("fluid")) cfFatal("twoFluid: case has no fluid{} block");
    const Dictionary& fd = caseDict.getDict("fluid");
    c.fluid.rhoL = fd.getScalar("rhoL");
    c.fluid.rhoG = fd.getScalar("rhoG");
    c.fluid.muL = fd.getScalar("muL");
    c.fluid.muG = fd.getScalar("muG");
    c.fluid.cL = fd.getScalarOr("cL", 0.0);
    c.fluid.cG = fd.getScalarOr("cG", 0.0);
    c.fluid.gasR = fd.getScalarOr("gasR", 0.0);
    // 5.75 / 5.77 (A5): the liquid phase may declare a power-law rheology
    // and/or a yield stress. The declarations that would be silently
    // ignored are refused instead -- the worst kind of silent number is
    // one the file seems to ask for. A zero declaration is the undeclared
    // law verbatim (the held set is empty and the wall keeps its old
    // Churchill expression bit for bit).
    c.fluid.flowIndex = fd.getScalarOr("flowIndex", 0.0);
    c.fluid.consistency = fd.getScalarOr("consistency", 0.0);
    c.fluid.yieldStress = fd.getScalarOr("yieldStress", 0.0);
    if (c.fluid.flowIndex < 0)
        cfFatal("twoFluid: fluid has a negative flowIndex (0 means a"
                " Newtonian liquid)");
    if (c.fluid.flowIndex > 0 && !(c.fluid.consistency > 0))
        cfFatal("twoFluid: fluid declares a flowIndex without a positive"
                " consistency: at n != 1 the fluid's muL is not the same"
                " quantity as K (Pa.s^n)");
    if (c.fluid.flowIndex == 0 && fd.has("consistency"))
        cfFatal("twoFluid: fluid declares a consistency without a"
                " flowIndex: without n the consistency is not a quantity"
                " the wall law reads");
    if (c.fluid.yieldStress < 0)
        cfFatal("twoFluid: fluid has a negative yieldStress: the static"
                " wall cannot hold less than zero (a negative threshold"
                " would make the held set empty)");

    if (!caseDict.has("closure"))
        cfFatal("twoFluid: case has no closure{} block");
    const Dictionary& cd = caseDict.getDict("closure");
    const std::string model = cd.getWordOr("model", "wallis");
    if (model == "wallis") c.closure.model = 0;
    else if (model == "zuber") c.closure.model = 1;
    else if (model == "dynamic") c.closure.model = 2;
    else if (model == "full") c.closure.model = 3;
    else cfFatal("twoFluid: unknown closure model '" + model + "'");
    c.closure.slipVelocity = cd.getScalarOr("slipVelocity", 0.0);
    c.closure.c0 = cd.getScalarOr("C0", 1.0);
    c.closure.driftVelocity = cd.getScalarOr("driftVelocity", 0.0);
    c.closure.interfacialDrag = cd.getScalarOr("interfacialDrag", 0.0);
    // shut-in-safe algebraic variant (opt-in): gravity-projected drift +
    // conservative receive-capacity thermal assembly. Default false keeps
    // the historical scheme bit-for-bit (all standing checkers).
    c.closure.gravitySlip = cd.has("gravitySlip") && cd.getBool("gravitySlip");
    // the projection is meaningful only where an inclination exists; on a
    // horizontal algebraic-closure line it degenerates the drift to zero,
    // which is exactly the shut-in-safe intent, so no refusal - but a
    // gravitySlip declaration alongside a dynamic/full closure is a
    // category error: those closures already drive slip from buoyancy
    if (c.closure.gravitySlip && (c.closure.model == 2 || c.closure.model == 3))
        cfFatal("twoFluid: closure gravitySlip belongs to the algebraic"
                " closures (wallis/zuber); the dynamic/full closures"
                " already drive the slip from buoyancy");

    if (!caseDict.has("initial"))
        cfFatal("twoFluid: case has no initial{} block");
    const Dictionary& id = caseDict.getDict("initial");
    c.alpha0 = id.getScalarOr("alpha", 0.0);
    c.water0 = id.getScalarOr("water", 0.0);
    c.temperature0 = id.getScalarOr("temperature", 300.0);
    // 5.77: optional starting velocity of a pressure-driven pipe
    c.vm0 = id.getScalarOr("vm", 0.0);
    if (c.alpha0 < 0 || c.alpha0 > 1 || c.water0 < 0 || c.water0 > 1)
        cfFatal("twoFluid: initial alpha/water must be in [0,1]");
    if (id.has("fills"))
    {
        const Dictionary& fills = id.getDict("fills");
        for (const std::string& name : fills.keys())
        {
            const Dictionary& f = fills.getDict(name);
            TwoFluidFill fill;
            fill.from = f.getScalar("from");
            fill.to = f.getScalar("to");
            fill.alpha = f.getScalar("alpha");
            c.fills.push_back(fill);
        }
    }

    if (caseDict.has("thermal"))
    {
        c.thermalOn = true;
        const Dictionary& th = caseDict.getDict("thermal");
        c.thermal.wallHeatTransfer = th.getScalar("wallHeatTransfer");
        const std::string wm = th.getWordOr("wallModel", "constant");
        if (wm == "constant") c.thermal.wallModel = 0;
        else if (wm == "geothermal") c.thermal.wallModel = 1;
        else cfFatal("twoFluid: unknown thermal wallModel '" + wm + "'");
        c.thermal.wallTemperature = th.getScalarOr("wallTemperature", 0.0);
        c.thermal.surfaceTemperature
            = th.getScalarOr("surfaceTemperature", 0.0);
        c.thermal.geothermalGradient
            = th.getScalarOr("gradient", 0.0);
        c.thermal.inletDepth = th.getScalarOr("inletDepth", 0.0);
        if (c.thermal.wallHeatTransfer < 0)
            cfFatal("twoFluid: negative wallHeatTransfer");
    }

    if (caseDict.has("formation"))
    {
        c.formationOn = true;
        const Dictionary& fm = caseDict.getDict("formation");
        c.formation.conductivity = fm.getScalar("conductivity");
        c.formation.outerRadius = fm.getScalar("outerRadius");
        c.formation.rhoC = fm.getScalar("rhoC");
        c.formation.shells = static_cast<int>(fm.getScalar("shells"));
        const std::string si = fm.getWordOr("semiInfinite", "no");
        if (si == "yes") c.formation.semiInfinite = true;
        else if (si != "no")
            cfFatal("twoFluid: formation semiInfinite must be 'yes' or 'no'"
                    " (got '" + si + "')");
    }

    if (caseDict.has("hydrate"))
    {
        c.hydrateOn = true;
        const Dictionary& hd = caseDict.getDict("hydrate");
        c.hydrate.curveA = hd.getScalar("curveA");
        c.hydrate.curveB = hd.getScalar("curveB");
        c.hydrate.growthRate = hd.getScalar("growthRate");
        c.hydrate.dissociationRate = hd.getScalar("dissociationRate");
        c.hydrate.gasContent = hd.getScalar("gasContent");
        c.hydrate.latentHeat = hd.getScalar("latentHeat");
        c.hydrate.initial = hd.getScalarOr("initial", 0.0);
    }

    if (caseDict.has("wax"))
    {
        c.waxOn = true;
        const Dictionary& wd = caseDict.getDict("wax");
        // Tf / dHf default to the n-C24H50 NIST constants; the rates and
        // wax fraction are case declarations (no built-in rate table).
        c.wax.meltingTemperature
            = wd.getScalarOr("meltingTemperature", 323.75);
        c.wax.fusionEnthalpy = wd.getScalarOr("fusionEnthalpy", 54890.0);
        c.wax.waxFraction = wd.getScalarOr("waxFraction", 0.0);
        c.wax.growthRate = wd.getScalarOr("growthRate", 0.0);
        c.wax.dissolutionRate = wd.getScalarOr("dissolutionRate", 0.0);
        c.wax.latentHeat = wd.getScalarOr("latentHeat", 0.0);
        c.wax.initial = wd.getScalarOr("initial", 0.0);
        // 5.98: both keys are absent by default -> the old hydraulic path.
        c.wax.feedback = wd.has("feedback") && wd.getBool("feedback");
        c.wax.depositRoughness = wd.getScalarOr("depositRoughness", 0.0);
    }

    if (caseDict.has("asphaltene"))
    {
        // 5.99: the AOP / width / fraction / rates are all case declarations
        // (an experimental onset, not a kernel correlation).
        c.asphalteneOn = true;
        const Dictionary& ad = caseDict.getDict("asphaltene");
        c.asphaltene.onsetPressure = ad.getScalar("onsetPressure");
        c.asphaltene.onsetWidth = ad.getScalarOr("onsetWidth", 0.0);
        c.asphaltene.asphalteneFraction
            = ad.getScalarOr("asphalteneFraction", 0.0);
        c.asphaltene.growthRate = ad.getScalarOr("growthRate", 0.0);
        c.asphaltene.dissolutionRate
            = ad.getScalarOr("dissolutionRate", 0.0);
        c.asphaltene.initial = ad.getScalarOr("initial", 0.0);
    }

    if (caseDict.has("scale"))
    {
        // 5.99: brine chemistry is declared; Ksp defaults to calcite at
        // 25 C (log Ksp = -8.480, Plummer & Busenberg 1982).
        c.scaleOn = true;
        const Dictionary& sd = caseDict.getDict("scale");
        c.scale.calciumActivity = sd.getScalarOr("calciumActivity", 0.0);
        c.scale.carbonateActivity
            = sd.getScalarOr("carbonateActivity", 0.0);
        c.scale.solubilityProduct
            = sd.getScalarOr("solubilityProduct", 3.311e-9);
        c.scale.scaleFraction = sd.getScalarOr("scaleFraction", 0.0);
        c.scale.growthRate = sd.getScalarOr("growthRate", 0.0);
        c.scale.dissolutionRate = sd.getScalarOr("dissolutionRate", 0.0);
        c.scale.initial = sd.getScalarOr("initial", 0.0);
    }

    if (caseDict.has("entrainment"))
    {
        c.entrainmentOn = true;
        const Dictionary& ed = caseDict.getDict("entrainment");
        c.entrainment.vCrit = ed.getScalarOr("vCrit", 0.0);
        c.entrainment.vRef = ed.getScalarOr("vRef", 1.0);
        c.entrainment.relaxTime = ed.getScalarOr("relaxTime", 0.0);
        c.entrainment.initial = ed.getScalarOr("initial", 0.0);
    }

    if (caseDict.has("emulsion"))
    {
        c.emulsionOn = true;
        const Dictionary& ed = caseDict.getDict("emulsion");
        const std::string m = ed.getWord("model");
        if (m == "brinkman") c.emulsion.model = 1;
        else if (m == "richardson") c.emulsion.model = 2;
        else cfFatal("twoFluid: unknown emulsion model '" + m + "'");
        c.emulsion.muWater = ed.getScalar("muWater");
    }

    if (caseDict.has("inhibitor"))
    {
        c.inhibitorOn = true;
        const Dictionary& id = caseDict.getDict("inhibitor");
        c.inhibitor.kDepression = id.getScalar("kDepression");
        c.inhibitor.initial = id.getScalarOr("initial", 0.0);
    }

    // O4 (5.163): rigorous inhibitor phase equilibrium (opt-in). The
    // activity constant is required whenever the block is declared (also
    // when disabled, so a re-enabled file cannot silently carry a typo);
    // enable false keeps the Hammerschmidt path of inhibitor{}.
    if (caseDict.has("inhibitorEquilibrium"))
    {
        const Dictionary& ied = caseDict.getDict("inhibitorEquilibrium");
        c.inhibitorEquilibriumOn = true;
        if (ied.has("enable"))
        {
            const std::string w = ied.getWord("enable");
            if (w == "true") c.inhibitorEquilibriumOn = true;
            else if (w == "false") c.inhibitorEquilibriumOn = false;
            else cfFatal("twoFluid: inhibitorEquilibrium enable must be"
                         " true or false");
        }
        c.inhibitorEquilibrium.activityK = ied.getScalar("activityK");
        if (ied.has("basis"))
        {
            const std::string b = ied.getWord("basis");
            if (b == "liquid") c.inhibitorEquilibrium.basisWater = false;
            else if (b == "water") c.inhibitorEquilibrium.basisWater = true;
            else cfFatal("twoFluid: unknown inhibitorEquilibrium basis '"
                         + b + "' (liquid or water)");
        }
    }

    if (caseDict.has("sand"))
    {
        c.sandOn = true;
        const Dictionary& sd = caseDict.getDict("sand");
        c.sand.terminalVelocity = sd.getScalar("terminalVelocity");
        c.sand.criticalVelocity = sd.getScalar("criticalVelocity");
        c.sand.erodeRate = sd.getScalarOr("erodeRate", 0.0);
        c.sand.initial = sd.getScalarOr("initial", 0.0);
        c.sand.bedInitial = sd.getScalarOr("bedInitial", 0.0);
        if (sd.has("feedback"))
        {
            const std::string w = sd.getWord("feedback");
            if (w == "true") c.sand.feedback = true;
            else if (w == "false") c.sand.feedback = false;
            else cfFatal("twoFluid: sand feedback must be true or false");
        }
    }

    if (caseDict.has("composition"))
    {
        // G1: named feed components, each with its mole fraction. Loud on
        // unknown names (requireComponent), count outside 2..8, negative
        // fractions and sums further than 1e-8 from 1; within that
        // tolerance the feed is auto-normalised to exactly 1 (documented
        // in TwoFluid.h).
        c.compositionOn = true;
        const Dictionary& cd = caseDict.getDict("composition");
        double sum = 0;
        for (const std::string& key : cd.keys())
        {
            const Dictionary::Node& node = cd.find(key);
            if (node.kind != Dictionary::Kind::scalar)
                cfFatal("twoFluid: composition entry '" + key + "' must be"
                        " a mole fraction number");
            pvt::requireComponent(key); // unknown names are loud
            if (!(node.scalar >= 0))
                cfFatal("twoFluid: composition entry '" + key + "' has a"
                        " negative mole fraction");
            c.compositionNames.push_back(key);
            c.compositionZ0.push_back(node.scalar);
            sum += node.scalar;
        }
        if (c.compositionNames.size() < 2)
            cfFatal("twoFluid: composition needs at least 2 components");
        if (c.compositionNames.size() > 8)
            cfFatal("twoFluid: composition supports at most 8 components");
        if (std::fabs(sum - 1.0) > 1e-8)
            cfFatal("twoFluid: composition mole fractions must sum to 1"
                    " within 1e-8");
        for (double& z : c.compositionZ0) z /= sum;
    }

    if (caseDict.has("slug"))
    {
        c.slugOn = true;
        const Dictionary& sd = caseDict.getDict("slug");
        c.slug.alphaThreshold = sd.getScalarOr("alphaThreshold", 0.5);
        c.slug.minLength = sd.getScalarOr("minLength", 0.0);
        if (c.slug.alphaThreshold <= 0 || c.slug.alphaThreshold >= 1)
            cfFatal("twoFluid: slug alphaThreshold must be in (0, 1)");
        if (c.slug.minLength < 0)
            cfFatal("twoFluid: negative slug minLength");
    }

    if (caseDict.has("terrain"))
    {
        c.terrainOn = true;
        const Dictionary& td = caseDict.getDict("terrain");
        c.terrain.minDepth = td.getScalarOr("minDepth", 0.0);
        if (c.terrain.minDepth < 0)
            cfFatal("twoFluid: negative terrain minDepth");
        c.terrain.slugVelocity = td.getScalarOr("slugVelocity", 0.0);
        if (c.terrain.slugVelocity < 0)
            cfFatal("twoFluid: negative terrain slugVelocity");
        if (td.has("dynamicUs"))
        {
            const std::string w = td.getWord("dynamicUs");
            if (w == "true") c.terrain.dynamicUs = true;
            else if (w == "false") c.terrain.dynamicUs = false;
            else cfFatal("twoFluid: terrain dynamicUs must be true or false");
        }
        c.terrain.gasCapTime = td.getScalarOr("gasCapTime", 0.0);
        c.terrain.maxSlugVelocity = td.getScalarOr("maxSlugVelocity", 0.0);
        if (c.terrain.dynamicUs)
        {
            if (c.terrain.slugVelocity <= 0)
                cfFatal("twoFluid: terrain dynamicUs needs a positive"
                        " slugVelocity gate (the cycle switch)");
            if (!(c.terrain.gasCapTime > 0))
                cfFatal("twoFluid: terrain dynamicUs needs a positive"
                        " gasCapTime");
            if (c.terrain.maxSlugVelocity < 0)
                cfFatal("twoFluid: negative terrain maxSlugVelocity");
            if (!(c.fluid.gasR > 0))
                cfFatal("twoFluid: terrain dynamicUs reads the gas-cap EOS"
                        " pressure but no gasR is declared");
        }
    }

    if (caseDict.has("slugTracking"))
    {
        // O2 (5.161): Lagrangian slug tracking. The block is opt-in and
        // enable false (or no block at all) leaves the run byte-identical.
        const Dictionary& sd = caseDict.getDict("slugTracking");
        bool enable = true;
        if (sd.has("enable"))
        {
            const std::string w = sd.getWord("enable");
            if (w == "true") enable = true;
            else if (w == "false") enable = false;
            else cfFatal("twoFluid: slugTracking enable must be true or false");
        }
        c.slugTracking.alphaThreshold = sd.getScalarOr("alphaThreshold", 0.5);
        c.slugTracking.minLength = sd.getScalarOr("minLength", 0.0);
        if (c.slugTracking.alphaThreshold <= 0
            || c.slugTracking.alphaThreshold >= 1)
            cfFatal("twoFluid: slugTracking alphaThreshold must be in (0, 1)");
        if (c.slugTracking.minLength < 0)
            cfFatal("twoFluid: negative slugTracking minLength");
        const std::string mode = sd.getWordOr("velocityMode", "constant");
        if (mode == "constant") c.slugTracking.velocityMode = 0;
        else if (mode == "gasCap") c.slugTracking.velocityMode = 1;
        else cfFatal("twoFluid: unknown slugTracking velocityMode '" + mode
                     + "' (constant or gasCap)");
        // E7b (5.182): opt-in Bendiksen front closure. Declaring
        // frontClosure bendiksen replaces the declared velocityMode law
        // with v_front = C0 v_m + V_d (Bendiksen 1984, turbulent branch;
        // source and applicability in TwoFluidSlugTracking). The two laws
        // are mutually exclusive: any legacy front-velocity declaration
        // alongside bendiksen is rejected loudly, so a run either rides a
        // fully declared migration law or the closure, never a mixture.
        if (sd.has("frontClosure"))
        {
            const std::string fc = sd.getWord("frontClosure");
            if (fc == "bendiksen")
            {
                if (sd.has("velocityMode"))
                    cfFatal("twoFluid: slugTracking frontClosure bendiksen"
                            " is mutually exclusive with velocityMode (the"
                            " Bendiksen closure replaces the declared"
                            " migration law)");
                if (sd.has("slugVelocity"))
                    cfFatal("twoFluid: slugTracking frontClosure bendiksen"
                            " is mutually exclusive with slugVelocity");
                if (sd.has("gasCapTime"))
                    cfFatal("twoFluid: slugTracking frontClosure bendiksen"
                            " is mutually exclusive with gasCapTime (no"
                            " gas-cap state rides a Bendiksen front)");
                if (sd.has("maxSlugVelocity"))
                    cfFatal("twoFluid: slugTracking frontClosure bendiksen"
                            " is mutually exclusive with maxSlugVelocity"
                            " (the clamp belongs to the gasCap law)");
                c.slugTracking.bendiksenC0 = sd.getScalarOr("C0", 1.2);
                if (!(c.slugTracking.bendiksenC0 > 0))
                    cfFatal("twoFluid: slugTracking Bendiksen C0 must be"
                            " positive");
                if (sd.has("driftVelocity"))
                {
                    c.slugTracking.bendiksenVd
                        = sd.getScalarOr("driftVelocity", 0.0);
                    c.slugTracking.bendiksenVdDeclared = true;
                }
                if (!(c.fluid.rhoL > c.fluid.rhoG))
                    cfFatal("twoFluid: slugTracking Bendiksen drift needs"
                            " rhoL > rhoG");
                c.slugTracking.frontMode = 1;
            }
            else
                cfFatal("twoFluid: unknown slugTracking frontClosure '"
                        + fc + "' (bendiksen is the only closure this"
                               " build ships)");
        }
        c.slugTracking.slugVelocity = sd.getScalarOr("slugVelocity", 0.0);
        if (c.slugTracking.slugVelocity < 0)
            cfFatal("twoFluid: negative slugTracking slugVelocity");
        c.slugTracking.gasCapTime = sd.getScalarOr("gasCapTime", 0.0);
        c.slugTracking.maxSlugVelocity = sd.getScalarOr("maxSlugVelocity", 0.0);
        if (c.slugTracking.frontMode == 1)
        {
            // the Bendiksen closure owns the front; the legacy modes'
            // required declarations do not apply (and were rejected above
            // if present)
        }
        else if (c.slugTracking.velocityMode == 0)
        {
            if (!(c.slugTracking.slugVelocity > 0))
                cfFatal("twoFluid: slugTracking constant mode needs a"
                        " positive slugVelocity");
        }
        else
        {
            if (!(c.slugTracking.gasCapTime > 0))
                cfFatal("twoFluid: slugTracking gasCap mode needs a positive"
                        " gasCapTime (the cap-pressure relaxation)");
            if (c.slugTracking.maxSlugVelocity < 0)
                cfFatal("twoFluid: negative slugTracking maxSlugVelocity");
            if (!(c.fluid.gasR > 0))
                cfFatal("twoFluid: slugTracking gasCap mode reads the gas-cap"
                        " EOS pressure but no gasR is declared");
        }
        // E7b (5.182): opt-in SLUG/FILM cell partition, reported as a
        // `zone` key on every cell row. Purely diagnostic over the
        // tracking state; enable false keeps the whole block (zoning
        // included) inert and the run byte-identical.
        if (sd.has("zoning"))
        {
            const std::string zw = sd.getWord("zoning");
            if (zw == "true") c.slugTracking.zoning = true;
            else if (zw == "false") c.slugTracking.zoning = false;
            else cfFatal("twoFluid: slugTracking zoning must be true or"
                         " false");
        }
        c.slugTrackingOn = enable;
    }

    if (caseDict.has("purgeWave"))
    {
        // O3 (5.162): the distributed riser purge (blowdown) wave. The
        // block is opt-in and enable false (or no block at all) leaves the
        // run byte-identical.
        const Dictionary& pd = caseDict.getDict("purgeWave");
        bool enable = true;
        if (pd.has("enable"))
        {
            const std::string w = pd.getWord("enable");
            if (w == "true") enable = true;
            else if (w == "false") enable = false;
            else cfFatal("twoFluid: purgeWave enable must be true or false");
        }
        c.purgeWave.triggerFactor = pd.getScalarOr("triggerFactor", 1.0);
        if (!(c.purgeWave.triggerFactor > 0))
            cfFatal("twoFluid: purgeWave triggerFactor must be positive");
        c.purgeWave.continuityAlpha = pd.getScalarOr("continuityAlpha", 0.5);
        if (c.purgeWave.continuityAlpha <= 0
            || c.purgeWave.continuityAlpha >= 1)
            cfFatal("twoFluid: purgeWave continuityAlpha must be in (0, 1)");
        c.purgeWave.alphaResidual = pd.getScalarOr("alphaResidual", 0.05);
        if (c.purgeWave.alphaResidual < 0 || c.purgeWave.alphaResidual >= 1)
            cfFatal("twoFluid: purgeWave alphaResidual must be in [0, 1)");
        if (c.purgeWave.alphaResidual >= c.purgeWave.continuityAlpha)
            cfFatal("twoFluid: purgeWave alphaResidual must stay below"
                    " continuityAlpha (a blown-dry column must not"
                    " re-trigger as continuous)");
        const std::string mode = pd.getWordOr("waveMode", "closed");
        if (mode == "closed") c.purgeWave.waveMode = 0;
        else if (mode == "constant") c.purgeWave.waveMode = 1;
        else cfFatal("twoFluid: unknown purgeWave waveMode '" + mode
                     + "' (closed or constant)");
        c.purgeWave.uWave = pd.getScalarOr("uWave", 0.0);
        if (c.purgeWave.uWave < 0)
            cfFatal("twoFluid: negative purgeWave uWave");
        if (c.purgeWave.waveMode == 1 && !(c.purgeWave.uWave > 0))
            cfFatal("twoFluid: purgeWave constant mode needs a positive"
                    " uWave");
        c.purgeWave.gasCapTime = pd.getScalarOr("gasCapTime", 0.0);
        if (!(c.purgeWave.gasCapTime > 0))
            cfFatal("twoFluid: purgeWave needs a positive gasCapTime"
                    " (the cap-pressure relaxation)");
        if (!(c.fluid.gasR > 0))
            cfFatal("twoFluid: purgeWave reads the gas-cap EOS pressure"
                    " but no gasR is declared");
        c.purgeWave.maxCycles =
            static_cast<int>(pd.getScalarOr("maxCycles", 0.0));
        if (c.purgeWave.maxCycles < 0)
            cfFatal("twoFluid: negative purgeWave maxCycles");
        c.purgeWave.riserStart = pd.getScalarOr("riserStart", 0.0);
        c.purgeWave.riserEnd = pd.getScalarOr("riserEnd", 0.0);
        if (c.purgeWave.riserStart < 0)
            cfFatal("twoFluid: negative purgeWave riserStart");
        if (c.purgeWave.riserEnd < 0)
            cfFatal("twoFluid: negative purgeWave riserEnd");
        if (c.purgeWave.riserEnd > 0
            && c.purgeWave.riserEnd <= c.purgeWave.riserStart)
            cfFatal("twoFluid: purgeWave riserEnd must exceed riserStart");
        c.purgeWaveOn = enable;
    }

    if (caseDict.has("pigging"))
    {
        // E6 (5.176): reduced 1D pig tracking. The block is opt-in and
        // enable false (or no block at all) leaves the run byte-identical.
        const Dictionary& gd = caseDict.getDict("pigging");
        bool enable = true;
        if (gd.has("enable"))
        {
            const std::string w = gd.getWord("enable");
            if (w == "true") enable = true;
            else if (w == "false") enable = false;
            else cfFatal("twoFluid: pigging enable must be true or false");
        }
        c.pigging.startX = gd.getScalarOr("startX", 0.0);
        c.pigging.endX = gd.getScalarOr("endX", 0.0);
        c.pigging.pigVelocity = gd.getScalarOr("pigVelocity", 0.0);
        c.pigging.blockage = gd.getScalarOr("blockage", 0.0);
        c.pigging.pushEfficiency = gd.getScalarOr("pushEfficiency", 0.0);
        c.pigging.frontCapacity = gd.getScalarOr("frontCapacity", 0.0);
        if (c.pigging.startX < 0)
            cfFatal("twoFluid: negative pigging startX");
        if (c.pigging.endX < 0)
            cfFatal("twoFluid: negative pigging endX");
        if (c.pigging.endX > 0 && c.pigging.endX <= c.pigging.startX)
            cfFatal("twoFluid: pigging endX must exceed startX");
        if (c.pigging.pigVelocity < 0)
            cfFatal("twoFluid: negative pigging pigVelocity");
        if (c.pigging.blockage < 0 || c.pigging.blockage >= 1.0)
            cfFatal("twoFluid: pigging blockage must be in [0, 1)");
        if (c.pigging.pushEfficiency < 0 || c.pigging.pushEfficiency > 1.0)
            cfFatal("twoFluid: pigging pushEfficiency must be in [0, 1]");
        if (!(c.pigging.frontCapacity > 0))
            cfFatal("twoFluid: pigging needs a positive frontCapacity");
        if (enable && !(c.pigging.pigVelocity > 0))
            cfFatal("twoFluid: pigging needs a positive pigVelocity");
        c.piggingOn = enable;
    }

    if (caseDict.has("acoustics"))
    {
        // E7c (5.183): linear along-the-line pressure-wave acoustics for
        // the EOS mode. The block is opt-in and enable false (or no block
        // at all) leaves the run byte-identical. Model, source and the
        // honest applicability boundary live in the TwoFluid.h class
        // documentation; the declared values are validated here, loudly.
        const Dictionary& ad = caseDict.getDict("acoustics");
        bool enable = true;
        if (ad.has("enable"))
        {
            const std::string w = ad.getWord("enable");
            if (w == "true") enable = true;
            else if (w == "false") enable = false;
            else cfFatal("twoFluid: acoustics enable must be true or"
                         " false");
        }
        c.acoustics.initialGauge = ad.getScalarOr("initialGauge", 0.0);
        if (c.acoustics.initialGauge < 0)
            cfFatal("twoFluid: acoustics initialGauge must be positive"
                    " (a negative acoustic overpressure is not a"
                    " depressurization)");
        c.acoustics.soundSpeed = ad.getScalarOr("soundSpeed", 0.0);
        if (c.acoustics.soundSpeed < 0)
            cfFatal("twoFluid: acoustics soundSpeed must be positive"
                    " (0 selects the model path of aModel: the default"
                    " twoPhase Wood mixture speed, the legacy pure-gas"
                    " eos speed, or Wood's full mixture speed)");
        // wave-speed model: the default twoPhase rides Wood's mixture law
        // with an incompressible liquid when aL is not declared; "eos" is
        // the legacy pure-gas-everywhere path, "wood" the full mixture law
        // with a required declared aL. A declared word is read loudly;
        // unknown words are refused, not ignored.
        if (ad.has("aModel"))
        {
            const std::string w = ad.getWord("aModel");
            if (w == "eos") c.acoustics.aModel = "eos";
            else if (w == "wood") c.acoustics.aModel = "wood";
            else if (w == "twoPhase") c.acoustics.aModel = "twoPhase";
            else cfFatal("twoFluid: unknown acoustics aModel '" + w
                         + "' (twoPhase, eos or wood)");
        }
        if (ad.has("liquidSoundSpeed"))
        {
            c.acoustics.liquidSoundSpeed
                = ad.getScalarOr("liquidSoundSpeed", 0.0);
            if (!(c.acoustics.liquidSoundSpeed > 0))
                cfFatal("twoFluid: acoustics liquidSoundSpeed must be"
                        " positive (it is the LIQUID phase sound speed"
                        " aL of Wood's equation; the kernel invents no"
                        " liquid EOS to derive it)");
        }
        // the declared wholesale override and the full Wood mixture model
        // are mutually exclusive: soundSpeed replaces the wave speed
        // wholesale, aModel wood derives it - declaring both is a
        // contradiction, and a wood declaration without aL has nothing
        // to derive from (the liquid sound speed is not invented). The
        // default twoPhase model needs no aL (incompressible liquid) and
        // a declared soundSpeed simply overrides it.
        if (c.acoustics.aModel == "wood")
        {
            if (c.acoustics.soundSpeed > 0)
                cfFatal("twoFluid: acoustics soundSpeed and aModel wood"
                        " are mutually exclusive (a declared speed"
                        " replaces the wave speed wholesale; wood"
                        " derives it from the phase speeds)");
            if (!(c.acoustics.liquidSoundSpeed > 0))
                cfFatal("twoFluid: aModel wood needs a positive"
                        " liquidSoundSpeed (Wood's equation needs the"
                        " LIQUID phase sound speed aL; the kernel"
                        " invents no liquid EOS)");
        }
        c.acoustics.cfl = ad.getScalarOr("cfl", 0.5);
        if (c.acoustics.cfl <= 0 || c.acoustics.cfl > 0.5)
            cfFatal("twoFluid: acoustics cfl must be in (0, 0.5] (the"
                    " pressure-release end updates its face on a half"
                    " control volume, which halves the explicit"
                    " limit)");
        if (ad.has("openEnd"))
        {
            const std::string oe = ad.getWord("openEnd");
            if (oe == "outlet") c.acoustics.openAtInlet = false;
            else if (oe == "inlet") c.acoustics.openAtInlet = true;
            else cfFatal("twoFluid: unknown acoustics openEnd '" + oe
                         + "' (outlet or inlet)");
        }
        if (enable)
        {
            // the wave speed is the density derivative of the kernel's
            // own gas law; without the EOS there is no a to take
            if (!(c.fluid.gasR > 0))
                cfFatal("twoFluid: acoustics needs the gas EOS (fluid"
                        " gasR > 0): the wave speed is the EOS density"
                        " derivative");
            if (!(c.acoustics.initialGauge > 0))
                cfFatal("twoFluid: enabled acoustics needs a positive"
                        " initialGauge (the linear field vents from"
                        " it)");
        }
        c.acousticsOn = enable;
    }
    // ---- 5.187 (F-class): the transient risk report -----------------------
    // Three read-only channels, the transient counterpart of the steady P22
    // `riskReport{}` (5.159). Every key is read loudly; nothing here writes a
    // model field, so declaring the block cannot move the hydraulics.
    if (caseDict.has("riskReport"))
    {
        const Dictionary& rr = caseDict.getDict("riskReport");
        static const char* kRiskKeys[] = {
            "erosionC", "hydrate", "corrosion",
            "hydrateLowK", "hydrateMediumK"};
        for (const std::string& k : rr.keys())
        {
            bool known = false;
            for (const char* w : kRiskKeys)
                if (k == w) { known = true; break; }
            if (!known)
                cfFatal("twoFluid: riskReport declares unknown key '" + k
                        + "': the report reads erosionC, the hydrate line,"
                          " the sour-gas fractions and their bands");
        }
        c.riskOn = true;
        c.risk.declared = true;
        if (rr.has("erosionC"))
        {
            c.risk.erosionC = rr.getScalar("erosionC");
            if (!(c.risk.erosionC > 0))
                cfFatal("twoFluid: riskReport erosionC must be positive: it"
                        " divides the API RP 14E velocity limit");
        }
        c.risk.hydrateLowK = rr.getScalarOr("hydrateLowK", 3.0);
        c.risk.hydrateMediumK = rr.getScalarOr("hydrateMediumK", 8.0);
        if (!(c.risk.hydrateLowK > 0)
            || !(c.risk.hydrateMediumK > c.risk.hydrateLowK))
            cfFatal("twoFluid: riskReport hydrate bands must be positive and"
                    " strictly ascending (low < medium)");
        if (rr.has("hydrate"))
        {
            const Dictionary& hb = rr.getDict("hydrate");
            static const char* kHydKeys[] = {"Tref", "pRef", "dH"};
            for (const std::string& k : hb.keys())
            {
                bool known = false;
                for (const char* w : kHydKeys)
                    if (k == w) { known = true; break; }
                if (!known)
                    cfFatal("twoFluid: riskReport hydrate declares unknown"
                            " key '" + k + "' (Tref, pRef, dH)");
            }
            c.risk.hydrateTref = hb.getScalar("Tref");
            c.risk.hydratePRef = hb.getScalar("pRef");
            c.risk.hydrateDH = hb.getScalar("dH");
            if (!(c.risk.hydrateTref > 0) || !(c.risk.hydratePRef > 0)
                || !(c.risk.hydrateDH > 0))
                cfFatal("twoFluid: riskReport hydrate Tref, pRef and dH must"
                        " all be positive (one dissociation curve)");
            if (!(c.fluid.gasR > 0))
                cfFatal("twoFluid: riskReport hydrate prices the dissociation"
                        " line at each section's own EOS pressure: declare the"
                        " gas EOS (fluid gasR > 0)");
            c.risk.hydrateOn = true;
        }
        if (rr.has("corrosion"))
        {
            const Dictionary& cb = rr.getDict("corrosion");
            static const char* kCorKeys[] = {"co2", "h2s", "co2Low", "co2High"};
            for (const std::string& k : cb.keys())
            {
                bool known = false;
                for (const char* w : kCorKeys)
                    if (k == w) { known = true; break; }
                if (!known)
                    cfFatal("twoFluid: riskReport corrosion declares unknown"
                            " key '" + k + "' (co2, h2s, co2Low, co2High)");
            }
            c.risk.co2Fraction = cb.getScalarOr("co2", 0.0);
            c.risk.h2sFraction = cb.getScalarOr("h2s", 0.0);
            if (!(c.risk.co2Fraction >= 0) || !(c.risk.co2Fraction <= 1)
                || !(c.risk.h2sFraction >= 0) || !(c.risk.h2sFraction <= 1))
                cfFatal("twoFluid: riskReport corrosion co2/h2s are gas-phase"
                        " mole fractions and must be in [0,1]");
            c.risk.co2LowPa = cb.getScalarOr("co2Low", 1e4);
            c.risk.co2HighPa = cb.getScalarOr("co2High", 1e5);
            if (!(c.risk.co2LowPa > 0)
                || !(c.risk.co2HighPa > c.risk.co2LowPa))
                cfFatal("twoFluid: riskReport corrosion CO2 bands must be"
                        " positive and strictly ascending (co2Low < co2High)");
            if (!(c.fluid.gasR > 0))
                cfFatal("twoFluid: riskReport corrosion prices the CO2/H2S"
                        " partial pressures at the highest section pressure:"
                        " declare the gas EOS (fluid gasR > 0)");
            c.risk.corrosionOn = true;
        }
    }

    // ---- 5.219 (F-class): hydrate blockage --------------------------------
    // The opt-in counterpart of the 5.186 non-blocking reading: a declared
    // blockage{} charges h to the bore (layer) and to the line (bridging
    // plug). Undeclared, no field moves.
    if (caseDict.has("blockage"))
    {
        const Dictionary& bk = caseDict.getDict("blockage");
        static const char* kBlockageKeys[] = {
            "packingFraction", "depositionRate", "erosionRate"};
        for (const std::string& k : bk.keys())
        {
            bool known = false;
            for (const char* w : kBlockageKeys)
                if (k == w) { known = true; break; }
            if (!known)
                cfFatal("twoFluid: blockage declares unknown key '" + k
                        + "' (packingFraction, depositionRate, erosionRate)");
        }
        if (!c.hydrateOn)
            cfFatal("twoFluid: blockage{} charges the hydrate fraction to the"
                    " bore, but no hydrate{} is declared: the blockage would"
                    " act on nothing");
        c.blockageOn = true;
        c.blockage.packingFraction = bk.getScalarOr("packingFraction", 0.6);
        c.blockage.depositionRate = bk.getScalarOr("depositionRate", 0.0);
        c.blockage.erosionRate = bk.getScalarOr("erosionRate", 0.0);
        if (!(c.blockage.packingFraction > 0)
            || c.blockage.packingFraction > 1.0)
            cfFatal("twoFluid: blockage packingFraction must be in (0, 1]:"
                    " it is the bridging limit on the hydrate fraction");
        if (c.blockage.depositionRate < 0)
            cfFatal("twoFluid: negative blockage depositionRate: the layer"
                    " grows by deposition, stripping is the erosionRate");
        if (c.blockage.erosionRate < 0)
            cfFatal("twoFluid: negative blockage erosionRate");
    }

    // ---- 5.220 (F-class): erosion / corrosion hydraulic feedback ----------
    // The opt-in counterpart of the 5.187 read-only risk report: the risk
    // moves the bore. Undeclared, no diameter term reads a velocity.
    if (caseDict.has("hydraulicFeedback"))
    {
        const Dictionary& hf = caseDict.getDict("hydraulicFeedback");
        static const char* kFeedbackKeys[] = {
            "erosionC", "erosionRate", "h2s", "corrosionRate"};
        for (const std::string& k : hf.keys())
        {
            bool known = false;
            for (const char* w : kFeedbackKeys)
                if (k == w) { known = true; break; }
            if (!known)
                cfFatal("twoFluid: hydraulicFeedback declares unknown key '"
                        + k + "' (erosionC, erosionRate, h2s, corrosionRate)");
        }
        c.hydraulicFeedbackOn = true;
        c.hydraulicFeedback.erosionC = hf.getScalarOr("erosionC", 100.0);
        c.hydraulicFeedback.erosionRate
            = hf.getScalarOr("erosionRate", 0.0);
        c.hydraulicFeedback.h2sFraction = hf.getScalarOr("h2s", 0.0);
        c.hydraulicFeedback.corrosionRate
            = hf.getScalarOr("corrosionRate", 0.0);
        if (!(c.hydraulicFeedback.erosionC > 0))
            cfFatal("twoFluid: hydraulicFeedback erosionC must be positive:"
                    " it divides the API RP 14E velocity limit");
        if (c.hydraulicFeedback.erosionRate < 0)
            cfFatal("twoFluid: negative hydraulicFeedback erosionRate: wear"
                    " grows the bore, it cannot shrink it");
        if (c.hydraulicFeedback.corrosionRate < 0)
            cfFatal("twoFluid: negative hydraulicFeedback corrosionRate");
        if (!(c.hydraulicFeedback.h2sFraction >= 0)
            || !(c.hydraulicFeedback.h2sFraction <= 1))
            cfFatal("twoFluid: hydraulicFeedback h2s is a gas-phase mole"
                    " fraction and must be in [0,1]");
        if (c.hydraulicFeedback.corrosionRate > 0 && !(c.fluid.gasR > 0))
            cfFatal("twoFluid: hydraulicFeedback corrosion prices the H2S"
                    " partial pressure at the highest section EOS pressure:"
                    " declare the gas EOS (fluid gasR > 0)");
    }

    // ---- 5.222 (F-class): time-varying ambient temperature ----------------
    // The wall's ambient level follows a declared piecewise-linear schedule.
    // Undeclared, wallTemperature(i) is the constant/geothermal expression.
    if (caseDict.has("ambientSchedule"))
    {
        const Dictionary& as = caseDict.getDict("ambientSchedule");
        static const char* kScheduleKeys[] = {"times", "ambientT"};
        for (const std::string& k : as.keys())
        {
            bool known = false;
            for (const char* w : kScheduleKeys)
                if (k == w) { known = true; break; }
            if (!known)
                cfFatal("twoFluid: ambientSchedule declares unknown key '" + k
                        + "' (times, ambientT)");
        }
        if (!c.thermalOn)
            cfFatal("twoFluid: ambientSchedule{} moves the thermal wall, but"
                    " no thermal{} is declared: the schedule would feed"
                    " nothing");
        c.ambientScheduleOn = true;
        c.ambientSchedule.times = as.getScalarList("times");
        c.ambientSchedule.ambientT = as.getScalarList("ambientT");
        if (c.ambientSchedule.times.size()
            != c.ambientSchedule.ambientT.size())
            cfFatal("twoFluid: ambientSchedule times and ambientT need the"
                    " same length");
        if (c.ambientSchedule.times.size() < 2)
            cfFatal("twoFluid: ambientSchedule needs at least two points"
                    " (times, ambientT)");
        for (std::size_t k = 0; k < c.ambientSchedule.times.size(); ++k)
        {
            if (c.ambientSchedule.times[k] < 0)
                cfFatal("twoFluid: ambientSchedule times must be non-negative");
            if (k && !(c.ambientSchedule.times[k - 1]
                       < c.ambientSchedule.times[k]))
                cfFatal("twoFluid: ambientSchedule times must be strictly"
                        " ascending");
        }
    }
}

} // namespace

TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict)
{
    TwoFluidCase c;
    loadTwoFluidCommon(caseDict, c);

    if (caseDict.has("pipe"))
    {
        const Dictionary& pd = caseDict.getDict("pipe");
        TwoFluidGeo g;
        g.name = pd.getWordOr("name", "pipe");
        g.length = pd.getScalar("length");
        g.diameter = pd.getScalar("diameter");
        if (pd.has("dInner"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares dInner: the "
                    "annulus section (5.73) is a steady-network declaration");
        if (pd.has("ductHeight"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares ductHeight: the "
                    "rectangular duct (5.80) is a steady-network declaration");
        if (pd.has("eccentricity"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares eccentricity: the "
                    "eccentric annulus (5.87) is a steady-network declaration");
        if (pd.has("flowPattern"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares flowPattern: the "
                    "writable flow-pattern dispatch (5.88) overrides the "
                    "Beggs-Brill map of the steady black-oil network; the "
                    "two-fluid row classifies its own regimes");
        if (pd.has("liftDomePressure") || pd.has("liftAreaRatio")
            || pd.has("liftSpringPressure") || pd.has("liftDischargeCoeff")
            || pd.has("liftKappa"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares a lift valve: the "
                    "gas-lift valve (5.89) is a steady-network injection-string "
                    "element priced on node pressures; the two-fluid row moves "
                    "its own phases and has no bellows to open");
        if (pd.has("rodPlungerDiameter") || pd.has("rodStrokeLength")
            || pd.has("rodStrokesPerMinute") || pd.has("rodDiameter")
            || pd.has("rodLength") || pd.has("rodYoungModulus")
            || pd.has("rodDensity"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares a rod pump: the "
                    "subsurface rod pump (5.90) is a steady-network element "
                    "whose rate is its stroke times its barrel; the two-fluid "
                    "row moves its own phases through the bore and has no "
                    "polished rod to pull a card against");
        g.roughness = pd.getScalarOr("roughness", 0.0);
        g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
        g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
        c.pipes.push_back(g);
    }
    else if (caseDict.has("pipes"))
    {
        // named sub-blocks in first-occurrence order form the serial chain
        const Dictionary& ps = caseDict.getDict("pipes");
        for (const std::string& k : ps.keys())
        {
            const Dictionary& pd = ps.getDict(k);
            TwoFluidGeo g;
            g.name = k;
            g.length = pd.getScalar("length");
            g.diameter = pd.getScalar("diameter");
            g.roughness = pd.getScalarOr("roughness", 0.0);
            g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
            g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
            c.pipes.push_back(g);
        }
    }
    else
    {
        cfFatal("twoFluid: case has no pipe{} or pipes{} block");
    }

    if (!caseDict.has("boundaries"))
        cfFatal("twoFluid: case has no boundaries{} block");
    const Dictionary& bd = caseDict.getDict("boundaries");
    if (!bd.has("inlet") || !bd.has("outlet"))
        cfFatal("twoFluid: boundaries{} needs inlet{} and outlet{}");
    c.inlet = loadTwoFluidBoundary(bd.getDict("inlet"));
    c.outlet = loadTwoFluidBoundary(bd.getDict("outlet"));

    if (c.outlet.schedule.size() > 0)
        cfFatal("twoFluid: rate schedules are only supported on the inlet");
    if (c.inlet.schedule.size() > 0
        && c.inlet.type != TwoFluidBoundary::Type::phaseRates)
        cfFatal("twoFluid: a rate schedule needs a phaseRates inlet");

    return c;
}

TwoFluidNetworkCase loadTwoFluidNetworkCase(const Dictionary& caseDict)
{
    // E7c (5.183): the acoustic layer is a serial-row feature - the
    // network report has no pA/aunode columns yet, so a network
    // declaration would run unreported physics. Loud, not silent.
    if (caseDict.has("acoustics"))
        cfFatal("twoFluid: acoustics is a serial-row feature; the branch"
                " network row carries no acoustic report columns yet");
    // 5.187: the risk report is a serial-row feature for the same reason - the
    // network frame has no `risk` columns, and a declaration that runs
    // unreported physics is refused rather than swallowed.
    if (caseDict.has("riskReport"))
        cfFatal("twoFluid: riskReport is a serial-row feature; the branch"
                " network row carries no risk report columns yet");
    // 5.219/5.220/5.222: same rule for the three opt-in blockage, wear and
    // ambient-schedule blocks - the network report carries no blk/wallT
    // columns, so a declaration there would run unreported physics.
    if (caseDict.has("blockage"))
        cfFatal("twoFluid: blockage is a serial-row feature; the branch"
                " network row carries no blockage report columns yet");
    if (caseDict.has("hydraulicFeedback"))
        cfFatal("twoFluid: hydraulicFeedback is a serial-row feature; the"
                " branch network row carries no wear report columns yet");
    if (caseDict.has("ambientSchedule"))
        cfFatal("twoFluid: ambientSchedule is a serial-row feature; the"
                " branch network row carries no wallT report columns yet");
    TwoFluidNetworkCase nc;
    nc.perPipe.networkMode = true;
    loadTwoFluidCommon(caseDict, nc.perPipe);

    if (!caseDict.has("network"))
        cfFatal("twoFluid: network case has no network{} block");
    const Dictionary& nd = caseDict.getDict("network");
    if (!nd.has("pipes")) cfFatal("twoFluid: network{} needs pipes{}");
    const Dictionary& ps = nd.getDict("pipes");
    for (const std::string& k : ps.keys())
    {
        const Dictionary& pd = ps.getDict(k);
        TwoFluidGeo g;
        g.name = k;
        g.length = pd.getScalar("length");
        g.diameter = pd.getScalar("diameter");
        if (pd.has("dInner"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares dInner: the "
                    "annulus section (5.73) is a steady-network declaration");
        if (pd.has("ductHeight"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares ductHeight: the "
                    "rectangular duct (5.80) is a steady-network declaration");
        if (pd.has("eccentricity"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares eccentricity: the "
                    "eccentric annulus (5.87) is a steady-network declaration");
        if (pd.has("flowPattern"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares flowPattern: the "
                    "writable flow-pattern dispatch (5.88) overrides the "
                    "Beggs-Brill map of the steady black-oil network; the "
                    "two-fluid row classifies its own regimes");
        if (pd.has("liftDomePressure") || pd.has("liftAreaRatio")
            || pd.has("liftSpringPressure") || pd.has("liftDischargeCoeff")
            || pd.has("liftKappa"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares a lift valve: the "
                    "gas-lift valve (5.89) is a steady-network injection-string "
                    "element priced on node pressures; the two-fluid row moves "
                    "its own phases and has no bellows to open");
        if (pd.has("rodPlungerDiameter") || pd.has("rodStrokeLength")
            || pd.has("rodStrokesPerMinute") || pd.has("rodDiameter")
            || pd.has("rodLength") || pd.has("rodYoungModulus")
            || pd.has("rodDensity"))
            cfFatal("twoFluid: pipe '" + g.name + "' declares a rod pump: the "
                    "subsurface rod pump (5.90) is a steady-network element "
                    "whose rate is its stroke times its barrel; the two-fluid "
                    "row moves its own phases through the bore and has no "
                    "polished rod to pull a card against");
        g.roughness = pd.getScalarOr("roughness", 0.0);
        g.angleDeg = pd.getScalarOr("angleDeg", 0.0);
        g.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));
        nc.pipes.push_back(g);
    }
    if (!nd.has("junctions")) cfFatal("twoFluid: network{} needs junctions{}");
    const Dictionary& js = nd.getDict("junctions");
    for (const std::string& k : js.keys())
    {
        const Dictionary& jd = js.getDict(k);
        TwoFluidJunctionSpec s;
        s.name = k;
        s.ends = jd.getWordList("ends");
        if (s.ends.size() < 2)
            cfFatal("twoFluid: junction '" + k + "' needs at least two ends");
        nc.junctions.push_back(s);
    }
    if (nd.has("boundaries"))
    {
        const Dictionary& bd = nd.getDict("boundaries");
        for (const std::string& k : bd.keys())
        {
            if (k.rfind('.') == std::string::npos)
                cfFatal("twoFluid: network boundary key '" + k
                        + "' must be pipe.inlet or pipe.outlet");
            nc.boundaries[k] = loadTwoFluidBoundary(bd.getDict(k));
        }
    }
    for (const auto& entry : nc.boundaries)
        if (entry.second.schedule.size() > 0)
            cfFatal("twoFluid: rate schedules do not support network"
                    " junctions (boundary '" + entry.first + "')");

    // topology validation: every pipe joins one or two junctions (a span
    // carries none), that/those ends boundary-free, a single-junction
    // pipe's other end carrying a boundary; a rate-driven far end must
    // inject toward the junction; the dynamic/full closures do not
    // support junctions yet
    if (nc.perPipe.closure.model == 2 || nc.perPipe.closure.model == 3)
        cfFatal("twoFluid: the dynamic/full closures do not support"
                " junctions yet");
    std::map<std::string, char> junctionEnd; // "pipe.inlet" -> i/o
    for (const TwoFluidJunctionSpec& j : nc.junctions)
    {
        for (const std::string& ref : j.ends)
        {
            const auto dot = ref.rfind('.');
            if (dot == std::string::npos)
                cfFatal("twoFluid: junction end '" + ref
                        + "' must be pipe.inlet or pipe.outlet");
            const std::string p = ref.substr(0, dot);
            const std::string e = ref.substr(dot + 1);
            if (e != "inlet" && e != "outlet")
                cfFatal("twoFluid: junction end '" + ref + "' must end in"
                        " .inlet or .outlet");
            bool found = false;
            for (const TwoFluidGeo& g : nc.pipes) found |= g.name == p;
            if (!found)
                cfFatal("twoFluid: junction '" + j.name + "' references"
                        " unknown pipe '" + p + "'");
            const std::string key = p + "." + e;
            if (junctionEnd.count(key))
                cfFatal("twoFluid: pipe end '" + key + "' joins two"
                        " junctions");
            junctionEnd[key] = e == "outlet" ? 'o' : 'i';
        }
    }
    for (const TwoFluidGeo& g : nc.pipes)
    {
        const bool inJ = junctionEnd.count(g.name + ".inlet") > 0;
        const bool outJ = junctionEnd.count(g.name + ".outlet") > 0;
        if (!inJ && !outJ)
            cfFatal("twoFluid: pipe '" + g.name
                    + "' needs one or two junction ends");
        if (inJ && outJ) continue; // span between two junctions: both ends
                                   // boundary-free, driven by pJ at each
        const std::string farKey = g.name + (inJ ? ".outlet" : ".inlet");
        if (!nc.boundaries.count(farKey))
            cfFatal("twoFluid: pipe '" + g.name + "' end '" + farKey
                    + "' has no boundary and no junction");
        const std::string jKey = g.name + (inJ ? ".inlet" : ".outlet");
        if (nc.boundaries.count(jKey))
            cfFatal("twoFluid: junction end '" + jKey
                    + "' must not also carry a boundary");
        const TwoFluidBoundary& far = nc.boundaries.at(farKey);
        if (far.type == TwoFluidBoundary::Type::phaseRates)
        {
            if (far.vsl + far.vsg < 0)
                cfFatal("twoFluid: phaseRates injection with negative total"
                        " rate is not supported");
            if (!outJ)
                cfFatal("twoFluid: pipe '" + g.name + "': a rate-driven far"
                        " end needs the junction at its outlet");
        }
        if (nc.perPipe.thermalOn
            && far.type != TwoFluidBoundary::Type::sealed
            && !(far.temperature > 0))
            cfFatal("twoFluid: pipe '" + g.name
                    + "': thermal on needs a far-end temperature");
    }

    return nc;
}

TwoFluidNetwork::TwoFluidNetwork(const TwoFluidNetworkCase& c)
    : cfg_(c)
{
    for (int i = 0; i < static_cast<int>(cfg_.pipes.size()); ++i)
        pipeIndex_[cfg_.pipes[i].name] = i;

    for (const TwoFluidJunctionSpec& js : cfg_.junctions)
    {
        Junction jn;
        jn.name = js.name;
        jn.alphaJ = cfg_.perPipe.alpha0;
        jn.waterJ = cfg_.perPipe.water0;
        jn.hydrateJ = cfg_.perPipe.hydrateOn ? cfg_.perPipe.hydrate.initial
                                             : 0.0;
        jn.waxJ = cfg_.perPipe.waxOn ? cfg_.perPipe.wax.initial : 0.0;
        jn.entrJ = cfg_.perPipe.entrainmentOn
                       ? cfg_.perPipe.entrainment.initial
                       : 0.0;
        jn.inhibJ = cfg_.perPipe.inhibitorOn
                        ? cfg_.perPipe.inhibitor.initial
                        : 0.0;
        jn.sandJ = cfg_.perPipe.sandOn ? cfg_.perPipe.sand.initial : 0.0;
        jn.zJ = cfg_.perPipe.compositionOn ? cfg_.perPipe.compositionZ0
                                           : std::vector<double>(); // G1
        jn.tJ = cfg_.perPipe.temperature0;
        jn.rhoJ = cfg_.perPipe.fluid.rhoG;
        bool pressureRef = false;
        for (const std::string& ref : js.ends)
        {
            const auto dot = ref.rfind('.');
            EndRef e;
            e.pipe = pipeIndex_.at(ref.substr(0, dot));
            e.atOutlet = ref.substr(dot + 1) == "outlet";
            const std::string farKey = cfg_.pipes[e.pipe].name
                                     + (e.atOutlet ? ".inlet" : ".outlet");
            // a span pipe's far end is another junction end: no boundary,
            // keep the default-constructed (sealed) EndRef.far - the warm
            // start skips it and the per-step global solve drives the span
            // from its two junction pressures
            const auto farIt = cfg_.boundaries.find(farKey);
            if (farIt != cfg_.boundaries.end()) e.far = farIt->second;
            pressureRef |= e.far.type
                        == TwoFluidBoundary::Type::fixedPressure;
            jn.ends.push_back(e);
        }
        if (!pressureRef)
            cfFatal("twoFluid: junction '" + jn.name + "' has no pressure"
                    " reference (needs a fixed-pressure far end)");
        // initial junction pressure: the steady conductance-weighted mean
        // of the far fixed pressures (the t = 0 report and the Newton
        // warm start)
        const TwoFluidFluid& f = cfg_.perPipe.fluid;
        const double a0 = cfg_.perPipe.alpha0;
        const double muBar = a0 * f.muL + (1.0 - a0) * f.muG;
        double num = 0, den = 0;
        for (const EndRef& e : jn.ends)
        {
            if (e.far.type != TwoFluidBoundary::Type::fixedPressure)
                continue;
            const TwoFluidGeo& g = cfg_.pipes[e.pipe];
            const double cond = M_PI * g.diameter * g.diameter / 4.0
                              * g.diameter * g.diameter
                              / (32.0 * muBar * g.length);
            num += cond * e.far.pressure;
            den += cond;
        }
        jn.p = num / den;
        junctions_.push_back(jn);
    }

    for (int i = 0; i < static_cast<int>(cfg_.pipes.size()); ++i)
    {
        const TwoFluidGeo& g = cfg_.pipes[i];
        TwoFluidCase pc = cfg_.perPipe;
        pc.pipes.assign(1, g);
        // the junction end gets a placeholder ghost (fixed pressure at the
        // initial state) that the per-step network solve overwrites
        TwoFluidBoundary junctionGhost;
        junctionGhost.type = TwoFluidBoundary::Type::fixedPressure;
        junctionGhost.pressure = 0;
        junctionGhost.alpha = pc.alpha0;
        junctionGhost.waterFraction = pc.water0;
        junctionGhost.hydrateFraction
            = pc.hydrateOn ? pc.hydrate.initial : 0.0;
        junctionGhost.waxDepositFraction
            = pc.waxOn ? pc.wax.initial : 0.0;
        junctionGhost.entrainedFraction
            = pc.entrainmentOn ? pc.entrainment.initial : 0.0;
        junctionGhost.inhibitorFraction
            = pc.inhibitorOn ? pc.inhibitor.initial : 0.0;
        junctionGhost.sandFraction
            = pc.sandOn ? pc.sand.initial : 0.0;
        junctionGhost.zInlet = pc.compositionOn ? pc.compositionZ0
                                                : std::vector<double>();
        junctionGhost.temperature = pc.temperature0;
        bool jIn = false, jOut = false;
        for (const Junction& jn : junctions_)
            for (const EndRef& e : jn.ends)
                if (e.pipe == i) (e.atOutlet ? jOut : jIn) = true;
        // junction ends get the placeholder ghost the per-step network
        // solve overwrites; the loader guarantees every non-junction end
        // carries a boundary
        pc.inlet = jIn ? junctionGhost
                       : cfg_.boundaries.at(g.name + ".inlet");
        pc.outlet = jOut ? junctionGhost
                         : cfg_.boundaries.at(g.name + ".outlet");
        pipes_.push_back(std::make_unique<TwoFluidPipe>(pc));
    }
}

void TwoFluidNetwork::step(double dt)
{
    if (!(dt > 0)) cfFatal("twoFluid: dt must be positive");
    solveNetwork(dt);
    for (auto& p : pipes_) p->step(dt);
    t_ += dt;
}

void TwoFluidNetwork::solveNetwork(double dt)
{
    const int np = static_cast<int>(pipes_.size());
    const int nj = static_cast<int>(junctions_.size());

    // per-pipe junction ends (-1 = no junction at that end)
    std::vector<int> jIn(np, -1), jOut(np, -1);
    for (int j = 0; j < nj; ++j)
        for (const EndRef& e : junctions_[j].ends)
            (e.atOutlet ? jOut : jIn)[e.pipe] = j;

    auto endBoundary = [&](int i, bool atInlet) -> const TwoFluidBoundary*
    {
        const std::string key = cfg_.pipes[i].name
                              + (atInlet ? ".inlet" : ".outlet");
        const auto it = cfg_.boundaries.find(key);
        return it == cfg_.boundaries.end() ? nullptr : &it->second;
    };

    // pipe classification: free = pressure-driven (unknown velocity) -
    // either a fixed-pressure far end or a span between two junctions;
    // otherwise the velocity is imposed (rate far end, or sealed dead leg)
    std::vector<bool> freePipe(np, false);
    std::vector<double> pFarIn(np, 0.0), pFarOut(np, 0.0), uImposed(np, 0.0);
    std::vector<int> pos(np, -1); // row/column of the free pipes
    int nf = 0;
    for (int i = 0; i < np; ++i)
    {
        const bool span = jIn[i] >= 0 && jOut[i] >= 0;
        const TwoFluidBoundary* far
            = span ? nullptr : endBoundary(i, jIn[i] < 0);
        if (span
            || (far && far->type == TwoFluidBoundary::Type::fixedPressure))
        {
            freePipe[i] = true;
            pos[i] = nf++;
            if (!span)
            {
                if (jIn[i] < 0) pFarIn[i] = far->pressure;
                else pFarOut[i] = far->pressure;
            }
        }
        else if (far && far->type == TwoFluidBoundary::Type::phaseRates)
            uImposed[i] = far->vsl + far->vsg;
        // sealed far end: dead leg, uImposed stays 0
    }

    // unknowns: the free pipes' velocities (pipe order), then the junction
    // pressures; the momentum rows are pipe-diagonal in the velocities but
    // couple to the junction pressures at the pipe's ends (a span couples
    // to both), the balance rows couple the free velocities per junction
    const int pBase = nf;
    std::vector<double> x(nf + nj, 0.0);
    for (int i = 0; i < np; ++i)
        if (freePipe[i]) x[pos[i]] = pipes_[i]->mixtureVelocity();
    for (int j = 0; j < nj; ++j) x[pBase + j] = junctions_[j].p;

    // the pressure at a pipe end: the junction pressure there, or the far
    // fixed pressure
    auto pEnd = [&](int i, bool atInlet, const std::vector<double>& xv)
    {
        const int j = atInlet ? jIn[i] : jOut[i];
        if (j >= 0) return xv[pBase + j];
        return atInlet ? pFarIn[i] : pFarOut[i];
    };
    auto uOf = [&](int i, const std::vector<double>& xv)
    {
        return freePipe[i] ? xv[pos[i]] : uImposed[i];
    };

    // 5.77: latch each pressure-driven pipe ONCE from the step's STARTING
    // junction pressures (the same anti-oscillation rule as 5.49 -- the
    // drive a held pipe sees depends on node pressures this very solve is
    // moving). Rate-driven and sealed legs never hold.
    for (int i = 0; i < np; ++i)
    {
        if (freePipe[i])
            pipes_[i]->latchStaticHold(pEnd(i, true, x) - pEnd(i, false, x));
        else
            pipes_[i]->clearStaticHold();
    }

    std::vector<double> r(nf + nj, 0.0), scale(nf + nj, 0.0);
    // 5.77: a held (or crossing) pipe's row becomes `u = 0`. Its own
    // columns into the junction pressures go away with it -- a held pipe
    // does not answer to the drive, it refuses it -- while the junction
    // balance rows keep their velocity column, so continuity still sees
    // the forced zero and the volume the pipe occupies stays sealed.
    std::vector<bool> heldRow(nf, false);
    auto residuals = [&](const std::vector<double>& xv)
    {
        for (int i = 0; i < np; ++i)
        {
            if (!freePipe[i]) continue;
            if (heldRow[pos[i]])
                r[pos[i]] = xv[pos[i]];
            else
                r[pos[i]] = pipes_[i]->momentumResidual(
                    dt, xv[pos[i]], pEnd(i, true, xv) - pEnd(i, false, xv));
        }
        for (int j = 0; j < nj; ++j)
        {
            double q = 0;
            for (const EndRef& e : junctions_[j].ends)
                q += (e.atOutlet ? 1.0 : -1.0) * pipes_[e.pipe]->area()
                   * uOf(e.pipe, xv);
            r[pBase + j] = q;
        }
    };
    for (int i = 0; i < np; ++i)
    {
        if (!freePipe[i]) continue;
        scale[pos[i]] = pipes_[i]->momentumScale(dt)
                      + std::abs(pEnd(i, true, x))
                      + std::abs(pEnd(i, false, x));
    }
    for (int j = 0; j < nj; ++j)
    {
        double qs = 1e-12;
        for (const EndRef& e : junctions_[j].ends)
            qs += pipes_[e.pipe]->area() * (std::abs(uOf(e.pipe, x)) + 1.0);
        scale[pBase + j] = qs;
    }

    // 5.77: unlike a serial pipe, a branch's far-end pressure is a Newton
    // UNKNOWN (the junction pressure). A latch taken from the step's
    // starting junction state can be stale during start-up: the solve
    // itself may lift that drive past the wall. Releasing is one-way
    // within the step (unlatched[]), which is what keeps the switch from
    // oscillating -- the ratchet 5.49 buys with the latch, preserved on
    // the branch that actually needs it.
    std::vector<bool> unlatched(np, false);
    bool ok = false;
    for (int it = 0; it < 60; ++it)
    {
        // 5.77: re-read the crossing witness on the current iterate. A
        // moving column that has gone through zero under a sub-threshold
        // drive parks at u = 0 for this step, exactly like 5.49's
        // `pd.held || crossing` row; the latched held pipes are parked
        // with it. Once parked, u * vMOld_ == 0 keeps it parked.
        for (int i = 0; i < np; ++i)
        {
            if (!freePipe[i]) continue;
            TwoFluidPipe& p = *pipes_[i];
            const double dp = pEnd(i, true, x) - pEnd(i, false, x);
            if (p.staticHold() && !unlatched[i]
                && std::abs(dp - p.gravityHeadPa())
                       > p.yieldThresholdPa())
            {
                p.clearStaticHold();
                unlatched[i] = true;
            }
            const bool crossing = p.yieldCrossing(dt, x[pos[i]], dp);
            if (crossing || p.staticHold())
            {
                x[pos[i]] = 0.0;
                heldRow[pos[i]] = true;
            }
            else heldRow[pos[i]] = false;
        }
        residuals(x);
        ok = true;
        for (int q = 0; q < nf + nj; ++q)
            if (!(std::abs(r[q]) <= 1e-12 * scale[q])) { ok = false; break; }
        if (ok) break;

        std::vector<std::vector<double>> J(
            nf + nj, std::vector<double>(nf + nj, 0.0));
        for (int i = 0; i < np; ++i)
        {
            if (!freePipe[i]) continue;
            TwoFluidPipe& p = *pipes_[i];
            if (heldRow[pos[i]])
            {
                // the forced row has slope 1 in its own velocity and no
                // sensitivity to the junction pressures it refuses
                J[pos[i]][pos[i]] = 1.0;
                continue;
            }
            const double u = x[pos[i]];
            const double dp = pEnd(i, true, x) - pEnd(i, false, x);
            const double h
                = std::max(1e-6 * std::max(std::abs(u), 1e-3), 1e-12);
            double slope;
            if (std::abs(u) <= h && p.yieldThresholdPa() > 0)
            {
                // 5.77: one-sided slope at rest (see momentumSolve). The
                // network Newton's centred stencil would otherwise feed
                // the threshold jump into the global Jacobian and make
                // the junction pressure column singular.
                const double u2 = u + (u >= 0 ? h : -h);
                slope = (p.momentumResidual(dt, u2, dp) - r[pos[i]])
                      / (u2 - u);
            }
            else
            {
                slope = (p.momentumResidual(dt, u + h, dp)
                       - p.momentumResidual(dt, u - h, dp)) / (2.0 * h);
            }
            if (!(slope > 0))
                cfFatal("twoFluid: non-monotone pipe momentum at pipe '"
                        + cfg_.pipes[i].name + "'");
            J[pos[i]][pos[i]] = slope;
            // dR/dpEnd = -1 each way, so dR/dpJ = -d(dpDrive)/dpJ:
            // -1 when the junction sits at the pipe's inlet, +1 at its
            // outlet (a span pipe carries both)
            if (jIn[i] >= 0) J[pos[i]][pBase + jIn[i]] = -1.0;
            if (jOut[i] >= 0) J[pos[i]][pBase + jOut[i]] = 1.0;
        }
        for (int j = 0; j < nj; ++j)
            for (const EndRef& e : junctions_[j].ends)
                if (freePipe[e.pipe])
                    J[pBase + j][pos[e.pipe]]
                        = (e.atOutlet ? 1.0 : -1.0)
                        * pipes_[e.pipe]->area();
        std::vector<double> dxv = r;
        for (double& v : dxv) v = -v;
        detail::solveDenseLinear(J, dxv);
        for (int q = 0; q < nf + nj; ++q) x[q] += dxv[q];
    }
    if (!ok)
        cfFatal("twoFluid: network Newton did not converge");
    for (int j = 0; j < nj; ++j) junctions_[j].p = x[pBase + j];

    // every pipe takes its velocity (free pipes the converged one), then
    // per junction the inflow streams mix: holdup/water/hydrate volume-flow
    // weighted, temperature heat-capacity-flow weighted, gas density
    // gas-mass-flux weighted (with the EOS this closes the junction gas
    // mass balance exactly)
    for (int i = 0; i < np; ++i)
        pipes_[i]->setExternalMomentum(uOf(i, x));

    const TwoFluidFluid& f = cfg_.perPipe.fluid;
    for (int j = 0; j < nj; ++j)
    {
        Junction& jn = junctions_[j];
        double qL = 0, qG = 0, wNum = 0, hNum = 0, waxNum = 0, enNum = 0,
               iNum = 0, sNum = 0, eNum = 0, eDen = 0, gNum = 0;
        // G1: per-component liquid-flow-weighted mixing (empty when no
        // composition{} feed was declared)
        const int nz = static_cast<int>(cfg_.perPipe.compositionZ0.size());
        std::vector<double> zNum(nz, 0.0);
        for (const EndRef& e : jn.ends)
        {
            const double s = e.atOutlet ? 1.0 : -1.0;
            if (s * uOf(e.pipe, x) <= 0) continue; // fed or static: no inflow
            const TwoFluidPipe& p = *pipes_[e.pipe];
            const bool atInlet = !e.atOutlet;
            const double fl = p.endLiquidFlux(atInlet);
            const double fg = p.endGasVolumeFlux(atInlet);
            if (fl * s < 0 || fg * s < 0)
                cfFatal("twoFluid: countercurrent phase flux at junction '"
                        + jn.name + "' not supported");
            const double aL = p.area() * std::abs(fl);
            const double aG = p.area() * std::abs(fg);
            const double rhoG = p.endGasDensity(atInlet);
            qL += aL;
            qG += aG;
            wNum += aL * p.endWater(atInlet);
            hNum += aL * p.endHydrate(atInlet);
            waxNum += aL * p.endWax(atInlet);
            enNum += aL * p.endEntrained(atInlet);
            iNum += aL * p.endInhibitor(atInlet);
            sNum += aL * p.endSand(atInlet);
            for (int k = 0; k < nz; ++k)
                zNum[k] += aL * p.endComposition(k, atInlet);
            eNum += (aL * f.rhoL * f.cL + aG * rhoG * f.cG)
                  * p.endTemperature(atInlet);
            eDen += aL * f.rhoL * f.cL + aG * rhoG * f.cG;
            gNum += aG * rhoG;
        }
        const double qM = qL + qG;
        if (qM > 0)
        {
            jn.alphaJ = qL / qM;
            if (qL > 0)
            {
                jn.waterJ = wNum / qL;
                jn.hydrateJ = hNum / qL;
                jn.waxJ = waxNum / qL;
                jn.entrJ = enNum / qL;
                jn.inhibJ = iNum / qL;
                jn.sandJ = sNum / qL;
                if (nz > 0)
                    for (int k = 0; k < nz; ++k) jn.zJ[k] = zNum[k] / qL;
            }
            if (eDen > 0) jn.tJ = eNum / eDen;
            if (qG > 0) jn.rhoJ = gNum / qG;
        }

        // junction faces: feeding branches get a fixed-pressure ghost at
        // the mixed state (their face flux upwinds from the boundary cell);
        // supplied branches get the junction split imposed as phaseRates at
        // the pipe's own flux functions at the mixed holdup (F + G = vM
        // keeps the volume balance exact for every closure; the phase split
        // is exact for the homogeneous closure, closure-approximate with
        // slip). A span pipe can feed one junction and be supplied by the
        // other - each junction mutates only its own end.
        for (const EndRef& e : jn.ends)
        {
            TwoFluidPipe& p = *pipes_[e.pipe];
            const double s = e.atOutlet ? 1.0 : -1.0;
            TwoFluidBoundary b;
            b.pressure = jn.p;
            b.alpha = jn.alphaJ;
            b.waterFraction = jn.waterJ;
            b.hydrateFraction = jn.hydrateJ;
            b.waxDepositFraction = jn.waxJ;
            b.entrainedFraction = jn.entrJ;
            b.inhibitorFraction = jn.inhibJ;
            b.sandFraction = jn.sandJ;
            if (nz > 0) b.zInlet = jn.zJ; // G1
            b.temperature = jn.tJ;
            b.rhoGInjected = jn.rhoJ;
            if (s * uOf(e.pipe, x) > 0)
            {
                b.type = TwoFluidBoundary::Type::fixedPressure;
            }
            else
            {
                // Conservative split: each supplied branch takes the
                // junction's mixed phase fractions applied to its own
                // mixture rate. Summing the branch rates then reproduces
                // the junction's phase rates exactly, so a slip closure can
                // no longer create or destroy phase volume at the node.
                const double vMag = std::abs(uOf(e.pipe, x));
                b.type = TwoFluidBoundary::Type::phaseRates;
                b.vsl = vMag * jn.alphaJ;
                b.vsg = vMag * (1.0 - jn.alphaJ);
            }
            if (e.atOutlet) p.setOutletBoundary(b);
            else p.setInletBoundary(b);
        }
    }
}

std::vector<std::string> TwoFluidNetwork::pipeNames() const
{
    std::vector<std::string> names;
    for (const TwoFluidGeo& g : cfg_.pipes) names.push_back(g.name);
    return names;
}

const TwoFluidPipe& TwoFluidNetwork::pipe(const std::string& name) const
{
    return *pipes_.at(pipeIndex_.at(name));
}

std::vector<std::string> TwoFluidNetwork::junctionNames() const
{
    std::vector<std::string> names;
    for (const Junction& jn : junctions_) names.push_back(jn.name);
    return names;
}

double TwoFluidNetwork::junctionPressure(const std::string& name) const
{
    for (const Junction& jn : junctions_)
        if (jn.name == name) return jn.p;
    cfFatal("twoFluid: unknown junction '" + name + "'");
    return 0;
}

double TwoFluidNetwork::junctionNetRate(const std::string& name) const
{
    double net = 0;
    for (double q : junctionPipeRates(name)) net += q;
    return net;
}

std::vector<double> TwoFluidNetwork::junctionPipeRates(
    const std::string& name) const
{
    for (const Junction& jn : junctions_)
        if (jn.name == name)
        {
            std::vector<double> rates;
            for (const EndRef& e : jn.ends)
                rates.push_back((e.atOutlet ? 1.0 : -1.0)
                                * pipes_[e.pipe]->area()
                                * pipes_[e.pipe]->mixtureVelocity());
            return rates;
        }
    cfFatal("twoFluid: unknown junction '" + name + "'");
    return {};
}

} // namespace cf
