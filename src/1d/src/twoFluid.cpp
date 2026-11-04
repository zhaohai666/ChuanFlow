#include "cf/1d/TwoFluid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "cf/common/Error.h"

#include "cf/1d/Network.h" // churchillFriction
#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double mixtureRho(const TwoFluidFluid& f, double a)
{
    return a * f.rhoL + (1.0 - a) * f.rhoG;
}

double mixtureMu(const TwoFluidFluid& f, double a)
{
    return a * f.muL + (1.0 - a) * f.muG;
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

TwoFluidPipe::TwoFluidPipe(const TwoFluidCase& c)
    : cfg_(c)
{
    const TwoFluidFluid& f = cfg_.fluid;
    if (!(f.rhoL > 0) || !(f.rhoG > 0) || !(f.muL > 0) || !(f.muG > 0))
        cfFatal("twoFluid: fluid needs positive densities and viscosities");
    const TwoFluidGeo& g = cfg_.geo;
    if (!(g.length > 0) || !(g.diameter > 0) || g.sections < 1)
        cfFatal("twoFluid: pipe needs positive length/diameter and >= 1"
                " sections");
    if (g.roughness < 0)
        cfFatal("twoFluid: negative roughness");

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
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
        cfFatal("twoFluid: phaseRates injection with negative total rate is"
                " not supported");

    n_ = g.sections;
    dx_ = g.length / n_;
    sinTheta_ = std::sin(g.angleDeg * M_PI / 180.0);
    alpha_.assign(n_, cfg_.alpha0);
    water_.assign(n_, cfg_.water0);
    for (const TwoFluidFill& fill : cfg_.fills)
    {
        if (fill.from > fill.to || fill.alpha < 0 || fill.alpha > 1)
            cfFatal("twoFluid: fill block needs from <= to and alpha in"
                    " [0,1]");
        for (int i = 0; i < n_; ++i)
        {
            const double x = (i + 0.5) * dx_;
            if (x >= fill.from && x <= fill.to) alpha_[i] = fill.alpha;
        }
    }
    // incompressible phases: an imposed rate is felt along the whole pipe
    // from t = 0; pressure-driven and sealed pipes start from rest
    if (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates)
        vM_ = vMOld_ = cfg_.inlet.vsl + cfg_.inlet.vsg;

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
            const double dExit = th.inletDepth
                               - cfg_.geo.length * sinTheta_;
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
}

void TwoFluidPipe::fluxCoeffs(double& qa, double& qb, double& qc) const
{
    // F(a) = a_l * vL(a): liquid volume flux per unit area
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

double TwoFluidPipe::liquidFlux(double alpha) const
{
    double qa, qb, qc;
    fluxCoeffs(qa, qb, qc);
    return qc + alpha * (qb + qa * alpha);
}

double TwoFluidPipe::gasVelocity(double alpha) const
{
    if (cfg_.closure.model == 0)
        return vM_ + alpha * cfg_.closure.slipVelocity;
    return cfg_.closure.c0 * vM_ + cfg_.closure.driftVelocity;
}

double TwoFluidPipe::cellFriction(int i, double vM) const
{
    // homogeneous Churchill law on mixture properties; |dp| per cell
    const TwoFluidFluid& f = cfg_.fluid;
    const double a = alpha_[i];
    const double rho = mixtureRho(f, a);
    const double mu = mixtureMu(f, a);
    const double v = std::abs(vM);
    if (v == 0) return 0.0;
    const double re = rho * v * cfg_.geo.diameter / mu;
    const double fr = churchillFriction(re, cfg_.geo.roughness
                                               / cfg_.geo.diameter)
                    * (dx_ / cfg_.geo.diameter) * 0.5 * rho * v * v;
    return fr;
}

std::vector<double> TwoFluidPipe::cellPressureLosses(double vM) const
{
    // signed pressure loss per cell: friction (odd in vM) + gravity head
    std::vector<double> loss(n_);
    for (int i = 0; i < n_; ++i)
    {
        const double grav = mixtureRho(cfg_.fluid, alpha_[i]) * kGravity
                          * sinTheta_ * dx_;
        loss[i] = (vM >= 0 ? cellFriction(i, vM) : -cellFriction(i, vM)) + grav;
    }
    return loss;
}

void TwoFluidPipe::momentumSolve(double dt)
{
    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;

    if (in.type == TwoFluidBoundary::Type::phaseRates)
    {
        vM_ = in.vsl + in.vsg;
        return;
    }
    if (in.type == TwoFluidBoundary::Type::sealed)
    {
        // sealed both ends: the mixture volume flux vanishes identically
        vM_ = 0.0;
        return;
    }

    // fixed-pressure inlet: implicit Euler on the mixture momentum with the
    // holdup frozen at the old time level
    double rhoBar = 0;
    for (int i = 0; i < n_; ++i)
        rhoBar += mixtureRho(cfg_.fluid, alpha_[i]);
    rhoBar /= n_;
    const double area = M_PI * cfg_.geo.diameter * cfg_.geo.diameter / 4.0;
    const double inertia = rhoBar * cfg_.geo.length / area;
    double gravSum = 0;
    for (int i = 0; i < n_; ++i)
        gravSum += mixtureRho(cfg_.fluid, alpha_[i]) * kGravity * sinTheta_
                 * dx_;
    const double dpDrive = in.pressure - out.pressure;

    auto residual = [&](double u)
    {
        double fr = 0;
        for (int i = 0; i < n_; ++i)
            fr += (u >= 0 ? cellFriction(i, u) : -cellFriction(i, u));
        return inertia * (u - vMOld_) / dt + fr + gravSum - dpDrive;
    };

    double u = vMOld_;
    const double scale = inertia * (std::abs(vMOld_) + 1.0) / dt
                       + std::abs(dpDrive) + std::abs(gravSum) + 1.0;
    for (int it = 0; it < 60; ++it)
    {
        const double r = residual(u);
        if (std::abs(r) <= 1e-12 * scale) break;
        const double h = std::max(1e-6 * std::max(std::abs(u), 1e-3), 1e-12);
        const double slope = (residual(u + h) - residual(u - h)) / (2.0 * h);
        if (!(slope > 0))
            cfFatal("twoFluid: non-monotone momentum residual");
        u -= r / slope;
    }
    if (!std::isfinite(u)) cfFatal("twoFluid: momentum solve diverged");
    vM_ = u;
}

void TwoFluidPipe::transportStep(double dt)
{
    double qa, qb, qc;
    fluxCoeffs(qa, qb, qc);

    // CFL guard for the explicit Godunov transport (monotone under
    // dt * max|F'| / dx <= 1; F' is monotone so the extremes sit at 0 and 1)
    const double cmax = std::max(std::abs(qb), std::abs(2.0 * qa + qb));
    if (dt * cmax / dx_ > 1.0 + 1e-12)
        cfFatal("twoFluid: CFL violated (dt * |dF/dalpha| / dx = "
                + std::to_string(dt * cmax / dx_)
                + ") - reduce dt");

    const TwoFluidBoundary& in = cfg_.inlet;
    const TwoFluidBoundary& out = cfg_.outlet;
    std::vector<double> flux(n_ + 1), wflux(n_ + 1);

    for (int j = 0; j <= n_; ++j)
    {
        if (j == 0)
        {
            if (in.type == TwoFluidBoundary::Type::sealed) flux[0] = 0;
            else if (in.type == TwoFluidBoundary::Type::phaseRates)
                flux[0] = in.vsl; // injected liquid superficial velocity
            else
                flux[0] = godunovQuadratic(qa, qb, qc, in.alpha, alpha_[0]);
        }
        else if (j == n_)
        {
            if (out.type == TwoFluidBoundary::Type::sealed) flux[n_] = 0;
            else flux[n_] = liquidFlux(alpha_[n_ - 1]); // zero gradient
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
            wflux[n_] = water_[n_ - 1] * flux[n_];
        }
        else
        {
            wflux[j] = (flux[j] >= 0 ? water_[j - 1] : water_[j]) * flux[j];
        }
    }

    static constexpr double kTol = 1e-9;
    for (int i = 0; i < n_; ++i)
    {
        const double aNew = alpha_[i] + dt / dx_ * (flux[i] - flux[i + 1]);
        if (aNew < -kTol || aNew > 1.0 + kTol)
            cfFatal("twoFluid: holdup left [0,1] - reduce dt");
        const double afOld = alpha_[i] * water_[i];
        double afNew = afOld + dt / dx_ * (wflux[i] - wflux[i + 1]);
        const double aCap = std::max(aNew, 0.0);
        if (afNew < -kTol || afNew > aCap + kTol)
            cfFatal("twoFluid: water cut left [0, holdup] - reduce dt");
        alpha_[i] = std::min(std::max(aNew, 0.0), 1.0);
        afNew = std::min(std::max(afNew, 0.0), aCap);
        water_[i] = aCap > kTol ? afNew / aCap : 0.0;
    }
}

void TwoFluidPipe::step(double dt)
{
    if (!(dt > 0)) cfFatal("twoFluid: dt must be positive");
    momentumSolve(dt);
    transportStep(dt);
    thermalStep(dt);
    vMOld_ = vM_;
    t_ += dt;
}

double TwoFluidPipe::wallTemperature(int i) const
{
    const TwoFluidThermal& th = cfg_.thermal;
    if (th.wallModel == 0) return th.wallTemperature;
    const double depth = th.inletDepth - (i + 0.5) * dx_ * sinTheta_;
    return th.surfaceTemperature + th.geothermalGradient * depth;
}

void TwoFluidPipe::thermalStep(double dt)
{
    if (!cfg_.thermalOn) return;
    const TwoFluidFluid& f = cfg_.fluid;
    const TwoFluidThermal& th = cfg_.thermal;
    const double area = M_PI * cfg_.geo.diameter * cfg_.geo.diameter / 4.0;
    auto rhoc = [&](double a)
    {
        return a * f.rhoL * f.cL + (1.0 - a) * f.rhoG * f.cG;
    };

    // |K| per face: enthalpy conductance W/K of the mixture stream, with the
    // upwind composition of each face (ghost cells at the open boundaries)
    std::vector<double> kFace(n_ + 1, 0.0);
    if (vM_ > 0)
    {
        double rhocUp;
        if (cfg_.inlet.type == TwoFluidBoundary::Type::phaseRates)
            rhocUp = (cfg_.inlet.vsl * f.rhoL * f.cL
                      + cfg_.inlet.vsg * f.rhoG * f.cG) / vM_;
        else
            rhocUp = rhoc(cfg_.inlet.alpha);
        kFace[0] = area * vM_ * rhocUp;
        for (int j = 1; j <= n_; ++j)
            kFace[j] = area * vM_ * rhoc(alpha_[j - 1]);
    }
    else if (vM_ < 0)
    {
        for (int j = 0; j < n_; ++j)
            kFace[j] = area * (-vM_) * rhoc(alpha_[j]);
        kFace[n_] = area * (-vM_) * rhoc(cfg_.outlet.alpha);
    }

    // implicit Euler: diagonal = capacity + outflow + wall; the inflow face
    // couples the upwind neighbour (or the ghost temperature at the boundary)
    std::vector<std::vector<double>> mat(n_, std::vector<double>(n_, 0.0));
    std::vector<double> rhs(n_, 0.0);
    for (int i = 0; i < n_; ++i)
    {
        const double cap = area * dx_ * rhoc(alpha_[i]);
        const double kw = th.wallHeatTransfer * dx_;
        mat[i][i] = cap / dt + kFace[vM_ >= 0 ? i + 1 : i] + kw;
        rhs[i] = cap / dt * T_[i] + kw * wallTemperature(i);
        if (vM_ > 0)
        {
            if (i > 0) mat[i][i - 1] = -kFace[i];
            else rhs[i] += kFace[0] * cfg_.inlet.temperature;
        }
        else if (vM_ < 0)
        {
            if (i < n_ - 1) mat[i][i + 1] = -kFace[i + 1];
            else rhs[i] += kFace[n_] * cfg_.outlet.temperature;
        }
    }
    detail::solveDenseLinear(mat, rhs);

    inletHeatRate_ = 0;
    if (vM_ > 0) inletHeatRate_ = kFace[0] * cfg_.inlet.temperature;
    else if (vM_ < 0) inletHeatRate_ = kFace[n_] * cfg_.outlet.temperature;
    wallHeatRate_ = 0;
    for (int i = 0; i < n_; ++i)
    {
        wallHeatRate_ += th.wallHeatTransfer * dx_
                       * (wallTemperature(i) - rhs[i]);
        T_[i] = rhs[i];
    }
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
    for (int i = 0; i < n_; ++i) v[i] = liquidFlux(alpha_[i]);
    return v;
}

std::vector<double> TwoFluidPipe::cellVsg() const
{
    std::vector<double> v(n_);
    for (int i = 0; i < n_; ++i)
        v[i] = (1.0 - alpha_[i]) * gasVelocity(alpha_[i]);
    return v;
}

TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict)
{
    TwoFluidCase c;
    if (!caseDict.has("fluid")) cfFatal("twoFluid: case has no fluid{} block");
    const Dictionary& fd = caseDict.getDict("fluid");
    c.fluid.rhoL = fd.getScalar("rhoL");
    c.fluid.rhoG = fd.getScalar("rhoG");
    c.fluid.muL = fd.getScalar("muL");
    c.fluid.muG = fd.getScalar("muG");
    c.fluid.cL = fd.getScalarOr("cL", 0.0);
    c.fluid.cG = fd.getScalarOr("cG", 0.0);

    if (!caseDict.has("pipe")) cfFatal("twoFluid: case has no pipe{} block");
    const Dictionary& pd = caseDict.getDict("pipe");
    c.geo.name = pd.getWordOr("name", "pipe");
    c.geo.length = pd.getScalar("length");
    c.geo.diameter = pd.getScalar("diameter");
    c.geo.roughness = pd.getScalarOr("roughness", 0.0);
    c.geo.angleDeg = pd.getScalarOr("angleDeg", 0.0);
    c.geo.sections = static_cast<int>(pd.getScalarOr("sections", 50.0));

    if (!caseDict.has("closure"))
        cfFatal("twoFluid: case has no closure{} block");
    const Dictionary& cd = caseDict.getDict("closure");
    const std::string model = cd.getWordOr("model", "wallis");
    if (model == "wallis") c.closure.model = 0;
    else if (model == "zuber") c.closure.model = 1;
    else cfFatal("twoFluid: unknown closure model '" + model + "'");
    c.closure.slipVelocity = cd.getScalarOr("slipVelocity", 0.0);
    c.closure.c0 = cd.getScalarOr("C0", 1.0);
    c.closure.driftVelocity = cd.getScalarOr("driftVelocity", 0.0);

    if (!caseDict.has("initial"))
        cfFatal("twoFluid: case has no initial{} block");
    const Dictionary& id = caseDict.getDict("initial");
    c.alpha0 = id.getScalarOr("alpha", 0.0);
    c.water0 = id.getScalarOr("water", 0.0);
    c.temperature0 = id.getScalarOr("temperature", 300.0);
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

    if (!caseDict.has("boundaries"))
        cfFatal("twoFluid: case has no boundaries{} block");
    const Dictionary& bd = caseDict.getDict("boundaries");
    auto loadBoundary = [](const Dictionary& d) {
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
        b.pressure = d.getScalarOr("p", 0.0);
        b.alpha = d.getScalarOr("alpha", 0.0);
        b.waterFraction = d.getScalarOr("waterFraction", 0.0);
        b.temperature = d.getScalarOr("temperature", 0.0);
        return b;
    };
    if (!bd.has("inlet") || !bd.has("outlet"))
        cfFatal("twoFluid: boundaries{} needs inlet{} and outlet{}");
    c.inlet = loadBoundary(bd.getDict("inlet"));
    c.outlet = loadBoundary(bd.getDict("outlet"));

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

    return c;
}

} // namespace cf
