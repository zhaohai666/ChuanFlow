#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"

namespace cf
{

/// M2 second increment: distributed 1D drift-flux two-fluid pipe model (the
/// OLGA-style three-equation core in its incompressible-phase limit). The
/// state along the pipe is the per-section liquid volume fraction alpha
/// (gas holds 1 - alpha) plus a passive water cut transported with the
/// liquid:
///   d(alpha)/dt + d(alpha vL)/dx = 0                 [liquid mass]
///   alpha_g = 1 - alpha (volume conservation)        [gas mass]
///   rho_m L du/dt = dp - sum(friction) - sum(rho_m g sin(theta) dx)
///                                                    [mixture momentum]
/// with a slip closure turning alpha into phase velocities. The algebraic
/// closures are the drift-flux form vG = C0 vM + vD (Zuber-Findlay) and the
/// constant-slip form vG = vM + alpha us (Wallis). The dynamic closure
/// carries a per-section slip state s = vG - vL derived from the two phase
/// momenta with the shared dp/dx eliminated (the reduced single-pressure
/// form of the two-momentum model):
///   K_s ds/dt = (1/rhoG - 1/rhoL)(rho_m du/dt + loss/dx) - D sigma,
/// K_s = 1 + alpha alpha_g (rhoG - rhoL)^2/(rhoL rhoG) the slip-inertia
/// factor, D = 1/(rhoL alpha) + 1/(rhoG alpha_g), sigma = C s |s| the
/// interfacial drag per unit volume (C = interfacialDrag, kg/m4) and loss
/// the per-cell friction + gravity pressure loss. The drag is treated
/// semi-implicitly (secant at the old slip), which is unconditionally
/// stable and relaxes s towards the terminal slip s* =
/// sign(drive) sqrt(|drive|/(C D)) with the geometric ratio
/// 1/(1 + dt C |s*| D / K_s). The mixture momentum feels the slip inertia
/// alpha alpha_g (rhoG - rhoL) ds/dt lagged one step; it vanishes at steady
/// state, so the steady holdups of the algebraic closures are reproduced.
/// Phases are incompressible and the area is constant, so the mixture
/// velocity is uniform along the pipe. With the gas EOS on (fluid gasR > 0)
/// the gas density becomes a per-section state evolved by gas mass
/// continuity on the fixed gas volume (1 - alpha, the sections stay full):
///   d(rhoG (1 - alpha))/dt + d(rhoG alpha_g vG)/dx = 0   [gas mass]
///   p_i = rhoG_i gasR T_i                                [ideal-gas EOS]
/// The gas volume flux alpha_g vG is a quadratic in alpha for all three
/// closures, so the same Godunov machinery transports gas mass; the
/// thermodynamic pressure p_i follows from the transported density and
/// feeds back (lagged one step) into every mixture density: the momentum
/// inertia/gravity/friction weights, the slip drive and the heat capacity.
/// The mixture momentum itself stays lumped (one velocity per pipe) - the
/// thermodynamic p_i is not a local momentum driver, so acoustic dynamics
/// and the friction pressure profile inside the pipe are beyond this
/// reduced form (distributed pressure momentum is the follow-up).
/// Holdup transport advances with the exact Godunov flux of the
/// quadratic liquid-flux function F(alpha) = alpha vL(alpha), so shocks
/// propagate at the Rankine-Hugoniot speed and the CFL-limited scheme is
/// monotone (alpha stays in [0, 1]).
/// Wall friction is the homogeneous Churchill law on mixture properties;
/// gravity uses the local mixture density.
/// With a thermal{} block the sections also carry a mixture temperature (the
/// phases stay in thermal equilibrium): implicit Euler on
///   C_i dT_i/dt = H_in,i - H_out,i + k' dx (T_wall,i - T_i),
/// C_i the section's phase-weighted rho*c capacity and the enthalpy flux
/// riding the phase volume fluxes of the step's transport, each phase
/// upwinded on its own flux direction (with slip the liquid no longer moves
/// at the mixture velocity; the slip-0 limit reproduces the homogeneous
/// mixture-weighted form). Implicit in T, so no CFL on top
/// of the transport one.
/// A case may chain several pipe segments (pipes{} block of named
/// sub-blocks in flow order - e.g. a vertical riser followed by a horizontal
/// flowline, a well trajectory). Each segment carries its own
/// length/diameter/roughness/angle/sections; the cells are concatenated and
/// every junction face is an ordinary interior face of the same Godunov
/// transport (the junction slip averages the two neighbouring states). The
/// chain needs a uniform bore: the mixture momentum is lumped on ONE
/// velocity, which requires a constant area. Per-section dx/slope/relative
/// roughness drive the friction, gravity head, thermal capacity and CFL per
/// cell; the fills' from/to coordinates measure the arc position from the
/// chain inlet.
/// With a hydrate{} block the sections also carry a hydrate volume fraction
/// h of the liquid+hydrate slurry (riding the liquid like the water cut).
/// Growth is driven by the subcooling dT = T_eq(p) - T with the equilibrium
/// curve ln(p[Pa]) = curveA + curveB / T_eq (curveB < 0), first order in
/// the free-water fraction and implicit in h:
///   dh/dt =  growthRate * dT * (water + h)   (dT > 0, growth)
///   dh/dt = -dissociationRate * |dT| * h     (dT < 0, dissociation)
/// (the implicit growth law converges to the full conversion h = water + h,
/// no separate nucleation threshold in this form). Formation consumes
/// gasContent kg of gas per m3 of hydrate - a real gas mass sink on the EOS
/// inventory, capped by the available gas - and moves the free-water volume
/// into hydrate 1:1 (the volume-preserving convention of the transient
/// network kernel). Formation releases latentHeat J/m3 (dissociation
/// absorbs it) as a source in the temperature equation. Requires the gas
/// EOS (fluid gasR > 0) and a thermal{} block: both p_i and T_i drive the
/// kinetics.

struct TwoFluidFluid
{
    double rhoL = 1000; // kg/m3
    double rhoG = 20;   // kg/m3 (with the EOS on: the gas density the
                        // injected rates and the initial state refer to)
    double muL = 0.01;  // Pa s
    double muG = 1e-5;  // Pa s
    double cL = 0;      // J/(kg K), needed when the thermal block is present
    double cG = 0;      // J/(kg K)
    double gasR = 0;    // J/(kg K) specific gas constant; > 0 switches the
                        // ideal-gas EOS on and rhoG becomes the per-section
                        // state evolved by gas mass continuity
};

struct TwoFluidGeo
{
    std::string name;
    double length = 0;   // m
    double diameter = 0; // m
    double roughness = 0;
    double angleDeg = 0; // from horizontal, positive = uphill inlet -> outlet
    int sections = 50;   // number of holdup cells along the pipe
};

struct TwoFluidClosure
{
    int model = 0;             // 0 = wallis (constant slip velocity),
                               // 1 = zuber (vG = C0 vM + vD),
                               // 2 = dynamic (transported slip state)
    double slipVelocity = 0;   // wallis us, m/s (positive = gas runs ahead)
    double c0 = 1.0;           // zuber profile parameter
    double driftVelocity = 0;  // zuber vD, m/s
    double interfacialDrag = 0;// dynamic: kg/m4, drag law sigma = C s |s|
                               // per unit volume; required > 0 for model 2
};

struct TwoFluidBoundary
{
    enum class Type { sealed, phaseRates, fixedPressure };
    Type type = Type::sealed;
    double vsl = 0;          // phaseRates: imposed liquid superficial velocity
    double vsg = 0;          // phaseRates: imposed gas superficial velocity
    double pressure = 0;     // fixedPressure: Pa
    double alpha = 0;        // fixedPressure inlet: composition of the
                             // inflow ghost cell
    double waterFraction = 0;// water cut of the inflow liquid
    double hydrateFraction = 0; // hydrate fraction of the inflow liquid
    double temperature = 0;  // K, temperature of the inflow ghost cell
                             // (required on both ends when thermal is on,
                             // the unused one for potential backflow)
};

/// Wall/formation heat exchange for the mixture temperature equation.
/// wallHeatTransfer is an overall coefficient per unit pipe length (W/(m K)):
/// the film, wall and cement resistances folded into one number, so the
/// exchange per section is k' dx (T_wall - T). The wall temperature is either
/// constant or a geothermal profile T = surfaceTemperature + gradient * d,
/// with d the true vertical depth below surface (inletDepth at s = 0,
/// decreasing uphill along the pipe).
struct TwoFluidThermal
{
    double wallHeatTransfer = 0;   // W/(m K), 0 = insulated wall
    int wallModel = 0;             // 0 = constant, 1 = geothermal
    double wallTemperature = 0;    // K (constant model)
    double surfaceTemperature = 0; // K (geothermal at depth 0)
    double geothermalGradient = 0; // K/m
    double inletDepth = 0;         // m below surface at s = 0
};

struct TwoFluidFill
{
    double from = 0;  // m along the pipe (cells whose centre is inside)
    double to = 0;
    double alpha = 0;
};

/// Gas-hydrate kinetics (see the class documentation above). The equilibrium
/// curve ln(p[Pa]) = curveA + curveB / T_eq is the first-order Clausius
/// form of the dissociation line; the defaults anchor roughly on methane
/// hydrate (T_eq ~ 295 K at 3 MPa). gasContent is the stoichiometric gas
/// mass per m3 of hydrate (CH4 x 5.75 H2O ~ 122 kg/m3) and latentHeat the
/// formation heat release per m3 (~ 4.2e8 J/m3).
struct TwoFluidHydrate
{
    double curveA = 35.2;            // ln(p[Pa]) = curveA + curveB / T_eq
    double curveB = -5684.0;         // K (negative: p up -> T_eq up)
    double growthRate = 0.01;        // 1/(s K), first order in free water
    double dissociationRate = 0.01;  // 1/(s K), first order in h
    double gasContent = 122.0;       // kg gas consumed per m3 hydrate
    double latentHeat = 4.2e8;       // J released per m3 hydrate formed
    double initial = 0;              // initial hydrate fraction of the slurry
};

struct TwoFluidCase
{
    TwoFluidFluid fluid;
    std::vector<TwoFluidGeo> pipes; // serial chain in flow order
    TwoFluidClosure closure;
    TwoFluidBoundary inlet;
    TwoFluidBoundary outlet;
    double alpha0 = 0;    // uniform initial liquid fraction
    double water0 = 0;    // uniform initial water cut of the liquid
    std::vector<TwoFluidFill> fills; // sub-range overrides of alpha0
    bool thermalOn = false;      // energy equation solved (thermal{} present)
    double temperature0 = 300;   // K, uniform initial temperature
    TwoFluidThermal thermal;
    bool hydrateOn = false;      // hydrate kinetics solved (hydrate{} present)
    TwoFluidHydrate hydrate;
};

class TwoFluidPipe
{
public:
    explicit TwoFluidPipe(const TwoFluidCase& c);

    /// One time step: implicit-Euler momentum for the mixture velocity
    /// (frozen holdup), the slip update, then explicit Godunov transport of
    /// alpha and the water cut. Fatals when dt violates the CFL condition.
    void step(double dt);

    double time() const { return t_; }
    int sections() const { return n_; }
    /// Per-cell arc step, m (one entry per section along the chain).
    std::vector<double> cellSizes() const { return dxs_; }

    /// Signed mixture (total volume-flux) velocity, m/s.
    double mixtureVelocity() const { return vM_; }
    /// p(inlet) - p(outlet) from the momentum balance, Pa.
    double pressureDrop() const;
    /// Node pressures at the cell faces, Pa (size sections + 1). Gauge
    /// relative to the outlet for rate-driven and sealed cases.
    std::vector<double> nodePressure() const;

    std::vector<double> cellAlpha() const { return alpha_; }
    std::vector<double> cellWater() const;  // water cut of the liquid
    std::vector<double> cellVsl() const;    // superficial liquid velocity
    std::vector<double> cellVsg() const;    // superficial gas velocity
    /// Slip velocity vG - vL per section, m/s; zeros unless the dynamic
    /// closure is on (the algebraic closures carry no slip state).
    std::vector<double> cellSlip() const { return slip_; }
    /// Thermodynamic section pressure from the gas EOS (p = rhoG gasR T),
    /// Pa; zeros unless the EOS is on. This is the transported density's
    /// pressure, not the hydraulic loss profile (see nodePressure()).
    std::vector<double> cellPressure() const;
    std::vector<double> cellTemperature() const { return T_; }
    /// Hydrate volume fraction of the liquid+hydrate slurry per section;
    /// zeros unless the hydrate{} block is present.
    std::vector<double> cellHydrate() const { return hydrate_; }
    /// Wall heat the sections absorb this step, W (positive = into the fluid);
    /// zero unless thermal is on.
    double wallHeatRate() const { return wallHeatRate_; }
    /// Enthalpy the inflow boundaries inject, W (transport + sources only).
    double inletHeatRate() const { return inletHeatRate_; }

private:
    void momentumSolve(double dt);
    void slipStep(double dt);
    void transportStep(double dt);
    void hydrateStep(double dt);
    void thermalStep(double dt);
    double wallTemperature(int i) const;
    // hydrate equilibrium temperature T_eq(p) from ln(p) = A + B / T_eq
    double equilibriumTemperature(double p) const
    {
        return cfg_.hydrate.curveB / (std::log(p) - cfg_.hydrate.curveA);
    }
    // liquid-flux function coefficients F(a) = qa a^2 + qb a + qc at the
    // current mixture velocity (algebraic closures; the dynamic closure
    // builds face-local coefficients from the slip state)
    void fluxCoeffs(double& qa, double& qb, double& qc) const;
    void gasFluxCoeffs(double& qa, double& qb, double& qc) const;
    double liquidFlux(double alpha, double s) const;
    double gasVolumeFlux(double alpha, double s) const;
    double gasVelocity(double alpha, double s) const;
    bool eosOn() const { return cfg_.fluid.gasR > 0; }
    double gasDensityOf(double p, double temperature) const
    {
        return p / (cfg_.fluid.gasR * temperature);
    }
    double sectionTemperature(int i) const
    {
        return cfg_.thermalOn ? T_[i] : cfg_.temperature0;
    }
    double cellFriction(int i, double vM) const; // sgn-agnostic |dp/dx| fric
    std::vector<double> cellPressureLosses(double vM) const; // fr + grav, Pa
    double mixtureRhoCell(int i) const; // with the per-section gas density
    double cellGasDensity(int i) const; // guarded: constant rhoG when sealed
    double ghostGasDensity(const TwoFluidBoundary& b) const;

    TwoFluidCase cfg_;
    int n_ = 0;
    double diameter_ = 0, length_ = 0; // uniform bore and total arc length
    std::vector<double> dxs_;    // per-cell arc step
    std::vector<double> sinT_;   // per-cell slope
    std::vector<double> relrough_;
    std::vector<double> centerDrop_; // vertical drop of each cell centre
                                     // from the chain inlet
    std::vector<double> alpha_, water_, T_;
    std::vector<double> slip_;    // dynamic closure: vG - vL per section
    std::vector<double> sDotOld_; // last step's ds/dt, the lagged slip
                                  // inertia the mixture momentum feels
    std::vector<double> gasMass_; // gas inventory per unit area kg/m2 =
                                  // rhoG (1 - alpha) dx; the transported
                                  // state when the EOS is on
    std::vector<double> liquidFaceFlux_; // this step's face liquid volume
                                         // flux (superficial, m/s); the
                                         // enthalpy carrier of the thermal
                                         // step
    std::vector<double> gasFaceFlux_;    // this step's face gas volume flux
    std::vector<double> hydrate_; // hydrate volume fraction of the slurry
    std::vector<double> hydrateDelta_; // this step's hydrate volume change
                                       // per bulk m3 (the latent-heat source
                                       // the thermal step consumes)
    double wallHeatRate_ = 0;
    double inletHeatRate_ = 0;
    double vM_ = 0;   // current mixture velocity
    double vMOld_ = 0;
    double t_ = 0;
};

/// Case grammar (see tests/cases/twoFluidFront/case.dict):
///   fluid { rhoL 1000; rhoG 20; muL 0.01; muG 1e-5; cL 4180; cG 2150;
///           gasR 287; }   // gasR > 0: ideal-gas EOS on, rhoG is the
///                         // gas density the injected rates refer to
///   pipe { length 100; diameter 0.1; roughness 0; angleDeg 0; sections 200; }
///         // or a serial chain of segments in flow order (uniform bore):
///         // pipes { riser { length 400; diameter 0.02; roughness 0;
///         //                 angleDeg 90; sections 40; }
///         //         flowline { length 400; diameter 0.02; roughness 0;
///         //                    angleDeg 0; sections 40; } }
///   closure { model wallis; slipVelocity 0.3; }   // or model zuber;
///                                                 //      C0 1.2; driftVelocity 0.15;
///                                                 // or model dynamic;
///                                                 //      interfacialDrag 500;
///   initial { alpha 0; water 0; temperature 300; fills { ... } }
///   boundaries
///   {
///       inlet { type phaseRates; vsl 1; vsg 0; waterFraction 0.25;
///               temperature 330; }
///       outlet { type fixedPressure; p 0; temperature 300; }
///   }
///   thermal      // optional: presence switches the temperature state on
///   {
///       wallHeatTransfer 5;      // W/(m K), 0 = insulated
///       wallModel constant;      // or geothermal;
///       wallTemperature 280;     // constant model
///       surfaceTemperature 288;  // geothermal: T(d) = surf + gradient*d
///       gradient 0.03;           // K/m of true vertical depth
///       inletDepth 2000;         // m below surface at s = 0
///   }
///   hydrate      // optional: presence switches the hydrate kinetics on;
///                // requires the gas EOS (gasR > 0) and a thermal{} block
///   {
///       curveA 35.2;             // ln(p[Pa]) = curveA + curveB / T_eq
///       curveB -5684;            // K
///       growthRate 0.01;         // 1/(s K), first order in free water
///       dissociationRate 0.01;   // 1/(s K), first order in h
///       gasContent 122;          // kg gas consumed per m3 hydrate
///       latentHeat 4.2e8;        // J released per m3 hydrate formed
///       initial 0;               // initial hydrate fraction of the slurry
///   }
///   time { dt 0.25; endTime 30; reportEvery 4; }
TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict);

} // namespace cf
