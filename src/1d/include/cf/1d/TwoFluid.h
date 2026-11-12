#pragma once

#include <cmath>
#include <map>
#include <memory>
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
/// The full closure (model full) solves the gas momentum per section
/// instead - the distributed pressure momentum form: the gas share of the
/// mixture row alpha_g dp (with the proportional wall split alpha_g fr
/// cancelling against the gas friction share) minus gas gravity minus the
/// interfacial drag drives the gas material acceleration
/// dvG/dt = du/dt + alpha ds/dt + vG dvG/dx, which with vG = vM + alpha s
/// reduces to
///   ds/dt = [(rho_m - rhoG)(g sinT + du/dt) - rhoG vG dvG/dx
///            - sigma/alpha_g]/(alpha rhoG),
/// terminal slip s* = sign(R) sqrt(alpha_g |R| / C) with
/// R = (rho_m - rhoG)(g sinT + du/dt): buoyancy, mixture acceleration and
/// the gas convection drive the slip, wall friction rides the pressure
/// field. At a uniform steady state the convective term vanishes identically
/// (steady fluxes are uniform, so rhoG and vG are).
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
/// reduced form (the model-full slip ODE carries the mixture-row pressure
/// gradient as its drive; a per-section mixture momentum is the follow-up).
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
/// With an entrainment{} block the sections also carry a droplet fraction d
/// of the liquid (the film-droplet split of annular flow; transport level
/// only in this increment - the closure, momentum and properties stay
/// mixture-based). The droplet content qd = alpha d rides the gas at its
/// velocity: the face flux qd vG upwinds on the gas direction, the film
/// flux F - qd vG completes the liquid flux F, so the volume balance
/// survives the split. d relaxes towards the equilibrium entrained
/// fraction E_eq = clamp((vG - vCrit)/vRef, 0, 1) with implicit Euler at
/// the exchange time relaxTime (vG from the closure's gas velocity). Inflow
/// boundaries carry entrainedFraction of the inflow liquid (rate faces
/// impose that share of the liquid rate; pressure ghosts inject
/// qd_ghost vG_ghost); junction mixing weights the droplet fraction with
/// the liquid volume flow like the water cut.

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
    // 5.75 (A5): the LIQUID phase's power-law rheology. The mixture's wall
    // stress becomes tau_w = alpha K' gamma^n + (1 - alpha) muG gamma on the
    // Newtonian wall rate gamma = 8 vm/d, with K' = K ((3n+1)/4n)^n the
    // Rabinowitsch-corrected consistency -- so at alpha = 1 the wall reads
    // back the exact single-phase power-law pipe law, at n = 1 (and with no
    // declaration at all) the Newtonian mixture law this row has always
    // had. The apparent mixture viscosity tau_w/gamma drives the same
    // Churchill blend. flowIndex 0 (the default) means Newtonian.
    // 5.77 (A5 increment 3): the liquid share of the wall can also carry a
    // yield stress tau_0, so tau_w gains + alpha tau_0 and the distributed
    // static-wall threshold becomes sum_i (4 dx_i/d) alpha_i tau_0. A
    // pressure-driven column offered less than that is HELD at u = 0 by
    // the latch 5.49 gave the single-phase row: the wall law is not a
    // function at u = 0, so its row is replaced instead of inverted. Zero
    // declaration keeps the old law verbatim (bit for bit).
    double flowIndex = 0;    // n; 0 = Newtonian liquid
    double consistency = 0;  // K, Pa s^n; must be positive when n > 0
    double yieldStress = 0;  // tau_0, Pa; 0 = no static wall
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
                               // 2 = dynamic (transported slip state),
                               // 3 = full gas momentum (transported slip
                               // state driven by the per-cell gas momentum
                               // balance: buoyancy + mixture acceleration
                               // against the interfacial drag; the wall
                               // friction rides the pressure field through
                               // the proportional split)
    double slipVelocity = 0;   // wallis us, m/s (positive = gas runs ahead)
    double c0 = 1.0;           // zuber profile parameter
    double driftVelocity = 0;  // zuber vD, m/s
    double interfacialDrag = 0;// dynamic/full: kg/m4, drag law
                               // sigma = C s |s| per unit volume;
                               // required > 0 for models 2 and 3
};

/// One piecewise-constant inlet-rate point: from t on the phaseRates
/// boundary imposes vsl/vsg until the next point; the boundary's scalar
/// rates govern before the first point. The loader resolves missing keys
/// to "hold the previous rate" and sorts the points by t (t must be > 0,
/// the constructor reads the scalar rates at t = 0).
struct TwoFluidSchedulePoint
{
    double t = 0;   // s, when this rate pair takes over
    double vsl = 0; // m/s, imposed liquid superficial velocity
    double vsg = 0; // m/s, imposed gas superficial velocity
};

struct TwoFluidBoundary
{
    enum class Type { sealed, phaseRates, fixedPressure };
    Type type = Type::sealed;
    double vsl = 0;          // phaseRates: imposed liquid superficial velocity
    double vsg = 0;          // phaseRates: imposed gas superficial velocity
    std::vector<TwoFluidSchedulePoint> schedule; // inlet rate schedule
                                                 // (piecewise constant,
                                                 // see TwoFluidSchedulePoint)
    double pressure = 0;     // fixedPressure: Pa
    double alpha = 0;        // fixedPressure inlet: composition of the
                             // inflow ghost cell
    double waterFraction = 0;// water cut of the inflow liquid
    double hydrateFraction = 0; // hydrate fraction of the inflow liquid
    double entrainedFraction = 0; // droplet fraction of the inflow liquid
                                  // (entrainment on): rate faces impose it
                                  // as that share of the liquid rate,
                                  // pressure ghosts inject it at the ghost
                                  // gas velocity
    double temperature = 0;  // K, temperature of the inflow ghost cell
                             // (required on both ends when thermal is on,
                             // the unused one for potential backflow)
    double rhoGInjected = 0; // gas density the boundary injects, kg/m3
                             // (junction mixing); 0 = the fluid reference
                             // rhoG. Consulted at rate-type faces by the
                             // gas mass flux and the enthalpy advection.
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

/// Slug-body detection over the transported holdup field: a statistics
/// layer, not a two-regime model. A body is a maximal run of consecutive
/// sections with alpha >= alphaThreshold whose total length reaches
/// minLength; slugBodies() reports each run's face positions, length and
/// liquid inventory. The detection rides the holdup the transport already
/// solves and adds no physics - it is the bookkeeping a slug-tracking
/// capability reports and verifies its models against.
struct TwoFluidSlug
{
    double alphaThreshold = 0.5; // liquid fraction a body needs
    double minLength = 0;        // m, shortest run that counts
};

struct TwoFluidSlugBody
{
    double rear = 0;      // m along the pipe (face coordinate)
    double front = 0;     // m
    double length = 0;    // m
    double inventory = 0; // liquid volume per m2 of bore, sum alpha dx
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

/// Film-droplet split of the liquid (annular flow, transport level only:
/// the closure, momentum and mixture properties stay untouched - the
/// droplet fraction rides the gas but does not feed back yet). E_eq is the
/// equilibrium entrained fraction the per-cell droplet fraction relaxes
/// towards: a linear ramp of the closure gas velocity between vCrit (no
/// entrainment below) and vCrit + vRef (full entrainment above), relaxed at
/// the exchange time constant relaxTime.
struct TwoFluidEntrainment
{
    double vCrit = 0;      // m/s, gas velocity where entrainment starts
    double vRef = 1;       // m/s, ramp width to full entrainment
    double relaxTime = 0;  // s, film <-> droplet exchange time constant
    double initial = 0;    // initial droplet fraction of the liquid
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
    // 5.77: a pressure-driven pipe may start already moving, so the
    // moving -> held crossing is observable without a time-varying
    // boundary. 0 (the default) is "start from rest", exactly as before.
    double vm0 = 0;       // initial mixture velocity, m/s
    std::vector<TwoFluidFill> fills; // sub-range overrides of alpha0
    bool thermalOn = false;      // energy equation solved (thermal{} present)
    double temperature0 = 300;   // K, uniform initial temperature
    TwoFluidThermal thermal;
    bool hydrateOn = false;      // hydrate kinetics solved (hydrate{} present)
    TwoFluidHydrate hydrate;
    bool entrainmentOn = false;  // film-droplet split solved (entrainment{}
                                 // present)
    TwoFluidEntrainment entrainment;
    bool slugOn = false;         // slug-body statistics solved (slug{} present)
    TwoFluidSlug slug;
    bool networkMode = false;    // branch-network member: the junction end
                                 // boundary is owned by the TwoFluidNetwork
                                 // (mutated every step), so the constructor
                                 // skips the serial boundary-combination
                                 // checks
};

/// One junction of a branch network: the named pipe ends it joins, each
/// "pipeName.inlet" or "pipeName.outlet".
struct TwoFluidJunctionSpec
{
    std::string name;
    std::vector<std::string> ends;
};

/// Branch-network case: several TwoFluidPipe objects sharing junction
/// nodes. The per-pipe physics (fluid/closure/initial/thermal/hydrate) is
/// shared; geometry and the non-junction boundary of each pipe are per
/// pipe. Every pipe joins one or two junctions (a two-junction pipe is a
/// span between nodes and carries no boundary of its own); the per-step
/// hydraulics solve as one global Newton over the whole network.
struct TwoFluidNetworkCase
{
    TwoFluidCase perPipe; // networkMode = true; pipes/inlet/outlet unused
    std::vector<TwoFluidGeo> pipes;
    std::vector<TwoFluidJunctionSpec> junctions;
    std::map<std::string, TwoFluidBoundary> boundaries; // "pipe.inlet"
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
    /// 5.77 (A5 increment 3): 1 when the static wall held this pipe at
    /// rest on the step just taken (its momentum row was replaced by
    /// u = 0), 0 when the moving wall law was solved. Only pressure-driven
    /// pipes can hold: a rate boundary imposes its velocity, a sealed pipe
    /// is static by construction. At t = 0 no step has been taken, so the
    /// flag reads 0 (the constructor state) even for a held case.
    bool staticHold() const { return held_; }
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
    /// Droplet fraction of the liquid per section (entrainment on);
    /// zeros unless the entrainment{} block is present.
    std::vector<double> cellEntrained() const { return entrained_; }
    /// Slug bodies in the current holdup field: maximal runs of sections
    /// with alpha >= alphaThreshold and length >= minLength, in flow
    /// order (empty unless the slug{} block is present). See TwoFluidSlug.
    std::vector<TwoFluidSlugBody> slugBodies() const;
    /// Wall heat the sections absorb this step, W (positive = into the fluid);
    /// zero unless thermal is on.
    double wallHeatRate() const { return wallHeatRate_; }
    /// Enthalpy the inflow boundaries inject, W (transport + sources only).
    double inletHeatRate() const { return inletHeatRate_; }

    // ---- junction support (used by TwoFluidNetwork) ----

    /// Bore area, m2.
    double area() const { return M_PI * diameter_ * diameter_ / 4.0; }
    /// Impose the mixture velocity for the next step (the network's
    /// coupled junction Newton owns the momentum equation of junction
    /// branches; the flag is consumed by the next step()).
    void setExternalMomentum(double vM)
    {
        externalMomentum_ = true;
        vM_ = vM;
    }
    /// Overwrite one end's boundary. The network rewrites the junction
    /// end every step: a fixed-pressure ghost at the mixed state for
    /// branches feeding the junction, or the imposed junction split
    /// (phaseRates with the mixed composition) for branches it supplies.
    void setInletBoundary(const TwoFluidBoundary& b) { cfg_.inlet = b; }
    void setOutletBoundary(const TwoFluidBoundary& b) { cfg_.outlet = b; }
    /// Fixed-pressure-style momentum residual at trial velocity u with
    /// the caller's signed driving difference folded in:
    ///   inertia (u - uOld)/dt - slipInertia + fr(u) + grav - dpDrive
    /// on the frozen holdup. The junction Newton assembles one row per
    /// pressure-driven branch from this.
    double momentumResidual(double dt, double u, double dpDrive) const;
    /// Natural magnitude of that residual (Newton tolerance anchor):
    /// inertia (|uOld| + 1)/dt + |grav| + 1.
    double momentumScale(double dt) const;

    // ---- distributed static wall (5.77 / A5 increment 3) ----

    /// The drop the yielded wall can hold without moving, Pa:
    ///   sum_i (4 dx_i / d) alpha_i tau_0.
    /// Only the liquid share of each cell's wall yields, so a gas cell
    /// contributes nothing and a non-uniform holdup field weights each
    /// section by its own arc length. At uniform alpha this is exactly the
    /// single-phase 4 tau_0 L / d of 5.49 scaled by alpha. Geometry and
    /// holdup only -- the threshold knows neither n nor K.
    double yieldThresholdPa() const;
    /// Static gravity head sum_i rho_m,i g sinT_i dx_i, Pa. The drive the
    /// static wall has to resist is dpDrive - this, the same quantity the
    /// momentum residual moves.
    double gravityHeadPa() const;
    /// Latch the held set ONCE at the start of a step (5.49): the decision
    /// reads the old velocity and the drive offered THEN, not the node
    /// pressures the same solve is moving. dpDrive = p(inlet) - p(outlet).
    void latchStaticHold(double dpDrive);
    /// A moving column under a sub-threshold drive has NO root on its own
    /// branch: friction is bounded below by the threshold there, so once
    /// the old velocity is small enough even u -> 0 cannot balance the
    /// drive. It parks at u = 0 for this step; the next step's latch
    /// (old velocity now zero) owns it from there. This is 5.49's
    /// `crossing` stated on the RESIDUAL, because the yielded wall law is
    /// not a function at u = 0 and a Newton differencing through zero
    /// would straddle the threshold jump (5.49's own frictionSlope note).
    bool yieldCrossing(double dt, double u, double dpDrive) const;
    /// Rate-driven and sealed pipes never hold: clear the flag so a pipe
    /// whose boundary is rewritten by a network cannot inherit a stale
    /// latch.
    void clearStaticHold() { held_ = false; }
    /// State at the given end (atInlet = face 0, else face n): the
    /// boundary-cell holdup/water/hydrate/temperature/gas density the
    /// junction mixing draws from, and the pipe's own phase volume
    /// fluxes there (feeding branches deliver these to the junction).
    double endAlpha(bool atInlet) const
    {
        return atInlet ? alpha_[0] : alpha_[n_ - 1];
    }
    double endWater(bool atInlet) const
    {
        return atInlet ? water_[0] : water_[n_ - 1];
    }
    double endHydrate(bool atInlet) const
    {
        return cfg_.hydrateOn ? (atInlet ? hydrate_[0] : hydrate_[n_ - 1])
                              : 0.0;
    }
    double endEntrained(bool atInlet) const
    {
        return cfg_.entrainmentOn
                   ? (atInlet ? entrained_[0] : entrained_[n_ - 1])
                   : 0.0;
    }
    double endTemperature(bool atInlet) const
    {
        return cfg_.thermalOn ? (atInlet ? T_[0] : T_[n_ - 1])
                              : cfg_.temperature0;
    }
    double endGasDensity(bool atInlet) const
    {
        return atInlet ? cellGasDensity(0) : cellGasDensity(n_ - 1);
    }
    double endLiquidFlux(bool atInlet) const
    {
        return liquidFlux(atInlet ? alpha_[0] : alpha_[n_ - 1], 0.0);
    }
    double endGasVolumeFlux(bool atInlet) const
    {
        return gasVolumeFlux(atInlet ? alpha_[0] : alpha_[n_ - 1], 0.0);
    }
    /// The pipe's phase volume fluxes at an arbitrary holdup with its own
    /// (signed) mixture velocity - the junction split a supplied branch
    /// receives: F(aJ), G(aJ) with F + G = vM for every closure, so the
    /// volume balance survives the split exactly.
    double liquidFluxAt(double alpha) const { return liquidFlux(alpha, 0.0); }
    double gasVolumeFluxAt(double alpha) const
    {
        return gasVolumeFlux(alpha, 0.0);
    }

private:
    void momentumSolve(double dt);
    void slipStep(double dt);
    void transportStep(double dt);
    void entrainmentStep(double dt);
    void hydrateStep(double dt);
    void thermalStep(double dt);
    /// Evaluate the inlet rate schedule at the step midpoint (robust to
    /// t_ accumulation drift) into the boundary's scalar rates (no-op
    /// without a schedule).
    void applySchedule(double dt);
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
    // closures with a transported slip state share the face-local flux
    // machinery (models 2 and 3 differ only in the slip ODE)
    bool slipClosure() const
    {
        return cfg_.closure.model == 2 || cfg_.closure.model == 3;
    }
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
    /// Momentum residual in the u -> 0 limit on the branch the column
    /// came FROM: friction tends to the signed distributed threshold, so
    /// this is the exact rest value of the discontinuous wall law (5.77).
    double restResidualPa(double dt, double dpDrive) const;
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
    std::vector<double> entrained_; // droplet fraction of the liquid
                                    // (entrainment on)
    double wallHeatRate_ = 0;
    double inletHeatRate_ = 0;
    double vM_ = 0;   // current mixture velocity
    double vMOld_ = 0;
    double t_ = 0;
    bool held_ = false;       // 5.77: the static wall owns this step's row
    double driveLatch_ = 0;   // 5.77: dpDrive - gravityHead at latch time
    bool externalMomentum_ = false; // one-shot: the next step() skips
                                    // momentumSolve (junction branches)
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
///                                                 // or model full;
///                                                 //      interfacialDrag 500;
///   initial { alpha 0; water 0; temperature 300; vm 0; fills { ... } }
///         // vm: optional initial mixture velocity for a pressure-driven
///         // pipe (0 = start from rest); useful to observe a moving column
///         // stop under a sub-threshold drive (5.77)
///   boundaries
///   {
///       inlet { type phaseRates; vsl 1; vsg 0; waterFraction 0.25;
///               temperature 330;
///               schedule                    // optional inlet rate
///               {                           // schedule (phaseRates only):
///                   p30 { t 30; vsl 0.4; vsg 0; }   // from t on these
///                   p35 { t 35; vsl 1;   vsg 0; }   // rates hold until
///               }                                   // the next point;
///       }                                           // missing keys hold
///       outlet { type fixedPressure; p 0; temperature 300; }  // the
///   }                                               // previous rate
///   slug         // optional: slug-body statistics over the holdup field
///   {
///       alphaThreshold 0.3;      // liquid fraction a body needs
///       minLength 1;             // m, shortest run that counts
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
///   entrainment  // optional: presence switches the film-droplet split on
///   {
///       vCrit 0.2;               // m/s, gas velocity where entrainment
///                                // starts
///       vRef 0.6;                // m/s, ramp width to full entrainment
///       relaxTime 8;             // s, film<->droplet exchange time
///       initial 0;               // initial droplet fraction of the liquid
///   }
///   time { dt 0.25; endTime 30; reportEvery 4; }
TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict);

/// Branch network of TwoFluidPipe objects joined at junction nodes. Each
/// pipe joins one or two junctions: a single-junction pipe has a real
/// boundary at its other end, a two-junction pipe is a span between nodes.
/// Every step the hydraulics solve as ONE global Newton: the unknowns are
/// the pressure-driven pipes' mixture velocities (fixed-pressure far end,
/// or span) plus every junction pressure; the rows are the implicit-Euler
/// momentum residual per driven pipe (analytic sensitivity to the junction
/// pressures at its ends) and the volume balance sum(s A vM) = 0 per
/// junction. Rate-driven far ends fix their branch velocity, sealed far
/// ends are dead legs; each junction needs at least one fixed-pressure far
/// end as its pressure reference. After the solve every junction mixes its
/// inflows (holdup/water/hydrate/droplet volume-flow-weighted, temperature
/// heat-capacity-flow-weighted, gas density mass-flux-weighted) and every
/// junction face is set: feeding branches get a fixed-pressure ghost at
/// the mixed state, supplied branches get the junction split imposed as
/// phaseRates with the pipe's own flux functions at the mixed holdup
/// (F + G = vM keeps the volume balance exact; the phase split is exact
/// for the homogeneous closure and closure-approximate with slip). A span
/// pipe can feed one junction and be supplied by the other. The dynamic
/// closure is not supported yet.
class TwoFluidNetwork
{
public:
    explicit TwoFluidNetwork(const TwoFluidNetworkCase& c);

    /// One network step: solve every junction (branch velocities + node
    /// pressure + face state), then advance every pipe.
    void step(double dt);

    double time() const { return t_; }
    std::vector<std::string> pipeNames() const;
    const TwoFluidPipe& pipe(const std::string& name) const;
    std::vector<std::string> junctionNames() const;
    /// Last solved junction pressure, Pa.
    double junctionPressure(const std::string& name) const;
    /// Signed net mixture volume rate into the junction, m3/s (the
    /// balance the Newton closed; ~0 at machine precision).
    double junctionNetRate(const std::string& name) const;
    /// Signed mixture volume rate of each connected pipe INTO the
    /// junction, m3/s (connection order: positive = the pipe feeds the
    /// junction).
    std::vector<double> junctionPipeRates(const std::string& name) const;

private:
    struct EndRef
    {
        int pipe = 0;            // index into pipes_
        bool atOutlet = false;   // the junction sits at the pipe's outlet
        TwoFluidBoundary far;    // the pipe's non-junction boundary
                                 // (sealed for a span pipe's junction end:
                                 // the other end is another junction)
    };
    struct Junction
    {
        std::string name;
        std::vector<EndRef> ends;
        double p = 0;            // last solved junction pressure
        double net = 0;          // last balance residual (diagnostic)
        std::vector<double> rates; // signed per-pipe rates into the node
        // last mixed junction-face state (kept when nothing flows in)
        double alphaJ = 0, waterJ = 0, hydrateJ = 0, entrJ = 0, tJ = 0,
               rhoJ = 0;
    };
    void solveNetwork(double dt);

    TwoFluidNetworkCase cfg_;
    std::vector<std::unique_ptr<TwoFluidPipe>> pipes_;
    std::map<std::string, int> pipeIndex_;
    std::vector<Junction> junctions_;
    double t_ = 0;
};

/// Network case grammar: the shared physics blocks are the serial ones
/// (fluid/closure/initial/thermal/hydrate); the topology lives in a
/// network{} block:
///   network
///   {
///       pipes { branchA { length 10; diameter 0.02; angleDeg 0;
///                         sections 10; } ... }
///       junctions { J { ends ( branchA.outlet branchB.outlet
///                             outlet.inlet ); } }
///       boundaries
///       {
///           branchA.inlet { type phaseRates; vsl 0.5; vsg 0;
///                           temperature 300; }
///           outlet.outlet { type fixedPressure; p 0; temperature 300; }
///       }
///   }
/// (boundary keys name the pipe end; ends claimed by a junction must not
/// carry a boundary, every other end must. A pipe may join two junctions -
/// a span between nodes with both ends boundary-free:
///   J1 { ends ( srcA.outlet srcB.outlet mid.inlet ); }
///   J2 { ends ( mid.outlet sink1.inlet sink2.inlet ); } )
TwoFluidNetworkCase loadTwoFluidNetworkCase(const Dictionary& caseDict);

} // namespace cf
