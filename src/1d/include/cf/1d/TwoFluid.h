#pragma once

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
///   rho L/A du/dt = dp - sum(friction) - sum(rho_m g sin(theta) dx)
///                                                    [mixture momentum]
/// with the drift-flux slip closure vG = C0 vM + vD (Zuber-Findlay) or the
/// constant-slip form vG = vM + alpha us (Wallis). Phases are incompressible
/// and the area is constant, so the mixture velocity is uniform along the
/// pipe; per-section gas compressibility (which makes the momentum equation
/// itself distributed) is the next increment. Holdup transport advances with
/// the exact Godunov flux of the quadratic liquid-flux function
/// F(alpha) = alpha vL(alpha), so shocks propagate at the Rankine-Hugoniot
/// speed and the CFL-limited scheme is monotone (alpha stays in [0, 1]).
/// Wall friction is the homogeneous Churchill law on mixture properties;
/// gravity uses the local mixture density.
/// With a thermal{} block the sections also carry a mixture temperature (the
/// phases stay in thermal equilibrium): implicit Euler on
///   C_i dT_i/dt = H_in,i - H_out,i + k' dx (T_wall,i - T_i),
/// C_i the section's phase-weighted rho*c capacity and the enthalpy flux
/// riding the mixture velocity with the upwind cell's composition (the same
/// homogeneous closure as the friction law). Implicit in T, so no CFL on top
/// of the transport one.

struct TwoFluidFluid
{
    double rhoL = 1000; // kg/m3
    double rhoG = 20;   // kg/m3
    double muL = 0.01;  // Pa s
    double muG = 1e-5;  // Pa s
    double cL = 0;      // J/(kg K), needed when the thermal block is present
    double cG = 0;      // J/(kg K)
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
                               // 1 = zuber (vG = C0 vM + vD)
    double slipVelocity = 0;   // wallis us, m/s (positive = gas runs ahead)
    double c0 = 1.0;           // zuber profile parameter
    double driftVelocity = 0;  // zuber vD, m/s
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

struct TwoFluidCase
{
    TwoFluidFluid fluid;
    TwoFluidGeo geo;
    TwoFluidClosure closure;
    TwoFluidBoundary inlet;
    TwoFluidBoundary outlet;
    double alpha0 = 0;    // uniform initial liquid fraction
    double water0 = 0;    // uniform initial water cut of the liquid
    std::vector<TwoFluidFill> fills; // sub-range overrides of alpha0
    bool thermalOn = false;      // energy equation solved (thermal{} present)
    double temperature0 = 300;   // K, uniform initial temperature
    TwoFluidThermal thermal;
};

class TwoFluidPipe
{
public:
    explicit TwoFluidPipe(const TwoFluidCase& c);

    /// One time step: implicit-Euler momentum for the mixture velocity
    /// (frozen holdup), then explicit Godunov transport of alpha and the
    /// water cut. Fatals when dt violates the CFL condition.
    void step(double dt);

    double time() const { return t_; }
    int sections() const { return n_; }
    double cellSize() const { return dx_; }

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
    std::vector<double> cellTemperature() const { return T_; }
    /// Wall heat the sections absorb this step, W (positive = into the fluid);
    /// zero unless thermal is on.
    double wallHeatRate() const { return wallHeatRate_; }
    /// Enthalpy the inflow boundaries inject, W (transport + sources only).
    double inletHeatRate() const { return inletHeatRate_; }

private:
    void momentumSolve(double dt);
    void transportStep(double dt);
    void thermalStep(double dt);
    double wallTemperature(int i) const;
    // liquid-flux function coefficients F(a) = qa a^2 + qb a + qc at the
    // current mixture velocity
    void fluxCoeffs(double& qa, double& qb, double& qc) const;
    double liquidFlux(double alpha) const;
    double gasVelocity(double alpha) const;
    double cellFriction(int i, double vM) const; // sgn-agnostic |dp/dx| fric
    std::vector<double> cellPressureLosses(double vM) const; // fr + grav, Pa

    TwoFluidCase cfg_;
    int n_ = 0;
    double dx_ = 0, sinTheta_ = 0;
    std::vector<double> alpha_, water_, T_;
    double wallHeatRate_ = 0;
    double inletHeatRate_ = 0;
    double vM_ = 0;   // current mixture velocity
    double vMOld_ = 0;
    double t_ = 0;
};

/// Case grammar (see tests/cases/twoFluidFront/case.dict):
///   fluid { rhoL 1000; rhoG 20; muL 0.01; muG 1e-5; cL 4180; cG 2150; }
///   pipe { length 100; diameter 0.1; roughness 0; angleDeg 0; sections 200; }
///   closure { model wallis; slipVelocity 0.3; }   // or model zuber;
///                                                 //      C0 1.2; driftVelocity 0.15;
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
///   time { dt 0.25; endTime 30; reportEvery 4; }
TwoFluidCase loadTwoFluidCase(const Dictionary& caseDict);

} // namespace cf
