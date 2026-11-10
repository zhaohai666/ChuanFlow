#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"
#include "cf/pvt/BlackOil.h"

namespace cf
{

/// Constant-property single-phase fluid for the 1D network kernel.
/// Property-varying fills use PvtFluid below.
struct FluidProperties
{
    double rho = 850.0; // kg/m3
    double mu = 1e-3;   // Pa.s
};

/// Tabulated single-phase fluid: rho(p) and mu(p) on a strictly ascending
/// pressure grid, piecewise-linear between nodes and clamped outside the
/// table ends. The steady pipe law marches the local mechanical-energy
/// balance with properties evaluated at the local pressure, so a
/// compressible fill (e.g. an isothermal gas, rho = p/(R T)) contributes
/// its shrinking density to the gravity head exactly where it acts.
struct PvtFluid
{
    std::vector<double> pressure; // Pa, strictly ascending
    std::vector<double> rho;      // kg/m3 at those pressures
    std::vector<double> mu;       // Pa.s at those pressures

    double density(double p) const;
    double viscosity(double p) const;
};

struct NetworkPipe
{
    std::string name;
    std::string fromName;
    std::string toName;
    label from = -1; // resolved by OneDNetwork::addPipe
    label to = -1;
    double length = 0;     // m, along the pipe
    double diameter = 0;   // m
    double roughness = 0;  // m, absolute roughness
    double elevation = 0;  // m, z(to) - z(from); adds rho*g*dz to the drop

    // inline equipment (single-phase pipes; the black-oil march ignores it)
    double valveK = 0;         // valve loss coefficient: dp = K rho v^2/2,
                               // always opposing the flow
    double pumpHead = 0;       // m of head at pumpRate (> 0 switches the
                               // pump on); constant head when pumpRate = 0
    double pumpShutoffHead = 0;// m at zero rate (defaults to pumpHead)
    double pumpRate = 0;       // m3/s where pumpHead is quoted

    // Power-law rheology (Metzner & Reed 1955): when `flowIndex` > 0 the
    // pipe's fluid is tau = K gdot^n instead of Newtonian, and the single
    // change to the Darcy-Weisbach law is the Reynolds number in the
    // friction-factor call -- the generalized one
    //     Re_MR = rho u^(2-n) d^n / (K' 8^(n-1)),  K' = K ((3n+1)/(4n))^n
    // which reduces to rho u d / mu *exactly* at n = 1, K = mu, so a
    // declared-but-Newtonian pipe is bit-identical to an undeclared one.
    // With this substitution Churchill's laminar branch reads 64/Re_MR,
    // which is the exact Hagen-Poiseuille analogue for a power-law fluid,
    // so the laminar drop is not a correlation but a closed form. The
    // exponent on (3n+1)/(4n) is what makes it exact: that factor is
    // Rabinowitsch's wall shear rate over the apparent one, and the stress
    // the profile holds is K times it to the n, not to the 1. The
    // textbook's exponent-1 constant is a turbulent-collapse convention and
    // sits ((3n+1)/(4n))^(1-1/n) off the laminar solution (25% in the rate at
    // n = 0.5); 5.47 measures that and fixes it.
    double flowIndex = 0;      // n; 0 = Newtonian (the fluid's own mu)
    double consistency = 0;    // K, Pa.s^n; required when flowIndex > 0

    // Bingham plastic (tau = tau_0 + mu_p gdot): `yieldStress` is tau_0 in Pa
    // and mu_p is the fluid's mu. This is NOT another substitution into the
    // Reynolds number, because the drop no longer vanishes with the rate: at
    // q = 0 the wall still holds dpStatic = 4 tau_0 L / d, and below that the
    // exact answer is q = 0 (the pipe is stuck, and it takes whatever drop is
    // asked of it up to the threshold). Above it the law is Buckingham's 1921
    // closed form, inverted for the wall stress, and the reported friction
    // factor is the exact f = 8 tau_w / (rho v^2) rather than a Churchill
    // number -- which makes f Re g(xi) = 64 a restatement of the solution.
    // tau_0 -> 0 recovers Hagen-Poiseuille term for term, so the undeclared
    // value is a sentinel, not a branch.
    //
    // Declaring `yieldStress` *together with* a `flowIndex` is
    // Herschel-Bulkley, tau = tau_0 + K gdot^n, and the same two keywords are
    // read as its threshold and its consistency. Its laminar rate is one
    // termwise integrable expression in the wall stress,
    //     Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n)
    //         [ tau_0^2/(n+1) + 2 tau_0 D/(2n+1) + D^2/(3n+1) ],  D = tw - tau_0,
    // which is Buckingham's bracket at n = 1, K = mu and the power-law closed
    // form above at tau_0 = 0 -- so both shipped laws are limits of this one,
    // and each is asserted against its own closed form. mu does not enter it,
    // which is why a Herschel-Bulkley pipe reports the *apparent* Reynolds
    // number rho u d / (tw d / 8u) rather than a material one.
    double yieldStress = 0;    // Pa; 0 = no threshold

    /// Coupled element (the 1D side of a 1D<->3D interface): when > 0 this
    /// pipe carries NO 1D resistance of its own – its whole loss is
    /// `coupledResistance * |Q|`, supplied by the partner model, and the
    /// friction/valve terms are skipped so nothing is counted twice. The
    /// coupling layer re-solves the network with an updated value until the
    /// drop the network sees equals the drop the detailed model reports.
    /// Length and diameter stay meaningful: they are the detailed region's own
    /// geometry, used for the reported velocity and Reynolds number.
    double coupledResistance = 0; // Pa/(m3/s), i.e. Pa.s/m3
};

/// Steady single-phase node-model pipeline network (the M1 "SteadySolver"
/// core): unknowns are node pressures, each pipe carries a flow-rate/pressure
/// drop relation with Darcy-Weisbach friction (Churchill all-regime friction
/// factor) plus the gravity head. Solved by damped Newton.
class OneDNetwork
{
public:
    std::vector<std::string> nodeNames;
    std::vector<double> nodePressureFixed; // NaN for unknown (free) nodes
    std::vector<double> nodeExternalRate;  // m3/s injected at the node
    FluidProperties fluid;
    std::vector<NetworkPipe> pipes;
    // when set, every pipe is solved as an isothermal black-oil two-phase
    // Beggs-Brill pipe and `sources` rates are stock-tank oil rates
    // (sm3/s); `fluid` is then unused
    bool hasBlackOil = false;
    pvt::BlackOilFluid blackOil;
    // tabulated-property single-phase mode (ignored when hasBlackOil is
    // set): `sources` rates are volumetric m3/s at the upstream end of
    // each pipe
    bool hasPvtFluid = false;
    PvtFluid pvtFluid;

    label addNode(const std::string& name);
    label nodeIndex(const std::string& name) const; // cfFatal if missing
    label addPipe(const NetworkPipe& pipe);         // cfFatal on bad nodes

    void setNodePressure(const std::string& node, double p)
    {
        nodePressureFixed[nodeIndex(node)] = p;
    }
    void setNodeRate(const std::string& node, double q)
    {
        nodeExternalRate[nodeIndex(node)] += q;
    }

    std::size_t nNodes() const { return nodeNames.size(); }
};

struct SteadyControls
{
    int maxIter = 100;
    double tol = 1e-9;   // relative to total injected rate
    double relax = 1.0;  // Newton damping
};

/// The published laminar-turbulent transition locus of a yielded pipe. Wilson &
/// Thomas (Can. J. Chem. Eng. 84:520, 2006) show that for a Bingham plastic the
/// condition at transition depends on the Hedstrom number alone,
///     He = tau_0 d^2 rho / mu^2,      Re_B = rho |v| d / mu,
/// with Re_B,t = 2100 at low He, 80 He^0.40 in the mid band (their Trel form is
/// 80 He^-0.10) and 25 He^0.5 above He = 1e5, where
///     Trel = V / (tau_0/rho)^0.5 = Re_B / He^0.5
/// is an identity, not a correlation. `band` is which of the three printed
/// ranges the pipe's He falls in (0: He <= 1700, 1: <= 1e5, 2: above), and stays
/// -1 for any pipe whose law is not Bingham's -- a Herschel-Bulkley pipe
/// declares no plastic viscosity, so neither He nor Re_B is a statement about
/// it. `dpAtTransition` is what makes the locus usable: it is the drop the
/// pipe's OWN wall law (5.46/5.47) asks for at the locus rate, i.e. how far the
/// injection end has to move before the laminar statement stops being the
/// one that applies. Nothing here re-prices the friction.
struct YieldedTransition
{
    int band = -1;
    double hedstrom = 0;            // He
    double rePlastic = 0;           // Re_B at the solved rate
    double reTransition = 0;        // Re_B,t(He)
    double rateAtTransition = 0;    // m3/s corresponding to Re_B,t
    double dpAtTransition = 0;      // Pa the wall law asks at that rate
    double relTransitionSpeed = 0;  // Trel = V_t / (tau_0/rho)^0.5
    bool pastLaminarLimit = false;  // the solved rate is past the locus
};

struct PipeFlowResult
{
    std::string name;
    double dpTotal = 0;     // p(from) - p(to), Pa
    double dpFriction = 0;  // dpTotal - rho*g*dz (equipment included; the
                            // pipe friction alone is dpFriction - dpValve
                            // - dpPump)
    double dpValve = 0;     // valve loss, K rho |v| v / 2 (opposes the flow)
    double dpPump = 0;      // pump contribution, -rho*g*H(q): negative when
                            // the pump adds head in the from -> to direction
    double velocity = 0;    // m/s, signed from -> to
    double flowRate = 0;    // m3/s
    double reynolds = 0;
    double frictionFactor = 0;
    double conductance = 0; // dQ/d(dpTotal), s
    YieldedTransition transition; // the locus, for a Bingham pipe only
    double holdup = 0;      // two-phase: outlet liquid volume fraction
    int regime = -1;        // two-phase: Beggs-Brill regime of the last
                            // segment (0 segregated, 1 transition,
                            // 2 intermittent, 3 distributed; -1 single-phase)
    double massRate = 0;    // two-phase: kg/s, invariant along the pipe
    bool choked = false;    // two-phase: the march hit the pressure floor at
                            // the imposed end pressures
};

struct SteadyNetworkResult
{
    std::vector<double> nodePressure;
    std::vector<PipeFlowResult> pipeFlows;
    int iterations = 0;
    double residual = 0;      // max |node residual| in m3/s
    double relativeResidual = 0;
    bool converged = false;
};

/// Churchill (1977) explicit all-regime friction factor. Re must be > 0;
/// returns the laminar limit 64/Re for Re -> 0.
double churchillFriction(double reynolds, double relativeRoughness);

/// The friction drop (Pa, positive) the declared wall law resists at a rate
/// magnitude |q|: Newtonian or power-law through Churchill on the
/// Metzner-Reed number, a yielded pipe through its own wall stress. Declared
/// here because the transient momentum row and the steady pipe loss are two
/// callers of one law; a rheology read differently in the two kernels would be
/// a third set of assertions nobody could reconcile.
double pipeFrictionPa(const NetworkPipe& pipe, const FluidProperties& fluid,
                      double q);

/// The transition locus of one yielded pipe at a rate magnitude |q| (see
/// YieldedTransition). Reported, never charged: the friction a yielded pipe
/// resists is still its laminar wall law, because the published turbulent-yield
/// closure is not landed yet.
YieldedTransition yieldedTransitionOf(const NetworkPipe& pipe,
                                      const FluidProperties& fluid, double q);

/// Signed flow through one pipe for a given total pressure difference
/// dp = p(from) - p(to). Q/dpTotal conductance via central differencing.
void solvePipeFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                   double dpTotal, PipeFlowResult& out);

/// Signed flow through one pipe for given absolute end pressures with
/// tabulated rho(p)/mu(p): the rate is solved by a bracketed secant on a
/// segment-wise march of the local energy balance (mass rate fixed by the
/// upstream density, so the fill expands as the pressure falls), and the
/// conductance comes from a one-sided probe into the live branch — the
/// static-head knee leaves a Q = 0 plateau, as in the two-phase law.
/// Flow always runs from the high-pressure end to the low-pressure end;
/// Q is signed from -> to.
bool solvePipeFlowPvt(const NetworkPipe& pipe, const PvtFluid& fluid,
                      double pFrom, double pTo, PipeFlowResult& out);

SteadyNetworkResult solveSteadyNetwork(const OneDNetwork& net,
                                       const SteadyControls& ctl);

/// Beggs-Brill pipe law Q(dpTotal) with absolute end pressures (the black-oil
/// PVT needs the pressure level, not just the difference): flow runs from the
/// high-pressure end to the low-pressure end, the rate is solved by a
/// bracketed secant on the march, and the conductance dQ/d(dpTotal) comes
/// from a central difference. Q is signed from -> to. Returns false when the
/// pipe chokes (the march hits the pressure floor before reaching pLow).
bool solvePipeFlowTwoPhase(const NetworkPipe& pipe,
                           const pvt::BlackOilFluid& fluid, double pFrom,
                           double pTo, PipeFlowResult& out);

/// Case dictionary grammar (see tests/cases/steadyNetwork/case.dict):
///   fluid { rho 850; mu 0.1; }
///   nodes { W1 {}; Outlet { pressure 1.0e6; } }
///   pipes { P1 { from W1; to Outlet; length 1000; diameter 0.05;
///                roughness 0; elevation 5;
///                valve 2.5;                 // loss coefficient K
///                pumpHead 25;               // m of head at pumpRate
///                pumpRate 0.02;             // 0 = constant-head pump
///                pumpShutoffHead 30;        // defaults to pumpHead
///                coupledResistance 3.8e6;   // Pa/(m3/s): the 1D side of a
///                                           // 1D<->3D interface (see above)
///              } }
///   sources { s1 { node W1; rate 0.002; } }
///   solve { maxIter 100; tol 1e-10; relax 1.0; }
/// Tabulated-property variant: fluid { pvt { pressure ( 1e5 1e6 4e6 );
///   rho ( 5 25 100 ); mu ( 1e-5 1.15e-5 1.5e-5 ); } } — piecewise-linear
/// rho(p)/mu(p), clamped at the table ends.
OneDNetwork loadNetworkFromDict(const Dictionary& caseDict);
SteadyControls loadSteadyControls(const Dictionary& caseDict);

} // namespace cf
