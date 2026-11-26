#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"
#include "cf/common/Rheology.h"
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

/// One station of a live-flash PVT march (D10/D7, §5.121): the state a true
/// tpFlash + LBC + Weinaug-Katz evaluation produced at the local pressure.
/// Empty on every table/constant/black-oil path, so no existing readout
/// changes. Populated only on the opt-in per-station flash march.
struct PvtFlashStation
{
    double p = 0;      // Pa, station pressure
    double rho = 0;    // kg/m3, homogeneous flash-mixture density
    double mu = 0;     // Pa s, mass-fraction-weighted phase LBC viscosity
    double sigma = 0;  // N/m, per-station Weinaug-Katz surface tension
};

/// Opaque live-flash state (composition, T, parachors, sigma switch). Defined
/// in network.cpp so this header does not need the pvt flash headers; held by
/// shared_ptr so PvtFluid keeps value semantics.
struct LiveFlashState;

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
    // 5.144 (G6): the stream's inlet temperature (K) for the PASSIVE
    // thermal march. 0 = undeclared: the table path never reads it and
    // every existing number is bit for bit; a pipe with thermal physics
    // on a pvt fluid without this declaration is refused.
    double temperature = 0;

    double density(double p) const;
    double viscosity(double p) const;
    /// D10/D7: per-station Weinaug-Katz sigma. Returns 0 on the table path
    /// and whenever the case did not opt into sigmaFromP, so an existing
    /// case cannot observe this method.
    double surfaceTension(double p) const;
    /// Evaluate rho/mu/sigma at p from ONE flash when liveFlash is set (the
    /// march calls this once per station); otherwise table lookup + sigma 0.
    /// 5.194 (G6d): the optional T is the STATE temperature -- the flash, and
    /// the LBC viscosity that rides on the flash's own molar volume, are
    /// evaluated at (T, p) instead of at the feed's single temperature. T <= 0
    /// means "the feed's temperature", which is the pre-feedback path bit for
    /// bit and the only state a single-temperature table has.
    void evaluateAt(double p, double& rhoOut, double& muOut,
                    double& sigmaOut, double T = 0.0) const;

    /// 5.189: the compositional feed, when the case declared one. Held even
    /// on the table path (where liveFlash is null) so an EOS-derived
    /// Joule-Thomson has a composition to evaluate against; null on every
    /// black-oil / plain-table case.
    std::shared_ptr<LiveFlashState> composition;

    /// 5.189: the Joule-Thomson coefficient from the fluid's own EOS at the
    /// local (T, p), in K/Pa. `srk` selects Soave-Redlich-Kwong; otherwise
    /// Peng-Robinson. Returns 0 when no composition was declared, which the
    /// parser refuses up front so a case cannot silently price zero.
    double eosJouleThomson(double T, double p, double cp, bool srk) const;

    /// Null (the default) => original piecewise-linear table path, every
    /// existing number bit-for-bit unchanged. Non-null => each call flashes
    /// the held composition at the local pressure.
    std::shared_ptr<LiveFlashState> liveFlash;
};

/// 5.158 (P18): one named segment of a DECLARED terrain trajectory. The
/// trajectory is the pipe's inclination laid out along its arc: segment k
/// covers the arc interval [sum of the earlier lengths, + its own length)
/// at angleDeg from horizontal (positive uphill). It is a REPORTING
/// geometry only -- the Beggs-Brill march keeps pricing the whole pipe at
/// the single declared `elevation` slope -- so the trajectory can describe
/// the dips the lumped slope hides without moving a hydraulic number.
struct TerrainTrajectorySegment
{
    double length = 0;   // m of arc along the pipe (> 0)
    double angleDeg = 0; // from horizontal, positive = uphill, in [-90, 90]
};

struct NetworkPipe
{
    std::string name;
    std::string fromName;
    std::string toName;
    label from = -1; // resolved by OneDNetwork::addPipe
    label to = -1;
    double length = 0;     // m, along the pipe
    double diameter = 0;   // m; the OUTER bore when annulusInner > 0, the
                           // WIDTH when ductHeight > 0 (5.80)
    double annulusInner = 0; // m; > 0 makes the bore a concentric annulus of
                             // this inner diameter (5.73); 0 = circular pipe,
                             // and every number below is bit-identical to the
                             // circular branch
    // Eccentric annulus (5.87, C3 increment 2): when annulusInner > 0 this is
    // the axis offset (m) of the inner bore from the outer one, 0 <= e <
    // (diameter - annulusInner)/2. The bore is still the gap between two
    // circles, so the area pi(Do^2-Di^2)/4, the wetted perimeter pi(Do+Di) and
    // Dh = Do - Di are all eccentricity-independent -- an offset closes the gap
    // on one side exactly as far as it opens it on the other, and both bores
    // stay fully wetted. What moves is the laminar profile: the wide side
    // carries disproportionately more of the flow, so the conductance RISES
    // above Fredrickson-Bird's and C = f.Re_Dh falls (94.713 concentric -> 71.2
    // at e = 0.5 gap, 43.5 at 0.98). e absent or 0 resolves nothing and every
    // concentric number is bit-identical. Newtonian single-phase steady only
    // this increment: the power-law annulus (5.74) is refused because its
    // zero-shear-radius integrals are concentric by construction, and so are
    // yield, a coupled edge and a duct clash.
    double eccentricity = 0; // m; axis offset of the inner bore, 0 = concentric
    // Rectangular duct (5.80): when ductHeight > 0 the bore is a rectangle of
    // width `diameter` and height `ductHeight` (both m). The laminar anchor is
    // the exact Fourier-series solution of the cross-section Poisson equation
    // grad^2 u = -G/mu: C(aspect) = f.Re_Dh is 96 for parallel plates
    // (aspect -> 0) and 56.91 for a square. 0 = not a duct, and every circular
    // / annulus number is bit-identical when this is absent. Newtonian
    // single-phase steady only this increment (power-law / yield / coupled /
    // two-phase / PVT / transient / two-fluid refused loudly).
    double ductHeight = 0;   // m; > 0 makes the bore a rectangular duct
    double roughness = 0;  // m, absolute roughness
    double elevation = 0;  // m, z(to) - z(from); adds rho*g*dz to the drop

    // inline equipment (single-phase pipes; the black-oil march ignores it)
    double valveK = 0;         // valve loss coefficient: dp = K rho v^2/2,
                               // always opposing the flow
    double pumpHead = 0;       // m of head at pumpRate (> 0 switches the
                               // pump on); constant head when pumpRate = 0
    double pumpShutoffHead = 0;// m at zero rate (defaults to pumpHead)
    double pumpRate = 0;       // m3/s where pumpHead is quoted
    // Multi-stage ESP + affinity laws (5.81, D1 increment 1): pumpStages N is
    // a count of identical stages in series (head adds at the same flow);
    // pumpSpeed s = N/N_rated rescales the curve by the Euler similarity laws
    // Q ~ s, H ~ s^2, i.e. H_s(q) = s^2 H_1(q/s) = s^2 h0 - (h0-h1)(q/q1)^2
    // (the quadratic term is speed-invariant, only shutoff scales). Combined
    // H(q) = stages * (s^2 h0 - (h0-h1)(q/q1)^2). stages 1 / speed 1 (or
    // absent) is the pump law above bit for bit. Single-phase steady only;
    // gas interference is a later increment.
    double pumpStages = 1;     // identical stages in series (>= 1)
    double pumpSpeed = 1;      // relative speed N/N_rated (> 0)
    // ESP gas interference (5.86, D1 increment 2): the catalog head is
    // degraded by the published Turpin (1986) two-phase correlation, read
    // through the shared pumpHeadOf so it multiplies the whole curve above
    // (stages, speed, constant head included). pumpGasFraction g is the free
    // gas volume fraction at pump intake in [0,1) (=> gas/liquid ratio
    // q_G/q_L = g/(1-g)); pumpIntakePressure is the ABSOLUTE intake pressure
    // in Pa (converted to psia for the correlation, p_i). With g > 0 the
    // returned head is H_tp = H_sp * exp(-a1 * (q_G/q_L)),
    // a1 = (346430/p_i^2)(q_G/q_L) - 410/p_i, p_i in psia. g absent or 0
    // short-circuits the multiplier, so the zero-declaration path stays IEEE
    // bit-identical to 5.81. Turpin is a curve-fit to Lea & Bearden air/water
    // and diesel/CO2 ESP data (25-500 psia, intake gas fraction to ~50%), a
    // pump-data match, not a universal predictor.
    double pumpGasFraction = 0;      // free gas volume fraction at intake [0,1)
    double pumpIntakePressure = 0;   // Pa, absolute pump intake pressure (p_i)

    // Beggs-Brill flow-pattern OVERRIDE (5.88, D4 increment 2): 5.85 made the
    // pattern map a queryable single source of truth but left the user unable
    // to say anything back to it. This is that channel: 0 segregated,
    // 2 intermittent, 3 distributed forces WHICH published holdup law prices
    // every segment of this two-phase pipe, while the map keeps being evaluated
    // and keeps reporting its own verdict (regimeFromMap) beside the one used.
    // The override is legitimate because Beggs-Brill is a correlation and a
    // field operator routinely knows the pattern; it is contained because it
    // touches only the two correlations that read a regime -- the horizontal
    // holdup and the inclination coefficient set -- leaving the published
    // friction multiplier, the acceleration term and the march itself alone.
    // 1 is refused: transition is not a holdup law, it is the map's own
    // interpolation between two of them, and forcing it outside the band that
    // defines its blend weight has no published meaning. -1 (absent) asks the
    // map, which is exactly the 5.85 path, bit for bit.
    int flowPattern = -1;

    // 5.143 (G4): WHICH published holdup law prices a two-phase pipe. 0 =
    // beggsBrill (the shipped law, bit for bit), 1 = hagedornBrown (Whitson
    // et al. 2017 restatement), 2 = dunsRos (API Bulletin 34 region
    // skeleton), 3 = olgasSteady (5.155, the OLGAS-steady public mechanistic
    // skeleton - see the provenance caveats in beggsBrill.cpp). The law
    // replaces ONLY the holdup block of the gradient; the flow-pattern map
    // stays evaluated and reported, and flowPattern (which selects one of
    // Beggs-Brill's own regime laws) is refused alongside a non-zero value.
    int correlation = 0;

    // 2026-11-25: with correlation hagedornBrown, opt into the Griffith
    // bubble / mist holdup branches of the classic vertical method (the
    // Whitson restatement alone carries only the middle chart). Absent or
    // false keeps every existing Hagedorn-Brown number bit-for-bit.
    bool hbGriffith = false;

    // 5.144 (G6): the PASSIVE steady temperature march. None of these
    // numbers feeds the pressure solution - they price the outlet
    // temperature along the march the hydraulics already ran - so a
    // declaration here cannot move a hydraulic answer. jouleThomson is
    // K/Pa applied to the march's own dp (signed; positive = cooling on
    // expansion, the usual gas behaviour); wallHeatTransfer is the overall
    // coefficient per unit length (W/(m K), film + wall + cement folded
    // into one number, as in the transient thermal block), or one value
    // per march segment through wallHeatTransferSegments (exactly 24
    // entries); the wall is constant (wallTemperature) or geothermal from
    // the inlet datum (surfaceTemperature + geothermalGradient * depth
    // below the inlet); heatCapacity is the stream's c_p (J/(kg K)),
    // needed by the wall term through m c_p dT/ds = U (T_wall - T). All
    // zero is the undeclared isothermal march, bit for bit.
    double jouleThomson = 0;
    /// 5.189: `jouleThomson eos` / `jouleThomson eosSrk` - the coefficient is
    /// the fluid's own EOS value at each segment's state instead of the
    /// constant above. Needs a compositional fluid; the parser refuses it
    /// anywhere else.
    bool jouleThomsonFromEos = false;
    bool jouleThomsonEosIsSrk = false;
    double wallHeatTransfer = 0;
    std::vector<double> wallHeatTransferSegments;
    int wallModel = 0;              // 0 constant, 1 geothermal
    double wallTemperature = 0;     // K
    double surfaceTemperature = 0;  // K (geothermal at the inlet datum)
    double geothermalGradient = 0;  // K/m
    double heatCapacity = 0;        // J/(kg K)

    /// 5.194 (G6d): TEMPERATURE FEEDBACK -- the march's own T reprices the
    /// fluid. Until now the temperature march was PASSIVE: it priced T along
    /// the pressure profile and nothing read the result back, so rho and mu
    /// were always the feed's single-temperature snapshot. With this switch
    /// (thermal{ temperatureFeedback true; }) every segment evaluates its
    /// rho/mu/sigma at (p, T) with the T the march has reached, and the mass
    /// rate is booked at the entry state (pUp, T_in). That is a COUPLED
    /// march: T now moves dp, and dp moves T, which is why this needs a
    /// statement in the case rather than being always on. Only the
    /// compositional live-flash path can honour it -- a tabulated rho(p)/mu(p)
    /// table is one temperature's snapshot, and the black-oil layer is
    /// isothermal -- so the parser refuses it there instead of pretending.
    bool pvtTemperatureFeedback = false;

    // 5.191 (G7): a positive-displacement machine (PCP / screw pump) on this
    // pipe. Its law is the machine's own, not the centrifugal curve above:
    //   q = n V - C dp    (displacement minus slip),  n rev/s, V m3/rev,
    //   C m3/(s Pa) the slip conductance,
    // so the head it adds is H = (nV - q) / (C rho g) and the row's slope
    // dH/dq = -1/(C rho g) is analytic. A rate above nV makes it a brake
    // (H < 0) rather than a pump; that is the law's own statement, not a
    // case error. Refused alongside a centrifugal `pumpHead` (one row, one
    // machine curve), a compressor, and on a coupled edge.
    bool pcpOn = false;
    double pcpDisplacement = 0; // m3/rev
    double pcpSpeed = 0;        // rev/s
    double pcpSlip = 0;         // m3/(s Pa)

    // 5.192 (G8): the jet pump as a row -- the slot that carries TWO streams.
    // A constant-area jet pump hands part of the motive stream's momentum to
    // an entrained stream, so a row that carries only one rate has nowhere to
    // put the second. The reduction this row takes is the EDUCATOR reading:
    // both streams are drawn from the row's own upstream node, so the row's
    // rate is the combined q_d = q_n + q_s and the machine splits it
    // internally into the nozzle stream q_n and the entrained q_s = q - q_n.
    // Mass then closes at both nodes with the one rate the network solves,
    // which is what a lumped nodal balance can carry; a surface power-fluid
    // line is a separate mass source, and the steady network has no slot for
    // one (that is a boundary, not an approximation).
    //
    // The rise is the element's own throat balance read as a head
    //   H(q) = rho v_n^2 / m^2 * [ (m-1-r)^2/(m-1) - F (1+r)^2/2 ],
    //   r = (q - q_n)/q_n,  m = A_t/A_n,  v_n = q_n/A_n,  F = f L_t/D_t,
    // which is the same algebra jetPumpSolve prices (mode B: the declared
    // suction rate) and increasing in q, so the rate bracket keeps its
    // monotone branch. Outside the element's envelope the head is clamped to
    // its own two states: r <= 0 is NO-SUCTION (the stall head, the machine
    // passing its motive stream alone) and H <= 0 is FREE-JET (zero rise at
    // the entrainment ceiling).
    bool jetOn = false;
    double jetMotiveRate = 0;    // m3/s through the nozzle q_n
    double jetNozzleArea = 0;    // m2, A_n
    double jetThroatArea = 0;    // m2, A_t > A_n
    double jetThroatLength = 0;  // m, L_t
    double jetFrictionFactor = 0;// Darcy f on the throat wall

    /// 5.188 (G6b): the layered wall, when `insulation{}` priced U instead of
    /// a declared wallHeatTransfer. The stack is radial from innerRadius
    /// outward, thickness/conductivity one entry per layer, plus an optional
    /// outer film coefficient (W/(m2 K), 0 = no film). The per-unit-length
    /// overall coefficient is the series sum
    ///   U = 1 / ( sum_i ln(r_{i+1}/r_i)/(2 pi k_i) + [1/(2 pi r_out h_o)] ),
    /// which is what the thermal march consumes; outerRadius is kept only so
    /// the report can show the stack it priced.
    bool insulatedLayered = false;
    std::vector<double> layerThickness;   // m, outer radius step per layer
    std::vector<double> layerConductivity; // W/(m K)
    double insulationInnerRadius = 0;     // m
    double insulationOuterRadius = 0;     // m
    double insulationOuterFilm = 0;       // W/(m2 K), 0 = no film

    /// 5.144 (G6): does this pipe carry any PRICED thermal physics? A
    /// Joule-Thomson coefficient or a positive wall coefficient anywhere.
    /// Declared-but-inert keys are refused by the parser; this is what the
    /// marches and addPipe scope checks act on.
    bool thermalOn() const
    {
        if (jouleThomson != 0.0 || jouleThomsonFromEos
            || wallHeatTransfer > 0.0) return true;
        for (double u : wallHeatTransferSegments)
            if (u > 0.0) return true;
        return false;
    }

    /// 5.144 (G6): does the pipe declare wall exchange (as opposed to a
    /// pure Joule-Thomson readout)? The relaxation term needs c_p and a
    /// wall temperature.
    bool wallExchangeOn() const
    {
        if (wallHeatTransfer > 0.0 || !wallHeatTransferSegments.empty())
            return true;
        return false;
    }

    // Gas-lift valve (5.89, D1 increment 3): this pipe IS an injection-pressure-
    // operated valve between the injection (annulus) node and the production
    // (tubing) node. `liftDomePressure` > 0 switches it on; `diameter` is then
    // the PORT bore, not a pipe bore (the gas has no other path), and the pipe
    // carries no length, roughness, elevation, rheology or other equipment.
    // Two published closures price it, both in GasLiftValve.h: the bellows
    // force balance decides WHETHER it flows (Guo et al. 2007 Eq 13.42/13.43),
    // the isentropic nozzle decides HOW MUCH. It is a tabulated-PVT-only
    // element because the nozzle needs the gas density at the upstream pressure
    // the valve actually sees; on a constant-property or black-oil network the
    // declaration is refused rather than quietly priced with the wrong rho.
    // Open/closed is discontinuous, so the closing band is not guessed: an
    // injection pressure between Pvc and Pvo is reported as indeterminate and
    // priced closed, because the steady solve has no memory of where the valve
    // came from.
    double liftDomePressure = 0;   // Pa absolute, nitrogen charge at valve T (Pd)
    double liftAreaRatio = 0;      // R = Ap/Ab, in (0,1)
    double liftSpringPressure = 0; // Pa, spring term St (0 = no spring)
    double liftDischargeCoeff = 0; // Cd, in (0,1]
    double liftKappa = 0;          // cp/cv of the lift gas (> 1)

    // Rod pump (5.90, D1 increment 4): this pipe IS a subsurface sucker-rod
    // pump between the intake node it draws from and the tubing node it
    // discharges into, in the declared direction. `rodPlungerDiameter` > 0
    // switches it on; it is the only pump here that is not a head curve,
    // because a plunger does not add head at a rate the network chooses -- it
    // moves a volume the machine dictates, A_p S_eff N, and the pressure rise
    // it works against is whatever the columns above and below it add up to.
    // The element is therefore read backwards from every other pipe: its rate
    // is the input and its `dpTotal` (negative while it lifts) is the answer.
    // What comes out is the dynamometer card -- the buoyed rod weight, the
    // fluid load the plunger actually holds, the stretch it costs, the stroke
    // left over, the work per stroke and the power -- all of it statics plus
    // Hooke's law, priced at the network's own two port pressures. See
    // RodPump.h for the four statements. It runs on a constant-property fill
    // only this increment: the barrel's mass balance is exact when the fluid's
    // density does not ride with its pressure, and the node model books volume
    // per supplying row, so a compressible fill would double-count the
    // expansion on the way to the surface. The wave equation that turns the
    // polished-rod card into the downhole one (Gibbs 1953), partial fillage,
    // slippage past the plunger and gas in the barrel are later increments.
    double rodPlungerDiameter = 0; // m, barrel bore A_p (the switch)
    double rodStrokeLength = 0;    // m, polished-rod stroke S
    double rodStrokesPerMinute = 0;// spm
    double rodDiameter = 0;        // m, rod string gauge A_r
    double rodLength = 0;          // m, pump setting depth L
    double rodYoungModulus = 0;    // Pa, rod string E
    double rodDensity = 0;         // kg/m3, rod string density

    // Compressor (5.156, P10): this pipe's line carries a steady pressure
    // booster at its `from` end. `compressorRatio` (the key `pressureRatio`)
    // > 0 switches it on; the ratio falls with the square of the rate through
    // the same quadratic family as pumpHeadOf (see compressorRatioOf), and
    // the effective statement to the pipe law is that the from-end pressure
    // the march sees is pFrom * PR(|q|). Direction is decided at the
    // zero-rate (shutoff) boost, so a machine can push a line against a
    // node-to-node gradient the bare pipe would run backwards under. It runs
    // on the tabulated-PVT fill only, like the lift valve: the boost is a
    // pressure-level statement and the readouts (power, discharge
    // temperature) need the stream's declared inlet temperature. The
    // thermodynamics are PASSIVE (§5.144's convention): the machine's kappa,
    // R_spec, Z and efficiency price its power and T2 readouts from the
    // hydraulic state, never into it. A head-mode catalog entry
    // (polytropicHead) is the point element's grammar (Compressor.h), not a
    // pipe key: on a line the curve is what the network needs.
    double compressorRatio = 0;         // PR at ratedRate (> 1); the switch
    double compressorShutoffRatio = 0;  // PR at zero rate (defaults to ratio)
    double compressorRatedRate = 0;     // m3/s where pressureRatio is quoted
    double compressorEfficiency = 1;    // polytropic eta, in (0,1]
    double compressorKappa = 0;         // cp/cv, > 1
    double compressorGasConstant = 0;   // R_spec, J/(kg K)
    double compressorZ = 1;             // declared constant compressibility

    // 5.158 (P18): the slug-START forecast. A declared trajectory lays the
    // pipe's dips out along its arc; at every dip (a strict downhill -> uphill
    // angle change) the liquid the march's own holdup profile says is flowing
    // down the upstream run is pooled on the uphill run downstream of the
    // dip, and the pooled level against the declared fill threshold, together
    // with the upstream segment's flow pattern, is the verdict. Pure report:
    // the march never reads any of this back. Absent, every field below is
    // inert and no line is printed.
    bool slugDeclared = false;
    std::vector<TerrainTrajectorySegment> slugTrajectory;
    double slugFillThreshold = 1.0; // pooled level / uphill length, >= this
                                    // at a segregated upstream -> slug start

    // 5.159 (P22): the erosion / hydrate / corrosion risk report, priced
    // from the march's own per-segment profile. Erosion: the API RP 14E
    // velocity limit Ve = C / sqrt(rho_m) in field units (C in the published
    // 100 ft/s-(lb/ft3)^0.5 convention, rho_m converted kg/m3 -> lb/ft3),
    // read at the segment of worst vm/Ve. Hydrate: the declared methane-
    // hydrate Clausius-Clapeyron line (Tref/pRef/dH, §5.122 family) priced
    // at each segment's march pressure, subcooling Teq - T read at the worst
    // segment (T is the march temperature when the pipe carries thermal
    // physics, else the fluid's isothermal temperature). Corrosion: CO2/H2S
    // partial pressures at the highest march pressure against the declared
    // bands. Every one is a report, not a charge. Absent, inert.
    bool riskDeclared = false;
    double erosionC = 100.0;            // API RP 14E C factor (> 0)
    bool hydrateDeclared = false;       // all three hydrate keys present
    double hydrateTref = 0;             // K, dissociation reference temperature
    double hydratePRef = 0;             // Pa, dissociation pressure at Tref
    double hydrateDH = 0;               // J/mol, molar dissociation enthalpy
    double hydrateLowK = 3;             // K subcooling bands: < low NONE..LOW,
    double hydrateMediumK = 8;          //   < medium LOW..MEDIUM, else HIGH
    double co2Fraction = 0;             // gas-phase CO2 mole fraction [0,1]
    double h2sFraction = 0;             // gas-phase H2S mole fraction [0,1]
    double co2LowPa = 1e4;              // CO2 partial-pressure bands (Pa)
    double co2HighPa = 1e5;

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

    // Fredrickson-Bird (1958) power-law annulus, resolved once by addPipe
    // when annulusInner > 0 and a flowIndex is declared: the zero-shear
    // radius (m) and the dimensionless profile integral with
    //     Q = (dp / (2 K L))^(1/n) . 2 pi . annulusDimQ.
    // Both depend only on (kappa = dInner/dOuter, n), never on the rate, so
    // they are solved once instead of inside every rate solve.
    double annulusLambda = 0;
    double annulusDimQ = 0;

    // Rectangular-duct laminar bracket (5.80), resolved once by addPipe when
    // ductHeight > 0: bracket(aspect) = 1 - (192 aspect/pi^5) sum_{n odd}
    // tanh(n pi/(2 aspect))/n^5, aspect = short/long. The rate is
    // Q = (dp/L) . long . short^3 . bracket/(12 mu) and C = 96/((1+aspect)^2
    // bracket); both depend only on the aspect, never on the rate, so the
    // series is summed once instead of inside every rate solve.
    double rectBracket = 0;

    // Eccentric-annulus cross-section integral (5.87), resolved once by addPipe
    // when eccentricity > 0: K_geom = integral_A w dA of the laminar Poisson
    // solution grad^2 w = -1 on the bore, m^4. The whole flow law needs only
    // this one number -- conductance = K_geom/(mu L) and C = 2 Dh^2 A/K_geom --
    // and it depends on (Ro, Ri, e) alone, never on the rate, the fluid or the
    // pressure, so the bipolar-coordinate solve runs at load time instead of
    // inside every rate solve. 0 = concentric, which stays on the closed form.
    double annulusEccK = 0;

    /// Coupled element (the 1D side of a 1D<->3D interface): when > 0 this
    /// pipe carries NO 1D resistance of its own – its whole loss is
    /// `coupledResistance * |Q|`, supplied by the partner model, and the
    /// friction/valve terms are skipped so nothing is counted twice. The
    /// coupling layer re-solves the network with an updated value until the
    /// drop the network sees equals the drop the detailed model reports.
    /// Length and diameter stay meaningful: they are the detailed region's own
    /// geometry, used for the reported velocity and Reynolds number.
    double coupledResistance = 0; // Pa/(m3/s), i.e. Pa.s/m3

    /// The fluid of a coupled edge's PARTNER region. A case declares it with
    /// the same `flowIndex` / `consistency` / `yieldStress` keys a wall-law
    /// pipe uses, but the 1D row never evaluates it: this pipe's loss is the
    /// partner's own statement (§5.43), so the parser moves the declaration
    /// here and leaves the three wall-law fields above at zero. The coupling
    /// apps read it back and hand it to the 3D detailed block, which turns it
    /// into a per-face apparent viscosity. Absent on an uncoupled pipe.
    Rheology partnerRheology;

    // 5.193 (G9): a WELL-TRAJECTORY SURVEY. Until now a pipe could only state
    // one slope for its whole arc (`elevation`/length, or the terrain
    // `trajectory{}` the slug forecast reads), so a build section and a
    // horizontal lateral were both priced as if the fluid fell down one
    // straight line. A survey states inclination against MEASURED DEPTH, the
    // way a directional well is actually delivered:
    //     survey { md (0 500 1000 ...); inclination (0 42 88 ...); }
    // md is metres along the hole, monotone, from 0 to the pipe's length;
    // inclination is degrees from horizontal, positive uphill. Between
    // stations the inclination is taken LINEAR IN MD (a telescoped survey:
    // azimuth is not carried, because hydraulics never read it -- only the
    // vertical projection sin(theta) does), and the march then prices its own
    // gravity segment by segment instead of one sin for the whole line.
    // `elevation` stops being a declaration and becomes a consequence: absent
    // it, the survey's own integral supplies it, and a declared one that
    // disagrees is refused rather than quietly overridden.
    bool surveyOn = false;
    std::vector<double> surveyMd;       // m, monotone, 0 .. length
    std::vector<double> surveyInclDeg;  // deg from horizontal, + uphill
};

/// 5.193 (G9): the vertical drop a survey's own integral states,
///     dz = integral_0^L sin(theta(s)) ds,
/// exact per survey interval (theta linear in s, so sin(theta(s)) integrates
/// to (cos(th1) - cos(th2))/(th2 - th1) and the theta = 0 interval to
/// sin(th) ds). Absent a survey this is `p.elevation`, bit for bit.
double surveyElevationOf(const NetworkPipe& p);

/// 5.193 (G9): sin(theta) at the MIDPOINT of each of nSeg equal march
/// segments, i.e. the gravity term the march should charge for segment k.
/// Absent a survey this is nSeg copies of elevation/length -- the constant
/// path the marches had before, bit for bit.
std::vector<double> surveySegSinOf(const NetworkPipe& p, int nSeg);

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
    // D10/D7: an optional per-station surface-tension provider published by
    // the loader when a live-flash case declares sigmaFromP. The PVT march
    // consumes its own PvtFluid::surfaceTension directly; this hook is the
    // same WK callback exposed to the (separate) Beggs-Brill two-phase march
    // via solvePipeFlowTwoPhase. Empty on every existing case.
    std::function<double(double)> sigmaOfP;

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
/// pipe's OWN laminar wall law (5.46/5.47) asks for at the locus rate, i.e. how
/// far the injection end has to move before the laminar statement stops being
/// the one that applies. The locus numbers themselves are not a charge; past
/// Re_B,t a smooth Bingham pipe is priced by the Wilson-Thomas turbulent
/// closure (5.83), and dpAtTransition stays laminar because it is read AT the
/// locus, not past it.
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

/// The dynamometer card of one rod pump, as the steady network reads it (see
/// RodPump.h for the four statements behind every number). `state` is -1 on
/// every row that is not a pump; 0 full card, 1 no fluid load across the
/// plunger, 2 the load stretches the string past the whole stroke. Loads are
/// newtons of rod force at the pump's depth, `maxLoad`/`minLoad` at the
/// polished rod and `downholeMaxLoad`/`downholeMinLoad` at the plunger, so the
/// two pairs differ by exactly the buoyed string weight. Travels are metres of
/// polished-rod stroke: `stretch` is what the fluid load costs and
/// `effectiveStroke` is what the plunger actually travels, the parallelogram's
/// width being `stroke - stretch` and its height exactly `fluidLoad`. The rate
/// is barrel volume at the intake conditions and the power is the card's area
/// per second, so `workPerCycle == fluidLoad * effectiveStroke` and
/// `hydraulicPower == dpPump * displacement` are both restatements of the same
/// parallelogram, not a second law.
struct RodPumpReadout
{
    int state = -1;
    double plungerArea = 0;      // A_p, m2
    double rodArea = 0;          // A_r, m2
    double rodWeightInAir = 0;   // N
    double buoyancy = 0;         // N
    double rodWeightInFluid = 0; // N
    double fluidLoad = 0;        // Fo, N
    double fluidLoadMax = 0;     // rho g L A_p, N
    double upthrust = 0;         // fluidLoadMax - fluidLoad, N (signed)
    double maxLoad = 0;          // W_max at the polished rod, N
    double minLoad = 0;          // W_min, N
    double downholeMaxLoad = 0;  // N
    double downholeMinLoad = 0;  // N
    double counterbalance = 0;   // (W_max + W_min)/2, N
    double stretch = 0;          // delta_f, m
    double stroke = 0;           // S, m
    double effectiveStroke = 0;  // S_eff, m
    double efficiency = 0;       // S_eff / S
    double displacement = 0;     // A_p S_eff N, m3/s
    double massRate = 0;         // rho * displacement, kg/s
    double workPerCycle = 0;     // Fo S_eff, J
    double hydraulicPower = 0;   // work per cycle x strokes per second, W
};

/// 5.158 (P18): one terrain dip's slug-start readout (see slugReadoutOf).
struct SlugDipReadout
{
    double s = 0;           // m, arc position of the dip
    double volume = 0;      // m3, liquid pooled on the uphill run
    double criterion = 0;   // volume / (A * uphill length), a level fraction
    int regime = -1;        // march regime at the upstream run's midpoint
    bool regimeOk = false;  // the upstream run is segregated
    bool slugStart = false; // regimeOk and criterion >= the fill threshold
};

struct SlugReadout
{
    std::string pipeName;
    double fillThreshold = 1.0;
    std::vector<SlugDipReadout> dips; // in arc order
};

/// 5.159 (P22): the erosion / hydrate / corrosion readout of one pipe (see
/// riskReadoutOf). Each `*On` says that section has a line to print; the
/// seg/p/T numbers are read at the WORST march segment of that section.
struct RiskReadout
{
    std::string pipeName;
    bool erosionOn = false;
    int erosionSeg = -1;
    double erosionVm = 0;     // m/s at the worst segment
    double erosionVelim = 0;  // API RP 14E Ve, m/s
    double erosionRatio = 0;  // vm / Ve
    int erosionOverCount = 0; // segments with ratio > 1
    bool hydrateOn = false;
    int hydrateSeg = -1;
    double hydrateP = 0;      // Pa at the worst segment
    double hydrateT = 0;      // K at the worst segment
    double hydrateTeq = 0;    // K, dissociation temperature at hydrateP
    double hydrateSubcool = 0; // Teq - T, K (negative = no risk)
    int hydrateGrade = 0;     // 0 NONE, 1 LOW, 2 MEDIUM, 3 HIGH
    bool corrosionOn = false;
    double co2Partial = 0;    // Pa at the highest march pressure
    double h2sPartial = 0;    // Pa
    int co2Band = 0;          // 0 LOW, 1 MEDIUM, 2 HIGH
    bool sour = false;        // H2S partial pressure >= 0.5 psia (3447 Pa)
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
    RodPumpReadout card;          // the dynamometer card, for a rod pump only
    double holdup = 0;      // two-phase: outlet liquid volume fraction
    int regime = -1;        // two-phase: Beggs-Brill regime of the last
                            // segment (0 segregated, 1 transition,
                            // 2 intermittent, 3 distributed; -1 single-phase)
    // 5.85 (D4 increment 1): the flow-pattern dispatch surface read out at
    // the outlet segment -- its own regime, no-slip liquid content, Froude
    // number and the four L boundaries that selected its holdup. `regime`
    // above stays the along-pipe envelope (regimeMax); these describe the
    // last segment alone, so the map is readable where the march ended.
    int regimeOutlet = -1;
    // 5.88 (D4 increment 2): what the map said at the outlet segment. Equal to
    // regimeOutlet unless pipe.flowPattern overrode it, in which case the pair
    // is the override's footprint on the report: `regimeOutlet` is what was
    // priced, this is what the correlation would have picked on its own.
    int regimeFromMapOutlet = -1;
    double lambdaOutlet = 0;
    double nfrOutlet = 0;
    double patternL1 = 0, patternL2 = 0, patternL3 = 0, patternL4 = 0;
    double massRate = 0;    // two-phase: kg/s, invariant along the pipe
    // 5.144 (G6): outlet temperature (K) of a marched pipe (black-oil or
    // pvt-table) whose thermal physics was declared; 0 on every undeclared
    // or constant-property path, where the march is isothermal.
    double temperature = 0;
    // 5.89 (D1 increment 3): the gas-lift valve's own statement, on a valve
    // pipe only (disposition -1 everywhere else). `flowRate` is the volume rate
    // at the UPSTREAM (injection-side) density the nozzle was priced on, and
    // `massRate` is Cd times the isentropic flux, so massRate/flowRate is that
    // density by construction. `valveDisposition` separates the two ways a
    // valve can pass nothing: below its closing pressure it is shut, between
    // the closing and opening pressures it is in the spread and the steady
    // solve has no basis to pick, so it is priced shut and says so.
    int valveDisposition = -1; // 0 closed, 1 open, 2 in the spread band
    double valveOpeningPressure = 0; // Pvo at the solved tubing pressure, Pa abs
    double valveClosingPressure = 0; // Pvc, Pa abs
    double valveSpread = 0;    // Pvo - Pvc, Pa; negative once the tubing alone
                               // holds the port off its seat
    double valveProductionEffect = 0; // R/(1-R) Pt: how far the tubing has
                                      // pulled Pvo down from its test-rack value
    double valveTestRackOpening = 0; // the same valve's Pvo at 1 atm tubing side
    bool valveChoked = false;        // sonic at the port: downstream is unheard
    bool choked = false;    // two-phase: the march hit the pressure floor at
                            // the imposed end pressures
    // D10/D7: station-by-station live-flash state, populated only when the
    // pipe ran with a true-flash PvtFluid. sigmaOutlet is the WK surface
    // tension at the outlet station (0 on every other path).
    std::vector<PvtFlashStation> pvtFlashStations;
    double sigmaOutlet = 0;
    // 5.156 (P10): the compressor readout on a pvt-table pipe carrying a
    // pressureRatio. Ratio 1 and zeros on every other path, so no existing
    // line gains a token. `compressorBoost` is the absolute pressure the
    // machine adds at its `from` end, pFrom (PR - 1); `compressorPower` is
    // the shaft power m H / eta at the curve's ratio and the march's own
    // mass rate; `compressorTemperatureOut` is the passive discharge
    // temperature T1 PR^((k-1)/(k eta)).
    double compressorRatio = 1;
    double compressorBoost = 0;
    double compressorPower = 0;
    double compressorTemperatureOut = 0;
    // 5.158 / 5.159: the final march's per-segment profile, copied from the
    // black-oil march's readout vectors (empty on every single-phase / pvt /
    // constant-property path). The slug and risk reports read these; nothing
    // that prices a hydraulic number does.
    std::vector<double> marchSegPressure;
    std::vector<double> marchSegHoldup;
    std::vector<double> marchSegLambda;
    std::vector<double> marchSegVm;
    std::vector<double> marchSegRhoM;
    std::vector<int> marchSegRegime;
    std::vector<double> marchSegTemperature;
    // 5.188 (G6b): the per-segment ABSOLUTE position of the march, filled only
    // on a thermal pipe: x in metres from the pipe inlet, one entry per
    // point (the N segment inlets plus the outlet), so a report can print
    // "p and T at this distance" without re-deriving ds.
    std::vector<double> marchSegX;
    double marchThermalU = 0;   // W/(m K) the march actually priced
    double marchOutletPressure = 0; // Pa, the march's own exit pressure
    // 5.189: the Joule-Thomson coefficient the march used at its INLET state
    // and where it came from (0 = the declared constant, 1 = pengRobinson
    // state difference, 2 = SRK's analytic dZ/dT). 0.0 / 0 on every pipe
    // without thermal physics.
    double dpPcp = 0;           // Pa, the PCP's own contribution (5.191)
    double dpJet = 0;           // Pa, the jet pump's own rise (5.192)
    double jetFlowRatio = 0;    // r = q_s/q_n the row settled on (5.192)
    int jetState = -1;          // -1 none, 0 PUMPING, 1 NO-SUCTION, 2 FREE-JET
    double marchJouleThomson = 0;
    int marchJtMode = 0;
    /// 5.194 (G6d): did the march price properties at its own T, and at what
    /// entry state did it book the mass rate? Both are zero / false on every
    /// pipe without the switch, so no existing readout moves.
    bool marchTemperatureFeedback = false;
    double marchEntryDensity = 0; // kg/m3 at (pUp, T_in)
    double marchMassRate = 0;     // kg/s, the march's conserved quantity
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
/// Metzner-Reed number, a yielded pipe through its own laminar wall stress, or
/// -- for a smooth Bingham plastic strictly past its locus -- the Wilson-Thomas
/// turbulent closure (5.83). Declared here because the transient momentum row
/// and the steady pipe loss are two callers of one law; a rheology read
/// differently in the two kernels would be a third set of assertions nobody
/// could reconcile.
double pipeFrictionPa(const NetworkPipe& pipe, const FluidProperties& fluid,
                      double q);

/// The transition locus of one yielded pipe at a rate magnitude |q| (see
/// YieldedTransition). The locus numbers are a report, not a charge; but a
/// smooth Bingham pipe strictly past Re_B,t is now priced by the Wilson-Thomas
/// turbulent closure (5.83) inside pipeFrictionPa, while dpAtTransition --
/// evaluated AT the locus -- stays on the laminar wall law.
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
                           double pTo, PipeFlowResult& out,
                           const std::function<double(double)>* sigmaOfP = nullptr);

/// 5.158 (P18): the slug-start forecast of one pipe. Reads the pipe's
/// declared trajectory and fill threshold, and the final march's per-segment
/// profile from `flow.marchSeg*`. Empty result on every undeclared pipe.
SlugReadout slugReadoutOf(const NetworkPipe& pipe, const PipeFlowResult& flow);

/// 5.159 (P22): the erosion / hydrate / corrosion readout of one pipe.
/// Reads the pipe's declared risk fields and the final march's per-segment
/// profile from `flow.marchSeg*`. Empty result on every undeclared pipe.
RiskReadout riskReadoutOf(const NetworkPipe& pipe,
                          const pvt::BlackOilFluid& fluid,
                          const PipeFlowResult& flow);

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
/// Live-flash variant (D10/D7): fluid { pvt { pressure(...);
///   compositional { composition {...}; temperature 320;
///                   liveFlash true; sigmaFromP true;
///                   parachors fromTable; } } } — rho/mu/sigma are re-flashed
/// at every march station. liveFlash / sigmaFromP absent or false keeps the
/// 5.110 precomputed-table path and every existing number bit-identical.
/// Gas-lift valve (5.89, D1 increment 3), on a fluid { pvt } network only:
///   pipes { VLV { from ANN; to T; length 0; diameter 0.0047625;
///                 liftDomePressure 7.0e6;     // Pa abs; the switch, 0 = no valve
///                 liftAreaRatio 0.10;         // R = Ap/Ab, in (0,1)
///                 liftSpringPressure 0;       // Pa, the St term
///                 liftDischargeCoeff 0.665;   // Cd, in (0,1]
///                 liftKappa 1.3; } }          // cp/cv of the lift gas
/// `diameter` is the port bore, `length` is required and must read 0, and
/// `elevation` must be absent or 0: a valve is a point element, and a wall law,
/// bore shape, pump, coupled edge or flow-pattern override co-declared with it
/// is refused rather than unread.
/// Rod pump (5.90, D1 increment 4), on a constant-property `fluid` network only:
///   pipes { PUMP { from ANN; to TUB; length 0;
///                  rodPlungerDiameter 0.04445;  // barrel bore; the switch
///                  rodStrokeLength 1.5;         // polished-rod stroke S
///                  rodStrokesPerMinute 8;       // spm
///                  rodDiameter 0.019050;        // rod string gauge
///                  rodLength 1500;              // pump setting depth L
///                  rodYoungModulus 2.0e11;      // string E
///                  rodDensity 7850; } }         // string rho
/// `from` is the intake the barrel fills from and `to` the tubing it delivers
/// into; the pump is a point element like the valve, so `diameter`, `length`
/// and `elevation` all read 0, and a wall law, bore shape, head curve, lift
/// valve, coupled edge or flow-pattern override co-declared with it is refused.
OneDNetwork loadNetworkFromDict(const Dictionary& caseDict);
SteadyControls loadSteadyControls(const Dictionary& caseDict);

} // namespace cf
