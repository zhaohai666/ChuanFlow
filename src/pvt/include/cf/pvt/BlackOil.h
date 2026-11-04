#pragma once

#include <string>
#include <vector>

namespace cf::pvt
{

/// Black-oil property layer: the standard petroleum-engineering correlation
/// set (Standing, Vazquez-Beggs, Beggs-Robinson, Lee, Dranchuk-Abou-Kassem),
/// implemented from the published defining equations. All primitive
/// correlations below use their native oilfield units (psia, degF, degR,
/// scf/STB, rb/STB, cp, lb/ft3) exactly as published; the SI adapter is the
/// BlackOilFluid/BlackOilState pair at the bottom, which is what the 1D
/// kernels consume.
///
/// Standard conditions for all SI surface rates: 101325 Pa, 288.15 K.
/// Gas specific gravity gamma_g = M_gas / M_air, M_air = 28.9647 g/mol.
/// API gravity is derived from the stock-tank oil density against water at
/// 60 degF (999.014 kg/m3): API = 141.5/gamma_o - 131.5.

// ---- oilfield-unit primitives (published form, oilfield units) ----

/// Standing (1947) bubble point, psia:
/// pb = 18.2 [(Rs/gamma_g)^0.83 10^a - 1.4], a = 0.00091(T-460) - 0.0125 API.
/// Clamped at 0 (dead oil -> pb = 0).
double standingBubblePointPsia(double gasGravity, double api, double TF,
                               double rsScfStb);

/// Vazquez-Beggs (1980) solution GOR, scf/STB:
/// Rs = C1 gamma_g p^C2 exp(C3 API / T), T in degR; C1/C2/C3 switch at
/// API 30 (0.0362/1.0937/25.724 and 0.0178/1.187/23.931).
double vazquezBeggsRsScfStb(double pPsia, double gasGravity, double api,
                            double TF);

/// Vazquez-Beggs (1980) oil formation volume factor, rb/STB. Below the
/// supplied saturation pressure the saturated form
/// Bo = 1 + C1 Rs + C2 (T-60) API/gamma_g + C3 Rs (T-60) API/gamma_g
/// applies; above it the undersaturated branch Bo = Bob exp[Co (pb - p)]
/// with the Vazquez-Beggs compressibility Co. `rsScfStb` is the solution GOR
/// at `pPsia`, `pbPsia` the saturation pressure.
double vazquezBeggsBoRbStb(double pPsia, double rsScfStb, double pbPsia,
                           double gasGravity, double api, double TF);

/// Beggs-Robinson (1975) dead oil viscosity, cp:
/// mu_od = 10^X - 1, X = 10^(3.0324 - 0.02023 API) T^-1.163 (T in degF).
double beggsRobinsonDeadOilCp(double api, double TF);

/// Beggs-Robinson (1975) live oil viscosity, cp: mu_o = A mu_od^B,
/// A = 10.715 (Rs + 100)^-0.515, B = 5.44 (Rs + 150)^-0.338.
double beggsRobinsonLiveOilCp(double deadOilCp, double rsScfStb);

/// Lee-Gonzalez-Eakin (1966) gas viscosity, cp:
/// mu = 1e-4 K exp(X rho^Y) with rho in g/cm3, T in degR,
/// K = (9.4 + 0.02 M) T^1.5 / (209 + 19 M + T),
/// X = 3.5 + 986/T + 0.01 M, Y = 2.4 - 0.2 X.
double leeGasViscosityCp(double molecularWeight, double TR, double rhoGPerCm3);

/// Standing pseudo-criticals of a gas of specific gravity gamma_g (degR /
/// psia): low-gravity form below 0.75, the heavier-gas form above.
void standingPseudoCritical(double gasGravity, double& ppcPsia,
                            double& tpcRankine);

/// Dranchuk-Abou-Kassem (1975) real-gas deviation factor: 11-constant fit of
/// the Standing-Katz chart, solved by damped Newton on the reduced density.
/// p in psia, T in degR, pseudo-criticals in the same units.
double dranchukAbouKassemZ(double pPsia, double TR, double ppcPsia,
                           double tpcRankine);

/// Jain (1976) explicit friction factor seed / Colebrook-White Darcy friction
/// factor solved by fixed-point Newton to machine precision.
/// Re < 2040 returns the laminar 64/Re.
double colebrookDarcyFriction(double reynolds, double relativeRoughness);

// ---- SI-facing adapter ----

/// Isothermal black-oil stream in SI. Per sm3 of stock-tank oil the stream
/// carries gorStd sm3 of gas at standard conditions and, when waterCut > 0,
/// water at the volume ratio WOR = waterCut/(1 - waterCut) to oil.
struct BlackOilFluid
{
    double gasGravity = 0.65;     // air = 1
    double oilDensityStd = 850.0; // kg/m3, stock-tank oil
    double gorStd = 0.0;          // sm3/sm3 producing GOR
    double waterCut = 0.0;        // fraction of std liquid volume
    double waterDensity = 1000.0; // kg/m3 (input constant)
    double waterViscosity = 1e-3; // Pa s (input constant)
    double temperature = 300.0;   // K, isothermal
    double surfaceTension = 0.03; // N/m gas-liquid (input constant)
};

/// Black-oil state at pressure p: phase volumes and properties per sm3 of
/// stock-tank oil. Above the stream's own bubble point all gas is dissolved
/// (rs = gorStd) and the undersaturated Bo branch applies; below it rs
/// follows Vazquez-Beggs and the excess gas flows free with FVF bg.
struct BlackOilState
{
    double pb = 0;       // Pa, stream bubble point
    double rs = 0;       // sm3/sm3 dissolved gas
    double bo = 1;       // rm3/sm3 oil volume factor
    double bg = 1;       // rm3/sm3 free-gas volume factor
    double freeGas = 0;  // sm3/sm3 free gas at std
    double rhoOil = 0;   // kg/m3 at p
    double rhoGas = 0;   // kg/m3 at p
    double rhoWater = 0; // kg/m3 at p
    double muOil = 0;    // Pa s at p
    double muGas = 0;    // Pa s at p
    double z = 1;        // free-gas deviation factor
};

BlackOilState evaluateBlackOil(const BlackOilFluid& fluid, double pPa);

/// Invariant mass carried by the stream per sm3 of stock-tank oil, computed
/// from standard-condition densities only: stock-tank oil + produced gas +
/// water. Phase volume changes along a pipe redistribute this mass between
/// phases but never change it, which makes it the marching conserved quantity.
double blackOilMassPerOilSm3(const BlackOilFluid& fluid);

/// Liquid-phase properties at p: volume-weighted density of reservoir oil and
/// water (water FVF = 1) and Arrhenius log-mix viscosity at the in-situ
/// volume fractions. Water properties are input constants for now.
void blackOilLiquid(const BlackOilFluid& fluid, const BlackOilState& st,
                    double& rhoL, double& muL);

} // namespace cf::pvt
