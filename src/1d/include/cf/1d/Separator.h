#pragma once

#include "cf/pvt/BlackOil.h"

namespace cf
{

/// P10 (5.156): a single-stage black-oil separator as a 1D steady point
/// element -- the same physics the pvt cascade app (M6) prices, read at the
/// separator's own (p, T) and reported as the two outlet STREAMS of a feed
/// whose stock-tank oil rate is imposed. It is not a network row: a separator
/// takes one inlet and gives two outlets at (almost) one pressure, which the
/// two-node pipe model of §5.140 cannot carry, so it lives beside the jet
/// pump and the plunger gas lift as a closed-form element the solver reports
/// and never iterates on.
///
/// The allocation is component-wise and exact by construction. Per sm3 of
/// stock-tank oil the feed carries oilDensityStd kg of oil, gorStd sm3 of gas
/// at standard conditions and, at waterCut w, w/(1-w) sm3 of water. Flashing
/// the stream at (p, T) through the shipped evaluateBlackOil leaves rs sm3 of
/// gas still dissolved; the rest, gorStd - rs, is liberated. The liquid
/// outlet therefore holds oil + dissolved gas + water, the gas outlet holds
/// the liberated gas, and an imperfect machine (efficiency < 1) carries
/// (1 - efficiency) of the condensable mass into the gas outlet as entrained
/// liquid. Both outlets leave at p - dP: the drop prices the internals and is
/// reported, not marched.
///
/// Declared simplifications, all of them the point of the element rather
/// than accidents: one stage at one (p, T) (the cascade app owns multi-stage
/// pressure letdown); the liberated gas and the liquid share the stage
/// pressure (no flash-drum elevation); entrainment is a mass split, not a
/// droplet-size statement; the water is at its declared std density (the
/// black-oil layer books water incompressible).
struct SeparatorSolution
{
    double pressure = 0;            // Pa, stage pressure
    double temperature = 0;         // K, stage temperature
    double efficiency = 1;          // liquid recovery in [0,1]
    double outletPressure = 0;      // Pa, both outlets: pressure - dP
    double rs = 0;                  // sm3/sm3 solution GOR of the outlet liquid
    double gasLiberated = 0;        // sm3/sm3 free gas at standard conditions
    double bo = 1;                  // rm3/sm3 oil FVF at stage p,T
    double bg = 1;                  // rm3/sm3 free-gas FVF at stage p,T
    double z = 1;                   // gas deviation factor (DAK) at stage p,T
    double rhoGasStd = 0;           // kg/sm3 gas at standard conditions
    double feedMassRate = 0;        // kg/s in
    double liquidMassRateIdeal = 0; // kg/s condensables at efficiency = 1
    double entrainedMassRate = 0;   // kg/s (1-eff) of the condensables
    double liquidMassRate = 0;      // kg/s the liquid outlet actually carries
    double gasMassRate = 0;         // kg/s liberated gas + entrainment
    double closure = 0;             // kg/s liquid + gas - feed (0 by design)
};

/// Price one separator stage. Every input is a declared number; nothing here
/// iterates. Out-of-envelope declarations are refused loudly: a feed rate
/// that is not positive has no outlet stream, a stage pressure or
/// temperature that is not positive cannot flash anything, an efficiency
/// outside [0,1] is not a recovery, and a negative dP would deliver the
/// outlets above the stage. The fluid itself is validated by the same
/// statements the pvt cascade makes (positive gasGravity and oilDensity,
/// non-negative gor, waterCut in [0,1)).
SeparatorSolution separatorSolve(const pvt::BlackOilFluid& fluid,
                                 double feedRateSm3PerS,
                                 double pressurePa,
                                 double temperatureK,
                                 double efficiency,
                                 double dPPa);

} // namespace cf
