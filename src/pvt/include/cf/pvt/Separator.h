#pragma once

#include <vector>

#include "cf/pvt/BlackOil.h"

namespace cf::pvt
{

/// One stage of a black-oil separator cascade, per sm3 of stock-tank oil.
///
/// All gas volumes are at standard conditions (sm3/sm3 oil); all oil volumes
/// are at the stage's p,T (rm3/sm3). The gas liberated at a stage is the drop
/// in solution GOR across that stage: the liquid entering carries
/// `availableGas`, the liquid leaving carries `rs`, and the difference leaves
/// as free gas. Every stage re-flashes the *remaining* liquid, so its bubble
/// point falls as gas is removed -- that is the only state handed onward.
struct SeparatorStage
{
    int index = 0;
    double pressure = 0;      // Pa, stage pressure
    double temperature = 0;   // K, stage temperature
    double availableGas = 0;  // sm3/sm3 solution GOR of the incoming liquid
    double rs = 0;            // sm3/sm3 solution GOR of the outgoing liquid
    double gasLiberated = 0;  // sm3/sm3 free gas released here (std conditions)
    double bo = 1;            // rm3/sm3 oil FVF at stage p,T
    double bg = 1;            // rm3/sm3 free-gas FVF at stage p,T
    double z = 1;             // gas deviation factor (DAK) at stage p,T
    double rhoOil = 0;        // kg/m3 live oil at stage p,T
    double gasMass = 0;       // kg/sm3 mass of gas released here
    double liquidMass = 0;    // kg/sm3 liquid leaving (oil + dissolved gas + water)
    double cumGasMass = 0;    // kg/sm3 cumulative gas removed through this stage
    double liveApi = 0;       // API gravity of the stage liquid (from rhoOil)
};

/// Result of a separator cascade, per sm3 of stock-tank oil.
struct SeparatorResult
{
    std::vector<SeparatorStage> stages;
    double stockTankApi = 0; // deg API from oilDensityStd (constant)
    double gasMassPerSm3 = 0;// kg/sm3 gas at standard conditions
    double massIn = 0;       // kg/sm3 feed mass (blackOilMassPerOilSm3)
    double massOut = 0;      // kg/sm3 final liquid + all gas removed
    double massImbalance = 0;// kg/sm3 = massOut - massIn
    double totalGasSm3 = 0;  // sm3/sm3 gas removed over all stages
};

/// Stateless N-stage black-oil separator cascade.
///
/// `pressures` are the stage pressures in Pa (any order is accepted; a gas
/// release only happens where solution GOR falls). `temperatures[i]` is the
/// temperature of stage i in K. No hidden state is carried between calls:
/// the only thing passed from one stage to the next is the outgoing solution
/// GOR `rs`, computed by evaluateBlackOil on a fluid whose gorStd has been
/// reset to the gas still in solution.
///
/// Mass closure is exact by construction and is reported for verification:
///   liquidMass_i + cumGasMass_i == blackOilMassPerOilSm3(fluid)
/// because liquidMass_i = m_oil + rs_i*m_gas + wor*m_water and
/// cumGasMass_i = (gorStd - rs_i)*m_gas.
SeparatorResult separateBlackOil(const BlackOilFluid& fluid,
                                 const std::vector<double>& pressures,
                                 const std::vector<double>& temperatures);

/// Isothermal convenience overload: every stage runs at the same T.
SeparatorResult separateBlackOil(const BlackOilFluid& fluid,
                                 const std::vector<double>& pressures,
                                 double temperature);

} // namespace cf::pvt
