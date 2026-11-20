#pragma once

#include <string>
#include <vector>

#include "cf/1d/BeggsBrill.h"
#include "cf/common/Dictionary.h"

namespace cf
{

/// Inflow performance relationship of a producer: bottom-hole flowing
/// pressure p_wf as a function of the stock-tank oil rate q. Three classic
/// closed-form models, SI throughout (Pa, sm3/s of stock-tank oil):
///   linear (Darcy inflow, single-phase or undersaturated oil, p_wf > pb):
///     q = J (p_res - p_wf),  J = productivityIndex
///   vogel (solution-gas drive, Vogel 1968):
///     q/qmax = 1 - 0.2 x - 0.8 x^2 with x = p_wf / p_res, qmax the AOF
///   jones (Rawlins-Schellhardt style deliverability):
///     q = C (p_res^2 - p_wf^2)^n,  C = deliverability
/// Each form carries its own closed-form inverse; pWf() clamps the rate
/// into the form's domain [0, rateMax()] and rate() the pressure into
/// [0, p_res], so the curves are total functions on their bracket.
struct IpR
{
    enum class Form { linear, vogel, jones };

    Form form = Form::linear;
    double reservoirPressure = 0; // Pa, average reservoir pressure
    double productivityIndex = 0; // linear: sm3/(s Pa)
    double maxRate = 0;           // vogel: absolute open flow, sm3/s
    double deliverability = 0;    // jones: C, sm3/(s Pa^(2n))
    double exponent = 1;          // jones: n

    /// Rate at which p_wf reaches 0: J p_res / qmax / C p_res^(2n).
    double rateMax() const;
    double pWf(double rate) const;  // Pa; the closed-form inverse
    double rate(double pWf) const;  // sm3/s; the curve itself
};

/// The wellbore as one Beggs-Brill pipe: measured depth along the flow path
/// plus the elevation rise from bottom hole to wellhead.
struct WellLift
{
    double length = 0;    // m
    double diameter = 0;  // m
    double roughness = 0; // m
    double elevation = 0; // m, z(wellhead) - z(bottomhole)
};

/// Vertical lift performance: the bottom-hole pressure that carries `rate`
/// of stock-tank oil (with the fluid's GOR / water cut / temperature) from
/// bottom hole up to a wellhead held at `wellheadPressure`. The producer's
/// stream flows upward, so the march runs bottom-up: the returned pressure
/// is the inlet whose Beggs-Brill march lands exactly on the wellhead
/// pressure (expanding bracket + Illinois secant, the same solver family as
/// the network pipe law - pOutlet rises sublinearly in pIn because higher
/// pressure packs the free gas smaller and eases the friction).
double vlpBottomHolePressure(const pvt::BlackOilFluid& fluid,
                             const WellLift& well, double rate,
                             double wellheadPressure);

struct NodalControls
{
    int maxIter = 100;   // secant iterations
    double tol = 1e-9;   // relative rate tolerance of the final bracket
};

struct NodalResult
{
    double rate = 0;               // sm3/s stock-tank oil at the operating point
    double bottomHolePressure = 0; // Pa
    double staticBottomPressure = 0; // Pa, the rate-0 column, for diagnosis
    double holdup = 0;             // outlet liquid holdup at the operating point
    int regime = -1;               // Beggs-Brill regime of the last segment
    int iterations = 0;
    bool flows = false;  // false when the reservoir cannot lift the static
                         // column: the rate stays 0 and the result is the
                         // static column itself
    bool converged = false;
};

/// Nodal analysis (the Pipesim well planner): the operating point is where
/// the IPR and the VLP cross, i.e. the largest root of
///   f(q) = p_vlp(q) - p_ipr(q).
/// p_vlp is U-shaped in q - holdup relief lowers the column before friction
/// raises it - and p_ipr falls monotonically, so f can vanish more than
/// once; the stable operating point is the largest root, bracketed from the
/// high-rate end exactly like the network pipe law. A well whose reservoir
/// cannot overcome the static column reports flows=false with rate 0.
NodalResult solveNodalAnalysis(const IpR& ipr, const pvt::BlackOilFluid& fluid,
                               const WellLift& well, double wellheadPressure,
                               const NodalControls& ctl = {});

/// 5.157 (P5): one lateral of a multi-lateral well. A lateral carries its own
/// IPR (any of the three closed forms, independently of the main bore and of
/// its siblings) and its own VLP bore (the same Beggs-Brill family as the
/// main wellbore). Its top opens into the main bore's bottom hole, so the
/// junction pressure - the main bore's flowing bottom-hole pressure - is the
/// lateral's wellhead back-pressure.
struct LateralSpec
{
    std::string name;
    IpR ipr;
    WellLift well;
};

/// 5.157: one lateral's operating point at the solved junction pressure.
struct LateralResult
{
    std::string name;
    double rate = 0;                 // sm3/s stock-tank oil
    double bottomHolePressure = 0;   // Pa, the lateral's own bottom hole
    double staticBottomPressure = 0; // Pa, its rate-0 column at the junction
    double holdup = 0;
    int regime = -1;
    bool flows = false;
    bool converged = false;
};

/// 5.157: the multi-lateral nodal solve. `main` carries the total rate
/// Q = sum(q_i) and the junction (main bottom-hole) pressure; `laterals`
/// holds one entry per declared lateral, in declaration order.
struct MultiLateralResult
{
    NodalResult main;
    std::vector<LateralResult> laterals;
    int iterations = 0;
    bool converged = false;
};

/// Multi-lateral nodal analysis (the Pipesim multi-lateral completion): the
/// junction pressure p_j is the fixed point of
///   F(p_j) = p_j - p_vlp_main( sum_i q_i(p_j) ),
/// where q_i(p_j) is lateral i's own IPR x VLP operating point at back-pressure
/// p_j (solveNodalAnalysis on the lateral's bore). Each q_i falls with p_j and
/// the main VLP rises with the total, so F rises and the bracket has exactly
/// one root; the solve is the same expanding-bracket + Illinois secant family
/// as the single-well nodal and the network pipe law. A lateral whose
/// reservoir cannot lift its own column against p_j contributes rate 0.
MultiLateralResult solveNodalLaterals(const pvt::BlackOilFluid& fluid,
                                      const WellLift& mainWell,
                                      double wellheadPressure,
                                      const std::vector<LateralSpec>& laterals,
                                      const NodalControls& ctl = {});

/// One nodal-analysis case parsed from a dictionary:
///   fluid { blackoil { gasGravity 0.65; oilDensity 850; gor 90;
///                      waterCut 0; temperature 340; surfaceTension 0.03; } }
///   ipr { form vogel; reservoirPressure 2.2e7; maxRate 0.006; }
///        // linear: form linear; reservoirPressure; productivityIndex
///        // jones:  form jones; reservoirPressure; deliverability; exponent
///   well { length 2400; diameter 0.062; roughness 1e-5; elevation 2200; }
///   wellheadPressure 2.5e6;
///   solve { maxIter 100; tol 1e-9; }
///   laterals { L1 { ipr { ... } well { ... } } L2 { ... } }  // 5.157, optional
/// When laterals{} is declared the case is a multi-lateral completion: the
/// main bore is a pure lift conduit, the total rate is the sum of the lateral
/// rates and the main ipr{} no longer feeds the solve (the header still
/// reports it, and the lateral lines announce the switch).
struct NodalCase
{
    pvt::BlackOilFluid fluid;
    IpR ipr;
    WellLift well;
    double wellheadPressure = 0;
    NodalControls controls;
    std::vector<LateralSpec> laterals; // 5.157: empty = the single-well case
};

NodalCase loadNodalCaseFromDict(const Dictionary& caseDict);

} // namespace cf
