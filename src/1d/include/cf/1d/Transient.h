#pragma once

#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// M2 first increment: transient 1D network hydraulics with
///   * pipe inertance  rho*L/A * dQ/dt,
///   * node compliance V/K (K = mixture bulk modulus; K = 0 -> incompressible
///     port enforcing sum(q) = 0),
///   * Darcy-Weisbach friction (same Churchill law as the steady kernel),
///   * valve local loss dp = K rho |v| v / 2 (same K convention as the
///     steady kernel) with optional linear close/open schedules K(t),
///   * gravity head,
/// advanced by implicit Euler with a Newton solve on (Q, p), plus drift-flux
/// transport of the mobile liquid volume per node (slug tracking) and an
/// exact integrating-factor first-order hydrate kinetics substep.
/// Phase volumes per node: V = V_L (mobile liquid) + V_h (hydrated,
/// immobile) + V_G (gas fills the rest).

struct TransientNodeConfig
{
    std::string name;
    double volume = 0;           // m3, fixed geometric node volume
    double bulkModulus = 0;      // Pa, mixture bulk modulus (0 = incompressible port)
    bool pressureFixed = false;  // Dirichlet pressure
    double pressure = 0;         // fixed value, or initial value for free nodes
    double liquidFraction = 1.0; // initial mobile-liquid fraction of volume
    double hydrateFraction = 0.0;// initial hydrate fraction of total liquid volume
    double hydrateRate = 0.0;    // 1/s, first-order liquid -> hydrate rate
};

struct TransientPipeConfig
{
    std::string name;
    std::string fromName;
    std::string toName;
    double length = 0;
    double diameter = 0;
    double roughness = 0;
    double elevation = 0; // z(to) - z(from)
    double initialFlow = 0;
    double valveK = 0;           // open-valve loss coefficient dp = K rho v^2/2
    // linear schedule: K ramps valveK -> valveKClosed over
    // [valveCloseStart, valveCloseEnd], holds closed, then ramps back over
    // [valveOpenStart, valveOpenEnd]; negative start times disable the ramps
    double valveCloseStart = -1;
    double valveCloseEnd = -1;
    double valveKClosed = 0;
    double valveOpenStart = -1;
    double valveOpenEnd = -1;
};

struct TransientSource
{
    std::string nodeName;
    double rate = 0;     // m3/s total injected volume
    double fraction = 1; // liquid fraction of the injected volume
};

struct TransientCase
{
    FluidProperties fluid;
    std::vector<TransientNodeConfig> nodes;
    std::vector<TransientPipeConfig> pipes;
    std::vector<TransientSource> sources;
};

struct TransientSample
{
    double t = 0;
    std::vector<double> pressure;        // per node
    std::vector<double> liquidFraction;  // mobile V_L / V per node
    std::vector<double> hydrateFraction; // V_h / (V_L + V_h) per node
    std::vector<double> flowRate;        // per pipe
};

struct TransientSeries
{
    std::vector<std::string> nodeNames;
    std::vector<std::string> pipeNames;
    std::vector<TransientSample> samples;
};

/// dh/dt = k (1 - h) integrated exactly: h -> 1 - (1 - h0) exp(-k dt).
/// Hydrate volume grows at the expense of mobile liquid (volume preserving).
void hydrateKineticsStep(double& liquidVolume, double& hydrateVolume,
                         double rate, double dt);

class TransientNetwork
{
public:
    explicit TransientNetwork(const TransientCase& c);

    /// One implicit-Euler step of size dt; returns max scaled Newton residual.
    double step(double dt);

    double time() const { return t_; }
    const std::vector<std::string>& nodeNames() const { return net_.nodeNames; }
    const std::vector<std::string>& pipeNames() const;

    const std::vector<double>& pressure() const { return p_; }
    std::vector<double> pipeFlows() const { return q_; }
    std::vector<double> liquidFractions() const;
    std::vector<double> hydrateFractions() const;
    std::vector<double> liquidVolumes() const;
    std::vector<double> hydrateVolumes() const;

private:
    struct PipeDyn
    {
        NetworkPipe geo; // reuse the steady pipe geometry + friction law
        double qPrev = 0;
        double inertia = 0; // rho L / A
        // valve schedule: K ramps geo.valveK -> kClosed over
        // [closeStart, closeEnd] and back over [openStart, openEnd];
        // negative start times disable the ramps
        double closeStart = -1, closeEnd = -1, kClosed = 0;
        double openStart = -1, openEnd = -1;
    };
    struct NodeDyn
    {
        double volume = 0, bulk = 0;
        bool fixed = false;
        double pPrev = 0;
        double unk = -1;
        double vL = 0, vH = 0;
        double source = 0, sourceFraction = 0;
        double rate = 0;
    };

    OneDNetwork net_; // node names + resolved pipe endpoints (from/to)
    std::vector<std::string> pipeNames_;
    std::vector<PipeDyn> pipes_;
    std::vector<NodeDyn> nodes_;
    std::vector<double> p_;  // current node pressures (all nodes)
    std::vector<double> q_;  // current pipe flows
    FluidProperties fluid_;
    double t_ = 0;
};

/// Advance the case by endTime/dt steps, sampling every reportEvery steps
/// (sample 0 is the initial state).
TransientSeries runTransient(const TransientCase& c, double dt,
                             double endTime, int reportEvery);

/// Case grammar (see tests/cases/transientBlowdown/case.dict and
/// tests/cases/transientValve/case.dict):
///   fluid { rho 850; mu 0.1; }
///   nodes { TANK { volume 0.5; bulkModulus 2.0e8; pressure 5.0e5;
///                  liquidFraction 1; hydrateRate 0; }
///            OUT { fixedPressure 0; } }
///   pipes { P1 { from TANK; to OUT; length 100; diameter 0.05;
///                valve 20; valveClose { start 2; end 6; closedK 1.0e12; }
///                valveOpen { start 10; end 14; } } }
///   sources { s1 { node TANK; rate 0.0; fraction 1; } }
///   time { dt 0.01; endTime 10; reportEvery 1; }
TransientCase loadTransientFromDict(const Dictionary& caseDict);
struct TransientTimeControls
{
    double dt = 0.01;
    double endTime = 1.0;
    int reportEvery = 1;
};
TransientTimeControls loadTransientControls(const Dictionary& caseDict);

} // namespace cf
