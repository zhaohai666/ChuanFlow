#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// Constant-property single-phase fluid for the 1D network kernel.
/// Full PVT flash models will live in src/pvt and replace this later.
struct FluidProperties
{
    double rho = 850.0; // kg/m3
    double mu = 1e-3;   // Pa.s
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

struct PipeFlowResult
{
    std::string name;
    double dpTotal = 0;     // p(from) - p(to), Pa
    double dpFriction = 0;  // dpTotal - rho*g*dz
    double velocity = 0;    // m/s, positive from -> to
    double flowRate = 0;    // m3/s
    double reynolds = 0;
    double frictionFactor = 0;
    double conductance = 0; // dQ/d(dpTotal), s
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

/// Signed flow through one pipe for a given total pressure difference
/// dp = p(from) - p(to). Q/dpTotal conductance via central differencing.
void solvePipeFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                   double dpTotal, PipeFlowResult& out);

SteadyNetworkResult solveSteadyNetwork(const OneDNetwork& net,
                                       const SteadyControls& ctl);

/// Case dictionary grammar (see tests/cases/steadyNetwork/case.dict):
///   fluid { rho 850; mu 0.1; }
///   nodes { W1 {}; Outlet { pressure 1.0e6; } }
///   pipes { P1 { from W1; to Outlet; length 1000; diameter 0.05;
///                roughness 0; elevation 5; } }
///   sources { s1 { node W1; rate 0.002; } }
///   solve { maxIter 100; tol 1e-10; relax 1.0; }
OneDNetwork loadNetworkFromDict(const Dictionary& caseDict);
SteadyControls loadSteadyControls(const Dictionary& caseDict);

} // namespace cf
