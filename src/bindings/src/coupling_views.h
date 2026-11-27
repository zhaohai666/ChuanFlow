// Shared read-back views for the steady and transient 1D<->3D coupling
// families. The kernel structs (`EdgeExchange`, `BlockMeasure`) are copied into
// these plain POD views so Python owns the numbers; field names track the
// kernel struct (camelCase -> snake_case on the pybind side), not the app's log
// wording.
#pragma once

#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/coupling/CoupledNetwork.h"

namespace cfbind
{

/// A snapshot of one edge's last interface exchange, copied out of the kernel
/// `EdgeExchange` so Python owns the numbers.
struct EdgeExchangeView
{
    cf::label edge = -1;
    std::string name;
    double resistance = 0;
    double q = 0;
    double dpPa = 0;
    double residual = 0;
    cf::label sweeps = 0;
    bool blockConverged = false;
    double inertance = 0;
    bool unsteady = false;
    double dpFrictionPa = 0;
    double dpInertialPa = 0;
    double qPrev = 0;
    double dq = 0;
    double rowDropPa = 0;
    double rowError = 0;
    double chordPa = 0;
    std::string inletShape = "flat";
    double fluxImposed = 0;
    double betaImposed = 1;
    double betaIn = 1;
    double betaOut = 1;
    double inletArea = 0;
    double pRisePa = 0;
    double pRiseX = 0;
    double dpConvectivePa = 0;
    double convectiveRowPa = 0;
    bool thermalOn = false;
    double inletTemperature = 0;
    double wallHeatRate = 0;
    double inletHeatRate = 0;
    double enthalpyOut = 0;
    double lineHeatFlux = 0;
    double t1d = 0;
    double thermalResidual = 0;
    double thermalBalance = 0;
};

inline EdgeExchangeView toView(const cf::coupling::EdgeExchange& ex)
{
    EdgeExchangeView v;
    v.edge = ex.edge;
    v.name = ex.name;
    v.resistance = ex.resistance;
    v.q = ex.q;
    v.dpPa = ex.dpPa;
    v.residual = ex.residual;
    v.sweeps = ex.sweeps;
    v.blockConverged = ex.blockConverged;
    v.inertance = ex.inertance;
    v.unsteady = ex.unsteady;
    v.dpFrictionPa = ex.dpFrictionPa;
    v.dpInertialPa = ex.dpInertialPa;
    v.qPrev = ex.qPrev;
    v.dq = ex.dq;
    v.rowDropPa = ex.rowDropPa;
    v.rowError = ex.rowError;
    v.chordPa = ex.chordPa;
    v.inletShape = ex.inletShape;
    v.fluxImposed = ex.fluxImposed;
    v.betaImposed = ex.betaImposed;
    v.betaIn = ex.betaIn;
    v.betaOut = ex.betaOut;
    v.inletArea = ex.inletArea;
    v.pRisePa = ex.pRisePa;
    v.pRiseX = ex.pRiseX;
    v.dpConvectivePa = ex.dpConvectivePa;
    v.convectiveRowPa = ex.convectiveRowPa;
    v.thermalOn = ex.thermalOn;
    v.inletTemperature = ex.inletTemperature;
    v.wallHeatRate = ex.wallHeatRate;
    v.inletHeatRate = ex.inletHeatRate;
    v.enthalpyOut = ex.enthalpyOut;
    v.lineHeatFlux = ex.lineHeatFlux;
    v.t1d = ex.t1d;
    v.thermalResidual = ex.thermalResidual;
    v.thermalBalance = ex.thermalBalance;
    return v;
}

/// A snapshot of one block's own last solve (`BlockMeasure`), the numbers the
/// app prints on the `block seed` line.
struct BlockMeasureView
{
    cf::label cells = 0;
    cf::label sweeps = 0;
    bool converged = false;
    double dpPa = 0;
    double imposedQ = 0;
    double fluxIn = 0;
    double fluxOut = 0;
    double planeFluxMin = 0;
    double planeFluxMax = 0;
    double probeSlope = 0;
    double probeSpread = 0;
    cf::label probeSections = 0;
    double inletArea = 0;
    double plugSpeed = 0;
    double uxMax = 0;
    double uxCentre = 0;
    std::string turbModel;
    double nutMax = 0;
    double wallYPlusMax = 0;
    double wallUTau = 0;
    std::string inletShape;
    double fluxImposed = 0;
    double betaImposed = 0;
    double betaIn = 0;
    double betaOut = 0;
    std::vector<std::pair<double, double>> betaStations;
    std::vector<std::pair<double, double>> pressureStations;
    // Unsteady-partner diagnostics (zero/default on quasi-steady blocks).
    bool unsteady = false;
    double dpFrictionPa = 0;
    double dpInertialPa = 0;
    double inertanceDynamic = 0;
    double inertanceProbe = 0;
    double inertance = 0;
    double betaProfile = 0;
    double momentumFluxImposedPa = 0;
    double dtLevel = 0;
};

inline BlockMeasureView toView(const cf::coupling::BlockMeasure& m)
{
    BlockMeasureView v;
    v.cells = m.cells;
    v.sweeps = m.sweeps;
    v.converged = m.converged;
    v.dpPa = m.dpPa;
    v.imposedQ = m.imposedQ;
    v.fluxIn = m.fluxIn;
    v.fluxOut = m.fluxOut;
    v.planeFluxMin = m.planeFluxMin;
    v.planeFluxMax = m.planeFluxMax;
    v.probeSlope = m.probeSlope;
    v.probeSpread = m.probeSpread;
    v.probeSections = m.probeSections;
    v.inletArea = m.inletArea;
    v.plugSpeed = m.plugSpeed;
    v.uxMax = m.uxMax;
    v.uxCentre = m.uxCentre;
    v.turbModel = m.turbModel;
    v.nutMax = m.nutMax;
    v.wallYPlusMax = m.wallYPlusMax;
    v.wallUTau = m.wallUTau;
    v.inletShape = m.inletShape;
    v.fluxImposed = m.fluxImposed;
    v.betaImposed = m.betaImposed;
    v.betaIn = m.betaIn;
    v.betaOut = m.betaOut;
    v.betaStations = m.betaStations;
    v.pressureStations = m.pressureStations;
    v.unsteady = m.unsteady;
    v.dpFrictionPa = m.dpFrictionPa;
    v.dpInertialPa = m.dpInertialPa;
    v.inertanceDynamic = m.inertanceDynamic;
    v.inertanceProbe = m.inertanceProbe;
    v.inertance = m.inertance;
    v.betaProfile = m.betaProfile;
    v.momentumFluxImposedPa = m.momentumFluxImposedPa;
    v.dtLevel = m.dtLevel;
    return v;
}

} // namespace cfbind
