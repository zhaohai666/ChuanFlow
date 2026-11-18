#pragma once

#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// D8 increment 6: quadrature-based moment method (QMOM) for particle
/// coagulation. The population balance is tracked through its first four
/// moments M_k = integral n(v) v^k dv (k=0..3), and a two-node Gaussian
/// quadrature is reconstructed from them (the standard QMOM closure, McGraw
/// 1997; Marchisio & Fox "Computational Models for Polydisperse Particulate
/// and Multiphase Systems", ch. 5).
///
/// The coagulation closure is the Smoluchowski constant-kernel beta:
///   dM_0/dt = -(beta/2) M_0^2              (number falls: 2 -> 1 per event)
///   dM_1/dt = 0                            (mass conserved exactly)
///   dM_2/dt = beta M_1^2
///   dM_3/dt = 3 beta M_1 M_2
/// These follow from dM_k/dt = (beta/2) sum_{p=1}^{k-1} C(k,p) M_p M_{k-p}
/// plus the separate number loss -beta M_0^2/2 (Friedlander, "Smoke, Dust and
/// Haze"; Ramabhadran & Seinfeld). With a constant kernel the system has the
/// closed form the checker reconstructs independently:
///   M_0(t) = M_00 / (1 + beta M_00 t / 2)
///   M_1(t) = M_10
///   M_2(t) = M_20 + beta M_10^2 t
///   M_3(t) = M_30 + 3 beta M_10 M_20 t + (3/2) beta^2 M_10^2 t^2
///
/// The two-node quadrature is the Gaussian rule matching M_0..M_3: its nodes
/// are the eigenvalues of the Jacobi matrix built from the moments and its
/// weights reproduce all four moments. The defining identity sum_i w_i v_i^k
/// = M_k (k=0..3) is what the checker grades, so it does not trust the
/// reconstruction algorithm.
struct QmomConfig
{
    bool enabled = false;
    scalar beta = 0;       // coagulation kernel [m3/s] (constant)
    scalar m0 = 0;         // initial M_0 (number density, 1/m3)
    scalar m1 = 0;         // initial M_1 (volume fraction, dimensionless)
    scalar m2 = 0;         // initial M_2 (m3/m3)
    scalar m3 = 0;         // initial M_3 (m6/m3)
    scalar dt = 0;         // s
    int steps = 0;
    int reportEvery = 1;

    void validate() const;
};

QmomConfig readQmomConfig(const Dictionary& caseDict);

struct QmomStep
{
    int step = 0;
    scalar t = 0;
    scalar m[4] = {0, 0, 0, 0};
    scalar w[2] = {0, 0};  // quadrature weights
    scalar v[2] = {0, 0};  // quadrature nodes (particle volumes)
    scalar massError = 0;  // |M_1(t) - M_1(0)| / M_1(0)
};

struct QmomResult
{
    std::vector<QmomStep> history;
    bool valid = false;
};

/// Constant-kernel Smoluchowski moment source, exposed so a checker can
/// reconstruct every step with the same published law.
void coagulationSource(const scalar m[4], scalar beta, scalar dm[4]);

/// Reconstruct a two-node Gaussian quadrature from M_0..M_3. Returns false
/// when the moments are not a valid measure (e.g. negative variance).
bool quadratureFromMoments(const scalar m[4], scalar w[2], scalar v[2]);

QmomResult runQmom(const QmomConfig& cfg);

} // namespace cf
