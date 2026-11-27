"""Transient 1D<->3D coupling family Python binding tests (CTest
`python.coupling_transient`).

Covers `chflow_core.solve_transient_coupled_case`:

  * single-edge transient coupled case (`transientCoupled1d3d`): the whole
    per-step series, the 1D state at every level, the interface exchanges and
    the partner block's fields;
  * multi-edge transient coupled case (`transientCoupledMulti1d3d`): one
    exchange per edge at every level, series-edge topology pins its rate;
  * failure paths (missing case file, seed index out of range) raise loudly.

The numeric 口径 is the one from `tests/check_transient_coupled_1d3d.py` and
`tests/check_transient_coupled_multi.py` (same closed forms, same tolerances),
not a looser home-grown threshold.
"""
import math
import os
import re
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

import chflow_core

ROOT = Path(__file__).resolve().parents[2]
SINGLE_CASE = ROOT / "tests" / "cases" / "transientCoupled1d3d" / "case.dict"
MULTI_CASE = ROOT / "tests" / "cases" / "transientCoupledMulti1d3d" / "case.dict"

# Geometry / fluid of tests/cases/transientCoupled1d3d/case.dict, mirrored from
# check_transient_coupled_1d3d.py so the closed form and tolerances stay the
# same source.
MU, RHO = 0.1, 800.0
H, W, L = 0.02, 0.02, 0.5
Q_IN = 2.4e-4
P_REF = 1.0e6
DT, N_STEPS = 0.02, 60
TOL_INTERFACE = 1e-4

SEC_X = 0.25 * math.pi * 0.02 ** 2
SEC_BY = 0.25 * math.pi * 0.03 ** 2
SEC_LEAD = 0.25 * math.pi * 0.1 ** 2
M_X = RHO * 0.5 / SEC_X
M_BY = RHO * 1.5 / SEC_BY
M_LEAD = RHO * 5.0 / SEC_LEAD
R_BY = 128.0 * MU * 1.5 / (math.pi * 0.03 ** 4)

# Multi-edge case constants (check_transient_coupled_multi.py).
Q_IN_M = 3.0e-5
DT_M, N_STEPS_M = 0.05, 12


class TransientCoupledSingleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = chflow_core.solve_transient_coupled_case(str(SINGLE_CASE), 0)

    def test_series_shape_and_clock(self):
        self.assertEqual(self.result.node_names, ["Src", "J1", "J2", "Snk"])
        self.assertEqual(self.result.pipe_names, ["A", "X", "By", "B"])
        self.assertEqual(len(self.result.steps), N_STEPS)
        for n, st in enumerate(self.result.steps, start=1):
            self.assertLess(abs(st.t - n * DT), 1e-12)
            self.assertEqual(len(st.pressure), 4)
            self.assertEqual(len(st.flow_rate), 4)
            self.assertEqual(len(st.edges), 1)
            self.assertTrue(math.isfinite(st.thermal_residual))

    def test_per_step_continuity_and_interface(self):
        for st in self.result.steps:
            f = dict(zip(self.result.pipe_names, st.flow_rate))
            p = dict(zip(self.result.node_names, st.pressure))
            ex = st.edges[0]
            self.assertLessEqual(ex.residual, TOL_INTERFACE)
            self.assertTrue(ex.block_converged)
            self.assertLess(abs(f["A"] - Q_IN), 1e-12 * Q_IN)
            self.assertLess(abs(f["A"] - f["X"] - f["By"]), 1e-12 * Q_IN)
            self.assertLess(abs(f["B"] - f["X"] - f["By"]), 1e-12 * Q_IN)
            self.assertLess(abs(p["Snk"] - P_REF), 1e-9)
            # the reported rate is the coupled pipe's own flow
            self.assertAlmostEqual(ex.q, f["X"], delta=1e-12 * Q_IN)

    def test_ramp_monotonic_and_settled(self):
        qx = [st.edges[0].q for st in self.result.steps]
        qby = [dict(zip(self.result.pipe_names, st.flow_rate))["By"]
               for st in self.result.steps]
        self.assertTrue(all(qx[i] < qx[i + 1] for i in range(N_STEPS - 1)))
        self.assertTrue(all(qby[i] > qby[i + 1] for i in range(N_STEPS - 1)))

        r_x = self.result.state[0].resistance
        q_inf = R_BY * Q_IN / (r_x + R_BY)
        gap = (q_inf - qx[-1]) / q_inf
        self.assertGreater(gap, 0.0)
        self.assertLessEqual(gap, 0.005)

        # first step: the at-rest network splits the sudden rate by inertance
        q1 = Q_IN * (R_BY + M_BY / DT) / ((M_X + M_BY) / DT + r_x + R_BY)
        self.assertLess(abs(qx[0] / q1 - 1.0), 1e-3)
        self.assertLess(qx[0], 0.95 * q_inf)

        # geometric ramp: q^n = A + B q^(n-1), destination and decay per step
        xs, ys = qx[1:-1], qx[2:]
        k = len(xs)
        sx, sy = sum(xs), sum(ys)
        b = (k * sum(p * f for p, f in zip(xs, ys)) - sx * sy) / (
            k * sum(v * v for v in xs) - sx * sx)
        a = (sy - b * sx) / k
        q_inf_fit = a / (1.0 - b)
        self.assertGreater(b, 0.9)
        self.assertLess(b, 1.0)
        self.assertLess(abs(q_inf_fit / q_inf - 1.0), 0.01)

    def test_block_measure_and_fields(self):
        measures = self.result.block_measures
        self.assertEqual(len(measures), 1)
        m = measures[0]
        self.assertEqual(m.cells, 1000)
        self.assertTrue(m.converged)
        self.assertLess(abs(m.flux_out - m.flux_in), 1e-7 * m.flux_in)
        self.assertLess(abs(m.ux_centre / m.plug_speed - 1.5), 0.03)

        u = self.result.block_velocity(0)
        p = self.result.block_pressure(0)
        self.assertEqual(u.shape, (1000, 3))
        self.assertEqual(p.shape, (1000,))
        self.assertEqual(u.dtype, np.float64)
        self.assertEqual(p.dtype, np.float64)
        self.assertTrue(u.flags["C_CONTIGUOUS"])
        self.assertTrue(p.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(u).all())
        self.assertTrue(np.isfinite(p).all())

    def test_state_matches_last_step(self):
        last = self.result.steps[-1]
        self.assertAlmostEqual(self.result.state[0].q, last.edges[0].q, places=12)
        self.assertAlmostEqual(self.result.state[0].dp_pa, last.edges[0].dp_pa,
                               places=12)

    def test_out_of_range_raises(self):
        with self.assertRaises(IndexError):
            chflow_core.solve_transient_coupled_case(str(SINGLE_CASE), 99)

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_transient_coupled_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


class TransientCoupledMultiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = chflow_core.solve_transient_coupled_case(str(MULTI_CASE), 0)

    def test_two_edges_at_every_level(self):
        self.assertEqual(self.result.edge_names, ["X1", "X2"])
        self.assertEqual(len(self.result.steps), N_STEPS_M)
        for st in self.result.steps:
            self.assertEqual(len(st.edges), 2)
            f = dict(zip(self.result.pipe_names, st.flow_rate))
            self.assertLess(abs(f["A"] - Q_IN_M), 1e-12 * Q_IN_M)
            self.assertLess(abs(f["X2"] - Q_IN_M), 1e-12 * Q_IN_M)
            self.assertLess(abs(f["B"] - Q_IN_M), 1e-12 * Q_IN_M)
            self.assertLess(abs(f["A"] - f["X1"] - f["By"]), 1e-12 * Q_IN_M)
            for ex in st.edges:
                self.assertLessEqual(ex.residual, TOL_INTERFACE)
                self.assertTrue(ex.block_converged)

    def test_block_measures_and_fields(self):
        measures = self.result.block_measures
        self.assertEqual(len(measures), 2)
        self.assertEqual([m.cells for m in measures], [1000, 800])
        for i, n in enumerate((1000, 800)):
            u = self.result.block_velocity(i)
            p = self.result.block_pressure(i)
            self.assertEqual(u.shape, (n, 3))
            self.assertEqual(p.shape, (n,))
            self.assertEqual(u.dtype, np.float64)
            self.assertEqual(p.dtype, np.float64)
            self.assertTrue(np.isfinite(u).all())
            self.assertTrue(np.isfinite(p).all())

    def test_cli_binary_parity_when_available(self):
        """Compare seed-0 (`lo`) last-step rates against the multi-edge CLI
        solver when its built binary is present; skipped otherwise."""
        candidates = [
            ROOT / "build-wt" / "src" / "coupling" / "apps" /
            "chflow-coupled-1d3d-transient-multi",
            ROOT / "build" / "src" / "coupling" / "apps" /
            "chflow-coupled-1d3d-transient-multi",
        ]
        binary = next((p for p in candidates if p.exists()), None)
        if binary is None:
            self.skipTest("transient multi-edge CLI binary not built")

        r = subprocess.run([str(binary), str(MULTI_CASE)],
                           capture_output=True, text=True,
                           cwd=str(MULTI_CASE.parent))
        self.assertEqual(r.returncode, 0, (r.stdout + r.stderr)[-400:])
        cli_q = {}
        for ln in r.stdout.splitlines():
            m = re.search(r"transient coupled start lo edge (\S+) name (\S+) q "
                          r"([-+0-9.eE]+)", ln)
            if m:
                cli_q[m.group(2)] = float(m.group(3))
        for ex in self.result.state:
            self.assertIn(ex.name, cli_q)
            self.assertLess(abs(ex.q - cli_q[ex.name]) / cli_q[ex.name], 1e-6)


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
