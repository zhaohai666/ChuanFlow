"""Coupling-family Python binding tests (CTest `python.coupling`).

Covers the 1D<->3D coupled steady solve exposed through `chflow.solve_coupled_case`:

  * read case.dict -> run the partitioned interface solve -> structured result;
  * the operating point matches the same analytic closed form the CLI checker
    `tests/check_coupled_1d3d.py` uses (planar Poiseuille block resistance in
    parallel with the 1D bypass), with the checker's own tolerances;
  * per-edge exchange, network node pressures and pipe flows, and per-block
    measures are retrievable independently;
  * multi-edge cases expose one exchange per edge;
  * failure paths (missing case file, seed index out of range) raise loudly.

The numeric 口径 is deliberately the one from `check_coupled_1d3d.py`, not a
looser home-grown threshold.
"""
import math
import os
import re
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

import chflow

ROOT = Path(__file__).resolve().parents[2]
SINGLE_CASE = ROOT / "tests" / "cases" / "coupled1d3d" / "case.dict"
MULTI_CASE = ROOT / "tests" / "cases" / "coupled1d3dMulti" / "case.dict"

# Geometry / fluid of tests/cases/coupled1d3d/case.dict (mirrored from
# check_coupled_1d3d.py so the closed form and tolerances stay the same source).
MU, RHO = 0.1, 800.0
H, W, L = 0.02, 0.02, 0.5
Q_IN = 2.4e-4
L_BY, D_BY = 1.5, 0.03
TOL_Q = 0.04          # checker: 4% on the split flow
TOL_DP = 0.08         # checker: 8% on the block drop
TOL_RESIDUAL = 1e-4   # couple{} tolerance the case asks for
TOL_CONTINUITY = 0.02 # checker: 2% block in/out imbalance


def analytic_operating_point():
    r_block = 12.0 * MU * L / (H ** 3 * W)
    r_by = 128.0 * MU * L_BY / (math.pi * D_BY ** 4)
    q_block = Q_IN * r_by / (r_block + r_by)
    dp_block = r_block * q_block
    return q_block, dp_block


class CoupledSingleEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = chflow.solve_coupled_case(str(SINGLE_CASE), 0)

    def test_run_converged(self):
        self.assertTrue(self.result.converged)
        self.assertGreaterEqual(self.result.iterations, 1)
        self.assertLessEqual(self.result.residual, TOL_RESIDUAL)

    def test_exchange_shape_matches_kernel(self):
        self.assertEqual(self.result.edge_names, ["X"])
        exchanges = self.result.exchanges
        self.assertEqual(len(exchanges), 1)
        ex = exchanges[0]
        self.assertGreaterEqual(ex.edge, 0)
        self.assertEqual(ex.name, "X")
        # Exchange carries the same fields the kernel EdgeExchange does.
        for attr in ("resistance", "q", "dp_pa", "residual", "sweeps",
                     "block_converged", "dp_friction_pa", "dp_inertial_pa",
                     "inlet_area", "beta_in", "beta_out"):
            self.assertTrue(hasattr(ex, attr), attr)

    def test_operating_point_matches_closed_form(self):
        q_exact, dp_exact = analytic_operating_point()
        ex = self.result.exchanges[0]
        self.assertLess(abs(ex.q - q_exact) / q_exact, TOL_Q)
        self.assertLess(abs(ex.dp_pa - dp_exact) / dp_exact, TOL_DP)

    def test_network_results_structured(self):
        self.assertEqual(len(self.result.node_pressure), 4)
        self.assertEqual(len(self.result.pipe_flows), 4)
        by_name = {f.name: f for f in self.result.pipe_flows}
        self.assertIn("X", by_name)
        self.assertIn("By", by_name)
        # The coupled pipe carries the exchange rate the block was measured at.
        self.assertAlmostEqual(by_name["X"].flow_rate, self.result.exchanges[0].q,
                               places=12)
        # Parallel branches share their drop; leads split at the inlet rate.
        self.assertAlmostEqual(by_name["X"].dp_total, by_name["By"].dp_total,
                               places=3)
        self.assertAlmostEqual(by_name["A"].flow_rate, Q_IN, places=12)

    def test_block_measures_and_fields(self):
        measures = self.result.block_measures
        self.assertEqual(len(measures), 1)
        m = measures[0]
        self.assertEqual(m.cells, 1000)
        self.assertTrue(m.converged)
        # Block continuity: what enters leaves the detailed region.
        self.assertLess(abs(m.flux_in - m.flux_out) / m.flux_in,
                        TOL_CONTINUITY)

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

    def test_history_per_iteration(self):
        self.assertGreaterEqual(len(self.result.history), 1)
        for step in self.result.history:
            self.assertIn("iteration", step)
            self.assertEqual(len(step["edges"]), 1)

    def test_out_of_range_raises(self):
        with self.assertRaises(IndexError):
            chflow.solve_coupled_case(str(SINGLE_CASE), 99)

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow.solve_coupled_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


class CoupledMultiEdgeTests(unittest.TestCase):
    def test_multi_edge_exchanges_indexable(self):
        result = chflow.solve_coupled_case(str(MULTI_CASE), 0)
        self.assertTrue(result.converged)
        self.assertEqual(result.edge_names, ["X1", "X2"])
        self.assertEqual(len(result.exchanges), 2)
        self.assertEqual([e.name for e in result.exchanges], ["X1", "X2"])
        self.assertEqual(len(result.block_measures), 2)
        for ex in result.exchanges:
            self.assertGreater(ex.q, 0.0)
            self.assertTrue(math.isfinite(ex.dp_pa))
            self.assertTrue(math.isfinite(ex.residual))

    def test_cli_binary_parity_when_available(self):
        """Compare the python numbers against the multi-edge CLI solver when
        its built binary is present; skipped otherwise. The CLI runs all three
        starting sets; the python binding takes seed 0 (`lo`), so only the
        `settled lo` lines are compared."""
        candidates = [
            ROOT / "build-wt" / "src" / "coupling" / "apps" / "chflow-coupled-1d3d-multi",
            ROOT / "build" / "src" / "coupling" / "apps" / "chflow-coupled-1d3d-multi",
        ]
        binary = next((p for p in candidates if p.exists()), None)
        if binary is None:
            self.skipTest("coupled multi-edge CLI binary not built")

        r = subprocess.run([str(binary), str(MULTI_CASE)],
                           capture_output=True, text=True,
                           cwd=str(MULTI_CASE.parent))
        self.assertEqual(r.returncode, 0, (r.stdout + r.stderr)[-400:])
        cli_q = {}
        for ln in r.stdout.splitlines():
            m = re.search(r"settled lo edge (\S+) q ([-+0-9.eE]+)", ln)
            if m:
                cli_q[m.group(1)] = float(m.group(2))
        result = chflow.solve_coupled_case(str(MULTI_CASE), 0)
        for ex in result.exchanges:
            self.assertIn(ex.name, cli_q)
            self.assertLess(abs(ex.q - cli_q[ex.name]) / cli_q[ex.name], TOL_Q)


class CoupledCliParitySingleTests(unittest.TestCase):
    def test_single_edge_cli_parity(self):
        """Compare seed-0 flow against the single-edge CLI solver when built."""
        candidates = [
            ROOT / "build-wt" / "src" / "coupling" / "apps" / "chflow-coupled-1d3d",
            ROOT / "build" / "src" / "coupling" / "apps" / "chflow-coupled-1d3d",
        ]
        binary = next((p for p in candidates if p.exists()), None)
        if binary is None:
            self.skipTest("coupled single-edge CLI binary not built")

        r = subprocess.run([str(binary), str(SINGLE_CASE)],
                           capture_output=True, text=True,
                           cwd=str(SINGLE_CASE.parent))
        self.assertEqual(r.returncode, 0, (r.stdout + r.stderr)[-400:])
        m = re.search(r"coupled seed 0: iterations \d+ q ([-+0-9.eE]+)",
                      r.stdout)
        self.assertIsNotNone(m, "no seed-0 coupled line in CLI output")
        cli_q = float(m.group(1))

        result = chflow.solve_coupled_case(str(SINGLE_CASE), 0)
        ex = result.exchanges[0]
        self.assertLess(abs(ex.q - cli_q) / cli_q, TOL_Q)


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
