"""Scalar transport family Python binding tests (CTest
`python.scalar_transport`).

Covers `chflow_core.solve_scalar_transport_case`:

  * advection-diffusion duct (`advectionDiffusion`): solved scalar vs the
    1D steady convection-diffusion analytic profile;
  * pure diffusion (`pureDiffusion`): solved scalar vs the linear z profile;
  * field shape / dtype / copy semantics and the missing-case failure path.

The closed forms and tolerances are the ones from
`tests/check_scalar_profile.py` (advection 5e-2, linearZ 1e-8), not looser
home-grown thresholds.
"""
import math
import sys
import unittest
from pathlib import Path

import numpy as np

import chflow_core

ROOT = Path(__file__).resolve().parents[2]
ADV_CASE = ROOT / "tests" / "cases" / "advectionDiffusion" / "case.dict"
DIFF_CASE = ROOT / "tests" / "cases" / "pureDiffusion" / "case.dict"

# check_scalar_profile.py: analytic_advection(x, U, L, D)
U_ADV, L_ADV, D_ADV = 0.1, 1.0, 0.01
TOL_ADV = 5e-2
TOL_LINEAR_Z = 1e-8


def analytic_advection(x):
    pe = U_ADV * L_ADV / D_ADV
    if abs(pe) < 1e-12:
        return 1.0 - x / L_ADV
    e1 = math.exp(pe * x / L_ADV)
    e2 = math.exp(pe)
    return (e1 - e2) / (1.0 - e2)


class ScalarAdvectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = chflow_core.solve_scalar_transport_case(str(ADV_CASE))

    def test_solved_and_structured(self):
        self.assertTrue(self.result.converged)
        self.assertEqual(self.result.n_cells(), 100)
        T = self.result.T()
        self.assertEqual(T.shape, (100,))
        self.assertEqual(T.dtype, np.float64)
        self.assertTrue(T.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(T).all())
        centres = self.result.cell_centres()
        self.assertEqual(centres.shape, (100, 3))
        self.assertTrue(np.isfinite(centres).all())

    def test_profile_matches_closed_form(self):
        T = self.result.T()
        x = self.result.cell_centres()[:, 0]
        err = max(abs(t - analytic_advection(xi)) for t, xi in zip(T, x))
        self.assertLessEqual(err, TOL_ADV, f"max abs error {err:.4e}")


class ScalarDiffusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = chflow_core.solve_scalar_transport_case(str(DIFF_CASE))

    def test_solved_and_structured(self):
        self.assertTrue(self.result.converged)
        self.assertEqual(self.result.n_cells(), 50)
        T = self.result.T()
        self.assertEqual(T.shape, (50,))
        self.assertEqual(T.dtype, np.float64)
        self.assertTrue(T.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(T).all())
        pts = self.result.points()
        self.assertEqual(pts.ndim, 2)
        self.assertEqual(pts.shape[1], 3)

    def test_linear_profile(self):
        T = self.result.T()
        z = self.result.cell_centres()[:, 2]
        zs = self.result.points()[:, 2]
        zmin, zmax = float(zs.min()), float(zs.max())
        err = max(abs(t - (zi - zmin) / (zmax - zmin))
                  for t, zi in zip(T, z))
        self.assertLessEqual(err, TOL_LINEAR_Z, f"max abs error {err:.4e}")


class ScalarFailurePathTests(unittest.TestCase):
    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_scalar_transport_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
