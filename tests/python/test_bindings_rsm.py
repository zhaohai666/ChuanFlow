"""RSM family Python binding tests (CTest `python.rsm`).

Covers `chflow_core.solve_rsm_case`: second-moment Reynolds-stress closure on
the seeded linear shear, advanced by the transient PISO coupling.

Note: unlike the LES family, the repository ships no independent closed-form or
stdout regression checker for the RSM app, so the anchors here are the kernel's
own declared invariants rather than an external closed solution:
  * isotropic start: R11=R22=R33=2 k0/3, off-diagonals 0, k=k0, eps=eps0;
  * realizability held after transport: diagonals >= 0, |Rij| <= sqrt(Rii Rjj),
    k and epsilon above their floors;
  * six stress components are readable as NumPy fields with the expected shapes.
"""
import math
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

import chflow_core

ROOT = Path(__file__).resolve().parents[2]

DICT_TEMPLATE = """// RSM shear case, generated in a temp dir by the Python binding test.
block
{{
    min   (0 0 0);
    max   (1 0.5 0.1);
    cells (4 8 1);

    xmin  {{ name left;   type symmetry; }}
    xmax  {{ name right;  type symmetry; }}
    ymin  {{ name bottom; type wall;     }}
    ymax  {{ name top;    type wall;     }}
    zmin  {{ name front;  type symmetry; }}
    zmax  {{ name back;   type symmetry; }}
}}

nu 0.01;

patches
{{
    left   {{ UType zeroGradient; p 0; pType fixedValue; }}
    right  {{ UType zeroGradient; pType zeroGradient;    }}
    bottom {{ U (0 0 0); UType fixedValue; pType zeroGradient; }}
    top    {{ U (1 0 0); UType fixedValue; pType zeroGradient; }}
    front  {{ UType zeroGradient; pType zeroGradient; }}
    back   {{ UType zeroGradient; pType zeroGradient; }}
}}

schemes
{{
    ddt       Euler;
    div       upwind;
    laplacian corrected;
    grad      linear;
}}

rsm
{{
    C1 {C1}; C1eps 1.44; C2eps 1.92; sigmaR 1.0; sigmaEps 1.3;
    relaxation {relax}; k0 {k0}; eps0 {eps0};
}}

profile
{{
    type linearShear;
    rate 2;
}}

piso
{{
    dt 0.001; endTime 0.005; correctors 2; reportEvery 1;
}}
"""

DEFAULTS = {"C1": 1.8, "relax": 0.7, "k0": 0.01, "eps0": 0.001}


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def solve_text(text):
    with tempfile.TemporaryDirectory(prefix="cf_rsm_py_") as td:
        p = os.path.join(td, "case.dict")
        open(p, "w").write(text)
        return chflow_core.solve_rsm_case(p)


def make_case(**overrides):
    cfg = dict(DEFAULTS)
    cfg.update(overrides)
    return solve_text(DICT_TEMPLATE.format(**cfg))


class RsmAnchorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = make_case()

    def test_initial_isotropic_state(self):
        k0, eps0 = DEFAULTS["k0"], DEFAULTS["eps0"]
        r = self.res.history[0]
        diag = 2.0 * k0 / 3.0
        for i in range(3):
            self.assertTrue(close(r.stress[i], diag, 1e-9, 1e-12),
                            f"initial R{i+1}{i+1}")
        for i in range(3, 6):
            self.assertTrue(close(r.stress[i], 0.0, 1e-9, 1e-12),
                            f"initial off-diagonal {i}")
        # k_ is rebuilt inside correctTurbulence (the first PISO step), so at
        # the app's pre-step report the trace-derived k, not the k_ field, is
        # the declared isotropic value.
        trace_k = 0.5 * (r.stress[0] + r.stress[1] + r.stress[2])
        self.assertTrue(close(trace_k, k0, 1e-9, 1e-12))
        self.assertTrue(close(r.epsilon, eps0, 1e-9, 1e-12))
        self.assertEqual(r.nut, 0.0)

    def test_realizability_after_transport(self):
        n = self.res.n_cells()
        r11, r22, r33 = self.res.stress_11, self.res.stress_22, self.res.stress_33
        r12, r13, r23 = self.res.stress_12, self.res.stress_13, self.res.stress_23
        for i in range(n):
            self.assertGreaterEqual(r11[i], 0.0)
            self.assertGreaterEqual(r22[i], 0.0)
            self.assertGreaterEqual(r33[i], 0.0)
            self.assertLessEqual(abs(r12[i]), math.sqrt(r11[i] * r22[i]) + 1e-12)
            self.assertLessEqual(abs(r13[i]), math.sqrt(r11[i] * r33[i]) + 1e-12)
            self.assertLessEqual(abs(r23[i]), math.sqrt(r22[i] * r33[i]) + 1e-12)
            self.assertGreaterEqual(self.res.k[i], 1e-12)
            self.assertGreaterEqual(self.res.epsilon[i], 1e-12)
            self.assertTrue(math.isfinite(self.res.nut[i]))
        self.assertTrue(all(math.isfinite(r.energy) for r in self.res.history))
        self.assertTrue(close(self.res.history[-1].t, 0.005, 1e-6, 1e-12))

    def test_field_shapes(self):
        n = self.res.n_cells()
        self.assertEqual(self.res.U.shape, (n, 3))
        self.assertEqual(self.res.p.shape, (n,))
        for f in (self.res.stress_11, self.res.stress_22, self.res.stress_33,
                  self.res.stress_12, self.res.stress_13, self.res.stress_23,
                  self.res.epsilon, self.res.k, self.res.nut,
                  self.res.production):
            self.assertEqual(f.shape, (n,))
        self.assertGreaterEqual(len(self.res.history), 2)
        self.assertEqual(self.res.history[0].stress.__len__(), 6)


class RsmRefusalTests(unittest.TestCase):
    def _reject(self, overrides, needle):
        with self.assertRaises(RuntimeError) as cm:
            make_case(**overrides)
        self.assertIn(needle, str(cm.exception))

    def test_negative_c1(self):
        self._reject({"C1": 0.0}, "rsm: C1 must be positive")

    def test_relaxation_out_of_range(self):
        self._reject({"relax": 1.5}, "rsm: relaxation must be in (0, 1]")

    def test_nonpositive_k0(self):
        self._reject({"k0": 0.0}, "rsm: k0 and eps0 must be positive")

    def test_nonpositive_eps0(self):
        self._reject({"eps0": 0.0}, "rsm: k0 and eps0 must be positive")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_rsm_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
