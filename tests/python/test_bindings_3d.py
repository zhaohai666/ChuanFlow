"""3D-family Python binding tests (CTest `python.solver3d`).

Minimum real chain: read the block mesh -> build fields -> run the built-in
steady P1 radiation solve -> read the fields back as NumPy arrays.

  * the P1 field is checked against the same double-exponential closed form the
    CLI checker `tests/check_p1_slab.py` rebuilds from the case constants, with
    that checker's own tolerances;
  * the Field <-> ndarray conversion contract is pinned: dtype float64,
    C-order copies, scalar (n,) / vector (n, 3), strict shape and dtype checks
    on write-back;
  * loud refusals stay loud.
"""
import math
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

import chflow

ROOT = Path(__file__).resolve().parents[2]
P1_CASE = ROOT / "tests" / "cases" / "radiationP1Slab" / "case.dict"

SIGMA = 5.670374419e-8
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def _scalar(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    assert m, f"no scalar {key}"
    return float(m.group(1))


def _vector(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", text)
    assert m, f"no vector {key}"
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def analytic_profile(case_text):
    """Rebuild the P1 closed form, mirroring check_p1_slab.py exactly."""
    ka = _scalar(case_text, "absorption")
    ks = _scalar(case_text, "scattering")
    tm = _scalar(case_text, "mediumTemperature")
    mn = _vector(case_text, "min")
    mx = _vector(case_text, "max")
    nc = _vector(case_text, "cells")
    tl = _scalar(re.search(r"left\s*\{[^}]*\}", case_text).group(0),
                 "wallTemperature")
    tr = _scalar(re.search(r"right\s*\{[^}]*\}", case_text).group(0),
                 "wallTemperature")

    lx = mx[0] - mn[0]
    nx = int(nc[0])
    dx = lx / nx
    beta = ka + ks
    gamma = 1.0 / (3.0 * beta)
    mm = math.sqrt(ka / gamma)

    gm = 4 * SIGMA * tm ** 4
    g0 = 4 * SIGMA * tl ** 4
    gl = 4 * SIGMA * tr ** 4
    c1 = g0 - gm
    c2 = (gl - gm - c1 * math.cosh(mm * lx)) / math.sinh(mm * lx)

    def g_at(x):
        return gm + c1 * math.cosh(mm * x) + c2 * math.sinh(mm * x)

    return {"ka": ka, "gm": gm, "mn": mn, "dx": dx, "nx": nx,
            "jump": abs(g0 - gl), "g_at": g_at}


class MeshAndFieldContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh = chflow.load_mesh(str(P1_CASE))
        cls.n = cls.mesh.n_cells()

    def test_mesh_geometry_arrays(self):
        centres = self.mesh.cell_centres()
        volumes = self.mesh.cell_volumes()
        self.assertEqual(centres.shape, (self.n, 3))
        self.assertEqual(volumes.shape, (self.n,))
        self.assertEqual(centres.dtype, np.float64)
        self.assertEqual(volumes.dtype, np.float64)
        self.assertTrue(centres.flags["C_CONTIGUOUS"])
        self.assertEqual(len(self.mesh.patch_names()), 6)
        self.assertTrue(np.isfinite(centres).all())
        self.assertTrue((volumes > 0).all())

    def test_scalar_field_round_trip_and_copy(self):
        f = chflow.ScalarField3d(self.mesh, "s", 1.0)
        self.assertEqual(f.n_cells(), self.n)
        arr = f.to_numpy()
        self.assertEqual(arr.shape, (self.n,))
        self.assertEqual(arr.dtype, np.float64)
        self.assertTrue(arr.flags["C_CONTIGUOUS"])
        np.testing.assert_allclose(arr, 1.0)

        # to_numpy() is a copy: mutating the array must not touch the field.
        arr[0] = 999.0
        self.assertEqual(f.to_numpy()[0], 1.0)

        values = np.full(self.n, 2.5, dtype=np.float64)
        f.from_numpy(values)
        np.testing.assert_allclose(f.to_numpy(), 2.5)

    def test_scalar_field_rejects_bad_input(self):
        f = chflow.ScalarField3d(self.mesh, "s")
        with self.assertRaises(ValueError):  # size mismatch
            f.from_numpy(np.full(self.n - 1, 1.0))
        with self.assertRaises(ValueError):  # wrong rank
            f.from_numpy(np.zeros((self.n, 3)))
        with self.assertRaises(TypeError):  # dtype must be float64
            f.from_numpy(np.full(self.n, 1.0, dtype=np.float32))

    def test_vector_field_round_trip(self):
        f = chflow.VectorField3d(self.mesh, "v")
        arr = f.to_numpy()
        self.assertEqual(arr.shape, (self.n, 3))
        self.assertEqual(arr.dtype, np.float64)
        self.assertTrue(arr.flags["C_CONTIGUOUS"])

        values = np.tile(np.array([1.0, 2.0, 3.0]), (self.n, 1))
        f.from_numpy(values)
        back = f.to_numpy()
        np.testing.assert_allclose(back[:, 0], 1.0)
        np.testing.assert_allclose(back[:, 1], 2.0)
        np.testing.assert_allclose(back[:, 2], 3.0)

    def test_vector_field_rejects_bad_input(self):
        f = chflow.VectorField3d(self.mesh, "v")
        with self.assertRaises(ValueError):  # wrong rank
            f.from_numpy(np.zeros(self.n))
        with self.assertRaises(ValueError):  # wrong last dimension
            f.from_numpy(np.zeros((self.n, 2)))
        with self.assertRaises(TypeError):  # dtype must be float64
            f.from_numpy(np.zeros((self.n, 3), dtype=np.float32))


class P1SolveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case_text = P1_CASE.read_text()
        cls.ana = analytic_profile(cls.case_text)
        cls.result = chflow.solve_p1_case(str(P1_CASE))

    def test_solve_converged_and_balanced(self):
        self.assertTrue(self.result.solved)
        self.assertGreaterEqual(self.result.iterations, 1)
        scale = max(abs(self.result.volume_source),
                    abs(self.result.boundary_flux), 1.0)
        self.assertLess(abs(self.result.balance) / scale, 1e-6)
        self.assertLess(self.result.volume_source * self.result.boundary_flux, 0)

    def test_field_shapes_and_dtypes(self):
        g = self.result.G()
        q = self.result.q_rad()
        s = self.result.net_source()
        self.assertEqual(g.shape, (self.result.n_cells(),))
        self.assertEqual(q.shape, (self.result.n_cells(), 3))
        self.assertEqual(s.shape, (self.result.n_cells(),))
        for arr in (g, q, s):
            self.assertEqual(arr.dtype, np.float64)
            self.assertTrue(arr.flags["C_CONTIGUOUS"])
            self.assertTrue(np.isfinite(arr).all())

    def test_profile_matches_closed_form(self):
        centres = self.result.mesh.cell_centres()
        g = self.result.G()
        order = np.argsort(centres[:, 0])
        tol_abs = 3e-3 * self.ana["jump"]  # checker's O(dx^2) tolerance
        for rank, idx in enumerate(order):
            xc = self.ana["mn"][0] + (rank + 0.5) * self.ana["dx"]
            self.assertLessEqual(abs(g[idx] - self.ana["g_at"](xc)), tol_abs,
                                 f"cell {idx} at x={xc}")

    def test_net_source_law(self):
        g = self.result.G()
        s = self.result.net_source()
        expected = self.ana["ka"] * (g - self.ana["gm"])
        np.testing.assert_allclose(s, expected, rtol=1e-9, atol=1.0)

    def test_monotonic_and_transverse_flux(self):
        centres = self.result.mesh.cell_centres()
        g = self.result.G()
        q = self.result.q_rad()
        order = np.argsort(centres[:, 0])
        gseq = g[order]
        self.assertTrue(np.all(gseq[:-1] > gseq[1:]))
        self.assertLess(np.max(np.abs(q[:, 1])), 1e-4)
        self.assertLess(np.max(np.abs(q[:, 2])), 1e-4)

    def test_missing_radiation_refused(self):
        bad = re.sub(r"radiation\s*\{.*?^\}", "", self.case_text,
                     count=1, flags=re.M | re.S)
        with tempfile.TemporaryDirectory(prefix="cf_p1_py_") as td:
            p = Path(td) / "bad.dict"
            p.write_text(bad)
            with self.assertRaises(RuntimeError):
                chflow.solve_p1_case(str(p))

    def test_cli_binary_parity_when_available(self):
        candidates = [
            ROOT / "build-wt" / "src" / "3d" / "apps" / "chflow-3d-p1",
            ROOT / "build" / "src" / "3d" / "apps" / "chflow-3d-p1",
        ]
        binary = next((p for p in candidates if p.exists()), None)
        if binary is None:
            self.skipTest("P1 CLI binary not built")

        r = subprocess.run([str(binary), str(P1_CASE)],
                           capture_output=True, text=True,
                           cwd=str(P1_CASE.parent))
        self.assertEqual(r.returncode, 0, (r.stdout + r.stderr)[-400:])
        cli = {}
        for ln in r.stdout.splitlines():
            m = re.search(r"p1 cell (\d+): G (" + NUM + r")", ln)
            if m:
                cli[int(m.group(1))] = float(m.group(2))
        g = self.result.G()
        tol_abs = 3e-3 * self.ana["jump"]
        for i, gv in cli.items():
            self.assertLessEqual(abs(g[i] - gv), tol_abs, f"cell {i}")


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
