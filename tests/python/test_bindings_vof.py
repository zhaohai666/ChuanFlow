"""VOF family Python binding tests (CTest `python.vof`).

Covers `chflow_core.solve_vof_case` on the two closed-form transport cases:

  * vofTranslation -- a plane interface carried by a uniform velocity: phase
    volume conserves against U*A*t, alpha stays in [0,1], the crossing moves at
    U; the plane reads kappa ~ 0.
  * vofPlic -- a cell-sharp disk carried by solid-body rotation, reconstructed
    with PLIC: the band matches alpha's own band, planes are unit planes lying
    on their offset, interface area spans the perimeter, kappaPlic reads -1/R
    on the band and the normal stays radial.

Closed forms and tolerances follow check_vof.py; the momentum / surface-tension
branch is out of this family's scope and refuses loudly.
"""
import math
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

import chflow_core

ROOT = Path(__file__).resolve().parents[2]
TRANS = ROOT / "tests" / "cases" / "vofTranslation" / "case.dict"
PLIC = ROOT / "tests" / "cases" / "vofPlic" / "case.dict"

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError(f"no scalar {key}")
    return float(m.group(1))


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError(f"no {name} block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced " + name)


def vectors(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", block_text)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


class TranslationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(TRANS).read()
        blk = block(cls.txt, "block")
        cls.mn = vectors(blk, "min")
        cls.mx = vectors(blk, "max")
        cls.nc = vectors(blk, "cells")
        vblk = block(cls.txt, "velocity")
        cls.U = vectors(vblk, "U")
        ablk = block(cls.txt, "alpha")
        cls.plane = vectors(ablk, "point")
        cls.endTime = scalar(cls.txt, "endTime")
        cls.dx = (cls.mx[0] - cls.mn[0]) / int(cls.nc[0])
        cls.area = (cls.mx[1] - cls.mn[1]) * (cls.mx[2] - cls.mn[2])
        cls.result = chflow_core.solve_vof_case(str(TRANS))

    def test_fields_structured(self):
        a, n, k = self.result.alpha, self.result.normal, self.result.kappa
        n_cells = self.result.n_cells()
        self.assertEqual(a.shape, (n_cells,))
        self.assertEqual(n.shape, (n_cells, 3))
        self.assertEqual(k.shape, (n_cells,))
        self.assertEqual(a.dtype, np.float64)
        self.assertTrue(a.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(a).all())
        self.assertTrue(np.isfinite(n).all())
        self.assertTrue(np.isfinite(k).all())
        self.assertFalse(self.result.plic_enabled)

    def test_boundedness(self):
        a = self.result.alpha
        self.assertGreaterEqual(a.min(), -1e-12)
        self.assertLessEqual(a.max(), 1.0 + 1e-12)

    def test_conservation(self):
        rows = self.result.reports
        self.assertTrue(rows)
        last = rows[-1]
        fed = self.U[0] * self.area * last.t
        a_out = self.result.alpha[-int(self.nc[1] * self.nc[2]):].max()
        deficit = self.result.mass0 + fed - last.mass
        self.assertGreaterEqual(deficit, -1e-12)
        self.assertLessEqual(deficit, fed * a_out + 1e-12)
        self.assertLessEqual(deficit / self.result.mass0, 1e-5)

    def test_interface_speed(self):
        rows = sorted(self.result.reports, key=lambda r: r.t)
        first, last = rows[0], rows[-1]
        speed = (last.x50 - first.x50) / (last.t - first.t)
        self.assertLessEqual(abs(speed / self.U[0] - 1.0), 0.01)
        bias = last.x50 - (self.plane[0] + self.U[0] * last.t)
        self.assertLessEqual(abs(bias), 0.4 * self.dx)

    def test_plane_curvature(self):
        band = np.linalg.norm(self.result.normal, axis=1) > 1e-12
        if band.any():
            self.assertLessEqual(np.abs(self.result.kappa[band]).max(), 1e-6)


class PlicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(PLIC).read()
        blk = block(cls.txt, "block")
        cls.mn = vectors(blk, "min")
        cls.mx = vectors(blk, "max")
        cls.nc = vectors(blk, "cells")
        ablk = block(cls.txt, "alpha")
        cblk = block(ablk, "circle")
        cls.center = vectors(cblk, "center")
        cls.R = scalar(cblk, "radius")
        cls.dz = cls.mx[2] - cls.mn[2]
        cls.tol = 1e-6
        cls.result = chflow_core.solve_vof_case(str(PLIC))

    def test_fields_structured(self):
        r = self.result
        n_cells = r.n_cells()
        self.assertEqual(r.alpha.shape, (n_cells,))
        self.assertEqual(r.alpha_initial.shape, (n_cells,))
        self.assertEqual(r.plic_normal.shape, (n_cells, 3))
        self.assertEqual(r.plic_offset.shape, (n_cells,))
        self.assertEqual(r.plic_area.shape, (n_cells,))
        self.assertEqual(r.plic_centroid.shape, (n_cells, 3))
        self.assertEqual(r.plic_fraction.shape, (n_cells,))
        self.assertEqual(r.kappa_plic.shape, (n_cells,))
        self.assertEqual(r.alpha.dtype, np.float64)
        self.assertTrue(r.plic_enabled)
        self.assertTrue(np.isfinite(r.alpha).all())

    def test_exact_filling(self):
        exact = math.pi * self.R * self.R * self.dz
        self.assertLessEqual(abs(self.result.mass0 / exact - 1.0), 1e-8)

    def test_band_matches_alpha(self):
        a = self.result.alpha_initial
        alpha_band = {i for i, v in enumerate(a)
                      if self.tol < v < 1.0 - self.tol}
        plic_band = {i for i, v in enumerate(self.result.plic_normal)
                     if math.hypot(v[0], v[1], v[2]) > 0.5}
        self.assertEqual(plic_band, alpha_band)
        self.assertEqual(self.result.band_cells, len(alpha_band))

    def test_planes_are_unit_planes(self):
        band = [i for i, v in enumerate(self.result.plic_normal)
                if math.hypot(v[0], v[1], v[2]) > 0.5]
        nrm, off, area, ctr = (self.result.plic_normal, self.result.plic_offset,
                               self.result.plic_area, self.result.plic_centroid)
        unit = max(abs(math.hypot(*nrm[i]) - 1.0) for i in band)
        planar = max(abs(nrm[i][0] * ctr[i][0] + nrm[i][1] * ctr[i][1]
                         + nrm[i][2] * ctr[i][2] - off[i]) for i in band)
        self.assertLessEqual(unit, 5e-10)
        self.assertLessEqual(planar, 5e-10)
        self.assertGreater(min(area[i] for i in band), 0.0)
        outside = [i for i in range(len(area)) if i not in set(band)]
        self.assertEqual(max(abs(area[i]) for i in outside), 0.0)

    def test_area_span_and_curvature(self):
        band = [i for i, v in enumerate(self.result.plic_normal)
                if math.hypot(v[0], v[1], v[2]) > 0.5]
        perimeter = 2.0 * math.pi * self.R * self.dz
        a_sum = self.result.plic_area[band].sum()
        self.assertTrue(0.85 < a_sum / perimeter < 1.05)
        kmean = self.result.kappa_plic[band].mean()
        self.assertLess(kmean, 0.0)
        self.assertLessEqual(abs(kmean * self.R + 1.0), 0.1)

    def test_normals_radial(self):
        nrm, ctr = self.result.plic_normal, self.result.plic_centroid
        angles = []
        for i, v in enumerate(nrm):
            if math.hypot(v[0], v[1], v[2]) <= 0.5:
                continue
            rx = ctr[i][0] - self.center[0]
            ry = ctr[i][1] - self.center[1]
            r = math.hypot(rx, ry)
            if r < 1e-9:
                continue
            d = -(v[0] * rx + v[1] * ry) / r
            angles.append(math.degrees(math.acos(max(-1.0, min(1.0, d)))))
        self.assertTrue(angles)
        self.assertLessEqual(sum(angles) / len(angles), 12.0)
        self.assertLessEqual(max(angles), 20.0)

    def test_fraction_readback(self):
        band = [i for i, v in enumerate(self.result.plic_normal)
                if math.hypot(v[0], v[1], v[2]) > 0.5]
        err = max(abs(self.result.plic_fraction[i] - self.result.alpha_initial[i])
                  for i in band)
        self.assertLessEqual(err, 1e-6)


class RefusalTests(unittest.TestCase):
    def test_momentum_mode_refused(self):
        # vofDroplet declares surfaceTension{} only; the momentum branch is
        # out of this family's scope, so the binding refuses loudly instead.
        with self.assertRaises(RuntimeError) as cm:
            chflow_core.solve_vof_case(
                str(ROOT / "tests" / "cases" / "vofDroplet" / "case.dict"))
        self.assertIn("momentum", str(cm.exception))

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_vof_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
