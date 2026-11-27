"""ALE moving-mesh family Python binding tests (CTest `python.ale`).

Covers `chflow_core.solve_ale_case`: rigid translation and uniform expansion
vertex motion with the Geometric Conservation Law. Anchors follow check_ale.py:

  * translate: closed-cell face area vectors are unchanged, so the mesh-volume
    flux sum(U_m.Sf) is exactly zero: dV/dt = 0 to roundoff; the mesh translates
    rigidly and keeps its total volume;
  * expand: vertices scale by (1+rate dt) per step, so volume grows as the cube
    of the stretch factor and xMax as the factor; the forward-Euler GCL residual
    is O(dt) and halves when dt halves at fixed end time;
  * loud refusals of bad motion, non-positive dt/steps/rate and a missing block.
"""
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

import chflow_core

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "tests" / "cases"
TCASE = CASES / "aleTranslation" / "case.dict"
ECASE = CASES / "aleExpansion" / "case.dict"

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
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


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError(f"no scalar {key}")
    return float(m.group(1))


def vectors(bt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", bt)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def parse(txt):
    ad = block(txt, "ale")
    blk = block(txt, "block")
    mn, mx, nc = vectors(blk, "min"), vectors(blk, "max"), vectors(blk, "cells")
    return {
        "motion": re.search(r"\bmotion\s+(\w+);", ad).group(1),
        "translation": vectors(ad, "translation") if "translation" in ad else (0, 0, 0),
        "rate": scalar(ad, "expansionRate") if "expansionRate" in ad else 0.0,
        "dt": scalar(ad, "dt"), "steps": int(scalar(ad, "steps")),
        "Lx": mx[0] - mn[0], "Ly": mx[1] - mn[1], "Lz": mx[2] - mn[2],
        "xMax0": mx[0],
    }


def solve_text(text):
    with tempfile.TemporaryDirectory(prefix="cf_ale_py_") as td:
        p = os.path.join(td, "case.dict")
        open(p, "w").write(text)
        return chflow_core.solve_ale_case(p)


class TranslateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(TCASE).read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_ale_case(str(TCASE))

    def test_roundoff_gcl_rigid_motion(self):
        cfg, rows, V0 = self.cfg, self.res.history, self.res.initial_volume
        self.assertEqual(len(rows), cfg["steps"])
        vx = cfg["translation"][0]
        for rec in rows:
            self.assertLess(rec.max_gcl, 1e-12)
            self.assertTrue(close(rec.total_volume, V0, 1e-12, 1e-15))
            self.assertTrue(close(rec.x_max, cfg["xMax0"] + vx * rec.t, 1e-9, 1e-9))
            self.assertTrue(close(rec.sum_gcl, 0.0, 0.0, 1e-12))

    def test_final_geometry_shapes(self):
        n = self.res.n_cells()
        self.assertEqual(self.res.volumes.shape, (n,))
        self.assertEqual(self.res.cell_centres.shape, (n, 3))
        self.assertEqual(self.res.points.shape, (self.res.n_points(), 3))


class ExpandTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(ECASE).read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_ale_case(str(ECASE))

    def test_volume_and_length_law(self):
        cfg, rows = self.cfg, self.res.history
        V0 = cfg["Lx"] * cfg["Ly"] * cfg["Lz"]
        self.assertEqual(len(rows), cfg["steps"])
        for rec in rows:
            factor = (1.0 + cfg["rate"] * cfg["dt"]) ** rec.step
            self.assertTrue(close(rec.total_volume, V0 * factor ** 3, 1e-9, 1e-12))
            self.assertTrue(close(rec.x_max, cfg["xMax0"] * factor, 1e-9, 1e-9))
            self.assertLess(rec.max_gcl, 1e-3)

    def test_gcl_error_first_order_in_dt(self):
        def final_gcl(dt, steps):
            txt = self.txt.replace("dt %s;" % self.cfg["dt"], "dt %s;" % dt) \
                          .replace("steps %d;" % self.cfg["steps"], "steps %d;" % steps)
            return solve_text(txt).history[-1].max_gcl

        g1 = final_gcl("0.01", 10)
        g2 = final_gcl("0.005", 20)
        g4 = final_gcl("0.0025", 40)
        self.assertTrue(1.7 < g1 / g2 < 2.3, f"ratio {g1 / g2}")
        self.assertTrue(1.7 < g2 / g4 < 2.3, f"ratio {g2 / g4}")


class RefusalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tb = open(TCASE).read()
        cls.eb = open(ECASE).read()

    def _reject(self, text, needle):
        with self.assertRaises(RuntimeError) as cm:
            solve_text(text)
        self.assertIn(needle, str(cm.exception))

    def test_negative_dt(self):
        self._reject(self.tb.replace("dt 0.1;", "dt -0.1;"),
                     "dt must be finite and positive")

    def test_zero_steps(self):
        self._reject(self.tb.replace("steps 10;", "steps 0;"),
                     "steps must be at least 1")

    def test_bad_motion(self):
        self._reject(re.sub(r"\bmotion\s+\w+\s*;", "motion spin;", self.tb, count=1),
                     "motion must be 'translate' or 'expand'")

    def test_zero_expansion_rate(self):
        self._reject(self.eb.replace("expansionRate 1.0;", "expansionRate 0;"),
                     "expansionRate must be finite and positive")

    def test_missing_block(self):
        txt = re.sub(r"ale\s*\{.*?^\}", "", self.tb, count=1, flags=re.M | re.S)
        self._reject(txt, "needs an ale")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_ale_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
