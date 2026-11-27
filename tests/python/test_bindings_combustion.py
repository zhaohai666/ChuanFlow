"""Combustion family Python binding tests (CTest `python.combustion`).

Covers `chflow_core.solve_combustion_case` on the fixed-oxidiser stagnant slab:

  * Arrhenius k(T) independently graded against math.exp at the three scan
    temperatures;
  * fuel mass fraction against the sinh closed form D Y'' - k Yo Y = 0;
  * rate law rate = k Yf Yo; pinned oxidiser stays at fixedOxidiser;
  * discrete fuel balance inflow - outflow = integrated sink;
  * loud refusals (negative diffusivity / preExp / activationE, no block);
  * field shape / dtype / copy semantics.

Closed forms and tolerances are the ones from check_combustion.py, not looser
home-grown thresholds.
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
CASE = ROOT / "tests" / "cases" / "combustionFixedOxidiser" / "case.dict"

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


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


class CombustionFixedOxidiserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = open(CASE).read()
        cb = block(cls.base, "combustion")
        cls.D = scalar(cb, "diffusivity")
        cls.A = scalar(cb, "preExp")
        cls.Ea = scalar(cb, "activationE")
        cls.R = scalar(cb, "gasConstant")
        cls.T = scalar(cb, "temperature")
        cls.Y0 = scalar(cb, "fuelInlet")
        cls.Yo = scalar(cb, "fixedOxidiser")
        blk = block(cls.base, "block")
        cls.mn = vectors(blk, "min")
        cls.mx = vectors(blk, "max")
        cls.nc = vectors(blk, "cells")
        cls.Lx = cls.mx[0] - cls.mn[0]
        cls.Nx = int(cls.nc[0])
        cls.dx = cls.Lx / cls.Nx
        cls.k = cls.A * math.exp(-cls.Ea / (cls.R * cls.T))
        cls.result = chflow_core.solve_combustion_case(str(CASE))

    def test_arrhenius_independent(self):
        for fac in (0.8, 1.0, 1.25):
            Tp = self.T * fac
            kp = self.result.arrhenius_rate(Tp)
            kexp = self.A * math.exp(-self.Ea / (self.R * Tp))
            self.assertTrue(close(kp, kexp, 1e-9),
                            f"k(T={Tp})={kp} != {kexp}")

    def test_fields_structured(self):
        Yf, Yo, rate = self.result.Yf, self.result.Yo, self.result.reaction_rate
        self.assertEqual(Yf.shape, (self.Nx,))
        self.assertEqual(Yo.shape, (self.Nx,))
        self.assertEqual(rate.shape, (self.Nx,))
        self.assertEqual(Yf.dtype, np.float64)
        self.assertTrue(Yf.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(Yf).all())
        self.assertTrue(np.isfinite(Yo).all())
        self.assertTrue(np.isfinite(rate).all())
        self.assertTrue(self.result.solved)
        self.assertTrue(self.result.fixed_oxidiser)

    def test_analytic_fuel_profile(self):
        p = math.sqrt(self.k * self.Yo / self.D)
        tol = 3e-3 * self.Y0
        centres = self.result.cell_centres()
        for i in range(self.Nx):
            xc = self.mn[0] + (i + 0.5) * self.dx
            ya = self.Y0 * math.sinh(p * (self.Lx - xc)) / math.sinh(p * self.Lx)
            self.assertLessEqual(abs(self.result.Yf[i] - ya), tol,
                                 f"cell {i}: Yf {self.result.Yf[i]} vs "
                                 f"analytic {ya}")
            self.assertAlmostEqual(centres[i, 0], xc, delta=1e-12)

    def test_rate_law_and_pinned_oxidiser(self):
        for i in range(self.Nx):
            yf, yo, rate = self.result.Yf[i], self.result.Yo[i], self.result.reaction_rate[i]
            self.assertTrue(close(rate, self.k * yf * yo, 1e-9, 1e-12),
                            f"cell {i}: rate {rate} != k Yf Yo")
            self.assertLess(abs(yo - self.Yo), 1e-12)

    def test_fuel_balance(self):
        resid = self.result.fuel_in - self.result.fuel_out - self.result.fuel_consumed
        p = math.sqrt(self.k * self.Yo / self.D)
        analytic_sink = sum(
            self.k * self.Y0 * math.sinh(
                p * (self.Lx - (self.mn[0] + (i + 0.5) * self.dx)))
            / math.sinh(p * self.Lx) * self.dx
            for i in range(self.Nx)) * self.Yo
        self.assertLessEqual(abs(resid), 1e-6 * max(analytic_sink, 1e-3),
                             f"balance residual {resid} (sink ~{analytic_sink})")


class CombustionRefusalTests(unittest.TestCase):
    def _run_text(self, text):
        with tempfile.TemporaryDirectory(prefix="cf_comb_") as td:
            p = os.path.join(td, "case.dict")
            open(p, "w").write(text)
            with self.assertRaises(RuntimeError) as cm:
                chflow_core.solve_combustion_case(p)
            return str(cm.exception)

    def test_refusals_loud(self):
        base = open(CASE).read()
        cases = [
            (base.replace("diffusivity 1.0;", "diffusivity -1.0;"),
             "diffusivity must be finite and positive"),
            (base.replace("preExp 100.0;", "preExp -100.0;"),
             "preExp A must be finite and positive"),
            (base.replace("activationE 50000.0;", "activationE -1.0;"),
             "activationE Ea cannot be negative"),
            (re.sub(r"combustion\s*\{.*?^\}", "", base, count=1,
                    flags=re.M | re.S),
             "needs a combustion"),
        ]
        for text, needle in cases:
            msg = self._run_text(text)
            self.assertIn(needle, msg, f"missing refusal: {needle}")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_combustion_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
