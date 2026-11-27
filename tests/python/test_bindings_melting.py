"""Phase-change melting family Python binding tests (CTest
`python.melting`).

Covers `chflow_core.solve_melting_case`:

  * latent-step cases a/b/c: one step that enters the mushy band from below,
    goes deeper, and crosses beyond the liquidus. The per-cell enthalpy
    inversion must match the closed form, the liquid fraction must be that
    temperature's own f_l, the sensible-only answer must differ (the
    "one step across the interval must not skip the latent heat" regression),
    and the enthalpy gain must equal the declared source;
  * melting front case: enthalpy balance, front motion, f_l bounds, Darcy
    liquid/solid separation, and the two loud refusals;
  * field shape / dtype / copy semantics and the missing-case failure path.

Closed forms and tolerances are the ones from check_melting_latent_step.py and
check_melting_front.py, not looser home-grown thresholds.
"""
import math
import re
import sys
import unittest
from pathlib import Path

import numpy as np

import chflow_core

ROOT = Path(__file__).resolve().parents[2]
FRONT_CASE = ROOT / "tests" / "cases" / "meltingFront" / "case.dict"
BASE_CASE = ROOT / "tests" / "cases" / "meltingFront" / "caseBase.dict"
REJECT_CASE = ROOT / "tests" / "cases" / "meltingFront" / "caseReject.dict"
LATENT_CASES = {
    "a": ROOT / "tests" / "cases" / "meltingLatentStep" / "a" / "case.dict",
    "b": ROOT / "tests" / "cases" / "meltingLatentStep" / "b" / "case.dict",
    "c": ROOT / "tests" / "cases" / "meltingLatentStep" / "c" / "case.dict",
}

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def clean(txt):
    out = []
    for line in txt.splitlines():
        line = re.sub(r"//.*$", "", line)
        line = re.sub(r"#.*$", "", line)
        out.append(line)
    return "\n".join(out)


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        return None
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


def vec(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", txt)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def f_of(t, ts, tl):
    return min(max((t - ts) / (tl - ts), 0.0), 1.0)


class LatentStepTests(unittest.TestCase):
    """5.196: latent heat settled by enthalpy inversion; a one-shot interval
    crossing must neither skip nor double-count the latent heat."""

    def _check_one(self, tag, path):
        with open(path) as fh:
            text = clean(fh.read())
        cp = scalar(text, "cp")
        rho = scalar(text, "rho")
        L = scalar(block(text, "phaseChange"), "latentHeat")
        ts = scalar(block(text, "phaseChange"), "solidusT")
        tl = scalar(block(text, "phaseChange"), "liquidusT")
        T0 = scalar(text, "initialT")
        q = scalar(text, "volumetricSource")
        dt = scalar(block(text, "piso"), "dt")
        lc = L / cp

        res = chflow_core.solve_melting_case(str(path))
        T = res.T()
        fl = res.liquid_fraction()
        self.assertEqual(T.shape, fl.shape)
        self.assertEqual(T.dtype, np.float64)
        self.assertTrue(T.flags["C_CONTIGUOUS"])
        self.assertEqual(len(T), 1)
        self.assertTrue(np.isfinite(T).all())
        self.assertTrue(np.isfinite(fl).all())

        T1, fl1 = float(T[0]), float(fl[0])
        dT = q / (rho * cp) * dt
        h = T0 + lc * f_of(T0, ts, tl) + dT
        if h <= ts:
            want = h
        elif h >= tl + lc:
            want = h - lc
        else:
            want = (h + lc * ts / (tl - ts)) / (1.0 + lc / (tl - ts))

        # M1: T is the enthalpy inversion of the step (1e-10 rel floor).
        self.assertLessEqual(abs(T1 - want), 1e-10 * max(abs(want), 1.0),
                             f"{tag}: printed {T1:.12g} vs closed form "
                             f"{want:.12g} (h {h:.6g})")
        # M2: f_l is that temperature's own liquid fraction.
        self.assertLessEqual(abs(fl1 - f_of(T1, ts, tl)), 1e-10)
        # M3: the sensible-only answer is a different number (regression gate).
        self.assertGreater(abs((T0 + dT) - T1), 1.0,
                           f"{tag}: sensible {T0 + dT:.6g} vs printed {T1:.6g}")
        # M4: enthalpy gain equals the declared source (latent neither skipped
        # nor double counted).
        def H(t):
            return cp * (t - T0) + L * (f_of(t, ts, tl) - f_of(T0, ts, tl))

        self.assertLessEqual(abs(H(T1) - cp * dT),
                             1e-9 * max(abs(cp * dT), 1.0),
                             f"{tag}: H {H(T1):.12g} vs source {cp * dT:.12g}")

    def test_case_a_enters_band_from_below(self):
        self._check_one("a", LATENT_CASES["a"])

    def test_case_b_deeper_in_band(self):
        self._check_one("b", LATENT_CASES["b"])

    def test_case_c_crosses_beyond_liquidus(self):
        self._check_one("c", LATENT_CASES["c"])


class MeltingFrontTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(FRONT_CASE) as fh:
            cls.text = clean(fh.read())
        cls.result = chflow_core.solve_melting_case(str(FRONT_CASE))

    def test_report_series_and_clock(self):
        rows = self.result.reports
        self.assertGreaterEqual(len(rows), 3)
        self.assertAlmostEqual(rows[0].t, 0.0, delta=1e-12)
        self.assertTrue(all(rows[i].t >= rows[i - 1].t - 1e-12
                            for i in range(1, len(rows))))
        p = block(self.text, "piso")
        endTime = scalar(p, "endTime")
        self.assertAlmostEqual(rows[-1].t, endTime, delta=1e-9)
        self.assertAlmostEqual(rows[0].H, 0.0, delta=1e-9)
        self.assertAlmostEqual(rows[0].Qin, 0.0, delta=1e-9)

    def test_enthalpy_balance(self):
        e = block(self.text, "energy")
        rho, cp = scalar(e, "rho"), scalar(e, "cp")
        ph = block(e, "phaseChange")
        L, ts, tl = (scalar(ph, "latentHeat"), scalar(ph, "solidusT"),
                     scalar(ph, "liquidusT"))
        blk = block(self.text, "block")
        mn, mx = vec(blk, "min"), vec(blk, "max")
        Vtot = (mx[0] - mn[0]) * (mx[1] - mn[1]) * (mx[2] - mn[2])
        Hscale = max(rho * cp * Vtot * (tl - ts), rho * Vtot * L, 1e-30)

        he, qe = self.result.reports[-1].H, self.result.reports[-1].Qin
        self.assertLessEqual(abs(he - qe),
                             0.02 * max(abs(he), abs(qe), 1e-30))
        self.assertLessEqual(abs(he - qe), 0.02 * Hscale)

    def test_front_and_liquid_bounds(self):
        blk = block(self.text, "block")
        mn, mx, cells = vec(blk, "min"), vec(blk, "max"), vec(blk, "cells")
        ylo, yhi = mn[1], mx[1]
        dy = (yhi - ylo) / int(cells[1])
        body = self.result.reports[1:]

        ys = [r.yF for r in body]
        self.assertTrue(all(ylo - 1e-12 <= y <= yhi + 1e-12 for y in ys))
        dys = [ys[i] - ys[i - 1] for i in range(1, len(ys))]
        self.assertTrue(all(d >= -1e-9 for d in dys) or
                        all(d <= 1e-9 for d in dys))
        self.assertGreaterEqual(max(ys) - min(ys), 0.5 * dy)

        for r in body:
            self.assertTrue(0.0 <= r.flMin <= r.flMax <= 1.0)
        self.assertLess(body[-1].flMin, 0.1)
        self.assertGreater(body[-1].flMax, 0.9)

    def test_darcy_liquid_solid_separation(self):
        p = block(self.text, "piso")
        dt = scalar(p, "dt")
        F = vec(p, "force") if re.search(r"\bforce\b", p) else (0.0, 0.0, 0.0)
        Fmag = math.sqrt(F[0] ** 2 + F[1] ** 2 + F[2] ** 2)
        ph = block(block(self.text, "energy"), "phaseChange")
        A, eps = scalar(ph, "mushyConstant"), scalar(ph, "mushyEpsilon")
        damp_at_solid = A * (1.0 - 0.1) ** 2 / (0.1 ** 3 + eps)
        usolid_cap = (6.0 * dt * Fmag + 1e-9
                      if damp_at_solid >= 0.5 / dt else float("inf"))

        pairs = [r for r in self.result.reports[1:]
                 if r.flMin < 0.1 and r.flMax > 0.9]
        self.assertGreaterEqual(len(pairs), 2)
        for r in pairs:
            self.assertGreater(r.uLiquid, 0.0)
            self.assertLessEqual(r.uSolid, 0.1 * r.uLiquid)
            self.assertLessEqual(r.uSolid, usolid_cap)

    def test_final_fields(self):
        T = self.result.T()
        fl = self.result.liquid_fraction()
        U = self.result.U()
        p = self.result.p()
        n = self.result.n_cells()
        self.assertEqual(T.shape, (n,))
        self.assertEqual(fl.shape, (n,))
        self.assertEqual(U.shape, (n, 3))
        self.assertEqual(p.shape, (n,))
        self.assertEqual(T.dtype, np.float64)
        self.assertTrue(np.isfinite(T).all())
        self.assertTrue(np.isfinite(fl).all())
        self.assertTrue(np.isfinite(U).all())

    def test_refusals_raise(self):
        with self.assertRaises(RuntimeError) as cm:
            chflow_core.solve_melting_case(str(BASE_CASE))
        self.assertIn("phaseChange", str(cm.exception))

        with self.assertRaises(RuntimeError) as cm:
            chflow_core.solve_melting_case(str(REJECT_CASE))
        msg = str(cm.exception)
        self.assertTrue(any(k in msg for k in
                            ("liquidusT must exceed solidusT",
                             "latentHeat must be positive",
                             "mushyConstant must be positive",
                             "mushyEpsilon must be positive")))

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_melting_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
