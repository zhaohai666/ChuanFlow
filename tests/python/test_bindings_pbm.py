"""Mixture PBM family Python binding tests (CTest `python.pbm`).

Covers `chflow_core.solve_mixture_pbm_case`: the monodisperse moment balance
with drift-flux slip. The trajectory is replayed independently by the same
explicit-Euler rule from the case constants; closed forms, physicality, the
slip=0 redistribution and the no-growth twin follow check_mixture_pbm.py.
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
CASE = ROOT / "tests" / "cases" / "mixturePbm" / "case.dict"
BASE = ROOT / "tests" / "cases" / "mixturePbm" / "caseBase.dict"

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError(f"no scalar {key}")
    return float(m.group(1))


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


def parse(txt):
    d = block(txt, "mixturePbm")
    return {
        "alpha0": scalar(d, "alpha0"), "d0": scalar(d, "d0"),
        "rhoL": scalar(d, "rhoL"), "rhoG": scalar(d, "rhoG"),
        "slip": scalar(d, "slip"), "j": scalar(d, "j"),
        "J": scalar(d, "J"), "breakage": scalar(d, "breakage"),
        "nFrag": scalar(d, "nFrag"), "dt": scalar(d, "dt"),
        "steps": int(scalar(d, "steps")),
        "reportEvery": int(scalar(d, "reportEvery")),
    }


def derive(cfg):
    v = math.pi * cfg["d0"] ** 3 / 6.0
    k = cfg["j"] / (cfg["j"] + cfg["slip"])
    m0 = cfg["alpha0"] / v
    return {"v": v, "k": k, "m0": [m0, cfg["alpha0"], m0 * v * v, m0 * v * v * v]}


def source(cfg, d, m):
    v, n = d["v"], cfg["nFrag"]
    f2 = 1.0 - 1.0 / n
    f3 = 1.0 - 1.0 / (n * n)
    b = cfg["breakage"]
    return [
        cfg["J"] + (n - 1.0) * b * m[0],
        cfg["J"] * v,
        cfg["J"] * v * v - f2 * b * v * v * m[0],
        cfg["J"] * v ** 3 - f3 * b * v ** 3 * m[0],
    ]


def replay(cfg, d):
    m = list(d["m0"])
    out = {0: list(m)}
    for step in range(1, cfg["steps"] + 1):
        dm = source(cfg, d, m)
        m = [m[i] + cfg["dt"] * dm[i] for i in range(4)]
        out[step] = list(m)
    return out


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


class PbmMainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASE).read()
        cls.cfg = parse(cls.txt)
        cls.d = derive(cls.cfg)
        cls.res = chflow_core.solve_mixture_pbm_case(str(CASE))
        cls.by_step = {s.step: s for s in cls.res.history}

    def test_records_match_replay(self):
        cfg, d, hist = self.cfg, self.d, replay(self.cfg, self.d)
        expected = list(range(0, cfg["steps"] + 1, cfg["reportEvery"]))
        if cfg["steps"] not in expected:
            expected.append(cfg["steps"])
        self.assertEqual(sorted(self.by_step), sorted(expected))
        for step, rec in self.by_step.items():
            self.assertLessEqual(abs(rec.t - step * cfg["dt"]), 1e-12)
            for i in range(4):
                self.assertTrue(close(rec.moments[i], hist[step][i], 1e-8, 1e-30),
                                f"step {step} M{i}")
            alpha = hist[step][1] * d["k"]
            self.assertTrue(close(rec.alpha, alpha, 1e-8, 1e-30))
            uG = cfg["j"] + (1.0 - alpha) * cfg["slip"]
            uL = cfg["j"] - alpha * cfg["slip"]
            self.assertTrue(close(rec.u_gas, uG, 1e-8, 1e-12))
            self.assertTrue(close(rec.u_liquid, uL, 1e-8, 1e-12))

    def test_closed_forms(self):
        cfg, d = self.cfg, self.d
        t = cfg["steps"] * cfg["dt"]
        last = self.by_step[cfg["steps"]]
        m1_cf = cfg["alpha0"] + cfg["J"] * d["v"] * t
        self.assertTrue(close(last.moments[1], m1_cf, 1e-8, 1e-15))
        A = (cfg["nFrag"] - 1.0) * cfg["breakage"]
        if A > 0:
            m0_cf = (d["m0"][0] + cfg["J"] / A) * math.exp(A * t) - cfg["J"] / A
            self.assertTrue(close(last.moments[0], m0_cf, 2e-2, 1e-12))
        self.assertTrue(close(last.u_gas - last.u_liquid, cfg["slip"], 1e-8, 1e-12))

    def test_physicality_and_growth(self):
        for rec in self.res.history:
            for mi in rec.moments:
                self.assertTrue(math.isfinite(mi) and mi > 0)
            self.assertTrue(0.0 < rec.alpha < 1.0)
        first, last = self.res.history[0], self.res.history[-1]
        self.assertGreater(last.moments[0], first.moments[0])
        self.assertGreater(last.moments[1], first.moments[1])
        self.assertGreater(last.moments[2], 0)
        self.assertGreater(last.moments[3], 0)

    def test_source_law(self):
        cfg, d = self.cfg, self.d
        m = d["m0"]
        want = source(cfg, d, m)
        got = self.res.source(m)
        for i in range(4):
            self.assertTrue(close(got[i], want[i], 1e-9, 1e-12))

    def test_closure_constants(self):
        self.assertTrue(close(self.res.v_part, self.d["v"], 1e-12))
        self.assertTrue(close(self.res.holdup_factor, self.d["k"], 1e-12))
        self.assertTrue(close(self.res.m0_initial, self.d["m0"][0], 1e-12))
        self.assertTrue(close(self.res.m1_initial, self.d["m0"][1], 1e-12))


class SlipAndTwinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASE).read()
        cls.cfg = parse(cls.txt)
        cls.main = chflow_core.solve_mixture_pbm_case(str(CASE))
        cls.base = chflow_core.solve_mixture_pbm_case(str(BASE))

    def test_slip_zero_redistributes_holdup(self):
        with tempfile.TemporaryDirectory(prefix="cf_pbm_") as td:
            p = os.path.join(td, "slip0.dict")
            open(p, "w").write(self.txt.replace("slip 0.2;", "slip 0.0;"))
            twin = chflow_core.solve_mixture_pbm_case(p)
        s0, m0 = twin.history[0], self.main.history[0]
        self.assertTrue(close(s0.alpha, self.cfg["alpha0"], 1e-8, 1e-12))
        self.assertGreater(s0.alpha, m0.alpha)
        self.assertTrue(close(s0.u_gas - s0.u_liquid, 0.0, 1e-8, 1e-12))

    def test_no_growth_twin_constant(self):
        cfg = parse(open(BASE).read())
        d = derive(cfg)
        for rec in self.base.history:
            for i in range(4):
                self.assertTrue(close(rec.moments[i], d["m0"][i], 1e-8, 1e-30))
            self.assertTrue(close(rec.alpha, d["m0"][1] * d["k"], 1e-8, 1e-30))


class RefusalTests(unittest.TestCase):
    def _run_text(self, text):
        with tempfile.TemporaryDirectory(prefix="cf_pbm_") as td:
            p = os.path.join(td, "case.dict")
            open(p, "w").write(text)
            with self.assertRaises(RuntimeError) as cm:
                chflow_core.solve_mixture_pbm_case(p)
            return str(cm.exception)

    def test_declaration_rejections(self):
        base = open(CASE).read()
        cases = [
            (base.replace("J 1000000;", "J -1000000;"),
             "J (nucleation rate) must be finite and non-negative"),
            (base.replace("breakage 0.5;", "breakage -0.5;"),
             "breakage frequency must be finite and non-negative"),
            (base.replace("alpha0 0.05;", "alpha0 1.5;"), "alpha0 too large"),
            (base.replace("rhoL 1000;", "rhoL -1000;"),
             "rhoL must be finite and positive"),
            (re.sub(r"mixturePbm\s*\{.*?^\}", "", base, count=1, flags=re.M | re.S),
             "needs a mixturePbm"),
        ]
        for text, needle in cases:
            msg = self._run_text(text)
            self.assertIn(needle, msg, f"missing refusal: {needle}")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_mixture_pbm_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
