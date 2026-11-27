"""QMOM family Python binding tests (CTest `python.qmom`).

Covers `chflow_core.solve_qmom_case`: two-node quadrature moment closure for
constant-kernel coagulation. Moments are replayed independently by the published
explicit-Euler Smoluchowski source; the two-node quadrature must reproduce
M0..M3. Closed forms, mass conservation, node growth and refusals follow
check_qmom.py.
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
CASE = ROOT / "tests" / "cases" / "qmomCoagulation" / "case.dict"

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
    d = block(txt, "qmom")
    return {
        "beta": scalar(d, "beta"), "m0": scalar(d, "m0"), "m1": scalar(d, "m1"),
        "m2": scalar(d, "m2"), "m3": scalar(d, "m3"),
        "dt": scalar(d, "dt"), "steps": int(scalar(d, "steps")),
        "reportEvery": int(scalar(d, "reportEvery")),
    }


def source(cfg, m):
    return [-0.5 * cfg["beta"] * m[0] ** 2,
            0.0,
            cfg["beta"] * m[1] ** 2,
            3.0 * cfg["beta"] * m[1] * m[2]]


def replay(cfg):
    m = [cfg["m0"], cfg["m1"], cfg["m2"], cfg["m3"]]
    out = {0: list(m)}
    for step in range(1, cfg["steps"] + 1):
        dm = source(cfg, m)
        m = [m[k] + cfg["dt"] * dm[k] for k in range(4)]
        out[step] = list(m)
    return out


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


class QmomTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASE).read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_qmom_case(str(CASE))
        cls.by_step = {s.step: s for s in cls.res.history}

    def test_records_match_replay(self):
        cfg, hist = self.cfg, replay(self.cfg)
        expected = list(range(0, cfg["steps"] + 1, cfg["reportEvery"]))
        if cfg["steps"] not in expected:
            expected.append(cfg["steps"])
        self.assertEqual(sorted(self.by_step), sorted(expected))
        for step, rec in self.by_step.items():
            for k in range(4):
                self.assertTrue(close(rec.moments[k], hist[step][k], 1e-9, 1e-9),
                                f"step {step} M{k}")

    def test_closed_form_final(self):
        cfg = self.cfg
        t = cfg["steps"] * cfg["dt"]
        last = self.by_step[cfg["steps"]]
        cf = [cfg["m0"] / (1.0 + cfg["beta"] * cfg["m0"] * t / 2.0),
              cfg["m1"],
              cfg["m2"] + cfg["beta"] * cfg["m1"] ** 2 * t,
              cfg["m3"] + 3.0 * cfg["beta"] * cfg["m1"] * cfg["m2"] * t
              + 1.5 * cfg["beta"] ** 2 * cfg["m1"] ** 2 * t ** 2]
        for k in range(4):
            self.assertTrue(close(last.moments[k], cf[k], 2e-2, 1e-6),
                            f"final M{k}")

    def test_quadrature_reproduces_moments(self):
        for rec in self.res.history:
            w, v = rec.weights, rec.nodes
            for k in range(4):
                rebuilt = w[0] * v[0] ** k + w[1] * v[1] ** k
                self.assertTrue(close(rebuilt, rec.moments[k], 1e-8, 1e-8),
                                f"step {rec.step} moment {k}")

    def test_mass_conserved_nodes_grow(self):
        for rec in self.res.history:
            self.assertLess(rec.mass_error, 1e-12)
            self.assertGreater(rec.nodes[0], 0)
            self.assertGreater(rec.nodes[1], 0)
            self.assertGreaterEqual(rec.weights[0], 0)
            self.assertGreaterEqual(rec.weights[1], 0)
        first, last = self.res.history[0], self.res.history[-1]
        self.assertGreater(last.nodes[1], first.nodes[1])
        self.assertLess(last.moments[0], first.moments[0])

    def test_source_and_quadrature_helpers(self):
        cfg = self.cfg
        m = [cfg["m0"], cfg["m1"], cfg["m2"], cfg["m3"]]
        want = source(cfg, m)
        got = self.res.source(m)
        for k in range(4):
            self.assertTrue(close(got[k], want[k], 1e-9, 1e-12))
        quad = self.res.quadrature(m)
        self.assertIsNotNone(quad)
        w, v = quad
        for k in range(4):
            rebuilt = w[0] * v[0] ** k + w[1] * v[1] ** k
            self.assertTrue(close(rebuilt, m[k], 1e-8, 1e-8))


class RefusalTests(unittest.TestCase):
    def _run_text(self, text):
        with tempfile.TemporaryDirectory(prefix="cf_qmom_") as td:
            p = os.path.join(td, "case.dict")
            open(p, "w").write(text)
            with self.assertRaises(RuntimeError) as cm:
                chflow_core.solve_qmom_case(p)
            return str(cm.exception)

    def test_declaration_rejections(self):
        base = open(CASE).read()
        cases = [
            (base.replace("beta 0.05;", "beta -0.05;"),
             "beta must be a finite positive"),
            (base.replace("m1 1.0;", "m1 -1.0;"),
             "M_1 must be finite and positive"),
            (base.replace("m2 1.2;", "m2 0.5;"),
             "violate M_0 M_2 >= M_1^2"),
            (base.replace("steps 20;", "steps 0;"), "steps must be at least 1"),
            (re.sub(r"qmom\s*\{.*?^\}", "", base, count=1, flags=re.M | re.S),
             "needs a qmom"),
        ]
        for text, needle in cases:
            msg = self._run_text(text)
            self.assertIn(needle, msg, f"missing refusal: {needle}")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_qmom_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
