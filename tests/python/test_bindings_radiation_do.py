"""Grey radiation DO / S2S / P1 family binding tests (CTest `python.radiation_do`).

Covers `chflow_core.solve_radiation_case`. Anchors are rebuilt analytically from
the case constants, mirroring check_radiation_do.py:

  * DO: per-ordinate exponential attenuation on a non-scattering black-wall slab;
    cell-averaged G, cell-centre q, net source S, closed-form wall fluxes and the
    discrete RTE balance, at 1e-6 relative;
  * S2S: transparent two-plate limit q = sigma(TL^4 - TR^4),
    G = 2 sigma(TL^4 + TR^4), uniform, zero volumetric source, F_LR = 1;
  * P1: qualitative agreement trend with the DO transport solution;
  * loud refusals for unknown model / coefficients / quadrature / missing block.
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
CASES = ROOT / "tests" / "cases" / "radiationDo"
SIGMA = 5.670374419e-8
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

HALF = {
    4: [(0.8611363115940526, 0.3478548451374538),
        (0.3399810435848563, 0.6521451548625461)],
    8: [(0.9602898564975363, 0.1012285362903763),
        (0.7966664774136267, 0.2223810344533745),
        (0.5255324099163290, 0.3137066458778873),
        (0.1834346424956498, 0.3626837833783620)],
}


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
    rad = block(txt, "radiation")
    blk = block(txt, "block")
    patches = block(txt, "patches")
    mn, mx, nc = vectors(blk, "min"), vectors(blk, "max"), vectors(blk, "cells")
    qm = re.search(r"\bquadrature\s+[Ss]?(\d+)", txt)
    return {
        "ka": scalar(rad, "absorption"), "ks": scalar(rad, "scattering"),
        "Tm": scalar(rad, "mediumTemperature"),
        "Lx": mx[0] - mn[0], "Nx": int(nc[0]), "mn": mn,
        "TL": scalar(block(patches, "left"), "wallTemperature"),
        "TR": scalar(block(patches, "right"), "wallTemperature"),
        "order": int(qm.group(1)) if qm else 8,
    }


def analytic(cfg, order):
    N, Lx, dx = cfg["Nx"], cfg["Lx"], cfg["Lx"] / cfg["Nx"]
    beta = cfg["ka"] + cfg["ks"]
    pi = math.pi
    Ib = SIGMA * cfg["Tm"] ** 4 / pi
    IL = SIGMA * cfg["TL"] ** 4 / pi
    IR = SIGMA * cfg["TR"] ** 4 / pi
    half = HALF[order]
    ords = [(mu, w) for mu, w in half] + [(-mu, w) for mu, w in half]

    def face_I(mu, x):
        if mu > 0:
            return Ib + (IL - Ib) * math.exp(-beta * x / mu)
        return Ib + (IR - Ib) * math.exp(-beta * (Lx - x) / abs(mu))

    def avg_I(mu, xc):
        a = abs(mu)
        i_in = face_I(mu, xc - dx / 2) if mu > 0 else face_I(mu, xc + dx / 2)
        return Ib + (i_in - Ib) * (a / (beta * dx)) * (1 - math.exp(-beta * dx / a))

    G, q, S = {}, {}, {}
    for i in range(N):
        xc = cfg["mn"][0] + (i + 0.5) * dx
        G[i] = 2 * pi * sum(w * avg_I(mu, xc) for mu, w in ords)
        q[i] = 2 * pi * sum(w * mu * face_I(mu, xc) for mu, w in ords)
        S[i] = cfg["ka"] * (G[i] - 4 * pi * Ib)
    q0 = 2 * pi * sum(w * mu * face_I(mu, 0.0) for mu, w in ords)
    qL = 2 * pi * sum(w * mu * face_I(mu, Lx) for mu, w in ords)
    vol = dx * sum(S.values())
    return G, q, S, q0, qL, vol, q0 - qL - vol


def solve_text(text):
    with tempfile.TemporaryDirectory(prefix="cf_rad_py_") as td:
        p = os.path.join(td, "case.dict")
        open(p, "w").write(text)
        return chflow_core.solve_radiation_case(p)


class DoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASES / "case.dict").read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_radiation_case(str(CASES / "case.dict"))

    def test_fields_match_closed_form(self):
        G, q, S, q0, qL, vol, bal = analytic(self.cfg, self.cfg["order"])
        n = self.cfg["Nx"]
        self.assertEqual(self.res.G.shape, (n,))
        self.assertEqual(self.res.q.shape, (n, 3))
        self.assertEqual(self.res.S.shape, (n,))
        for i in range(n):
            self.assertTrue(close(self.res.G[i], G[i], 1e-6, 1e-6), f"G {i}")
            self.assertTrue(close(self.res.S[i], S[i], 1e-6, 1e-6), f"S {i}")
            self.assertTrue(close(self.res.q[i, 0], q[i], 1e-6, 1e-6), f"qx {i}")
            self.assertLess(abs(self.res.q[i, 1]), 1e-6)
            self.assertLess(abs(self.res.q[i, 2]), 1e-6)
        self.assertTrue(close(self.res.q_left, q0, 1e-6, 1e-6))
        self.assertTrue(close(self.res.q_right, -qL, 1e-6, 1e-6))
        self.assertTrue(close(self.res.vol_source, vol, 1e-6, 1e-6))
        gmean = sum(G.values()) / n
        self.assertTrue(close(self.res.G.mean(), gmean, 1e-6, 1e-6))
        scale = max(abs(q0), abs(qL), abs(vol), 1.0)
        self.assertLess(abs(self.res.balance) / scale, 1e-8)
        self.assertTrue(close(self.res.balance, bal, 1e-6, 1e-9))

    def test_physics_signs_and_ordering(self):
        self.assertEqual(self.res.model, "do")
        self.assertGreater(self.res.q_left, 0)
        self.assertLess(self.res.q_right, 0)
        g = self.res.G
        self.assertLess(g.min(), g.mean())
        self.assertLess(g.mean(), g.max())
        self.assertTrue(all(g[i] > g[i + 1] for i in range(len(g) - 1)))

    def test_s4_quadrature(self):
        txt = self.txt.replace("quadrature S8;", "quadrature S4;")
        res = solve_text(txt)
        cfg = dict(self.cfg, order=4)
        G, q, S, q0, qL, vol, _ = analytic(cfg, 4)
        n = cfg["Nx"]
        for i in range(n):
            self.assertTrue(close(res.G[i], G[i], 1e-6, 1e-6))
            self.assertTrue(close(res.S[i], S[i], 1e-6, 1e-6))
            self.assertTrue(close(res.q[i, 0], q[i], 1e-6, 1e-6))
        self.assertTrue(close(res.q_left, q0, 1e-6, 1e-6))
        self.assertTrue(close(res.q_right, -qL, 1e-6, 1e-6))
        self.assertTrue(close(res.vol_source, vol, 1e-6, 1e-6))


class S2sTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASES / "case.dict").read()
        cls.cfg = parse(cls.txt)
        cls.res = solve_text(cls.txt.replace("model do;", "model s2s;"))

    def test_transparent_plate_limit(self):
        c = self.cfg
        q = SIGMA * (c["TL"] ** 4 - c["TR"] ** 4)
        G = 2 * SIGMA * (c["TL"] ** 4 + c["TR"] ** 4)
        self.assertEqual(self.res.model, "s2s")
        self.assertTrue(close(self.res.q_left, q, 1e-6, 1e-6))
        self.assertTrue(close(self.res.q_right, -q, 1e-6, 1e-6))
        self.assertTrue(close(self.res.G.mean(), G, 1e-6, 1e-6))
        self.assertTrue(close(self.res.G.min(), self.res.G.max(), 1e-9, 1e-9))
        self.assertLess(abs(self.res.vol_source), 1e-6 * max(abs(q), 1.0))
        self.assertLess(abs(self.res.balance) / max(abs(q), 1.0), 1e-9)
        self.assertTrue(close(self.res.view_factor, 1.0, 1e-12, 1e-12))
        for i in range(c["Nx"]):
            self.assertTrue(close(self.res.q[i, 0], q, 1e-6, 1e-6))
            self.assertTrue(close(self.res.G[i], G, 1e-6, 1e-6))
            self.assertLess(abs(self.res.S[i]), 1e-6)


class P1TrendTests(unittest.TestCase):
    def test_p1_trend_agrees_with_do(self):
        do = chflow_core.solve_radiation_case(str(CASES / "case.dict"))
        p1 = chflow_core.solve_radiation_case(str(CASES / "caseP1.dict"))
        self.assertEqual(p1.model, "p1")
        self.assertTrue(close(p1.G.mean(), do.G.mean(), 0.25))
        self.assertGreater(p1.q_left, 0)
        self.assertLess(p1.q_right, 0)
        self.assertLess(p1.G.min(), p1.G.max())


class RefusalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = open(CASES / "case.dict").read()

    def _reject(self, text, needle):
        with self.assertRaises(RuntimeError) as cm:
            solve_text(text)
        self.assertIn(needle, str(cm.exception))

    def test_unknown_model(self):
        self._reject(re.sub(r"\bmodel\s+\w+\s*;", "model bogus;", self.base, count=1),
                     "unknown radiation model")

    def test_negative_absorption(self):
        self._reject(self.base.replace("absorption 1.0;", "absorption -1.0;"),
                     "absorption kappa_a must be a finite positive")

    def test_negative_scattering(self):
        self._reject(self.base.replace("scattering 0;", "scattering -0.5;"),
                     "scattering kappa_s cannot be negative")

    def test_zero_medium_temperature(self):
        self._reject(self.base.replace("mediumTemperature 1200;", "mediumTemperature 0;"),
                     "mediumTemperature must be a finite positive")

    def test_zero_wall_temperature(self):
        self._reject(self.base.replace("wallTemperature 1500;", "wallTemperature 0;"),
                     "wallTemperature must be finite and positive")

    def test_bad_quadrature(self):
        self._reject(self.base.replace("quadrature S8;", "quadrature S16;"),
                     "unknown quadrature")

    def test_missing_radiation_block(self):
        txt = re.sub(r"radiation\s*\{.*?^\}", "", self.base, count=1, flags=re.M | re.S)
        self._reject(txt, "needs a radiation")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_radiation_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
