"""LES SGS family Python binding tests (CTest `python.les`).

Covers `chflow_core.solve_les_case`: algebraic Smagorinsky / WALE sub-grid
viscosity on the seeded linear shear, advanced by the transient PISO coupling.
Anchors are rebuilt analytically from the case text (mirroring check_les_sgs.py):
uniform block mesh filter width Delta = V^(1/3), unidirectional shear |S|=rate,
Smagorinsky nut = (Cs Delta)^2 |S|, WALE nut ~ 0 on pure shear, and the analytic
volume-averaged kinetic energy of the seeded profile.
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
CASES = ROOT / "tests" / "cases" / "lesSgs"

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


def scalar(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        if default is not None:
            return default
        raise RuntimeError(f"no scalar {key}")
    return float(m.group(1))


def word(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_]\w*)\s*;", txt)
    if not m:
        if default is not None:
            return default
        raise RuntimeError(f"no word {key}")
    return m.group(1)


def vector(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", txt)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def parse(txt):
    blk = block(txt, "block")
    lo, hi, nc = vector(blk, "min"), vector(blk, "max"), vector(blk, "cells")
    Nx, Ny, Nz = int(nc[0]), int(nc[1]), int(nc[2])
    Lx, Ly, Lz = hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]
    V = Lx * Ly * Lz / (Nx * Ny * Nz)
    les = block(txt, "les") if re.search(r"\bles\s*\{", txt) else ""
    prof = block(txt, "profile") if re.search(r"\bprofile\s*\{", txt) else ""
    piso = block(txt, "piso")
    return {
        "Nx": Nx, "Ny": Ny, "Nz": Nz, "Ly": Ly, "V": V, "delta": V ** (1 / 3),
        "model": word(les, "model", "smagorinsky"),
        "Cs": scalar(les, "Cs", 0.17), "Cw": scalar(les, "Cw", 0.325),
        "rate": scalar(prof, "rate", 1.0),
        "endTime": scalar(piso, "endTime"),
    }


def solve_text(text):
    with tempfile.TemporaryDirectory(prefix="cf_les_py_") as td:
        p = os.path.join(td, "case.dict")
        open(p, "w").write(text)
        return chflow_core.solve_les_case(p)


class SmagorinskyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASES / "case.dict").read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_les_case(str(CASES / "case.dict"))

    def test_t0_reports_match_analytic(self):
        c, r = self.cfg, self.res.history[0]
        self.assertTrue(close(r.t, 0.0, 0.0, 1e-12))
        self.assertTrue(close(r.strain_max, abs(c["rate"]), 1e-6, 1e-12))
        self.assertGreater(r.strain_max, 0)
        nut_true = c["Cs"] ** 2 * c["delta"] ** 2 * abs(c["rate"])
        self.assertTrue(close(r.nut_max, nut_true, 1e-6, 1e-15))
        self.assertTrue(close(r.nut_min, nut_true, 1e-6, 1e-15))
        self.assertTrue(close(r.nut_mean, nut_true, 1e-6, 1e-15))
        self.assertGreater(r.nut_max, 0)
        self.assertTrue(close(self.res.filter_width[0], c["delta"], 1e-6, 1e-15))

    def test_energy_analytic_and_decay(self):
        c, rows = self.cfg, self.res.history
        Ny = c["Ny"]
        E0 = 0.5 * c["rate"] ** 2 * c["Ly"] ** 2 * (4 - 1 / Ny ** 2) / 12
        self.assertTrue(close(rows[0].energy, E0, 1e-9, 1e-15))
        self.assertTrue(all(math.isfinite(r.energy) for r in rows))
        self.assertTrue(all(r.energy >= 0.0 for r in rows))
        self.assertLessEqual(rows[-1].energy, rows[0].energy * (1 + 1e-6) + 1e-15)
        self.assertTrue(close(rows[-1].t, c["endTime"], 1e-6, 1e-12))
        self.assertGreaterEqual(len(rows), 2)
        self.assertTrue(all(rows[i + 1].t >= rows[i].t - 1e-12
                            for i in range(len(rows) - 1)))

    def test_field_shapes(self):
        n = self.res.n_cells()
        self.assertEqual(self.res.U.shape, (n, 3))
        self.assertEqual(self.res.p.shape, (n,))
        self.assertEqual(self.res.nut.shape, (n,))
        self.assertEqual(self.res.strain_rate.shape, (n,))
        self.assertEqual(self.res.filter_width.shape, (n,))
        self.assertEqual(self.res.model, "les-smagorinsky")


class WaleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(CASES / "caseWale.dict").read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_les_case(str(CASES / "caseWale.dict"))

    def test_wale_nut_zero_on_pure_shear(self):
        c, rows = self.cfg, self.res.history
        nut_smag = c["Cs"] ** 2 * c["delta"] ** 2 * abs(c["rate"])
        tol = 1e-6 * max(nut_smag, 1e-12)
        self.assertTrue(close(rows[0].strain_max, abs(c["rate"]), 1e-6, 1e-12))
        for r in rows:
            self.assertGreaterEqual(r.nut_max, 0.0)
            self.assertLessEqual(r.nut_max, tol)
            self.assertLessEqual(r.nut_min, tol)

    def test_smagorinsky_twin_is_nonzero(self):
        c = self.cfg
        twin = re.sub(r"\bmodel\s+wale\s*;", "model smagorinsky;", self.txt, count=1)
        res = solve_text(twin)
        nut_smag = c["Cs"] ** 2 * c["delta"] ** 2 * abs(c["rate"])
        self.assertGreater(res.history[0].nut_max, 0)
        self.assertTrue(close(res.history[0].nut_max, nut_smag, 1e-6, 1e-15))


class ZeroCsTests(unittest.TestCase):
    def test_cs_zero_nut_off_strain_kept(self):
        res = chflow_core.solve_les_case(str(CASES / "caseZero.dict"))
        r = res.history[0]
        self.assertEqual(r.nut_max, 0.0)
        self.assertEqual(r.nut_min, 0.0)
        self.assertEqual(r.nut_mean, 0.0)
        self.assertGreater(r.strain_max, 0)


class RefusalTests(unittest.TestCase):
    def test_unknown_model(self):
        base = open(CASES / "case.dict").read()
        bad = re.sub(r"\bmodel\s+\w+\s*;", "model bogus;", base, count=1)
        with self.assertRaises(RuntimeError) as cm:
            solve_text(bad)
        self.assertIn("les: unknown model 'bogus'", str(cm.exception))

    def test_negative_cs(self):
        base = open(CASES / "case.dict").read()
        bad = re.sub(r"\b(Cs\s+)" + NUM + r"\s*;", r"\g<1>-0.1;", base, count=1)
        with self.assertRaises(RuntimeError) as cm:
            solve_text(bad)
        self.assertIn("les: Cs must be finite and non-negative", str(cm.exception))

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_les_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
