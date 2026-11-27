"""DPM family Python binding tests (CTest `python.dpm`).

Covers `chflow_core.solve_dpm_case` on the two verified ducts:

  * dpmChannel -- one-way: steady carrier solved once, particles tracked through
    the frozen field. Schiller-Naumann Cd/Re, drag acceleration, explicit-Euler
    velocity/position updates, launch state and terminal status are reconciled
    against the closed forms in check_dpm_channel.py.
  * dpmTwoWay -- Picard carrier <-> trajectory sweep; outer count, stopping
    residual, per-cell drag impulse (particle drag integral) are checked.

Closed forms and tolerances follow the checkers; no looser home-grown thresholds.
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
CHANNEL = ROOT / "tests" / "cases" / "dpmChannel" / "case.dict"
TWOWAY = ROOT / "tests" / "cases" / "dpmTwoWay" / "case.dict"

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VEC = r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)"


def scalar(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    if not m:
        raise RuntimeError(f"no scalar {key}")
    return float(m.group(1))


def block(text, name):
    m = re.search(r"^\s*" + re.escape(name) + r"\s*\{(.*?)^\s*\}", text, re.M | re.S)
    if not m:
        raise RuntimeError(f"no {name} block")
    return m.group(1)


def key_vec(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+" + VEC, block_text)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def mag(v):
    return math.sqrt(sum(x * x for x in v))


def cd_of(reynolds):
    if reynolds <= 1000:
        return (24.0 / reynolds) * (1.0 + 0.15 * reynolds ** 0.687)
    return 0.44


def diff(a, b):
    return tuple(x - y for x, y in zip(a, b))


class ChannelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = open(CHANNEL).read()
        dpm = block(cls.text, "dpm")
        cls.dt = scalar(dpm, "dt")
        cls.rhoF = scalar(dpm, "carrierRho")
        cls.mu = cls.rhoF * scalar(cls.text, "nu")
        cls.gravity = key_vec(dpm, "gravity")
        cls.traceSteps = int(scalar(dpm, "traceSteps"))
        bb = block(cls.text, "block")
        cls.xMax = key_vec(bb, "max")[0]
        cls.yMin = key_vec(bb, "min")[1]
        cls.result = chflow_core.solve_dpm_case(str(CHANNEL))

    def test_carrier_and_fields_structured(self):
        U, p = self.result.U, self.result.p
        n = self.result.n_cells()
        self.assertEqual(U.shape, (n, 3))
        self.assertEqual(p.shape, (n,))
        self.assertEqual(U.dtype, np.float64)
        self.assertTrue(U.flags["C_CONTIGUOUS"])
        self.assertTrue(np.isfinite(U).all())
        self.assertTrue(np.isfinite(p).all())
        self.assertTrue(self.result.carrier_converged)
        self.assertFalse(self.result.two_way)

    def test_launch_state_and_terminal_status(self):
        parts = self.result.particles
        self.assertEqual(len(parts), 2)
        self.assertEqual(parts[0].diameter, 1e-3)
        self.assertEqual(parts[0].rho_particle, 2500.0)
        self.assertEqual(parts[1].diameter, 3e-3)
        self.assertEqual(parts[1].rho_particle, 8000.0)
        self.assertEqual(parts[0].status, "outlet")
        self.assertEqual(parts[1].status, "wall")
        self.assertLess(abs(parts[0].x[0] - (self.xMax - 0.5e-3)), 1e-5)
        self.assertLess(abs(parts[1].x[1] - (self.yMin + 1.5e-3)), 1e-4)
        self.assertEqual(self.result.outlet_count, 1)
        self.assertEqual(self.result.wall_count, 1)
        self.assertEqual(self.result.expired_count, 0)

    def test_schiller_naumann_and_reynolds(self):
        by_id = {}
        for s in self.result.samples:
            by_id.setdefault(s.id, []).append(s)
        samples = by_id[0]
        rhoP, d = 2500.0, 1e-3
        for s in samples:
            rel = diff(s.u_fluid, s.u)
            wantRe = self.rhoF * mag(rel) * d / self.mu
            self.assertLess(abs(s.reynolds - wantRe) / wantRe, 1e-9)
            self.assertLess(abs(s.cd - cd_of(s.reynolds)) / s.cd, 1e-9)
            rate = 3.0 * self.rhoF * s.cd / (4.0 * rhoP * d)
            wantA = tuple(rate * r for r in rel)
            self.assertLess(mag(diff(s.drag_acceleration, wantA)), 1e-7)

    def test_explicit_euler_updates(self):
        by_id = {}
        for s in self.result.samples:
            by_id.setdefault(s.id, []).append(s)
        seq = [s for s in sorted(by_id[0], key=lambda s: s.step)
               if s.step < self.traceSteps]
        body = tuple((1.0 - self.rhoF / 2500.0) * g for g in self.gravity)
        for a, b in zip(seq, seq[1:]):
            want_u = tuple(a.u[i] + self.dt * (a.drag_acceleration[i] + body[i])
                            for i in range(3))
            self.assertLess(mag(diff(b.u, want_u)), 1e-10)
            want_x = tuple(a.x[i] + self.dt * b.u[i] for i in range(3))
            self.assertLess(mag(diff(b.x, want_x)), 1e-12)

    def test_particle_approaches_carrier(self):
        by_id = {}
        for s in self.result.samples:
            by_id.setdefault(s.id, []).append(s)
        first = sorted(by_id[0], key=lambda s: s.step)[0]
        last = self.result.particles[0]
        self.assertLessEqual(first.u[0], 1e-12)
        self.assertGreater(last.u[0], 0.05)
        self.assertLess(abs(last.u[0] - 0.12), 0.05)

    def test_standalone_drag_law(self):
        for re in (0.1, 1.0, 50.0, 500.0, 2000.0):
            self.assertLess(abs(chflow_core.dpm_sphere_cd(re) - cd_of(re)) / cd_of(re),
                            1e-12)
        uf, up = (0.1, 0.0, 0.0), (0.0, 0.0, 0.0)
        cd = cd_of(self.rhoF * 0.1 * 1e-3 / self.mu)
        rate = 3.0 * self.rhoF * cd / (4.0 * 2500.0 * 1e-3)
        want = tuple(rate * r for r in diff(uf, up))
        got = chflow_core.dpm_drag_acceleration(uf, up, self.rhoF, 2500.0, 1e-3, self.mu)
        self.assertLess(mag(diff(got, want)), 1e-10)


class TwoWayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = open(TWOWAY).read()
        dpm = block(cls.text, "dpm")
        cls.tol = scalar(dpm, "twoWayTol")
        cls.outer = int(scalar(dpm, "outerIterations"))
        cls.rhoF = scalar(dpm, "carrierRho")
        cls.mu = cls.rhoF * scalar(cls.text, "nu")
        cls.gravity = key_vec(dpm, "gravity")
        cls.result = chflow_core.solve_dpm_case(str(TWOWAY))

    def test_picard_converged_in_budget(self):
        self.assertTrue(self.result.two_way)
        self.assertGreaterEqual(self.result.outer_used, 1)
        self.assertLessEqual(self.result.outer_used, self.outer)
        self.assertLess(self.result.two_way_residual, self.tol)
        self.assertTrue(self.result.carrier_converged)
        self.assertEqual(self.result.expired_count, 0)

    def test_total_impulse_is_drag_integral(self):
        want = [0.0, 0.0, 0.0]
        for p in self.result.particles:
            mp = p.rho_particle * math.pi * p.diameter ** 3 / 6.0
            body = tuple((1.0 - self.rhoF / p.rho_particle) * g for g in self.gravity)
            du = diff(p.u, p.u0)
            for i in range(3):
                want[i] += mp * (du[i] - body[i] * p.age)
        got = self.result.cell_drag_impulse.sum(axis=0)
        self.assertGreater(mag(got), 1e-8)
        for i in range(3):
            self.assertLess(abs(got[i] - want[i]), 1e-8 * max(1.0, abs(want[i])))


class RefusalTests(unittest.TestCase):
    def _run_text(self, text):
        with tempfile.TemporaryDirectory(prefix="cf_dpm_") as td:
            p = os.path.join(td, "case.dict")
            open(p, "w").write(text)
            with self.assertRaises(RuntimeError) as cm:
                chflow_core.solve_dpm_case(p)
            return str(cm.exception)

    def test_declaration_rejections(self):
        base = open(CHANNEL).read()
        cases = [
            (base.replace("drag schillerNaumann", "drag unknownLaw"),
             "drag law 'unknownLaw' is not implemented"),
            (re.sub(r"(stream0\s*\{.*?diameter\s+)-?[\d.eE+-]+", r"\g<1>-1",
                    base, count=1, flags=re.S),
             "diameter must be finite and positive"),
            (re.sub(r"\ndpm\s*\{.*?\n\}\n", "\n", base, count=1, flags=re.S),
             "needs a dpm{} declaration"),
        ]
        for text, needle in cases:
            msg = self._run_text(text)
            self.assertIn(needle, msg, f"missing refusal: {needle}")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_dpm_case(
                str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
