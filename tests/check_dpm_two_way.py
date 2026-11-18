"""D8 increment 3: two-way DPM momentum feedback in a steady 3D duct.

The checker does not trust the trajectory report. It rebuilds the total drag
impulse from the printed end states (u_end - u_0 - body*t, times particle
mass), converts every printed per-cell impulse into a cell acceleration from
the case geometry, and requires the carrier checksum to move while the
trajectory stage itself still leaves U untouched.
"""

import math
import os
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(ROOT, "cases", "dpmTwoWay", "case.dict")
BASE = os.path.join(ROOT, "cases", "dpmChannel", "case.dict")
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VEC = r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)"


def scalar(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def block(text, name):
    m = re.search(r"^\s*" + re.escape(name) + r"\s*\{(.*?)^\}", text, re.M | re.S)
    if not m:
        raise RuntimeError("no " + name + " block")
    return m.group(1)


def key_vec(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+" + VEC, block_text)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def vec3(m, first):
    return tuple(float(m.group(i)) for i in (first, first + 1, first + 2))


def launches_of(output):
    out = {}
    for line in output.splitlines():
        m = re.search(r"dpm launch (\d+)\s*: stream (\d+) x " + VEC
                      + r" u " + VEC + r" diameter (" + NUM + r") rho (" + NUM + r")", line)
        if m:
            out[int(m.group(1))] = {"u": vec3(m, 6), "d": float(m.group(9)),
                                    "rho": float(m.group(10))}
    return out


def ends_of(output):
    out = {}
    for line in output.splitlines():
        m = re.search(r"dpm particle (\d+) end: status (\w+) t (" + NUM + r") x "
                      + VEC + r" u " + VEC, line)
        if m:
            out[int(m.group(1))] = {"status": m.group(2), "t": float(m.group(3)),
                                    "u": vec3(m, 7)}
    return out


def sources_of(output):
    out = {}
    for line in output.splitlines():
        m = re.search(r"dpm source cell (\d+) impulse " + VEC + r" force " + VEC, line)
        if m:
            out[int(m.group(1))] = {"impulse": vec3(m, 2), "force": vec3(m, 5)}
    return out


def diff(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mag(v):
    return math.sqrt(sum(x * x for x in v))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


class TestDpmTwoWay(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if len(sys.argv) < 3:
            raise SystemExit("usage: check_dpm_two_way.py <exe> <case.dict>")
        cls.exe = os.path.abspath(sys.argv[1])
        cls.tmp = tempfile.TemporaryDirectory(prefix="cf_twoway_")
        with open(CASE) as fh:
            cls.text = fh.read()
        with open(BASE) as fh:
            cls.baseText = fh.read()
        dpm = block(cls.text, "dpm")
        cls.tol = scalar(dpm, "twoWayTol")
        cls.relax = scalar(dpm, "twoWayRelax")
        cls.outer = int(scalar(dpm, "outerIterations"))
        cls.rhoF = scalar(dpm, "carrierRho")
        cls.mu = cls.rhoF * scalar(cls.text, "nu")
        cls.gravity = key_vec(dpm, "gravity")
        bb = block(cls.text, "block")
        lo, hi, cells = (key_vec(bb, "min"), key_vec(bb, "max"), key_vec(bb, "cells"))
        cls.vol = ((hi[0] - lo[0]) / cells[0]) * ((hi[1] - lo[1]) / cells[1]) \
            * ((hi[2] - lo[2]) / cells[2])

        cls.completed = subprocess.run([cls.exe, CASE], cwd=cls.tmp.name, text=True,
                                       capture_output=True, timeout=600)
        if cls.completed.returncode != 0:
            raise RuntimeError("two-way case failed:\n" + cls.completed.stdout[-3000:]
                               + cls.completed.stderr[-3000:])
        out = cls.completed.stdout
        cls.launches = launches_of(out)
        cls.ends = ends_of(out)
        cls.sources = sources_of(out)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_picard_converges_in_declared_budget(self):
        rounds = re.findall(r"dpm twoWay outer (\d+) residual (" + NUM + r")", self.completed.stdout)
        self.assertGreaterEqual(len(rounds), 2)
        self.assertLessEqual(len(rounds), self.outer)
        used = re.search(r"dpm twoWay outerUsed (\d+) residual (" + NUM + r")",
                         self.completed.stdout)
        self.assertIsNotNone(used)
        self.assertLessEqual(int(used.group(1)), self.outer)
        self.assertLess(float(used.group(2)), self.tol)
        line = next(l for l in self.completed.stdout.splitlines() if "dpm carrier frozen" in l)
        c = dict((m.group(1), float(m.group(2)))
                 for m in re.finditer(r"\b(uMaxBefore|uMaxAfter|sumUBefore|sumUAfter) ("
                                      + NUM + r")", line))
        self.assertEqual(c["uMaxBefore"], c["uMaxAfter"])
        self.assertEqual(c["sumUBefore"], c["sumUAfter"])

    def test_carrier_actually_moved_and_tracking_is_frozen(self):
        m = re.search(r"dpm twoWay carrier baseline uMax (" + NUM + r") sumU (" + NUM
                      + r") final uMax (" + NUM + r") final sumU (" + NUM + r")",
                      self.completed.stdout)
        self.assertIsNotNone(m)
        u0, s0, u1, s1 = (float(m.group(i)) for i in range(1, 5))
        # not round-off: the one-way baseline is a different printed number
        self.assertGreater(abs(u1 - u0), 1e-9)
        self.assertGreater(abs(s1 - s0), 1e-8)
        line = next(l for l in self.completed.stdout.splitlines() if "dpm carrier frozen" in l)
        c = dict((mm.group(1), float(mm.group(2)))
                 for mm in re.finditer(r"\b(uMaxBefore|uMaxAfter|sumUBefore|sumUAfter) ("
                                       + NUM + r")", line))
        self.assertEqual(c["uMaxBefore"], c["uMaxAfter"])
        self.assertEqual(c["sumUBefore"], c["sumUAfter"])

    def test_total_impulse_is_particle_drag_integral(self):
        # du/dt = aDrag + body  =>  integral(aDrag dt) = (u_end - u_0) - body*t
        want = [0.0, 0.0, 0.0]
        for pid, end in self.ends.items():
            launch = self.launches[pid]
            mp = launch["rho"] * math.pi * launch["d"] ** 3 / 6.0
            body = tuple((1.0 - self.rhoF / launch["rho"]) * g for g in self.gravity)
            du = diff(end["u"], launch["u"])
            for i in range(3):
                want[i] += mp * (du[i] - body[i] * end["t"])
        got = [sum(s["impulse"][i] for s in self.sources.values()) for i in range(3)]
        self.assertGreater(mag(got), 1e-8)
        for i in range(3):
            self.assertLess(abs(got[i] - want[i]), 1e-8 * max(1.0, abs(want[i])))

    def test_cell_force_is_impulse_over_window_and_volume(self):
        window = max(e["t"] for e in self.ends.values())
        self.assertGreater(window, 0)
        for cell, s in self.sources.items():
            self.assertLessEqual(cell, 999)
            for i in range(3):
                expected = -s["impulse"][i] / (self.rhoF * self.vol * window)
                self.assertLess(abs(s["force"][i] - expected),
                                1e-8 * max(1e-12, abs(expected)))

    def test_explicit_no_declaration_is_byte_identical(self):
        twin = re.sub(r"(\n\s*dpm\s*\{)", r"\1\n    twoWay no;", self.baseText, count=1)
        self.assertNotEqual(twin, self.baseText)
        path = os.path.join(self.tmp.name, "explicit_no.dict")
        with open(path, "w") as fh:
            fh.write(twin)
        r = subprocess.run([self.exe, path], cwd=self.tmp.name, text=True,
                           capture_output=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        self.assertNotIn("dpm twoWay", r.stdout)
        ref = subprocess.run([self.exe, BASE], cwd=self.tmp.name, text=True,
                             capture_output=True, timeout=600)
        self.assertEqual(ref.returncode, 0)
        self.assertEqual(r.stdout, ref.stdout)

    def test_declaration_rejections(self):
        cases = [
            (re.sub(r"outerIterations \d+;", "outerIterations 0;", self.text, count=1),
             "outerIterations must be at least 1"),
            (re.sub(r"twoWayRelax [\d.eE+-]+;", "twoWayRelax 1.5;", self.text, count=1),
             "twoWayRelax must be in (0, 1]"),
            (re.sub(r"twoWayTol [\d.eE+-]+;", "twoWayTol 0;", self.text, count=1),
             "twoWayTol must be finite and positive"),
        ]
        for text, needle in cases:
            path = os.path.join(self.tmp.name, "reject.dict")
            with open(path, "w") as fh:
                fh.write(text)
            r = subprocess.run([self.exe, path], cwd=self.tmp.name, text=True,
                               capture_output=True, timeout=600)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn(needle, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
