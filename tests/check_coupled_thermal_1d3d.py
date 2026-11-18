"""D8 increment 1: 1D<->3D conjugate-heat interface.

The 3D block owns an EnergyModel and reports its wall heat / enthalpy account;
the 1D side carries one lumped thermal node per coupled edge and hands its
temperature back as the 3D inlet condition.  Every printed number is read from
the app's full-precision log; the checker rebuilds the 1D implicit-Euler
recurrence and the 3D heat balance independently from the case dict.
"""

import math
import os
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.abspath(__file__))
THERMAL_CASE = os.path.join(ROOT, "cases", "coupledThermal1d3d", "case.dict")
PARENT_CASE = os.path.join(ROOT, "cases", "transientCoupled1d3d", "case.dict")
TOKEN = re.compile(r"([A-Za-z][A-Za-z0-9]*)\s+([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")


def block(text, name):
    m = re.search(r"^\s*" + re.escape(name) + r"\s*\{([^{}]*)\}", text, re.M | re.S)
    if not m:
        raise RuntimeError("no " + name + " block in case")
    return m.group(1)


def scalar(text, key):
    m = re.search(r"\b" + re.escape(key)
                  + r"\s+([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*;", text)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def word(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_]\w*)\s*;", text)
    if not m:
        raise RuntimeError("no word " + key)
    return m.group(1)


def set_key(text, key, value):
    pat = re.compile(r"(\b" + re.escape(key) + r"\s+)[-+]?[\d.eE+]+(\s*;)")
    out, n = pat.subn(lambda m: m.group(1) + str(value) + m.group(2), text, 1)
    if n != 1:
        raise RuntimeError("set_key: marker not found: " + key)
    return out


def drop_block(text, name):
    pat = re.compile(r"\n\s*" + re.escape(name) + r"\s*\{[^{}]*\}", re.S)
    out, n = pat.subn("\n", text, 1)
    if n != 1:
        raise RuntimeError("drop_block: marker not found: " + name)
    return out


def run(exe, case_path, cwd):
    return subprocess.run([exe, case_path], cwd=cwd, text=True,
                          capture_output=True, timeout=300)


def steps_of(output):
    steps = []
    for line in output.splitlines():
        if "transient step " not in line or ":" not in line:
            continue
        body = line.split(":", 1)[1]
        steps.append({m.group(1): float(m.group(2)) for m in TOKEN.finditer(body)})
    return steps


class TestCoupledThermal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if len(sys.argv) < 3:
            raise SystemExit("usage: check_coupled_thermal_1d3d.py <exe> <case.dict>")
        cls.exe = os.path.abspath(sys.argv[1])
        cls.tmp = tempfile.TemporaryDirectory(prefix="cf_thermal_")
        with open(THERMAL_CASE) as fh:
            cls.case_text = fh.read()
        cls.energy = block(cls.case_text, "energy")
        cls.thermal = block(cls.case_text, "thermal")
        cls.rho = scalar(cls.energy, "rho")
        cls.cp = scalar(cls.energy, "cp")
        cls.capacity = scalar(cls.thermal, "capacity1d")
        cls.relax = scalar(cls.thermal, "relax")
        cls.sourceT = scalar(cls.thermal, "sourceTemperature")
        cls.length = scalar(cls.thermal, "length")
        cls.tol = scalar(cls.thermal, "tolerance")
        cls.dt = scalar(block(cls.case_text, "time"), "dt")
        cls.wall = word(cls.thermal, "wallPatch")

        cls.live = run(cls.exe, THERMAL_CASE, cls.tmp.name)
        if cls.live.returncode != 0:
            raise RuntimeError("thermal case failed:\n" + cls.live.stdout[-3000:]
                               + cls.live.stderr[-3000:])
        cls.steps = steps_of(cls.live.stdout)

        twin = drop_block(drop_block(cls.case_text, "energy"), "thermal")
        twin_path = os.path.join(cls.tmp.name, "twin.dict")
        with open(twin_path, "w") as fh:
            fh.write(twin)
        cls.twin = run(cls.exe, twin_path, cls.tmp.name)

        with open(PARENT_CASE) as fh:
            parent = fh.read()
        parent = re.sub(r"(\bendTime\s+)[-+]?[\d.eE+]+(\s*;)",
                        lambda m: m.group(1) + str(cls.dt * len(cls.steps)) + m.group(2),
                        parent, count=1)
        ref_path = os.path.join(cls.tmp.name, "reference.dict")
        with open(ref_path, "w") as fh:
            fh.write(parent)
        cls.reference = run(cls.exe, ref_path, cls.tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_thermal_fields_on_every_step(self):
        expected = int(round(0.2 / self.dt))
        self.assertEqual(len(self.steps), expected)
        for k, s in enumerate(self.steps, 1):
            self.assertEqual(s["thermal"], 1)
            self.assertAlmostEqual(s["t"], k * self.dt, delta=1e-12)
            self.assertGreater(s["qWall"], 0)
            self.assertGreater(s["T1"], self.sourceT)

    def test_temperature_handoff_and_line_heat(self):
        for k, s in enumerate(self.steps):
            self.assertEqual(s["Tin"], s["T1"])
            self.assertLess(abs(s["qPrime"] - s["qWall"] / self.length)
                            / max(abs(s["qPrime"]), 1e-300), 1e-12)
            self.assertLessEqual(s["tres"], self.tol)
            if k:
                self.assertGreater(s["T1"], self.steps[k - 1]["T1"])

    def test_3d_heat_balance_is_not_a_printed_claim(self):
        for s in self.steps:
            recomputed = s["qWall"] + s["qIn"] - s["hOut"]
            self.assertLess(abs(recomputed - s["hbal"]), 1e-6)
            self.assertLess(abs(recomputed), 1e-5)

    def test_1d_lumped_implicit_euler_recurrence(self):
        prevT = self.sourceT
        rhoCp = self.rho * self.cp
        worst = 0.0
        for s in self.steps:
            predicted = prevT + self.dt / self.capacity * (
                rhoCp * s["q"] * (self.sourceT - prevT) + s["qWall"])
            worst = max(worst, abs(s["T1"] - predicted))
            prevT = s["T1"]
        # The accepted pass still carries the case's residual (K); the printed
        # node is the iterate handed to the 3D inlet, not the post-update value.
        self.assertLess(worst, 2.0 * self.tol)

    def test_zero_declaration_twin_is_byte_identical(self):
        self.assertEqual(self.twin.returncode, 0, self.twin.stdout[-2000:])
        self.assertEqual(self.reference.returncode, 0, self.reference.stdout[-2000:])
        self.assertEqual(self.twin.stdout, self.reference.stdout)
        self.assertFalse(any("thermal 1" in l for l in self.twin.stdout.splitlines()))
        self.assertTrue(any("thermal 1" in l for l in self.live.stdout.splitlines()))

    def test_hydraulics_are_not_moved_by_the_energy_pass(self):
        ref = steps_of(self.reference.stdout)
        self.assertEqual(len(ref), len(self.steps))
        for a, b in zip(self.steps, ref):
            self.assertLess(abs(a["q"] - b["q"]) / b["q"], 1e-5)
            self.assertLess(abs(a["dpBlock"] - b["dpBlock"]) / b["dpBlock"], 1e-5)

    def test_declaration_rejections(self):
        cases = [
            (drop_block(self.case_text, "thermal"), "energy{} is declared without thermal{}"),
            (re.sub(r"\bwallPatch\s+bottom\s*;", "", self.case_text),
             "thermal{} needs a wallPatch"),
            (set_key(self.case_text, "capacity1d", 0), "capacity1d must be finite and positive"),
        ]
        for text, needle in cases:
            path = os.path.join(self.tmp.name, "reject.dict")
            with open(path, "w") as fh:
                fh.write(text)
            r = run(self.exe, path, self.tmp.name)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn(needle, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
