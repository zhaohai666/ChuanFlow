"""D8 increment 2: one-way Lagrangian particles in a solved 3D duct.

The checker never trusts the app's trajectory: it parses the frozen carrier
checksum, the launch state, and the first explicit-Euler samples, then rebuilds
the Schiller-Naumann drag coefficient, particle Reynolds number, drag
acceleration, buoyancy and the velocity/position updates from the case dict.
"""

import math
import os
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(ROOT, "cases", "dpmChannel", "case.dict")
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VEC = r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)"


def scalar(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def block(text, name):
    m = re.search(r"^\s*" + re.escape(name) + r"\s*\{(.*?)^\s*\}", text, re.M | re.S)
    if not m:
        raise RuntimeError("no " + name + " block")
    return m.group(1)


def key_vec(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+" + VEC, block_text)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def launch_of(output):
    launches = {}
    for line in output.splitlines():
        m = re.search(r"dpm launch (\d+)\s*: stream (\d+) x " + VEC
                      + r" u " + VEC + r" diameter (" + NUM + r") rho (" + NUM + r")",
                      line)
        if m:
            g = [float(m.group(i)) for i in range(1, 11)]
            launches[int(g[0])] = {
                "stream": int(g[1]), "x": tuple(g[2:5]), "u": tuple(g[5:8]),
                "d": g[8], "rho": g[9],
            }
    return launches


def samples_of(output):
    samples = {}
    for line in output.splitlines():
        m = re.search(r"dpm sample (\d+) (\d+)\s*: t (" + NUM + r") x " + VEC
                      + r" u " + VEC + r" uf " + VEC + r" Re (" + NUM + r")"
                      + r" Cd (" + NUM + r") aDrag " + VEC, line)
        if not m:
            continue
        g = [float(m.group(i)) for i in range(1, 18)]
        key = (int(g[0]), int(g[1]))
        samples[key] = {
            "t": g[2], "x": tuple(g[3:6]), "u": tuple(g[6:9]),
            "uf": tuple(g[9:12]), "re": g[12], "cd": g[13],
            "aDrag": tuple(g[14:17]),
        }
    return samples


def ends_of(output):
    ends = {}
    for line in output.splitlines():
        m = re.search(r"dpm particle (\d+) end: status (\w+) t (" + NUM + r") x "
                      + VEC + r" u " + VEC, line)
        if m:
            ends[int(m.group(1))] = {
                "status": m.group(2), "t": float(m.group(3)),
                "x": tuple(float(m.group(i)) for i in (4, 5, 6)),
                "u": tuple(float(m.group(i)) for i in (7, 8, 9)),
            }
    return ends


def diff(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mag(v):
    return math.sqrt(sum(x * x for x in v))


def cd_of(reynolds):
    if reynolds <= 1000:
        return (24.0 / reynolds) * (1.0 + 0.15 * reynolds ** 0.687)
    return 0.44


class TestDpmChannel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if len(sys.argv) < 3:
            raise SystemExit("usage: check_dpm_channel.py <exe> <case.dict>")
        cls.exe = os.path.abspath(sys.argv[1])
        cls.tmp = tempfile.TemporaryDirectory(prefix="cf_dpm_")
        with open(CASE) as fh:
            cls.text = fh.read()
        cls.dpm = block(cls.text, "dpm")
        cls.dt = scalar(cls.dpm, "dt")
        cls.rhoF = scalar(cls.dpm, "carrierRho")
        cls.mu = cls.rhoF * scalar(cls.text, "nu")
        cls.gravity = key_vec(cls.dpm, "gravity")
        cls.traceSteps = int(scalar(cls.dpm, "traceSteps"))
        cls.blockDict = block(cls.text, "block")
        cls.xMax = key_vec(cls.blockDict, "max")[0]
        cls.yMin = key_vec(cls.blockDict, "min")[1]

        cls.completed = subprocess.run([cls.exe, CASE], cwd=cls.tmp.name, text=True,
                                 capture_output=True, timeout=300)
        if cls.completed.returncode != 0:
            raise RuntimeError("dpm case failed:\n" + cls.completed.stdout[-3000:]
                               + cls.completed.stderr[-3000:])
        cls.launches = launch_of(cls.completed.stdout)
        cls.samples = samples_of(cls.completed.stdout)
        cls.ends = ends_of(cls.completed.stdout)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_carrier_is_frozen_to_the_bit(self):
        line = next(l for l in self.completed.stdout.splitlines()
                    if "dpm carrier frozen" in l)
        before = dict((m.group(1), float(m.group(2)))
                      for m in re.finditer(r"\b(uMaxBefore|uMaxAfter|sumUBefore|sumUAfter) ("
                                           + NUM + r")", line))
        self.assertEqual(before["uMaxBefore"], before["uMaxAfter"])
        self.assertEqual(before["sumUBefore"], before["sumUAfter"])

    def test_launch_state_and_terminal_status(self):
        self.assertEqual(len(self.launches), 2)
        self.assertEqual(self.launches[0]["d"], 1e-3)
        self.assertEqual(self.launches[0]["rho"], 2500.0)
        self.assertEqual(self.launches[1]["d"], 3e-3)
        self.assertEqual(self.launches[1]["rho"], 8000.0)
        self.assertEqual(self.ends[0]["status"], "outlet")
        self.assertEqual(self.ends[1]["status"], "wall")
        # the outlet crossing is detected at a discrete step, not on the face
        self.assertLess(abs(self.ends[0]["x"][0] - (self.xMax - 0.5e-3)), 1e-5)
        # likewise the wall hit is recorded at the step that crossed the radius
        self.assertLess(abs(self.ends[1]["x"][1] - (self.yMin + 1.5e-3)), 1e-4)
        line = next(l for l in self.completed.stdout.splitlines() if "dpm summary" in l)
        self.assertIn("particles 2 outlet 1 wall 1 expired 0", line)

    def test_schiller_naumann_and_reynolds_from_case(self):
        for (pid, step), s in self.samples.items():
            if pid != 0:
                continue
            rel = diff(s["uf"], s["u"])
            wantRe = self.rhoF * mag(rel) * self.launches[0]["d"] / self.mu
            self.assertLess(abs(s["re"] - wantRe) / wantRe, 1e-9)
            self.assertLess(abs(s["cd"] - cd_of(s["re"])) / s["cd"], 1e-9)
            rate = 3.0 * self.rhoF * s["cd"] / (4.0 * self.launches[0]["rho"]
                                                * self.launches[0]["d"])
            wantA = tuple(rate * r for r in rel)
            self.assertLess(mag(diff(s["aDrag"], wantA)), 1e-7)

    def test_explicit_euler_velocity_and_position(self):
        rhoP = self.launches[0]["rho"]
        body = tuple((1.0 - self.rhoF / rhoP) * g for g in self.gravity)
        for step in range(self.traceSteps - 1):
            a, b = self.samples[(0, step)], self.samples[(0, step + 1)]
            want_u = tuple(a["u"][i] + self.dt * (a["aDrag"][i] + body[i])
                           for i in range(3))
            self.assertLess(mag(diff(b["u"], want_u)), 1e-10)
            # the app advances velocity first, then position with the new one
            want_x = tuple(a["x"][i] + self.dt * b["u"][i] for i in range(3))
            self.assertLess(mag(diff(b["x"], want_x)), 1e-12)

    def test_particle_approaches_carrier(self):
        first, last = self.samples[(0, 0)], self.ends[0]
        self.assertLess(first["u"][0], 1e-12)
        self.assertGreater(last["u"][0], 0.05)
        self.assertLess(abs(last["u"][0] - 0.12), 0.05)

    def test_declaration_rejections(self):
        cases = [
            (self.text.replace("drag schillerNaumann", "drag unknownLaw"),
             "drag law 'unknownLaw' is not implemented"),
            (re.sub(r"(stream0\s*\{.*?diameter\s+)-?[\d.eE+-]+", r"\g<1>-1",
                    self.text, count=1, flags=re.S),
             "diameter must be finite and positive"),
            (re.sub(r"\ndpm\s*\{.*?\n\}\n", "\n", self.text, count=1, flags=re.S),
             "needs a dpm{} declaration"),
        ]
        for text, needle in cases:
            path = os.path.join(self.tmp.name, "reject.dict")
            with open(path, "w") as fh:
                fh.write(text)
            r = subprocess.run([self.exe, path], cwd=self.tmp.name, text=True,
                               capture_output=True, timeout=300)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn(needle, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
