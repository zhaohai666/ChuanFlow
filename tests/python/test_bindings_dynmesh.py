"""Dynamic-mesh 6DoF family Python binding tests (CTest `python.dynmesh`).

Covers `chflow_core.solve_dynmesh_case`: GCL-corrected relative-flux scalar
transport on a rigidly translating mesh, plus the translational (mass + linear
drag + constant force) and rotational (inertia + rotational drag + stiffness)
rigid-body ODEs. Closed forms follow check_dynamic_mesh_6dof.py:

  * translation:  m v' = F - k v, v -> F/k;
  * rotation:     I theta'' + c theta' + K theta = 0, damped harmonic oscillator;
  * a uniform scalar is invariant under the GCL-corrected relative flux;
  * the static twin (zero mesh velocity, zero force) leaves the body at rest;
  * loud refusals of negative mass, drag and inertia.
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
CASES = ROOT / "tests" / "cases" / "dynamicMesh6Dof"
MAIN = CASES / "case.dict"
STATIC = CASES / "caseStatic.dict"
REJECT = CASES / "caseReject.dict"

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
    dm = block(txt, "dynamicMesh")
    return {
        "meshV": vectors(dm, "meshVelocity"), "F": vectors(dm, "bodyForce"),
        "m": scalar(dm, "mass"), "k": scalar(dm, "drag"),
        "I": scalar(dm, "inertia"), "c": scalar(dm, "rotDrag"),
        "K": scalar(dm, "rotStiffness"), "theta0": scalar(dm, "theta0"),
        "omega0": scalar(dm, "omega0"), "phi0": scalar(dm, "initialValue"),
        "dt": scalar(dm, "dt"), "steps": int(scalar(dm, "steps")),
        "reportEvery": int(scalar(dm, "reportEvery")),
    }


def body_exact(t, m, k, F, I, c, K, theta0, omega0):
    if k > 0:
        tau = m / k
        vt = F[0] / k
        decay = math.exp(-t / tau)
        v = vt * (1.0 - decay)
        x = vt * (t - tau * (1.0 - decay))
    else:
        v = F[0] * t / m
        x = F[0] * t * t / (2.0 * m)

    if K > 0:
        alpha = c / (2.0 * I)
        disc = K / I - alpha * alpha
        if disc > 0:
            wd = math.sqrt(disc)
            A = theta0
            B = (omega0 + alpha * theta0) / wd
            env = math.exp(-alpha * t)
            th = env * (A * math.cos(wd * t) + B * math.sin(wd * t))
            om = env * (omega0 * math.cos(wd * t) - (alpha * B + A * wd) * math.sin(wd * t))
        elif disc == 0:
            A = theta0
            B = omega0 + alpha * theta0
            env = math.exp(-alpha * t)
            th = env * (A + B * t)
            om = env * (B - alpha * (A + B * t))
        else:
            r = math.sqrt(-disc)
            l1, l2 = -alpha + r, -alpha - r
            B = (omega0 - l1 * theta0) / (l2 - l1)
            A = theta0 - B
            th = A * math.exp(l1 * t) + B * math.exp(l2 * t)
            om = A * l1 * math.exp(l1 * t) + B * l2 * math.exp(l2 * t)
    elif c > 0:
        om = omega0 * math.exp(-c * t / I)
        th = theta0 + omega0 * I / c * (1.0 - math.exp(-c * t / I))
    else:
        om = omega0
        th = theta0 + omega0 * t
    return x, v, th, om


def solve_text(text):
    with tempfile.TemporaryDirectory(prefix="cf_dynmesh_py_") as td:
        p = os.path.join(td, "case.dict")
        open(p, "w").write(text)
        return chflow_core.solve_dynmesh_case(p)


class MainCaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(MAIN).read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_dynmesh_case(str(MAIN))

    def test_records_count(self):
        c = self.cfg
        expected = c["steps"] + 1 if c["reportEvery"] == 1 else \
            1 + (c["steps"] + c["reportEvery"] - 1) // c["reportEvery"]
        self.assertEqual(len(self.res.history), expected)

    def test_uniform_scalar_and_gcl(self):
        for rec in self.res.history:
            self.assertTrue(close(rec.scalar_mean, self.cfg["phi0"], 0.0, 1e-9))
            self.assertLess(rec.mass_balance, 1e-9)

    def test_body_closed_forms(self):
        c = self.cfg
        for rec in self.res.history:
            ex, ev, eth, eom = body_exact(
                rec.t, c["m"], c["k"], c["F"], c["I"], c["c"], c["K"],
                c["theta0"], c["omega0"])
            self.assertTrue(close(rec.x, ex, 1e-6, 1e-7), f"t={rec.t} x")
            self.assertTrue(close(rec.v, ev, 1e-6, 1e-7), f"t={rec.t} v")
            self.assertTrue(close(rec.theta, eth, 1e-6, 1e-7), f"t={rec.t} theta")
            self.assertTrue(close(rec.omega, eom, 1e-6, 1e-7), f"t={rec.t} omega")

    def test_phi_field_uniform(self):
        phi = self.res.phi
        self.assertEqual(phi.shape, (self.res.n_cells(),))
        self.assertTrue(all(close(v, self.cfg["phi0"], 0.0, 1e-9) for v in phi))


class StaticTwinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.txt = open(STATIC).read()
        cls.cfg = parse(cls.txt)
        cls.res = chflow_core.solve_dynmesh_case(str(STATIC))

    def test_body_at_rest(self):
        self.assertEqual(len(self.res.history), self.cfg["steps"] + 1)
        for rec in self.res.history:
            self.assertTrue(close(rec.x, 0.0, 0.0, 1e-12))
            self.assertTrue(close(rec.v, 0.0, 0.0, 1e-12))
            self.assertTrue(close(rec.theta, 0.0, 0.0, 1e-12))
            self.assertTrue(close(rec.omega, 0.0, 0.0, 1e-12))
            self.assertTrue(close(rec.scalar_mean, self.cfg["phi0"], 0.0, 1e-9))
            self.assertLess(rec.mass_balance, 1e-12)


class RefusalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = open(MAIN).read()

    def _reject(self, text, needle):
        with self.assertRaises(RuntimeError) as cm:
            solve_text(text)
        self.assertIn(needle, str(cm.exception))

    def test_negative_mass(self):
        self._reject(open(REJECT).read(), "mass must be finite and positive")

    def test_negative_drag(self):
        txt = re.sub(r"(\bdrag\s+)" + NUM + r"(\s*;)", r"\g<1>-1\g<2>",
                     self.text, count=1)
        self._reject(txt, "drag must be finite and non-negative")

    def test_negative_inertia(self):
        txt = re.sub(r"(\binertia\s+)" + NUM + r"(\s*;)", r"\g<1>-0.5\g<2>",
                     self.text, count=1)
        self._reject(txt, "inertia must be finite and positive")

    def test_missing_case_raises(self):
        with self.assertRaises(RuntimeError):
            chflow_core.solve_dynmesh_case(str(ROOT / "tests" / "cases" / "no-such.dict"))


if __name__ == "__main__":
    sys.exit(unittest.main(verbosity=2))
