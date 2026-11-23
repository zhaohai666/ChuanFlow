#!/usr/bin/env python3
"""Batch D / F10: independent grading of the two-step chain-branching 0D
reactor (ignition and extinction of detailed chemistry).

Nothing is trusted from the solver's own arithmetic. Every constant is read
from the case text, the three ODEs

    dT/dt  = q r1 - cooling (T - Tenv)
    dYf/dt = -r1                       r1 = A1 exp(-E1/(R T)) Yf Yr
    dYr/dt = +r1 - r2                  r2 = A2 exp(-E2/(R T)) Yr

are re-integrated with this script's own RK4 on a four-times finer substep,
and every printed (t, T, Yf, Yr) sample is compared at 1e-3 relative. The
solver's reported ignition delay and extinction Damkoehler threshold are
re-derived here by an independent bisection on the same ODEs. The adiabatic
invariant T + q Yf, species bounds, ignition-before-threshold, radical decay
for the extinguish variant and loud rejections are all asserted.

usage: check_combustion_detailed.py <solver> <case> <ignite> <extinguish> <reject>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
nchecks = 0


def check(cond, msg):
    global nchecks
    nchecks += 1
    if not cond:
        fails.append(msg)


def close(a, b, rel, ab=0.0):
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), ab)


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("no " + name + " block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced " + name + " block")


def parse_case(path):
    txt = open(path).read()
    cb = block(txt, "detailedChemistry")
    c = {k: scalar(cb, k) for k in (
        "A1", "E1", "A2", "E2", "q", "T0", "gasConstant", "Yf0", "Yr0",
        "cooling", "Tenv", "dt", "ignitionRise", "tIgnMax")}
    c["nSteps"] = int(scalar(cb, "nSteps"))
    c["nSub"] = int(scalar(cb, "nSub"))
    return c, txt


def kv(T, A, E, R):
    return A if E == 0 else A * math.exp(-E / (R * T))


def rhs(T, Yf, Yr, c, beta):
    a = kv(T, c["A1"], c["E1"], c["gasConstant"])
    b = kv(T, c["A2"], c["E2"], c["gasConstant"])
    r1 = a * Yf * Yr
    r2 = b * Yr
    return (c["q"] * r1 - beta * (T - c["Tenv"]), -r1, r1 - r2)


def rk4(s, h, c, beta):
    k1 = rhs(*s, c, beta)
    s2 = tuple(s[i] + 0.5 * h * k1[i] for i in range(3))
    k2 = rhs(*s2, c, beta)
    s3 = tuple(s[i] + 0.5 * h * k2[i] for i in range(3))
    k3 = rhs(*s3, c, beta)
    s4 = tuple(s[i] + h * k3[i] for i in range(3))
    k4 = rhs(*s4, c, beta)
    return tuple(s[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                 for i in range(3))


def run(c, beta, nsub, collect=True):
    s = (c["T0"], c["Yf0"], c["Yr0"])
    h = c["dt"] / nsub
    trig = c["T0"] + c["ignitionRise"] * c["q"] * c["Yf0"]
    peak = s[0]
    tign = -1.0
    hist = [(0.0,) + s] if collect else None
    for i in range(c["nSteps"]):
        for _ in range(nsub):
            s = rk4(s, h, c, beta)
        t = (i + 1) * c["dt"]
        if s[0] > peak:
            peak = s[0]
        if tign < 0 and s[0] >= trig:
            tign = t
        if collect:
            hist.append((t,) + s)
    return hist, peak, tign


def diagnose(c):
    a0 = kv(c["T0"], c["A1"], c["E1"], c["gasConstant"])
    b0 = kv(c["T0"], c["A2"], c["E2"], c["gasConstant"])
    branch_da = a0 * c["Yf0"] / b0 if b0 > 0 else math.inf
    trig = c["T0"] + c["ignitionRise"] * c["q"] * c["Yf0"]

    def ignites(beta):
        return run(c, beta, c["nSub"], collect=False)[1] >= trig

    if not ignites(0):
        beta_crit = 0.0  # even adiabatic chemistry cannot ignite
    else:
        lo, hi = 0.0, 1.0
        while ignites(hi) and hi < 1e6:
            hi *= 2
        for _ in range(36):
            mid = 0.5 * (lo + hi)
            if ignites(mid):
                lo = mid
            else:
                hi = mid
        beta_crit = 0.5 * (lo + hi)
    loss = b0 + beta_crit
    da_crit = a0 * c["Yf0"] / loss if loss > 0 else math.inf
    _, peak, tign = run(c, c["cooling"], c["nSub"], collect=False)
    return dict(branchDa=branch_da, betaCrit=beta_crit, DaCrit=da_crit,
                ignitionDelay=tign, peakT=peak, ignited=1 if tign >= 0 else 0)


def parse_output(out):
    hist = []
    params = {}
    numbers = {}
    crow = re.compile(r"^COMBUST\s+(" + NUM + r")\s+(" + NUM + r")\s+("
                      + NUM + r")\s+(" + NUM + r")\s*$")
    for line in out.splitlines():
        line = line.strip()
        m = crow.match(line)
        if m:
            hist.append(tuple(float(m.group(i)) for i in range(1, 5)))
        elif line.startswith("DETAILED params") or line.startswith("DETAILED number"):
            toks = line.split()
            target = params if "params" in line else numbers
            for i in range(2, len(toks) - 1, 2):
                try:
                    target[toks[i]] = float(toks[i + 1])
                except ValueError:
                    pass
    return hist, params, numbers


def solve(path):
    r = subprocess.run([BINARY, path], capture_output=True, text=True)
    return r


if len(sys.argv) != 6:
    raise SystemExit("usage: check_combustion_detailed.py <solver> <case> "
                     "<ignite> <extinguish> <reject>")
BINARY = sys.argv[1]
P_CASE, P_IGN, P_EXT, P_REJ = sys.argv[2:6]

CONFIGS = []
for path, tag in ((P_CASE, "case"), (P_IGN, "ignite"), (P_EXT, "extinguish")):
    cfg, text = parse_case(path)
    r = solve(path)
    check(r.returncode == 0,
          "%s: solver returned %d: %s" % (tag, r.returncode,
                                          (r.stdout + r.stderr)[-400:]))
    hist, params, numbers = parse_output(r.stdout)
    CONFIGS.append((tag, cfg, text, hist, params, numbers))

    # The solver must echo the case-declared constants, not a canned history.
    for key in ("A1", "E1", "A2", "E2", "q", "T0", "gasConstant", "Yf0",
                "Yr0", "cooling", "Tenv", "dt", "ignitionRise"):
        check(key in params and close(params[key], cfg[key], 1e-12, 0.0),
              "%s: params echo %s=%s != case %g" %
              (tag, key, params.get(key), cfg[key]))
    for key in ("nSteps", "nSub"):
        check(key in params and params[key] == cfg[key],
              "%s: params echo %s=%s != case %d" %
              (tag, key, params.get(key), cfg[key]))

    # ---- independently integrated history (4x finer RK4) ----------------
    fine = run(cfg, cfg["cooling"], 4 * cfg["nSub"], collect=True)[0]
    check(len(hist) == cfg["nSteps"] + 1,
          "%s: expected %d history rows, got %d" %
          (tag, cfg["nSteps"] + 1, len(hist)))
    check(len(fine) == len(hist),
          "%s: independent history length mismatch" % tag)
    for (t, T, Yf, Yr), (tf, Tf, Yff, Yrf) in zip(hist, fine):
        check(close(t, tf, 1e-12, 1e-12), "%s: t %g != %g" % (tag, t, tf))
        check(close(T, Tf, 1e-3, 1e-6),
              "%s: T(t=%g) %g != independent %g" % (tag, t, T, Tf))
        check(close(Yf, Yff, 1e-3, 1e-9),
              "%s: Yf(t=%g) %g != independent %g" % (tag, t, Yf, Yff))
        check(close(Yr, Yrf, 1e-3, 1e-9),
              "%s: Yr(t=%g) %g != independent %g" % (tag, t, Yr, Yrf))

    # ---- species / energy bounds ----------------------------------------
    inv0 = cfg["T0"] + cfg["q"] * cfg["Yf0"]
    for (t, T, Yf, Yr) in hist:
        check(-1e-10 <= Yf <= cfg["Yf0"] + 1e-10,
              "%s: Yf=%g outside [0,%g]" % (tag, Yf, cfg["Yf0"]))
        check(-1e-14 <= Yr <= cfg["Yf0"] + cfg["Yr0"] + 1e-6,
              "%s: Yr=%g outside species bounds" % (tag, Yr))
        check(T > 0 and T <= inv0 + 1e-3,
              "%s: T=%g outside (0,%g]" % (tag, T, inv0))
        if cfg["cooling"] == 0:
            check(abs((T + cfg["q"] * Yf) - inv0)
                  <= 1e-6 * max(1.0, cfg["q"] * cfg["Yf0"]),
                  "%s: adiabatic invariant T+qYf broken at t=%g (%g vs %g)" %
                  (tag, t, T + cfg["q"] * Yf, inv0))

    # ---- reported diagnostic numbers are real computed outcomes ---------
    d = diagnose(cfg)
    for key, val in d.items():
        check(key in numbers and close(numbers.get(key, math.nan), val, 1e-3, 1e-9),
              "%s: reported %s=%s != independent %g" %
              (tag, key, numbers.get(key), val))

# ---- ignition cases: before the declared threshold -----------------------
for tag, cfg, text, hist, params, numbers in CONFIGS:
    if tag == "extinguish":
        continue
    trig = cfg["T0"] + cfg["ignitionRise"] * cfg["q"] * cfg["Yf0"]
    d = diagnose(cfg)
    check(d["branchDa"] > 1, "%s: branchDa %g <= 1 (branching must win)" %
          (tag, d["branchDa"]))
    check(numbers.get("ignited") == 1, "%s: solver did not report ignition" % tag)
    check(0 <= d["ignitionDelay"] < cfg["tIgnMax"],
          "%s: ignition delay %g not before declared tIgnMax %g" %
          (tag, d["ignitionDelay"], cfg["tIgnMax"]))
    check(d["peakT"] >= trig and hist[-1][1] >= trig,
          "%s: temperature never crosses ignition trigger %g" % (tag, trig))

# ---- extinguish case: radical decays, no runaway --------------------------
ext = CONFIGS[2]
cfg, hist, numbers = ext[1], ext[3], ext[5]
d = diagnose(cfg)
trig = cfg["T0"] + cfg["ignitionRise"] * cfg["q"] * cfg["Yf0"]
check(d["branchDa"] < 1, "extinguish: branchDa %g >= 1 (termination must win)" %
      d["branchDa"])
check(numbers.get("ignited") == 0, "extinguish: solver wrongly reports ignition")
check(all(hist[i + 1][3] <= hist[i][3] * (1 + 1e-6)
          for i in range(len(hist) - 1)),
      "extinguish: radical fraction is not monotonically decaying")
check(hist[-1][3] < cfg["Yr0"],
      "extinguish: final radical %g not below seed %g" % (hist[-1][3], cfg["Yr0"]))
check(max(s[1] for s in hist) < trig,
      "extinguish: temperature reached the ignition trigger (runaway)")
check(abs(hist[-1][1] - cfg["Tenv"]) < 1.0,
      "extinguish: final T %g did not settle at Tenv %g" %
      (hist[-1][1], cfg["Tenv"]))

# ---- loud rejections ------------------------------------------------------
BASE = CONFIGS[0][2]
repls = [
    (open(P_REJ).read(), "A1", "reject_case"),
    (BASE.replace("A1 2000;", "A1 -2000;"), "A1", "neg_A1"),
    (BASE.replace("A2 40;", "A2 -40;"), "A2", "neg_A2"),
    (BASE.replace("E1 60000;", "E1 -1;"), "E1", "neg_E1"),
    (BASE.replace("E2 30000;", "E2 -1;"), "E2", "neg_E2"),
    (BASE.replace("Yf0 1;", "Yf0 1.5;"), "Yf0", "bad_Yf0"),
    (BASE.replace("Yr0 1e-4;", "Yr0 -1e-4;"), "Yr0", "bad_Yr0"),
    (re.sub(r"detailedChemistry\s*\{.*?^\}", "", BASE, count=1,
            flags=re.M | re.S), "detailedChemistry", "no_block"),
]
for txt, needle, tag in repls:
    with tempfile.TemporaryDirectory(prefix="cf_f10_") as td:
        p = os.path.join(td, tag + ".dict")
        open(p, "w").write(txt)
        rr = solve(p)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d): %s" %
              (needle, rr.returncode, (rr.stdout + rr.stderr)[-200:]))

if fails:
    print("combustionDetailed FAILED %d of %d checks:" % (len(fails), nchecks))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("combustionDetailed PASS: %d assertions; independent fine-RK4 history "
      "(1e-3), adiabatic invariant, species bounds, computed ignition delay "
      "and extinction Damkoehler threshold, two ignitions before threshold, "
      "radical-decay extinction, 8 loud refusals" % nchecks)
