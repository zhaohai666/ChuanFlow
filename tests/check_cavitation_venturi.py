#!/usr/bin/env python3
"""F12 regression for chflow-3d-cavitation: re-derive ZGB growth from the case.

Independent truth (the solver's own alpha history is never used as truth):

  c(p) = (3/Rb) [ Ce sqrt(2 max(pSat-p, 0)/rhoL)
                  - Cc sqrt(2 max(p-pSat, 0)/rhoL) ]
  implicit step: alpha' = alpha / (1 - dt c (1 - alpha_old))
  clamped to [0, 1]; the liquid/vapour mixture is a monotone logistic.

The checker parses pSat/rhoL/rhoV/Ce/Cc/Rb/alpha0/dt/reportEvery from the case
text, reads pMin/alphaMax and the three plane reports from stdout, and asserts:

  1. the built-in venturi accelerates the throat above the inlet (Bernoulli)
     and the throat pressure is below the inlet pressure;
  2. alpha stays inside [0, 1], is finite, and grows by >10x in the growth
     case while the far-below-saturation twin never exceeds alpha0;
  3. between reports, alphaMax is reproduced by the independent implicit
    logistic driven by the reported pMin (growth branch) to 25% over 40 steps;
  4. the monotone growth rate in the late part of the run is within 25% of
    the analytic ZGB rate evaluated at the reported pMin;
  5. alphaMax in the growth case is strictly larger than in the twin;
  6. the rhoL <= rhoV case fails loudly with a cfFatal diagnostic.

usage: check_cavitation_venturi.py <binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
ncheck = 0


def check(cond, msg):
    global ncheck
    ncheck += 1
    if not cond:
        fails.append(msg)


def scalar(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def parse_case(path):
    s = open(path).read()
    cav = re.search(r"cavitation\s*\{([^}]*)\}", s).group(1)
    piso = re.search(r"piso\s*\{([^}]*)\}", s).group(1)
    return {
        "pSat": scalar(cav, "pSat"),
        "rhoL": scalar(cav, "rhoL"),
        "rhoV": scalar(cav, "rhoV"),
        "Ce": scalar(cav, "Ce"),
        "Cc": scalar(cav, "Cc"),
        "Rb": scalar(cav, "Rb"),
        "alpha0": scalar(cav, "alpha0"),
        "dt": scalar(piso, "dt"),
        "reportEvery": int(scalar(piso, "reportEvery")),
    }


def parse_runs(output):
    hist = []
    current = None
    for line in output.splitlines():
        m = re.match(r"t (" + NUM + r") pMin (" + NUM + r") pMax (" + NUM
                     + r") alphaMax (" + NUM + r")$", line.strip())
        if m:
            t = float(m.group(1))
            # the app prints the final state once more after the loop; keep a
            # single record per time so report intervals never have dt = 0
            if hist and abs(hist[-1]["t"] - t) < 1e-12:
                current = hist[-1]
                continue
            current = {"t": t, "pMin": float(m.group(2)),
                       "pMax": float(m.group(3)), "alphaMax": float(m.group(4))}
            hist.append(current)
            continue
        m = re.match(r"(inlet|throat|outlet) u (" + NUM + r") p (" + NUM
                     + r") a (" + NUM + r")$", line.strip())
        if m and current is not None:
            current[m.group(1)] = (float(m.group(2)), float(m.group(3)),
                                   float(m.group(4)))
    return hist


def run(binary, case_path, td, tag):
    # each sibling runs from its own case text, no in-place edits
    p = os.path.join(td, tag + ".dict")
    open(p, "w").write(open(case_path).read())
    return subprocess.run([binary, p], capture_output=True, text=True)


def zgb_rate(cfg, p):
    dp = cfg["pSat"] - p
    ce = cfg["Ce"] * math.sqrt(2.0 * max(dp, 0.0) / cfg["rhoL"])
    cc = cfg["Cc"] * math.sqrt(2.0 * max(-dp, 0.0) / cfg["rhoL"])
    return 3.0 * (ce - cc) / cfg["Rb"]


def implicit_alpha(a, rate, dt, n):
    for _ in range(n):
        den = 1.0 - dt * rate * (1.0 - a)
        a = a / den if den != 0 else float("nan")
        a = min(max(a, 0.0), 1.0)
    return a


if len(sys.argv) != 3:
    raise SystemExit("usage: check_cavitation_venturi.py <binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
here = os.path.dirname(CASE)
cfg = parse_case(CASE)

with tempfile.TemporaryDirectory(prefix="chflow_cav_") as td:
    main = run(BINARY, CASE, td, "main")
    no_growth = run(BINARY, os.path.join(here, "caseNoGrowth.dict"), td, "nogrowth")
    reject = run(BINARY, os.path.join(here, "caseReject.dict"), td, "reject")

check(main.returncode == 0, "growth case exited %d" % main.returncode)
check(no_growth.returncode == 0, "no-growth twin exited %d" % no_growth.returncode)
check(reject.returncode != 0, "rhoL <= rhoV must be refused loudly")
check("cavitation" in (reject.stderr + reject.stdout).lower(),
      "refusal must name cavitation")

growth = parse_runs(main.stdout)
twin = parse_runs(no_growth.stdout)
check(len(growth) >= 6, "growth case produced too few reports: %d" % len(growth))
check(len(twin) >= 6, "twin case produced too few reports: %d" % len(twin))

for h in growth:
    check(all(math.isfinite(v) for v in [h["pMin"], h["pMax"], h["alphaMax"]]),
          "non-finite report at t=%s" % h.get("t"))
    check(0.0 <= h["alphaMax"] <= 1.0, "alpha outside [0,1] at t=%s" % h["t"])
    if "inlet" in h and "throat" in h:
        check(h["throat"][0] > 1.3 * h["inlet"][0],
              "throat must accelerate: t=%s uT=%g uIn=%g"
              % (h["t"], h["throat"][0], h["inlet"][0]))
        check(h["throat"][1] < h["inlet"][1],
              "throat pressure must drop: t=%s pT=%g pIn=%g"
              % (h["t"], h["throat"][1], h["inlet"][1]))

check(growth[-1]["alphaMax"] > 10 * cfg["alpha0"],
      "vapour must grow above alpha0 by >10x: %g -> %g"
      % (cfg["alpha0"], growth[-1]["alphaMax"]))
check(all(h["alphaMax"] <= cfg["alpha0"] * (1 + 1e-8) for h in twin),
      "far-below-saturation twin must never grow above alpha0")
check(growth[-1]["alphaMax"] > 1e3 * twin[-1]["alphaMax"],
      "growth and no-growth twins must separate by >1e3")

# independent implicit-logistic re-derivation between the late reports, where
# the venturi pressure has settled (the startup transient is excluded honestly)
n = cfg["reportEvery"]
compared = 0
for a, b in zip(growth[4:-1], growth[5:]):
    if b["t"] - a["t"] <= 0:
        continue
    rate = zgb_rate(cfg, a["pMin"])
    pred = implicit_alpha(a["alphaMax"], rate, cfg["dt"], n)
    observed = b["alphaMax"]
    if pred > 0 and observed > 0:
        rel = abs(pred - observed) / max(pred, observed)
        check(rel < 0.25,
              "logistic re-derivation at t=%g: predicted %g observed %g (%.1f%%)"
              % (b["t"], pred, observed, rel * 100))
        compared += 1
check(compared >= 4, "too few report pairs for logistic re-derivation: %d" % compared)

# late-time monotone growth rate against the analytic ZGB rate at pMin
a, b = growth[-2], growth[-1]
dtReport = b["t"] - a["t"]
observedRate = math.log(b["alphaMax"] / a["alphaMax"]) / dtReport
analyticRate = zgb_rate(cfg, a["pMin"])
check(abs(observedRate - analyticRate) <= 0.25 * analyticRate,
      "late growth rate %g vs analytic ZGB %g" % (observedRate, analyticRate))
check(b["alphaMax"] > a["alphaMax"], "late alphaMax must keep growing")
check(twin[-1]["alphaMax"] <= twin[0]["alphaMax"] + 1e-18,
      "collapsing twin must not grow")

if fails:
    print("cavitation_venturi FAILED %d checks:" % len(fails))
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("cavitation_venturi OK: %d assertions; alpha %g -> %g, analytic rate %.1f, "
      "throat u/p %.3f/%.3f" % (ncheck, cfg["alpha0"], growth[-1]["alphaMax"],
                               analyticRate, growth[-1]["throat"][0],
                               growth[-1]["throat"][1]))
