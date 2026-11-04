#!/usr/bin/env python3
"""Verification of nodal analysis (chflow-1d-nodal).

The shipped cases cross a Vogel inflow curve with the Beggs-Brill
vertical-lift curve of one producing well; the operating point is the largest
root of p_vlp(q) - p_ipr(q). The correlation constants are anchored at the
unit level (tests/unit/nodalTests.cpp); this script checks the assembled
solve, over a living gas well, its dead-oil twin and a dead well:

  * every case converges and the flowing ones flow (flows=true, iterations>0),
  * the operating rate reproduces the Vogel identity recomputed here in
    closed form from the printed bottom-hole pressure (the solve's own
    residual is invisible at 12 digits beyond ~1e-6 relative),
  * the operating point is a crossing, not a tangent: the printed q and pwf
    sit inside the IPR capacity, above the wellhead pressure and below the
    reservoir, and the flowing bottom hole sits UNDER the static column
    (free gas lightens the column),
  * golden anchors: the measured operating points are asserted against the
    values recorded at calibration,
  * the dead-oil twin (GOR 0, everything else identical) produces LESS at a
    HIGHER bottom-hole pressure: gas lifts the well,
  * the dead well (reservoir below the static column) reports rate 0,
    flows=false and a bottom hole exactly on the static column.

usage: check_nodal.py <path-to-chflow-1d-nodal> <vogel-case.dict> <deadoil-case.dict> <deadwell-case.dict>
"""
import re
import subprocess
import sys
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(__doc__)

exe = sys.argv[1]
cases = {name: sys.argv[i] for i, name in
         enumerate(("vogel", "deadoil", "deadwell"), start=2)}

failures = []
notes = []


def check(cond, msg):
    if cond:
        notes.append("ok:   " + msg)
    else:
        failures.append(msg)


def run(name):
    path = cases[name]
    proc = subprocess.run([exe, str(path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed on {name} rc={proc.returncode}\n"
                         f"{proc.stdout}\n{proc.stderr}")
    out = {}
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        if toks[0] == "column":
            out["static"] = float(toks[2])
        elif toks[0] == "well" and toks[1] == "q":
            out["q"], out["pwf"], out["pth"] = float(toks[2]), float(toks[4]), float(toks[6])
        elif toks[0] == "well" and toks[1] == "holdup":
            out["holdup"], out["regime"] = float(toks[2]), toks[4]
        elif toks[0] == "nodal":
            out["iterations"] = int(toks[2].split("=")[1])
            out["flows"] = toks[3].split("=")[1] == "true"
    return out


def ipr_constants(path):
    """The IPR numbers the case dict asks for (vogel form only)."""
    txt = Path(path).read_text()
    pres = float(re.search(r"reservoirPressure\s+([0-9.eE+-]+)", txt).group(1))
    qmax = float(re.search(r"maxRate\s+([0-9.eE+-]+)", txt).group(1))
    return pres, qmax


def vogel_rate(pwf, pres, qmax):
    x = pwf / pres
    return qmax * (1.0 - 0.2 * x - 0.8 * x * x)


# the operating points recorded when this checker was calibrated; the solver
# is deterministic, so any drift means the physics moved
GOLDEN = {
    "vogel": {"q": 0.00290785579167, "pwf": 15120468.4282,
              "static": 18642604.4027},
    "deadoil": {"q": 0.00057633047675, "pwf": 20796742.5522,
                "static": 20769607.9235},
    "deadwell": {"static": 18642604.4027},
}

results = {}
for name in cases:
    results[name] = run(name)
    notes.append(f"{name}: q {results[name]['q']:.9g} pwf "
                 f"{results[name]['pwf']:.9g} static "
                 f"{results[name]['static']:.9g} holdup "
                 f"{results[name]['holdup']:.4g} regime "
                 f"{results[name]['regime']} iterations "
                 f"{results[name]['iterations']} flows "
                 f"{results[name]['flows']}")

# ---- the living gas well: a crossing that satisfies both curves ----
v = results["vogel"]
pres, qmax = ipr_constants(cases["vogel"])
check(v["flows"] and v["iterations"] > 0,
      "the gas well flows after a real secant solve, not by fiat")
check(0.0 < v["q"] < qmax,
      "the operating rate sits inside the IPR capacity (0, qmax)")
check(v["pth"] < v["pwf"] < pres,
      "the bottom hole is between the wellhead and the reservoir")
q_vogel = vogel_rate(v["pwf"], pres, qmax)
check(abs(v["q"] / q_vogel - 1.0) < 1e-5,
      f"the operating rate satisfies the Vogel identity recomputed here "
      f"(q {v['q']:.9g} vs {q_vogel:.9g}, {100.0 * (v['q'] / q_vogel - 1.0):+.3g}%)")
check(v["pwf"] < v["static"],
      f"the flowing bottom hole sits under the static column: free gas "
      f"lightens it ({v['static'] - v['pwf']:.3g} Pa lighter)")
check(0.0 < v["holdup"] < 1.0 and v["regime"] in
      {"segregated", "transition", "intermittent", "distributed"},
      "the gas-bearing column reports a two-phase holdup and a named regime")

# ---- the dead-oil twin: gas lifts the well ----
d = results["deadoil"]
check(d["flows"] and d["iterations"] > 0, "the dead-oil twin flows too")
check(d["q"] < v["q"],
      f"without gas the operating rate falls ({d['q']:.4g} < {v['q']:.4g} "
      f"sm3/s): gas lifts the well")
check(d["pwf"] > v["pwf"],
      "and the dead-oil twin needs more bottom-hole pressure at its own point")
q_vogel_d = vogel_rate(d["pwf"], pres, qmax)
check(abs(d["q"] / q_vogel_d - 1.0) < 1e-5,
      "the dead-oil twin's point satisfies the same Vogel identity")
check(d["holdup"] >= 0.999,
      "the dead-oil column is liquid full (GOR 0 leaves no free gas)")

# ---- the dead well: honesty about the static column ----
w = results["deadwell"]
check(not w["flows"] and w["q"] == 0.0 and w["iterations"] == 0,
      "the dead well reports rate 0 with flows=false and no secant work")
check(abs(w["pwf"] - w["static"]) < 1e-6 * w["static"],
      "the dead well's bottom hole is exactly the static column")

# ---- golden anchors ----
for name, gold in GOLDEN.items():
    for key, ref in gold.items():
        got = results[name][key]
        check(abs(got / ref - 1.0) < 1e-9,
              f"{name} {key} {got:.12g} holds the golden anchor {ref:.12g}")

print("\n".join(notes))
if failures:
    print("\nFAILURES:")
    print("\n".join(failures))
    raise SystemExit(1)
print(f"all {sum(1 for _ in notes)} nodal checks passed")
