#!/usr/bin/env python3
"""Analytic verification of the Weinaug-Katz compositional surface tension.

usage: check_pvt_surface_tension.py <exe> <case.dict>

D7 increment 2 (5.94). The Weinaug & Katz (1943) parachor surface tension
  sigma[N/m] = [ sum_i P_i (x_i/Vml - y_i/Vmg) ]^4
is checked two independent ways:

  formula    a pure-Python transcription reproduces the two worked examples
             published by the `chemicals` library (Weinaug_Katz) to 1e-12 --
             this pins the SI parachor units, the molar-volume form, the fourth
             power and the vapour-neglected (Vmg -> large) limit.
  kernel     for every two-phase flash the checker re-derives sigma from the
             DECLARED parachors and the reported PR phase state (Vml = ML/rhoL,
             Vmg = MV/rhoV, x, y) and weighs it against the C++ report.
  physics    single-phase flashes have no interface (sigma == 0); two-phase
             sigma > 0; sigma falls as the phases converge toward the bubble
             point; neglecting the vapour density can only raise sigma.
  gate       removing the `surfaceTension` block drops the report and leaves
             every other pvt line bit-identical (untouched path).

The parachors are declared test inputs, not a tabulated dataset: this increment
transcribes only the verified equation.
"""
import os
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]


def weinaug_katz(Ps, Vml, Vmg, xs, ys):
    rhoL = 1.0 / Vml
    rhoV = 1.0 / Vmg
    tot = sum(P * (x * rhoL - y * rhoV) for P, x, y in zip(Ps, xs, ys))
    return tot ** 4


fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


# ---- formula self-test against the published worked examples ----
close("Weinaug-Katz reproduces the chemicals published example",
      weinaug_katz([5.1e-5, 7.2e-5], 0.000125, 0.02011, [.4, .6], [.6, .4]),
      0.06547479150776776, 0.0, 1.0e-12)
close("Weinaug-Katz reproduces the vapour-neglected example",
      weinaug_katz([5.1e-5, 7.2e-5], 0.000125, 1.0e100, [.4, .6], [.6, .4]),
      0.06701752894095361, 0.0, 1.0e-12)


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"pvt failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


with open(case) as fh:
    src = fh.read()
mp = re.search(r"parachors\s*\(([^)]*)\)", src)
if not mp:
    raise SystemExit("case declares no surfaceTension parachors")
parachors = [float(v) for v in mp.group(1).split()]


def parse(out):
    flashes, cur = [], None
    for line in out.splitlines():
        t = line.replace("[chflow]", "").split()
        if not t or t[0] != "pvt":
            continue
        k = t[1]
        if k == "flash":
            kv = dict(zip(t[2::2], t[3::2]))
            cur = {"T": float(kv["T"]), "p": float(kv["p"]),
                   "beta": float(kv["beta"]), "single": kv["single"] == "1",
                   "x": {}, "y": {}}
            flashes.append(cur)
        elif k == "liquid" and cur is not None:
            cur["rhoL"], cur["ML"] = float(t[3]), float(t[5])
        elif k == "vapour" and cur is not None:
            cur["rhoV"], cur["MV"] = float(t[3]), float(t[5])
        elif k in ("x", "y") and cur is not None:
            cur[k][t[2]] = float(t[3])
        elif k == "surfaceTension" and cur is not None:
            cur["sigma"] = float(t[2])
    return flashes


out = run(case)
flashes = parse(out)
check("case produced flashes", len(flashes) >= 3, f"{len(flashes)} flashes")
check("case reports surfaceTension", all("sigma" in f for f in flashes),
      "no pvt surfaceTension lines (is the surfaceTension block present?)")

two_phase = []
for f in flashes:
    label = f"p={f['p']:.6g}"
    if "sigma" not in f:
        continue
    if f["single"]:
        check(f"{label} single-phase sigma is zero", f["sigma"] == 0.0, f"{f['sigma']}")
        continue
    names = list(f["x"].keys())
    check(f"{label} one parachor per component", len(parachors) == len(names),
          f"{len(parachors)} parachors vs {len(names)} components")
    xs = [f["x"][n] for n in names]
    ys = [f["y"][n] for n in names]
    Vml = f["ML"] / f["rhoL"]
    Vmg = f["MV"] / f["rhoV"]
    ref = weinaug_katz(parachors, Vml, Vmg, xs, ys)
    close(f"{label} sigma vs independent Weinaug-Katz", f["sigma"], ref, 1.0e-9, 1.0e-15)
    check(f"{label} two-phase sigma positive", f["sigma"] > 0.0, f"{f['sigma']}")
    neg = weinaug_katz(parachors, Vml, 1.0e100, xs, ys)
    check(f"{label} vapour-neglected sigma >= full", neg >= f["sigma"] - 1.0e-15,
          f"neg {neg:.6g} full {f['sigma']:.6g}")
    two_phase.append((f["p"], f["sigma"]))

# ---- physics: sigma falls as the phases converge toward the bubble point ----
two_phase.sort()
for (p0, s0), (p1, s1) in zip(two_phase, two_phase[1:]):
    check(f"sigma decreases from p={p0:.6g} to p={p1:.6g}", s0 >= s1,
          f"{s0:.6g} vs {s1:.6g}")
check("case has a two-phase branch", len(two_phase) >= 2, f"{len(two_phase)} two-phase flashes")

# ---- gate: removing the surfaceTension block drops the report, rest unchanged ----
twin = re.sub(r"surfaceTension\s*\{[^}]*\}", "", src, count=1)
twin_path = case + ".nosurf"
with open(twin_path, "w") as fh:
    fh.write(twin)
twin_out = run(twin_path)
os.remove(twin_path)
check("untouched case reports no surfaceTension", "pvt surfaceTension" not in twin_out,
      "surfaceTension leaked into the gated-off report")
base_other = [l for l in out.splitlines() if "pvt surfaceTension" not in l]
twin_other = [l for l in twin_out.splitlines() if "pvt surfaceTension" not in l]
check("gated-off pvt report bit-identical", base_other == twin_other,
      "a non-surfaceTension pvt line changed when the block was removed")

if fails:
    print("\n".join(f"FAIL: {m}" for m in fails))
    raise SystemExit(f"{len(fails)} pvt-surface-tension check(s) failed")

print(f"pvt-surface-tension: {len(flashes)} flashes ({len(two_phase)} two-phase), "
      f"published examples matched, kernel == independent Weinaug-Katz reference, "
      f"gated report clean - all checks passed")