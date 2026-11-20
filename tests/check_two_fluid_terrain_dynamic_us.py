#!/usr/bin/env python3
"""G2 gas-cap-driven terrain slug velocity (dynamicUs) for the distributed
two-fluid kernel (chflow-1d-twofluid).

usage: check_two_fluid_terrain_dynamic_us.py <exe> <case.dict>

The captured lump travels at the energy-balance ejection value
    u_s = sqrt(2 max(p_cap - p_front, 0) / rhoL)
with p_cap relaxing (implicit Euler) towards the EOS pressure of the cell
at the slug rear and p_front the front cell's EOS pressure plus the
hydrostatic head rhoL g (z(front) - z(dip)) still to be lifted. The
checker re-derives that algebraic relation from the REPORTED p profile,
pcap and front every frame (1e-6), witnesses the maxSlugVelocity clamp,
the terrain-dependent velocity span and mean advance, the exact kinematic
twin completing a full capture/discharge cycle (deposits counter), the
dynamicUs-false kinematics, and the loud rejections.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASE = sys.argv[1:3]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def expect_fatal(name, path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


# geometry (mirrors the case's pipes block: 30/-10, 5/-60, 5/+60, 10/+5,
# one metre per section)
PIPES = [(30, -10.0), (5, -60.0), (5, 60.0), (10, 5.0)]
DX = 1.0
N = sum(n for n, _ in PIPES)
Z = [0.0]
for n, ang in PIPES:
    for _ in range(n):
        Z.append(Z[-1] + DX * math.sin(math.radians(ang)))
ARC = [0.0]
for n, _ in PIPES:
    for _ in range(n):
        ARC.append(ARC[-1] + DX)
TOTAL = ARC[-1]
DIP = min(range(1, N), key=lambda j: Z[j])  # interior low point (face 35)
RHO_L = 1000.0
G = 9.80665
US_GATE = 0.5
MAXU = 1.5
DT = 0.01
STEPS_PER_FRAME = 100


def z_at(pos):
    if pos <= 0:
        return Z[0]
    if pos >= TOTAL:
        return Z[N]
    i = int(pos / DX)
    if i >= N:
        i = N - 1
    f = (pos - ARC[i]) / DX
    return Z[i] + f * (Z[i + 1] - Z[i])


def cell_of(pos):
    i = int(pos / DX)
    return min(max(i, 0), N - 1)


def run(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    frames = []
    for line in r.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "cells": [], "slugs": []})
        elif tok[0] == "cell" and frames:
            c, key = {}, None
            for w in tok[2:]:
                if key is None:
                    key = w
                else:
                    c[key] = float(w)
                    key = None
            frames[-1]["cells"].append(c)
        elif tok[0] == "terrain" and tok[1] == "slug" and frames:
            s = {"front": float(tok[tok.index("front") + 1]),
                 "vol": float(tok[tok.index("vol") + 1]),
                 "active": "ACTIVE" in tok}
            if "us" in tok:
                s["us"] = float(tok[tok.index("us") + 1])
                s["pcap"] = float(tok[tok.index("pcap") + 1])
            frames[-1]["slugs"].append(s)
        elif tok[0] == "terrain" and tok[1] == "deposits" and frames:
            frames[-1]["deposits"] = int(tok[2])
    return frames


text = open(CASE).read()
frames = run(CASE)

# ---- algebraic us relation every active frame ----
n_checked = 0
us_seen = []
capture_t = None
first_front = None
last_front = None
last_t = None
for fr in frames:
    for s in fr["slugs"]:
        if not s["active"]:
            continue
        if capture_t is None:
            capture_t = fr["t"]
            first_front = s["front"]
        last_front = s["front"]
        last_t = fr["t"]
        if "us" not in s:
            continue
        us_seen.append(s["us"])
        fi = cell_of(s["front"])
        p_front = fr["cells"][fi]["p"] + RHO_L * G * max(
            0.0, z_at(s["front"]) - Z[DIP])
        u_want = math.sqrt(2.0 * max(s["pcap"] - p_front, 0.0) / RHO_L)
        if MAXU > 0:
            u_want = min(u_want, MAXU)
        if abs(s["us"] - u_want) > 1e-6 * max(u_want, 1.0):
            fails.append(f"us relation broken at t={fr['t']}: "
                         f"got {s['us']} want {u_want}")
        n_checked += 1
if n_checked == 0:
    fails.append("no active dynamic slug frame found")
else:
    print(f"ok:   u_s = sqrt(2 max(p_cap - p_front, 0)/rho_L) re-derived "
          f"from the reported p profile on {n_checked} frames (1e-6)")

# ---- clamp + terrain-dependent dynamics ----
if us_seen:
    check("dynamic us never exceeds maxSlugVelocity",
          max(us_seen), MAXU, 0, 1e-9)
    span = max(us_seen) - min(us_seen)
    print(f"{'ok:  ' if span > 1e-3 else 'FAIL:'} travel velocity varies "
          f"over the run (span {span}) -- not the constant-us crawl")
    if span <= 1e-3:
        fails.append("dynamic us constant")
    if capture_t is None or last_front is None or last_front <= first_front:
        fails.append("slug never travelled")
    else:
        rate = (last_front - first_front) / (last_t - capture_t)
        print(f"{'ok:  ' if abs(rate - US_GATE) > 0.2 * US_GATE else 'FAIL:'} "
              f"mean advance {rate:.4f} m/s vs kinematic {US_GATE} m/s "
              f"(terrain/EOS-dependent travel)")
        if abs(rate - US_GATE) <= 0.2 * US_GATE:
            fails.append("dynamic travel rate equals the kinematic rate")

# ---- dynamicUs false: the exact kinematic twin ----
with tempfile.TemporaryDirectory(prefix="cf_tus_") as td:
    p = os.path.join(td, "twin.dict")
    open(p, "w").write(text.replace("dynamicUs true", "dynamicUs false"))
    fr_t = run(p)
    bad = []
    # kinematic relation: consecutive ACTIVE frames of the SAME lump
    # advance by us*dt*steps; a front that falls back is a re-capture
    # (the previous lump discharged) and starts a new witness
    prev = None
    for fr in fr_t:
        for s in fr["slugs"]:
            if s["active"]:
                if prev is not None and s["front"] >= prev[1]:
                    want = prev[1] + US_GATE * (fr["t"] - prev[0])
                    if abs(s["front"] - want) > 1e-6:
                        bad.append((fr["t"], s["front"], want))
                prev = (fr["t"], s["front"])
    check("dynamicUs false keeps the exact constant-us kinematics",
          float(not bad), 1, 0)
    # the full capture/travel/discharge cycle is witnessed on the twin:
    # the deposits counter increments when a lump discharges
    twin_dep = max(fr.get("deposits", 0) for fr in fr_t)
    print(f"{'ok:  ' if twin_dep >= 1 else 'FAIL:'} kinematic twin completes "
          f"a capture/discharge cycle (deposits {twin_dep})")
    if twin_dep < 1:
        fails.append("kinematic twin never discharged a lump")

    # ---- rejections ----
    p = os.path.join(td, "noreos.dict")
    open(p, "w").write(re.sub(r"\n\s*gasR 287;", "", text))
    expect_fatal("dynamicUs without the EOS gasR", p)
    p = os.path.join(td, "notau.dict")
    open(p, "w").write(re.sub(r"\n\s*gasCapTime 30;", "", text))
    expect_fatal("dynamicUs without gasCapTime", p)
    p = os.path.join(td, "badword.dict")
    open(p, "w").write(text.replace("dynamicUs true", "dynamicUs maybe"))
    expect_fatal("dynamicUs word outside true/false", p)
    p = os.path.join(td, "nogate.dict")
    open(p, "w").write(text.replace("slugVelocity 0.5", "slugVelocity 0"))
    expect_fatal("dynamicUs without the slugVelocity gate", p)

print()
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("all dynamicUs judgments passed")
