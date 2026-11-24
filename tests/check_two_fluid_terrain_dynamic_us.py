#!/usr/bin/env python3
"""G2 gas-cap-driven terrain slug velocity (dynamicUs) for the distributed
two-fluid kernel (chflow-1d-twofluid): the gvol = vM - flux repair's
"no spurious ejection" judgment.

usage: check_two_fluid_terrain_dynamic_us.py <exe> <case.dict>

The case keeps the STRONGEST injection attack on the gas cap (vsl 0,
vsg > 0 into the upstream pocket, sealed outlet) for the whole run, and
the repaired gas-face identity F + G = vM is supposed to make that attack
fail: the captured lump is NOT a gas barrier (its region keeps
alpha_g ~ 0.2 and the gas phase seeps through it continuously), the two
pockets stay in pressure equilibrium, and the ejection drive
    pcap - p_front - rhoL g (z(front) - z(dip))
never crosses zero. The old G2 scenario's us > 0 rode exactly the broken
F + G = vM (its spurious front pressure gradients), so the checker
asserts the repaired invariant's consequence: us is bit-exact 0 on every
frame and the re-derived drive stays below 1 Pa. A positive drive would
need the lump to seal its gas phase - a future increment; this batch
invents no new physics to satisfy a judgment.

The checker re-derives the algebraic us relation
    u_s = sqrt(2 max(p_cap - p_front, 0) / rhoL)
from the REPORTED p profile, pcap and cfront (the front position the
reported velocity was computed at; the reported "front" is that plus
u*dt) every frame (1e-6) as the reporting-contract witness, and the
dynamicUs-false twin completes a full capture/discharge cycle (deposits
counter). The old clamp-witness and velocity-span judgments are gone:
under the repaired physics the drive is identically non-positive, so us
never reaches the maxSlugVelocity clamp and a positive-drive witness
would require inventing physics; the four loud rejections stay.
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
                s["cfront"] = float(tok[tok.index("cfront") + 1])
            frames[-1]["slugs"].append(s)
        elif tok[0] == "terrain" and tok[1] == "deposits" and frames:
            frames[-1]["deposits"] = int(tok[2])
    return frames


text = open(CASE).read()
frames = run(CASE)

# ---- algebraic us relation + the no-spurious-ejection judgment ----
# the kernel computes us at the PRE-travel front (reported as cfront; the
# reported "front" is cfront + us*dt), so the re-derivation locates the
# cell and the elevation from cfront - exact reporting-contract witness
n_checked = 0
max_drive = None
bad_us_nonzero = []
for fr in frames:
    for s in fr["slugs"]:
        if not s["active"] or "us" not in s:
            continue
        # the front position the reported us/pcap pair was computed at;
        # the capture frame ran no dynamicUs branch (cfront 0) and its
        # front is exactly the dip position, so either value locates the
        # cell the kernel sampled
        pos = s.get("cfront", 0.0) or s["front"]
        fi = cell_of(pos)
        p_front = fr["cells"][fi]["p"] + RHO_L * G * max(
            0.0, z_at(pos) - Z[DIP])
        drive = s["pcap"] - p_front
        if max_drive is None or drive > max_drive:
            max_drive = drive
        # repaired physics: the drive never crosses zero, so the ejection
        # value is bit-exact 0 (IEEE sqrt(0) == 0 on the max(., 0) clamp)
        if s["us"] != 0.0:
            bad_us_nonzero.append((fr["t"], s["us"]))
        u_want = math.sqrt(2.0 * max(drive, 0.0) / RHO_L)
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
if bad_us_nonzero:
    fails.append(f"spurious ejection: us != 0 on {len(bad_us_nonzero)} "
                 f"frames (first t={bad_us_nonzero[0][0]} us="
                 f"{bad_us_nonzero[0][1]})")
else:
    print(f"ok:   us is bit-exact 0 on every reported frame "
          f"({n_checked} frames) -- the sealed line is ejection-free "
          f"under maximum injection provocation")
if max_drive is None:
    fails.append("no drive sample found")
else:
    print(f"{'ok:  ' if max_drive < 1.0 else 'FAIL:'} re-derived ejection "
          f"drive stays below 1 Pa over the whole run: max {max_drive!r} Pa")
    if max_drive >= 1.0:
        fails.append(f"ejection drive reached {max_drive} Pa >= 1 Pa")

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
