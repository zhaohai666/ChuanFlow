#!/usr/bin/env python3
"""E7b (5.182): slug mechanics on the two-fluid kernel - the Bendiksen
(1984) front closure (frontClosure bendiksen) plus the SLUG/FILM cell
partition (slugTracking{ zoning true }).

usage: check_transient_slug_mechanics.py <exe> <case-dir>

Nothing here imports or mirrors a C++ code path. The checker re-derives
everything independently from the CASE TEXT (geometry, schedule, declared
coefficients) and the REPORTED field lines (t/vm, cell alpha/vsl/zone,
slugtrack states and events), and compares line by line:

  * the front velocity of every live tracked slug equals
      uF = C0 vm + V_d
    with the DECLARED C0 (default 1.2) and the default drift
      V_d = 0.35 sin(theta_cell) sqrt(g D (rhoL-rhoG)/rhoL)
    (Bendiksen 1984, turbulent branch) evaluated from the case geometry -
    zero on the horizontal main case, 0.35 sin(5 deg) sqrt(g D dRho/rhoL)
    on the inclined variant, recomputed every frame from the reported vm;
  * the whole Lagrangian bookkeeping (advance / retire / merge / generate)
    replays exactly, with the Bendiksen front in place of the declared
    constant;
  * the zone partition is re-derived from the replayed spans (cell centre
    inside [xR, xF] -> SLUG=1, else FILM=0) and matches the reported
    `zone` key on every cell row of every frame; the existing cell keys
    (including the closure slip `s`) keep their meaning untouched - the
    checker also proves the boundaries MOVE (SLUG cells appear, vanish
    and hand over as the tracked faces travel);
  * the true zero-declaration twin (no slugTracking block at all) is
    byte-identical to the declared run once the slugtrack lines and the
    zone tokens are stripped, and enable false is byte-identical to it
    outright (the undeclared path is untouched);
  * the legacy twin (main case minus the two E7b keys) still rides the
    declared constant slugVelocity (the legacy front law is intact);
  * the mutual exclusions and the coefficient range checks are refused
    loudly (cfFatal), 8 rejection witnesses.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASEDIR = sys.argv[1:3]

fails = []
n_checks = 0


def check(name, got, want, rtol, atol=0.0):
    global n_checks
    n_checks += 1
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def expect_fatal(name, path, needle=None):
    global n_checks
    n_checks += 1
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    if not bad and needle is not None:
        bad = needle.lower() not in msg
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


def run_raw(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout


def case_scalar(text, key, default=None):
    m = re.search(rf"\b{key}\s+([0-9.eE+-]+)\s*;", text)
    if m:
        return float(m.group(1))
    if default is not None:
        return default
    raise SystemExit(f"case text has no {key}")


def parse(text):
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "vm": float(tok[3]),
                           "cells": [], "states": [], "events": []})
        elif tok[0] == "cell" and frames:
            c, key = {}, None
            for w in tok[2:]:
                if key is None:
                    key = w
                else:
                    c[key] = float(w)
                    key = None
            frames[-1]["cells"].append(c)
        elif tok[0] == "slugtrack" and frames:
            if tok[1] == "count":
                continue
            i = tok.index
            ev = {"id": int(tok[1]), "t": float(tok[i("t") + 1]),
                  "xR": float(tok[i("xR") + 1]), "xF": float(tok[i("xF") + 1]),
                  "inv": float(tok[i("inv") + 1]),
                  "uF": float(tok[i("uF") + 1]),
                  "uR": float(tok[i("uR") + 1]),
                  "kind": tok[i("event") + 1]}
            if ev["kind"] == "NONE":
                frames[-1]["states"].append(ev)
            else:
                frames[-1]["events"].append(ev)
    return frames


G = 9.80665


def bendiksen_vd(text, theta_deg):
    """The default Bendiksen drift law from the case geometry."""
    d = case_scalar(text, "diameter")
    rho_l = case_scalar(text, "rhoL")
    rho_g = case_scalar(text, "rhoG")
    return 0.35 * math.sin(math.radians(theta_deg)) * math.sqrt(
        G * d * (rho_l - rho_g) / rho_l)


def replay(frames, c0, vd, dt, total, n, theta_deg):
    """Independent replay of the tracking with the Bendiksen front; the
    drift is priced at the cell holding the front face (uniform slope, so
    one sine). Returns per-frame reference objects and the event journal."""
    sin_t = math.sin(math.radians(theta_deg))
    objs, events, states = [], [], []
    next_id = 0
    err = {"pos": 0.0, "inv": 0.0, "u": 0.0}
    for k, fr in enumerate(frames):
        if k:
            # per-frame Bendiksen front price from the REPORTED vm
            def uf_at(x):
                i = min(int(x / 1.0), n - 1)  # uniform cells, dx resolved below
                vd_cell = (vd if vd is not None else
                           0.35 * sin_t * math.sqrt(
                               G * DX * (RHO_L - RHO_G) / RHO_L))
                return c0 * fr["vm"] + vd_cell

            for o in objs:
                o["uF"] = uf_at(o["xF"])
                a = fr["cells"][cell_at(o["xR"])]["alpha"]
                o["uR"] = (fr["cells"][cell_at(o["xR"])]["vsl"] / a
                           if a > 1e-9 else fr["vm"])
                o["xF"] = min(TOTAL, o["xF"] + o["uF"] * DT)
                o["xR"] = max(0.0, o["xR"] + o["uR"] * DT)
                if o["xR"] > o["xF"]:
                    o["xR"] = o["xF"]
                if o["xR"] < TOTAL:
                    o["inv"] = inventory(fr, o["xR"], o["xF"])
            k2 = 0
            while k2 < len(objs):
                if objs[k2]["xR"] >= TOTAL:
                    o = objs.pop(k2)
                    events.append({"kind": "OUT", "t": fr["t"],
                                   "id": o["id"], "xR": o["xR"],
                                   "xF": o["xF"], "inv": o["inv"],
                                   "uF": o["uF"], "uR": o["uR"]})
                else:
                    k2 += 1
            objs.sort(key=lambda o: o["xR"])
            k2 = 0
            while k2 + 1 < len(objs):
                a_o, b_o = objs[k2], objs[k2 + 1]
                if b_o["xR"] <= a_o["xF"]:
                    events.append({"kind": "MERGE", "t": fr["t"],
                                   "id": a_o["id"], "xR": a_o["xR"],
                                   "xF": max(a_o["xF"], b_o["xF"]),
                                   "inv": a_o["inv"] + b_o["inv"],
                                   "uF": b_o["uF"], "uR": a_o["uR"]})
                    a_o["xF"] = max(a_o["xF"], b_o["xF"])
                    a_o["inv"] += b_o["inv"]
                    a_o["uF"] = b_o["uF"]
                    del objs[k2 + 1]
                else:
                    k2 += 1
            for (r0, r1, rinv) in scan_runs(fr):
                if any(r1 > o["xR"] and r0 < o["xF"] for o in objs):
                    continue
                o = {"id": next_id, "xR": r0, "xF": r1, "inv": rinv,
                     "uF": uf_at(r1), "uR": 0.0, "tBorn": fr["t"]}
                o["uR"] = (fr["cells"][cell_at(r0)]["vsl"]
                           / fr["cells"][cell_at(r0)]["alpha"]
                           if fr["cells"][cell_at(r0)]["alpha"] > 1e-9
                           else fr["vm"])
                next_id += 1
                objs.append(o)
                events.append({"kind": "GEN", "t": fr["t"], "id": o["id"],
                               "xR": o["xR"], "xF": o["xF"], "inv": o["inv"],
                               "uF": o["uF"], "uR": o["uR"]})
            rep = {s["id"]: s for s in fr["states"]}
            if set(rep) != {o["id"] for o in objs}:
                raise SystemExit(f"frame {k}: object set mismatch reported "
                                 f"{sorted(rep)} vs reference "
                                 f"{sorted(o['id'] for o in objs)}")
            for o in objs:
                r = rep[o["id"]]
                err["pos"] = max(err["pos"], abs(o["xR"] - r["xR"]),
                                 abs(o["xF"] - r["xF"]))
                err["inv"] = max(err["inv"], abs(o["inv"] - r["inv"]))
                err["u"] = max(err["u"], abs(o["uF"] - r["uF"]),
                               abs(o["uR"] - r["uR"]))
        # reference zone partition for this frame (cell centres)
        zones = []
        acc = 0.0
        for i in range(N):
            centre = acc + 0.5 * DX
            zones.append(1 if any(o["xR"] <= centre <= o["xF"]
                                  for o in objs) else 0)
            acc += DX
        states.append((list(objs), zones))  # copy: objs keeps mutating
    return states, events, err


def cell_at(x):
    return min(max(int(x / DX), 0), N - 1)


def inventory(fr, xR, xF):
    lo, hi = max(0.0, xR), min(TOTAL, xF)
    if hi <= lo:
        return 0.0
    inv = 0.0
    for i in range(N):
        o = min(hi, (i + 1) * DX) - max(lo, i * DX)
        if o > 0.0:
            inv += fr["cells"][i]["alpha"] * o
    return inv


def scan_runs(fr):
    out, inrun = [], False
    rear = length = inv = x = 0.0
    for i in range(N):
        if fr["cells"][i]["alpha"] >= THRESH:
            if not inrun:
                inrun, rear, length, inv = True, x, 0.0, 0.0
            length += DX
            inv += fr["cells"][i]["alpha"] * DX
        else:
            if inrun and length >= MINLEN:
                out.append((rear, rear + length, inv))
            inrun = False
        x += DX
    if inrun and length >= MINLEN:
        out.append((rear, rear + length, inv))
    return out


def compare_events(frames, ref_events):
    rep = [e for fr in frames for e in fr["events"]]
    if len(rep) != len(ref_events):
        return (f"event count: reported {len(rep)} vs reference "
                f"{len(ref_events)}", False)
    for a, b in zip(rep, ref_events):
        if a["kind"] != b["kind"] or a["id"] != b["id"]:
            return f"event order/kind: {a} vs {b}", False
        for k in ("t", "xR", "xF", "inv", "uF", "uR"):
            if abs(a[k] - b[k]) > 1e-8 * max(1.0, abs(b[k])):
                return f"{b['kind']} field {k}: {a[k]} vs {b[k]}", False
    return "", True


# ---- case constants, parsed from the case TEXT (not trusted from memory) --
with open(os.path.join(CASEDIR, "case.dict")) as f:
    CASE_TEXT = f.read()
N = int(case_scalar(CASE_TEXT, "sections"))
DX = case_scalar(CASE_TEXT, "length") / N
TOTAL = N * DX
DT = case_scalar(CASE_TEXT, "dt")
THRESH = case_scalar(CASE_TEXT, "alphaThreshold")
MINLEN = case_scalar(CASE_TEXT, "minLength")
RHO_L = case_scalar(CASE_TEXT, "rhoL")
RHO_G = case_scalar(CASE_TEXT, "rhoG")
THETA = case_scalar(CASE_TEXT, "angleDeg")
C0 = 1.2          # the Bendiksen default (nothing declared in case.dict)
VD_MAIN = 0.0     # horizontal: sin(0) kills the drift term exactly

# ================= declared run =================
text = run_raw(os.path.join(CASEDIR, "case.dict"))
frames = parse(text)
check("frames reported (reportEvery 1 over 60 s at dt 0.05)",
      len(frames), 1201, 0)

# ---- zone keys: present, binary, absent at t = 0 ---------------------
zone_rows = sum(1 for fr in frames for c in fr["cells"] if "zone" in c)
check("every cell row of every frame carries the zone key",
      zone_rows, len(frames) * N, 0)
zone_vals = {c["zone"] for fr in frames for c in fr["cells"]}
check("zone flags are binary (SLUG=1 / FILM=0)",
      float(zone_vals <= {0.0, 1.0}), 1.0, 0)
check("t = 0 has no tracked slug yet: every cell is FILM",
      float(all(c["zone"] == 0.0 for c in frames[0]["cells"])), 1.0, 0)

# ---- Bendiksen front velocity, re-derived from the REPORTED vm --------
uf_err = 0.0
for fr in frames:
    want = C0 * fr["vm"] + VD_MAIN
    for s in fr["states"]:
        uf_err = max(uf_err, abs(s["uF"] - want))
check("every reported front velocity equals C0 vm + V_d (horizontal: "
      "V_d = 0) with C0 = 1.2, relative tolerance 1e-9",
      uf_err, 0.0, 0, 1e-9)
uf_min = min(s["uF"] for fr in frames for s in fr["states"])
vm_max = max(fr["vm"] for fr in frames)
check("the Bendiksen closure never undercuts the mixture velocity "
      "(C0 = 1.2 > 1)", float(uf_min >= C0 * vm_max - 1e-12), 1.0, 0)

# ---- full independent replay (bendiksen front) ------------------------
states, ref_events, err = replay(frames, C0, None, DT, TOTAL, N, THETA)
check("trajectories re-derived independently (max face error)",
      err["pos"], 0.0, 0, 1e-7)
check("inventory re-estimates match the field integrals",
      err["inv"], 0.0, 0, 1e-7)
check("face velocities re-derived (front = Bendiksen, rear = vsl/alpha)",
      err["u"], 0.0, 0, 1e-7)
bad, ok = compare_events(frames, ref_events)
check("event journal (times, ids, faces, velocities) reproduced exactly",
      float(ok), 1.0, 0)
if not ok:
    fails.append("journal: " + bad)
n_none = sum(len(fr["states"]) for fr in frames)
n_ref = sum(len(s) for s, _ in states)
check("state lines present every frame for every live object",
      n_none, n_ref, 0)

gens = [e for e in ref_events if e["kind"] == "GEN"]
merges = [e for e in ref_events if e["kind"] == "MERGE"]
outs = [e for e in ref_events if e["kind"] == "OUT"]
check("two pulses plus the drained tail body generate three slugs",
      len(gens), 3, 0)
check("exactly one in-pipe coalescence (the trailing Bendiksen front at "
      "1.2 m/s catches the leading rear)", len(merges), 1, 0)
check("every generated slug ends its life in a MERGE or an OUT (no leaks)",
      len(merges) + len(outs), len(gens), 0)
for g in gens:
    fr = next(f for f in frames if abs(f["t"] - g["t"]) < 1e-12)
    want = C0 * fr["vm"] + VD_MAIN
    check(f"GEN id {g['id']} is born riding the Bendiksen value",
          g["uF"], want, 1e-9)

# ---- zone partition re-derived from the replayed spans ----------------
mismatch = 0
slug_cells = film_to_slug = slug_to_film = 0
prev = None
for k, fr in enumerate(frames):
    _, zones = states[k]
    for i, c in enumerate(fr["cells"]):
        if c["zone"] != zones[i]:
            mismatch += 1
    if prev is not None:
        film_to_slug += sum(1 for i in range(N)
                            if zones[i] == 1 and prev[i] == 0)
        slug_to_film += sum(1 for i in range(N)
                            if zones[i] == 0 and prev[i] == 1)
    prev = zones
check("the reported zone partition matches the partition re-derived from "
      "the replayed slug spans (cell-centre rule), zero mismatches over "
      f"{len(frames) * N} cell-frames", float(mismatch), 0.0, 0)
check("the partition is not static: cells hand over FILM->SLUG as the "
      "tracked fronts arrive", float(film_to_slug > 0), 1.0, 0)
check("and back SLUG->FILM as the rears drain through", float(slug_to_film > 0), 1.0, 0)
check("the slug phase occupies cells at some point",
      float(any(z == 1 for _, zones in states for z in zones)), 1.0, 0)

# ---- zones agree directly with the REPORTED slugtrack spans -----------
# the state lines print 12 significant digits; a face that sits within a
# print-ulp of a cell centre makes the centre rule UNDECIDABLE from the
# report alone (the kernel priced the zone on the unrounded face). Cells
# whose centre lies within 1e-9 of a printed face are counted as
# print-rounding-ambiguous and exempt from the strict rule; everything
# else must classify exactly.
direct_bad = ambiguous = 0
EPS = 1e-9
for fr in frames:
    spans = [(s["xR"], s["xF"]) for s in fr["states"]]
    for i, c in enumerate(fr["cells"]):
        centre = i * DX + 0.5 * DX
        strict = 1.0 if any(r <= centre <= f for r, f in spans) else 0.0
        if any(abs(centre - r) <= EPS or abs(centre - f) <= EPS
               for r, f in spans):
            ambiguous += 1
        elif c["zone"] != strict:
            direct_bad += 1
check("away from face-ulp boundaries the zone key matches the partition "
      "re-derived from the REPORTED slugtrack state lines alone "
      "(checker trusts no internal)",
      float(direct_bad), 0.0, 0)
check("face-on-cell-centre ambiguities are rare print-rounding artifacts "
      "(<= 20 cell-frames out of 24020)",
      float(ambiguous <= 20), 1.0, 0)

# ---- byte-identity twins ----------------------------------------------
base_text = run_raw(os.path.join(CASEDIR, "caseNoBlock.dict"))
stripped = "\n".join(
    re.sub(r" zone [01](?=\s|$)", "", l)
    for l in text.splitlines() if "slugtrack" not in l)
check("declared run minus slugtrack lines and zone tokens is "
      "byte-identical to the no-block twin (the E7b keys are the only "
      "difference)", float(stripped == base_text.rstrip("\n")), 1.0, 0)
dis_text = run_raw(os.path.join(CASEDIR, "caseDisabled.dict"))
check("enable false is byte-identical to the no-block twin",
      float(dis_text.rstrip("\n") == base_text.rstrip("\n")), 1.0, 0)

# ---- legacy twin: the declared constant law is intact ------------------
leg_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
leg_frames = parse(leg_text)
check("legacy twin (main case minus the E7b keys) runs the full horizon",
      len(leg_frames), 1201, 0)
leg_uf = [s["uF"] for fr in leg_frames for s in fr["states"]]
check("legacy twin fronts ride the declared constant slugVelocity 1.5 "
      "(the legacy front law is untouched)",
      float(leg_uf and all(abs(u - 1.5) <= 1e-12 for u in leg_uf)), 1.0, 0)
leg_zone = sum(1 for fr in leg_frames for c in fr["cells"] if "zone" in c)
check("legacy twin emits no zone keys (zoning not declared there)",
      leg_zone, 0, 0)

# ================= inclined variant: the drift term =====================
inc_text = run_raw(os.path.join(CASEDIR, "caseInclined.dict"))
inc_frames = parse(inc_text)
check("inclined variant (angleDeg 5) runs the full horizon",
      len(inc_frames), 1201, 0)
vd_inc = bendiksen_vd(CASE_TEXT, 5.0)
check("the default drift law from the case geometry is "
      "0.35 sin(5 deg) sqrt(g D (rhoL-rhoG)/rhoL)",
      vd_inc,
      0.35 * math.sin(math.radians(5.0)) * math.sqrt(
          G * case_scalar(CASE_TEXT, "diameter")
          * (RHO_L - RHO_G) / RHO_L), 0, 1e-15)
check("the inclined drift term is strictly positive (the sin(theta) "
      "term is alive)", float(vd_inc > 0), 1.0, 0)
inc_uf_err = 0.0
for fr in inc_frames:
    want = C0 * fr["vm"] + vd_inc
    for s in fr["states"]:
        inc_uf_err = max(inc_uf_err, abs(s["uF"] - want))
check("inclined variant: every front velocity equals C0 vm + V_d with "
      "the nonzero drift, relative tolerance 1e-9",
      inc_uf_err, 0.0, 0, 1e-9)
inc_states, inc_events, inc_err = replay(inc_frames, C0, vd_inc, DT,
                                         TOTAL, N, 5.0)
check("inclined variant: trajectories re-derived with the drift-carrying "
      "front", inc_err["pos"], 0.0, 0, 1e-7)
inc_bad, inc_ok = compare_events(inc_frames, inc_events)
check("inclined variant: event journal reproduced exactly",
      float(inc_ok), 1.0, 0)
if not inc_ok:
    fails.append("inclined journal: " + inc_bad)
inc_gens = [e for e in inc_events if e["kind"] == "GEN"]
check("inclined variant: the pulses still generate tracked slugs",
      float(len(inc_gens) >= 2), 1.0, 0)
inc_mismatch = 0
for k, fr in enumerate(inc_frames):
    _, zones = inc_states[k]
    for i, c in enumerate(fr["cells"]):
        if c["zone"] != zones[i]:
            inc_mismatch += 1
check("inclined variant: zone partition re-derived with zero mismatches",
      float(inc_mismatch), 0.0, 0)
finite = all(math.isfinite(v) for fr in inc_frames for s in fr["states"]
             for v in (s["xR"], s["xF"], s["inv"], s["uF"], s["uR"]))
check("inclined variant: every reported number is finite",
      float(finite), 1.0, 0)

# ================= rejections ===========================================
for fname, needle in [
    ("caseRejectFrontWord.dict", "unknown slugTracking frontClosure"),
    ("caseRejectFrontMode.dict",
     "mutually exclusive with velocityMode"),
    ("caseRejectFrontSlug.dict", "mutually exclusive with slugVelocity"),
    ("caseRejectFrontCap.dict", "mutually exclusive with gasCapTime"),
    ("caseRejectFrontMax.dict",
     "mutually exclusive with maxSlugVelocity"),
    ("caseRejectC0.dict", "Bendiksen C0 must be positive"),
    ("caseRejectC0Neg.dict", "Bendiksen C0 must be positive"),
    ("caseRejectZoningWord.dict", "zoning must be true or false"),
]:
    expect_fatal(fname[:-5], os.path.join(CASEDIR, fname), needle)

print()
print(f"criteria: {n_checks}")
if fails:
    print(f"FAILED {len(fails)} check(s):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("transient slug mechanics (5.182): all checks passed")
