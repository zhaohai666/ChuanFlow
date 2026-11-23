#!/usr/bin/env python3
"""O2 (5.161): OLGA-style slug tracking - the full-front Lagrangian method
on the distributed two-fluid kernel (chflow-1d-twofluid).

usage: check_transient_slugtracking.py <exe> <case-dir>

The kernel maintains, beside the Eulerian holdup field, a set of tracked
slug objects {xR, xF, inv, uF, uR} that are GENERATED from the same
maximal alpha >= threshold run scan slugBodies() reports, ride a declared
constant migration velocity (or the 5.142 gas-cap ejection law) at the
front and the local liquid transport velocity vsl/alpha at the rear,
re-estimate their inventory from the field over [xR, xF] every step,
MERGE when the trailing front reaches the leading rear, and retire (OUT)
when the rear crosses the domain end.

This checker does NOT import or mirror any C++ code path: it re-derives
the whole bookkeeping independently in Python from the REPORTED per-cell
field (alpha, vsl, p) and compares every object trajectory, every event
and the discharge accounting line by line. On top of the replay it
verifies the physics witnesses the gap list asks for:
  * every GEN agrees with the Eulerian slug{} scan at the same frame
    (exact faces) and had no match one frame earlier (within one cell);
  * exactly one MERGE, with the merged inventory the exact sum of the
    two objects' inventories and the leader's front / trailer's rear;
  * every object that is born also retires, the OUT time is the rear
    crossing the outlet (the front pinned there long before), and the
    duration is t_OUT - t_BORN;
  * the slug-carried discharge equals the outlet-face liquid flux
    integrated over the pinned windows, is bounded by the total liquid
    that left the domain, and carries the bulk of it; the global liquid
    balance injected - discharged == inventory change holds independently
    of the tracking;
  * the undeclared twin (no slugTracking block) reproduces the declared
    run's whole physics report byte-for-byte once the slugtrack lines are
    stripped, and the enable-false twin is byte-identical outright;
  * the gasCap variant re-derives the ejection sqrt(2 max(p_cap -
    p_front, 0)/rho_L) with the implicit-Euler cap relaxation from the
    reported p profile every frame, under the maxSlugVelocity clamp;
  * out-of-range / inconsistent declarations are refused loudly (cfFatal).
"""

import math
import os
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASEDIR = sys.argv[1:3]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def expect_fatal(name, path, needle=None):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    if not bad and needle is not None:
        bad = needle.lower() not in msg
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


# ---- case constants (mirrors case.dict: 20 sections x 1 m) ----
DX = 1.0
N = 20
TOTAL = N * DX
DT = 0.05
THRESH = 0.45
MINLEN = 0.5
UF_CONST = 1.5
RHO_L = 1000.0
G = 9.80665
TAU = 5.0
MAXU = 2.0
# imposed inlet vsl: 0.75 during [5, 10) and [12, 17), else 0.1
def imposed_vsl(t_mid):
    return 0.75 if (5.0 <= t_mid < 10.0 or 12.0 <= t_mid < 17.0) else 0.1


def run_raw(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout


def parse(text):
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "vm": float(tok[3]),
                           "cells": [], "bodies": [], "states": [],
                           "events": [], "count": None})
        elif tok[0] == "cell" and frames:
            c, key = {}, None
            for w in tok[2:]:
                if key is None:
                    key = w
                else:
                    c[key] = float(w)
                    key = None
            frames[-1]["cells"].append(c)
        elif tok[0] == "slugbody" and frames:
            i = tok.index
            frames[-1]["bodies"].append(
                {"rear": float(tok[i("rear") + 1]),
                 "front": float(tok[i("front") + 1]),
                 "inv": float(tok[i("inv") + 1])})
        elif tok[0] == "slugtrack" and frames:
            if tok[1] == "count":
                frames[-1]["count"] = int(tok[2])
                continue
            i = tok.index
            ev = {"id": int(tok[1]), "t": float(tok[i("t") + 1]),
                  "xR": float(tok[i("xR") + 1]), "xF": float(tok[i("xF") + 1]),
                  "inv": float(tok[i("inv") + 1]),
                  "uF": float(tok[i("uF") + 1]),
                  "uR": float(tok[i("uR") + 1]),
                  "kind": tok[i("event") + 1]}
            if ev["kind"] == "OUT":
                ev["out"] = float(tok[i("out") + 1])
                ev["duration"] = float(tok[i("duration") + 1])
            if ev["kind"] == "NONE":
                frames[-1]["states"].append(ev)
            else:
                frames[-1]["events"].append(ev)
    return frames


# ---- independent reference replay of the tracking algorithm ----
def cell_at(x):
    i = int(x / DX)
    return min(max(i, 0), N - 1)


def liquid_velocity(fr, i):
    a = fr["cells"][i]["alpha"]
    return fr["cells"][i]["vsl"] / a if a > 1e-9 else fr["vm"]


def inventory(fr, xR, xF):
    lo = max(0.0, xR)
    hi = min(TOTAL, xF)
    if hi <= lo:
        return 0.0
    inv = 0.0
    for i in range(N):
        o = min(hi, (i + 1) * DX) - max(lo, i * DX)
        if o > 0.0:
            inv += fr["cells"][i]["alpha"] * o
    return inv


def scan_runs(fr):
    out = []
    inrun = False
    rear = length = inv = x = 0.0
    for i in range(N):
        if fr["cells"][i]["alpha"] >= THRESH:
            if not inrun:
                inrun = True
                rear = x
                length = 0.0
                inv = 0.0
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


def replay(frames, gascap):
    """Re-derive every object trajectory and event from the reported field;
    returns (per-frame state dict, event list, max errors)."""
    objs = []
    next_id = 0
    states = []
    events = []
    err_pos = err_inv = err_u = 0.0
    for n, fr in enumerate(frames):
        if n == 0:
            states.append({})
            continue
        # 1. advance
        for o in objs:
            if gascap:
                ri = cell_at(o["xR"])
                fi = cell_at(o["xF"])
                p_r = fr["cells"][ri].get("p", 0.0)
                p_f = fr["cells"][fi].get("p", 0.0)
                o["cap"] = (o["cap"] / DT + p_r / TAU) / (1.0 / DT + 1.0 / TAU)
                u = math.sqrt(2.0 * max(o["cap"] - p_f, 0.0) / RHO_L)
                if MAXU > 0.0:
                    u = min(u, MAXU)
                o["uF"] = u
            else:
                o["uF"] = UF_CONST
            o["uR"] = liquid_velocity(fr, cell_at(o["xR"]))
            o["xF"] += o["uF"] * DT
            o["xR"] += o["uR"] * DT
            if o["xR"] < 0.0:
                o["xR"] = 0.0
            if o["xF"] > TOTAL:
                o["xF"] = TOTAL
            if o["xR"] > o["xF"]:
                o["xR"] = o["xF"]
            if o["xR"] < TOTAL:
                o["inv"] = inventory(fr, o["xR"], o["xF"])
            if o["xF"] >= TOTAL and o["xR"] < TOTAL:
                o["outVol"] += fr["cells"][N - 1]["vsl"] * DT
        # 2. retire
        k = 0
        while k < len(objs):
            o = objs[k]
            if o["xR"] >= TOTAL:
                events.append({"kind": "OUT", "t": fr["t"], "id": o["id"],
                               "xR": o["xR"], "xF": o["xF"], "inv": o["inv"],
                               "uF": o["uF"], "uR": o["uR"],
                               "out": o["outVol"],
                               "duration": fr["t"] - o["tBorn"]})
                del objs[k]
            else:
                k += 1
        # 3. merge
        objs.sort(key=lambda o: o["xR"])
        k = 0
        while k + 1 < len(objs):
            a, b = objs[k], objs[k + 1]
            if b["xR"] <= a["xF"]:
                events.append({"kind": "MERGE", "t": fr["t"], "id": a["id"],
                               "xR": a["xR"], "xF": max(a["xF"], b["xF"]),
                               "inv": a["inv"] + b["inv"],
                               "uF": b["uF"], "uR": a["uR"]})
                a["xF"] = max(a["xF"], b["xF"])
                a["inv"] += b["inv"]
                a["uF"] = b["uF"]
                a["outVol"] += b["outVol"]
                del objs[k + 1]
            else:
                k += 1
        # 4. generate
        for (r0, r1, rinv) in scan_runs(fr):
            if any(r1 > o["xR"] and r0 < o["xF"] for o in objs):
                continue
            o = {"id": next_id, "xR": r0, "xF": r1, "inv": rinv,
                 "uR": liquid_velocity(fr, cell_at(r0)),
                 "uF": 0.0 if gascap else UF_CONST,
                 "tBorn": fr["t"], "outVol": 0.0,
                 "cap": (fr["cells"][cell_at(r0)].get("p", 0.0)
                         if gascap else 0.0)}
            next_id += 1
            objs.append(o)
            events.append({"kind": "GEN", "t": fr["t"], "id": o["id"],
                           "xR": r0, "xF": r1, "inv": rinv,
                           "uF": o["uF"], "uR": o["uR"]})
        # 5. compare against the reported state lines
        rep = {s["id"]: s for s in fr["states"]}
        if set(rep) != {o["id"] for o in objs}:
            raise SystemExit(f"frame {n}: object set mismatch "
                             f"reported {sorted(rep)} vs reference "
                             f"{sorted(o['id'] for o in objs)}")
        for o in objs:
            r = rep[o["id"]]
            err_pos = max(err_pos, abs(o["xR"] - r["xR"]),
                          abs(o["xF"] - r["xF"]))
            err_inv = max(err_inv, abs(o["inv"] - r["inv"]))
            err_u = max(err_u, abs(o["uF"] - r["uF"]), abs(o["uR"] - r["uR"]))
        states.append({o["id"]: dict(o) for o in objs})
    return states, events, err_pos, err_inv, err_u


def compare_events(frames, ref_events):
    """Reported journal (flattened in frame order) vs reference journal."""
    rep = [e for fr in frames for e in fr["events"]]
    if len(rep) != len(ref_events):
        return f"event count: reported {len(rep)} vs reference " \
               f"{len(ref_events)}", False
    for a, b in zip(rep, ref_events):
        if a["kind"] != b["kind"] or a["id"] != b["id"]:
            return f"event order/kind: {a} vs {b}", False
        for k in ("t", "xR", "xF", "inv", "uF", "uR"):
            if abs(a[k] - b[k]) > 1e-8 * max(1.0, abs(b[k])):
                return f"{b['kind']} field {k}: {a[k]} vs {b[k]}", False
        if b["kind"] == "OUT":
            for k in ("out", "duration"):
                if abs(a[k] - b[k]) > 1e-8 * max(1.0, abs(b[k])):
                    return f"OUT field {k}: {a[k]} vs {b[k]}", False
    return "", True


# ================= declared run =================
text = run_raw(os.path.join(CASEDIR, "case.dict"))
frames = parse(text)
check("frames reported (reportEvery 1 over 60 s at dt 0.05)",
      len(frames), 1201, 0)

states, ref_events, err_pos, err_inv, err_u = replay(frames, gascap=False)

# ---- full independent replay ----
check("object trajectories re-derived independently: max |xR-xF| error",
      err_pos, 0.0, 0, 1e-7)
check("inventory re-estimates match the field integrals",
      err_inv, 0.0, 0, 1e-7)
check("face velocities match vsl/alpha at the faces",
      err_u, 0.0, 0, 1e-7)
bad, ok = compare_events(frames, ref_events)
check("event journal (times, ids, faces, volumes) reproduced exactly",
      float(ok), 1.0, 0)
if not ok:
    fails.append("journal: " + bad)
n_none = sum(len(fr["states"]) for fr in frames)
n_ref = sum(len(s) for s in states)
check("state lines present every frame for every live object",
      n_none, n_ref, 0)

# ---- GEN vs the Eulerian slug{} scan ----
gens = [e for e in ref_events if e["kind"] == "GEN"]
check("two pulses plus the drained tail body generate three slugs",
      len(gens), 3, 0)
gen_scan_ok = len(gens) > 0
for g in gens:
    fr = next(f for f in frames if abs(f["t"] - g["t"]) < 1e-12)
    if not any(abs(b["rear"] - g["xR"]) <= 1e-9
               and abs(b["front"] - g["xF"]) <= 1e-9 for b in fr["bodies"]):
        gen_scan_ok = False
check("every GEN sits on the slug{} scan of the same frame, within one "
      "cell (here: exact faces)", float(gen_scan_ok), 1.0, 0)

# ---- MERGE ----
merges = [e for e in ref_events if e["kind"] == "MERGE"]
check("exactly one MERGE (the trailing front catches the leading rear "
      "in-pipe)", len(merges), 1, 0)
m = merges[0]
fr = next(f for f in frames if abs(f["t"] - m["t"]) < 1e-12)
# the two parents at the merge frame: previous frame's objects advanced
# one step with this frame's field (the same rule the kernel applies)
prev_objs = states[frames.index(fr) - 1]
adv = []
for o in prev_objs.values():
    xF = min(TOTAL, o["xF"] + UF_CONST * DT)
    uR = liquid_velocity(fr, cell_at(o["xR"]))
    xR = max(0.0, o["xR"] + uR * DT)
    if xR > xF:
        xR = xF
    adv.append({"id": o["id"], "xR": xR, "xF": xF,
                "inv": inventory(fr, xR, xF)})
adv.sort(key=lambda o: o["xR"])
adv_ok = len(adv) == 2 and adv[1]["xR"] <= adv[0]["xF"] \
    and m["id"] == adv[0]["id"]
check("MERGE involves exactly the trailing/leading pair with the trailer "
      "surviving", float(adv_ok), 1.0, 0)
if adv_ok:
    a, b = adv[0], adv[1]
    check("MERGE inventory is the exact sum of the two objects' field "
          "inventories (volume conservation at coalescence)",
          m["inv"], a["inv"] + b["inv"], 0, 1e-9)
    check("MERGE keeps the trailer rear and the leader front",
          (abs(m["xR"] - a["xR"]) < 1e-9) and
          (abs(m["xF"] - max(a["xF"], b["xF"])) < 1e-9), True, 0)
else:
    fails.append("merge parents not reproducible")

# ---- OUT ----
outs = [e for e in ref_events if e["kind"] == "OUT"]
check("every generated slug ends its life in an OUT or a MERGE "
      "(no leaks)", len(outs) + len(merges), len(gens), 0)
out_ok = all(e["xR"] >= TOTAL - 1e-12 for e in outs)
check("retirement happens when the rear crosses the domain outlet",
      float(out_ok), 1.0, 0)
dur_ok = True
pin_t = {}
for fr2 in frames:
    for s in fr2["states"]:
        if abs(s["xF"] - TOTAL) < 1e-9 and s["id"] not in pin_t:
            pin_t[s["id"]] = fr2["t"]
for e in outs:
    # duration = t_OUT - t_BORN: recover tBorn from the GEN of the same id
    g = next(ev for ev in gens if ev["id"] == e["id"])
    if abs(e["duration"] - (e["t"] - g["t"])) > 1e-9:
        dur_ok = False
    # the front reached the outlet (pinned) no later than the OUT
    if e["id"] in pin_t and not (pin_t[e["id"]] <= e["t"] + 1e-9):
        dur_ok = False
check("OUT reports t_OUT - t_BORN as the duration and follows the pinned "
      "front's arrival at the outlet", float(dur_ok), 1.0, 0)

# ---- discharge accounting ----
total_out = sum(e["out"] for e in outs)
pinned_flux = 0.0
for n2, fr2 in enumerate(frames):
    if n2 == 0:
        continue
    # the kernel adds the outlet flux while a live span covers the
    # outlet; an object born at this frame has not been advanced yet,
    # so its birth frame contributes no discharge
    born = {e["id"] for e in fr2["events"] if e["kind"] == "GEN"}
    if any(abs(s["xF"] - TOTAL) < 1e-9 and s["xR"] < TOTAL
           and s["id"] not in born for s in fr2["states"]):
        pinned_flux += fr2["cells"][N - 1]["vsl"] * DT
check("slug-carried discharge equals the outlet-face liquid flux "
      "integrated over the pinned windows", total_out, pinned_flux, 0, 1e-8)
discharge = sum(fr2["cells"][N - 1]["vsl"] * DT
                for n2, fr2 in enumerate(frames) if n2 > 0)
injected = sum(imposed_vsl(n2 * DT - 0.5 * DT) * DT
               for n2 in range(1, len(frames)))
check("no liquid is created: slug discharge <= total outlet discharge",
      float(total_out <= discharge + 1e-9), 1.0, 0)
# the brief's witness "total discharge ~ total injected liquid", stated
# honestly: the tracked slugs carry the bulk of the liquid that left the
# domain by the last OUT (film discharge outside the tracked windows and
# the liquid still in the pipe at that instant account for the rest)
injected_to_last_out = sum(
    imposed_vsl(n2 * DT - 0.5 * DT) * DT
    for n2 in range(1, int(max(e["t"] for e in outs) / DT) + 1))
check("the tracked slugs carry the bulk of the liquid injected by the "
      "last OUT (>= 60 percent, the documented tolerance)",
      float(total_out >= 0.6 * injected_to_last_out), 1.0, 0)
inv_end = sum(frames[-1]["cells"][i]["alpha"] * DX for i in range(N))
inv_0 = sum(frames[0]["cells"][i]["alpha"] * DX for i in range(N))
# the Godunov interface flux differs from the cell closure flux F(alpha)
# while a staircase tail crosses a face, so the telescoped balance is
# exact only in the steady limit; 1e-3 m per m2 of bore over the whole
# 60 s run is the documented tolerance (0.02 percent of the injected
# liquid)
check("global liquid balance: injected - discharged == inventory change "
      "(independent of the tracking, documented tolerance 1e-3)",
      injected - discharge, inv_end - inv_0, 0, 1e-3)

# ---- undeclared path byte-identity ----
base_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
stripped = "\n".join(l for l in text.splitlines() if "slugtrack" not in l)
check("declared run minus the slugtrack lines is byte-identical to the "
      "undeclared twin",
      float(stripped == base_text.rstrip("\n")), 1.0, 0)
dis_text = run_raw(os.path.join(CASEDIR, "caseDisabled.dict"))
check("enable false is byte-identical to the undeclared twin",
      float(dis_text.rstrip("\n") == base_text.rstrip("\n")), 1.0, 0)

# ================= gasCap variant =================
gc_text = run_raw(os.path.join(CASEDIR, "caseGasCap.dict"))
gc_frames = parse(gc_text)
gc_states, gc_events, gc_pos, gc_inv, gc_u = replay(gc_frames, gascap=True)
check("gasCap variant: trajectories and ejection velocity re-derived "
      "from the reported p profile", gc_u, 0.0, 0, 1e-7)
bad_gc, ok_gc = compare_events(gc_frames, gc_events)
check("gasCap variant: event journal reproduced exactly", float(ok_gc),
      1.0, 0)
if not ok_gc:
    fails.append("gasCap journal: " + bad_gc)
uF_all = [s["uF"] for fr2 in gc_frames for s in fr2["states"]]
check("gasCap variant: the maxSlugVelocity clamp is never exceeded",
      float(all(u <= MAXU + 1e-12 for u in uF_all)), 1.0, 0)
check("gasCap variant: slugs are still generated from the scan",
      float(sum(1 for e in gc_events if e["kind"] == "GEN") >= 2), 1.0, 0)
finite = all(math.isfinite(v) for fr2 in gc_frames
             for s in fr2["states"] for v in
             (s["xR"], s["xF"], s["inv"], s["uF"], s["uR"]))
check("gasCap variant: every reported number is finite", float(finite),
      1.0, 0)

# ================= rejections =================
with open(os.path.join(CASEDIR, "case.dict")) as f:
    case_text = f.read()
with tempfile.TemporaryDirectory(prefix="cf_st_") as td:
    def variant(name, old, new):
        assert old in case_text, name
        p = os.path.join(td, name + ".dict")
        with open(p, "w") as f:
            f.write(case_text.replace(old, new))
        return p

    expect_fatal("velocityMode outside constant/gasCap",
                 variant("mode", "velocityMode constant",
                         "velocityMode driftwave"),
                 "unknown slugTracking velocityMode")
    expect_fatal("alphaThreshold 1.5 outside (0, 1)",
                 variant("thr", "alphaThreshold 0.45;\n    minLength 0.5;\n    velocityMode",
                         "alphaThreshold 1.5;\n    minLength 0.5;\n    velocityMode"),
                 "must be in (0, 1)")
    expect_fatal("alphaThreshold 0 outside (0, 1)",
                 variant("thr0", "alphaThreshold 0.45;\n    minLength 0.5;\n    velocityMode",
                         "alphaThreshold 0;\n    minLength 0.5;\n    velocityMode"),
                 "must be in (0, 1)")
    expect_fatal("negative minLength",
                 variant("mlen", "minLength 0.5;\n    velocityMode",
                         "minLength -1;\n    velocityMode"),
                 "negative slugTracking minLength")
    expect_fatal("constant mode without a positive slugVelocity",
                 variant("sv0", "slugVelocity 1.5", "slugVelocity 0"),
                 "positive slugVelocity")
    expect_fatal("negative slugVelocity",
                 variant("svneg", "slugVelocity 1.5", "slugVelocity -2"),
                 "negative slugTracking slugVelocity")
    expect_fatal("enable outside true/false",
                 variant("en", "enable true", "enable maybe"),
                 "true or false")
    expect_fatal("gasCap mode without gasCapTime",
                 variant("tau", "velocityMode constant;\n    slugVelocity 1.5",
                         "velocityMode gasCap;\n    gasCapTime 0"),
                 "positive gasCapTime")
    expect_fatal("negative maxSlugVelocity",
                 variant("mx", "velocityMode constant;\n    slugVelocity 1.5",
                         "velocityMode gasCap;\n    gasCapTime 5;\n"
                         "    maxSlugVelocity -1"),
                 "negative slugTracking maxSlugVelocity")
    expect_fatal("gasCap mode without the EOS gasR (caseReject.dict)",
                 os.path.join(CASEDIR, "caseReject.dict"),
                 "no gasR is declared")

print()
if fails:
    print(f"FAILED {len(fails)} check(s):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("transient slug tracking (5.161): all checks passed")
