#!/usr/bin/env python3
"""5.78 / A8 increment 1: the 3D partner block carries TWO fluids.

A single-path steady network (Src -> A -> X(3D) -> B -> Snk) injects a fixed
rate, so every arm is measured at EXACTLY the same q: the difference between a
declared arm and its undeclared twin is the phase field's own contribution,
with no branch re-splitting to fit around.

The block declares `detail { phases { rhoL; rhoG; muL; muG; alpha (...);
gravity (gx gy gz); viscosityMix (yes|no); } }`. Only the LIQUID share of each
cell is heavy: the body force is g * rho_m/rho_ref, and because the FVM kernel
keeps p kinematic, multiplying it back by rho_ref turns the measured drop into
the physical head  -L_meas <rho_m> g_flow  on top of the Poiseuille friction.
`viscosityMix yes` also lets alpha price the per-face viscosity band by band.

Every arm (including the undeclared twin and the refusal probes) is derived
from the one shipped case into the SYSTEM temp directory -- no scratch dict is
ever written into a case directory.

usage: check_coupled_multiphase.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASE = sys.argv[1], sys.argv[2]
BASE = open(CASE, encoding="utf-8").read()

fails = []


def check(what, ok, detail=""):
    print(f"{'ok:  ' if ok else 'FAIL:'} {what}" + (f" ({detail})" if detail else ""))
    if not ok:
        fails.append(what)


def close(what, got, want, tol, scale=None):
    s = abs(want) if scale is None else scale
    ok = abs(got - want) <= tol * max(s, 1e-30)
    check(what, ok, f"got={got:.12g} want={want:.12g} tol={tol:.3g}")
    return ok


def body(txt):
    """Strip // comments (the dict dialect's only comment marker)."""
    return "\n".join(line.split("//")[0] for line in txt.splitlines())


def scalar(txt, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)\s*;", txt)
    return float(m.group(1)) if m else default


def word(txt, key, default=""):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z0-9_]+)\s*;", txt)
    return m.group(1) if m else default


def numbers(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", txt)
    return [float(x) for x in m.group(1).split()] if m else []


def bracket(txt, key):
    """Body of a named { } block, with NESTED braces (block{ xmin{...} ... })."""
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", txt)
    if not m:
        return ""
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    return ""


# ---- independent closed form, rebuilt from the case text ------------------
FB = body(BASE)
RHO = scalar(bracket(FB, "fluid"), "rho")
MU = scalar(bracket(FB, "fluid"), "mu")
BLOCK = bracket(FB, "block")
LO, HI = numbers(BLOCK, "min"), numbers(BLOCK, "max")
CELLS = numbers(BLOCK, "cells")
Q_IN = scalar(bracket(FB, "feed"), "rate")

# The inlet patch names both the flow axis and its sign: its outward normal
# points out of the region, so xmin/ymin/zmin inlet => +axis, xmax... => -axis.
AXIS, SIGN = 0, 1
for k, name in enumerate(("x", "y", "z")):
    if re.search(name + r"min\s*\{[^}]*type\s+inlet", FB):
        AXIS, SIGN = k, 1
    if re.search(name + r"max\s*\{[^}]*type\s+inlet", FB):
        AXIS, SIGN = k, -1
L = abs(HI[AXIS] - LO[AXIS])
NX = CELLS[AXIS]
L_MEAS = L * (NX - 1) / NX      # first-to-last CELL CENTRE span

PH = re.compile(r"    phases \{.*?\n    \}\n", re.S)
PH_BASE = bracket(body(BASE), "phases")
RHOL = scalar(PH_BASE, "rhoL") or RHO
RHOG = scalar(PH_BASE, "rhoG") or RHO
# Derived from the shipped case, never hard-coded twice: editing rhoL/rhoG in
# case.dict moves the closed form AND every derived arm together.
RHO_DECL = f"rhoL {RHOL}; rhoG {RHOG};"


def arm(text):
    return PH.sub("    phases {\n" + text + "\n    }\n", BASE)


def closed(txt):
    """(alphaMean, rhoMean, muMean, bodyDropPa, fricRatio) from the dict."""
    b = body(txt)
    ph = bracket(b, "phases")
    if not ph.strip():
        return None
    rhoL = scalar(ph, "rhoL") or RHO
    rhoG = scalar(ph, "rhoG") or RHO
    muL = scalar(ph, "muL") or MU
    muG = scalar(ph, "muG") or MU
    bands = numbers(ph, "alpha")
    a = sum(bands) / len(bands) if bands else scalar(ph, "alpha0")
    rhoM = a * rhoL + (1.0 - a) * rhoG
    muM = a * muL + (1.0 - a) * muG
    g = numbers(ph, "gravity") or [0.0, 0.0, 0.0]
    gFlow = g[AXIS] * SIGN
    drop = -rhoM * gFlow * L_MEAS
    mix = word(ph, "viscosityMix", "no") == "yes"
    return a, rhoM, muM, drop, (muM / MU if mix else 1.0), tuple(bands), mix


def run(txt, tag):
    """One arm in the system temp; returns a dict of the numbers it printed."""
    p = os.path.join(TD, tag + ".dict")
    open(p, "w", encoding="utf-8").write(txt)
    r = subprocess.run([EXE, p], capture_output=True, text=True, cwd=TD)
    if r.returncode != 0:
        raise SystemExit(f"arm {tag} failed rc={r.returncode}\n"
                        f"{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    out = {}
    for line in r.stdout.splitlines():
        if "coupled solved" in line:
            out["dp"] = float(re.search(r"\bdp\s+([-+0-9.eE]+)", line).group(1))
            out["q"] = float(re.search(r"\bq\s+([-+0-9.eE]+)", line).group(1))
        elif "block seed" in line:
            for key in ("sweeps", "probeSlope", "probeSpread", "dpPa", "uxCentre"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    out[key] = float(m.group(1))
        elif "block phases seed" in line:
            for key in ("alphaMean", "rhoMean", "muMean", "bodyDrop"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    out[key] = float(m.group(1))
    return out


with tempfile.TemporaryDirectory(prefix="cf_multiphase_") as TD:
    ARMS = {
        "plain": PH.sub("", BASE),
        "zero": arm("        alpha0 0;"),
        "main": BASE,
        "uniform": arm(f"        {RHO_DECL} alpha ( 1 1 1 1 ); "
                       "gravity ( -2 0 0 );"),
        "gas": arm(f"        {RHO_DECL} alpha ( 0 0 0 0 ); "
                   "gravity ( -2 0 0 );"),
        "half": arm(f"        {RHO_DECL} alpha ( 1 1 0 0 ); "
                    "gravity ( -1 0 0 );"),
        "reverse": arm(f"        {RHO_DECL} alpha ( 1 1 0 0 ); "
                       "gravity ( 0.5 0 0 );"),
        "mix": arm("        muL 0.2; muG 0.05; alpha ( 1 1 1 1 ); "
                   "viscosityMix yes;"),
        "mix_bands": arm("        muL 0.2; muG 0.05; alpha ( 1 1 0 0 ); "
                         "viscosityMix yes;"),
    }
    R = {name: run(txt, name) for name, txt in ARMS.items()}

    # ---- (1) a zero declaration is the old solve, bit for bit -------------
    for key in ("dp", "q", "sweeps", "probeSlope", "dpPa"):
        close(f"undeclared twin: {key} bit-identical with a declared-but-empty"
              f" phases block", R["zero"][key], R["plain"][key], 0.0,
              scale=1.0)
    check("a declared-but-empty phases block still reports its (trivial) field",
          all(k in R["zero"] for k in
              ("alphaMean", "rhoMean", "muMean", "bodyDrop")))

    # ---- (2) every arm keeps the imposed rate -----------------------------
    for name in ARMS:
        close(f"arm {name}: the single path carries exactly the injected rate",
              R[name]["q"], Q_IN, 1e-12)

    # ---- (3) the printed field IS the closed form -------------------------
    for name, txt in ARMS.items():
        c = closed(txt)
        if c is None:
            continue
        a, rhoM, muM, drop, ratio, bands, mix = c
        for key, want in (("alphaMean", a), ("rhoMean", rhoM),
                          ("muMean", muM), ("bodyDrop", drop)):
            close(f"arm {name}: printed {key} = dict-only closed form",
                  R[name][key], want, 1e-9, scale=max(abs(want), 1.0))

    dp0 = R["plain"]["dp"]
    bodyMain = closed(BASE)[3]

    # ---- (4) the head is linear in gravity, and reverses with it ----------
    dMain = R["main"]["dp"] - dp0
    close("half gravity halves the measured head",
          (R["half"]["dp"] - dp0) / dMain, 0.5, 1e-3)
    close("reversing gravity reverses the head (-0.25x)",
          (R["reverse"]["dp"] - dp0) / dMain, -0.25, 1e-3)
    check("reverse < undeclared < main (the head has a direction)",
          R["reverse"]["dp"] < dp0 < R["main"]["dp"])

    # The difference is the analytic head to within the documented O(dx/L)
    # footprint of the discrete pressure gradient and the entry length
    # (measured 0.1 - 0.9% on this 50-cell line; the checker quantifies it).
    for name in ("main", "half", "uniform", "gas", "reverse"):
        want = closed(ARMS[name])[3]
        got = R[name]["dp"] - dp0
        close(f"arm {name}: dp - dp_plain = the analytic distributed head",
              got, want, 0.02, scale=max(abs(want), 1.0))

    # ---- (5) DISTRIBUTED witness: only the liquid share is heavy ----------
    check("stratified arm is neither the all-liquid nor the all-gas head",
          abs(bodyMain - closed(ARMS["uniform"])[3]) > 1.0
          and abs(bodyMain - closed(ARMS["gas"])[3]) > 1.0)
    close("the stratified arm's mean density blends the two bands",
          R["main"]["rhoMean"], 0.5 * RHOL + 0.5 * RHOG, 1e-9, scale=1000.0)
    # The shipped column is deliberately NOT at rho_ref: a block that ignored
    # the phase field and just used the network's rho would sit at rhoMean=800,
    # and a checker that only compared against 800 could never notice.
    check("a declared phase field really changes the block's density",
          abs(R["main"]["rhoMean"] - R["zero"]["rhoMean"]) > 50.0
          and R["gas"]["rhoMean"] < R["main"]["rhoMean"]
          < R["uniform"]["rhoMean"],
          f"gas {R['gas']['rhoMean']} < stratified {R['main']['rhoMean']} "
          f"< liquid {R['uniform']['rhoMean']}, "
          f"undeclared {R['zero']['rhoMean']}")
    check("a frozen interface leaves a visible kink in the probe window",
          R["main"]["probeSpread"] > 10.0 * R["plain"]["probeSpread"],
          f"spread {R['main']['probeSpread']:.3f} vs "
          f"{R['plain']['probeSpread']:.3f} Pa")

    # ---- (6) the same alpha also prices the wall law ---------------------
    close("uniform mixture doubles the developed friction gradient",
          R["mix"]["probeSlope"] / R["plain"]["probeSlope"], 2.0, 1e-2)
    close("stratified mixture prices the wall band by band (muMean)",
          R["mix_bands"]["muMean"], 0.125, 1e-9, scale=1.0)
    check("a stratified wall is neither the plain nor the uniform-mix drop",
          R["plain"]["dp"] < R["mix_bands"]["dp"] < R["mix"]["dp"],
          f"{R['plain']['dp']:.2f} < {R['mix_bands']['dp']:.2f} "
          f"< {R['mix']['dp']:.2f} Pa")
    close("stratified friction is the length-weighted series sum",
          R["mix_bands"]["dp"], dp0 * 1.25, 0.03)
    check("the stratified wall leaves a kink the uniform mix does not",
          R["mix_bands"]["probeSpread"] > 10.0 * R["mix"]["probeSpread"])

    # ---- (7) loud refusals: never a silently ignored declaration ---------
    probes = [
        ("negative alpha band",
         arm(f"        {RHO_DECL} alpha ( 1 -0.1 0 0 );"),
         "must be in [0,1]"),
        ("negative rhoL",
         arm("        rhoL -1; alpha ( 1 0 );"),
         "cannot be negative"),
        ("gravity with two components",
         arm("        alpha ( 1 0 ); gravity ( 1 2 );"),
         "three finite components"),
        ("unknown viscosityMix word",
         arm("        alpha ( 1 0 ); viscosityMix maybe;"),
         "understands 'yes' and 'no'"),
        ("alpha with a single band",
         arm("        alpha ( 1 );"),
         "at least two bands"),
        ("alpha bands together with a uniform alpha0",
         arm("        alpha ( 1 0 ); alpha0 0.5;"),
         "mutually exclusive"),
        ("viscosityMix with a non-Newtonian partner fluid",
         BASE.replace("coupledResistance 3.75e6;",
                      "coupledResistance 3.75e6;\n"
                      "         flowIndex 0.5; consistency 0.15;")
             .replace("gravity ( -2 0 0 );",
                      "gravity ( -2 0 0 ); viscosityMix yes;"),
         "cannot both own the per-face viscosity"),
    ]
    for tag, txt, needle in probes:
        p = os.path.join(TD, "probe_reject.dict")
        open(p, "w", encoding="utf-8").write(txt)
        rr = subprocess.run([EXE, p], capture_output=True, text=True, cwd=TD)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {tag}: non-zero exit naming the reason", ok,
              f"needle '{needle}'")

if fails:
    raise SystemExit(f"coupledMultiphase FAILED: {len(fails)} checks -> {fails}")
print("coupledMultiphase OK: distributed two-phase head + band-priced wall,"
      " zero-declaration twin bit-identical, 7 loud refusals")
