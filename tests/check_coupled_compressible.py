#!/usr/bin/env python3
"""5.79 / A8 increment 2: the 1D<->3D partner region is COMPRESSIBLE.

A single-path steady network (Src -> A -> X(3D) -> B -> Snk) injects a fixed
rate, so q_in is known exactly and dp(q) is priced against an isothermal
ideal-gas closed form with no branch re-splitting to fit around (5.78's lesson).

`detail { compressible { pRef; RsT; mu; } }` makes the partner a plane duct
marched along the mesh's own cross-sections:

    rho = p/(RsT),  mDot = rho_in q_in,  dp/dx = -12 mu u / H^2,
    u = mDot/(rho A)   =>   p dp/dx = -12 mu mDot RsT/(W H^3)
    =>  p_in^2 - p_out^2 = 24 mu L p_in q_in/(W H^3)

L is the first-to-last CELL CENTRE span (0.49 m on 50 x-cells), because that
is the span whose pressure difference the block reports -- the same anchor as
5.78. `qMag` is the INLET volume rate: the outlet rate is a RESULT, and
`fluxOut > fluxIn` is the visible signature of the EOS. RsT cancels out of
dp (it only sets the absolute density), which is itself a sharp criterion.

Every arm and every refusal probe is written into the SYSTEM temp directory;
no scratch dict is ever left in a case directory.

usage: check_coupled_compressible.py <solver-binary> <case.dict>
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


def bracket(txt, key):
    """Body of a named { } block, with NESTED braces."""
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


def scalar(txt, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)\s*;", txt)
    return float(m.group(1)) if m else default


def numbers(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", txt)
    return [float(x) for x in m.group(1).split()] if m else []


# ---- independent closed form, rebuilt from the case text ------------------
FB = body(BASE)
RHO = scalar(bracket(FB, "fluid"), "rho")
MU = scalar(bracket(FB, "fluid"), "mu")
BLOCKTXT = bracket(FB, "block")
LO, HI, CELLS = (numbers(BLOCKTXT, k) for k in ("min", "max", "cells"))
AXIS, SIGN = 0, 1
for k, name in enumerate(("x", "y", "z")):
    if re.search(name + r"min\s*\{[^}]*type\s+inlet", FB):
        AXIS, SIGN = k, 1
    if re.search(name + r"max\s*\{[^}]*type\s+inlet", FB):
        AXIS, SIGN = k, -1
L = abs(HI[AXIS] - LO[AXIS])
NX = CELLS[AXIS]
L_MEAS = L * (NX - 1) / NX          # first-to-last CELL CENTRE span
T = [abs(HI[k] - LO[k]) for k in range(3) if k != AXIS]
H, W = min(T), max(T)
AREA = W * H

CP = re.compile(r"    compressible \{.*?\n    \}\n", re.S)


def arm(text):
    return CP.sub("    compressible {\n" + text + "\n    }\n", BASE)


def closed(txt):
    """dict-only isothermal gas solution, or None for an undeclared case."""
    b = body(txt)
    # ONLY detail{} owns the declaration. An identically named block at the
    # top level is not a partner declaration, and the loader ignores it.
    cp = bracket(bracket(b, "detail"), "compressible")
    if not cp.strip():
        return None
    pIn = scalar(cp, "pRef")
    RsT = scalar(cp, "RsT")
    mu = scalar(cp, "mu") or MU
    qIn = scalar(bracket(b, "feed"), "rate")
    rhoIn = pIn / RsT
    mDot = rhoIn * qIn
    K = 24.0 * mu * L_MEAS * pIn * qIn / (W * H ** 3)
    p2Out = pIn ** 2 - K
    if p2Out <= 0:
        return {"choked": True}
    pOut = math.sqrt(p2Out)
    rhoOut = pOut / RsT
    qOut = mDot / rhoOut
    return {
        "choked": False,
        "pIn": pIn, "pOut": pOut, "dp": pIn - pOut,
        "rhoIn": rhoIn, "rhoOut": rhoOut,
        "qIn": qIn, "qOut": qOut, "massRate": mDot,
        "uIn": qIn / AREA, "uOut": qOut / AREA,
        "p2Slope": -K / L_MEAS,
        "length": L_MEAS, "height": H, "width": W, "area": AREA,
        "sections": NX, "RsT": RsT, "mu": mu,
    }


def run(txt, tag):
    """One arm in the system temp; returns the numbers it printed."""
    p = os.path.join(TD, tag + ".dict")
    open(p, "w", encoding="utf-8").write(txt)
    r = subprocess.run([EXE, p], capture_output=True, text=True, cwd=TD)
    if r.returncode != 0:
        raise SystemExit(f"arm {tag} failed rc={r.returncode}\n"
                        f"{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    out = {"stations": {}, "seeds": {}, "hasCompressible": False}
    for line in r.stdout.splitlines():
        if "coupled solved:" in line:
            for key in ("q", "dp", "R", "iterations", "residual"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    out[key] = float(m.group(1))
        elif "block seed" in line:
            for key in ("sweeps", "fluxIn", "fluxOut", "planeFluxMin",
                        "planeFluxMax", "probeSlope", "probeSpread",
                        "probeSections", "dpPa", "uxMax", "uxCentre"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    out[key] = float(m.group(1))
        elif "block compressible seed" in line:
            out["hasCompressible"] = True
            for key in ("pIn", "pOut", "dp", "rhoIn", "rhoOut", "qIn",
                        "qOut", "massRate", "uIn", "uOut", "p2Slope",
                        "sections", "length", "height", "width", "area"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    out[key] = float(m.group(1))
        elif "block compressible station" in line:
            k = int(re.search(r"station\s+(\d+)", line).group(1))
            st = {}
            for key in ("x", "p", "rho", "q", "u"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    st[key] = float(m.group(1))
            out["stations"][k] = st
        elif "coupled seed " in line:
            s = int(re.search(r"coupled seed\s+(\d+)", line).group(1))
            sd = {}
            for key in ("q", "dp", "R"):
                m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", line)
                if m:
                    sd[key] = float(m.group(1))
            out["seeds"][s] = sd
    return out


# A `compressible{}` OUTSIDE detail{} is not a partner declaration: it is the
# zero-declaration twin, differing from `plain` only by text the loader ignores.
TOPLEVEL = (CP.sub("", BASE) + "\n"
            + "// same words, but outside detail{}: not a partner declaration\n"
            + CP.search(BASE).group(0).replace("    compressible {",
                                               "compressible {")
                                      .replace("\n    }\n", "\n}\n"))

with tempfile.TemporaryDirectory(prefix="cf_compressible_") as TD:
    ARMS = {
        "plain": CP.sub("", BASE),
        "toplevel": TOPLEVEL,
        "main": BASE,
        "pRef_low": arm("        pRef 5000; RsT 300000;"),
        "pRef_high": arm("        pRef 20000; RsT 300000;"),
        "rate_half": BASE.replace("rate 2.4e-4;", "rate 1.2e-4;"),
        "rate_double": BASE.replace("rate 2.4e-4;", "rate 4.8e-4;"),
        "mu_decl": arm("        pRef 10000; RsT 300000; mu 0.2;"),
        "RsT_high": arm("        pRef 10000; RsT 600000;"),
    }
    R = {name: run(txt, name) for name, txt in ARMS.items()}
    DECLARED = [n for n in ARMS if closed(ARMS[n]) is not None]

    # ---- (1) a zero declaration is the old SIMPLE solve, bit for bit ------
    for key in ("dp", "q", "R", "sweeps", "fluxIn", "fluxOut",
                "probeSlope", "probeSpread", "uxMax", "uxCentre", "dpPa"):
        close(f"zero declaration: {key} identical when compressible{{}} sits"
              f" outside detail{{}}",
              R["toplevel"][key], R["plain"][key], 0.0, scale=1.0)
    check("an undeclared partner prints no compressible state",
          not R["plain"]["hasCompressible"])
    # SIMPLE's inlet/outlet patch fluxes are discrete sums, so they agree to
    # the linear-solver tolerance rather than bit for bit.
    close("an undeclared partner conserves the VOLUME rate",
          R["plain"]["fluxIn"], R["plain"]["fluxOut"], 1e-9, scale=1.0)
    check("an undeclared partner still runs SIMPLE sweeps, not a march",
          R["plain"]["sweeps"] != R["main"]["sections"],
          f"sweeps {R['plain']['sweeps']} vs sections {R['main']['sections']}")

    # ---- (2) every arm keeps the imposed rate ----------------------------
    for name, txt in ARMS.items():
        close(f"arm {name}: the single path carries exactly the injected rate",
              R[name]["q"], scalar(bracket(body(txt), "feed"), "rate"), 1e-12)

    # ---- (3) the printed state IS the dict-only closed form ---------------
    FULL = ("pIn", "pOut", "dp", "rhoIn", "rhoOut", "qIn", "qOut",
            "massRate", "uIn", "uOut", "p2Slope", "length", "height",
            "width", "area", "sections")
    KEY6 = ("pOut", "dp", "rhoOut", "qOut", "massRate", "p2Slope")
    for name in DECLARED:
        c = closed(ARMS[name])
        keys = FULL if name == "main" else KEY6
        for key in keys:
            close(f"arm {name}: printed {key} = dict-only closed form",
                  R[name][key], c[key], 1e-9, scale=max(abs(c[key]), 1.0))
        close(f"arm {name}: the network is charged exactly the block drop",
              R[name]["dpPa"], c["dp"], 1e-9, scale=max(abs(c["dp"]), 1.0))

    # ---- (4) the march really is a field, not a lookup table -------------
    c = closed(BASE)
    for k, st in sorted(R["main"]["stations"].items()):
        close(f"station {k}: p^2 = pIn^2 + (dp2/dx) x",
              st["p"] ** 2, c["pIn"] ** 2 + c["p2Slope"] * st["x"],
              1e-9, scale=c["pIn"] ** 2)
        close(f"station {k}: rho = p / RsT",
              st["rho"], st["p"] / c["RsT"], 1e-9,
              scale=max(st["rho"], 1.0))
        close(f"station {k}: rho*q = the conserved mass rate",
              st["rho"] * st["q"], c["massRate"], 1e-9,
              scale=max(c["massRate"], 1.0))
        close(f"station {k}: u = q / A",
              st["u"], st["q"] / c["area"], 1e-9, scale=max(st["u"], 1.0))
    check("five independent stations were printed",
          len(R["main"]["stations"]) == 5,
          f"got {len(R['main']['stations'])}")

    # ---- (5) the EOS is visible in BOTH the field and the block ledger ----
    check("the gas expands: rhoOut < rhoIn and qOut > qIn",
          R["main"]["rhoOut"] < R["main"]["rhoIn"]
          and R["main"]["qOut"] > R["main"]["qIn"],
          f"rho {R['main']['rhoIn']:.6g} -> {R['main']['rhoOut']:.6g}, "
          f"q {R['main']['qIn']:.6g} -> {R['main']['qOut']:.6g}")
    check("the block ledger reports a NON-constant volume rate",
          R["main"]["fluxOut"] > R["main"]["fluxIn"],
          f"fluxIn {R['main']['fluxIn']:.8g} < fluxOut "
          f"{R['main']['fluxOut']:.8g}")
    close("planeFluxMin is the inlet rate", R["main"]["planeFluxMin"],
          R["main"]["qIn"], 1e-9, scale=1.0)
    close("planeFluxMax is the expanded outlet rate",
          R["main"]["planeFluxMax"], R["main"]["qOut"], 1e-9, scale=1.0)
    check("p(x) bows away from a straight line (p^2 is the linear one)",
          R["main"]["probeSpread"] > 10.0 * R["plain"]["probeSpread"],
          f"spread {R['main']['probeSpread']:.4f} Pa vs "
          f"{R['plain']['probeSpread']:.4f} Pa")
    check("the bow grows with the rate",
          R["rate_half"]["probeSpread"] < R["main"]["probeSpread"]
          < R["rate_double"]["probeSpread"])

    # RsT sets the absolute density but cancels out of dp: a partner that
    # folded it into the momentum law would fail these three.
    close("doubling RsT leaves dp bit-identical (RsT cancels in dp)",
          R["RsT_high"]["dp"], R["main"]["dp"], 0.0, scale=1.0)
    close("doubling RsT halves rhoIn",
          R["RsT_high"]["rhoIn"] / R["main"]["rhoIn"], 0.5, 1e-9)
    close("doubling RsT halves the mass rate",
          R["RsT_high"]["massRate"] / R["main"]["massRate"], 0.5, 1e-9)

    # A declared mu must override fluid.mu (K is linear in mu).
    close("a declared mu doubles d(p^2)/dx",
          R["mu_decl"]["p2Slope"] / R["main"]["p2Slope"], 2.0, 1e-9)

    # pRef is an ABSOLUTE pressure: the same rate at a different pressure
    # level gives a different drop. A constant-density look-up cannot pass.
    check("dp falls strictly as pRef rises",
          R["pRef_low"]["dp"] > R["main"]["dp"] > R["pRef_high"]["dp"],
          f"{R['pRef_low']['dp']:.2f} > {R['main']['dp']:.2f} > "
          f"{R['pRef_high']['dp']:.2f} Pa")
    check("the pRef steps are far outside the solver tolerance",
          R["pRef_low"]["dp"] - R["main"]["dp"] > 10.0
          and R["main"]["dp"] - R["pRef_high"]["dp"] > 10.0)

    # dp/q is NOT constant: the isothermal duct is a nonlinear partner.
    ratios = [R[n]["dp"] / R[n]["q"]
              for n in ("rate_half", "main", "rate_double")]
    check("dp/q grows with q: the partner is nonlinear, not R*q",
          ratios[0] < ratios[1] < ratios[2],
          " < ".join(f"{v:.6g}" for v in ratios))
    check("two independent seeds reach the same fixed point",
          len(R["main"]["seeds"]) == 2
          and abs(R["main"]["seeds"][0]["dp"]
                  - R["main"]["seeds"][1]["dp"]) < 1e-9
          and abs(R["main"]["seeds"][0]["R"]
                  - R["main"]["seeds"][1]["R"])
          < 1e-9 * R["main"]["seeds"][0]["R"])
    check("both seeds converged within the quasi-Newton budget",
          R["main"]["iterations"] <= 30 and R["main"]["residual"] < 1e-6,
          f"iterations {R['main']['iterations']}, "
          f"residual {R['main']['residual']:.3g}")

    # ---- (6) loud refusals: never a silently ignored declaration ---------
    probes = [
        ("negative pRef",
         arm("        pRef -1; RsT 300000;"),
         "finite positive absolute pressure"),
        ("zero RsT",
         arm("        pRef 10000; RsT 0;"),
         "finite and positive"),
        ("negative mu",
         arm("        pRef 10000; RsT 300000; mu -1;"),
         "cannot be negative"),
        ("a choked duct",
         arm("        pRef 1000; RsT 300000;"),
         "choked"),
        ("compressible together with phases",
         BASE.replace("    compressible {",
                      "    phases { alpha0 0.5; }\n    compressible {"),
         "cannot both own the partner's density"),
        ("compressible together with a non-Newtonian fluid",
         BASE.replace("coupledResistance 3.853e6;",
                      "coupledResistance 3.853e6;\n"
                      "         flowIndex 0.5; consistency 0.15;"),
         "cannot both own the partner's viscosity"),
        ("compressible together with a RANS wall model",
         BASE.replace("    probe { min 0.12; max 0.49; }",
                      "    turbulence { model kOmegaSST; }\n\n"
                      "    probe { min 0.12; max 0.49; }"),
         "turbulence model cannot both own the partner's viscosity"),
        ("compressible together with an unsteady level",
         BASE.replace("ddt       steady;", "ddt       Euler;"),
         "quasi-steady"),
    ]
    for tag, txt, needle in probes:
        p = os.path.join(TD, "probe_reject.dict")
        open(p, "w", encoding="utf-8").write(txt)
        rr = subprocess.run([EXE, p], capture_output=True, text=True, cwd=TD)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {tag}: non-zero exit naming the reason", ok,
              f"needle '{needle}'")

if fails:
    raise SystemExit(f"coupledCompressible FAILED: {len(fails)} checks -> {fails}")
print("coupledCompressible OK: isothermal EOS duct priced against a dict-only"
      " closed form, p^2(x) linear station by station, compressed outlet,"
      " nonlinear dp/q, zero-declaration twin bit-identical, 8 loud refusals")
