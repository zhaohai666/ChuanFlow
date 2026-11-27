#!/usr/bin/env python3
"""5.229 (H1b) verification: six-equation distributed phase momentum.

usage: check_two_fluid_six_equation.py <exe> <case.dict> <tmpdir>

Declaring sixEquation{} replaces the lumped mixture-momentum row by one
implicit-Euler momentum row PER PHASE PER CELL, driven on the cell's own EOS
pressure gradient (ghost faces carry the adjacent cell's pressure), with the
wall friction priced by the same Churchill law on each phase's own Reynolds
number and the dynamic closure's interfacial drag K = C|uG-uL| linearised
semi-implicit through the origin secant. The checker re-derives the exact
2x2 recursion frame by frame from the reported (alpha, p, uG, uL) columns:

    rhoG = p / (gasR T0)   (the declared ghost rhoG where the cell is liquid-full)
    dpdx = (pR - pL) / dx  (face-centred, ghost = adjacent cell at the ends)
    r_k  = churchill(Re_k, rr) 0.5 rho_k |u_k| u_k / d   (lagged at frame n)
    mG = ag rhoG / dt,  mL = a rhoL / dt,  K = C |uG - uL|
    bG = mG uG - ag (dpdx + rG),  bL = mL uL - a (dpdx + rL)
    det = mG mL + K (mG + mL)
    uG' = ((mL + K) bG + K bL) / det,  uL' = (K bG + (mG + K) bL) / det

and the mixture velocity the transport reads, as the arc-weighted
mass-weighted average of the phase velocities (priced on the state the
momentum solve marched from - the previous frame - since the transport
moves alpha and p after it; with writeBack the acoustic step re-derives it
on the final field). Arms:

  * MAIN (case.dict): the 2x2 recursion, the reported slip s == uG - uL, the
    vm identity, a sheared line that relaxes under a developing gradient.
  * UNIFORM (caseUniform.dict): dp/dx = 0 exactly - p stays at the
    EOS-consistent 2.115e6, alpha at the declared 0.5 - and the slip decays
    monotonically onto its lagged terminal value
        s* = sign(R) sqrt(|R| / (C (1/(ag rhoG) + 1/(a rhoL)))),
        R = rL/rhoL - rG/rhoG.
  * WRITE-BACK (caseAcoustic.dict): acoustics{} + writeBack true - the
    reported uG/uL equal the momentum result from the previous frame PLUS
    the average of the frame's own reported acoustic face velocities,
    weighted (1 - alpha) for the gas and alpha for the liquid.
  * PASSIVE twin (casePassive.dict) prints no uG/uL keys and no header; the
    enable-false variant of the main arm is byte-identical to it; the two
    standing decks (5.183 acoustics, 5.163 inhibitor equilibrium) witness
    the untouched undeclared path.
  * loud refusals: unknown key, half uG0/uL0 pair, writeBack with enable
    false, writeBack without acoustics, no gas EOS, non-dynamic closure,
    phaseRates / sealed boundary, sixEquation on a network row.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

EXE, CASE, TMP = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
CASE_DIR = CASE.parent

checks = 0
failures = []


def ok(name, cond, detail=""):
    global checks
    checks += 1
    if cond:
        print("ok:  ", name)
    else:
        msg = name + (" | " + detail if detail else "")
        print("FAIL:", msg)
        failures.append(msg)


def relerr(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def close(a, b, rel=1e-8, abs_=1e-9):
    return relerr(a, b) < rel or abs(a - b) < abs_


def run(path):
    out = subprocess.run([EXE, str(path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


def parse(text):
    frames = []
    cur = None
    for line in text.splitlines():
        m = re.match(r"\[chflow\] t (\S+) vm (\S+) dp (\S+)", line)
        if m:
            if cur is not None:
                frames.append(cur)
            cur = {"t": float(m.group(1)), "vm": float(m.group(2)),
                   "dp": float(m.group(3)), "cells": {}, "au": {}}
            continue
        if cur is None:
            continue
        m = re.match(r"\[chflow\] aunode (\d+) u (\S+)", line)
        if m:
            cur["au"][int(m.group(1))] = float(m.group(2))
            continue
        m = re.match(r"\[chflow\] cell (\d+) (.*)", line)
        if m:
            i = int(m.group(1))
            tok = m.group(2).split()
            d = {}
            j = 0
            while j + 1 < len(tok):
                try:
                    d[tok[j]] = float(tok[j + 1])
                except ValueError:
                    pass
                j += 2
            cur["cells"][i] = d
    if cur is not None:
        frames.append(cur)
    return frames


def scal(block, key):
    return float(re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", block).group(1))


def params(path):
    """the declared constants every arm's recursion is priced with"""
    s = "\n".join(l.split("//")[0] for l in Path(path).read_text().splitlines())
    pipe = re.search(r"pipe\s*\{(.*?)\}", s, re.S).group(1)
    fl = re.search(r"fluid\s*\{(.*?)\}", s, re.S).group(1)
    se = re.search(r"sixEquation\s*\{(.*?)\n\}", s, re.S).group(1)
    p = {
        "dt": scal(re.search(r"time\s*\{(.*?)\}", s, re.S).group(1), "dt"),
        "D": scal(pipe, "diameter"),
        "rr": scal(pipe, "roughness") / scal(pipe, "diameter"),
        "dx": scal(pipe, "length") / int(scal(pipe, "sections")),
        "rhoL": scal(fl, "rhoL"),
        "rhoGhost": scal(fl, "rhoG"),
        "muL": scal(fl, "muL"),
        "muG": scal(fl, "muG"),
        "gasR": scal(fl, "gasR"),
        "T0": scal(re.search(r"initial\s*\{(.*?)\}", s, re.S).group(1),
                   "temperature"),
        "C": scal(re.search(r"closure\s*\{(.*?)\}", s, re.S).group(1),
                  "interfacialDrag"),
        "uG0": scal(se, "uG0"),
        "uL0": scal(se, "uL0"),
    }
    return p


def churchill(re_, rr):
    re_ = max(re_, 1e-12)
    a1 = (7.0 / re_) ** 0.9 + 0.27 * rr
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re_) ** 16.0
    return 8.0 * ((8.0 / re_) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def rhoG_of(c, p):
    return c["p"] / (p["gasR"] * p["T0"]) \
        if 1.0 - c["alpha"] > 1e-12 else p["rhoGhost"]


def rphase(rho, mu, u, p):
    v = abs(u)
    if v == 0:
        return 0.0
    return churchill(rho * v * p["D"] / mu, p["rr"]) \
        * 0.5 * rho * v * u / p["D"]


def momentum_step(fr, p):
    """the 2x2 phase rows marched from a frame's reported state"""
    cells = fr["cells"]
    n = len(cells)
    res = {}
    for i in range(n):
        c = cells[i]
        a = c["alpha"]
        ag = 1.0 - a
        rg = rhoG_of(c, p)
        pFaceL = c["p"] if i == 0 else 0.5 * (cells[i - 1]["p"] + c["p"])
        pFaceR = (c["p"] if i == n - 1
                  else 0.5 * (c["p"] + cells[i + 1]["p"]))
        dpdx = (pFaceR - pFaceL) / p["dx"]
        rG = rphase(rg, p["muG"], c["uG"], p)
        rL = rphase(p["rhoL"], p["muL"], c["uL"], p)
        if a <= 1e-12:
            ug = c["uG"] + p["dt"] * (-dpdx - rG) / rg
            res[i] = (ug, ug)
        elif ag <= 1e-12:
            ul = c["uL"] + p["dt"] * (-dpdx - rL) / p["rhoL"]
            res[i] = (ul, ul)
        else:
            mG = ag * rg / p["dt"]
            mL = a * p["rhoL"] / p["dt"]
            K = p["C"] * abs(c["uG"] - c["uL"])
            bG = mG * c["uG"] - ag * (dpdx + rG)
            bL = mL * c["uL"] - a * (dpdx + rL)
            det = mG * mL + K * (mG + mL)
            res[i] = (((mL + K) * bG + K * bL) / det,
                      (K * bG + (mG + K) * bL) / det)
    return res


def vm_of(vel, state, p):
    """the arc-weighted mass-weighted mixture velocity the transport reads"""
    mom = mass = 0.0
    for i, c in vel["cells"].items():
        s = state["cells"][i]
        a = s["alpha"]
        ag = 1.0 - a
        rg = rhoG_of(s, p)
        mom += (ag * rg * c["uG"] + a * p["rhoL"] * c["uL"]) * p["dx"]
        mass += (ag * rg + a * p["rhoL"]) * p["dx"]
    return mom / mass if mass > 0 else 0.0


def slip_star(fr, p):
    """the lagged terminal slip of the pure slip-relaxation (dp/dx = 0) limit"""
    c = fr["cells"][0]
    a = c["alpha"]
    ag = 1.0 - a
    rg = rhoG_of(c, p)
    R = rphase(p["rhoL"], p["muL"], c["uL"], p) / p["rhoL"] \
        - rphase(rg, p["muG"], c["uG"], p) / rg
    dcoef = 1.0 / (ag * rg) + 1.0 / (a * p["rhoL"])
    return math.copysign(math.sqrt(abs(R) / (p["C"] * dcoef)), R)


def check_recursion(tag, fr, p, writeback=False):
    worst = worst_vm = worst_wb = 0.0
    for n in range(1, len(fr)):
        pred = momentum_step(fr[n - 1], p)
        cur = fr[n]
        for i, c in cur["cells"].items():
            pug, pul = pred[i]
            if writeback:
                up = 0.5 * (cur["au"][i] + cur["au"][i + 1])
                worst_wb = max(worst_wb, abs(c["uG"] - pug))
                pug += up * (1.0 - c["alpha"])
                pul += up * c["alpha"]
            worst = max(worst, relerr(c["uG"], pug), relerr(c["uL"], pul))
        # the mixture velocity the transport read was derived on the state
        # the momentum solve marched from; with writeBack the acoustic step
        # re-derives it on the frame's own final field
        state = cur if writeback else fr[n - 1]
        worst_vm = max(worst_vm, relerr(cur["vm"], vm_of(cur, state, p)))
    ok(tag + ": 2x2 recursion reproduced frame by frame (worst %.3g)" % worst,
       worst < 1e-8)
    ok(tag + ": reported vm = mass-weighted phase mixture (worst %.3g)"
       % worst_vm, worst_vm < 1e-8)
    return worst_wb


def main():
    base = CASE.read_text()
    pm = params(CASE)

    # ---- 1. the main arm ---------------------------------------------------
    rc, out = run(CASE)
    ok("main arm exits 0", rc == 0, f"rc={rc}")
    ok("main arm names its mode in the header",
       re.search(r"sixEquation: writeBack false uG0 3 uL0 1", out) is not None)
    fr = parse(out)
    ok("main arm reports uG/uL every frame",
       len(fr) > 2 and all("uG" in c and "uL" in c for f in fr
                           for c in f["cells"].values()))
    ok("main arm starts at the declared phase velocities",
       all(close(c["uG"], pm["uG0"], 1e-12, 1e-12)
           and close(c["uL"], pm["uL0"], 1e-12, 1e-12)
           for c in fr[0]["cells"].values()))
    check_recursion("main", fr, pm)
    ok("main arm: reported slip s = uG - uL",
       all(abs(c["s"] - (c["uG"] - c["uL"])) < 1e-9 for f in fr
           for c in f["cells"].values()))
    ok("main arm: the gradient really drives the line",
       abs(fr[-1]["vm"]) > abs(fr[0]["vm"]) and fr[-1]["dp"] > 0.0)

    # ---- 2. the uniform arm --------------------------------------------------
    uni = CASE_DIR / "caseUniform.dict"
    pu = params(uni)
    rc, out = run(uni)
    ok("uniform arm exits 0", rc == 0, f"rc={rc}")
    fu = parse(out)
    ok("uniform arm: p stays at the EOS-consistent 2.115e6 exactly",
       all(abs(c["p"] - 2115000.0) < 1e-6 for f in fu
           for c in f["cells"].values()))
    ok("uniform arm: alpha stays at the declared 0.5 exactly",
       all(abs(c["alpha"] - 0.5) < 1e-12 for f in fu
           for c in f["cells"].values()))
    check_recursion("uniform", fu, pu)

    worst_mono = 0.0
    decay = True
    for n in range(1, len(fu)):
        sstar = slip_star(fu[n - 1], pu)
        s_prev = fu[n - 1]["cells"][0]["uG"] - fu[n - 1]["cells"][0]["uL"]
        s_now = fu[n]["cells"][0]["uG"] - fu[n]["cells"][0]["uL"]
        worst_mono = max(worst_mono, abs(s_now - sstar) - abs(s_prev - sstar))
        decay = decay and abs(s_now) <= abs(s_prev) + 1e-12
    ok("uniform arm: slip approaches its lagged terminal s* monotonically "
       "(worst excess %.3g)" % worst_mono, worst_mono <= 1e-12)
    ok("uniform arm: the declared slip 2.0 decays every frame", decay)
    s_last = fu[-1]["cells"][0]["uG"] - fu[-1]["cells"][0]["uL"]
    sstar_last = slip_star(fu[-1], pu)
    ok("uniform arm: the final slip has decayed 95%% off the declared 2.0 "
       "and still rides above its terminal value (%.6g vs s* %.6g)"
       % (s_last, sstar_last),
       abs(s_last) < 0.05 * abs(pu["uG0"] - pu["uL0"])
       and sstar_last < s_last < pu["uG0"] - pu["uL0"])

    # ---- 3. the write-back arm -----------------------------------------------
    ac = CASE_DIR / "caseAcoustic.dict"
    pw = params(ac)
    rc, out = run(ac)
    ok("write-back arm exits 0", rc == 0, f"rc={rc}")
    ok("write-back arm names writeBack true",
       re.search(r"sixEquation: writeBack true uG0 2 uL0 1", out) is not None)
    fw = parse(out)
    ok("write-back arm reports one acoustic velocity per face",
       all(len(f["au"]) == len(f["cells"]) + 1 for f in fw))
    wb = check_recursion("write-back", fw, pw, writeback=True)
    ok("write-back arm: the acoustic step really moves the phase momenta "
       "(max increment %.3g)" % wb, wb > 1e-6)
    pa = [max(abs(c["pA"]) for c in f["cells"].values()) for f in fw]
    ok("write-back arm: the acoustic overpressure decays and vents "
       "(peak %.6g -> %.6g, min %.6g)" % (pa[0], pa[-1], min(pa)),
       pa[-1] < pa[0] and min(pa) < 0.5 * pa[0])

    # ---- 4. the passive twins --------------------------------------------------
    rc, outp = run(CASE_DIR / "casePassive.dict")
    ok("passive twin exits 0", rc == 0, f"rc={rc}")
    ok("undeclared sixEquation prints no uG/uL keys and no header",
       " uG " not in outp and "sixEquation" not in outp)

    TMP.mkdir(parents=True, exist_ok=True)
    inert = base.replace("sixEquation\n{\n",
                         "sixEquation\n{\n    enable false;\n")
    p = TMP / "case_inert.dict"
    p.write_text(inert)
    rc, outi = run(p)
    ok("enable-false arm exits 0", rc == 0, f"rc={rc}")
    ok("enable-false arm is byte-identical to the undeclared twin",
       outi == outp)

    # the standing decks are untouched by the new code
    for tag, deck in (("5.183 acoustics", "transientEosAcoustics"),
                      ("5.163 inhibitor equilibrium", "transientInhibitorEq")):
        rc, o = run(CASE_DIR.parent / deck / "case.dict")
        ok(tag + " deck still exits 0", rc == 0, f"rc={rc}")
        ok(tag + " deck prints no sixEquation report",
           " uG " not in o and "sixEquation" not in o)

    # ---- 5. loud refusals ------------------------------------------------------
    def variant(tag, new, old):
        q = TMP / ("case_reject_" + tag + ".dict")
        q.write_text(base.replace(old, new))
        rc, out = run(q)
        ok("refuses " + tag,
           rc != 0 and "sixEquation" in out,
           out.splitlines()[-1] if out else "no output")

    variant("unknownKey", "    uX 1;", "    uG0 3.0;")
    variant("halfVelocityPair", "    uG0 3.0;\n",
            "    uG0 3.0;\n    uL0 1.0;\n")
    variant("writeBackDisabled",
            "sixEquation\n{\n    enable false;\n    writeBack true;\n",
            "sixEquation\n{\n")
    variant("writeBackNoAcoustics", "sixEquation\n{\n    writeBack true;\n",
            "sixEquation\n{\n")
    variant("noGasEOS", "gasR 0;", "gasR 470;")
    variant("closureNotDynamic", "closure { model wallis; slipVelocity 0.5; }",
            "closure { model dynamic; interfacialDrag 2.0; }")
    variant("phaseRatesBoundary",
            "inlet  { type phaseRates; vsl 1.0; vsg 2.0; }",
            "inlet  { type fixedPressure; p 6.0e6; }")
    variant("sealedBoundary", "outlet { type sealed; }",
            "outlet { type fixedPressure; p 5.0e6; }")

    net = TMP / "case_reject_networkRow.dict"
    net.write_text("time { start 0; dt 0.1; endTime 0.1; reportEvery 1; }\n"
                   "network { pipes {} }\n"
                   "sixEquation { uG0 1; uL0 0; }\n")
    rc, out = run(net)
    ok("refuses sixEquationOnNetworkRow",
       rc != 0 and "sixEquation is a serial-row feature" in out,
       out.splitlines()[-1] if out else "no output")

    print()
    print("checks:", checks, "failures:", len(failures))
    if failures:
        for f in failures:
            print("FAILED:", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
