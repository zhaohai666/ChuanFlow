#!/usr/bin/env python3
"""5.220 verification: erosion / corrosion HYDRAULIC feedback on the bore.

usage: check_two_fluid_erosion_corrosion_feedback.py <exe> <case.dict> <tmpdir>

The 5.187 riskReport only reads the risk out; declaring hydraulicFeedback{}
lets it move the bore. The checker re-derives both arms independently from
the reported state of every frame:

  * EROSION arm (case.dict): the section's own API RP 14E limit
    Ve = 0.3048 * erosionC / sqrt(rho_m / 16.0184634), with rho_m priced on
    the frame's own (alpha, p) at the declared temperature; where |vm| > Ve
    the bore GROWS by exactly erosionRate * dt that step:
        dEff_n = dEff_{n-1} + dt * erosionRate * [ |vm_n| / Ve_n > 1 ]
  * CORROSION arm (caseCorrosion.dict): while the declared H2S fraction at
    the highest section pressure is >= 3447 Pa (0.5 psia) the bore grows by
    corrosionRate * dt on every section (the slow line keeps the erosion
    channel quiet - the checker asserts that too).
  * PASSIVE twin (casePassive.dict, hydraulicFeedback{} undeclared) prints no
    dEff column; the declared arm with both rates zeroed is byte-identical to
    it once the gated column is stripped; the 5.187 read-only deck is re-run
    as the witness that its d_eff still stays D sqrt(1 - wax).
  * loud refusals: erosionRate < 0, corrosionRate < 0, erosionC <= 0, h2s
    outside [0,1], unknown key, corrosion without the gas EOS, feedback on a
    network row.
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
                   "cells": {}}
            continue
        if cur is None:
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


def strip_deff(text):
    return "\n".join(re.sub(r" dEff [-+0-9.eE]+", "", l).rstrip()
                     for l in text.splitlines() if l.strip())


def main():
    base = CASE.read_text()
    src = "\n".join(l.split("//")[0] for l in base.splitlines())
    dt = scal(re.search(r"time\s*\{(.*?)\}", src, re.S).group(1), "dt")
    D = scal(re.search(r"pipe\s*\{(.*?)\}", src, re.S).group(1), "diameter")
    fl = re.search(r"fluid\s*\{(.*?)\}", src, re.S).group(1)
    rhoL = scal(fl, "rhoL")
    gasR = scal(fl, "gasR")
    T0 = scal(re.search(r"initial\s*\{(.*?)\}", src, re.S).group(1),
              "temperature")

    def arm_decls(path):
        s = "\n".join(l.split("//")[0] for l in Path(path).read_text()
                      .splitlines())
        hf = re.search(r"hydraulicFeedback\s*\{(.*?)\n\}", s, re.S).group(1)
        return (scal(hf, "erosionC"), scal(hf, "erosionRate"),
                scal(hf, "h2s"), scal(hf, "corrosionRate"))

    C, eroRate, h2s, corRate = arm_decls(CASE)

    def ve(rhoM):
        return 0.3048 * C / math.sqrt(rhoM / 16.0184634)

    def rho_m(c):
        return c["alpha"] * rhoL + (1.0 - c["alpha"]) * c["p"] / (gasR * T0)

    # ---- 1. the erosion arm ----------------------------------------------
    rc, out = run(CASE)
    ok("erosion arm exits 0", rc == 0, f"rc={rc}")
    fr = parse(out)
    ok("erosion arm reports dEff every frame",
       len(fr) > 2 and all("dEff" in c for f in fr
                           for c in f["cells"].values()))

    worst = 0.0
    grew = 0
    over_any = 0
    for n in range(1, len(fr)):
        for i, c in fr[n]["cells"].items():
            prev = fr[n - 1]["cells"][i]["dEff"]
            ratio = abs(fr[n]["vm"]) / ve(rho_m(c))
            if ratio > 1.0:
                over_any += 1
                grew += 1
            want = prev + dt * eroRate * (1.0 if ratio > 1.0 else 0.0)
            worst = max(worst, relerr(c["dEff"], want))
    ok("erosion: dEff_n = dEff_{n-1} + dt*rate above Ve (worst %.3g)" % worst,
       worst < 1e-9)
    ok("erosion: the run really crosses the limit", over_any > 0)
    ok("erosion: dEff grows every frame while over the limit",
       grew == (len(fr) - 1) * len(fr[0]["cells"]))
    ok("erosion: the bore ends wider than declared",
       all(c["dEff"] > fr[0]["cells"][i]["dEff"]
           for i, c in fr[-1]["cells"].items()))
    ok("erosion: frame 0 starts at the declared bore",
       all(abs(c["dEff"] - D) < 1e-12 for c in fr[0]["cells"].values()))

    # ---- 2. the corrosion arm --------------------------------------------
    Cc, eroRateC, h2sC, corRateC = arm_decls(CASE_DIR / "caseCorrosion.dict")

    def ve_c(rhoM):
        return 0.3048 * Cc / math.sqrt(rhoM / 16.0184634)

    rc, out = run(CASE_DIR / "caseCorrosion.dict")
    ok("corrosion arm exits 0", rc == 0, f"rc={rc}")
    fc = parse(out)
    worst = 0.0
    sour_any = 0
    over_any = 0
    for n in range(1, len(fc)):
        pmax = max(c["p"] for c in fc[n]["cells"].values())
        sour = h2sC * pmax >= 3447.0
        sour_any += 1 if sour else 0
        for i, c in fc[n]["cells"].items():
            ratio = abs(fc[n]["vm"]) / ve_c(rho_m(c))
            if ratio > 1.0:
                over_any += 1
            prev = fc[n - 1]["cells"][i]["dEff"]
            want = prev + dt * corRateC * (1.0 if sour else 0.0)
            worst = max(worst, relerr(c["dEff"], want))
    ok("corrosion: dEff_n = dEff_{n-1} + dt*rate while sour (worst %.3g)"
       % worst, worst < 1e-9)
    ok("corrosion: the sour line really holds (0.5 psia crossed)",
       sour_any == len(fc) - 1)
    ok("corrosion: the slow line never crosses Ve (erosion channel quiet)",
       over_any == 0)
    ok("corrosion: the bore ends wider than declared",
       all(c["dEff"] > fc[0]["cells"][i]["dEff"]
           for i, c in fc[-1]["cells"].items()))

    # ---- 3. the passive twins --------------------------------------------
    rc, outp = run(CASE_DIR / "casePassive.dict")
    ok("passive twin exits 0", rc == 0, f"rc={rc}")
    ok("undeclared feedback prints no dEff column", "dEff" not in outp)

    inert = re.sub(r"hydraulicFeedback\s*\{.*?\n\}",
                   "hydraulicFeedback { erosionC 100; erosionRate 0;"
                   " h2s 0; corrosionRate 0; }", src, flags=re.S)
    TMP.mkdir(parents=True, exist_ok=True)
    p = TMP / "case_inert.dict"
    p.write_text(inert)
    rc, outi = run(p)
    ok("inert-declared arm exits 0", rc == 0, f"rc={rc}")
    ok("inert feedback is byte-identical to the undeclared twin",
       strip_deff(outi) == strip_deff(outp))

    # the 5.187 read-only deck is untouched by the new code
    ref = CASE_DIR.parent / "twoFluidRiskReport" / "case.dict"
    rc, outr = run(ref)
    ok("5.187 deck still exits 0", rc == 0, f"rc={rc}")
    refTxt = "\n".join(l.split("//")[0] for l in ref.read_text().splitlines())
    refD = scal(re.search(r"pipe\s*\{(.*?)\}", refTxt, re.S).group(1),
                "diameter")
    refW = scal(re.search(r"wax\s*\{(.*?)\n\}", refTxt, re.S).group(1),
                "initial")
    fr5 = parse(outr)
    w5 = max(relerr(c["dEff"], refD * math.sqrt(max(1.0 - refW, 1e-12)))
             for f in fr5 for c in f["cells"].values() if "dEff" in c)
    ok("5.187 dEff stays D sqrt(1 - wax) (worst %.3g)" % w5, w5 < 1e-9)

    # ---- 4. loud refusals --------------------------------------------------
    def variant(tag, new, old):
        txt = base.replace(old, new)
        q = TMP / ("case_reject_" + tag + ".dict")
        q.write_text(txt)
        rc, out = run(q)
        ok("refuses " + tag, rc != 0 and "error" in out.lower(),
           out.splitlines()[-1] if out else "no output")

    variant("negativeErosionRate", "erosionRate -1;", "erosionRate 1.0e-5;")
    variant("negativeCorrosionRate", "corrosionRate -1;", "corrosionRate 0;")
    variant("zeroErosionC", "erosionC 0;", "erosionC 100;")
    variant("h2sOutOfRange", "h2s 1.5;", "h2s 0;")
    variant("unknownKey", "    erosionX 1;\n", "    erosionC 100;\n")
    noEos = src.replace("gasR 470; ", "")
    noEos = re.sub(r"hydraulicFeedback\s*\{.*?\n\}",
                   "hydraulicFeedback { erosionC 100; erosionRate 0;"
                   " h2s 0.01; corrosionRate 1e-6; }", noEos, flags=re.S)
    q = TMP / "case_reject_corroNoEos.dict"
    q.write_text(noEos)
    rc, out = run(q)
    ok("refuses corrosionWithoutEos",
       rc != 0 and "declare the gas EOS" in out,
       out.splitlines()[-1] if out else "no output")

    net = TMP / "case_reject_networkRow.dict"
    net.write_text("time { start 0; dt 0.1; endTime 0.1; reportEvery 1; }\n"
                   "network { pipes {} }\n"
                   "hydraulicFeedback { erosionRate 1e-5; }\n")
    rc, out = run(net)
    ok("refuses hydraulicFeedbackOnNetworkRow",
       rc != 0 and "hydraulicFeedback is a serial-row feature" in out,
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
