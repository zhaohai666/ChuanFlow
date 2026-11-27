#!/usr/bin/env python3
"""5.227 verification: hydrate BLOCKAGE (deposition layer + particle bridging).

usage: check_two_fluid_hydrate_blockage.py <exe> <case.dict> <tmpdir>

The 5.186 increment proved hydrate h NON-blocking; declaring blockage{}
charges it. The checker re-implements both arms independently:

  * DEPOSITION arm (case.dict): the wall layer obeys the implicit-Euler
    balance  delta_n = (delta_{n-1}/dt + depositionRate*h_n)/(1/dt +
    erosionRate)  re-derived frame by frame from the reported h and blk
    columns, and the bore loses the layer on both walls:
    dEff == D - 2*delta on every section of every frame (no wax is declared,
    so the thin-layer law contributes exactly D).
  * BRIDGING arm (caseBridge.dict): every frame whose final hydrate field
    reaches the packing limit (h >= 0.6 anywhere) must report vm == 0.
  * PASSIVE twin (casePassive.dict, blockage{} undeclared) prints no dEff/blk
    token, and the declared arm with the blockage made inert (depositionRate
    0, packingFraction 1) is byte-identical to it once the gated columns are
    stripped. The 5.186 deck itself is re-run as the witness that the
    non-blocking arm is untouched (dEff stays D sqrt(1 - wax), no blk).
  * loud refusals: packingFraction <= 0, depositionRate < 0, erosionRate < 0,
    unknown key, blockage without hydrate{}, blockage on a network row.
"""
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


def main():
    base = CASE.read_text()
    # comments carry the words we grep for - strip them before extraction
    src = "\n".join(l.split("//")[0] for l in base.splitlines())
    D = scal(re.search(r"pipe\s*\{(.*?)\}", src, re.S).group(1), "diameter")
    bk = re.search(r"blockage\s*\{(.*?)\n\}", src, re.S).group(1)
    dep = scal(bk, "depositionRate")
    ero = scal(bk, "erosionRate")
    pack = scal(bk, "packingFraction")
    dt = scal(re.search(r"time\s*\{(.*?)\}", src, re.S).group(1), "dt")

    # ---- 1. the deposition arm -------------------------------------------
    rc, out = run(CASE)
    ok("deposition arm exits 0", rc == 0, f"rc={rc}")
    fr = parse(out)
    ok("deposition arm has frames with dEff and blk",
       len(fr) > 2 and all("dEff" in c and "blk" in c
                           for f in fr for c in f["cells"].values()))

    # dEff == D - 2*delta on every section of every frame
    wde = 0.0
    for f in fr:
        for c in f["cells"].values():
            wde = max(wde, relerr(c["dEff"], D - 2.0 * c["blk"]))
    ok("dEff == D - 2*delta every frame (worst %.3g)" % wde, wde < 1e-9)

    # the layer balance re-derived frame by frame from the reported h
    wbl = 0.0
    for n in range(1, len(fr)):
        for i, c in fr[n]["cells"].items():
            prev = fr[n - 1]["cells"][i]["blk"]
            want = (prev / dt + dep * c["h"]) / (1.0 / dt + ero)
            wbl = max(wbl, relerr(c["blk"], want))
    ok("delta obeys the implicit deposition/erosion balance (worst %.3g)"
       % wbl, wbl < 1e-9)
    ok("frame 0 starts with no layer",
       all(c["blk"] == 0.0 for c in fr[0]["cells"].values()))
    ok("the layer really grows (final blk = %.3g m > 0)"
       % max(c["blk"] for c in fr[-1]["cells"].values()),
       max(c["blk"] for c in fr[-1]["cells"].values()) > 0)
    ok("the bore shrinks below D once the layer grows",
       all(c["dEff"] < D for f in fr[1:] for c in f["cells"].values()))
    # the packing limit is never reached here (water cut caps h), so the
    # deposition arm must not plug
    ok("no bridging in the deposition arm (h stays below the packing limit)",
       max(c["h"] for f in fr for c in f["cells"].values()) < pack)
    ok("the deposition arm keeps flowing (vm != 0 every frame)",
       all(f["vm"] != 0.0 for f in fr))

    # ---- 2. the bridging arm ---------------------------------------------
    rc, out = run(CASE_DIR / "caseBridge.dict")
    ok("bridging arm exits 0", rc == 0, f"rc={rc}")
    fb = parse(out)
    crossed = [f for f in fb
               if max(c["h"] for c in f["cells"].values()) >= pack]
    ok("the bridge really forms (h >= %.2g somewhere)" % pack, len(crossed) > 0)
    viol = [f["t"] for f in fb
            if max(c["h"] for c in f["cells"].values()) >= pack
            and f["vm"] != 0.0]
    ok("every bridged frame reports vm == 0 (violations: %s)" % viol[:3],
       not viol)
    ok("the line flows before the bridge (early vm != 0)",
       fb[0]["vm"] != 0.0 and fb[1]["vm"] != 0.0)
    firstBridge = next(n for n, f in enumerate(fb)
                       if max(c["h"] for c in f["cells"].values()) >= pack)
    ok("the plug holds once formed (all frames after the crossing are at rest)",
       all(f["vm"] == 0.0 for f in fb[firstBridge:]))

    # ---- 3. the passive twins --------------------------------------------
    rc, outp = run(CASE_DIR / "casePassive.dict")
    ok("passive twin exits 0", rc == 0, f"rc={rc}")
    ok("undeclared blockage prints no dEff/blk column",
       "dEff" not in outp and " blk " not in outp)

    # the declared arm with the blockage made inert is the twin verbatim
    inert = re.sub(r"blockage\s*\{.*?\n\}",
                   "blockage { packingFraction 1.0; depositionRate 0;"
                   " erosionRate 0; }", src, flags=re.S)
    p = TMP / "case_inert.dict"
    TMP.mkdir(parents=True, exist_ok=True)
    p.write_text(inert)
    rc, outi = run(p)
    ok("inert-declared arm exits 0", rc == 0, f"rc={rc}")

    def strip(text):
        return "\n".join(re.sub(r" (dEff|blk) [-+0-9.eE]+", "", l).rstrip()
                         for l in text.splitlines() if l.strip())

    ok("inert blockage is byte-identical to the undeclared twin",
       strip(outi) == strip(outp))

    # the 5.186 non-blocking deck is untouched by the new code
    ref = CASE_DIR.parent / "twoFluidHydrateNonBlocking" / "case.dict"
    rc, outr = run(ref)
    ok("5.186 deck still exits 0", rc == 0, f"rc={rc}")
    ok("5.186 deck prints no blockage column",
       " blk " not in outr)
    refTxt = "\n".join(l.split("//")[0] for l in ref.read_text().splitlines())
    refD = scal(re.search(r"pipe\s*\{(.*?)\}", refTxt, re.S).group(1),
                "diameter")
    refW = scal(re.search(r"wax\s*\{(.*?)\n\}", refTxt, re.S).group(1),
                "initial")
    fr5 = parse(outr)
    w5 = max(relerr(c["dEff"], refD * max(1.0 - refW, 1e-12) ** 0.5)
             for f in fr5 for c in f["cells"].values() if "dEff" in c)
    ok("5.186 dEff stays D sqrt(1 - wax) (worst %.3g)" % w5, w5 < 1e-9)

    # ---- 4. loud refusals --------------------------------------------------
    def variant(tag, new, old):
        txt = base.replace(old, new)
        q = TMP / ("case_reject_" + tag + ".dict")
        q.write_text(txt)
        rc, out = run(q)
        ok("refuses " + tag, rc != 0 and "error" in out.lower(),
           out.splitlines()[-1] if out else "no output")

    variant("packingZero", "packingFraction 0;", "packingFraction 0.6;")
    variant("packingNegative", "packingFraction -0.5;",
            "packingFraction 0.6;")
    variant("negativeDeposition", "depositionRate -1;",
            "depositionRate 2.0e-5;")
    variant("negativeErosion", "erosionRate -1;", "erosionRate 1.0e-4;")
    variant("unknownKey", "    packingX 1;\n", "    packingFraction 0.6;\n")
    noHyd = re.sub(r"hydrate\s*\{.*?\n\}\n", "", src, flags=re.S)
    q = TMP / "case_reject_noHydrate.dict"
    q.write_text(noHyd)
    rc, out = run(q)
    ok("refuses blockageWithoutHydrate",
       rc != 0 and "no hydrate{} is declared" in out,
       out.splitlines()[-1] if out else "no output")

    net = TMP / "case_reject_networkRow.dict"
    net.write_text("time { start 0; dt 0.1; endTime 0.1; reportEvery 1; }\n"
                   "network { pipes {} }\n"
                   "blockage { packingFraction 0.6; }\n")
    rc, out = run(net)
    ok("refuses blockageOnNetworkRow",
       rc != 0 and "blockage is a serial-row feature" in out,
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
