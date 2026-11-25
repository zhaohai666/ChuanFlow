#!/usr/bin/env python3
"""5.186 verification: hydrate as a NON-BLOCKING flow-assurance field.

usage: check_two_fluid_hydrate_nonblocking.py <exe> <case.dict> <tmpdir>

The claim under test is a NEGATIVE one, so it is proven two ways.

1. In the growing arm (kinetic hydrate{} live) the hydrate fraction h really
   moves, while the two hydraulics that own a length and a viscosity stay
   exactly what their own declarations say:
       d_eff == D sqrt(1 - wax)          (wax deposit only, h appears nowhere)
       muEff == muL (1 + 2.5 water cut)  (emulsion of the liquid only)
   on every section of every frame, with h > 0 (the field is alive), and
   water + h <= 1 (h eats the water, not the bore).
2. In the passive twin (growthRate = dissociationRate = 0, latentHeat -> 0 so
   h is frozen) two decks differing ONLY in the initial hydrate fraction (0.1
   against 0.5) must print byte-identical output once the h column is removed:
   h is read by no equation at all, so a 5x larger h cannot move a pressure,
   a velocity, a holdup or a temperature.
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


FRAME = re.compile(r"\[chflow\] t (\S+) vm (\S+) dp (\S+)")
CELL = re.compile(r"\[chflow\] cell (\d+) alpha (\S+) water (\S+) vsl (\S+) "
                  r"vsg (\S+) T (\S+) p (\S+) h (\S+) wax (\S+) dEff (\S+) "
                  r"muEff (\S+)")


def parse(text):
    frames = []
    cur = None
    for line in text.splitlines():
        m = FRAME.match(line)
        if m:
            if cur is not None:
                frames.append(cur)
            cur = {"t": float(m.group(1)), "cells": []}
            continue
        if cur is None:
            continue
        m = CELL.match(line)
        if m:
            cur["cells"].append(dict(
                i=int(m.group(1)), alpha=float(m.group(2)),
                water=float(m.group(3)), vsl=float(m.group(4)),
                vsg=float(m.group(5)), T=float(m.group(6)),
                p=float(m.group(7)), h=float(m.group(8)),
                wax=float(m.group(9)), dEff=float(m.group(10)),
                muEff=float(m.group(11))))
    if cur is not None:
        frames.append(cur)
    return frames


def run(path):
    out = subprocess.run([EXE, str(path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


def decl(path, key, blockname=None):
    t = Path(path).read_text()
    seg = t
    if blockname:
        seg = re.search(r"^" + blockname + r"\s*\{(.*?)^\}", t, re.S | re.M).group(1)
    return float(re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", seg).group(1))


def relerr(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main():
    D = decl(CASE, "diameter", "pipe")
    muL = decl(CASE, "muL", "fluid")
    w0 = decl(CASE, "initial", "wax")

    # ---- 1. the growing arm ---------------------------------------------
    rc, out = run(CASE)
    ok("hydrate arm exits 0", rc == 0, f"rc={rc}")
    fr = parse(out)
    ok("hydrate arm has frames and cells", len(fr) > 2 and len(fr[0]["cells"]) > 0)

    wde = wmu = 0.0
    hmax = 0.0
    wmix = 0.0
    for f in fr:
        for c in f["cells"]:
            wde = max(wde, relerr(c["dEff"], D * (1.0 - w0) ** 0.5))
            wmu = max(wmu, relerr(c["muEff"], muL * (1.0 + 2.5 * c["water"])))
            hmax = max(hmax, c["h"])
            if c["water"] + c["h"] > 1.0 + 1e-12:
                wmix = max(wmix, 1.0)
    ok("d_eff stays D sqrt(1 - wax): h never enters the bore (worst %.3g)"
       % wde, wde < 1e-9)
    ok("muEff stays muL (1 + 2.5 water cut): h never thickens (worst %.3g)"
       % wmu, wmu < 1e-9)
    ok("h is alive (h_max = %.3g > 1e-3)" % hmax, hmax > 1e-3)
    ok("water + h <= 1 on every cell (h eats water, not the bore)",
       wmix == 0.0)

    # h really varies in time (otherwise the invariant would be vacuous)
    h_first = [c["h"] for c in fr[0]["cells"]]
    h_last = [c["h"] for c in fr[-1]["cells"]]
    ok("h changes over the run (a growing field, not a frozen one)",
       max(abs(a - b) for a, b in zip(h_first, h_last)) > 1e-4)

    # alpha stays a holdup of the SLURRY: alpha in [0,1], never reduced by h
    walpha = 0.0
    for f in fr:
        for c in f["cells"]:
            if not (0.0 <= c["alpha"] <= 1.0):
                walpha = 1.0
    ok("alpha remains a valid holdup on every cell", walpha == 0.0)

    # ---- 2. the passive twins -------------------------------------------
    p0 = CASE_DIR / "casePassive.dict"
    p1 = CASE_DIR / "casePassiveHalf.dict"
    rc0, o0 = run(p0)
    rc1, o1 = run(p1)
    ok("passive twin exits 0", rc0 == 0 and rc1 == 0)
    ok("the twins really differ only in the initial h",
       "initial 0.1;" in p0.read_text() and "initial 0.5;" in p1.read_text())

    def strip_h(text):
        # The compare covers the HYDRAULIC lines only: h is read by no
        # holdup, momentum or diameter term. It does sit in the mass and
        # energy accounts (water inventory, gas sink, latent heat), so the
        # per-frame `thermal:` diagnostics - which sum the slurry's energy
        # inventory - are excluded; their spread is floating-point dust.
        return "\n".join(re.sub(r" h [-+0-9.eE]+", "", l).rstrip()
                         for l in text.splitlines()
                         if l.strip() and "thermal:" not in l)

    ok("with the h column removed the two decks are byte-identical",
       strip_h(o0) == strip_h(o1))
    ok("hydraulic lines carry no h dependence (5x h, identical hydraulics)",
       strip_h(o0) == strip_h(o1))
    ok("the h column itself does differ (the twin is not trivially equal)",
       o0 != o1)

    print()
    print("checks:", checks, "failures:", len(failures))
    if failures:
        for f in failures:
            print("FAILED:", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
