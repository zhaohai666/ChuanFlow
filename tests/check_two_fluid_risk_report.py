#!/usr/bin/env python3
"""5.187 (F-class) verification: the transient flow-assurance risk report of
the 1D two-fluid kernel (erosion / non-blocking hydrate / corrosion).

usage: check_two_fluid_risk_report.py <exe> <case.dict> <tmpdir>

The checker reads every declaration out of the case text, re-implements the
three laws locally and compares them to the report of every frame:

  * EROSION: rho_m from the frame's own printed (alpha, p, T) and the EOS
    rho_g = p/(gasR T); Ve = 0.3048 C / sqrt(rho_m/16.0184634); the worst
    vm/Ve section, its rho_m, the count of sections over the limit;
  * HYDRATE (the NON-BLOCKING reading): Teq from the declared
    Clausius-Clapeyron line (Tref, pRef, dH) at the reported section's own p
    and T, the subcooling Teq - T and its K band; and the invariant that the
    hydrate fraction h is in NO hydraulic term - the effective diameter stays
    D sqrt(1 - wax) (wax only) and the liquid viscosity stays
    muL (1 + 2.5 water-cut) (water only) on every section of every frame,
    while h itself really is non-zero (the kinetics are alive);
  * CORROSION: CO2/H2S partials at the highest section pressure, the CO2
    band and the 0.5 psia sour line.
A twin deck with the riskReport{} block deleted must reproduce every non-risk
line byte for byte (the report is read by no equation), and seven loud
refusals cover the grammar.
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
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


FRAME = re.compile(r"\[chflow\] t (\S+) vm (\S+) dp (\S+)")
CELL = re.compile(r"\[chflow\] cell (\d+) alpha (\S+) water (\S+) vsl (\S+) "
                  r"vsg (\S+) T (\S+) p (\S+) h (\S+) wax (\S+) dEff (\S+) "
                  r"muEff (\S+)")
EROS = re.compile(r"\[chflow\] risk erosion seg (-?\d+) vm (\S+) rhom (\S+) "
                  r"velim (\S+) ratio (\S+) over (\d+)")
HYD = re.compile(r"\[chflow\] risk hydrate seg (-?\d+) p (\S+) T (\S+) "
                 r"Teq (\S+) dTsub (\S+) grade (\d+)")
COR = re.compile(r"\[chflow\] risk corrosion co2 (\S+) h2s (\S+) band (\d+) "
                 r"sour (\d+)")


def parse(text):
    frames = []
    cur = None
    for line in text.splitlines():
        m = FRAME.match(line)
        if m:
            if cur is not None:
                frames.append(cur)
            cur = {"t": float(m.group(1)), "vm": float(m.group(2)),
                   "cells": {}, "erosion": None, "hydrate": None,
                   "corrosion": None}
            continue
        if cur is None:
            continue
        m = CELL.match(line)
        if m:
            cur["cells"][int(m.group(1))] = dict(
                alpha=float(m.group(2)), water=float(m.group(3)),
                vsl=float(m.group(4)), vsg=float(m.group(5)),
                T=float(m.group(6)), p=float(m.group(7)),
                h=float(m.group(8)), wax=float(m.group(9)),
                dEff=float(m.group(10)), muEff=float(m.group(11)))
            continue
        m = EROS.match(line)
        if m:
            cur["erosion"] = dict(seg=int(m.group(1)), vm=float(m.group(2)),
                                  rhom=float(m.group(3)),
                                  velim=float(m.group(4)),
                                  ratio=float(m.group(5)),
                                  over=int(m.group(6)))
            continue
        m = HYD.match(line)
        if m:
            cur["hydrate"] = dict(seg=int(m.group(1)), p=float(m.group(2)),
                                  T=float(m.group(3)), Teq=float(m.group(4)),
                                  dTsub=float(m.group(5)),
                                  grade=int(m.group(6)))
            continue
        m = COR.match(line)
        if m:
            cur["corrosion"] = dict(co2=float(m.group(1)),
                                    h2s=float(m.group(2)),
                                    band=int(m.group(3)),
                                    sour=int(m.group(4)))
    if cur is not None:
        frames.append(cur)
    return frames


def run(path):
    out = subprocess.run([EXE, str(path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


def scal(text, key):
    return float(re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", text).group(1))


def read_decls(path):
    t = Path(path).read_text()

    def block(name):
        m = re.search(r"^" + name + r"\s*\{(.*?)^\}", t, re.S | re.M)
        return m.group(1) if m else ""

    d = dict(
        gasR=scal(block("fluid"), "gasR"),
        rhoL=scal(block("fluid"), "rhoL"),
        muL=scal(block("fluid"), "muL"),
        D=scal(block("pipe"), "diameter"),
        C=scal(block("riskReport"), "erosionC"),
        Tref=scal(block("riskReport") and re.search(r"hydrate\s*\{(.*?)\}",
                                                    block("riskReport"),
                                                    re.S).group(1), "Tref"),
        pRef=scal(re.search(r"hydrate\s*\{(.*?)\}", block("riskReport"),
                            re.S).group(1), "pRef"),
        dH=scal(re.search(r"hydrate\s*\{(.*?)\}", block("riskReport"),
                          re.S).group(1), "dH"),
        co2=scal(re.search(r"corrosion\s*\{(.*?)\}", block("riskReport"),
                           re.S).group(1), "co2"),
        h2s=scal(re.search(r"corrosion\s*\{(.*?)\}", block("riskReport"),
                           re.S).group(1), "h2s"),
        co2Low=scal(re.search(r"corrosion\s*\{(.*?)\}", block("riskReport"),
                              re.S).group(1), "co2Low"),
        co2High=scal(re.search(r"corrosion\s*\{(.*?)\}", block("riskReport"),
                               re.S).group(1), "co2High"),
        wax0=scal(block("wax"), "initial"),
    )
    for key in ("hydrateLowK", "hydrateMediumK"):
        mm = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", block("riskReport"))
        d[key] = float(mm.group(1)) if mm else (3.0 if key == "hydrateLowK"
                                                 else 8.0)
    d["grade"] = lambda dts: 0 if dts < 0 else (1 if dts < d["hydrateLowK"]
                                                else (2 if dts < d["hydrateMediumK"]
                                                      else 3))
    d["ve"] = lambda rho_m: (0.3048 * d["C"]
                             / math.sqrt(rho_m / 16.0184634))
    d["teq"] = lambda p: 1.0 / (1.0 / d["Tref"]
                                - 8.314462618 * math.log(p / d["pRef"]) / d["dH"])
    return d


def main():
    d = read_decls(CASE)
    rc, out = run(CASE)
    ok("risk arm exits 0", rc == 0, f"rc={rc}")
    fr = parse(out)
    ok("every frame carries all three report lines",
       all(f["erosion"] and f["hydrate"] and f["corrosion"] for f in fr))

    def rho_m(cell):
        rg = cell["p"] / (d["gasR"] * cell["T"])
        return cell["alpha"] * d["rhoL"] + (1.0 - cell["alpha"]) * rg

    # ---- erosion --------------------------------------------------------
    we_seg = we_rho = we_ve = we_ratio = 0.0
    we_over = 0.0
    for f in fr:
        e = f["erosion"]
        cells = f["cells"]
        over = 0
        for i, c in cells.items():
            rm = rho_m(c)
            ratio = abs(f["vm"]) / d["ve"](rm)
            if ratio > 1.0:
                over += 1
        we_over = max(we_over, 1.0 if over != e["over"] else 0.0)
        # the worst section: the largest vm/Ve (vm is global, so the smallest
        # Ve); ties resolve to the first section, as the kernel scans upward
        i, c = min(cells.items(), key=lambda kv: (d["ve"](rho_m(kv[1])), kv[0]))
        rm = rho_m(c)
        we_seg = max(we_seg, 1.0 if i != e["seg"] else 0.0)
        we_rho = max(we_rho, relerr(e["rhom"], rm))
        we_ve = max(we_ve, relerr(e["velim"], d["ve"](rm)))
        we_ratio = max(we_ratio, relerr(e["ratio"], abs(f["vm"]) / d["ve"](rm)))
    ok("erosion worst section is the smallest-Ve section", we_seg == 0.0)
    ok("erosion rhom is the mixture density there (worst %.3g)" % we_rho,
       we_rho < 1e-9)
    ok("erosion Ve == 0.3048 C / sqrt(rho/16.0184634) (worst %.3g)" % we_ve,
       we_ve < 1e-9)
    ok("erosion ratio == vm/Ve (worst %.3g)" % we_ratio, we_ratio < 1e-9)
    ok("erosion over-count matches the sections above Ve", we_over == 0.0)
    ok("the run really crosses the limit (over-count > 0)",
       max(f["erosion"]["over"] for f in fr) > 0)
    ok("erosion prices a MIXTURE density, not rhoL",
       all(f["erosion"]["rhom"] < d["rhoL"] - 1e-9 for f in fr))

    # ---- hydrate: the readout, and the non-blocking invariants -----------
    wt = wgrade = wpres = 0.0
    wseg = 0.0
    for f in fr:
        h = f["hydrate"]
        c = f["cells"][h["seg"]]
        wseg = max(wseg, 1.0 if abs(c["p"] - h["p"]) > 1e-6
                   or abs(c["T"] - h["T"]) > 1e-9 else 0.0)
        teq = d["teq"](c["p"])
        wt = max(wt, relerr(h["Teq"], teq))
        wgrade = max(wgrade, 1.0 if h["grade"] != d["grade"](teq - c["T"])
                     else 0.0)
        # the reported pair must be the frame's worst (largest) subcooling
        worst = max(cells_teq(f["cells"], d) , key=lambda kv: kv[1])
        wpres = max(wpres, 1.0 if abs(worst[0] - h["seg"]) > 0 else 0.0)
    ok("hydrate line seeds the reported section's own p and T", wseg == 0.0)
    ok("hydrate Teq == the declared CC line (worst %.3g)" % wt, wt < 1e-9)
    ok("hydrate grade follows the K bands", wgrade == 0.0)
    ok("hydrate reports the largest-subcooling section", wpres == 0.0)
    ok("the hydrate channel reaches a risk band (grade >= 2 somewhere)",
       max(f["hydrate"]["grade"] for f in fr) >= 2)

    # NON-BLOCKING: h is in neither the diameter nor the viscosity, on every
    # cell of every frame, while h itself is alive.
    wde = wmu = 0.0
    hmax = 0.0
    for f in fr:
        for c in f["cells"].values():
            d_expect = d["D"] * math.sqrt(max(1.0 - d["wax0"], 1e-12))
            wde = max(wde, relerr(c["dEff"], d_expect))
            wmu = max(wmu, relerr(c["muEff"], d["muL"] * (1.0 + 2.5 * c["water"])))
            hmax = max(hmax, c["h"])
    ok("d_eff == D sqrt(1 - wax) with hydrate declared (worst %.3g)" % wde,
       wde < 1e-9)
    ok("muEff == muL (1 + 2.5 water cut) with hydrate declared (worst %.3g)"
       % wmu, wmu < 1e-9)
    ok("the hydrate fraction itself is alive (h_max = %.3g)" % hmax,
       hmax > 1e-3)
    ok("h never reaches the diameter or the viscosity (d_eff, muEff invariant)",
       wde < 1e-9 and wmu < 1e-9)

    # ---- corrosion -------------------------------------------------------
    wco2 = wh2s = 0.0
    wband = wps = 0.0
    for f in fr:
        c = f["corrosion"]
        pmax = max(cell["p"] for cell in f["cells"].values())
        wps = max(wps, 1.0 if pmax <= 0 else 0.0)
        wco2 = max(wco2, relerr(c["co2"], d["co2"] * pmax))
        wh2s = max(wh2s, relerr(c["h2s"], d["h2s"] * pmax))
        band = 0 if d["co2"] * pmax < d["co2Low"] else (
            1 if d["co2"] * pmax < d["co2High"] else 2)
        wband = max(wband, 1.0 if band != c["band"] else 0.0)
        wband = max(wband, 1.0 if c["sour"] != (d["h2s"] * pmax >= 3447.0)
                    else 0.0)
    ok("corrosion prices the highest section pressure", wps == 0.0)
    ok("CO2 partial == co2 * p_max (worst %.3g)" % wco2, wco2 < 1e-9)
    ok("H2S partial == h2s * p_max (worst %.3g)" % wh2s, wh2s < 1e-9)
    ok("CO2 band and the 0.5 psia sour line", wband == 0.0)

    # ---- the report is read by NO equation -------------------------------
    twin = CASE_DIR / "caseNoRisk.dict"
    rc2, out2 = run(twin)
    ok("twin exits 0", rc2 == 0, f"rc={rc2}")
    ok("twin prints no risk token", "risk " not in out2 and "risk\t" not in out2)
    stripped = "\n".join(l for l in out.splitlines() if "risk " not in l)
    ok("declared report is byte-identical to the undeclared twin",
       stripped.rstrip() == out2.rstrip())

    # ---- loud refusals ---------------------------------------------------
    TMP.mkdir(parents=True, exist_ok=True)
    base = CASE.read_text()

    def variant(tag, new, old):
        txt = base.replace(old, new)
        p = TMP / ("case_reject_" + tag + ".dict")
        p.write_text(txt)
        rc, out = run(p)
        ok("refuses " + tag, rc != 0 and "error" in out.lower(),
           out.splitlines()[-1] if out else "no output")

    variant("unknownKey", "    erosionX 100;\n",
            "    erosionC 100;\n")
    variant("partialHydrate", "hydrate { Tref 313.4; pRef 3.0e6; }",
            "hydrate { Tref 313.4; pRef 3.0e6; dH 47256.4; }")
    variant("hydrateWithoutEos", "", "gasR 470; ")
    variant("co2OutOfRange", "co2 1.5;", "co2 0.08;")
    variant("bandsNotAscending", "co2Low 1.0e5; co2High 1.0e4;",
            "co2Low 1.0e4; co2High 1.0e5;")
    variant("zeroErosionC", "erosionC 0;", "erosionC 100;")
    variant("hydrateBandsReversed", "hydrateLowK 9; hydrateMediumK 4;",
            "erosionC 100;")

    # 5.187: the report is a serial-row feature. Declaring it on a branch
    # network row must be refused by name (a declaration that runs unreported
    # physics is a bug, not a default), so the probe is a network deck whose
    # only content past the time block is the offending block.
    net = (TMP / "case_reject_networkRow.dict")
    net.write_text("time { start 0; dt 0.1; endTime 0.1; reportEvery 1; }\n"
                   "network { pipes {} }\n"
                   "riskReport { erosionC 100; }\n")
    rc, out = run(net)
    ok("refuses riskReportOnNetworkRow",
       rc != 0 and "riskReport is a serial-row feature" in out,
       out.splitlines()[-1] if out else "no output")

    print()
    print("checks:", checks, "failures:", len(failures))
    if failures:
        for f in failures:
            print("FAILED:", f)
        return 1
    return 0


def cells_teq(cells, d):
    out = []
    for i, c in cells.items():
        if c["p"] > 0:
            out.append((i, d["teq"](c["p"]) - c["T"]))
    return out


if __name__ == "__main__":
    sys.exit(main())
