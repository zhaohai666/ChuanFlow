#!/usr/bin/env python3
"""5.230 verification: time-varying ambient temperature for the wall.

usage: check_two_fluid_ambient_schedule.py <exe> <case.dict> <tmpdir>

The 5.144/5.188 wall is a constant or a fixed geothermal profile; declaring
ambientSchedule{ makes the ambient level follow a piecewise-linear schedule.
The checker re-implements both halves independently:

  * the REPORTED wallT column must equal the closed-form interpolation of
    (times, ambientT) at the frame time (held outside the window, the exact
    segment search of the kernel: last node at or before t, linear blend);
  * the temperature field must relax onto that MOVING wall: the sealed,
    liquid-full line is one lumped capacity per section,
        cap = A dx rhoL cL,  kw = k' dx
        T_n = (T_{n-1} cap/dt + kw wallT_n) / (cap/dt + kw)
    (implicit Euler, no flow, no EOS) - marched frame by frame to 1e-9, so
    the schedule is proven to drive the physics, not just to print.
  * PASSIVE twin (casePassive.dict, ambientSchedule{} undeclared) prints no
    wallT column and obeys the old constant-wall closed form; a constant
    schedule (ambientT all equal to the declared wallTemperature) is
    byte-identical to the twin once the wallT column is stripped; the 5.144
    relax deck itself is re-run as the untouched witness.
  * loud refusals: times not strictly ascending, times/ambientT length
    mismatch, one point only, unknown key, schedule without thermal{},
    schedule on a network row.
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


def nums(block, key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", block)
    return [float(x) for x in m.group(1).split()]


def interp(t, ts, tv):
    # the kernel's exact piecewise-linear ambient, held outside the window
    if t <= ts[0]:
        return tv[0]
    if t >= ts[-1]:
        return tv[-1]
    k = 0
    while k + 1 < len(ts) - 1 and ts[k + 1] <= t:
        k += 1
    f = (t - ts[k]) / (ts[k + 1] - ts[k])
    return tv[k] + f * (tv[k + 1] - tv[k])


def strip_wallt(text):
    return "\n".join(re.sub(r" wallT [-+0-9.eE]+", "", l).rstrip()
                     for l in text.splitlines() if l.strip())


def main():
    base = CASE.read_text()
    src = "\n".join(l.split("//")[0] for l in base.splitlines())
    asrc = re.search(r"ambientSchedule\s*\{(.*?)\n\}", src, re.S).group(1)
    ts = nums(asrc, "times")
    tv = nums(asrc, "ambientT")
    dt = scal(re.search(r"time\s*\{(.*?)\}", src, re.S).group(1), "dt")
    pipe = re.search(r"pipe\s*\{(.*?)\}", src, re.S).group(1)
    L = scal(pipe, "length")
    D = scal(pipe, "diameter")
    nsec = int(scal(pipe, "sections"))
    fl = re.search(r"fluid\s*\{(.*?)\}", src, re.S).group(1)
    rhoL = scal(fl, "rhoL")
    cL = scal(fl, "cL")
    th = re.search(r"thermal\s*\{(.*?)\}", src, re.S).group(1)
    kwall = scal(th, "wallHeatTransfer")
    wallT0 = scal(th, "wallTemperature")
    T0 = scal(re.search(r"initial\s*\{(.*?)\}", src, re.S).group(1),
              "temperature")
    alpha0 = scal(re.search(r"initial\s*\{(.*?)\}", src, re.S).group(1),
                  "alpha")

    # ---- 1. the schedule arm ---------------------------------------------
    rc, out = run(CASE)
    ok("schedule arm exits 0", rc == 0, f"rc={rc}")
    fr = parse(out)
    ok("schedule arm reports wallT every frame",
       len(fr) > 2 and all("wallT" in c for f in fr
                           for c in f["cells"].values()))

    # the reported wallT is the closed-form interpolation at the frame time
    wt = 0.0
    for f in fr:
        want = interp(f["t"], ts, tv)
        for c in f["cells"].values():
            wt = max(wt, relerr(c["wallT"], want))
    ok("wallT == the piecewise-linear schedule (worst %.3g)" % wt, wt < 1e-9)
    ok("the ambient really moves (wallT spans %.3g K)"
       % (max(c["wallT"] for f in fr for c in f["cells"].values())
          - min(c["wallT"] for f in fr for c in f["cells"].values())),
       max(c["wallT"] for f in fr for c in f["cells"].values())
       - min(c["wallT"] for f in fr for c in f["cells"].values()) > 1.0)

    # the temperature field relaxes onto the moving wall: the sealed,
    # liquid-full line is an independent lumped capacity per section
    dx = L / nsec
    area = math.pi * D * D / 4.0
    rhoG = scal(fl, "rhoG")
    cG = scal(fl, "cG")
    cap = area * dx * (alpha0 * rhoL * cL + (1.0 - alpha0) * rhoG * cG)
    kw = kwall * dx
    wT = 0.0
    for n in range(1, len(fr)):
        for i, c in fr[n]["cells"].items():
            prev = fr[n - 1]["cells"][i]["T"]
            tw = interp(fr[n]["t"], ts, tv)
            want = (prev * cap / dt + kw * tw) / (cap / dt + kw)
            wT = max(wT, relerr(c["T"], want))
    ok("T obeys the implicit-Euler march onto the moving wall (worst %.3g)"
       % wT, wT < 1e-9)
    ok("the schedule drives the field (final T != the constant-wall answer)",
       abs(fr[-1]["cells"][0]["T"] - (wallT0 + (T0 - wallT0)
                                      / (1.0 + dt * kw / cap)
                                      ** (len(fr) - 1))) > 1.0)

    # ---- 2. the passive twins --------------------------------------------
    rc, outp = run(CASE_DIR / "casePassive.dict")
    ok("passive twin exits 0", rc == 0, f"rc={rc}")
    ok("undeclared schedule prints no wallT column", "wallT" not in outp)
    fp = parse(outp)
    wP = 0.0
    for n in range(1, len(fp)):
        for i, c in fp[n]["cells"].items():
            prev = fp[n - 1]["cells"][i]["T"]
            want = (prev * cap / dt + kw * wallT0) / (cap / dt + kw)
            wP = max(wP, relerr(c["T"], want))
    ok("passive twin keeps the 5.144 constant-wall closed form (worst %.3g)"
       % wP, wP < 1e-9)

    # a constant schedule at the declared wall value is the old path verbatim
    const = re.sub(r"ambientT\s*\([^)]*\)",
                   "ambientT ( %g %g %g )" % (wallT0, wallT0, wallT0),
                   src)
    TMP.mkdir(parents=True, exist_ok=True)
    p = TMP / "case_const.dict"
    p.write_text(const)
    rc, outc = run(p)
    ok("constant-schedule arm exits 0", rc == 0, f"rc={rc}")
    ok("a constant schedule is byte-identical to the undeclared twin",
       strip_wallt(outc) == strip_wallt(outp))

    # the 5.144 relax deck is untouched by the new code
    ref = CASE_DIR.parent / "twoFluidThermalRelax" / "case.dict"
    rc, outr = run(ref)
    ok("5.144 deck still exits 0", rc == 0, f"rc={rc}")
    ok("5.144 deck prints no wallT column", "wallT" not in outr)

    # ---- 3. loud refusals --------------------------------------------------
    def variant(tag, new, old):
        txt = base.replace(old, new)
        q = TMP / ("case_reject_" + tag + ".dict")
        q.write_text(txt)
        rc, out = run(q)
        ok("refuses " + tag, rc != 0 and "error" in out.lower(),
           out.splitlines()[-1] if out else "no output")

    variant("timesNotAscending", "times    ( 0 30 20 );",
            "times    ( 0 30 60 );")
    variant("timesEqual", "times    ( 0 30 30 );",
            "times    ( 0 30 60 );")
    variant("sizeMismatch", "ambientT ( 280 310 );",
            "ambientT ( 280 310 290 );")
    variant("onePoint", "times    ( 0 );\n    ambientT ( 280 );",
            "times    ( 0 30 60 );\n    ambientT ( 280 310 290 );")
    variant("unknownKey", "    ambientX ( 1 );\n", "    times    ( 0 30 60 );\n")
    noTh = re.sub(r"thermal\s*\{.*?\n\}\n", "", src, flags=re.S)
    q = TMP / "case_reject_noThermal.dict"
    q.write_text(noTh)
    rc, out = run(q)
    ok("refuses scheduleWithoutThermal",
       rc != 0 and "no thermal{} is declared" in out,
       out.splitlines()[-1] if out else "no output")

    net = TMP / "case_reject_networkRow.dict"
    net.write_text("time { start 0; dt 0.1; endTime 0.1; reportEvery 1; }\n"
                   "network { pipes {} }\n"
                   "ambientSchedule { times ( 0 1 ); ambientT ( 1 2 ); }\n")
    rc, out = run(net)
    ok("refuses ambientScheduleOnNetworkRow",
       rc != 0 and "ambientSchedule is a serial-row feature" in out,
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
