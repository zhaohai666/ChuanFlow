#!/usr/bin/env python3
"""D2 increment 1 (5.96): wax precipitation / dissolution, independently
re-derived.

The kernel's waxStep integrates, per section, the first-order relaxation of
the deposited wax w onto the van't Hoff equilibrium w_eq = max(0, z - xsat):
    xsat(T) = exp[-(dHf/R)(1/T - 1/Tf)]
    w_new   = (w/dt + k w_eq) / (1/dt + k),  k = kg if w<w_eq else kd
waxStep runs BEFORE thermalStep in each step, so the wax increment printed at
frame n was integrated with the temperature reported at frame n-1. This
checker reads the emitted per-frame T and wax from ALL frames (the case uses
reportEvery 1) and reproduces that recursion with nothing but math.exp and the
declared constants -- it never calls a kernel constitutive function.

Checks:
  1. per-frame, per-cell implicit-Euler wax recursion from the previous
     frame's reported T and w (the core independent witness);
  2. the closed-form equilibrium xsat and w_eq = max(0, z - xsat);
  3. the two arms: warm cells (xsat >= z) keep w = 0, cold cells deposit;
     the precipitation boundary sits where xsat(T) = z, not at Tf;
  4. the late frames converge to the closed-form equilibrium deposit;
  5. bounded w in [0, z] every frame;
  6. the latent-heat twin raises the deposit-zone temperature, while a no-wax
     twin is bit-for-bit identical in temperature to the latentHeat=0 run;
  7. six loud rejections (bad Tf, dHf, waxFraction, rate, initial, no thermal).

Usage: check_two_fluid_wax.py <exe> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

R = 8.314462618
FAILS = []


def check(name, ok, detail=""):
    if not ok:
        FAILS.append(f"{name} {detail}".strip())


def near(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def parse_case(path):
    src = open(path).read()

    def key(block, name, default=None):
        m = re.search(r"\b" + re.escape(name) + r"\s+([-\d.eE+]+)\s*;", src)
        return float(m.group(1)) if m else default

    wax = {}
    wb = re.search(r"\bwax\s*\{(.*?)\n\}", src, re.S)
    if wb:
        for k in ("meltingTemperature", "fusionEnthalpy", "waxFraction",
                  "growthRate", "dissolutionRate", "latentHeat", "initial"):
            m = re.search(r"\b" + k + r"\s+([-\d.eE+]+)\s*;", wb.group(1))
            if m:
                wax[k] = float(m.group(1))
    dt = key(None, "dt")
    # wall profile: Tw(drop) = surface + gradient * (inletDepth - drop)
    tb = re.search(r"\bthermal\s*\{(.*?)\n\}", src, re.S)
    surface = grad = depth = None
    if tb:
        for k, dst in (("surfaceTemperature", "s"), ("gradient", "g"),
                       ("inletDepth", "d")):
            m = re.search(r"\b" + k + r"\s+([-\d.eE+]+)\s*;", tb.group(1))
            if m:
                v = float(m.group(1))
                if dst == "s":
                    surface = v
                elif dst == "g":
                    grad = v
                else:
                    depth = v
    # section-centre drops from the pipe geometry (vertical line)
    pb = re.search(r"\bpipe\s*\{(.*?)\n\}", src, re.S)
    length = float(re.search(r"\blength\s+([-\d.eE+]+)", pb.group(1)).group(1))
    nsec = int(float(re.search(r"\bsections\s+([-\d.eE+]+)", pb.group(1)).group(1)))
    dx = length / nsec
    drops = [(i + 0.5) * dx for i in range(nsec)]  # inlet at the bottom, drop grows
    return src, wax, dt, (surface, grad, depth), drops, nsec


def run(exe, case):
    r = subprocess.run([exe, case], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"case failed ({r.returncode}): " + (r.stdout + r.stderr)[-400:])
    frames, cur = [], None
    for line in r.stdout.splitlines():
        line = line.replace("[chflow] ", "")
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 2:
            cur = {"t": float(tok[1]), "cells": {}}
            frames.append(cur)
        elif tok[0] == "cell" and cur is not None:
            i = int(tok[1])
            d = {}
            for j in range(2, len(tok) - 1):
                if tok[j] in ("T", "wax"):
                    d[tok[j]] = float(tok[j + 1])
            cur["cells"][i] = d
    return frames


def xsat(T, Tf, dHf):
    return math.exp(-(dHf / R) * (1.0 / T - 1.0 / Tf))


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_two_fluid_wax.py <exe> <case.dict>")
    exe, case = sys.argv[1], sys.argv[2]
    src, wx, dt, wall, drops, nsec = parse_case(case)
    Tf, dHf = wx["meltingTemperature"], wx["fusionEnthalpy"]
    z, kg, kd = wx["waxFraction"], wx["growthRate"], wx["dissolutionRate"]
    check("wax case parsed", bool(wx) and dt > 0 and Tf > 0 and dHf > 0,
          f"wx={wx} dt={dt}")

    frames = run(exe, case)
    check("all frames emitted", len(frames) >= 2, f"{len(frames)} frames")
    check("every frame has all cells",
          all(len(f["cells"]) == nsec for f in frames), "cell count mismatch")

    # per-frame wall temperatures (closed form, independent of the kernel)
    surface, grad, depth = wall
    check("geothermal wall parsed", None not in (surface, grad, depth), str(wall))
    wallT = [surface + grad * (depth - d) for d in drops]

    # ---- 1 & 5: implicit-Euler recursion from the previous frame's T ----
    worst_rec = 0.0
    for n in range(1, len(frames)):
        prev, now = frames[n - 1], frames[n]
        for i in range(nsec):
            w0, w1 = prev["cells"][i]["wax"], now["cells"][i]["wax"]
            T0 = prev["cells"][i]["T"]
            xs = xsat(T0, Tf, dHf)
            drive = z - xs
            weq = max(drive, 0.0)
            k = kg if drive > 0 else kd
            pred = min(max((w0 / dt + k * weq) / (1.0 / dt + k), 0.0), z)
            worst_rec = max(worst_rec, abs(w1 - pred))
            check(f"frame {n} cell {i} wax bounded", -1e-12 <= w1 <= z + 1e-12,
                  f"w={w1}")
    check("wax recursion reproduced (1e-9)", worst_rec < 1e-9,
          f"worst {worst_rec:.3e}")

    # ---- 2 & 3: closed-form equilibrium and the two arms ----
    last = frames[-1]
    warm, cold = [], []
    for i in range(nsec):
        T = last["cells"][i]["T"]
        xs = xsat(T, Tf, dHf)
        weq = max(z - xs, 0.0)
        check(f"cell {i} equilibrium xsat positive", xs > 0, f"xsat={xs}")
        if xs >= z:
            warm.append(i)
            check(f"warm cell {i} stays clean", abs(last["cells"][i]["wax"]) < 1e-10,
                  f"w={last['cells'][i]['wax']}")
        else:
            cold.append(i)
    check("both arms visited", len(warm) > 0 and len(cold) > 0,
          f"warm={len(warm)} cold={len(cold)}")
    # precipitation boundary is xsat = z, i.e. T* = 1/(1/Tf - R ln z / dHf)
    Tstar = 1.0 / (1.0 / Tf - (R * math.log(z)) / dHf)
    check("precipitation boundary is xsat=z not Tf", abs(Tstar - Tf) > 10.0,
          f"T*={Tstar:.3f} Tf={Tf}")

    # cold cells deposit monotonically and toward w_eq
    for i in cold:
        seq = [f["cells"][i]["wax"] for f in frames]
        mono = all(seq[j + 1] >= seq[j] - 1e-12 for j in range(len(seq) - 1))
        check(f"cold cell {i} monotone deposit", mono, "non-monotone wax")
        T = last["cells"][i]["T"]
        weq = max(z - xsat(T, Tf, dHf), 0.0)
        near(f"cold cell {i} reaches equilibrium", seq[-1], weq, 1e-3, 1e-6)

    # wall profile: late-frame temperature equals the geothermal line
    worst_T = max(abs(last["cells"][i]["T"] - wallT[i]) for i in range(nsec))
    check("late temperature matches geothermal wall (1e-6)", worst_T < 1e-6,
          f"worst {worst_T:.3e}")

    # ---- 6: latent heat twin and no-wax twin ----
    base_dir = os.path.dirname(case)
    latent_case = os.path.join(base_dir, "case_latent.dict")
    nowax_case = os.path.join(base_dir, "case_nowax.dict")
    latent = run(exe, latent_case)
    nowax = run(exe, nowax_case)
    # latentHeat=0 base vs no-wax twin: temperature bit-for-bit identical
    mx_diff = 0.0
    for fb, fn in zip(frames, nowax):
        for i in range(nsec):
            mx_diff = max(mx_diff, abs(fb["cells"][i]["T"] - fn["cells"][i]["T"]))
    check("no-wax twin temperature is bit-identical", mx_diff == 0.0,
          f"max diff {mx_diff:.3e}")
    # latent case raises temperature where wax precipitates
    rise = 0.0
    for fb, fl in zip(frames, latent):
        for i in cold:
            rise = max(rise, fl["cells"][i]["T"] - fb["cells"][i]["T"])
    check("latent heat raises deposit-zone temperature", rise > 0.1,
          f"max rise {rise:.4f} K")

    # ---- 7: six loud rejections ----
    rejects = [
        ("negative meltingTemperature",
         src.replace("meltingTemperature 323.75", "meltingTemperature -1"),
         "positive meltingTemperature"),
        ("negative fusionEnthalpy",
         src.replace("fusionEnthalpy 54890", "fusionEnthalpy -1"),
         "positive meltingTemperature"),
        ("waxFraction above 1",
         src.replace("waxFraction 0.5", "waxFraction 1.5"),
         "waxFraction outside"),
        ("negative growthRate",
         src.replace("growthRate 0.01", "growthRate -0.01"),
         "negative wax rate"),
        ("initial above 1",
         src.replace("initial 0;", "initial 1.5;"),
         "initial wax fraction outside"),
        ("wax without thermal",
         re.sub(r"\nthermal\s*\{.*?\n\}\n", "\n", src, count=1, flags=re.S),
         "needs a thermal"),
    ]
    for what, txt, needle in rejects:
        fd, p = tempfile_path()
        with open(p, "w") as fh:
            fh.write(txt)
        rr = subprocess.run([exe, p], capture_output=True, text=True)
        os.remove(p)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"reject: {what}", ok, f"rc={rr.returncode}")

    if FAILS:
        print("\n".join("FAIL: " + m for m in FAILS))
        raise SystemExit(f"{len(FAILS)} wax check(s) failed")
    print(f"wax: {len(frames)} frames x {nsec} cells, {len(warm)} warm / "
          f"{len(cold)} cold, worst recursion {worst_rec:.3e}, "
          f"T* {Tstar:.3f} K, latent rise {rise:.3f} K - all checks passed")


def tempfile_path():
    import tempfile
    fd, p = tempfile.mkstemp(suffix=".dict")
    os.close(fd)
    return fd, p


if __name__ == "__main__":
    main()
