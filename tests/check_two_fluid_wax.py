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
    if len(sys.argv) not in (3, 4):
        raise SystemExit("usage: check_two_fluid_wax.py <exe> <case.dict> [mode]")
    exe, case = sys.argv[1], sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) == 4 else "local"
    if mode == "transport":
        return check_transport(exe, case)
    if mode != "local":
        raise SystemExit(f"unknown wax mode {mode}")
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


def check_transport(exe, case):
    """5.97 (D2 increment 2): deposited-wax transport in a rate-driven line.

    wallis slip 0 and alpha = 1 make the liquid flux the constant F = vsl;
    the deposit alpha*w then upwinds linearly. Each step runs transport FIRST,
    then the local van't Hoff reaction, so the frame-n deposit is reproduced
    from frame n-1 by:  w_trans = w + dt/dx*F*(w_up - w);
                        w_new   = (w_trans/dt + k*w_eq)/(1/dt + k).
    w_eq = max(0, z - xsat(T)); the case is isothermal so T is constant.
    """
    src = open(case).read()

    def scalar(name):
        m = re.search(r"\b" + name + r"\s+([-\d.eE+]+)\s*;", src)
        return float(m.group(1)) if m else None

    dt = scalar("dt")
    db = re.search(r"\binlet\s*\{(.*?)\}", src, re.S)
    vsl = float(re.search(r"\bvsl\s+([-\d.eE+]+)", db.group(1)).group(1))
    win = float(re.search(r"\bwaxDepositFraction\s+([-\d.eE+]+)", db.group(1)).group(1))
    pb = re.search(r"\bpipe\s*\{(.*?)\n\}", src, re.S)
    length = float(re.search(r"\blength\s+([-\d.eE+]+)", pb.group(1)).group(1))
    nsec = int(float(re.search(r"\bsections\s+([-\d.eE+]+)", pb.group(1)).group(1)))
    dx = length / nsec
    wb = re.search(r"\bwax\s*\{(.*?)\n\}", src, re.S)
    Tf = float(re.search(r"\bmeltingTemperature\s+([-\d.eE+]+)", wb.group(1)).group(1))
    dHf = float(re.search(r"\bfusionEnthalpy\s+([-\d.eE+]+)", wb.group(1)).group(1))
    z = float(re.search(r"\bwaxFraction\s+([-\d.eE+]+)", wb.group(1)).group(1))
    kg = float(re.search(r"\bgrowthRate\s+([-\d.eE+]+)", wb.group(1)).group(1))
    kd = float(re.search(r"\bdissolutionRate\s+([-\d.eE+]+)", wb.group(1)).group(1))
    w0 = float(re.search(r"\binitial\s+([-\d.eE+]+)\s*;", wb.group(1)).group(1))
    tb = re.search(r"\bthermal\s*\{(.*?)\n\}", src, re.S)
    Twall = float(re.search(r"\bwallTemperature\s+([-\d.eE+]+)", tb.group(1)).group(1))

    check("transport case parsed",
          None not in (dt, vsl, win, Tf, dHf, z, kg, kd) and dx > 0,
          f"dt={dt} vsl={vsl} nsec={nsec}")
    xs = xsat(Twall, Tf, dHf)
    weq = max(z - xs, 0.0)
    # the isothermal reaction drive direction is fixed
    k = kg if z - xs > 0 else kd
    F = vsl  # alpha = 1, wallis slip 0
    check("transport CFL", dt * F / dx <= 1.0 + 1e-12,
          f"dt*F/dx={dt*F/dx}")

    frames = run(exe, case)
    check("all frames emitted", len(frames) >= 2, f"{len(frames)} frames")
    check("every frame has all cells",
          all(len(f["cells"]) == nsec for f in frames), "cell count")
    # frame 0 is the initial state (no transport/reaction step yet)
    for i in range(nsec):
        near(f"frame 0 cell {i} initial deposit",
             frames[0]["cells"][i]["wax"], w0, 0, 1e-12)

    # reproduce transport-then-reaction from every frame to the next
    worst = 0.0
    for n in range(1, len(frames)):
        prev, now = frames[n - 1], frames[n]
        w = [prev["cells"][i]["wax"] for i in range(nsec)]
        trans = []
        for i in range(nsec):
            up = win if i == 0 else w[i - 1]
            trans.append(w[i] + dt / dx * F * (up - w[i]))
        for i in range(nsec):
            pred = (trans[i] / dt + k * weq) / (1.0 / dt + k)
            worst = max(worst, abs(now["cells"][i]["wax"] - pred))
            check(f"frame {n} cell {i} deposit bounded",
                  -1e-12 <= now["cells"][i]["wax"] <= z + 1e-12, "")
    check("transport+reaction recursion reproduced (1e-9)", worst < 1e-9,
          f"worst {worst:.3e}")

    # late frames approach the closed-form equilibrium from below everywhere
    last = frames[-1]
    for i in range(nsec):
        check(f"late cell {i} approaches equilibrium",
              w0 - 1e-9 <= last["cells"][i]["wax"] <= weq + 1e-9,
              f"w={last['cells'][i]['wax']} weq={weq}")
    check("equilibrium is above the inlet injection", weq > win,
          f"weq={weq} win={win}")

    # sealed twin: same wax{} but zero flux -> pure 5.96 local reaction, no
    # transport, no inflow. Every cell is identical and matches the standalone
    # implicit-Euler reaction w=(w/dt+k weq)/(1/dt+k). This also witnesses
    # that adding transport left the old local path bit-for-bit unchanged.
    sealed = os.path.join(os.path.dirname(case), "case_sealed.dict")
    twin = src.replace("type phaseRates; vsl 0.5; vsg 0; waterFraction 0;",
                       "type sealed;")
    twin = twin.replace("type fixedPressure; p 1.0e5;", "type sealed;")
    with open(sealed, "w") as fh:
        fh.write(twin)
    stw = run(exe, sealed)
    os.remove(sealed)
    mx_spread = 0.0
    mx_react = 0.0
    for fs in stw:
        vals = [fs["cells"][i]["wax"] for i in range(nsec)]
        mx_spread = max(mx_spread, max(vals) - min(vals))
    for n in range(1, len(stw)):
        prev_w = stw[n - 1]["cells"][0]["wax"]
        pred = (prev_w / dt + k * weq) / (1.0 / dt + k)
        mx_react = max(mx_react, abs(stw[n]["cells"][0]["wax"] - pred))
    check("sealed twin cells are all identical (pure local reaction)",
          mx_spread < 1e-12, f"max spread {mx_spread:.3e}")
    check("sealed twin matches standalone reaction (5.96 path unchanged)",
          mx_react < 1e-9, f"max {mx_react:.3e}")

    # front causality: compare the transported case to the pure-reaction twin.
    # Upwind with cF = dt*F/dx = 0.5 carries the inlet disturbance exactly one
    # cell per step, so at frame n the affected cells are precisely 0..n-1;
    # every cell >= n must still be bit-for-bit identical to pure reaction.
    cF = dt * F / dx
    front_n = 6
    fb, fs = frames[front_n], stw[front_n]
    reached = [i for i in range(nsec)
               if abs(fb["cells"][i]["wax"] - fs["cells"][i]["wax"]) > 1e-12]
    check("front reaches exactly the upstream cells (one cell per step)",
          reached == list(range(min(front_n, nsec))),
          f"reached={reached[:12]} cF={cF}")
    check("front is causal (cells ahead of the disturbance unchanged)",
          all(abs(fb["cells"][i]["wax"] - fs["cells"][i]["wax"]) < 1e-12
              for i in range(front_n, nsec)),
          "a cell ahead of the front differs from pure reaction")

    # at that early frame the inflow fraction 0.1 is still above the transient
    # deposit, so advection has ADDED wax at the inlet relative to pure
    # reaction. (Later the sign flips: the inflow holds w below w_eq.)
    early_gain = fb["cells"][0]["wax"] - fs["cells"][0]["wax"]
    check("inlet inflow adds deposit at the front (early frame)",
          early_gain > 1e-3, f"gain {early_gain:.4f}")

    # late-time inlet cell settles to the closed-form advection/reaction
    # steady state (cF*win + k*dt*weq)/(cF + k*dt). Since win < weq this sits
    # strictly below weq, and is far from the sealed twin's weq: a direct
    # witness that the wxf[0] boundary term is active.
    w_ss = (cF * win + k * dt * weq) / (cF + k * dt)
    w_inlet = frames[-1]["cells"][0]["wax"]
    check("inlet cell approaches advection-reaction steady state",
          abs(w_inlet - w_ss) < 1e-3, f"got {w_inlet:.6f} want {w_ss:.6f}")
    check("inlet steady state sits between win and weq",
          win < w_inlet < weq, f"win={win} w={w_inlet:.6f} weq={weq:.6f}")

    if FAILS:
        print("\n".join("FAIL: " + m for m in FAILS))
        raise SystemExit(f"{len(FAILS)} wax-transport check(s) failed")
    print(f"wax-transport: {len(frames)} frames x {nsec} cells, F={F} dx={dx}, "
          f"weq={weq:.6f}, worst recursion {worst:.3e} - all checks passed")


def tempfile_path():
    import tempfile
    fd, p = tempfile.mkstemp(suffix=".dict")
    os.close(fd)
    return fd, p


if __name__ == "__main__":
    main()
