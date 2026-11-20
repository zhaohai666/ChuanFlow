#!/usr/bin/env python3
"""Closed-form verification of the three G9 point elements (5.140):

  * jet pump (constant-area mixing momentum theory, Thornhcser 1955):
    p_d - p_s = rho vn^2 (m-1-r)^2/(m^2(m-1)) - f (L/D) rho v3^2/2, with
    continuity q_d = q_n + q_s and v3 = v_n (1+r)/m. Every arm's printed
    rates must satisfy that momentum balance to 1e-9; the ideal arm also
    pins the perfect-square root r = (m-1) - m sqrt(dp/(rho vn^2)) pointwise,
    and the imposed-rate arm round-trips through the same equation.
  * progressive cavity pump: q = cavitiesPerRev * cavityVolume * rpm/60
    - slip * dp recomputed exactly, power = q * dp, ev = q/q_th.
  * plunger gas lift: the released-gas schedule integrates to the declared
    batch (q_w * w == V_g, q_avg * T == V_g), the cycle window is exactly
    (t mod T) < w, the thresholds are the 5.89 bellows balance at depth and
    the nozzle cap is the 5.89 isentropic flux, recomputed here.

Plus the loud rejections: negative areas, rpm 0 with a declared machine,
negative slip, a window longer than its cycle, a back-pressured jet.

usage: check_g9_point_elements.py <exe> <case.dict>
"""
import math
import os
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

G = 9.80665


def run(exe_path, case_path):
    return subprocess.run([exe_path, case_path], capture_output=True, text=True)


def parse(stdout):
    jets, pcps, pgls, windows = {}, {}, {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "jetPump":
            d = dict(zip(tok[4::2], map(float, tok[5::2])))
            d["state"] = tok[3]
            jets[tok[1]] = d
        elif tok[0] == "pcp":
            d = dict(zip(tok[4::2], map(float, tok[5::2])))
            d["state"] = tok[3]
            pcps[tok[1]] = d
        elif tok[0] == "pgl":
            d = dict(zip(tok[4::2], map(float, tok[5::2])))
            d["state"] = tok[3]
            pgls[tok[1]] = d
        elif tok[0] == "pglWindow":
            windows.setdefault(tok[1], []).append(
                (float(tok[3]), float(tok[5])))
    return jets, pcps, pgls, windows


proc = run(exe, case)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
jets, pcps, pgls, windows = parse(proc.stdout)

fails = []
ncrit = 0


def check(ok, msg):
    global ncrit
    ncrit += 1
    if not ok:
        fails.append(msg)


def close(a, b, tol=1e-9, floor=1e-12):
    return abs(a - b) <= tol * max(abs(a), abs(b), floor)


check("point elements: 5 jet pumps, 5 pcp pumps, 4 plunger gas lifts"
      in proc.stdout, "element count line wrong")

# ---------------------------------------------------------------- jet pumps --
# geometry shared by the arms (mirrors the case dict; the checker re-derives
# every answer from it, trusting only the printed numbers as the unknowns)
RHO, AN, AT, QN = 850.0, 1.0e-4, 4.0e-4, 2.0e-3
M = AT / AN
VN = QN / AN
DT = math.sqrt(4.0 * AT / math.pi)


def gmatch(r, fline):
    s = M - 1.0 - r
    return s * s / (M - 1.0) - 0.5 * fline * (1.0 + r) ** 2


def solve_r(dp, fline):
    """root of g(r) = m^2 dp/(rho vn^2) on [0, m-1], or None when stalled."""
    target = dp * M * M / (RHO * VN * VN)
    if target >= gmatch(0.0, fline):
        return None
    a = 1.0 / (M - 1.0) - 0.5 * fline
    b = -(2.0 + fline)
    c = (M - 1.0) - 0.5 * fline - target
    if abs(a) < 1e-12:
        return -c / b
    disc = b * b - 4.0 * a * c
    return (-b - math.sqrt(disc)) / (2.0 * a)


for name, fline, ps, pd, state in [
        ("JP_IDEAL", 0.0, 2.0e6, 2.05e6, "PUMPING"),
        ("JP_FRICTION", 0.02 * 0.2 / DT, 2.0e6, 2.05e6, "PUMPING"),
        ("JP_STALL", 0.0, 2.0e6, 2.07e6, "NO-SUCTION"),
        ("JP_FREEJET", 0.0, 2.0e6, 2.0e6, "FREE-JET")]:
    p = jets.get(name)
    if p is None:
        check(False, f"{name}: no printed row")
        continue
    check(p["state"] == state, f"{name}: state {p['state']} != {state}")
    dp = pd - ps
    r = solve_r(dp, fline)
    r_exp = 0.0 if r is None else min(max(r, 0.0), M - 1.0)
    check(close(p["r"], r_exp), f"{name}: r {p['r']} != closed form {r_exp}")
    qs = r_exp * QN
    check(close(p["qs"], qs), f"{name}: qs {p['qs']} != {qs}")
    check(close(p["qd"], QN + qs), f"{name}: qd != qn + qs (continuity)")
    # the momentum balance itself, recomputed from the PRINTED rates; a
    # clamped NO-SUCTION state does not satisfy it (that is the stall) and
    # instead must print the exact overshoot as its residual
    r_pr = p["qs"] / p["qn"]
    v3_pr = p["v3"]
    lhs = dp
    rhs = (RHO * VN * VN * (M - 1.0 - r_pr) ** 2 / (M * M * (M - 1.0))
           - fline * RHO * v3_pr * v3_pr / 2.0)
    if state == "NO-SUCTION":
        check(abs(p["resid"] - (rhs - lhs)) <= 1e-9 * max(abs(lhs), 1.0),
              f"{name}: stall residual {p['resid']} != shutoff overshoot "
              f"{rhs - lhs}")
    else:
        check(abs(lhs - rhs) <= 1e-9 * max(abs(lhs), 1.0),
              f"{name}: momentum balance off by {lhs - rhs} Pa")
    check(close(p["psupply"], ps + 0.5 * RHO * VN * VN),
          f"{name}: supply pressure != ps + rho vn^2/2")
    if state == "PUMPING":
        check(abs(p["resid"]) <= 1e-9 * dp, f"{name}: residual {p['resid']}")
        N = dp / (p["psupply"] - pd)
        check(close(p["N"], N), f"{name}: head ratio {p['N']} != {N}")
        check(close(p["eff"], p["r"] * N), f"{name}: eff != r N")
        check(close(p["pu"], qs * dp), f"{name}: useful power != qs dp")
        check(close(p["pm"], QN * (p["psupply"] - pd)),
              f"{name}: motive power != qn (psupply - pd)")
        check(close(p["eff"], p["pu"] / p["pm"]), f"{name}: eff != pu/pm")
    if state == "NO-SUCTION":
        check(p["qs"] == 0.0 and p["qd"] == QN,
              "NO-SUCTION must pass the motive stream alone")
        check(p["resid"] < 0.0, "stall must report the overshoot as a number")
    if state == "FREE-JET":
        check(close(p["r"], M - 1.0), "free jet must sit at r = m-1")

# imposed-rate mode: r declared, chamber pressure is the answer; round-trip
pi = jets.get("JP_IMPOSED")
if pi is None:
    check(False, "JP_IMPOSED: no printed row")
else:
    fline = 0.02 * 0.2 / DT
    r = 1.0e-3 / QN
    dp = RHO * VN * VN * gmatch(r, fline) / (M * M)
    ps_exp = 2.05e6 - dp
    check(pi["state"] == "PUMPING", f"JP_IMPOSED state {pi['state']}")
    check(close(pi["ps"], ps_exp), f"JP_IMPOSED ps {pi['ps']} != {ps_exp}")
    # round trip: the solved chamber pressure fed back through mode A gives
    # the same suction rate
    r_rt = solve_r(2.05e6 - ps_exp, fline)
    check(r_rt is not None and close(r_rt, r),
          f"JP_IMPOSED round trip {r_rt} != {r}")
    check(close(pi["qs"], 1.0e-3), "JP_IMPOSED must honor the declared rate")

# --------------------------------------------------------------------- pcp --
for name, cav, vol, rpm, slip, dp, state in [
        ("PCP_BASE", 2, 1.2e-5, 120, 2.0e-12, 8.0e6, "PUMPING"),
        ("PCP_NOSLIP", 2, 1.2e-5, 120, 0.0, 8.0e6, "PUMPING"),
        ("PCP_STALL", 2, 1.2e-5, 120, 1.0e-10, 8.0e6, "STALLED"),
        ("PCP_NOLOAD", 2, 1.2e-5, 120, 2.0e-12, 0.0, "NO-LOAD"),
        ("PCP_PORT", 2, 1.2e-5, 120, 2.0e-12, 8.0e6, "PUMPING")]:
    p = pcps.get(name)
    if p is None:
        check(False, f"{name}: no printed row")
        continue
    qth = cav * vol * rpm / 60.0
    q = qth - slip * dp
    check(p["state"] == state, f"{name}: state {p['state']} != {state}")
    check(close(p["qth"], qth), f"{name}: qth {p['qth']} != {qth}")
    check(close(p["slip"], slip * dp), f"{name}: slip rate wrong")
    check(close(p["q"], q), f"{name}: q {p['q']} != qth - slip dp {q}")
    check(close(p["m"], 900.0 * q), f"{name}: mass rate != rho q")
    check(close(p["dp"], dp), f"{name}: clamped rise {p['dp']} != {dp}")
    check(close(p["ev"], q / qth), f"{name}: ev != q/qth")
    check(close(p["power"], q * dp), f"{name}: power != q dp")
check(close(pcps["PCP_BASE"]["q"], pcps["PCP_PORT"]["q"]),
      "PCP_PORT must reproduce PCP_BASE from its two ports")

# ------------------------------------------------------- plunger gas lift --
V_G, PREF, RHO_REF, T, W, Z, RHOL = 0.5, 1.0e5, 1.2, 600.0, 30.0, 1500.0, 900.0
PC_S, PT_S, R, KD, CD, AP = 8.0e6, 4.0e6, 0.10, 1.3, 0.665, \
    0.25 * math.pi * 0.01 ** 2
PC = PC_S + RHOL * G * Z
PT = PT_S + RHOL * G * Z


def pvo(pdm, pt):
    return pdm / (1.0 - R) - R / (1.0 - R) * pt


def nozzle(p_up, rho_up, p_down):
    rc = (2.0 / (KD + 1.0)) ** (KD / (KD - 1.0))
    rr = p_down / p_up
    if rr >= 1.0:
        return 0.0, False
    rs = rr if rr > rc else rc  # choked flow saturates the bracket at rc
    br = rs ** (2.0 / KD) - rs ** ((KD + 1.0) / KD)
    return CD * AP * math.sqrt(2.0 * rho_up * p_up * KD / (KD - 1.0) * br), \
        rr <= rc


for name, pdm, w, state in [("PGL_OPEN", 9.0e6, 30.0, "OPEN"),
                            ("PGL_LIMITED", 9.0e6, 0.001, "VALVE-LIMITED"),
                            ("PGL_INSPREAD", 21.0e6, 30.0, "IN-SPREAD"),
                            ("PGL_CLOSED", 22.0e6, 30.0, "CLOSED")]:
    p = pgls.get(name)
    if p is None:
        check(False, f"{name}: no printed row")
        continue
    check(p["state"] == state, f"{name}: state {p['state']} != {state}")
    check(close(p["pc"], PC), f"{name}: Pc(z) != casing + rho g z")
    check(close(p["pt"], PT), f"{name}: Pt(z) != tubing + rho g z")
    check(close(p["pvo"], pvo(pdm, PT)), f"{name}: Pvo != bellows balance")
    check(close(p["pvc"], pdm), f"{name}: Pvc != dome (St = 0)")
    mreq = RHO_REF * V_G / w
    check(close(p["mreq"], mreq), f"{name}: required rate != rho V / w")
    if state in ("OPEN", "VALVE-LIMITED"):
        mcap, choked = nozzle(PC, RHO_REF * PC / PREF, PT)
        check(close(p["mcap"], mcap), f"{name}: nozzle cap != isentropic flux")
        check(p["choked"] == (1.0 if choked else 0.0), f"{name}: choked flag")
        mdel = min(mcap, mreq)
        vol = mdel * w / RHO_REF
        check(close(p["mw"], mdel), f"{name}: window mass rate")
        check(close(p["qw"] * w, vol), f"{name}: schedule integral != batch")
        check(close(p["volref"], vol), f"{name}: delivered volume")
        check(close(p["qavg"] * T, vol), f"{name}: average x period != batch")
        check(close(p["voldh"], vol * PREF / PC), f"{name}: Boyle volume")
        check(close(p["qwd"], p["qw"] * PREF / PC), f"{name}: Boyle rate")
        if state == "OPEN":
            check(close(p["volref"], V_G),
                  "the OPEN batch must deliver the declared volume exactly")
        else:
            check(close(p["volref"], mcap * w / RHO_REF) and p["volref"] < V_G,
                  "VALVE-LIMITED must report the shortfall, not the promise")
    else:
        check(p["volref"] == 0.0 and p["qw"] == 0.0 and p["mcap"] == 0.0,
              f"{state} delivers nothing")

# the cycle window: open iff (t mod T) < w, edge probes included
for t, o in windows.get("PGL_OPEN", []):
    check(o == (1.0 if (t % T) < W else 0.0),
          f"window probe t={t}: {o} != expected")
check(len(windows.get("PGL_OPEN", [])) == 7, "window probes missing")

# ------------------------------------------------------------- rejections --
with open(case) as f:
    base_text = f.read()


def probe(mutate, needle, label):
    with tempfile.TemporaryDirectory() as td:
        text = mutate(base_text)
        cp = os.path.join(td, "case.dict")
        with open(cp, "w") as f:
            f.write(text)
        pr = run(exe, cp)
        check(pr.returncode != 0 and needle in (pr.stdout + pr.stderr),
              f"rejection '{label}' not loud (needle '{needle}')")


probe(lambda t: t.replace("nozzleArea 1.0e-4; throatArea 4.0e-4;\n                  throatLength 0; frictionFactor 0;\n                  motiveRate 2.0e-3; suctionPressure 2.0e6;\n                  dischargePressure 2.05e6",
                         "nozzleArea -1.0e-4; throatArea 4.0e-4; motiveRate 2.0e-3; suctionPressure 2.0e6; dischargePressure 2.05e6"),
      "nozzle area is positive", "negative nozzle area")
probe(lambda t: t.replace("nozzleArea 1.0e-4; throatArea 4.0e-4;\n                  throatLength 0; frictionFactor 0;\n                  motiveRate 2.0e-3; suctionPressure 2.0e6;\n                  dischargePressure 2.07e6",
                         "nozzleArea 5.0e-4; throatArea 4.0e-4; motiveRate 2.0e-3; suctionPressure 2.0e6; dischargePressure 2.07e6"),
      "throat area must exceed", "throat smaller than nozzle")
probe(lambda t: t.replace("dischargePressure 2.05e6; }\n    JP_FRICTION",
                         "dischargePressure 1.9e6; }\n    JP_FRICTION"),
      "below the suction", "discharge below suction")
probe(lambda t: t.replace("suctionRate 1.0e-3", "suctionRate 7.0e-3"),
      "outside [0, (m-1)", "imposed suction past the free-jet ceiling")
probe(lambda t: t.replace("rpm 120;\n                 slip 2.0e-12; pressureRise 8.0e6",
                         "rpm 0;\n                 slip 2.0e-12; pressureRise 8.0e6", 1),
      "rpm is positive", "rpm 0 with a declared machine")
probe(lambda t: t.replace("slip 2.0e-12; pressureRise 8.0e6",
                         "slip -2.0e-12; pressureRise 8.0e6", 1),
      "slip line is non-negative", "negative slip")
probe(lambda t: t.replace("cavityVolume 1.2e-5", "cavityVolume -1.2e-5", 1),
      "cavity volume is positive", "negative cavity volume")
probe(lambda t: t.replace("openWindow 30;\n                   releaseDepth 1500; casingSurfacePressure 8.0e6;\n                   tubingSurfacePressure 4.0e6; liquidDensity 900;\n                   liftDomePressure 9.0e6",
                         "openWindow 700;\n                   releaseDepth 1500; casingSurfacePressure 8.0e6;\n                   tubingSurfacePressure 4.0e6; liquidDensity 900;\n                   liftDomePressure 9.0e6"),
      "open window must lie", "window longer than the cycle")
probe(lambda t: t.replace("gasVolumePerCycle 0.5; referencePressure 1.0e5;\n                   gasDensityAtReference 1.2; cyclePeriod 600; openWindow 0.001",
                         "gasVolumePerCycle -0.5; referencePressure 1.0e5; gasDensityAtReference 1.2; cyclePeriod 600; openWindow 0.001"),
      "released gas volume per cycle is positive", "negative batch volume")
probe(lambda t: t.replace("releaseDepth 1500; casingSurfacePressure 8.0e6;\n                   tubingSurfacePressure 4.0e6; liquidDensity 900;\n                   liftDomePressure 22.0e6",
                         "releaseDepth -10; casingSurfacePressure 8.0e6; tubingSurfacePressure 4.0e6; liquidDensity 900; liftDomePressure 22.0e6"),
      "non-negative vertical", "negative release depth")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"g9-point-elements: {len(fails)} of {ncrit} checks FAILED")
print(f"g9-point-elements: all {ncrit} checks passed")
