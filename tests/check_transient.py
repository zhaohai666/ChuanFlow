#!/usr/bin/env python3
"""Analytic verification of the 1D transient kernel (chflow-1d-transient).

usage: check_transient.py <exe> <case.dict> <mode>
  mode = blowdown : tank draining through a laminar line. The implicit-Euler
            (Q, p) system is 2x2 linear, so we replicate the exact discrete
            update (tight tolerance) AND the physical slow decay
            exp(-t / (R*C)) (2 %).
  mode = rlc      : inertial sloshing between two compliant nodes, compared
            against the exact underdamped RLC solution for Q(t).
  mode = slug     : liquid fraction advected through incompressible ports at
            constant flow: exact cascade of volume delays.
  mode = valve    : tank blowdown through a line whose valve closes and
            reopens on a schedule. Substituting the node compliance into the
            pipe momentum row collapses each step to one monotone scalar
            equation (with K(t) evaluated at the implicit new time), solved
            here by bisection: exact (Q, p, alpha) replication plus anchors
            for the closed-valve leak, the frozen tank pressure and the
            resumed discharge after reopen.
"""
import math
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

exe, case, mode = sys.argv[1], sys.argv[2], sys.argv[3]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")

# ---- parse the time series ----
times, nodes_series, pipes_series = [], [], []
cur_t, cur_nodes, cur_pipes = None, {}, {}


def flush():
    global cur_t
    if cur_t is not None:
        times.append(cur_t)
        nodes_series.append(cur_nodes)
        pipes_series.append(cur_pipes)
        cur_t = None


for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "t":
        flush()
        cur_t, cur_nodes, cur_pipes = float(tok[1]), {}, {}
    elif tok[0] == "node":
        cur_nodes[tok[1]] = {
            "p": float(tok[3]), "alpha": float(tok[5]), "h": float(tok[7])}
    elif tok[0] == "pipe":
        cur_pipes[tok[1]] = float(tok[3])
flush()
if not times:
    raise SystemExit("no time series parsed")

fails = []


def check(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got={got:.6e} want={want:.6e}")
    if not ok:
        fails.append(name)


MU, RHO = 0.1, 850.0


def resistance(length, d):
    return 128.0 * MU * length / (math.pi * d ** 4)


if mode == "blowdown":
    V, K, L, DI = 1.0, 2.0e8, 1000.0, 0.02
    Rr = resistance(L, DI)
    C = V / K
    inertia = RHO * L / (0.25 * math.pi * DI * DI)
    p, q = 5.0e5, 0.0
    dt = times[1] - times[0]
    # exact replication of the simultaneous implicit-Euler 2x2 system
    # for unknowns (q, p) at the new time level:
    #   (I/dt + R) q - p           = (I/dt) q_prev   [pipe momentum]
    #   (dt K / V) q + p           = p_prev          [node compliance]
    ref = []
    t = 0.0
    for _ in times:
        ref.append((t, p, q))
        qn = (p + (inertia / dt) * q) / (inertia / dt + Rr + dt * K / V)
        pn = p - dt * qn * K / V
        q, p, t = qn, pn, t + dt
    for i, t in enumerate(times):
        check(f"p(t={t:.2f})", nodes_series[i]["TANK"]["p"], ref[i][1], 1e-9)
        check(f"Q(t={t:.2f})", pipes_series[i]["LINE"], ref[i][2], 1e-9)
    # physics: slow eigenvalue ~ exp(-t/(R C)); last sample
    tau = Rr * C
    check("continuous envelope p(end)",
          nodes_series[-1]["TANK"]["p"], 5.0e5 * math.exp(-times[-1] / tau), 0.02)

elif mode == "rlc":
    V, K = 1.0, 1.0e8
    L, DI = 2.49, 0.0519
    Rr = resistance(L, DI)
    inertia = RHO * L / (0.25 * math.pi * DI * DI)
    ceq = (V / K) / 2.0
    w0 = 1.0 / math.sqrt(inertia * ceq)
    zeta = 0.5 * Rr * math.sqrt(ceq / inertia)
    dp0 = 1.0e5
    if zeta >= 1.0:
        raise SystemExit("case misconfigured: not underdamped")
    wd = w0 * math.sqrt(1.0 - zeta * zeta)
    t_peak = math.atan(wd / (zeta * w0)) / wd
    q_peak = dp0 / (inertia * wd) * math.exp(-zeta * w0 * t_peak) * math.sin(wd * t_peak)

    qs = [(t, s["CROSS"]) for t, s in zip(times, pipes_series)]
    t_max, q_max = max(qs, key=lambda x: x[1])
    check("first peak time", t_max, t_peak, 0.03, 2 * (times[1] - times[0]))
    check("first peak amplitude", q_max, q_peak, 0.03)
    # oscillatory: sign must reverse at least twice within 1.5 s
    signs = [1 if q > 1e-9 else (-1 if q < -1e-9 else 0) for _, q in qs]
    flips = sum(1 for a, b in zip(signs, signs[1:]) if a > 0 and b < 0)
    flips += sum(1 for a, b in zip(signs, signs[1:]) if a < 0 and b > 0)
    if flips < 4:
        fails.append(f"flow did not oscillate (flips={flips})")
    else:
        print(f"ok:   oscillation sign flips = {flips} (>= 4)")
    # nodes exchange liquid but total inventory is conserved (closed system)
    v_tot = sum(s["A"]["alpha"] + s["B"]["alpha"] for s in nodes_series) / len(nodes_series)
    check("closed-system liquid inventory", v_tot, 1.0, 1e-6)

elif mode == "slug":
    L, DI = 100.0, 0.02
    Rr = resistance(L, DI)
    q_want = 6.0e5 / (5.0 * Rr)
    v_node = 1.0
    dt = times[1] - times[0]
    # The nodes are well-mixed volumes; the kernel transports with a snapshot
    # of the old fractions, so the exact discrete recursion is a closed-form
    # reference: replicate it from the reported flows and compare to 1e-9.
    a2 = a3 = a4 = 0.0
    for i in range(1, len(times)):
        q = pipes_series[i]["PC"]
        check(f"Q(t={times[i]:.0f})", q, q_want, 1e-3)
        d = q * dt / v_node
        new2 = a2 + d * (1.0 - a2)
        new3 = a3 + d * (a2 - a3)
        new4 = a4 + d * (a3 - a4)
        a2, a3, a4 = new2, new3, new4
        s = nodes_series[i]
        check(f"alpha2(t={times[i]:.0f})", s["N2"]["alpha"], a2, 1e-9)
        check(f"alpha3(t={times[i]:.0f})", s["N3"]["alpha"], a3, 1e-9)
        check(f"alpha4(t={times[i]:.0f})", s["N4"]["alpha"], a4, 1e-9)
    # slug physically reaches the terminal long before endTime
    if a4 < 0.5:
        fails.append("slug never reached the terminal")
    else:
        print(f"ok:   slug reached terminal (alpha4={a4:.3f})")

elif mode == "valve":
    V, K, L, DI = 5.0, 4.0e7, 100.0, 0.05
    KV, KC = 20.0, 1.0e12
    CS, CE, OS, OE = 2.0, 6.0, 10.0, 14.0
    Rr = resistance(L, DI)
    area = math.pi * DI * DI / 4.0
    inertia = RHO * L / area
    ck = K / V
    dt = times[1] - times[0]

    def k_valve(t):
        if t < CS:
            return KV
        if t < CE:
            return KV + (t - CS) / (CE - CS) * (KC - KV)
        if t < OS:
            return KC
        if t < OE:
            return KC + (t - OS) / (OE - OS) * (KV - KC)
        return KV

    def valve_dp(qq, kv):
        return kv * RHO * abs(qq) * qq / (2.0 * area * area)

    # per step: p_new = p_prev - (K dt / V) q substituted into
    # inertia/dt (q - q_prev) + R q + valve_dp(q, K(t+dt)) = p_new.
    # The schedule time must be accumulated in lockstep with the kernel
    # (t += dt each step, same doubles): the printed sample time carries
    # only 12 significant digits, and on a ramp that dt error maps to
    # dK = slope * dt ~ 1e2 - enough to break the 1e-9 gate near K = 20.
    p, q, a, tc = 5.0e5, 0.0, 1.0, 0.0
    for i in range(1, len(times)):
        tc += dt
        kn = k_valve(tc)

        def resid(qq, kn=kn, q_old=q, p_old=p):
            return (inertia / dt * (qq - q_old) + Rr * qq
                    + valve_dp(qq, kn) + ck * dt * qq - p_old)

        lo, hi = 0.0, p / (Rr + ck * dt)
        while resid(hi) < 0:
            hi *= 2.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if resid(mid) > 0:
                hi = mid
            else:
                lo = mid
        q = 0.5 * (lo + hi)
        check(f"Q(t={times[i]:.2f})", pipes_series[i]["LINE"], q, 1e-9)
        p = p - ck * dt * q
        check(f"p(t={times[i]:.2f})", nodes_series[i]["TANK"]["p"], p, 1e-9)
        a = a * (1.0 - dt * q / V)
        check(f"alpha(t={times[i]:.2f})", nodes_series[i]["TANK"]["alpha"], a, 1e-9)

    # closed window [8, 10]: K = 1e12 strangles the line to a leak
    closed = [abs(pipes_series[i]["LINE"]) for i in range(round(8.0 / dt),
                                                          round(10.0 / dt) + 1)]
    check("closed-valve leak", max(closed), 0.0, 0.0, 1e-6)
    i8, i10 = round(8.0 / dt), round(10.0 / dt)
    p8, p10 = nodes_series[i8]["TANK"]["p"], nodes_series[i10]["TANK"]["p"]
    check("tank pressure frozen while closed", abs(p10 - p8), 0.0, 0.0, 10.0)
    # reopen resumes the discharge at the quasi-steady p/R level
    q_end, p_end = pipes_series[-1]["LINE"], nodes_series[-1]["TANK"]["p"]
    if q_end <= 1.0e-3:
        fails.append(f"reopen did not resume the discharge (q(end)={q_end:.3e})")
    else:
        print(f"ok:   reopen resumes discharge (q(end)={q_end:.3e} m3/s)")
    if not p_end < p10 - 1.0e4:
        fails.append("tank pressure did not resume dropping after reopen")
    else:
        print(f"ok:   pressure dropping again (p dropped {p10 - p_end:.0f} Pa after reopen)")
    # laminar throughout: Churchill reduces to 64/Re below ~2000, which the
    # exact replication above relies on
    re_max = max(RHO * abs(s["LINE"]) / area * DI / MU for s in pipes_series)
    if re_max >= 2000.0:
        fails.append(f"case left the laminar regime (max Re = {re_max:.0f})")
    else:
        print(f"ok:   laminar throughout (max Re = {re_max:.0f})")

else:
    raise SystemExit(f"unknown mode {mode}")

if fails:
    raise SystemExit(f"transient {mode} verification FAILED: {fails}")
print(f"transient {mode} case matches the analytic reference")
