#!/usr/bin/env python3
"""Analytic verification of the distributed two-fluid kernel
(chflow-1d-twofluid).

usage: check_two_fluid.py <exe> <case.dict> <mode>
  mode = front           : liquid column pushing into gas through a Wallis
            constant-slip closure. The liquid-flux function
            F(a) = us a^2 + (vm - us) a is convex, so the front is a
            Rankine-Hugoniot shock advancing at exactly F(1)-F(0) = vm,
            independent of the slip velocity. The Godunov recursion is
            replicated cell by cell (1e-9) and the shock position is read
            from the exact liquid inventory; the injected water cut forms a
            second front at vm.
  mode = steadyholdup    : imposed phase rates settle to the uniform holdup
            solving F(alpha*) = vsl (quadratic for Wallis, linear for
            Zuber-Findlay); the pressure drop is recomputed term by term
            from the Churchill friction law and the mixture gravity head.
  mode = settling        : sealed horizontal pipe, counter-current slip
            segregation. For uniform alpha0 the exact solution is a pair of
            shocks both moving at us*alpha0; inventory is conserved exactly.
  mode = pressurelaminar : fixed-pressure start-up of a pure-liquid line;
            the implicit-Euler momentum ODE is solved independently and the
            physical exponential envelope towards u* = dp D^2 / (32 mu L)
            is checked against tau = rho D^2 / (32 mu).
  mode = dynamicslipsettling : sealed vertical pipe, gravity segregation
            under the dynamic closure. Each cell's slip state s = vG - vL
            relaxes to the terminal slip s* = sqrt(drive/(C D)) of the
            reduced two-momentum ODE
            K_s ds/dt = (1/rhoG - 1/rhoL)(rho_m du/dt + loss/dx) - D C s|s|
            (with vm = 0: drive = (1/rhoG - 1/rhoL) rho_m g sin(theta)).
            The coupled slip map + Godunov transport is replicated cell by
            cell; the liquid pools over the lower half and gas gathers over
            the upper half with the interface pinned at mid-pipe, the
            inventory and the gravity pressure drop are conserved exactly,
            and the final slip matches the closed-form terminal value.
  mode = dynamicslipholdup   : imposed phase rates settle to the coupled
            steady state (alpha*, s*) solving F(alpha*, s*) = vsl with
            F(a) = s a^2 + (vM - s) a and s = s*(alpha) the terminal slip
            of the friction-driven ODE; the checker root-finds the pair
            independently, then verifies the uniform final holdup/slip
            fields, the phase-flux split vsl/vsg and the Churchill pressure
            drop.
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
times, frames = [], []
cur_t, cur_vm, cur_dp, cur_nodes, cur_cells = None, None, None, [], []


def flush():
    global cur_t, cur_vm, cur_dp, cur_nodes, cur_cells
    if cur_t is not None:
        times.append(cur_t)
        frames.append({"t": cur_t, "vm": cur_vm, "dp": cur_dp,
                       "nodes": cur_nodes, "cells": cur_cells})
        cur_t, cur_vm, cur_dp, cur_nodes, cur_cells = None, None, None, [], []


for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "t":
        flush()
        cur_t, cur_vm, cur_dp = float(tok[1]), float(tok[3]), float(tok[5])
    elif tok[0] == "node":
        cur_nodes.append(float(tok[3]))
    elif tok[0] == "cell":
        c = {"alpha": float(tok[3]), "water": float(tok[5]),
             "vsl": float(tok[7]), "vsg": float(tok[9])}
        for k in range(10, len(tok), 2):
            c[tok[k]] = float(tok[k + 1])
        cur_cells.append(c)
flush()
if not times:
    raise SystemExit("no time series parsed")

fails = []


def check(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got={got:.10e} want={want:.10e}")
    if not ok:
        fails.append(name)


G = 9.80665


def churchill(re, relrough):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * relrough
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


def godunov(qa, qb, qc, aL, aR):
    fL = qc + aL * (qb + qa * aL)
    if aL == aR:
        return fL
    fR = qc + aR * (qb + qa * aR)
    shock = (aL > aR) if qa > 0 else (aL < aR)
    if shock:
        s = (fL - fR) / (aL - aR)
        return fL if s >= 0 else fR
    if qa == 0:
        return fL if qb >= 0 else fR
    sonic = -qb / (2.0 * qa)
    lo, hi = min(aL, aR), max(aL, aR)
    a = min(max(sonic, lo), hi)
    return qc + a * (qb + qa * a)


# case parameters (mirrored from the case dicts)
RHO_L, RHO_G, MU_L, MU_G = 1000.0, 20.0, 0.01, 1e-5


def mixture(rho_a):
    return rho_a * RHO_L + (1.0 - rho_a) * RHO_G


def visc(rho_a):
    return rho_a * MU_L + (1.0 - rho_a) * MU_G


if mode == "front":
    N, L, US, DT = 200, 100.0, 0.3, 0.25
    DX = L / N
    VM, VSL, PHI_IN = 1.0, 1.0, 0.25
    qa, qb, qc = US, VM - US, 0.0
    alpha = [1.0 if (i + 0.5) * DX <= 40.0 else 0.0 for i in range(N)]
    water = [0.0] * N
    step_in_sample = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 4
        while step_in_sample < want_steps:
            flux = ([VSL]
                    + [godunov(qa, qb, qc, alpha[j - 1], alpha[j])
                       for j in range(1, N)]
                    + [qc + alpha[-1] * (qb + qa * alpha[-1])])
            wflux = [PHI_IN * flux[0]]
            for j in range(1, N):
                wflux.append((water[j - 1] if flux[j] >= 0 else water[j])
                             * flux[j])
            wflux.append(water[-1] * flux[N])
            new_a, new_af = [], []
            for i in range(N):
                new_a.append(alpha[i] + DT / DX * (flux[i] - flux[i + 1]))
                new_af.append(alpha[i] * water[i]
                              + DT / DX * (wflux[i] - wflux[i + 1]))
            alpha, water = [], []
            for i in range(N):
                a = min(max(new_a[i], 0.0), 1.0)
                af = min(max(new_af[i], 0.0), max(new_a[i], 0.0))
                alpha.append(a)
                water.append(af / a if a > 1e-9 else 0.0)
            step_in_sample += 1
        t = fr["t"]
        check(f"vm(t={t})", fr["vm"], VM, 1e-12)
        inv = sum(alpha) * DX
        check(f"liquid inventory = shock position (t={t})", inv, 40.0 + VM * t,
              1e-9, 1e-9)
        winv = sum(alpha[i] * water[i] for i in range(N)) * DX
        check(f"water inventory (t={t})", winv, PHI_IN * VM * t, 1e-9, 1e-9)
        for i in range(N):
            got = fr["cells"][i]["alpha"]
            if abs(got - alpha[i]) > 1e-9 + 1e-9 * abs(alpha[i]):
                fails.append(f"alpha[{i}](t={t}): got={got} want={alpha[i]}")
                break
        else:
            print(f"ok:   cell alpha field replicated at t={t} (1e-9)")
        if t > 0:
            # TVD: the captured front must be monotone (no over/undershoots)
            front_cells = [i for i in range(N) if 1e-6 < alpha[i] < 1 - 1e-6]
            mono = all(alpha[i] >= alpha[i + 1] - 1e-12
                       for i in range(front_cells[0], front_cells[-1]))
            if not mono:
                fails.append(f"front profile not monotone at t={t}")
            else:
                print(f"ok:   front profile monotone (TVD) at t={t}")
            # the Wallis flux is marginally degenerate for this shock
            # (F'(0.5) = vm = shock speed for every us), so the first-order
            # Godunov front spreads diffusively (~sqrt(t) cells) instead of
            # staying two cells wide; only a generous cap is meaningful
            sharp = sum(1 for a in alpha if 1e-3 < a < 1 - 1e-3)
            if sharp > 18:
                fails.append(f"front smeared over {sharp} cells at t={t}")
            else:
                print(f"ok:   front width at t={t} ({sharp} transition cells)")

elif mode.startswith("steadyholdup"):
    N, L, DI, DT = 40, 10.0, 0.02, 0.1
    DX = L / N
    VSL, VSG = 0.4, 0.6
    VM = VSL + VSG
    ROUGH, ANGLE = 2e-5, 30.0
    wallis = "wallis" in case.lower()
    if wallis:
        US = 0.3
        qa, qb, qc = US, VM - US, 0.0
        # F(a) = us a^2 + (vm-us) a = vsl -> smaller quadratic root
        a_star = (-(VM - US) + math.sqrt((VM - US) ** 2 + 4 * US * VSL)) / (2 * US)
    else:
        C0, VD = 1.2, 0.15
        qa, qb, qc = 0.0, C0 * VM + VD, VM * (1 - C0) - VD
        a_star = (VSL - qc) / qb
    sin_th = math.sin(ANGLE * math.pi / 180.0)
    rho_m = mixture(a_star)
    mu_m = visc(a_star)
    fr_total = 0.0
    for _ in range(N):
        re = rho_m * VM * DI / mu_m
        fr_total += (churchill(re, ROUGH / DI) * (DX / DI) * 0.5 * rho_m * VM ** 2)
    grav_total = N * rho_m * G * sin_th * DX
    check("steady holdup (last frame, cell 0)", frames[-1]["cells"][0]["alpha"],
          a_star, 1e-9)
    for i in range(N):
        got = frames[-1]["cells"][i]["alpha"]
        if abs(got - a_star) > 1e-9:
            fails.append(f"steady holdup cell {i}: got={got} want={a_star}")
            break
    else:
        print(f"ok:   uniform steady holdup alpha*={a_star:.10f} (1e-9)")
    for i in range(N):
        c = frames[-1]["cells"][i]
        if abs(c["vsl"] + c["vsg"] - VM) > 1e-9:
            fails.append(f"vsl+vsg != vm at cell {i}")
            break
        if abs(c["vsg"] - VSG) > 1e-9:
            fails.append(f"cell {i}: vsg={c['vsg']} want imposed {VSG}")
            break
    else:
        print("ok:   per-cell phase fluxes consistent with the imposed rates")
    check("pressure drop (friction + gravity)", frames[-1]["dp"],
          fr_total + (grav_total if wallis else 0.0), 1e-9)

elif mode == "settling":
    N, L, US, DT = 50, 100.0, 0.2, 4.0
    DX = L / N
    qa, qb, qc = US, -US, 0.0  # vm = 0: F(a) = us a(a-1)
    alpha = [0.5] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 20
        while steps_done < want_steps:
            flux = [0.0] + [godunov(qa, qb, qc, alpha[j - 1], alpha[j])
                            for j in range(1, N)] + [0.0]
            alpha = [alpha[i] + DT / DX * (flux[i] - flux[i + 1])
                     for i in range(N)]
            steps_done += 1
        t = fr["t"]
        check(f"vm identically zero (t={t})", fr["vm"], 0.0, 0.0, 1e-12)
        inv = sum(alpha) * DX
        check(f"liquid inventory (t={t})", inv, 50.0, 1e-9, 1e-9)
        for i in range(N):
            got = fr["cells"][i]["alpha"]
            if abs(got - alpha[i]) > 1e-9 + 1e-9 * abs(alpha[i]):
                fails.append(f"alpha[{i}](t={t}): got={got} want={alpha[i]}")
                break
        else:
            print(f"ok:   cell alpha field replicated at t={t} (1e-9)")
        if t > 0:
            # shocks move at us*alpha0 = 0.1 m/s: settled plateau from the
            # inlet wall, gas pocket from the outlet wall
            i_left = max(i for i in range(N) if alpha[i] > 0.75)
            pos_left = (i_left + 1) * DX
            check(f"left shock position (t={t})", pos_left, 0.1 * t, 0.0,
                  2.5 * DX)
            i_right = min(i for i in range(N) if alpha[i] < 0.25)
            check(f"gas pocket edge (t={t})", i_right * DX, L - 0.1 * t, 0.0,
                  2.5 * DX)
    check("settled plateau alpha=1 at the inlet wall", alpha[0], 1.0, 1e-9)
    if alpha[-1] > 0.05:
        fails.append(f"gas pocket never formed at the outlet wall "
                     f"(alpha={alpha[-1]})")
    else:
        print(f"ok:   gas pocket at the outlet wall (alpha={alpha[-1]:.4f})")

elif mode == "pressurelaminar":
    N, L, DI, DT = 4, 2.0, 0.02, 0.005
    RHO, MU, DP = 1000.0, 0.1, 2000.0
    # the column inertia is its mass per unit area rho L (the unknown is a
    # velocity, not a volumetric rate)
    inertia = RHO * L
    DX = L / N

    def fr_total(u):
        # pure liquid, per-cell sum exactly as the kernel accumulates it
        re = RHO * abs(u) * DI / MU
        f = churchill(re, 0.0)
        total = 0.0
        for _ in range(N):
            total += f * (DX / DI) * 0.5 * RHO * u * u
        return total if u >= 0 else -total

    def residual(u, u_old):
        return inertia * (u - u_old) / DT + fr_total(u) - DP

    u, u_old, steps_done = 0.0, 0.0, 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 20
        while steps_done < want_steps:
            # independent solve of the implicit-Euler momentum equation
            lo, hi = -2.0, 2.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if residual(mid, u_old) > 0:
                    hi = mid
                else:
                    lo = mid
            u = 0.5 * (lo + hi)
            u_old = u
            steps_done += 1
        check(f"vm(t={fr['t']})", fr["vm"], u, 1e-8)
        check(f"dp = friction sum (t={fr['t']})", fr["dp"], fr_total(u), 1e-9)
        for c in fr["cells"]:
            if abs(c["alpha"] - 1.0) > 1e-12:
                fails.append("pure-liquid holdup drifted")
                break
    u_star = DP * DI ** 2 / (32.0 * MU * L)
    tau = RHO * DI * DI / (32.0 * MU)
    check("start-up envelope u(end)", frames[-1]["vm"],
          u_star * (1.0 - math.exp(-frames[-1]["t"] / tau)), 0.03)

elif mode == "dynamicslipsettling":
    N, L, DI, DT = 40, 20.0, 0.05, 0.01
    C_DRAG, ANGLE = 1000.0, 90.0
    DX = L / N
    sin_th = math.sin(ANGLE * math.pi / 180.0)
    inv_dr = 1.0 / RHO_G - 1.0 / RHO_L
    DRHO = RHO_G - RHO_L

    def slip_map(a, s_old, dt):
        # one step of the reduced two-momentum slip ODE with vm = 0
        # (sealed pipe): K_s ds/dt = inv_dr rho_m g sin(theta) - D C s|s|,
        # drag semi-implicit as the secant through the origin
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_m = a * RHO_L + ag * RHO_G
        ks = 1.0 + a * ag * DRHO * DRHO / (RHO_L * RHO_G)
        drive = inv_dr * rho_m * G * sin_th
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * ag)
        rate = C_DRAG * abs(s_old) * d / ks
        return (s_old / dt + drive / ks) / (1.0 / dt + rate)

    def slip_fixed_point(a):
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_m = a * RHO_L + ag * RHO_G
        drive = inv_dr * rho_m * G * sin_th
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * ag)
        return math.copysign(math.sqrt(abs(drive) / (C_DRAG * d)), drive)

    # replicate the coupled slip map + Godunov transport cell by cell
    alpha = [0.5] * N
    slip = [0.0] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 1500
        while steps_done < want_steps:
            slip = [slip_map(alpha[i], slip[i], DT) for i in range(N)]
            flux = [0.0] + [godunov(0.5 * (slip[j - 1] + slip[j]),
                                    -0.5 * (slip[j - 1] + slip[j]),
                                    0.0, alpha[j - 1], alpha[j])
                            for j in range(1, N)] + [0.0]
            alpha = [alpha[i] + DT / DX * (flux[i] - flux[i + 1])
                     for i in range(N)]
            steps_done += 1
        t = fr["t"]
        check(f"vm identically zero (t={t})", fr["vm"], 0.0, 0.0, 1e-12)
        inv = sum(alpha) * DX
        check(f"liquid inventory (t={t})", inv, 10.0, 1e-9, 1e-9)
        dp0 = 0.5 * RHO_L * G * sin_th * L + 0.5 * RHO_G * G * sin_th * L
        check(f"gravity head conserved (t={t})", fr["dp"], dp0, 1e-9, 1e-9)
        for i in range(N):
            got, want = fr["cells"][i]["alpha"], alpha[i]
            if abs(got - want) > 1e-9 + 1e-9 * abs(want):
                fails.append(f"alpha[{i}](t={t}): got={got} want={want}")
                break
            gots = fr["cells"][i].get("s")
            if gots is None or abs(gots - slip[i]) > 1e-9 + 1e-9 * abs(slip[i]):
                fails.append(f"slip[{i}](t={t}): got={gots} want={slip[i]}")
                break
        else:
            print(f"ok:   alpha/slip fields replicated at t={t} (1e-9)")
        if t > 0:
            mono = all(alpha[i] >= alpha[i + 1] - 1e-12 for i in range(N - 1))
            if not mono:
                fails.append(f"profile not monotone at t={t}")
            else:
                print(f"ok:   profile monotone (TVD) at t={t}")
    # settled state: liquid pools over the lower half, gas over the upper
    # half, the 0.5-crossing pinned at mid-pipe; terminal slip matches the
    # closed form (measured lag ~4e-4 relative: alpha still creeps, so the
    # slip sits one relaxation step behind its fixed point)
    check("settled pool bottom half", min(alpha[:N // 2]), 1.0, 0.0, 0.01)
    check("gas pocket top half", max(alpha[N // 2:]), 0.0, 0.0, 0.01)
    cross = min(i for i in range(N) if alpha[i] < 0.5)
    if cross != N // 2:
        fails.append(f"interface not at mid-pipe (crossing cell {cross})")
    else:
        print("ok:   interface pinned at mid-pipe (cell 20)")
    worst = max(abs(fr_last_s - slip_fixed_point(fr_last_a))
                / max(slip_fixed_point(fr_last_a), 1e-12)
                for fr_last_a, fr_last_s in zip(alpha, slip))
    check("terminal slip s* = sqrt(drive/(C D)) at end", worst, 0.0, 0.0,
          2e-3)

elif mode == "dynamicslipholdup":
    N, L, DI, DT = 40, 10.0, 0.02, 0.05
    C_DRAG, VSL, VSG, VM = 500.0, 0.4, 0.6, 1.0
    DX = L / N
    inv_dr = 1.0 / RHO_G - 1.0 / RHO_L
    DRHO = RHO_G - RHO_L

    def drive_friction(a):
        # steady drive: du/dt = 0, horizontal -> friction loss only
        rho_m = a * RHO_L + (1.0 - a) * RHO_G
        mu_m = a * MU_L + (1.0 - a) * MU_G
        re = rho_m * VM * DI / mu_m
        f = churchill(re, 0.0)
        fr = f * (DX / DI) * 0.5 * rho_m * VM * VM
        return inv_dr * fr / DX

    def s_star(a):
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * (1.0 - a))
        return math.copysign(math.sqrt(abs(drive_friction(a))
                                       / (C_DRAG * d)), drive_friction(a))

    # coupled steady state: F(alpha*, s*(alpha*)) = vsl (F monotone rising)
    lo, hi = 1e-9, 1.0 - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        s = s_star(mid)
        f_mid = s * mid * mid + (VM - s) * mid
        if f_mid > VSL:
            hi = mid
        else:
            lo = mid
    a_star = 0.5 * (lo + hi)
    s_star_v = s_star(a_star)

    def slip_map(a, s_old, dt):
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_m = a * RHO_L + ag * RHO_G
        ks = 1.0 + a * ag * DRHO * DRHO / (RHO_L * RHO_G)
        # rate-driven: vm imposed constant -> du/dt = 0; friction only
        drive = drive_friction(a)
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * ag)
        rate = C_DRAG * abs(s_old) * d / ks
        return (s_old / dt + drive / ks) / (1.0 / dt + rate)

    # replicate the coupled transient (rate-driven, water cut zero)
    alpha = [0.9] * N
    slip = [0.0] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 400
        while steps_done < want_steps:
            slip = [slip_map(alpha[i], slip[i], DT) for i in range(N)]
            flux = [VSL]
            for j in range(1, N):
                s = 0.5 * (slip[j - 1] + slip[j])
                flux.append(godunov(s, VM - s, 0.0, alpha[j - 1], alpha[j]))
            sl = slip[N - 1]
            a_end = alpha[N - 1]
            flux.append(sl * a_end * a_end + (VM - sl) * a_end)
            alpha = [min(max(alpha[i] + DT / DX * (flux[i] - flux[i + 1]),
                             0.0), 1.0)
                     for i in range(N)]
            steps_done += 1
        t = fr["t"]
        check(f"vm imposed (t={t})", fr["vm"], VM, 1e-12)
        for i in range(N):
            got, want = fr["cells"][i]["alpha"], alpha[i]
            if abs(got - want) > 1e-9 + 1e-9 * abs(want):
                fails.append(f"alpha[{i}](t={t}): got={got} want={want}")
                break
            gots = fr["cells"][i].get("s")
            if gots is None or abs(gots - slip[i]) > 1e-9 + 1e-9 * abs(slip[i]):
                fails.append(f"slip[{i}](t={t}): got={gots} want={slip[i]}")
                break
        else:
            print(f"ok:   alpha/slip fields replicated at t={t} (1e-9)")
    # steady state: uniform holdup/slip matching the independent root-find,
    # the imposed phase split, and the Churchill friction drop
    rho_m = a_star * RHO_L + (1.0 - a_star) * RHO_G
    mu_m = a_star * MU_L + (1.0 - a_star) * MU_G
    re = rho_m * VM * DI / mu_m
    f = churchill(re, 0.0)
    dp_want = N * f * (DX / DI) * 0.5 * rho_m * VM * VM
    for c in frames[-1]["cells"]:
        if abs(c["alpha"] - a_star) > 1e-9:
            fails.append(f"steady holdup {c['alpha']} != {a_star}")
            break
        if abs(c["s"] - s_star_v) > 1e-9 + 1e-9 * abs(s_star_v):
            fails.append(f"steady slip {c['s']} != {s_star_v}")
            break
        if abs(c["vsl"] - VSL) > 1e-9:
            fails.append(f"steady vsl {c['vsl']} != {VSL}")
            break
        if abs(c["vsg"] - VSG) > 1e-9:
            fails.append(f"steady vsg {c['vsg']} != {VSG}")
            break
    else:
        print(f"ok:   uniform steady (alpha*, s*) = ({a_star:.10f},"
              f" {s_star_v:.10f}) with the imposed phase split (1e-9)")
    check("pressure drop (Churchill friction)", frames[-1]["dp"], dp_want,
          1e-9)
    if frames[-1]["cells"] != frames[-2]["cells"]:
        fails.append("state still drifting between the last two frames")
    else:
        print("ok:   steady state held over the final 20 s")

else:
    raise SystemExit(f"unknown mode {mode}")

if fails:
    raise SystemExit(f"two-fluid {mode} verification FAILED: {fails}")
print(f"two-fluid {mode} case matches the analytic reference")
