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
  mode = compressibleholdup  : the dynamicslipholdup line with the
            ideal-gas EOS on (gasR, rhoG the injected reference density).
            The per-section gas inventory mG = rhoG (1 - alpha) dx evolves
            by upwind gas-mass continuity on the Godunov gas volume flux
            G(a) = vm + (s - vm) a - s a^2; the checker replicates the
            coupled slip + holdup + gas-mass recursion cell by cell and
            p = rhoG gasR T0 per section, the global gas inventory balance,
            and the steady state rhoG* = rhoG (injection-pinned) with the
            incompressible (alpha*, s*) pair.
  mode = wellchain        : serial pipe chain (vertical riser + horizontal
            flowline, uniform bore, dynamic closure, imposed rates). Each
            cell carries its own slope; the slip-ODE drive gains the
            per-cell gravity head
            drive_i = (1/rhoG - 1/rhoL) [fr_i/dx + rho_m g sinT_i]
            (rate-driven: du/dt = 0), and the segment junction is an
            ordinary interior Godunov face on the averaged slip. Each
            segment settles to its own per-cell root
            F(alpha*, s*(alpha*)) = vsl, s* = sqrt(drive/(C D))
            (riser alpha* ~= 0.7581 / s* ~= 1.9523; flowline
            alpha* ~= 0.5529 / s* ~= 0.6186), joined across a few-cell
            riser halo by a stationary shock carrying the flux from the
            riser root to the flowline root. The checker replicates the
            full transient cell by cell (1e-9), root-finds both segments'
            (alpha*, s*) independently for the bulk fields, recomputes the
            pressure drop term by term from the replicated profile, and
            requires flatness of the final frames.
  mode = hydratesealed   : sealed line, wallis slip 0 - every transport
            flux vanishes identically and each cell runs the pure per-cell
            hydrate recursion. Subcooling dT = T_eq(p_i) - T_i with
            T_eq = curveB / (ln p - curveA) and p_i = mg_i gasR T_i
            / ((1-a) dx); growth implicit in the free-water sink
            dh/dt = kg dT (w + h - h) -> h_new = (h/dt + kg dT (w+h))
            / (1/dt + kg dT) capped by the gas availability
            dhMax = (mg - trace)/(gasContent a dx); dissociation
            h_new = h / (1 + kd |dT| dt); water = (w+h) - h_new (the 1:1
            volume-preserving move) and the stoichiometric sink
            mg -= gasContent a dh dx; then the diagonal thermal step and
            the EOS pressure report. The checker replicates the recursion
            cell by cell (1e-9) and asserts the closed-form end state:
            full conversion h* = w0 with the gas inventory
            mg* = mg0 - gasContent a h* dx = 200 - 61 = 139 kg/m2, T* = wall
            and p* = 139 * 287 * 280 = 11,170,040 Pa exactly.
  mode = hydrateflow     : rate-driven flowing line (wallis slip 0.3, EOS,
            thermal, hydrate kinetics all on). The checker replicates the
            full coupled step cell by cell (1e-9): Godunov transport of
            alpha on F(a) = us a^2 + (vm - us) a with the water cut and the
            hydrate fraction riding the same liquid flux under the same
            clamp semantics, upwind gas-mass continuity on
            G(a) = -us a^2 + (us - vm) a + vm, the implicit hydrate
            recursion with the stoichiometric gas sink, and the implicit
            temperature solve whose enthalpy conductances are the per-face
            phase volume fluxes (rhoL cL vsl_face + rhoG cG vsg_face), each
            phase upwinded on its own flux direction. Global identities
            from emitted fields only: the hydrate-able composition
            w + h = waterFraction everywhere (the 1:1 conversion), the
            phase split vsl + vsg = vM, the gas inventory
            sum mg = m0 + net in - gasContent sum a h dx, the Churchill
            pressure drop term by term, the wall-heat sum, and the steady
            energy balance qInlet + qWall + latent release = outflow
            enthalpy.
  mode = networkspan    : two-junction pipe span - two pressure-driven
            sources (200/300 Pa) merging at J1, a span pipe `mid` with
            junctions at BOTH ends (no boundary of its own), then a split
            at J2 into two pressure-driven sinks (0 Pa). The per-step
            hydraulics are one global Newton over [five free-pipe
            velocities] + [p_J1, p_J2]; with wallis slip 0 and a uniform
            holdup each pipe is exactly the laminar implicit-Euler ODE
            I_i du_i/dt + K_i u_i = p_endIn - p_endOut (the span is driven
            by p_J1 - p_J2), closed by sum s A u = 0 per node. The checker
            solves the coupled 7x7 system independently every step (dense
            Gauss with pivoting; rtol 1e-4 covers the ~1e-6
            Churchill-vs-laminar deviation) and replicates the per-junction
            mixing + Godunov water transport + implicit upwind thermal
            solve cell by cell for all five pipes (1e-9). The steady state
            anchors on the nodal conductance closed form
            p1* = 141.9281, p2* = 12.2418 with the front-flushed phi1* and
            T1* the span hands to both sinks.
  mode = junction        : symmetric Y junction - two pressure-driven
            branches (inlets 200/300 Pa) merging at junction J into a
            common outlet line (0 Pa). With wallis slip 0 and a uniform
            holdup the mixture properties are constant, so each branch is
            exactly the laminar implicit-Euler momentum ODE
            I_i du_i/dt + K_i u_i = s_i (p_i - p_J) with I_i = rho_m L_i
            and K_i = 32 mu_m L_i / D_i^2, closed by the volume balance
            sum s_i A_i u_i = 0 and p_J = p_C + I_C (u_C - u_C_old)/dt
            + K_C u_C. The checker integrates the eliminated 2x2 system
            independently every step (rtol 1e-4 covers the ~1e-6
            Churchill-vs-laminar deviation) and replicates the junction
            mixing + Godunov water transport + implicit upwind thermal
            solve cell by cell (1e-9): the feeding branches'
            liquid-volume-flow-weighted water cut and
            heat-capacity-flow-weighted temperature set the supplied
            line's inflow state. The steady state anchors on the closed
            form p_J* = sum (A_i/K_i) p_i / sum (A_i/K_i) with the
            front-flushed phi_J* and T_J*.

  mode = entrained      : film-droplet split (entrainment first increment)
            - rate-driven horizontal pipe, wallis slip 0.3, so the gas
            velocity rides the holdup. The per-cell droplet fraction d of
            the liquid transports its content qd = alpha d at the gas
            velocity (face flux qd vG upwinded on the gas direction; the
            film flux F - qd vG completes the liquid flux F) and relaxes
            towards E_eq = clamp((vG - vCrit)/vRef, 0, 1) with implicit
            Euler at relaxTime. The injected liquid arrives 10% entrained
            against an initial d0 = 0.9, so the run exercises the droplet
            front, exchange in both directions and the steady profile. The
            imposed rate keeps vm = vsl + vsg exact, so the checker
            replicates every step cell by cell (1e-9) and anchors the last
            frame on the closed forms alpha* (root of F(alpha*) = vsl),
            the developing steady droplet profile — no uniform d* is a
            fixed point, the exact discrete one combines the face balance
            with the implicit Euler exchange into the linear recursion
            d_i = (vG* d_{i-1} + (dx/tau) E_eq*) / (vG* + dx/tau) fed by
            d_0 = (E_in vsl + alpha* (dx/tau) E_eq*) /
            (alpha* (vG* + dx/tau)) — the flushed water cut and the
            per-cell Churchill pressure drop.
  mode = fullmomentum    : full gas momentum closure (distributed pressure
            momentum first increment) - rate-driven vertical riser with the
            ideal-gas EOS on. The per-section slip solves
            ds/dt = [(rho_m - rhoG)(g sinT + aM) - sigma/alpha_g]
                    / (alpha rhoG)
            with sigma = C s|s| semi-implicit through the origin secant:
            the mixture-row pressure gradient is substituted (the gas share
            alpha_g fr cancels the proportional gas wall split), so only
            buoyancy and the mixture acceleration drive the slip and wall
            friction rides the pressure field. Rate-driven: aM = 0, so the
            terminal slip is s* = sqrt(alpha_g alpha (rhoL - rhoG) g / C)
            and the coupled fixed point F(alpha*, s*) = vsl with
            F(a) = s a^2 + (vM - s) a; volume conservation pins rhoG* = rhoG
            exactly (EOS pressure p* = rhoG gasR T0 uniform). The checker
            replicates the coupled slip + Godunov holdup + water + gas-mass
            recursion cell by cell (1e-9), verifies vsl + vsg = vm and the
            flushed water cut, the Churchill + gravity pressure drop and the
            pnode march from the outlet anchor, and the uniform steady state
            (alpha*, s*, rhoG* = rhoG).
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
if mode in ("junction", "networkspan"):
    # network emission: per frame a "t" line, then per junction a
    # "junction" line and per pipe "pipe" / "pnode" / "pcell" lines
    times, frames = [], []
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t":
            times.append(float(tok[1]))
            frames.append({"t": float(tok[1]), "junctions": {},
                           "pipes": {}, "pnode": {}, "pcell": {}})
            continue
        if not frames:
            continue
        f = frames[-1]
        if tok[0] == "junction":
            f["junctions"][tok[1]] = {"p": float(tok[3]),
                                      "net": float(tok[5])}
        elif tok[0] == "pipe":
            f["pipes"][tok[1]] = {"vm": float(tok[3]), "dp": float(tok[5])}
        elif tok[0] == "pnode":
            f["pnode"].setdefault(tok[1], []).append(float(tok[4]))
        elif tok[0] == "pcell":
            c = {"alpha": float(tok[4]), "water": float(tok[6]),
                 "vsl": float(tok[8]), "vsg": float(tok[10])}
            for k in range(11, len(tok), 2):
                c[tok[k]] = float(tok[k + 1])
            f["pcell"].setdefault(tok[1], []).append(c)
    if not times:
        raise SystemExit("no time series parsed")
else:
    times, frames = [], []
    cur_t, cur_vm, cur_dp, cur_nodes, cur_cells = None, None, None, [], []
    cur_qwall, cur_qinlet = 0.0, 0.0

    def flush():
        global cur_t, cur_vm, cur_dp, cur_nodes, cur_cells
        global cur_qwall, cur_qinlet
        if cur_t is not None:
            times.append(cur_t)
            frames.append({"t": cur_t, "vm": cur_vm, "dp": cur_dp,
                           "nodes": cur_nodes, "cells": cur_cells,
                           "qwall": cur_qwall, "qinlet": cur_qinlet})
            cur_t, cur_vm, cur_dp, cur_nodes, cur_cells = \
                None, None, None, [], []
            cur_qwall, cur_qinlet = 0.0, 0.0

    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t":
            flush()
            cur_t, cur_vm, cur_dp = float(tok[1]), float(tok[3]), float(tok[5])
        elif tok[0] == "node":
            cur_nodes.append(float(tok[3]))
        elif tok[0] == "thermal:":
            cur_qwall, cur_qinlet = float(tok[2]), float(tok[4])
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

elif mode == "compressibleholdup":
    N, L, DI, DT = 40, 10.0, 0.02, 0.05
    C_DRAG, VSL, VSG, VM = 500.0, 0.4, 0.6, 1.0
    GAS_R, T0 = 287.0, 300.0
    DX = L / N

    def rho_g_cell(mg, a):
        ag = 1.0 - a
        return mg / (ag * DX) if ag > 1e-12 else RHO_G

    def slip_map(a, s_old, mg, dt):
        # kernel slipStep with the per-section gas density of the EOS
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_g = rho_g_cell(mg, a)
        rho_m = a * RHO_L + ag * rho_g
        mu_m = a * MU_L + ag * MU_G
        re = rho_m * VM * DI / mu_m
        fr = churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM
        inv = 1.0 / rho_g - 1.0 / RHO_L
        ks = 1.0 + a * ag * (rho_g - RHO_L) ** 2 / (RHO_L * rho_g)
        drive = inv * fr / DX  # rate-driven: du/dt = 0; horizontal
        d = 1.0 / (RHO_L * a) + 1.0 / (rho_g * ag)
        rate = C_DRAG * abs(s_old) * d / ks
        return (s_old / dt + drive / ks) / (1.0 / dt + rate)

    def face_fluxes(alpha, slip):
        # liquid and gas volume fluxes of the dynamic closure: F convex,
        # G = vm - F concave, so the Godunov selections differ
        flux, gvol = [VSL], [VSG]
        for j in range(1, N):
            s = 0.5 * (slip[j - 1] + slip[j])
            flux.append(godunov(s, VM - s, 0.0, alpha[j - 1], alpha[j]))
            gvol.append(godunov(-s, s - VM, VM, alpha[j - 1], alpha[j]))
        sl = slip[N - 1]
        a_end = alpha[N - 1]
        flux.append(sl * a_end * a_end + (VM - sl) * a_end)
        gvol.append(VM + (sl - VM) * a_end - sl * a_end * a_end)
        return flux, gvol

    # independent steady state: gas mass in = out at the injected reference
    # density pins rhoG* = RHO_G, so the (alpha*, s*) pair is the
    # incompressible dynamicslipholdup one; root-find F(a, s*(a)) = vsl
    def drive_steady(a):
        rho_m = a * RHO_L + (1.0 - a) * RHO_G
        mu_m = a * MU_L + (1.0 - a) * MU_G
        re = rho_m * VM * DI / mu_m
        fr = churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM
        return (1.0 / RHO_G - 1.0 / RHO_L) * fr / DX

    def s_star_steady(a):
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * (1.0 - a))
        return math.copysign(math.sqrt(abs(drive_steady(a))
                                       / (C_DRAG * d)), drive_steady(a))

    lo, hi = 1e-9, 1.0 - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        s = s_star_steady(mid)
        if s * mid * mid + (VM - s) * mid > VSL:
            hi = mid
        else:
            lo = mid
    a_star = 0.5 * (lo + hi)
    s_star_v = s_star_steady(a_star)

    # coupled transient replication with the gas-mass state (kg/m2)
    alpha = [0.9] * N
    slip = [0.0] * N
    mg = [RHO_G * 0.1 * DX] * N
    steps_done = 0
    net_in = 0.0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 400
        while steps_done < want_steps:
            slip = [slip_map(alpha[i], slip[i], mg[i], DT) for i in range(N)]
            flux, gvol = face_fluxes(alpha, slip)
            gm = []
            for j in range(N + 1):
                if j == 0:
                    rho_up = RHO_G  # injected reference density
                elif j == N:
                    rho_up = (rho_g_cell(mg[N - 1], alpha[N - 1])
                              if gvol[j] >= 0 else 1.0e6 / (GAS_R * T0))
                elif gvol[j] >= 0:
                    rho_up = rho_g_cell(mg[j - 1], alpha[j - 1])
                else:
                    rho_up = rho_g_cell(mg[j], alpha[j])
                gm.append(rho_up * gvol[j])
            alpha = [min(max(alpha[i] + DT / DX * (flux[i] - flux[i + 1]),
                             0.0), 1.0) for i in range(N)]
            mg = [mg[i] - DT * (gm[i + 1] - gm[i]) for i in range(N)]
            net_in += DT * (gm[0] - gm[N])
            steps_done += 1
        t = fr["t"]
        check(f"vm imposed (t={t})", fr["vm"], VM, 1e-12)
        ok_frame = True
        for i in range(N):
            got = fr["cells"][i]
            want_p = rho_g_cell(mg[i], alpha[i]) * GAS_R * T0
            if (abs(got["alpha"] - alpha[i]) > 1e-9 + 1e-9 * abs(alpha[i])
                    or abs(got.get("s", 0.0) - slip[i])
                    > 1e-9 + 1e-9 * abs(slip[i])
                    or abs(got.get("p", 0.0) - want_p) > 1e-9 * want_p):
                fails.append(f"cell {i} (t={t}): got={got} want alpha="
                             f"{alpha[i]} s={slip[i]} p={want_p}")
                ok_frame = False
                break
        if ok_frame:
            print(f"ok:   alpha/slip/p fields replicated at t={t} (1e-9)")
    # global gas inventory: emitted p/alpha reconstruct the solver's
    # inventory, which must equal the replicated boundary-flux balance
    m0 = sum(RHO_G * 0.1 * DX for _ in range(N))
    m_end = sum(c["p"] * (1.0 - c["alpha"]) * DX / (GAS_R * T0)
                for c in frames[-1]["cells"])
    check("gas inventory balance", m_end, m0 + net_in, 1e-9)
    # steady state: the incompressible (alpha*, s*) pair with rhoG* = RHO_G
    rho_m = a_star * RHO_L + (1.0 - a_star) * RHO_G
    mu_m = a_star * MU_L + (1.0 - a_star) * MU_G
    re = rho_m * VM * DI / mu_m
    f = churchill(re, 0.0)
    dp_want = N * f * (DX / DI) * 0.5 * rho_m * VM * VM
    p_star = RHO_G * GAS_R * T0
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
        if abs(c["p"] - p_star) > 1e-9 * p_star:
            fails.append(f"steady p {c['p']} != {p_star} (rhoG* != rhoG)")
            break
    else:
        print(f"ok:   uniform steady (alpha*, s*) = ({a_star:.10f},"
              f" {s_star_v:.10f}) with rhoG* = rhoG and the EOS pressure"
              f" {p_star:.6e} (1e-9)")
    check("pressure drop (Churchill friction)", frames[-1]["dp"], dp_want,
          1e-9)
    if frames[-1]["cells"] != frames[-2]["cells"]:
        fails.append("state still drifting between the last two frames")
    else:
        print("ok:   steady state held over the final 20 s")

elif mode == "wellchain":
    # serial chain: vertical riser (cells 0-9, sin(theta) = 1) welded to a
    # horizontal flowline (cells 10-19, sin(theta) = 0), uniform bore,
    # dynamic closure, imposed phase rates (vM = vsl + vsg = 1 m/s). The
    # slip-ODE drive gains the per-cell gravity head
    #   drive_i = (1/rhoG - 1/rhoL) [fr_i/dx + rho_m g sinT_i]
    # (rate-driven: du/dt = 0), and the segment junction is an ordinary
    # interior Godunov face on the averaged slip, so each segment settles
    # to its own per-cell root F(alpha*, s*(alpha*)) = vsl with
    # F(a) = s a^2 + (vM - s) a, s* = sqrt(drive/(C D))
    # (riser alpha* ~= 0.7581 / s* ~= 1.9523; flowline alpha* ~= 0.5529 /
    # s* ~= 0.6186), joined across a few-cell riser halo by a stationary
    # shock: the junction face flux is exactly vsl. The checker replicates
    # the full transient cell by cell (1e-9), root-finds both segments'
    # (alpha*, s*) independently for the bulk fields, recomputes the
    # pressure drop term by term from the replicated steady profile, and
    # requires flatness of the final frames.
    N, DI, DT = 20, 0.02, 0.05
    C_DRAG, VSL, VSG, VM = 500.0, 0.4, 0.6, 1.0
    DX = 2.0
    SIN_T = [1.0] * 10 + [0.0] * 10
    STEPS_PER_FRAME = 1600
    INV_DR = 1.0 / RHO_G - 1.0 / RHO_L
    DRHO = RHO_G - RHO_L

    def drive_cell(a, sin_t):
        # loss/dx per cell of the slip ODE: Churchill friction + gravity
        rho_m = a * RHO_L + (1.0 - a) * RHO_G
        mu_m = a * MU_L + (1.0 - a) * MU_G
        re = rho_m * VM * DI / mu_m
        fr = churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM
        return INV_DR * (fr / DX + rho_m * G * sin_t)

    def slip_star(a, sin_t):
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * (1.0 - a))
        dr = drive_cell(a, sin_t)
        return math.copysign(math.sqrt(abs(dr) / (C_DRAG * d)), dr)

    def seg_root(sin_t):
        # coupled steady state: F(alpha*, s*(alpha*)) = vsl (F monotone
        # rising)
        lo, hi = 1e-9, 1.0 - 1e-9
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            s = slip_star(mid, sin_t)
            f_mid = s * mid * mid + (VM - s) * mid
            if f_mid > VSL:
                hi = mid
            else:
                lo = mid
        a_star = 0.5 * (lo + hi)
        return a_star, slip_star(a_star, sin_t)

    A_R, S_R = seg_root(1.0)
    A_F, S_F = seg_root(0.0)

    def slip_map(a, s_old, sin_t, dt):
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_m = a * RHO_L + ag * RHO_G
        ks = 1.0 + a * ag * DRHO * DRHO / (RHO_L * RHO_G)
        # rate-driven: vm imposed constant -> du/dt = 0
        drive = drive_cell(a, sin_t)
        d = 1.0 / (RHO_L * a) + 1.0 / (RHO_G * ag)
        rate = C_DRAG * abs(s_old) * d / ks
        return (s_old / dt + drive / ks) / (1.0 / dt + rate)

    # replicate the coupled transient over the concatenated chain
    alpha = [0.9] * N
    slip = [0.0] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * STEPS_PER_FRAME
        while steps_done < want_steps:
            slip = [slip_map(alpha[i], slip[i], SIN_T[i], DT)
                    for i in range(N)]
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
    # steady state: bulk cells of each segment sit at that segment's
    # independent root (riser halo cells 7-9 are the stationary-shock
    # boundary layer and are exempt from the root check)
    for i in range(7):
        c = frames[-1]["cells"][i]
        if (abs(c["alpha"] - A_R) > 1e-3) or (abs(c["s"] - S_R) > 1e-3):
            fails.append(f"riser bulk cell {i}: alpha={c['alpha']} "
                         f"s={c['s']} != root ({A_R:.10f}, {S_R:.10f})")
            break
    else:
        print(f"ok:   riser bulk cells 0-6 at the independent root "
              f"({A_R:.10f}, {S_R:.10f}) (1e-3)")
    for i in range(10, N):
        c = frames[-1]["cells"][i]
        if (abs(c["alpha"] - A_F) > 1e-9) or (abs(c["s"] - S_F) > 1e-9):
            fails.append(f"flowline cell {i}: alpha={c['alpha']} "
                         f"s={c['s']} != root ({A_F:.10f}, {S_F:.10f})")
            break
    else:
        print(f"ok:   flowline cells 10-19 at the independent root "
              f"({A_F:.10f}, {S_F:.10f}) (1e-9)")
    for i, c in enumerate(frames[-1]["cells"]):
        if abs(c["vsl"] + c["vsg"] - VM) > 1e-12:
            fails.append(f"phase-split identity cell {i}: vsl={c['vsl']} "
                         f"vsg={c['vsg']} sum={c['vsl'] + c['vsg']}")
            break
    else:
        print(f"ok:   per-cell phase split vsl + vsg = vM = {VM} (1e-12)")
    # closed-form pressure drop, term by term over the emitted profile
    dp_want = 0.0
    for i, c in enumerate(frames[-1]["cells"]):
        a = c["alpha"]
        rho_m = a * RHO_L + (1.0 - a) * RHO_G
        mu_m = a * MU_L + (1.0 - a) * MU_G
        re = rho_m * VM * DI / mu_m
        dp_want += churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM \
                 + rho_m * G * SIN_T[i] * DX
    check("pressure drop (Churchill friction + gravity, term by term)",
          frames[-1]["dp"], dp_want, 1e-9)
    if frames[-1]["cells"] != frames[-2]["cells"]:
        fails.append("state still drifting between the last two frames")
    else:
        print("ok:   steady state held over the final 80 s")

elif mode == "hydratesealed":
    # sealed line, wallis slip 0: every transport flux is identically zero,
    # so each cell runs the pure per-cell hydrate recursion (growth implicit
    # in the free-water sink, gas-availability cap, 1:1 water <-> hydrate
    # volume move, stoichiometric gas sink on the EOS inventory) followed by
    # the diagonal thermal step and the EOS pressure report
    RHO_G = 200.0  # case fluid rhoG (the shared default above is another case)
    N, DT = 10, 0.5
    DX, DI = 2.0, 0.02
    AREA = math.pi * DI * DI / 4.0
    GAS_R, T0 = 287.0, 300.0
    CL, CG = 4180.0, 2150.0
    A0, W0 = 0.5, 0.5
    KW, T_WALL = 50.0, 280.0
    CURVE_A, CURVE_B = 35.2, -5684.0
    K_G, K_D = 0.01, 0.01
    GAS_CONTENT, LATENT = 122.0, 4.2e8
    STEPS_PER_FRAME = 160

    def rho_g_cell_sealed(m, ag):
        return m / (ag * DX)

    alpha = [A0] * N
    water = [W0] * N
    hyd = [0.0] * N
    temp = [T0] * N
    mg = [RHO_G * (1.0 - A0) * DX] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * STEPS_PER_FRAME
        while steps_done < want_steps:
            for i in range(N):
                ag = 1.0 - alpha[i]
                p_i = mg[i] * GAS_R * temp[i] / (ag * DX)
                dT = CURVE_B / (math.log(p_i) - CURVE_A) - temp[i]
                total = water[i] + hyd[i]
                if dT >= 0:
                    r = K_G * dT
                    h_new = (hyd[i] / DT + r * total) / (1.0 / DT + r)
                    trace = 1e-9 * RHO_G * DX
                    dh_max = (mg[i] - trace) / (GAS_CONTENT * alpha[i] * DX)
                    h_new = min(h_new, hyd[i] + max(dh_max, 0.0))
                else:
                    h_new = hyd[i] / (1.0 + DT * K_D * (-dT))
                h_new = min(max(h_new, 0.0), total)
                dh = h_new - hyd[i]
                water[i] = min(max(total - h_new, 0.0), 1.0)
                hyd[i] = h_new
                mg[i] -= GAS_CONTENT * alpha[i] * dh * DX
                cap = AREA * DX * (alpha[i] * RHO_L * CL
                                   + ag * rho_g_cell_sealed(mg[i], ag) * CG)
                kw = KW * DX
                temp[i] = (cap / DT * temp[i] + kw * T_WALL
                           + LATENT * alpha[i] * dh * AREA * DX / DT) \
                          / (cap / DT + kw)
            steps_done += 1
        t = fr["t"]
        ok_frame = True
        for i in range(N):
            got = fr["cells"][i]
            want_p = rho_g_cell_sealed(mg[i], 1.0 - alpha[i]) * GAS_R * temp[i]
            if (abs(got["alpha"] - alpha[i]) > 1e-9
                    or abs(got["water"] - water[i]) > 1e-9 + 1e-9 * water[i]
                    or abs(got.get("h", 0.0) - hyd[i]) > 1e-9 + 1e-9 * hyd[i]
                    or abs(got.get("T", 0.0) - temp[i])
                    > 1e-9 + 1e-9 * temp[i]
                    or abs(got.get("p", 0.0) - want_p) > 1e-9 * want_p):
                fails.append(f"cell {i} (t={t}): got={got} want w={water[i]} "
                             f"h={hyd[i]} T={temp[i]} p={want_p}")
                ok_frame = False
                break
        if ok_frame:
            print(f"ok:   water/hydrate/T/p fields replicated at t={t} (1e-9)")
    # closed-form end state: full 1:1 conversion exhausts the water, the
    # stoichiometric sink leaves mg* = mg0 - gasContent a h* dx, the wall
    # pins T and the EOS pins p
    mg_end = RHO_G * (1.0 - A0) * DX - GAS_CONTENT * A0 * W0 * DX
    p_end = mg_end / ((1.0 - A0) * DX) * GAS_R * T_WALL
    c0 = frames[-1]["cells"][0]
    check("end hydrate fraction (full 1:1 conversion)", c0.get("h", 0.0),
          W0, 1e-9)
    check("end water cut (exhausted)", c0["water"], 0.0, 0.0, 1e-9)
    check("end temperature (pinned at wall)", c0.get("T", 0.0), T_WALL, 1e-9)
    check("end EOS pressure (mg* = 139 kg/m2)", c0.get("p", 0.0), p_end, 1e-9)
    if frames[-1]["cells"] != frames[-2]["cells"]:
        fails.append("state still drifting between the last two frames")
    else:
        print("ok:   steady state held over the final 80 s")

elif mode == "hydrateflow":
    # rate-driven flowing line with the EOS, the thermal block and the
    # hydrate kinetics all active: the checker replicates the full coupled
    # step cell by cell - Godunov transport of alpha (wallis F convex,
    # upwind from the left) with water and hydrate riding the same liquid
    # flux under the same clamp semantics, upwind gas-mass continuity on
    # the wallis gas flux G, the implicit hydrate recursion with the
    # stoichiometric gas sink, and the implicit temperature solve with the
    # per-face phase-volume-flux enthalpy conductances
    RHO_G = 200.0  # case fluid rhoG (the shared default above is another case)
    N, DI, DT = 20, 0.02, 0.2
    L = 10.0
    DX = L / N
    AREA = math.pi * DI * DI / 4.0
    VM, US, VSL, VSG, WF_IN = 1.0, 0.3, 0.4, 0.6, 0.25
    GAS_R, P_OUT, T_OUT = 287.0, 1.0e6, 280.0
    CL, CG = 4180.0, 2150.0
    T0, T_IN = 290.0, 290.0
    KW, T_WALL = 50.0, 280.0
    A0, W0 = 0.5, 0.25
    CURVE_A, CURVE_B = 35.2, -5684.0
    K_G, K_D = 0.01, 0.01
    GAS_CONTENT, LATENT = 122.0, 4.2e8
    STEPS_PER_FRAME = 100

    def solve_dense(mat, rhs):
        # port of src/1d/src/DenseSolve.h (Gauss elimination with partial
        # pivoting + back substitution)
        n = len(rhs)
        for k in range(n):
            piv, best = k, abs(mat[k][k])
            for i in range(k + 1, n):
                if abs(mat[i][k]) > best:
                    best, piv = abs(mat[i][k]), i
            if piv != k:
                mat[k], mat[piv] = mat[piv], mat[k]
                rhs[k], rhs[piv] = rhs[piv], rhs[k]
            for i in range(k + 1, n):
                m = mat[i][k] / mat[k][k]
                if m == 0:
                    continue
                for j in range(k, n):
                    mat[i][j] -= m * mat[k][j]
                rhs[i] -= m * rhs[k]
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            s = rhs[i] - sum(mat[i][j] * x[j] for j in range(i + 1, n))
            x[i] = s / mat[i][i]
        return x

    def rho_g_of(m, a):
        ag = 1.0 - a
        return m / (ag * DX) if ag > 1e-12 else RHO_G

    alpha = [A0] * N
    water = [W0] * N
    hyd = [0.0] * N
    temp = [T0] * N
    mg = [RHO_G * (1.0 - A0) * DX] * N
    steps_done = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * STEPS_PER_FRAME
        while steps_done < want_steps:
            # ---- transport ----
            flux = [VSL]
            for j in range(1, N):
                flux.append(godunov(US, VM - US, 0.0, alpha[j - 1], alpha[j]))
            a_end = alpha[N - 1]
            flux.append(US * a_end * a_end + (VM - US) * a_end)
            gvol = [VSG]
            for j in range(1, N):
                gvol.append(godunov(-US, US - VM, VM, alpha[j - 1], alpha[j]))
            gvol.append(VM + (US - VM) * a_end - US * a_end * a_end)
            wflux, hflux = [0.0] * (N + 1), [0.0] * (N + 1)
            for j in range(N + 1):
                if flux[j] == 0:
                    continue
                if j == 0:
                    wflux[j] = WF_IN * flux[j]      # injected composition
                    hflux[j] = 0.0 * flux[j]
                elif j == N:
                    wflux[j] = water[N - 1] * flux[j]
                    hflux[j] = hyd[N - 1] * flux[j]
                else:
                    wu = water[j - 1] if flux[j] >= 0 else water[j]
                    hu = hyd[j - 1] if flux[j] >= 0 else hyd[j]
                    wflux[j] = wu * flux[j]
                    hflux[j] = hu * flux[j]
            gmflux = []
            for j in range(N + 1):
                q = gvol[j]
                if q > 0:
                    rho_up = RHO_G if j == 0 else \
                        rho_g_of(mg[N - 1], alpha[N - 1]) if j == N else \
                        rho_g_of(mg[j - 1], alpha[j - 1])
                elif q < 0:
                    rho_up = P_OUT / (GAS_R * T_OUT) if j == N else \
                        rho_g_of(mg[j], alpha[j])
                else:
                    rho_up = 0.0
                gmflux.append(rho_up * q)
            for i in range(N):
                a_new = alpha[i] + DT / DX * (flux[i] - flux[i + 1])
                a_cap = max(a_new, 0.0)
                af = alpha[i] * water[i] + DT / DX * (wflux[i] - wflux[i + 1])
                ah = alpha[i] * hyd[i] + DT / DX * (hflux[i] - hflux[i + 1])
                alpha[i] = min(max(a_new, 0.0), 1.0)
                af = min(max(af, 0.0), a_cap)
                water[i] = af / a_cap if a_cap > 1e-9 else 0.0
                ah = min(max(ah, 0.0), a_cap)
                hyd[i] = ah / a_cap if a_cap > 1e-9 else 0.0
                mg[i] = max(mg[i] - DT * (gmflux[i + 1] - gmflux[i]), 0.0)
            # ---- hydrate kinetics ----
            h_delta = [0.0] * N
            for i in range(N):
                p_i = (mg[i] * GAS_R * temp[i] / ((1.0 - alpha[i]) * DX)
                       if 1.0 - alpha[i] > 1e-12 else 0.0)
                dT = CURVE_B / (math.log(p_i) - CURVE_A) - temp[i]
                total = water[i] + hyd[i]
                if dT >= 0:
                    r = K_G * dT
                    h_new = (hyd[i] / DT + r * total) / (1.0 / DT + r)
                    trace = 1e-9 * RHO_G * DX
                    dh_max = (mg[i] - trace) / (GAS_CONTENT * alpha[i] * DX)
                    h_new = min(h_new, hyd[i] + max(dh_max, 0.0))
                else:
                    h_new = hyd[i] / (1.0 + DT * K_D * (-dT))
                h_new = min(max(h_new, 0.0), total)
                dh = h_new - hyd[i]
                water[i] = min(max(total - h_new, 0.0), 1.0)
                hyd[i] = h_new
                mg[i] -= GAS_CONTENT * alpha[i] * dh * DX
                h_delta[i] = alpha[i] * dh
            # ---- thermal: implicit solve on phase-flux conductances ----
            mat = [[0.0] * N for _ in range(N)]
            rhs = [0.0] * N
            for i in range(N):
                cap = AREA * DX * (alpha[i] * RHO_L * CL
                                   + (1.0 - alpha[i])
                                   * rho_g_of(mg[i], alpha[i]) * CG)
                kw = KW * DX
                mat[i][i] += cap / DT + kw
                rhs[i] += cap / DT * temp[i] + kw * T_WALL
                rhs[i] += LATENT * h_delta[i] * AREA * DX / DT
            for j in range(N + 1):
                kl = AREA * RHO_L * CL * abs(flux[j])
                q = gvol[j]
                if q > 0:
                    rho_up = RHO_G if j == 0 else \
                        rho_g_of(mg[N - 1], alpha[N - 1]) if j == N else \
                        rho_g_of(mg[j - 1], alpha[j - 1])
                elif q < 0:
                    rho_up = P_OUT / (GAS_R * T_OUT) if j == N else \
                        rho_g_of(mg[j], alpha[j])
                else:
                    rho_up = 0.0
                kg = AREA * rho_up * CG * abs(q)
                for q_f, k_f in ((flux[j], kl), (q, kg)):
                    if q_f == 0:
                        continue
                    up = j - 1 if q_f > 0 else j
                    dn = j if q_f > 0 else j - 1
                    if 0 <= up < N:
                        mat[up][up] += k_f
                    if 0 <= dn < N:
                        if 0 <= up < N:
                            mat[dn][up] -= k_f
                        else:
                            rhs[dn] += k_f * (T_IN if j == 0 else T_OUT)
            temp = solve_dense(mat, rhs)
            steps_done += 1
        t = fr["t"]
        check(f"vm imposed (t={t})", fr["vm"], VM, 1e-12)
        ok_frame = True
        for i in range(N):
            got = fr["cells"][i]
            ag_i = 1.0 - alpha[i]
            want_p = (mg[i] * GAS_R * temp[i] / (ag_i * DX)
                      if ag_i > 1e-12 else 0.0)
            if (abs(got["alpha"] - alpha[i]) > 1e-9 + 1e-9 * alpha[i]
                    or abs(got["water"] - water[i]) > 1e-9 + 1e-9 * water[i]
                    or abs(got.get("h", 0.0) - hyd[i]) > 1e-9 + 1e-9 * hyd[i]
                    or abs(got.get("T", 0.0) - temp[i])
                    > 1e-9 + 1e-9 * temp[i]
                    or abs(got.get("p", 0.0) - want_p) > 1e-9 * want_p):
                fails.append(f"cell {i} (t={t}): got={got} want a={alpha[i]} "
                             f"w={water[i]} h={hyd[i]} T={temp[i]} "
                             f"p={want_p}")
                ok_frame = False
                break
        if ok_frame:
            print(f"ok:   alpha/water/hydrate/T/p replicated at t={t} (1e-9)")
    # identities from emitted fields only
    for i, c in enumerate(frames[-1]["cells"]):
        if abs(c["water"] + c.get("h", 0.0) - WF_IN) > 1e-12:
            fails.append(f"hydrate-able composition cell {i}: "
                         f"w+h={c['water'] + c.get('h', 0.0)} != {WF_IN}")
            break
    else:
        print(f"ok:   per-cell composition w + h = {WF_IN} everywhere "
              f"(1:1 conversion) (1e-12)")
    for i, c in enumerate(frames[-1]["cells"]):
        if abs(c["vsl"] + c["vsg"] - VM) > 1e-12:
            fails.append(f"phase-split identity cell {i}: vsl={c['vsl']} "
                         f"vsg={c['vsg']} sum={c['vsl'] + c['vsg']}")
            break
    else:
        print(f"ok:   per-cell phase split vsl + vsg = vM = {VM} (1e-12)")
    m_emit = sum(c["p"] * (1.0 - c["alpha"]) * DX / (GAS_R * c["T"])
                 for c in frames[-1]["cells"])
    c_out = frames[-1]["cells"][-1]
    rho_g_out = c_out["p"] / (GAS_R * c_out["T"])
    # emitted EOS inventory vs the replicated per-cell gas masses
    check("gas inventory (emitted EOS vs replicated mg)",
          m_emit, sum(mg), 1e-9)
    # steady rate balance from emitted fields only: the net gas inflow
    # (injected rhoG * vsg minus the outlet EOS density * vsg) is the gas
    # the kinetics consumes making hydrate, which then rides the liquid out
    # (formation rate = vsl_out * h_out - vsl_in * h_in, h_in = 0)
    check("steady gas balance (net inflow = hydrate consumption)",
          RHO_G * VSG - rho_g_out * c_out["vsg"],
          GAS_CONTENT * c_out["vsl"] * c_out["h"], 1e-9)
    dp_want = 0.0
    for c in frames[-1]["cells"]:
        rho_g = c["p"] / (GAS_R * c["T"])
        rho_m = c["alpha"] * RHO_L + (1.0 - c["alpha"]) * rho_g
        mu_m = c["alpha"] * MU_L + (1.0 - c["alpha"]) * MU_G
        re = rho_m * VM * DI / mu_m
        dp_want += churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM
    check("pressure drop (Churchill friction, term by term)",
          frames[-1]["dp"], dp_want, 1e-9)
    qwall_want = sum(KW * DX * (T_WALL - c["T"])
                     for c in frames[-1]["cells"])
    check("wall heat sum", frames[-1].get("qwall", 0.0), qwall_want, 1e-9)
    # steady energy balance: inflow enthalpy + wall + latent release equals
    # the outflow enthalpy carried by the outlet-face phase fluxes; the
    # latent source is the volumetric formation rate (vsl_out * h_out per
    # unit area) times the pipe area times the latent heat
    kl_out = AREA * RHO_L * CL * c_out["vsl"]
    kg_out = AREA * rho_g_out * CG * c_out["vsg"]
    h_out = (kl_out + kg_out) * c_out["T"]
    latent = LATENT * AREA * c_out["vsl"] * c_out["h"]
    check("steady energy balance (inlet + wall + latent = outflow)",
          h_out, frames[-1].get("qinlet", 0.0) + frames[-1].get("qwall", 0.0)
          + latent, 1e-9)
    # physics shape: water depletes and hydrate accumulates along x
    ws = [c["water"] for c in frames[-1]["cells"]]
    hs = [c.get("h", 0.0) for c in frames[-1]["cells"]]
    if any(ws[i] <= ws[i + 1] for i in range(N - 1)) or \
            any(hs[i] >= hs[i + 1] for i in range(N - 1)):
        fails.append("conversion profile not monotone along x")
    else:
        print("ok:   water depletes / hydrate accumulates monotonically")
    if frames[-1]["cells"] != frames[-2]["cells"]:
        fails.append("state still drifting between the last two frames")
    else:
        print("ok:   steady state held over the final 20 s")

elif mode == "junction":
    # symmetric Y junction (case.dict constants): two pressure-driven
    # branches merging at J into a common outlet line, uniform holdup and
    # wallis slip 0 - the mixture properties are constant, so the coupled
    # transient is the conductance network below
    RHO_L, RHO_G, MU_L, MU_G = 1000.0, 20.0, 0.01, 1e-5
    CL, CG = 4180.0, 2150.0
    ALPHA = 0.5
    RHO_M = ALPHA * RHO_L + (1.0 - ALPHA) * RHO_G
    MU_M = ALPHA * MU_L + (1.0 - ALPHA) * MU_G
    NA, N_OUT = 10, 8
    LA, LB, LC = 10.0, 10.0, 8.0
    DA, DB, DC = 0.02, 0.02, 0.03
    AA = math.pi * DA * DA / 4.0
    AB = math.pi * DB * DB / 4.0
    AC = math.pi * DC * DC / 4.0
    DXA, DXC = LA / NA, LC / N_OUT
    PA, PB, PC = 200.0, 300.0, 0.0
    WF_A, WF_B, T_A, T_B = 0.5, 0.1, 300.0, 320.0
    T0 = 300.0
    KA = 32.0 * MU_M * LA / DA ** 2
    KB = 32.0 * MU_M * LB / DB ** 2
    KC = 32.0 * MU_M * LC / DC ** 2
    IA, IB, IC = RHO_M * LA, RHO_M * LB, RHO_M * LC
    DT = 0.1
    N_STEPS = 12000
    REPORT = 100

    # steady closed form: conductance-weighted junction pressure
    cond_a, cond_b, cond_c = AA / KA, AB / KB, AC / KC
    pj_star = (cond_a * PA + cond_b * PB + cond_c * PC) \
        / (cond_a + cond_b + cond_c)
    ua_star = (PA - pj_star) / KA
    ub_star = (PB - pj_star) / KB
    uc_star = (AA * ua_star + AB * ub_star) / AC
    # inflow mixing: liquid-volume-flow-weighted water cut,
    # heat-capacity-flow-weighted temperature
    qla = AA * ua_star * ALPHA
    qlb = AB * ub_star * ALPHA
    qga = AA * ua_star * (1.0 - ALPHA)
    qgb = AB * ub_star * (1.0 - ALPHA)
    phi_star = (qla * WF_A + qlb * WF_B) / (qla + qlb)
    e_a = qla * RHO_L * CL + qga * RHO_G * CG
    e_b = qlb * RHO_L * CL + qgb * RHO_G * CG
    tj_star = (e_a * T_A + e_b * T_B) / (e_a + e_b)

    # independent trajectory: per step the eliminated 2x2 laminar
    # implicit-Euler momentum system, then the junction mixing and the
    # transport/thermal replication with the new velocities
    def advect_water(w, dx, f, w_in):
        # water rides the uniform liquid flux f (holdup stays ALPHA);
        # upwind Godunov with the ghost value feeding face 0 and the
        # outlet cell feeding face n, then the kernel's clamp semantics
        n = len(w)
        wf = [w_in * f] + [w[j - 1] * f for j in range(1, n)] \
            + [w[n - 1] * f]
        out = []
        for i in range(n):
            afn = ALPHA * w[i] + DT / dx * (wf[i] - wf[i + 1])
            afn = min(max(afn, 0.0), ALPHA)
            out.append(afn / ALPHA)
        return out

    def thermal(T, area, dx, u, t_in):
        # implicit upwind enthalpy advection: the donor cell drains at its
        # own (new) temperature, the receiver gains the donor's (new)
        # temperature, face 0 is fed by the ghost reservoir
        n = len(T)
        cap = area * dx * (ALPHA * RHO_L * CL
                           + (1.0 - ALPHA) * RHO_G * CG)
        kf = area * u * (RHO_L * CL * ALPHA + RHO_G * CG * (1.0 - ALPHA))
        diag = [cap / DT] * n
        rhs = [cap / DT * T[i] for i in range(n)]
        rhs[0] += kf * t_in
        for i in range(n - 1):
            diag[i] += kf
        diag[n - 1] += kf
        # lower-bidiagonal system: forward substitution
        out = [rhs[0] / diag[0]]
        for i in range(1, n):
            out.append((rhs[i] + kf * out[i - 1]) / diag[i])
        return out

    u_a_old = u_b_old = u_c_old = 0.0
    w_a = [0.5] * NA
    w_b = [0.5] * NA
    w_c = [0.5] * N_OUT
    t_a = [T0] * NA
    t_b = [T0] * NA
    t_c = [T0] * N_OUT
    traj = [{"t": 0.0, "uA": 0.0, "uB": 0.0, "uC": 0.0, "pJ": pj_star,
             "wA": w_a[:], "wB": w_b[:], "wC": w_c[:],
             "TA": t_a[:], "TB": t_b[:], "TC": t_c[:]}]
    for step in range(1, N_STEPS + 1):
        ac_i = IC / DT + KC
        k = ac_i / AC
        a11 = IA / DT + KA + k * AA
        a12 = k * AB
        a21 = k * AA
        a22 = IB / DT + KB + k * AB
        b1 = PA - PC + IA / DT * u_a_old + IC / DT * u_c_old
        b2 = PB - PC + IB / DT * u_b_old + IC / DT * u_c_old
        det = a11 * a22 - a12 * a21
        u_a = (b1 * a22 - b2 * a12) / det
        u_b = (a11 * b2 - a21 * b1) / det
        u_c = (AA * u_a + AB * u_b) / AC
        p_j = PC + IC * (u_c - u_c_old) / DT + KC * u_c

        # junction mix from the feeding branches' end cells at the new
        # velocities, before the transport advances the fields
        ala = AA * u_a * ALPHA
        alb = AB * u_b * ALPHA
        aga = AA * u_a * (1.0 - ALPHA)
        agb = AB * u_b * (1.0 - ALPHA)
        phi = (ala * w_a[-1] + alb * w_b[-1]) / (ala + alb)
        ea = ala * RHO_L * CL + aga * RHO_G * CG
        eb = alb * RHO_L * CL + agb * RHO_G * CG
        tj = (ea * t_a[-1] + eb * t_b[-1]) / (ea + eb)

        w_a = advect_water(w_a, DXA, u_a * ALPHA, WF_A)
        w_b = advect_water(w_b, DXA, u_b * ALPHA, WF_B)
        w_c = advect_water(w_c, DXC, u_c * ALPHA, phi)
        t_a = thermal(t_a, AA, DXA, u_a, T_A)
        t_b = thermal(t_b, AB, DXA, u_b, T_B)
        t_c = thermal(t_c, AC, DXC, u_c, tj)
        u_a_old, u_b_old, u_c_old = u_a, u_b, u_c
        if step % REPORT == 0:
            traj.append({"t": step * DT, "uA": u_a, "uB": u_b, "uC": u_c,
                         "pJ": p_j, "wA": w_a[:], "wB": w_b[:],
                         "wC": w_c[:], "TA": t_a[:], "TB": t_b[:],
                         "TC": t_c[:]})

    if len(frames) != len(traj):
        fails.append(f"frame count {len(frames)} != {len(traj)}")
    else:
        print(f"ok:   {len(frames)} frames at reportEvery {REPORT}")
    geo = {"branchA": (NA, DXA, DA, AA, KA), "branchB": (NA, DXA, DB, AB, KB),
           "outlet": (N_OUT, DXC, DC, AC, KC)}
    for fr, tr in zip(frames, traj):
        t = fr["t"]
        jn = fr["junctions"]["J"]
        # the Newton balance row is exact: net = sum s A u at the
        # converged velocities (tolerance 1e-12 * rate scale ~ 1e-15)
        check(f"junction net balance (t={t})", jn["net"], 0.0, 0.0, 1e-12)
        # junction pressure + branch velocities vs the independent
        # laminar trajectory (rtol covers the Churchill deviation)
        check(f"p_J (t={t})", jn["p"], tr["pJ"], 1e-4)
        for pn, ukey in (("branchA", "uA"), ("branchB", "uB"),
                         ("outlet", "uC")):
            n, dx, d, a, kk = geo[pn]
            u = tr[ukey]
            check(f"vm {pn} (t={t})", fr["pipes"][pn]["vm"], u, 1e-4)
            # the emitted dp is the friction sum alone (the inertia lives
            # in p_J, not in the loss-marched node pressures)
            re = RHO_M * abs(u) * d / MU_M
            dp_want = n * churchill(re, 0.0) * (dx / d) * 0.5 * RHO_M * u * u
            check(f"dp {pn} (t={t})", fr["pipes"][pn]["dp"], dp_want, 1e-9)
        for pn, wkey, tkey in (("branchA", "wA", "TA"), ("branchB", "wB", "TB"),
                               ("outlet", "wC", "TC")):
            cells = fr["pcell"][pn]
            for i, c in enumerate(cells):
                if abs(c["alpha"] - ALPHA) > 1e-9:
                    fails.append(f"{pn} alpha[{i}](t={t}) drifted")
                    break
                if abs(c["vsl"] + c["vsg"] - fr["pipes"][pn]["vm"]) > 1e-9:
                    fails.append(f"{pn} phase split [{i}](t={t})")
                    break
                if abs(c["water"] - tr[wkey][i]) > 1e-9:
                    fails.append(f"{pn} water[{i}](t={t}): got={c['water']}"
                                 f" want={tr[wkey][i]}")
                    break
                # the trajectory's velocities carry the ~1e-6
                # Churchill-vs-laminar deviation, so the implicit thermal
                # solve can only be matched to ~1e-9 (atol 1e-7)
                if abs(c["T"] - tr[tkey][i]) > 1e-7:
                    fails.append(f"{pn} T[{i}](t={t}): got={c['T']}"
                                 f" want={tr[tkey][i]}")
                    break
            else:
                print(f"ok:   {pn} alpha/water/T fields (t={t}) "
                      f"(1e-9 / 1e-9 / 1e-7)")

    # steady anchors on the closed form (last frame)
    fr = frames[-1]
    check("steady p_J*", fr["junctions"]["J"]["p"], pj_star, 1e-4)
    check("steady vm branchA*", fr["pipes"]["branchA"]["vm"], ua_star, 1e-4)
    check("steady vm branchB*", fr["pipes"]["branchB"]["vm"], ub_star, 1e-4)
    check("steady vm outlet*", fr["pipes"]["outlet"]["vm"], uc_star, 1e-4)
    # the loss-marched node pressures meet at the junction at steady state
    check("branchA far node = p_J*", fr["pnode"]["branchA"][-1], pj_star,
          1e-4)
    check("branchB far node = p_J*", fr["pnode"]["branchB"][-1], pj_star,
          1e-4)
    check("outlet head node = p_J*", fr["pnode"]["outlet"][0], pj_star, 1e-4)
    # front-flushed composition/temperature of the supplied line
    check("steady phi_J* (outlet water)", fr["pcell"]["outlet"][-1]["water"],
          phi_star, 1e-4)
    check("steady T_J* (outlet temperature)",
          fr["pcell"]["outlet"][-1]["T"], tj_star, 1e-4)
    check("branchA water undisturbed", fr["pcell"]["branchA"][-1]["water"],
          WF_A, 1e-9)
    check("branchB water flushed", fr["pcell"]["branchB"][-1]["water"],
          WF_B, 1e-9)
    check("branchA temperature undisturbed",
          fr["pcell"]["branchA"][-1]["T"], T_A, 1e-9)
    check("branchB temperature flushed", fr["pcell"]["branchB"][-1]["T"],
          T_B, 1e-9)
    drift = []
    for pn in ("branchA", "branchB", "outlet"):
        p1, p0 = fr["pipes"][pn], frames[-2]["pipes"][pn]
        if abs(p1["vm"] - p0["vm"]) > 1e-9 or abs(p1["dp"] - p0["dp"]) > 1e-9:
            drift.append(pn)
        for c1, c0 in zip(fr["pcell"][pn], frames[-2]["pcell"][pn]):
            for key in ("alpha", "water", "vsl", "vsg", "T"):
                if abs(c1[key] - c0[key]) > 1e-9:
                    drift.append(f"{pn}.{key}")
    if drift:
        fails.append(f"state still drifting between the last two frames: "
                     f"{drift[:4]}")
    else:
        print("ok:   steady state held over the final 10 s")

elif mode == "networkspan":
    # two-junction pipe span (case.dict constants): two pressure-driven
    # sources merging at J1, the span `mid` with junctions at BOTH ends
    # (no boundary of its own), then a split at J2 into two pressure-driven
    # sinks - uniform holdup and wallis slip 0 keep the mixture properties
    # constant, so the coupled transient is the resistance network below:
    # one global system over the five velocities and the two node pressures
    RHO_L, RHO_G, MU_L, MU_G = 1000.0, 20.0, 0.01, 1e-5
    CL, CG = 4180.0, 2150.0
    ALPHA = 0.5
    RHO_M = ALPHA * RHO_L + (1.0 - ALPHA) * RHO_G
    MU_M = ALPHA * MU_L + (1.0 - ALPHA) * MU_G
    DT = 0.1
    N_STEPS = 25000
    REPORT = 100
    WF_A, T_A = 0.5, 300.0
    WF_B, T_B = 0.1, 320.0
    T0 = 300.0
    PA, PB = 200.0, 300.0
    # pipe -> (sections, length, diameter); all bores horizontal
    geo = {"srcA": (10, 10.0, 0.02), "srcB": (10, 10.0, 0.02),
           "mid": (6, 6.0, 0.02), "sink1": (4, 4.0, 0.03),
           "sink2": (2, 2.0, 0.02)}
    area, kres, iner, dxc = {}, {}, {}, {}
    for pn, (n, length, d) in geo.items():
        a = math.pi * d * d / 4.0
        area[pn] = a
        kres[pn] = 32.0 * MU_M * length / d ** 2
        iner[pn] = RHO_M * length
        dxc[pn] = length / n
    AAS, AM = area["srcA"], area["mid"]
    A1, A2 = area["sink1"], area["sink2"]
    order = ["srcA", "srcB", "mid", "sink1", "sink2"]

    # steady closed form: nodal conductance solution (c = A/K)
    cond = {pn: area[pn] / kres[pn] for pn in geo}
    cs, cm = cond["srcA"], cond["mid"]
    c1, c2 = cond["sink1"], cond["sink2"]
    beta = cm / (cm + c1 + c2)
    p1_star = (cs * PA + cs * PB) / (2.0 * cs + cm * (1.0 - beta))
    p2_star = beta * p1_star
    u_star = {}
    u_star["srcA"] = (PA - p1_star) / kres["srcA"]
    u_star["srcB"] = (PB - p1_star) / kres["srcB"]
    u_star["mid"] = (p1_star - p2_star) / kres["mid"]
    u_star["sink1"] = p2_star / kres["sink1"]
    u_star["sink2"] = p2_star / kres["sink2"]
    # J1 mixes the inflows: liquid-volume-flow-weighted water cut,
    # heat-capacity-flow-weighted temperature (the span flushes to that
    # state and J2, fed by the span alone, hands it to both sinks)
    qla = area["srcA"] * u_star["srcA"] * ALPHA
    qlb = area["srcB"] * u_star["srcB"] * ALPHA
    qga = area["srcA"] * u_star["srcA"] * (1.0 - ALPHA)
    qgb = area["srcB"] * u_star["srcB"] * (1.0 - ALPHA)
    phi_star = (qla * WF_A + qlb * WF_B) / (qla + qlb)
    e_a = qla * RHO_L * CL + qga * RHO_G * CG
    e_b = qlb * RHO_L * CL + qgb * RHO_G * CG
    t_star = (e_a * T_A + e_b * T_B) / (e_a + e_b)

    def gauss(Amat, b):
        # dense Gaussian elimination with partial pivoting (no numpy)
        n = len(b)
        M = [list(Amat[i]) + [b[i]] for i in range(n)]
        for col in range(n):
            piv = max(range(col, n), key=lambda rr: abs(M[rr][col]))
            if abs(M[piv][col]) < 1e-300:
                raise SystemExit("singular network system")
            M[col], M[piv] = M[piv], M[col]
            pv = M[col][col]
            for rr in range(n):
                if rr != col and M[rr][col] != 0.0:
                    f = M[rr][col] / pv
                    for cc in range(col, n + 1):
                        M[rr][cc] -= f * M[col][cc]
        return [M[i][n] / M[i][i] for i in range(n)]

    def advect_water(w, dx, f, w_in):
        # water rides the uniform liquid flux f (holdup stays ALPHA);
        # upwind Godunov with the ghost value feeding face 0 and the
        # outlet cell feeding face n, then the kernel's clamp semantics
        n = len(w)
        wf = [w_in * f] + [w[j - 1] * f for j in range(1, n)] \
            + [w[n - 1] * f]
        out = []
        for i in range(n):
            afn = ALPHA * w[i] + DT / dx * (wf[i] - wf[i + 1])
            afn = min(max(afn, 0.0), ALPHA)
            out.append(afn / ALPHA)
        return out

    def thermal(T, areaa, dx, u, t_in):
        # implicit upwind enthalpy advection: the donor cell drains at its
        # own (new) temperature, the receiver gains the donor's (new)
        # temperature, face 0 is fed by the ghost reservoir
        n = len(T)
        cap = areaa * dx * (ALPHA * RHO_L * CL
                            + (1.0 - ALPHA) * RHO_G * CG)
        kf = areaa * u * (RHO_L * CL * ALPHA + RHO_G * CG * (1.0 - ALPHA))
        diag = [cap / DT] * n
        rhs = [cap / DT * T[i] for i in range(n)]
        rhs[0] += kf * t_in
        for i in range(n - 1):
            diag[i] += kf
        diag[n - 1] += kf
        # lower-bidiagonal system: forward substitution
        out = [rhs[0] / diag[0]]
        for i in range(1, n):
            out.append((rhs[i] + kf * out[i - 1]) / diag[i])
        return out

    # per step: the 7x7 laminar system [uA uB uC u1 u2 p1 p2] - one
    # implicit-Euler momentum row per pipe (the span driven by p1 - p2)
    # plus one volume-balance row per node; then the node mixes and the
    # transport/thermal replication with the new velocities
    u_old = {pn: 0.0 for pn in order}
    w = {pn: [0.5] * geo[pn][0] for pn in order}
    tfield = {pn: [T0] * geo[pn][0] for pn in order}
    traj = [{"t": 0.0, "u": {pn: 0.0 for pn in order},
             "p1": 250.0, "p2": 0.0,   # kernel warm start (conductance
                                       # mean of each node's far pressures)
             "w": {pn: w[pn][:] for pn in order},
             "T": {pn: tfield[pn][:] for pn in order}}]
    # momentum rows: (I/dt + K) u_i + s_in p1 + s_out p2 = p_far + I/dt u_old
    #   srcA: driven by PA - p1  -> +p1
    #   srcB: driven by PB - p1  -> +p1
    #   mid:  driven by p1 - p2  -> -p1 +p2   (the span couples to both)
    #   sink1/2: driven by p2 - 0 -> -p2
    rows = [("srcA", 1.0, 0.0, PA), ("srcB", 1.0, 0.0, PB),
            ("mid", -1.0, 1.0, 0.0), ("sink1", 0.0, -1.0, 0.0),
            ("sink2", 0.0, -1.0, 0.0)]
    for step in range(1, N_STEPS + 1):
        Amat = [[0.0] * 7 for _ in range(7)]
        b = [0.0] * 7
        for k, (pn, s1, s2, pfar) in enumerate(rows):
            Amat[k][k] = iner[pn] / DT + kres[pn]
            Amat[k][5] = s1
            Amat[k][6] = s2
            b[k] = pfar + iner[pn] / DT * u_old[pn]
        # balance rows: sum s A u = 0 (s = +1 when the node is at the
        # pipe's outlet)
        Amat[5] = [AAS, AAS, -AM, 0.0, 0.0, 0.0, 0.0]
        Amat[6] = [0.0, 0.0, AM, -A1, -A2, 0.0, 0.0]
        x = gauss(Amat, b)
        uA, uB, uC, u1, u2, p1, p2 = x
        unew = {"srcA": uA, "srcB": uB, "mid": uC,
                "sink1": u1, "sink2": u2}

        # node mixes from the feeding branches' end cells at the new
        # velocities, before the transport advances the fields
        ala = AAS * uA * ALPHA
        alb = AAS * uB * ALPHA
        aga = AAS * uA * (1.0 - ALPHA)
        agb = AAS * uB * (1.0 - ALPHA)
        phi1 = (ala * w["srcA"][-1] + alb * w["srcB"][-1]) / (ala + alb)
        ea = ala * RHO_L * CL + aga * RHO_G * CG
        eb = alb * RHO_L * CL + agb * RHO_G * CG
        t1 = (ea * tfield["srcA"][-1] + eb * tfield["srcB"][-1]) / (ea + eb)
        # J2 is fed by the span alone: the mix is the span's end state
        phi2 = w["mid"][-1]
        t2 = tfield["mid"][-1]

        w["srcA"] = advect_water(w["srcA"], dxc["srcA"], uA * ALPHA, WF_A)
        w["srcB"] = advect_water(w["srcB"], dxc["srcB"], uB * ALPHA, WF_B)
        w["mid"] = advect_water(w["mid"], dxc["mid"], uC * ALPHA, phi1)
        w["sink1"] = advect_water(w["sink1"], dxc["sink1"], u1 * ALPHA,
                                  phi2)
        w["sink2"] = advect_water(w["sink2"], dxc["sink2"], u2 * ALPHA,
                                  phi2)
        tfield["srcA"] = thermal(tfield["srcA"], AAS, dxc["srcA"], uA, T_A)
        tfield["srcB"] = thermal(tfield["srcB"], AAS, dxc["srcB"], uB, T_B)
        tfield["mid"] = thermal(tfield["mid"], AM, dxc["mid"], uC, t1)
        tfield["sink1"] = thermal(tfield["sink1"], A1, dxc["sink1"], u1,
                                  t2)
        tfield["sink2"] = thermal(tfield["sink2"], A2, dxc["sink2"], u2,
                                  t2)
        u_old = unew
        if step % REPORT == 0:
            traj.append({"t": step * DT, "u": dict(unew), "p1": p1,
                         "p2": p2, "w": {pn: w[pn][:] for pn in order},
                         "T": {pn: tfield[pn][:] for pn in order}})

    if len(frames) != len(traj):
        fails.append(f"frame count {len(frames)} != {len(traj)}")
    else:
        print(f"ok:   {len(frames)} frames at reportEvery {REPORT}")
    for fr, tr in zip(frames, traj):
        t = fr["t"]
        # the Newton balance rows are exact: net = sum s A u at the
        # converged velocities (tolerance 1e-12 * rate scale ~ 1e-15)
        for jn_name, want_p in (("J1", tr["p1"]), ("J2", tr["p2"])):
            jn = fr["junctions"][jn_name]
            check(f"{jn_name} net balance (t={t})", jn["net"], 0.0, 0.0,
                  1e-12)
            check(f"{jn_name} pressure (t={t})", jn["p"], want_p, 1e-4)
        for pn in order:
            n, length, d = geo[pn]
            u = tr["u"][pn]
            check(f"vm {pn} (t={t})", fr["pipes"][pn]["vm"], u, 1e-4)
            # the emitted dp is the friction sum alone (the inertia lives
            # in the node pressures, not in the loss-marched profiles);
            # evaluated at the KERNEL's converged vm: the checker's
            # laminar-trajectory velocity deviates from it by the ~1e-9
            # Churchill-vs-laminar difference, which the dp ~ u^2
            # amplification would smear past a 1e-9 march tolerance (the
            # velocity physics itself is covered by the 1e-4 vm check)
            vmk = fr["pipes"][pn]["vm"]
            re = RHO_M * abs(vmk) * d / MU_M
            dp_want = n * churchill(re, 0.0) * (dxc[pn] / d) * 0.5 \
                * RHO_M * vmk * vmk
            check(f"dp {pn} (t={t})", fr["pipes"][pn]["dp"], dp_want, 1e-9)
            cells = fr["pcell"][pn]
            for i, c in enumerate(cells):
                if abs(c["alpha"] - ALPHA) > 1e-9:
                    fails.append(f"{pn} alpha[{i}](t={t}) drifted")
                    break
                if abs(c["vsl"] + c["vsg"] - fr["pipes"][pn]["vm"]) > 1e-9:
                    fails.append(f"{pn} phase split [{i}](t={t})")
                    break
                if abs(c["water"] - tr["w"][pn][i]) > 1e-9:
                    fails.append(f"{pn} water[{i}](t={t}): got={c['water']}"
                                 f" want={tr['w'][pn][i]}")
                    break
                # the trajectory's velocities carry the ~1e-6
                # Churchill-vs-laminar deviation, so the implicit thermal
                # solve can only be matched to ~1e-9 (atol 1e-7)
                if abs(c["T"] - tr["T"][pn][i]) > 1e-7:
                    fails.append(f"{pn} T[{i}](t={t}): got={c['T']}"
                                 f" want={tr['T'][pn][i]}")
                    break
            else:
                print(f"ok:   {pn} alpha/water/T fields (t={t}) "
                      f"(1e-9 / 1e-9 / 1e-7)")

    # steady anchors on the closed form (last frame)
    fr = frames[-1]
    check("steady p1*", fr["junctions"]["J1"]["p"], p1_star, 1e-4)
    check("steady p2*", fr["junctions"]["J2"]["p"], p2_star, 1e-4)
    for pn in order:
        check(f"steady vm {pn}*", fr["pipes"][pn]["vm"], u_star[pn], 1e-4)
    # the loss-marched node pressures meet at the junctions at steady
    # state (src pipes anchored at their far inlet, the span and the sinks
    # at their far outlet / imposed rates)
    check("srcA far node = p1*", fr["pnode"]["srcA"][-1], p1_star, 1e-4)
    check("srcB far node = p1*", fr["pnode"]["srcB"][-1], p1_star, 1e-4)
    check("mid head node = p1*", fr["pnode"]["mid"][0], p1_star, 1e-4)
    check("mid tail node = p2*", fr["pnode"]["mid"][-1], p2_star, 1e-4)
    check("sink1 head node = p2*", fr["pnode"]["sink1"][0], p2_star, 1e-4)
    check("sink2 head node = p2*", fr["pnode"]["sink2"][0], p2_star, 1e-4)
    # front-flushed composition/temperature: srcA undisturbed, srcB
    # flushed, the span carries the J1 mix and hands it to both sinks
    check("steady phi1* (mid water)", fr["pcell"]["mid"][-1]["water"],
          phi_star, 1e-4)
    check("steady T1* (mid temperature)", fr["pcell"]["mid"][-1]["T"],
          t_star, 1e-4)
    check("sink1 water = phi1*", fr["pcell"]["sink1"][-1]["water"],
          phi_star, 1e-4)
    check("sink2 water = phi1*", fr["pcell"]["sink2"][-1]["water"],
          phi_star, 1e-4)
    check("sink1 temperature = T1*", fr["pcell"]["sink1"][-1]["T"],
          t_star, 1e-4)
    check("sink2 temperature = T1*", fr["pcell"]["sink2"][-1]["T"],
          t_star, 1e-4)
    check("srcA water undisturbed", fr["pcell"]["srcA"][-1]["water"],
          WF_A, 1e-9)
    check("srcB water flushed", fr["pcell"]["srcB"][-1]["water"],
          WF_B, 1e-9)
    check("srcA temperature undisturbed",
          fr["pcell"]["srcA"][-1]["T"], T_A, 1e-9)
    check("srcB temperature flushed", fr["pcell"]["srcB"][-1]["T"],
          T_B, 1e-9)
    drift = []
    for pn in order:
        pn1, pn0 = fr["pipes"][pn], frames[-2]["pipes"][pn]
        if abs(pn1["vm"] - pn0["vm"]) > 1e-9 \
                or abs(pn1["dp"] - pn0["dp"]) > 1e-9:
            drift.append(pn)
        for c1f, c0f in zip(fr["pcell"][pn], frames[-2]["pcell"][pn]):
            for key in ("alpha", "water", "vsl", "vsg", "T"):
                if abs(c1f[key] - c0f[key]) > 1e-9:
                    drift.append(f"{pn}.{key}")
    if drift:
        fails.append(f"state still drifting between the last two frames: "
                     f"{drift[:4]}")
    else:
        print("ok:   steady state held over the final 10 s")

elif mode == "entrained":
    # film-droplet split first increment: rate-driven wallis pipe, the
    # droplet fraction d of the liquid rides the gas velocity and relaxes
    # towards E_eq = clamp((vG - vCrit)/vRef, 0, 1) at the exchange time
    N, L, DI, DT = 20, 20.0, 0.02, 0.1
    DX = L / N
    VSL, VSG = 0.15, 0.25
    VM = VSL + VSG
    US = 0.3
    PHI_IN, PHI0 = 0.25, 0.8
    E_IN, E0 = 0.1, 0.9
    VCRIT, VREF, TAU = 0.2, 0.6, 8.0
    ROUGH = 0.0
    REPORT = 100
    qa, qb, qc = US, VM - US, 0.0           # F(a) = us a^2 + (vm-us) a
    ga, gb, gc = -US, US - VM, VM           # G(a) = -us a^2 + (us-vm) a + vm
    a_star = (-(VM - US) + math.sqrt((VM - US) ** 2 + 4 * US * VSL)) / (2 * US)
    vg_star = VM + US * a_star
    d_star = min(max((vg_star - VCRIT) / VREF, 0.0), 1.0)

    def vg_of(a):
        return VM + US * a

    alpha = [0.5] * N
    water = [PHI0] * N
    entr = [E0] * N
    step_in_sample = 0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * REPORT
        while step_in_sample < want_steps:
            flux = ([VSL]
                    + [godunov(qa, qb, qc, alpha[j - 1], alpha[j])
                       for j in range(1, N)]
                    + [qc + alpha[-1] * (qb + qa * alpha[-1])])
            gvol = ([VSG]
                    + [godunov(ga, gb, gc, alpha[j - 1], alpha[j])
                       for j in range(1, N)]
                    + [gc + alpha[-1] * (gb + ga * alpha[-1])])
            wflux = [PHI_IN * flux[0]]
            for j in range(1, N):
                wflux.append((water[j - 1] if flux[j] >= 0 else water[j])
                             * flux[j])
            wflux.append(water[-1] * flux[N])
            # droplet content flux: rate face imposes the inflow share of
            # the liquid rate, gas-flux faces upwind on the gas direction
            dflux = [E_IN * flux[0]]
            for j in range(1, N):
                if gvol[j] > 0:
                    dflux.append(alpha[j - 1] * entr[j - 1]
                                 * vg_of(alpha[j - 1]))
                elif gvol[j] < 0:
                    dflux.append(alpha[j] * entr[j] * vg_of(alpha[j]))
                else:
                    dflux.append(0.0)
            dflux.append(alpha[-1] * entr[-1] * vg_of(alpha[-1])
                         if gvol[N] > 0 else 0.0)
            new = []
            for i in range(N):
                na = alpha[i] + DT / DX * (flux[i] - flux[i + 1])
                naf = (alpha[i] * water[i]
                       + DT / DX * (wflux[i] - wflux[i + 1]))
                nqe = (alpha[i] * entr[i]
                       + DT / DX * (dflux[i] - dflux[i + 1]))
                cap = max(na, 0.0)
                na = min(max(na, 0.0), 1.0)
                naf = min(max(naf, 0.0), cap)
                nqe = min(max(nqe, 0.0), cap)
                new.append((na, naf / cap if cap > 1e-9 else 0.0,
                            nqe / cap if cap > 1e-9 else 0.0))
            alpha = [v[0] for v in new]
            water = [v[1] for v in new]
            entr = [v[2] for v in new]
            for i in range(N):
                if alpha[i] > 1e-9:
                    e_eq = min(max((vg_of(alpha[i]) - VCRIT) / VREF,
                                   0.0), 1.0)
                    entr[i] += (e_eq - entr[i]) * DT / (TAU + DT)
            step_in_sample += 1
        t = fr["t"]
        check(f"vm(t={t})", fr["vm"], VM, 1e-12)
        for i in range(N):
            c = fr["cells"][i]
            if abs(c["alpha"] - alpha[i]) > 1e-9 + 1e-9 * abs(alpha[i]):
                fails.append(f"alpha[{i}](t={t}): got={c['alpha']} "
                             f"want={alpha[i]}")
                break
            if abs(c["water"] - water[i]) > 1e-9 + 1e-9 * abs(water[i]):
                fails.append(f"water[{i}](t={t}): got={c['water']} "
                             f"want={water[i]}")
                break
            if abs(c["entr"] - entr[i]) > 1e-9 + 1e-9 * abs(entr[i]):
                fails.append(f"entr[{i}](t={t}): got={c['entr']} "
                             f"want={entr[i]}")
                break
            # the superficial phase fluxes are F(alpha), G(alpha): they sum
            # to the mixture velocity exactly and pin the split per cell
            if abs(c["vsl"] - (qc + alpha[i] * (qb + qa * alpha[i]))) > 1e-9:
                fails.append(f"vsl[{i}](t={t}): got={c['vsl']}")
                break
            if abs(c["vsg"] - (gc + alpha[i] * (gb + ga * alpha[i]))) > 1e-9:
                fails.append(f"vsg[{i}](t={t}): got={c['vsg']}")
                break
            if abs(c["vsl"] + c["vsg"] - VM) > 1e-9:
                fails.append(f"vsl+vsg != vm at cell {i} (t={t})")
                break
        else:
            print(f"ok:   alpha/water/entr fields replicated (t={t}) "
                  f"(1e-9)")
        # pressure drop: per-cell Churchill friction on the replicated
        # holdup (rate-driven anchors the gauge at the outlet)
        dp = 0.0
        for i in range(N):
            rho = alpha[i] * RHO_L + (1.0 - alpha[i]) * RHO_G
            mu = alpha[i] * MU_L + (1.0 - alpha[i]) * MU_G
            re = rho * VM * DI / mu
            dp += churchill(re, ROUGH / DI) * (DX / DI) * 0.5 * rho * VM ** 2
        check(f"dp = Churchill friction sum (t={t})", fr["dp"], dp, 1e-9)

    # steady anchors on the closed forms (last frame)
    fr = frames[-1]
    hold_bad = [i for i in range(N)
                if abs(fr["cells"][i]["alpha"] - a_star)
                > 1e-9 + 1e-9 * a_star]
    if hold_bad:
        fails.append(f"steady holdup cells {hold_bad[:4]} != alpha*")
    else:
        print(f"ok:   uniform steady holdup alpha*={a_star:.10f} (1e-9)")
    # no uniform droplet fraction is a fixed point: the inlet injects
    # E_IN vsl of droplet content while the interior carries alpha d vG,
    # so the steady state is the developing profile — the exact discrete
    # fixed point of transport+relaxation (clamps inactive there)
    dx_tau = DX / TAU
    e_eq_star = d_star
    d_prof = [(E_IN * VSL + a_star * dx_tau * e_eq_star)
              / (a_star * (vg_star + dx_tau))]
    for _ in range(1, N):
        d_prof.append((vg_star * d_prof[-1] + dx_tau * e_eq_star)
                      / (vg_star + dx_tau))
    ent_bad = [i for i in range(N)
               if abs(fr["cells"][i]["entr"] - d_prof[i])
               > 1e-9 + 1e-9 * d_prof[i]]
    if ent_bad:
        fails.append(f"steady entrained cells {ent_bad[:4]} != profile")
    else:
        print(f"ok:   steady droplet profile d_0={d_prof[0]:.10f} -> "
              f"d_{N - 1}={d_prof[-1]:.10f} "
              f"(E_eq*={e_eq_star:.10f} at vG*={vg_star:.10f}, 1e-9)")
    wat_bad = [i for i in range(N)
               if abs(fr["cells"][i]["water"] - PHI_IN) > 1e-9]
    if wat_bad:
        fails.append(f"steady water cells {wat_bad[:4]} != injected cut")
    else:
        print(f"ok:   water cut flushed to the injected {PHI_IN} (1e-9)")
    drift = []
    for key in ("alpha", "water", "entr"):
        for i in range(N):
            if abs(frames[-1]["cells"][i][key]
                   - frames[-2]["cells"][i][key]) > 1e-9:
                drift.append(f"{key}[{i}]")
    if drift:
        fails.append(f"state still drifting between the last two frames: "
                     f"{drift[:4]}")
    else:
        print("ok:   steady state held over the final 10 s")

elif mode == "fullmomentum":
    # full gas momentum closure (distributed pressure momentum): rate
    # driven vertical riser, EOS on. The per-section slip solves
    #   ds/dt = [(rho_m - rhoG)(g sinT + aM) - sigma/alpha_g]/(alpha rhoG)
    # with sigma = C s|s| semi-implicit through the origin secant: the
    # mixture-row pressure gradient is substituted (its alpha_g fr share
    # cancels the proportional gas wall split), so buoyancy and the
    # mixture acceleration drive the slip and wall friction rides the
    # pressure field. Rate-driven: vm is pinned at the injected rates from
    # the constructor on, so aM = 0 and the terminal slip is
    # s* = sqrt(alpha_g alpha (rhoL - rhoG) g / C) with the coupled
    # fixed point F(alpha*, s*) = vsl, F(a) = s a^2 + (vM - s) a
    N, L, DI, DT = 20, 20.0, 0.02, 0.05
    C_DRAG, VSL, VSG, VM = 1000.0, 0.15, 0.25, 0.4
    GAS_R, T0 = 1000.0, 300.0
    RHO_L, RHO_G, MU_L, MU_G = 1000.0, 500.0, 0.01, 1e-5
    SIN_T, PHI_IN, PHI0, A0 = 1.0, 0.3, 0.8, 0.5
    P_OUT = 1.5e8
    DX = L / N

    def rho_g_cell(mg, a):
        ag = 1.0 - a
        return mg / (ag * DX) if ag > 1e-12 else RHO_G

    def slip_map(a, s_old, mg, dt):
        # kernel slipStep model full: buoyancy drive (rate-driven aM = 0),
        # section-local gas density, origin-secant semi-implicit drag
        ag = 1.0 - a
        if a <= 1e-12 or ag <= 1e-12:
            return 0.0
        rho_g = rho_g_cell(mg, a)
        rho_m = a * RHO_L + ag * rho_g
        drive = (rho_m - rho_g) * G * SIN_T
        rate = C_DRAG * abs(s_old) / (ag * a * rho_g)
        return (s_old / dt + drive / (a * rho_g)) / (1.0 / dt + rate)

    def face_fluxes(alpha, slip):
        # liquid and gas volume fluxes of the full closure: the same face
        # machinery as the dynamic closure (interior godunov on the
        # averaged slip, outlet face the last cell's flux)
        flux, gvol = [VSL], [VSG]
        for j in range(1, N):
            s = 0.5 * (slip[j - 1] + slip[j])
            flux.append(godunov(s, VM - s, 0.0, alpha[j - 1], alpha[j]))
            gvol.append(godunov(-s, s - VM, VM, alpha[j - 1], alpha[j]))
        sl = slip[N - 1]
        a_end = alpha[N - 1]
        flux.append(sl * a_end * a_end + (VM - sl) * a_end)
        gvol.append(VM + (sl - VM) * a_end - sl * a_end * a_end)
        return flux, gvol

    # coupled steady state: s*(a) = sqrt(alpha_g a (rhoL - rhoG) g / C)
    # and the holdup root F(alpha*, s*(alpha*)) = vsl, iterated to machine
    # precision (the map contracts at ~0.3 per sweep)
    k_drive = (RHO_L - RHO_G) * G * SIN_T / C_DRAG

    def s_of(a):
        return math.sqrt((1.0 - a) * a * k_drive)

    a_it = 0.7
    for _ in range(200):
        s_it = s_of(a_it)
        disc = (VM - s_it) ** 2 + 4.0 * s_it * VSL
        a_new = ((s_it - VM) + math.sqrt(disc)) / (2.0 * s_it)
        if abs(a_new - a_it) < 1e-15:
            a_it = a_new
            break
        a_it = a_new
    a_star = a_it
    s_star_v = s_of(a_star)

    # coupled transient replication with the gas-mass state (kg/m2)
    alpha = [A0] * N
    slip = [0.0] * N
    water = [PHI0] * N
    mg = [RHO_G * (1.0 - A0) * DX] * N
    steps_done = 0
    net_in = 0.0
    for k, fr in enumerate(frames):
        want_steps = 0 if k == 0 else k * 400
        while steps_done < want_steps:
            slip = [slip_map(alpha[i], slip[i], mg[i], DT) for i in range(N)]
            flux, gvol = face_fluxes(alpha, slip)
            wf = [PHI_IN * flux[0]]
            for j in range(1, N):
                if flux[j] == 0:
                    wf.append(0.0)
                else:
                    wf.append((water[j - 1] if flux[j] >= 0 else water[j])
                              * flux[j])
            wf.append(water[N - 1] * flux[N])
            gm = []
            for j in range(N + 1):
                if j == 0:
                    rho_up = RHO_G  # injected reference density
                elif j == N:
                    rho_up = (rho_g_cell(mg[N - 1], alpha[N - 1])
                              if gvol[j] >= 0 else P_OUT / (GAS_R * T0))
                elif gvol[j] >= 0:
                    rho_up = rho_g_cell(mg[j - 1], alpha[j - 1])
                else:
                    rho_up = rho_g_cell(mg[j], alpha[j])
                gm.append(rho_up * gvol[j])
            new_a, new_w, new_mg = [], [], []
            for i in range(N):
                na = min(max(alpha[i] + DT / DX * (flux[i] - flux[i + 1]),
                             0.0), 1.0)
                af_new = (alpha[i] * water[i]
                          + DT / DX * (wf[i] - wf[i + 1]))
                cap = max(na, 0.0)
                af_new = min(max(af_new, 0.0), cap)
                new_a.append(na)
                new_w.append(af_new / cap if cap > 1e-9 else 0.0)
                new_mg.append(max(mg[i] - DT * (gm[i + 1] - gm[i]), 0.0))
            alpha, water, mg = new_a, new_w, new_mg
            net_in += DT * (gm[0] - gm[N])
            steps_done += 1
        t = fr["t"]
        check(f"vm imposed (t={t})", fr["vm"], VM, 1e-12)
        ok_frame = True
        for i in range(N):
            got = fr["cells"][i]
            want_p = rho_g_cell(mg[i], alpha[i]) * GAS_R * T0
            if (abs(got["alpha"] - alpha[i]) > 1e-9 + 1e-9 * abs(alpha[i])
                    or abs(got.get("s", 0.0) - slip[i])
                    > 1e-9 + 1e-9 * abs(slip[i])
                    or abs(got.get("p", 0.0) - want_p) > 1e-9 * want_p
                    or abs(got.get("water", 0.0) - water[i])
                    > 1e-9 + 1e-9 * abs(water[i])):
                fails.append(f"cell {i} (t={t}): got={got} want alpha="
                             f"{alpha[i]} s={slip[i]} p={want_p} "
                             f"water={water[i]}")
                ok_frame = False
                break
        if ok_frame:
            print(f"ok:   alpha/slip/p/water fields replicated at t={t}"
                  " (1e-9)")
    # global gas inventory: emitted p/alpha reconstruct the solver's
    # inventory, which must equal the replicated boundary-flux balance
    m0 = sum(RHO_G * (1.0 - A0) * DX for _ in range(N))
    m_end = sum(c["p"] * (1.0 - c["alpha"]) * DX / (GAS_R * T0)
                for c in frames[-1]["cells"])
    check("gas inventory balance", m_end, m0 + net_in, 1e-9)
    # steady state: the coupled (alpha*, s*) fixed point with rhoG* = RHO_G
    rho_m = a_star * RHO_L + (1.0 - a_star) * RHO_G
    mu_m = a_star * MU_L + (1.0 - a_star) * MU_G
    re = rho_m * VM * DI / mu_m
    fr_cell = churchill(re, 0.0) * (DX / DI) * 0.5 * rho_m * VM * VM
    grav_cell = rho_m * G * SIN_T * DX
    dp_want = N * (fr_cell + grav_cell)
    p_star = RHO_G * GAS_R * T0
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
        if abs(c["p"] - p_star) > 1e-9 * p_star:
            fails.append(f"steady p {c['p']} != {p_star} (rhoG* != rhoG)")
            break
    else:
        print(f"ok:   uniform steady (alpha*, s*) = ({a_star:.10f},"
              f" {s_star_v:.10f}) with rhoG* = rhoG and the EOS pressure"
              f" {p_star:.6e} (1e-9)")
    check("pressure drop (Churchill friction + gravity)", frames[-1]["dp"],
          dp_want, 1e-9)
    # node pressure march from the outlet anchor: p[j] = p[j+1] + loss_j
    node_bad = [j for j in range(N + 1)
                if abs(frames[-1]["nodes"][j]
                       - (P_OUT + (N - j) * (fr_cell + grav_cell)))
                > 1e-9 * (P_OUT + (N - j) * (fr_cell + grav_cell))]
    if node_bad:
        fails.append(f"pnode march off at nodes {node_bad[:4]}")
    else:
        print(f"ok:   pnode marches from the outlet anchor {P_OUT:.6e} by"
              f" the per-cell loss {fr_cell + grav_cell:.6e} (1e-9)")
    wat_bad = [i for i in range(N)
               if abs(frames[-1]["cells"][i]["water"] - PHI_IN) > 1e-9]
    if wat_bad:
        fails.append(f"steady water cells {wat_bad[:4]} != injected cut")
    else:
        print(f"ok:   water cut flushed to the injected {PHI_IN} (1e-9)")
    drift = []
    for key, tol in (("alpha", 1e-9), ("water", 1e-9), ("s", 1e-9),
                     ("p", 1e-9 * P_OUT)):
        for i in range(N):
            if (abs(frames[-1]["cells"][i][key]
                    - frames[-2]["cells"][i][key]) > tol):
                drift.append(f"{key}[{i}]")
    if drift:
        fails.append(f"state still drifting between the last two frames: "
                     f"{drift[:4]}")
    else:
        print("ok:   steady state held over the final 20 s")

else:
    raise SystemExit(f"unknown mode {mode}")

if fails:
    raise SystemExit(f"two-fluid {mode} verification FAILED: {fails}")
print(f"two-fluid {mode} case matches the analytic reference")
