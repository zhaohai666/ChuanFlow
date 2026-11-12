#!/usr/bin/env python3
r"""transientGelCreep (A2 / 5.74): gel strength and sub-threshold creep enter
the transient row.

  tau_static(t_rest) = tau0 + tau_g (1 - exp(-t_rest / t_g))   (clock at rest)
  held:  |drive| <= 4 tau_static(t_rest) L / d
  creep: q = pi d^4 drive / (128 mu_eff L),  mu_eff = mu + tau0/gamma_dot_c

The case is five independent source-fed chambers (V = 0.1, B = 1e8) draining
to 0 Pa pools through yielded pipes: each chamber ramps at exactly
(B/V) src dt = 102.5 Pa per step while held, so the opening LAYER of every arm
is an integer the checker can predict from the same latch recursion the kernel
runs, the creep arms answer a closed form at every one of the 600 steps, and
the opened arms park on the Buckingham wall law. The zero-declaration branch
is the pre-existing 5.49 row bit for bit -- the whole standing suite is the
witness, no case here carries it.

usage: check_transient_gel_creep.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

RHO, MU, D, L = 850.0, 0.1, 0.04, 20.0
TAU_PLAIN, TAU_HELD = 3.0, 40.0
GEL_S, GEL_T = 2.0, 50.0
CREEP_RATE = 0.5
V_CH, BULK, DT, T_END = 0.1, 1.0e8, 0.5, 300.0
SRC, SRC_HI = 2.05e-7, 2.5625e-7
STEP_PA = BULK / V_CH * SRC * DT              # 102.5 Pa per held step
G_CREEP = math.pi * D ** 4 / (128.0 * (MU + TAU_HELD / CREEP_RATE) * L)
G_PLAIN = math.pi * D ** 4 / (128.0 * MU * L)


def thr(t_rest):
    """The augmented static threshold: 4 (tau0 + tau_g(t_rest)) L / d."""
    tau_g = GEL_S * (1.0 - math.exp(-t_rest / GEL_T))
    return 4.0 * (TAU_PLAIN + tau_g) * L / D


def buckingham_dp(q):
    """dp of the 5.46 Buckingham wall law at rate q: solve
    F(tw) = tw - (4/3) tau0 + tau0^4 / (3 tw^3) = 32 mu q / (pi d^3)
    on the analytic bracket [tau0, ...] by bisection, dp = 4 L tw / d."""
    big = 32.0 * MU * q / (math.pi * D ** 3)

    def f(tw):
        return tw - 4.0 / 3.0 * TAU_PLAIN + TAU_PLAIN ** 4 / (3.0 * tw ** 3)

    lo, hi = TAU_PLAIN, TAU_PLAIN + big + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < big:
            lo = mid
        else:
            hi = mid
    return 4.0 * L * 0.5 * (lo + hi) / D


def parse(stdout):
    samples = []
    cur = None
    for ln in stdout.splitlines():
        m = re.match(r"\[chflow\] t (\S+)", ln)
        if m:
            cur = {"t": float(m.group(1)), "p": {}, "q": {}, "hold": {}}
            samples.append(cur)
            continue
        m = re.match(r"\[chflow\] node (\S+) p (\S+)", ln)
        if m and cur:
            cur["p"][m.group(1)] = float(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) Q (\S+)", ln)
        if m and cur:
            cur["q"][m.group(1)] = float(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) static hold (\S+)", ln)
        if m and cur:
            cur["hold"][m.group(1)] = int(m.group(2))
    return samples


def main():
    binary, case = sys.argv[1], sys.argv[2]
    r = subprocess.run([binary, case], capture_output=True, text=True,
                       cwd=os.path.dirname(case))
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit("transient gel/creep solver failed")
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    head = [ln for ln in r.stdout.splitlines()
            if "nodes" in ln and "pipes" in ln]
    require(bool(head) and "10 nodes, 5 pipes" in head[0],
            f"the case is the shipped 10-node / 5-pipe march ({head})")
    samples = parse(r.stdout)
    n_exp = int(round(T_END / DT)) + 1
    require(len(samples) == n_exp
            and abs(samples[-1]["t"] - T_END) < 1e-12,
            f"the trace is {n_exp} reported samples to t = {T_END:g} "
            f"(got {len(samples)})")

    def opening_layer(name):
        for i, s in enumerate(samples):
            if s["q"][name] != 0.0:
                return i          # sample i is the end of step i (t = i dt)
        return None

    # --- the opening layers are the predicted integers -------------------
    # While held the chamber ramps STEP_PA per step, so drive at step n is
    # STEP_PA (n-1); the gel clock has walked 0.5 (n-1) s of rest. The
    # predicted layer is the first n whose latched drive leaves
    # holdsStatic; the kernel must land exactly there, not one step near.
    def predicted_layer(gel):
        n = 1
        while n * DT <= T_END:
            t_rest = 0.0 if n == 1 else (n - 1) * DT
            limit = thr(t_rest) if gel else 4.0 * TAU_PLAIN * L / D
            if STEP_PA * (n - 1) > limit:
                return n
            n += 1
        return None

    for name, gel, chamber in (("Plain", False, "F1"), ("Gel", True, "F2")):
        pred = predicted_layer(gel)
        got = opening_layer(name)
        require(got == pred,
                f"{name}: opening layer = {got} (t = {got * DT:g} s), "
                f"exactly the latch recursion's prediction {pred} "
                f"(ramp {STEP_PA:g} Pa/step vs threshold "
                f"{'6000 + 4000 (1 - e^{-t_rest/50})' if gel else '6000'})")
        prev = samples[got - 1]
        require(prev["q"][name] == 0.0 and prev["hold"][name] == 1
                and abs(prev["p"][chamber] - STEP_PA * (got - 1))
                    <= 1e-9 * STEP_PA * (got - 1),
                f"{name}: the last held sample is pure ramp -- "
                f"p[{chamber}] = {prev['p'][chamber]:.4f} = "
                f"{STEP_PA:g} x {got - 1}, hold flag set, Q = 0")
    delay = opening_layer("Gel") - opening_layer("Plain")
    require(delay == 22,
            f"the gel arm opens {delay} layers ({delay * DT:g} s) after the "
            f"plain one: the 22-step gap IS the restart pressure the gel "
            f"adds, observable as calendar, not as a fitted number")
    n_open = opening_layer("Gel")
    t_rest_open = (n_open - 1) * DT
    drive_open = STEP_PA * (n_open - 1)
    require(drive_open > thr(t_rest_open)
            and STEP_PA * (n_open - 2) <= thr(t_rest_open - DT),
            f"the gel opening drive {drive_open:.1f} Pa sits against the "
            f"augmented threshold {thr(t_rest_open):.1f} Pa at t_rest = "
            f"{t_rest_open:g} s (the previous step's {STEP_PA * (n_open - 2):.1f} "
            f"still under {thr(t_rest_open - DT):.1f}): the wall held "
            f"{thr(t_rest_open) - 6000:.1f} Pa of rebuild, the closed form "
            f"4000 (1 - e^{{-t_rest/50}})")

    # --- the never-opening wall (5.49 witness inside the same march) -----
    # (the t = 0 sample prints before any latch has run: its hold flag is
    # the constructor's 0, so the flag is read from the solved samples on)
    bad = [s["t"] for s in samples if s["q"]["Held"] != 0.0]
    bad += [s["t"] for s in samples[1:] if s["hold"]["Held"] != 1]
    require(not bad,
            "Held: Q = 0 and the hold flag set at every one of the "
            f"{len(samples)} samples -- thr = 80000 Pa against a chamber "
            "that only ever reaches 61500")
    p3 = samples[-1]["p"]["F3"]
    require(abs(p3 - STEP_PA * 600) <= 1e-9 * STEP_PA * 600,
            f"Held's chamber ends at {p3:.2f} Pa = 102.5 x 600, the pure "
            "ramp with zero outflow, and 61500 < 80000 is why it never "
            "opened")

    # --- the creep arms: closed form at every step ------------------------
    worst = 0.0
    worst_hi = 0.0
    for i in range(1, len(samples)):
        q, q_hi = samples[i]["q"]["Creep"], samples[i]["q"]["CreepHi"]
        # the creep row is a forced value built on the LATCHED drive: the
        # pressure the step started from, i.e. the previous sample's report
        worst = max(worst, abs(q - G_CREEP * samples[i - 1]["p"]["F4"])
                    / max(q, 1e-300))
        worst_hi = max(worst_hi,
                       abs(q_hi - G_CREEP * samples[i - 1]["p"]["F5"])
                       / max(q_hi, 1e-300))
    require(worst <= 1e-9 and worst_hi <= 1e-9,
            f"Creep and CreepHi: q = pi d^4 p_latched / (128 mu_eff L) at "
            f"every step, mu_eff = {MU} + {TAU_HELD}/{CREEP_RATE} = "
            f"{MU + TAU_HELD / CREEP_RATE:g} -- worst relative error "
            f"{max(worst, worst_hi):.2e} over {len(samples) - 1} steps x 2 "
            "arms (the bi-viscosity branch IS that line)")
    ratios = [samples[i]["q"]["CreepHi"] / samples[i]["q"]["Creep"]
              for i in range(1, len(samples))
              if samples[i]["q"]["Creep"] != 0.0]
    require(all(abs(x - 1.25) <= 1e-9 * 1.25 for x in ratios),
            f"CreepHi / Creep = 1.25 at every step (worst deviation "
            f"{max(abs(x - 1.25) for x in ratios):.2e}): the creep branch "
            "is linear in the drive, so a 1.25x source is a 1.25x rate")
    p4, p5 = samples[-1]["p"]["F4"], samples[-1]["p"]["F5"]
    require(abs(p4 - SRC / G_CREEP) <= 1e-4 * SRC / G_CREEP
            and abs(p5 - SRC_HI / G_CREEP) <= 1e-4 * SRC_HI / G_CREEP,
            f"the creep chambers park on src/g = {SRC / G_CREEP:.1f} / "
            f"{SRC_HI / G_CREEP:.1f} Pa (got {p4:.1f} / {p5:.1f}): the "
            "creep rate carries exactly the injection at equilibrium")

    # --- the opened arms park on the Buckingham wall law ------------------
    dp_star = buckingham_dp(SRC)
    late = [s for s in samples if s["t"] >= 250.0]
    bad = []
    for s in late:
        for name, ch in (("Plain", "F1"), ("Gel", "F2")):
            if abs(s["q"][name] - SRC) > 1e-9 * SRC \
                    or abs(s["p"][ch] - dp_star) > 1e-9 * dp_star:
                bad.append((s["t"], name))
    require(not bad,
            f"Plain and Gel at t >= 250: Q = src = {SRC:g} and "
            f"p_up = dp* = {dp_star:.4f} Pa to 1e-9, the 5.46 Buckingham "
            f"root tw* = {dp_star * D / (4.0 * L):.7f} Pa -- the opened "
            "wall is the steady law, not a new one")
    require(dp_star > 4.0 * TAU_PLAIN * L / D,
            f"dp* = {dp_star:.1f} Pa exceeds the static threshold "
            f"{4.0 * TAU_PLAIN * L / D:g} Pa from above: the moving-wall "
            "law tends to its threshold, it never dips under it")

    # --- the rejections: loud, not silent defaults (probes mutate) --------
    text = open(case).read()
    probes = [
        ("negative gelStrength",
         text.replace("gelStrength 2;", "gelStrength -2;", 1),
         "negative"),
        ("gelStrength without a positive gelTime",
         text.replace("gelTime 50;", "gelTime 0;", 1),
         "positive gelTime"),
        ("gelStrength without a yieldStress",
         text.replace("gelStrength 2; gelTime 50;",
                      "gelStrength 2; gelTime 50; yieldStress 0;", 1),
         "no static wall"),
        ("negative creepRate",
         text.replace("creepRate 0.5;", "creepRate -0.5;", 1),
         "negative"),
        ("creepRate without a yieldStress",
         text.replace("creepRate 0.5; initialFlow 0;",
                      "creepRate 0.5; yieldStress 0; initialFlow 0;", 1),
         "no held branch"),
    ]
    for what, txt, needle in probes:
        p = os.path.join(os.path.dirname(case), "probe.dict")
        open(p, "w").write(txt)
        rr = subprocess.run([binary, p], capture_output=True, text=True,
                            cwd=os.path.dirname(case))
        os.remove(p)
        require(rr.returncode != 0 and needle in rr.stdout + rr.stderr,
                f"rejecting {what}: non-zero exit naming the reason "
                f"('{needle}')")

    print()
    if fails:
        print(f"FAILED: {len(fails)}")
        for f in fails:
            print("  -", f)
        raise SystemExit(1)
    print("OK: gel strength and sub-threshold creep enter the transient "
          "row as their closed forms -- the gel arm opens 22 layers late "
          "against the augmented threshold 4 (tau0 + tau_g (1 - "
          "e^{-t_rest/t_g})) L/d, the creep arms answer "
          "pi d^4 drive / (128 mu_eff L) at every step (and stay linear "
          "at 1.25x), the never-opening wall stays the 5.49 Q = 0, the "
          "opened walls park on the 5.46 Buckingham law, and the five "
          "declaration scope violations are refused loudly.")


if __name__ == "__main__":
    main()
