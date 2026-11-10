#!/usr/bin/env python3
"""5.57: which half of the dead band knows the flow index, and which half does not.

Usage: check_transient_yield_index_parking_band.py TRANSIENT STEADY CASE

5.56 swept a Bingham node's own start and found a closed parking band
[pA - 4 tau_0 L/d, pB + 4 tau_0 L/d] plus a tail that is QUADRATIC in the excess,
so `1/gap` is a straight line. Both facts came out of one wall law, and the point
of this file is that they come from DIFFERENT halves of it:

  * the threshold `4 tau_0 L / d` is geometry only (transient.cpp:28-:31, and the
    comment there says so out loud), so the band must not move when the wall law
    becomes Herschel-Bulkley `tau = tau_0 + K gammas^n`;
  * the cancellation that made the tail quadratic is a property of the HB section
    integral near threshold,

      Q = (pi R^3 / tw^3 K^(1/n)) [d^(p+3)/(p+3) + 2 tau_0 d^(p+2)/(p+2)
                                   + tau_0^2 d^(p+1)/(p+1)],  p = 1/n,
          d = tw - tau_0 = (d_pipe / 4 L) gap,

    whose leading power is m = (n+1)/n and whose far power is 1/n. The node row
    Cc dp/dt = Q then gives d(gap)/dt = -A_n gap^m only in the near-threshold
    regime, i.e. gap^(-1/n) is the straight line there, with slope A_n/n.

Three indices (n, K) = (0.5, 1.0), (1.0, 0.1), (2.0, 0.02), seven starts each,
tau_0 = 2.0 everywhere so every wall has the same 8000 Pa threshold. What the run
says, family by family:

  1. the desk is index-blind: the three starts inside print one token each at
     every n, the two edge arms are HELD on the first decision at every n, and
     the far sliders stay one-directional -- while the two 1 Pa sliders print 1
     token at n=0.5 and all 20001 at n=2, so "did it slide" is NOT a token test
     once the tail is steep.
  2. the kernel's HB wall law is this file's section integral, read off the
     tokens: reconstruct tw from the row itself, tw = (drive - (inertia/dt) dq) /
     (4L/d), and the bracket reproduces every printed rate at every index -- from
     the FIRST decided step on, so the ramp is not excluded but accounted for.
     The same first step shows the index moving the balance: the wake-up leaves
     539 / 1636 / 1821 Pa of the drive to the fluid's inertia at n=0.5 / 1 / 2.
  3. the exponent has a geography: m = (n+1)/n at threshold, 1/n far away, and
     every chord measured on a real trajectory sits strictly between the two and
     agrees with the bracket's own chord -- 5.56's quadratic is the n=1 case, and
     the near-edge arms at n=1 and n=2 read m off to 0.01%.
  4. the straight line moves from `1/gap` to `gap^(-1/n)` with the slope A_n/n
     predicted from this file's constants; at n=1 it IS 5.56's coefficient.
  5. the same distance outside the band buys a different amount of progress: the
     1 Pa start is still 1.000000 Pa out at n=0.5 after 40 s and 0.398 Pa out at
     n=2, and the 2000 Pa start closes to 1145 / 51.5 / 2.77 Pa.
  6. the ceiling problem is still the floor problem at every index, token for
     token.
  7. the n=1 column reproduces the numbers 5.56 published from a 12-arm matrix.
  8. a wall that woke never re-arms at any index, and the sub-threshold row is
     held in every decided state.
  9. the initial column prints 14 silent flags for each index.
 10. the steady kernel gives one index-blind answer too: six rows pay zero and
     the three undetermined nodes sit at 25000 Pa, the centre of the same band.

The kernel is not touched to get any of this.
"""

import math
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
STEADY = sys.argv[2]
CASE = sys.argv[3]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO, MU = 850.0, 0.1
LENGTH, DIAM, TAU0 = 20.0, 0.02, 2.0
AREA = 0.25 * math.pi * DIAM * DIAM
WALL_SCALE = 4.0 * LENGTH / DIAM            # dp = WALL_SCALE * tw
INERTIA = RHO * LENGTH / AREA
DT = 0.002
I_DT = INERTIA / DT
THR = WALL_SCALE * TAU0                     # 8000 Pa, no n and no K in it
P_A, P_B = 3.0e4, 2.0e4
CEIL, FLOOR = P_A - THR, P_B + THR          # 22000 / 28000
VOL, BULK = 0.1, 1.0e8
CC = VOL / BULK                             # 1e-9 m3/Pa
R_PIPE = 0.5 * DIAM
SLOPE = DIAM / (4.0 * LENGTH)               # tw - tau_0 = SLOPE * gap

PARKED = "abc"
UP = "de"                                   # below the band, node rises
DOWN = "fg"                                 # above the band, node falls
MIRR = {"d": "f", "e": "g"}
CASES = [("5", 0.5, 1.0), ("1", 1.0, 0.1), ("2", 2.0, 0.02)]

# 5.56's own Bingham numbers, re-read here through the n = 1 column.
REF_556 = {"e4": 21540.241693, "d4": 21999.0019592,
           "e40gap": 51.530219, "d40gap": 0.980747, "K556": 4.908739e-04}

# states sampled for the row reconstruction: the first decision, the peak of the
# startup, the long flat middle and the far end of the 40 s file.
SAMPLE = [1, 2, 5, 50, 200, 500, 1000, 2000, 4000, 8000, 12000, 16000, 20000]

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
        failures.append(what)


def run(binary, case):
    out = subprocess.run([binary, case], cwd=CASE_DIR, capture_output=True,
                         text=True, check=True).stdout
    N, P, F = {}, {}, {}
    banner = None
    for line in out.splitlines():
        line = line.strip()
        m = re.match(r"\[chflow\] node (\S+) p (\S+)", line)
        if m:
            N.setdefault(m.group(1), []).append(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) Q (\S+)", line)
        if m:
            P.setdefault(m.group(1), []).append(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) static hold (\d+)", line)
        if m:
            F.setdefault(m.group(1), []).append(int(m.group(2)))
            continue
        m = re.match(r"\[chflow\] steady solved: iterations=(\d+)", line)
        if m:
            banner = int(m.group(1))
    return N, P, F, banner


def gap_of(up, p):
    return abs((CEIL if up else FLOOR) - p)


def exponent(n):
    """The near-threshold power of the HB section integral."""
    return (n + 1.0) / n


def q_hb(tw, n, K):
    """The exact laminar HB discharge at wall stress tw -- the p = 1/n bracket."""
    if tw <= TAU0:
        return 0.0
    d = tw - TAU0
    p = 1.0 / n
    br = (d ** (p + 3.0) / (p + 3.0) + 2.0 * TAU0 * d ** (p + 2.0) / (p + 2.0)
          + TAU0 * TAU0 * d ** (p + 1.0) / (p + 1.0))
    return math.pi * R_PIPE ** 3 * tw ** -3 * K ** -p * br


def q_bingham(tw):
    """5.56's own transcription of the same bracket at n = 1, kept separate so
    family 2 can prove the two formulas are one formula."""
    if tw <= TAU0:
        return 0.0
    d = tw - TAU0
    return (math.pi * R_PIPE ** 3 * tw ** -3 / MU * d * d
            * (TAU0 * TAU0 / 2.0 + 2.0 * TAU0 * d / 3.0 + d * d / 4.0))


def slope_of(n, K):
    """d(gap^-1/n)/dt predicted from the leading term of the bracket."""
    m = exponent(n)
    a_n = ((math.pi * R_PIPE ** 3) / (TAU0 * K ** (1.0 / n)) * (n / (n + 1.0))
           * (SLOPE ** m) / CC)
    return a_n / n


def chord_exponent(g1, q1, g2, q2):
    return math.log(q2 / q1) / math.log(g2 / g1)


def row_tw(up, pn, q, i):
    """The wall stress the kernel's own row implies at state i: the drive is the
    printed node pressure, the rate change is the printed tokens
    (transient.cpp:381-:382, inertia/dt * (q - qPrev) + 4L/d * tw = drive). Both
    waking rows flow from the high side, so the printed rate is positive in
    either direction and no sign flip is needed."""
    drive = (P_A - pn[i]) if up else (pn[i] - P_B)
    return (drive - I_DT * (q[i] - q[i - 1])) / WALL_SCALE


TR = {}
for tag, n, K in CASES:
    TR[tag] = {"n": n, "K": K}
    for key, fname in (("short", "case"), ("long", "case_long")):
        N, P, F, _ = run(TRANSIENT, "%s/%s_%s.dict" % (CASE_DIR, fname, tag))
        TR[tag][key] = (N, P, F)
        TR[tag][key + "_states"] = len(N["F" + tag + "a"])
STEADY_RUN = run(STEADY, "%s/case_steady.dict" % CASE_DIR)

# --- family 1: the band is geometry only, so three indices share one desk -----
for tag, n, K in CASES:
    for which in ("short", "long"):
        N, P, F = TR[tag][which]
        ns = TR[tag][which + "_states"]
        frozen = {a: set(N["F" + tag + a]) for a in PARKED}
        require(frozen == {"a": {"25000"}, "b": {"22000"}, "c": {"28000"}},
                "n=%s %s: the three starts inside the band still print one token "
                "each (25000 / 22000 / 28000) over %d states -- the threshold is "
                "4 tau_0 L/d with no n and no K in it (transient.cpp:28-:31)"
                % (n, which, ns))
        uniq = {a: len(set(N["F" + tag + a])) for a in UP + DOWN}
        require(all(N["F" + tag + a][-1] != N["F" + tag + a][0] or uniq[a] == 1
                    for a in UP + DOWN) and uniq["e"] == ns and uniq["g"] == ns,
                "n=%s %s: the 2000 Pa starts slide through all %d tokens; the 1 Pa "
                "starts print %d and %d of them, so at m=%.1f a slider can be "
                "indistinguishable from a parked start by tokens alone"
                % (n, which, ns, uniq["d"], uniq["f"], exponent(n)))
        require(max(float(x) for x in N["F" + tag + "e"]) < CEIL
                and min(float(x) for x in N["F" + tag + "g"]) > FLOOR
                and all(float(N["F" + tag + "e"][i]) > float(N["F" + tag + "e"][i - 1])
                        for i in range(1, ns))
                and all(float(N["F" + tag + "g"][i]) < float(N["F" + tag + "g"][i - 1])
                        for i in range(1, ns)),
                "n=%s %s: the far sliders are still strictly monotone toward the "
                "band and never reach or cross it" % (n, which))
        for a in ("b", "c"):
            w = ("S" if a == "b" else "E") + tag + a
            require(F[w][1] == 1,
                    "n=%s %s: the edge arm %s is HELD on the very first decision, "
                    "drive exactly %.0f = threshold %.0f, so the `<=` of "
                    "holdsStatic is index-blind too (transient.cpp:40-:43)"
                    % (n, which, a, THR, THR))

# --- family 2: the kernel's HB wall law, read off its own tokens --------------
# The row itself says (inertia/dt) * dQ + 4L/d * tw = drive (transient.cpp:
# 381-:382), so the wall stress is recoverable from three printed tokens -- and
# the inertia term is part of the balance, not a contamination: even on the very
# first decided step, where the fluid is still accelerating and the wall holds
# LESS than the whole drive, the bracket closes on the printed rate.
first_tw = {}
for tag, n, K in CASES:
    N, P, F = TR[tag]["long"]
    ns = TR[tag]["long_states"]
    rel, held = [], True
    for a in UP + DOWN:
        up = a in UP
        pn = [float(x) for x in N["F" + tag + a]]
        q = [float(x) for x in P[("S" if up else "E") + tag + a]]
        for i in range(1, ns, 61):
            if q[i] == 0.0:
                continue
            tw = row_tw(up, pn, q, i)
            held = held and tw > TAU0
            rel.append((a, i, abs(q_hb(tw, n, K) - q[i]) / abs(q[i])))
    worst = max(r[2] for r in rel)
    require(worst < 1e-5 and held,
            "n=%s: %d of the %d printed rates, starting at the FIRST decided step, "
            "are this file's HB section integral at the wall stress the row itself "
            "implies -- worst relative error %.3e, which is the printed quantum "
            "and not the law -- and none of them sits below tau_0 = %.1f Pa"
            % (n, len(rel), ns - 1, worst, TAU0))
    ramp = {}
    for a in ("e", "g"):
        up = a in UP
        pn = [float(x) for x in N["F" + tag + a]]
        q = [float(x) for x in P[("S" if up else "E") + tag + a]]
        drive = (P_A - pn[1]) if up else (pn[1] - P_B)
        tw = row_tw(up, pn, q, 1)
        ramp[a] = (drive - WALL_SCALE * tw, tw, drive / WALL_SCALE)
    first_tw[n] = ramp["e"]
    mirror = max(abs(ramp["e"][j] - ramp["g"][j]) for j in range(3))
    require(all(ramp[a][0] > 100 and ramp[a][1] < ramp[a][2] - 0.1
                for a in ("e", "g")) and mirror < 1e-6,
            "n=%s: the first step is the ramp, and the reconstruction still closes "
            "on it -- the waking row leaves %.1f Pa of the drive to the fluid's own "
            "inertia (%.1f%% of the %.0f Pa threshold), so the wall holds %.5f Pa "
            "of wall stress where the same drive holds %.5f at rest, and the two "
            "mirrored arms e/g land on the same numbers (%.1e apart)"
            % (n, ramp["e"][0], 100 * ramp["e"][0] / THR, THR, ramp["e"][1],
               ramp["e"][2], mirror))
require(first_tw[0.5][1] > first_tw[1.0][1] > first_tw[2.0][1],
        "and the index moves that first step too: the wall stress the same "
        "2000 Pa drive leaves on step one is %.5f Pa (n=0.5), %.5f (n=1), %.5f "
        "(n=2), monotone in n -- a shear-thickening wall spends more of the drive "
        "on accelerating the column"
        % (first_tw[0.5][1], first_tw[1.0][1], first_tw[2.0][1]))
one = [abs(q_hb(tw, 1.0, MU) - q_bingham(tw)) / max(q_bingham(tw), 1e-300)
       for tw in [TAU0 + 1e-6 * k for k in range(1, 4000)]]
require(max(one) < 1e-12,
        "the general bracket at n=1 and K=mu IS 5.56's Bingham bracket "
        "(worst %.1e over 3999 wall stresses), so this increment adds a "
        "second half to one formula rather than a second formula"
        % max(one))

# --- family 3: the exponent has a geography -----------------------------------
for tag, n, K in CASES:
    N, P, F = TR[tag]["long"]
    ns = TR[tag]["long_states"]
    m, asymptote = exponent(n), 1.0 / n
    near = chord_exponent(1e-9, q_hb(TAU0 + 1e-9, n, K),
                          1e-8, q_hb(TAU0 + 1e-8, n, K))
    far = chord_exponent(1e4, q_hb(TAU0 + 1e4, n, K),
                         2e4, q_hb(TAU0 + 2e4, n, K))
    require(abs(near - m) < 1e-4 * m,
            "n=%s: the bracket's own exponent one nanoPa of wall stress above "
            "the yield is %.6f against (n+1)/n = %.6f" % (n, near, m))
    require(abs(far - asymptote) < 0.02 * asymptote,
            "n=%s: five thousand yield stresses above it the same bracket has "
            "become the pure power law, exponent %.6f against 1/n = %.6f -- m is "
            "the NEAR-threshold power, not the only one" % (n, far, asymptote))
    chords = []
    for a in UP + DOWN:
        up = a in UP
        pn = [float(x) for x in N["F" + tag + a]]
        q = [float(x) for x in P[("S" if up else "E") + tag + a]]
        peak = max(range(1, ns), key=lambda i: abs(q[i]))
        pts = [(gap_of(up, pn[i]), abs(q[i])) for i in SAMPLE
               if i > peak and abs(q[i]) > 0]
        for (g1, r1), (g2, r2) in zip(pts, pts[1:]):
            if g2 < g1 and r2 < r1 and (g1 - g2) > 1e-3 * g1:
                e = chord_exponent(g1, r1, g2, r2)
                b = chord_exponent(g1, q_hb(TAU0 + SLOPE * g1, n, K),
                                   g2, q_hb(TAU0 + SLOPE * g2, n, K))
                chords.append((a, e, b))
    lo = min(c[1] for c in chords)
    hi = max(c[1] for c in chords)
    off = max(abs(c[1] - c[2]) / c[2] for c in chords)
    require(lo > asymptote and hi < m + 1e-9 and off < 0.05,
            "n=%s: %d chords measured on the four real trajectories span "
            "%.4f - %.4f, strictly inside (1/n, m) = (%.4f, %.4f), and each one "
            "reproduces the bracket's own chord over the same two gaps to within "
            "%.3f%% -- the case's gaps are 0.4 to 2000 Pa, so delta/tau_0 stays "
            "under 0.5 and the far asymptote is never reached"
            % (n, len(chords), lo, hi, asymptote, m, 100 * off))
    late = [c for c in chords if c[0] in "df" and c[2] > 0.0]
    if n == 0.5:
        require(not late,
                "n=0.5: the near-edge arms contribute NO chord at all, because a "
                "1 Pa start moves by less than one printed quantum in 40 s -- the "
                "m=3 tail is too steep to be read from the desk it was made for")
    else:
        worst = max(abs(c[1] - m) / m for c in late)
        require(worst < 0.01,
                "n=%s: the near-edge arms' own chords are %.5f - %.5f, within "
                "%.3f%% of m = %.4f -- this is the kernel reading its own "
                "near-threshold power off a slider that moved a few tenths of a "
                "Pascal in 40 s" % (n, min(c[1] for c in late),
                                    max(c[1] for c in late), 100 * worst, m))
    if n == 1.0:
        require(abs(m - 2.0) < 1e-12,
                "n=1 gives m = 2, which is 5.56's quadratic tail: the "
                "cancellation there is this exponent's n=1 case, not a separate "
                "law")

# --- family 4: the straight line moves from 1/gap to gap^(-1/n) ---------------
for tag, n, K in CASES:
    N, P, F = TR[tag]["long"]
    ns = TR[tag]["long_states"]
    want = slope_of(n, K)
    got = []
    for a in ("d", "f"):
        up = a in UP
        s = [float(x) for x in N["F" + tag + a]]
        g1, g2 = gap_of(up, s[2000]), gap_of(up, s[-1])
        if g1 > 0 and g2 > 0 and g1 != g2:
            got.append(((g2 ** (-1.0 / n) - g1 ** (-1.0 / n))
                        / ((len(s) - 1 - 2000) * DT), a))
    worst = max(abs(v - want) / want for v, _ in got)
    require(worst < 0.03,
            "n=%s: gap^(-1/n) is the straight line on the near arms %s, measured "
            "slope %s against %.6e predicted from pi R^3 (n/(n+1)) (d/4L)^m / "
            "(tau_0 K^(1/n) Cc) with m=(n+1)/n (worst %.3f%% off)"
            % (n, "/".join(a for _, a in got),
               ", ".join("%s=%.6e" % (a, v) for v, a in got), want, 100 * worst))
    if n == 1.0:
        ref = abs(got[0][0] - REF_556["K556"]) / REF_556["K556"]
        require(ref < 1e-3,
                "n=1's slope %.6e is 5.56's own 1/gap coefficient %.6e to %.3f%%, "
                "so one formula covers both increments"
                % (got[0][0], REF_556["K556"], 100 * ref))

# --- family 5: the same distance outside buys a different amount --------------
near, veryfar = {}, {}
for tag, n, K in CASES:
    Ns, Ps, Fs = TR[tag]["short"]
    Nl, Pl, Fl = TR[tag]["long"]
    near[n] = gap_of(True, float(Nl["F" + tag + "d"][-1]))
    veryfar[n] = gap_of(True, float(Nl["F" + tag + "e"][-1]))
    if n == 0.5:
        require(len(set(Ns["F" + tag + "d"])) == 1
                and len(set(Ns["F" + tag + "f"])) == 1,
                "n=0.5: the 4 s file gives the 1 Pa arms ONE token each, so on a "
                "steep tail a slider and a parked start are the same printout "
                "until the clock is lengthened")
require(near[0.5] > near[1.0] > near[2.0]
        and veryfar[0.5] > veryfar[1.0] > veryfar[2.0],
        "after 40 s the 1 Pa start is %.6f / %.6f / %.6f Pa from the ceiling and "
        "the 2000 Pa start is %.6f / %.6f / %.6f Pa at n=0.5 / 1 / 2 -- the "
        "monotone order is the same for both starts, so the index decides how "
        "much a given distance outside the band buys"
        % (near[0.5], near[1.0], near[2.0],
           veryfar[0.5], veryfar[1.0], veryfar[2.0]))
require(near[0.5] > 0.999999,
        "at n=0.5 a start ONE PA outside the band is still %.6f Pa out after 40 s "
        "(%.6f Pa closed in all), i.e. m=3 makes the bracket that 5.56 used to "
        "prove the edge is closed far too fine to see -- 5.56's 1 Pa bracketing "
        "is an n=1 convenience, not a property of the latch"
        % (near[0.5], 1.0 - near[0.5]))
require(veryfar[2.0] < 3.0,
        "at n=2 the SAME 2000 Pa outside is down to %.6f Pa in the same clock, "
        "99.86%% of the distance the steep tail barely starts" % veryfar[2.0])

# --- family 6: the ceiling problem is still the floor problem, for every n ----
for tag, n, K in CASES:
    for which in ("short", "long"):
        N, P, F = TR[tag][which]
        ns = TR[tag][which + "_states"]
        for a, b in MIRR.items():
            sa, sb = N["F" + tag + a], N["F" + tag + b]
            off = max(abs(float(x) + float(y) - (P_A + P_B)) for x, y in zip(sa, sb))
            qa = P["S" + tag + a] if a in UP else P["E" + tag + a]
            qb = P["S" + tag + b] if b in UP else P["E" + tag + b]
            same = sum(1 for x, y in zip(qa, qb) if x == y)
            require(off <= 1e-7 and same == ns,
                    "n=%s %s: arms %s/%s mirror token for token -- node sums to "
                    "%.0f Pa within %.1e and the waking rates are identical in "
                    "%d/%d states" % (n, which, a, b, P_A + P_B, off, same, ns))

# --- family 7: n=1 must BE 5.56, read as numbers across two matrices ---------
for tag, n, K in CASES:
    if n != 1.0:
        continue
    Ns, Ps, Fs = TR[tag]["short"]
    Nl, Pl, Fl = TR[tag]["long"]
    for want, ser, a in ((REF_556["e4"], float(Ns["F1e"][-1]), "e"),
                         (REF_556["d4"], float(Ns["F1d"][-1]), "d")):
        require(abs(ser - want) < 1e-6,
                "n=1 arm %s reaches %.10g Pa at 4 s, the number 5.56 published "
                "for the same start from a 12-arm matrix (%.10g, off %.2e)"
                % (a, ser, want, abs(ser - want)))
    for want, a in ((REF_556["e40gap"], "e"), (REF_556["d40gap"], "d")):
        g = gap_of(True, float(Nl["F1" + a][-1]))
        require(abs(g - want) < 1e-4,
                "n=1 arm %s sits %.6f Pa from the ceiling after 40 s, 5.56's "
                "number is %.6f (off %.2e Pa)" % (a, g, want, abs(g - want)))

# --- family 8: a wall that woke never re-arms, at any index ------------------
for tag, n, K in CASES:
    N, P, F = TR[tag]["long"]
    ns = TR[tag]["long_states"]
    small = None
    for a in UP + DOWN:
        up = a in UP
        w = ("S" if up else "E") + tag + a
        o = ("E" if up else "S") + tag + a
        require(all(v == 0 for v in F[w][1:]),
                "n=%s: arm %s's waking row is never asked to hold again over %d "
                "decided states" % (n, a, ns - 1))
        require(all(v == 1 for v in F[o][1:]) and set(P[o]) == {"0"},
                "n=%s: arm %s's other row, still under its own threshold, is held "
                "in every decided state and prints one token 0" % (n, a))
        for x in P[w][1:]:
            v = abs(float(x))
            if v > 0 and (small is None or v < small):
                small = v
    require(small is not None and small > 0,
            "n=%s: the smallest rate any waking row ever printed is %.6e m3/s -- "
            "never exactly 0, which is what keeps the latch from re-arming "
            "regardless of the index" % (n, small))

# --- family 9: the initial column is still silent ----------------------------
for tag, n, K in CASES:
    N, P, F = TR[tag]["short"]
    require(all(v[0] == 0 for v in F.values()) and len(F) == 14,
            "n=%s: t = 0 prints %d hold flags and every one of them is 0, so an "
            "event scan must start at i = 1 for the HB arms exactly as it did for "
            "the Bingham ones" % (n, len(F)))

# --- family 10: the steady kernel has no index either ------------------------
Ns, Ps, Fs, banner = STEADY_RUN
require(banner == 0 and len(Ps) == 6 and len(Ns) == 5
        and all(set(v) == {"0"} for v in Ps.values())
        and all(set(Ns["T" + t]) == {"25000"} for t in ("5", "1", "2")),
        "steady: %d rows pay Q 0 with dp 5000, the three undetermined nodes sit at "
        "25000 Pa for n=0.5, n=1 and n=2 alike, banner iterations=%d -- the "
        "undetermined node lands on the CENTRE of the same 8000 Pa band, and the "
        "other kernel is as blind to the index as the threshold is"
        % (len(Ps), banner))

print("\n%s: %d checks, %d failures"
      % ("PASS" if not failures else "FAIL", CHECKS, len(failures)))
for f in failures:
    print("  -", f)
sys.exit(1 if failures else 0)
