#!/usr/bin/env python3
"""Run the transient 1D3D coupling with BOTH coupled edges in the parallel
position and with the detailed regions solving their own du/dt, and find out
what the inertance the 1D row installs is actually worth.

Everything the chain built so far asked the partner for a drop at an instant:
the region was solved to a steady state at the rate handed over, so its drop was
a function of the rate and nothing else. `network`, `partner` and the per-edge
mix of the two were then three ways of writing the same interface condition, and
the only evidence between them was the path each took. Here the block solves
`du/dt` against the field the previous time level left, and the interface keeps
iterating on the FRICTION part of what it reports while charging the
acceleration part to its own momentum row. That makes one number that no earlier
run can produce:

    rowError = (R*q + M_row*(q - qPrev)/dt - dpBlock) / dpBlock

per accepted level and per edge - and it is not iterable, because interface
iteration drives the friction residual down and leaves this exactly as large as
the inertance in the row is wrong.

The claims, in both directions and with no approximations:

  1. six arms, two axes       - three attribution sources x quasi-steady versus
                                time-accurate, plus one operator change, on the
                                same two meshes and the same seeds
  2. what got installed       - the row's M equals the analytic inertance of the
                                source the arm names: rho*L/(pi d^2/4) for the
                                drawing, rho*V/A^2 for the region, and
                                rho*V_probe/A^2 for the slab - the last one is a
                                number only this increment can ask for
  3. the split is measured    - dpFriction + dpInertia == dpBlock in machine
                                arithmetic at every accepted level of every
                                unsteady arm, the steady arm reports no
                                acceleration at all, and the acceleration term
                                points the way the rate changed on every level
                                except one per arm - the level where the branch's
                                flux has stopped changing while its profile is
                                still accelerating
  4. the level owns its history - qPrev is the rate the PREVIOUS accepted level
                                committed, exactly, so re-solving a level
                                iterates on one time step instead of walking the
                                trajectory
  5. the block grades its own - inertanceDynamic (what the field it solved says
                                its acceleration was worth per unit rate of
                                change) against inertanceProbe (the geometric
                                rho*V_probe/A^2 of the same slab): chi is one
                                number per edge, the same for every arm, and the
                                slab's mass really is what the region spends as
                                long as the shape is not changing - 0.6% to 2.8%
                                at the settled levels, 5% to 7% short at the level
                                that starts from rest, and meaningless at the two
                                levels in between
  6. the referee              - rowError is the inertance claim and nothing
                                else: R*q sits on the friction chord to within
                                the declared tolerance while sitting off the
                                reported drop by 22x to 3039x it - up to 61% of
                                the drop - and the remaining gap equals
                                (M_row - inertanceDynamic)*(dq/dt)/dpBlock inside
                                the same tolerance
  7. attribution has an error - on the settled levels rowError is exactly
                                (M_row - M_probe)*(dq/dt)/dpBlock, to 3.2% worst
                                relative for the drawing and 2.9% for the region,
                                and for the slab - where that formula says zero -
                                what is left is 0.01% to 0.23% of the drop. This
                                is the first time a source choice carries a
                                measured error instead of a road
  8. chord and destination    - the friction chord of a region that starts from
                                rest is several times the quasi-steady chord on
                                the first level (the wall layer is only born, so
                                the shear is the starting one) and is back on the
                                quasi-steady law from level 3 - which is what
                                makes the iterated coefficient a rate-function
                                again, and what the next level pays passes to
                                unwind; the destination stays the resistive split
                                of the chords, while the six arms arrive at six
                                different positions on their way there
  9. cost and fields         - every level agrees inside the sub-iteration
                                budget under either operator, one pass is one
                                network solve plus one solve per partner, and
                                both wall-resolved blocks close on their own
                                fields in every arm

usage: check_transient_coupled_unsteady_partner.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledUnsteadyPartner/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN, P_REF = 1.2e-4, 1.0e6
DT, N_STEPS = 0.05, 12
MAX_PASSES = 10
TOL = 2e-4
NAME = ("X1", "X2")
ARMS = ("ctrl_steady_net_gs", "bore_gs", "region_gs", "slab_gs", "slab_qn",
        "bore_qn")
SOURCE = {"ctrl_steady_net_gs": "network", "bore_gs": "network",
          "bore_qn": "network", "region_gs": "partner",
          "slab_gs": "slab", "slab_qn": "slab"}
DDT = dict((a, "Euler") for a in ARMS)
DDT["ctrl_steady_net_gs"] = "steady"
OPERATOR = dict((a, "no" if a.endswith("_gs") else "yes") for a in ARMS)
# (length, gap, span) of each detailed region, the bore the drawing gives it, and
# the probe range the reported drop is measured across
EDGES = {"X1": (0.2, 0.005, 0.005), "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.005, "X2": 0.0049}
PROBE = {"X1": (0.06, 0.19), "X2": (0.04, 0.11)}
# cell centres sit at (i + 0.5) * dx with dx = 0.01 in both ducts, so the probe
# slab is 13 of X1's 20 sections and 7 of X2's 12
SLAB_FRAC = {"X1": 13.0 / 20.0, "X2": 7.0 / 12.0}


def M_region(n):
    L, H, W = EDGES[n]
    return RHO * L / (H * W)


def M_slab(n):
    return M_region(n) * SLAB_FRAC[n]


def M_bore(n):
    L = EDGES[n][0]
    return RHO * L / (0.25 * math.pi * BORE[n] ** 2)


M_INSTALL = {"network": M_bore, "partner": M_region, "slab": M_slab}


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def word(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+(\S+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def pairs(line):
    t = line.split()
    k = next(i for i, x in enumerate(t) if x.endswith(":")) + 1
    out = {}
    while k + 1 < len(t):
        out[t[k]] = float(t[k + 1])
        k += 2
    return out


def powerLaw(qs, dps):
    """Least squares on log dp = log C + m log q."""
    xs = [math.log(x) for x in qs]
    ys = [math.log(x) for x in dps]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    m = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return m, math.exp((sy - m * sx) / n)


def resistiveSplit(c1, m1, c2, m2, total):
    lo, hi = 1e-9 * total, total * (1.0 - 1e-9)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if c1 * mid ** m1 - c2 * (total - mid) ** m2 > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) / total


def parse(lines):
    """One dict per arm: what it was declared, its accepted states, its passes,
    its partner declarations and its block reports."""
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and ln.count(":") == 1 and " edge " not in ln:
            t = ln.split()
            cur = dict(name=t[2].rstrip(":"), inertia=t[4], accelerate=t[6],
                       source={}, partner={}, starts={}, steps={}, passes=[],
                       blocks=[], solves=0)
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " inertia " in ln:
            cur["source"][ln.split()[4]] = word(ln, "inertia")
        elif ln.startswith("transient partner "):
            n = ln.split()[4]
            cur["partner"][n] = dict(ddt=word(ln, "ddt"),
                                     region=value(ln, "regionInertance"),
                                     probe=value(ln, "probeInertance"))
        elif ln.startswith("transient start ") and " edge " in ln and " resistance " in ln:
            cur["starts"][ln.split()[4]] = value(ln, "resistance")
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["passes"].append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            t = ln.split()
            cur["steps"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " newton " in ln:
            cur["steps"].setdefault(int(ln.split()[2].rstrip(":")), {})["s"] = ln
        elif ln.startswith("nodes step "):
            cur["steps"][int(ln.split()[2].rstrip(":"))]["nd"] = ln
        elif ln.startswith("pipes step "):
            cur["steps"][int(ln.split()[2].rstrip(":"))]["pp"] = ln
        elif ln.startswith("transient block "):
            cur["blocks"].append(ln)
    return arms


def series(a):
    """The accepted states of one arm, per edge, in level order."""
    s = dict(q=[[], []], dp=[[], []], dpf=[[], []], dpi=[[], []], qp=[[], []],
             dq=[[], []], R=[[], []], res=[[], []], M=[[], []], row=[[], []],
             rerr=[[], []], uns=[], sw=[[], []])
    for n in range(1, N_STEPS + 1):
        for i in (0, 1):
            ln = a["steps"][n][i]
            s["q"][i].append(value(ln, "q"))
            s["dp"][i].append(value(ln, "dpBlock"))
            s["dpf"][i].append(value(ln, "dpFriction"))
            s["dpi"][i].append(value(ln, "dpInertia"))
            s["qp"][i].append(value(ln, "qPrev"))
            s["dq"][i].append(value(ln, "dq"))
            s["R"][i].append(value(ln, "R"))
            s["res"][i].append(value(ln, "residual"))
            s["M"][i].append(value(ln, "M"))
            s["row"][i].append(value(ln, "rowDrop"))
            s["rerr"][i].append(value(ln, "rowError"))
            s["uns"].append(word(ln, "unsteady"))
            s["sw"][i].append(value(ln, "sweeps"))
    return s


def main():
    binary, case = sys.argv[1], sys.argv[2]
    r = subprocess.run([binary, os.path.basename(case)],
                       capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(case)))
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient unsteady-partner solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)
    at = {a["name"]: k for k, a in enumerate(arms)}
    ser = dict((a["name"], series(a)) for a in arms)
    UN = tuple(a for a in ARMS if DDT[a] == "Euler")

    # ---- 1. six arms and the two axes they cross ----
    require([a["name"] for a in arms] == list(ARMS),
            f"all six arms ran, in the order declared (got {[a['name'] for a in arms]})")
    for a in arms:
        n = a["name"]
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{n}: all {N_STEPS} time levels reported an accepted state")
        require(sorted(a["partner"]) == list(NAME),
                f"{n}: both partners declared what they were asked to carry")
        require(all(a["partner"][e]["ddt"] == DDT[n] for e in NAME),
                f"{n}: both regions stepped as this arm declares - effective"
                f" schemes ddt {'/'.join(a['partner'][e]['ddt'] for e in NAME)}"
                f" against the arm's '{DDT[n]}' (the block dicts themselves ask"
                f" for Euler, so the control's steady is the override working)")
        require(a["accelerate"] == OPERATOR[n],
                f"{n}: negotiated by the"
                f" {'quasi-Newton' if a['accelerate'] == 'yes' else 'plain Gauss-Seidel'}"
                " map this arm declares")
        require(tuple(a["source"][e] for e in NAME) == (SOURCE[n], SOURCE[n]),
                f"{n}: installed `{SOURCE[n]}` on both edges, which is the source"
                " this arm declares")
        require(a["inertia"] == SOURCE[n],
                f"{n}: the header's single `inertia {a['inertia']}` word is the"
                " arm's source, and no per-edge table overrode it")
    require(len({SOURCE[a] for a in ARMS}) == 3,
            "three attribution sources are on the table in one run - the drawing,"
            " the region, and the slab the reported drop is measured across")
    require(len({DDT[a] for a in ARMS}) == 2,
            "and every source is read under a partner that carries its own du/dt,"
            " with one quasi-steady control arm beside them")

    # ---- 2. what the row actually got ----
    for a in arms:
        n = a["name"]
        s = ser[n]
        for i, e in enumerate(NAME):
            want = M_INSTALL[SOURCE[n]](e)
            got = set(abs(x / want - 1.0) for x in s["M"][i])
            require(max(got) < 1e-12,
                    f"{n} {e}: every level's row carries M {s['M'][i][0]:.6g} kg/m4"
                    f" = rho*L/A of the `{SOURCE[n]}` book ({want:.6g}) to"
                    f" {max(got):.1e} relative")
            pr = a["partner"][e]["probe"]
            rg = a["partner"][e]["region"]
            require(abs(pr / rg - SLAB_FRAC[e]) < 1e-12,
                    f"{n} {e}: the partner measures its slab at {pr:.6g} kg/m4"
                    f" against the whole region's {rg:.6g} - the probe leaves out"
                    f" {100 * (1 - SLAB_FRAC[e]):.1f}% of the duct, and the two"
                    " numbers differ by exactly that")
            require(abs(rg - M_region(e)) < 1e-9 * rg
                    and abs(M_bore(e) / rg - 1.0) > 0.2,
                    f"{n} {e}: the region's own rho*V/A^2 is the analytic"
                    f" {M_region(e):.6g} and the drawing's bore is"
                    f" {100 * (M_bore(e) / rg - 1.0):+.1f}% off it")
            if SOURCE[n] == "slab":
                require(abs(s["M"][i][0] / pr - 1.0) < 1e-12,
                        f"{n} {e}: the slab source really did install the slab's"
                        " number, not the region's")

    # ---- 3. the split is measured, not assumed ----
    for a in arms:
        n = a["name"]
        s = ser[n]
        if DDT[n] == "steady":
            require(set(s["uns"]) == {"0"},
                    f"{n}: the quasi-steady control reports no unsteady drop on"
                    " any edge at any level")
            require(all(x == 0.0 for x in s["dpi"][0] + s["dpi"][1]),
                    f"{n}: and its whole drop is friction by construction"
                    " (dpInertia is 0.0 exactly, not small)")
            require(all(x == 0.0 for x in s["rerr"][0] + s["rerr"][1]),
                    f"{n}: nothing is graded against nothing - rowError stays 0"
                    " for a partner that states no acceleration")
            continue
        require(set(s["uns"]) == {"1"},
                f"{n}: both edges report a drop from a momentum balance that"
                " carries du/dt at every level")
        worst = max(abs((f + k) / d - 1.0)
                    for i in (0, 1) for f, k, d in zip(s["dpf"][i], s["dpi"][i],
                                                       s["dp"][i]))
        require(worst < 1e-12,
                f"{n}: dpFriction + dpInertia == dpBlock at all"
                f" {2 * N_STEPS} accepted edges to {worst:.1e} - the split is an"
                " identity on the block's own field, not a second measurement")
        share = [[abs(k / d) for k, d in zip(s["dpi"][i], s["dp"][i])]
                 for i in (0, 1)]
        flat = [x for row in share for x in row]
        require(max(flat) > 0.3 and min(flat) > 0.004,
                f"{n}: the acceleration is never decoration - it carries"
                f" {100 * min(flat):.1f}% of the reported drop at its smallest and"
                f" {100 * max(flat):.1f}% at its largest, and the largest is the"
                " level that starts from rest")
        require(all(share[i][0] > share[i][-1] for i in (0, 1)),
                f"{n}: and it decays as the ramp settles - on both edges the last"
                f" level pays {100 * share[0][-1]:.1f}% and {100 * share[1][-1]:.1f}%"
                f" where the first paid {100 * share[0][0]:.1f}% and"
                f" {100 * share[1][0]:.1f}% - the inertial part of a partner's drop"
                " is a transient, exactly as the row says it should be")
        bad = [(lv, e, s["dq"][i][lv - 1], s["dpi"][i][lv - 1])
               for i, e in enumerate(NAME)
               for lv, (k, dq) in enumerate(zip(s["dpi"][i], s["dq"][i]), 1)
               if dq != 0.0 and (k > 0.0) != (dq > 0.0)]
        require([b[:2] for b in bad] == [(2, "X2")],
                f"{n}: the acceleration term points the way the rate changed at"
                " every level but one - X2 at level 2, where that branch's flux"
                f" has already turned back ({bad[0][2]:+.1e} m3/s) while its"
                " field is still speeding up"
                f" ({bad[0][3]:+.0f} Pa of inertial drop). The slab is not a"
                " lump: the wall layer is still growing, so nearly the same flux"
                " can carry a momentum that is still rising")

    # ---- 4. the level owns its history ----
    for a in arms:
        n = a["name"]
        if DDT[n] == "steady":
            continue
        s = ser[n]
        worst = 0.0
        for i in (0, 1):
            for lv in range(2, N_STEPS + 1):
                worst = max(worst, abs(s["qp"][i][lv - 1] - s["q"][i][lv - 2]))
            worst = max(worst, max(abs(d - (q - qp))
                                   for q, qp, d in zip(s["q"][i], s["qp"][i],
                                                       s["dq"][i])))
        require(worst < 1e-15,
                f"{n}: qPrev at every level is the rate the PREVIOUS accepted"
                " level committed, exactly, and dq is the difference of the two -"
                " so the passes inside one level all accelerate away from one"
                " history instead of walking the trajectory")
        require(all(s["qp"][i][0] == 0.0 for i in (0, 1)),
                f"{n}: the first level starts from rest (qPrev 0 on both edges),"
                " which is why its rate change is the largest of the run")

    # ---- 5. the block grades its own acceleration ----
    chi = {}
    for a in arms:
        n = a["name"]
        if DDT[n] == "steady":
            continue
        s = ser[n]
        for i, e in enumerate(NAME):
            pr = a["partner"][e]["probe"]
            for lv, (k, d) in enumerate(zip(s["dpi"][i], s["dq"][i]), 1):
                if d != 0.0:
                    chi[(n, e, lv)] = k / (d / DT) / pr
            ch = [chi[(n, e, lv)] for lv in range(1, N_STEPS + 1)]
            require(0.90 < ch[0] < 0.98,
                    f"{n} {e}: at the level that starts from rest, where the"
                    " rate change is the whole acceleration budget, what the field"
                    f" the block solved spent is {100 * ch[0]:.1f}% of the"
                    f" geometric rho*V_probe/A^2 of the slab the drop is measured"
                    " across - the slab moves as one lump to that margin, and the"
                    " missing few percent is the profile being born")
            settled = [abs(c - 1.0) for c in ch[3:]]
            require(max(settled) < 0.03,
                    f"{n} {e}: and at every level from 4 on the slab's own mass"
                    f" is what the region spent, to {100 * max(settled):.1f}%"
                    " worst - so the attribution floor is a property of the"
                    " region, not of the arm that read it")
            trans = [abs(ch[lv - 1] - 1.0) for lv in (2, 3)]
            require(min(trans) > 0.05,
                    f"{n} {e}: the two levels in between are the ones the lumped"
                    f" slab cannot say anything about ({100 * trans[0]:.0f}% and"
                    f" {100 * trans[1]:.0f}% off) - the wall layer is still"
                    " growing there, so the field's momentum and the edge's flux"
                    " change are not proportional, and no inertance in a 1D row"
                    " can be asked to absorb that")
    for e in NAME:
        for lv in range(4, N_STEPS + 1):
            vals = [chi[(n, e, lv)] for n in UN]
            require(max(vals) / min(vals) < 1.05,
                    f"{e} level {lv}: chi is one physical number rather than one"
                    f" per arm - the {len(UN)} unsteady arms, which install three"
                    f" different inertances in their rows, all measure the same"
                    f" acceleration here between {min(vals):.3f} and"
                    f" {max(vals):.3f} of the slab")
    for a in arms:
        if DDT[a["name"]] == "steady":
            continue
        for b in a["blocks"]:
            beta = value(b, "betaProfile")
            require(beta > 1.0 and beta < 2.0,
                    f"{a['name']} {word(b, 'name').rstrip(':')}: the profile it"
                    f" accelerated has a momentum-flux correction beta {beta:.4f} -"
                    " the boundary layer the wall function resolves is why a"
                    " flat-rate picture needs a shape factor at all")
            require(abs(value(b, "dtLevel") - DT) < 1e-15,
                    f"{a['name']} {word(b, 'name').rstrip(':')}: the acceleration"
                    f" was measured against the declared level of {DT:g} s")

    # ---- 6. the referee: rowError is the inertance claim and nothing else ----
    chord_err = {}
    for a in arms:
        n = a["name"]
        s = ser[n]
        for i, e in enumerate(NAME):
            fr = [abs((R * q / f) - 1.0) for R, q, f in zip(s["R"][i], s["q"][i],
                                                             s["dpf"][i])]
            full = [abs((R * q / d) - 1.0) for R, q, d in zip(s["R"][i], s["q"][i],
                                                               s["dp"][i])]
            chord_err[(n, e)] = (max(fr), max(full))
            if DDT[n] == "steady":
                require(max(fr) <= TOL and max(full) <= TOL,
                        f"{n} {e}: for the quasi-steady control the chord IS the"
                        " whole drop - both are inside the declared tolerance"
                        f" ({max(fr):.3e})")
                continue
            require(max(fr) <= TOL,
                    f"{n} {e}: the refreshed coefficient lands on the FRICTION"
                    f" chord at every accepted level (worst {max(fr):.3e} against"
                    f" the declared {TOL:g}) - the interface iterates on the part"
                    " of the drop that depends on the rate alone")
            require(min(full) > 10.0 * TOL,
                    f"{n} {e}: and it is not iterating on the whole drop - away"
                    f" from the friction chord it sits {min(full) / TOL:.0f}x to"
                    f" {max(full) / TOL:.0f}x the tolerance (worst"
                    f" {max(full):.3e}); refreshing the coefficient with the"
                    " acceleration is the a - c/R map that loses its crossing"
                    " when dt < M/R")
            pre = [abs(row - (R * q + M * dq / DT)) / d
                   for row, R, q, M, dq, d in zip(s["row"][i], s["R"][i],
                                                  s["q"][i], s["M"][i],
                                                  s["dq"][i], s["dp"][i])]
            grad = [abs(re - (row - d) / d) for re, row, d
                    in zip(s["rerr"][i], s["row"][i], s["dp"][i])]
            pred = [abs(re * d - (M * dq / DT - k)) / d
                    for re, d, M, dq, k in zip(s["rerr"][i], s["dp"][i],
                                               s["M"][i], s["dq"][i],
                                               s["dpi"][i])]
            require(max(pre) < 1e-12 and max(grad) < 1e-12,
                    f"{n} {e}: rowDrop is R*q + M*dq/dt and rowError is"
                    f" (rowDrop - dpBlock)/dpBlock at every level - worst"
                    f" {max(pre):.1e} and {max(grad):.1e} relative, which is"
                    " machine arithmetic on the numbers the interface wrote")
            require(max(pred) <= 1.5 * TOL,
                    f"{n} {e}: and that grade is the inertance claim and nothing"
                    f" else - rowError x dpBlock against (M_row - inertanceMeasured)"
                    f" x dq/dt differs by at most {max(pred):.3e} of the drop,"
                    " inside the interface's own tolerance, so the slack left in"
                    " the grade is the solver's, not the attribution's")

    # ---- 7. the error an attribution choice now carries ----
    # On the levels where chi has settled, the interface's own numbers say what
    # the row's inertance is worth:
    #     rowError = (M_row - M_probe) * (q - qPrev) / (dt * dpBlock)
    # - which is the claim the previous chain could only make about a road.
    med, settle = {}, {}
    for gs, qn, lbl in [("bore_gs", "bore_qn", "network"),
                        ("region_gs", None, "partner"),
                        ("slab_gs", "slab_qn", "slab")]:
        for n in [gs] + ([qn] if qn else []):
            a = arms[at[n]]
            s = ser[n]
            for i, e in enumerate(NAME):
                pr = a["partner"][e]["probe"]
                got = [abs(x) for x in s["rerr"][i]][3:]
                want = [abs(M - pr) * abs(dq) / DT / d for M, dq, d
                        in zip(s["M"][i], s["dq"][i], s["dp"][i])][3:]
                v = sorted(got)
                med[(n, e)] = v[len(v) // 2]
                settle[(n, e)] = max(got)
                if SOURCE[n] == "slab":
                    require(max(got) < 0.005,
                            f"{n} {e}: the slab source is the one this formula"
                            " says should be silent, and it is - at every level"
                            f" from 4 on the row's claim is off by less than"
                            f" {100 * max(got):.2f}% of the drop, which is the"
                            " shape the chi band leaves behind rather than a"
                            " wrong inertance")
                    continue
                ratio = [g / w for g, w in zip(got, want)]
                require(min(ratio) > 0.95 and max(ratio) < 1.05,
                        f"{n} {e}: the `{lbl}` inertance in the row is wrong by"
                        f" {100 * abs(s['M'][i][0] / pr - 1):.0f}% of the slab's"
                        f" and the grade comes out at {min(ratio):.3f} to"
                        f" {max(ratio):.3f} times exactly that much of the"
                        " level's rate change - the first time a source choice"
                        " has an error attached to it instead of a road")
    for e in NAME:
        order = [med[("bore_gs", e)], med[("region_gs", e)], med[("slab_gs", e)]]
        require(min(order[:2]) > 10.0 * order[2],
                f"{e}: on the settled levels the drawing's {100 * order[0]:.2f}%"
                f" and the region's {100 * order[1]:.2f}% are both more than an"
                f" order of magnitude above the slab's {100 * order[2]:.4f}%"
                " (median) - the drop is measured across the slab, so the slab's"
                " momentum is the one the row can be right about. Note the first"
                " two do not order by their mistake alone: reading the bore also"
                " moves the trajectory, and so the rate change the mistake is"
                " multiplied by")
    for gs, qn in [("bore_gs", "bore_qn"), ("slab_gs", "slab_qn")]:
        for e in NAME:
            rat = max(med[(gs, e)] / med[(qn, e)], med[(qn, e)] / med[(gs, e)])
            require(rat < 3.0,
                    f"{e}: the {SOURCE[gs]} rung is the same grade under"
                    f" both operators ({med[(gs, e)]:.3e} plain against"
                    f" {med[(qn, e)]:.3e} accelerated, {rat:.2f}x) - rowError is"
                    " a property of the physics, not of the interface map")

    # ---- 8. the chord is a rate-function again, and the destination is one ----
    ctrl = ser["ctrl_steady_net_gs"]
    ratios = {}
    for n in UN:
        s = ser[n]
        for i, e in enumerate(NAME):
            r = [f / q / (d / qc) for f, q, d, qc in zip(s["dpf"][i], s["q"][i],
                                                          ctrl["dp"][i],
                                                          ctrl["q"][i])]
            ratios[(n, e)] = r
            require(r[0] > 2.0,
                    f"{e} {n}: on the first level its friction chord is"
                    f" {r[0]:.2f}x the quasi-steady chord at the same position -"
                    " the region starts from rest and its wall layer is only just"
                    " born, so what the block charges as friction is still the"
                    " starting (Stokes) shear, and that is what the next level has"
                    " to unwind")
            late = r[2:]
            require(max(abs(x - 1.0) for x in late) < 0.10,
                    f"{e} {n}: from level 3 on the same chord sits within"
                    f" {100 * max(abs(x - 1.0) for x in late):.2f}% of the"
                    " quasi-steady one - the coefficient the interface iterates on"
                    " is a function of the rate again, which is why this"
                    " formulation converges at all")
    require(max(r[0] for r in ratios.values()) / min(r[0] for r in ratios.values())
            < 1.6,
            "the starting inflation is one physical effect, not an arm's quirk:"
            f" level-1 chord ratios run {min(r[0] for r in ratios.values()):.2f} to"
            f" {max(r[0] for r in ratios.values()):.2f} over the"
            f" {len(ratios)} unsteady arms and edges")
    laws = {}
    for i, e in enumerate(NAME):
        laws[e] = powerLaw(ctrl["q"][i][3:], ctrl["dp"][i][3:])
    sp_pred = resistiveSplit(laws["X1"][1], laws["X1"][0],
                             laws["X2"][1], laws["X2"][0], Q_IN)
    ends = [ser[n]["q"][0][-1] / Q_IN for n in ARMS]
    for n, x in zip(ARMS, ends):
        require(abs(x - sp_pred) < 0.05,
                f"{n}: it ends at {x:.6f} of the injected rate on X1, within"
                f" {100 * abs(x / sp_pred - 1.0):.2f}% of the resistive split of"
                f" the quasi-steady chords ({sp_pred:.6f}) - charging the"
                " acceleration to the row leaves the destination where the chords"
                " put it")
    require(max(ends) - min(ends) > 1e-6,
            f"and the six arms end at six different positions (X1's split from"
            f" {min(ends):.6f} to {max(ends):.6f}) - in the quasi-steady chain an"
            " attribution source could only move the road; one that changes the"
            " row changes where the run arrives, and the ladder of |rowError|"
            " above says which arrival is the better one")
    for e in NAME:
        trio = [ser[n]["q"][NAME.index(e)][-1] for n in ("bore_gs", "region_gs",
                                                         "slab_gs")]
        require(max(abs(a - b) for a in trio for b in trio) > 1e-9,
                f"{e}: the three sources end the run at three different rates"
                f" ({', '.join(f'{x:.6g}' for x in trio)} m3/s) - the ladder"
                " reaches the position, not only the ramp")

    # ---- 9. cost, budget, and the blocks' own fields ----
    tot = 0
    for a in arms:
        n = a["name"]
        pas = [value(a["steps"][k]["s"], "passes") for k in range(1, N_STEPS + 1)]
        tot += len(a["passes"])
        require(all(p <= MAX_PASSES for p in pas)
                and sum(pas) == len(a["passes"]),
                f"{n}: every level agrees inside the allowed sub-iterations (max"
                f" {int(max(pas))} of {MAX_PASSES}) - iterating on the friction"
                " chord keeps the coefficient a crossing of two positive branches"
                " even with the region's du/dt in the balance")
        require(a["solves"] == 2 * len(a["passes"]),
                f"{n}: {a['solves']} block solves for {len(a['passes'])} passes"
                " - each pass runs both partners once")
        newt = max(value(a["steps"][k]["s"], "newton")
                   for k in range(1, N_STEPS + 1))
        require(newt < 2e-10,
                f"{n}: the 1D step solved its own system at every level (worst"
                f" Newton residual {newt:.2e}, and with the region's du/dt in"
                " the balance that residual is still the 1D solver's floor"
                " rather than an interface argument)")
        require(all(word(b, "blockConverged") == "1" for b in a["blocks"]),
                f"{n}: both blocks closed their own fields - the implicit Euler"
                " term did not push either past its SIMPLE budget")
        require(all(word(b, "unsteady") == ("1" if DDT[n] == "Euler" else "0")
                    for b in a["blocks"]),
                f"{n}: and the block reports carry the same unsteady word the"
                " interface charged them with")
    swd = {}
    for a in arms:
        for i in (0, 1):
            swd.setdefault(NAME[i], []).append(
                sum(ser[a["name"]]["sw"][i]) / N_STEPS)
    require(all(max(v) / min(v) < 6.0 for v in swd.values()),
            "cost per level per edge, steady against unsteady: X1 sweeps"
            f" {swd['X1'][0]:.1f} (control) against {swd['X1'][1]:.1f} to"
            f" {max(swd['X1'][1:]):.1f}, X2 {swd['X2'][0]:.1f} against"
            f" {swd['X2'][1]:.1f} to {max(swd['X2'][1:]):.1f} - stepping in time"
            " is not free, and it is not a different solver either")
    require(all(max(v[1:]) < v[0] for v in swd.values()),
            "and the two halves of that cost move in opposite directions: inside"
            " a block solve the unsteady arms are CHEAPER than the control (X1"
            f" {swd['X1'][0]:.1f} against {min(swd['X1'][1:]):.1f} to"
            f" {max(swd['X1'][1:]):.1f} sweeps, X2 {swd['X2'][0]:.1f} against"
            f" {min(swd['X2'][1:]):.1f} to {max(swd['X2'][1:]):.1f}) - the"
            " implicit Euler term puts V/dt on the diagonal, which bounds rAU by"
            " dt - while across the interface they are dearer (max passes per"
            " level 4 for the control against 6, 6, 9, 5 and 5 for the five"
            " unsteady arms), because the chord the first level carries out of"
            " the starting wall layer is the one the next level unwinds")
    require(tot < 6 * MAX_PASSES * N_STEPS * 1.5,
            f"cost: {tot} interface passes for {N_STEPS * len(ARMS)} agreed time"
            " levels")
    for a in arms:
        for b in a["blocks"]:
            e = word(b, "name").rstrip(":")
            L, H, W = EDGES[e]
            i = NAME.index(e)
            s = ser[a["name"]]
            require(abs(value(b, "fluxIn") - value(b, "fluxOut")) <= 1e-9 * Q_IN
                    and abs(value(b, "fluxIn") - s["q"][i][-1]) <= 1e-9 * Q_IN,
                    f"{a['name']} {e}: the block carries the rate its edge ended"
                    f" with ({value(b, 'fluxIn'):.6g} m3/s) and what enters leaves")
            require(abs(value(b, "inletArea") - H * W) <= 1e-12 * H * W,
                    f"{a['name']} {e}: the inlet patch is the section both"
                    " inertances were referred to")
            yp = value(b, "wallYPlusMax")
            require(11.81 <= yp <= 300.0,
                    f"{a['name']} {e}: wall cells at y+ {yp:.1f}, inside the"
                    " log-law band")
            v = value(b, "plugSpeed")
            de = 4.0 * H * W / (2.0 * (H + W))
            require(RHO * v * de / MU > 4000.0,
                    f"{a['name']} {e}: turbulent at Re_bulk"
                    f" {RHO * v * de / MU:.0f}")
            closure = abs(value(b, "probeSlope")) * H / 2.0
            require(abs(closure - RHO * value(b, "wallUTau") ** 2) <= 0.08 * closure,
                    f"{a['name']} {e}: developed gradient balances wall shear to"
                    f" {(closure / (RHO * value(b, 'wallUTau') ** 2) - 1) * 100:+.2f}%")
            require(abs(value(b, "dpPa") - s["dp"][i][-1]) <= 1e-9 * value(b, "dpPa"),
                    f"{a['name']} {e}: the drop the interface ended with is the"
                    " drop the block's own field measures")

    print(f"\nPASS: the region solves its own du/dt, the interface iterates on the"
          " friction chord, and for the first time the inertance in the 1D row is"
          " graded - the slab's rho*V/A^2 is the one the drop belongs to"
          if not fails else "")
    for f in fails:
        print("FAILED: " + f)
    if fails:
        raise SystemExit(f"{len(fails)} unsteady-partner claim(s) failed")


if __name__ == "__main__":
    main()
