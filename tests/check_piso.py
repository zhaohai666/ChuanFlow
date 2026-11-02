#!/usr/bin/env python3
"""Run chflow-piso-foam on the start-up channel and verify the transient
against the analytic series solution, plus the divergence-free corrected flux.

The case is a planar channel (H = 0.2, nu = 0.01) started from rest by a
constant body force a = 0.1, so u = (u(y,t), 0, 0) and

  u(y,t) = a H^2/(2 nu) [ xi(1-xi) - sum_n 4(1-(-1)^n)/(n pi)^3
                            sin(n pi xi) exp(-n^2 pi^2 nu t / H^2) ],  xi = y/H

The steady part is compared against the *discrete* steady state: a cell-centred
finite volume imposes the wall value half a cell away from the first node, so
the solved profile is the analytic parabola raised by (a/2nu)(dy/2)^2. The
checker asserts that shift explicitly instead of absorbing it into a loose
tolerance, which makes it a real regression test of the wall treatment.

usage: check_piso.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

H, NU, FORCE, CELLS_Y = 0.2, 0.01, 0.1, 20
U_STEADY = FORCE * H * H / (2.0 * NU)  # centreline value of the parabola
DY = H / CELLS_Y
# cell-centred finite volume imposes the wall value at the face one half-cell
# away, so the discrete steady state is the analytic parabola raised by the
# curvature over that half cell: (a/2nu)(dy/2)^2 = 1.25e-4 here.
U_EXACT_STEADY = FORCE / (2.0 * NU)


def steady(y):
    return U_EXACT_STEADY * (y * (H - y) + (DY / 2.0) ** 2)


def exact(xi, t, nmax=200):
    """Analytic start-up series with the discrete wall shift in the steady part."""
    s = steady(xi * H)
    for n in range(1, nmax + 1):
        b = 4.0 * (1.0 - (-1.0) ** n) / (n * math.pi) ** 3
        s -= U_STEADY * b * math.sin(n * math.pi * xi) * math.exp(
            -n * n * math.pi * math.pi * NU * t / (H * H))
    return s


def parse(stdout):
    hist = [(float(a), float(b), float(c)) for a, b, c, _ in
            re.findall(r"t (\S+) uMax (\S+) cont (\S+) correctors (\S+)", stdout)]
    prof = [(float(a), float(b)) for a, b in
            re.findall(r"profile y (\S+) u (\S+)", stdout)]
    done = re.search(r"steps (\d+) endTime (\S+) maxCorrectors (\d+) worstContinuity (\S+)",
                     stdout)
    return hist, prof, done


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr)
        raise SystemExit("solver failed")

    hist, prof, done = parse(r.stdout)
    if not hist or not prof or not done:
        raise SystemExit("could not parse the solver output")
    steps, endTime, maxCorr, worst = (int(done[1]), float(done[2]),
                                      int(done[3]), float(done[4]))
    dt = float(re.search(r"(?<!\w)dt\s+([0-9.eE+-]+);", open(case).read())[1])

    # uMax is the centre cell of the (symmetric, x-invariant) profile
    xiProbe = min((y / H for y, _ in prof), key=lambda x: abs(x - 0.5))
    uProbe = steady(xiProbe * H)

    rel = [(t, abs(u - exact(xiProbe, t)) / max(exact(xiProbe, t), 1e-12))
           for t, u, _ in hist]
    histErr = max(e for _, e in rel)
    histErrLate = max(e for t, e in rel if t > endTime / 2.0)

    profErr = max(abs(u - steady(y)) for y, u in prof)
    # the same profile measured against the *continuous* parabola carries the
    # half-cell wall offset (a/2nu)(dy/2)^2, first-order FVM wall treatment
    contErr = max(abs(u - U_EXACT_STEADY * y * (H - y)) for y, u in prof)

    # monotone start-up: the centre cell must rise to its discrete steady value
    # without overshoot
    monotone = all(b[1] >= a[1] - 1e-12 for a, b in zip(hist, hist[1:]))
    overshoot = max(u for _, u, _ in hist) - uProbe
    contMax = max(c for _, _, c in hist)

    print(f"steps={steps} endTime={endTime} maxCorrectors={maxCorr} "
          f"histErr={histErr:.3e} histErrLate={histErrLate:.3e} "
          f"profileErr={profErr:.3e} wallOffset={contErr:.3e} "
          f"contMax={contMax:.3e} worst={worst:.3e} "
          f"overshoot={overshoot:.3e} xiProbe={xiProbe:.4f}")

    checks = [
        # the first reported point carries the accumulated O(dt) Euler lag of a
        # rising solution plus the discrete-mode eigenvalue error; by half the
        # run time only the (tiny) residual transient is left
        (histErr < 2e-2, f"transient history off the series by {histErr:.3e}"),
        (histErrLate < 1e-5, f"late history off the series by {histErrLate:.3e}"),
        (profErr < 1e-6 * uProbe, f"final profile off the discrete steady "
                                  f"profile by {profErr:.3e}"),
        (abs(contErr - U_EXACT_STEADY * (DY / 2.0) ** 2) < 1e-6,
         f"wall offset {contErr:.3e} is not the expected half-cell shift"),
        (contMax < 1e-12, f"corrected flux not divergence-free: {contMax:.3e}"),
        (worst < 1e-12, f"worst continuity error {worst:.3e}"),
        (monotone, "start-up overshoots or oscillates"),
        (abs(overshoot) < 1e-6 * uProbe, f"steady limit mismatch {overshoot:.3e}"),
        (maxCorr <= 2, f"corrector loop not converging: {maxCorr} correctors"),
        (steps == round(endTime / dt), f"unexpected step count {steps}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    for msg in failed:
        print("FAIL:", msg)
    if failed:
        raise SystemExit("PISO verification FAILED")


if __name__ == "__main__":
    main()
