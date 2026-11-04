#!/usr/bin/env python3
"""Run chflow-piso-foam on the closed-box gravity case and verify the p_rgh
variable change does what the plain pressure formulation cannot: hold a closed
box at rest under a uniform body force.

The exact answer is U = 0 everywhere and p = rho g (H - y). The patch-uniform
pressure BCs cannot represent the linear hydrostatic field, so without p_rgh
the same box creeps to uMax 1.82663170862 m/s and still growing at t = 0.5.
With pRgh the uniform force is absorbed into p_rgh = p - force.x, the constant
per patch IS the right BC, and the solved state is exact:

  uMaxAll     0
  prghSpread  0
  pTrueRange  9.81 * 0.975   (the head over the cell-centre y extent)

usage: check_prgh.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

GRAV, H, END_TIME, DT, CELLS_Y = 9.81, 1.0, 0.5, 0.01, 40
# cell centres run from dy/2 to H - dy/2; the reconstructed true pressure
# spans g times exactly that extent
P_TRUE_RANGE = GRAV * (H - H / CELLS_Y)

GOLDEN = {
    "uMaxAll": 0.0,
    "prghSpread": 0.0,
    "pTrueRange": 9.56475,
    "worstContinuity": 0.0,
}


def parse(stdout):
    prgh = re.search(r"p_rgh: uMaxAll (\S+) prghSpread (\S+) pTrueRange (\S+)",
                     stdout)
    done = re.search(r"steps (\d+) endTime (\S+) maxCorrectors (\d+) "
                     r"worstContinuity (\S+)", stdout)
    prof = [(float(a), float(b)) for a, b in
            re.findall(r"profile y (\S+) u (\S+)", stdout)]
    return prgh, done, prof


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True,
                       cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr)
        raise SystemExit("solver failed")

    prgh, done, prof = parse(r.stdout)
    if not prgh or not done:
        raise SystemExit("could not parse the solver output")
    steps, endTime, maxCorr, worst = (int(done[1]), float(done[2]),
                                      int(done[3]), float(done[4]))
    uMaxAll, prghSpread, pTrueRange = (float(g) for g in prgh.groups())

    checks = [
        (uMaxAll == GOLDEN["uMaxAll"],
         f"closed box not at rest: uMaxAll {uMaxAll:.3e}"),
        (prghSpread == GOLDEN["prghSpread"],
         f"p_rgh not constant: spread {prghSpread:.3e}"),
        (abs(pTrueRange - P_TRUE_RANGE) < 1e-9,
         f"reconstructed head {pTrueRange:.6g} != {GRAV} * {H - H / CELLS_Y} "
         f"= {P_TRUE_RANGE:.6g}"),
        (abs(pTrueRange - GOLDEN["pTrueRange"]) < 1e-9,
         f"pTrueRange drifted from the golden anchor "
         f"{GOLDEN['pTrueRange']}"),
        (prof and max(abs(u) for _, u in prof) == 0.0,
         "final profile is not identically zero"),
        (worst == GOLDEN["worstContinuity"],
         f"worst continuity error {worst:.3e}"),
        (maxCorr <= 2, f"corrector loop not converging: {maxCorr} correctors"),
        (steps == round(END_TIME / DT), f"unexpected step count {steps}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    for msg in failed:
        print("FAIL:", msg)
    if failed:
        raise SystemExit("p_rgh verification FAILED")
    print(f"p_rgh verification OK: uMaxAll {uMaxAll:.3e} prghSpread "
          f"{prghSpread:.3e} pTrueRange {pTrueRange:.6g} "
          f"({P_TRUE_RANGE:.6g} expected) steps {steps} maxCorrectors "
          f"{maxCorr} worstContinuity {worst:.3e}")


if __name__ == "__main__":
    main()
