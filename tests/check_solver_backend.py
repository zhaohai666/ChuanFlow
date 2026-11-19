#!/usr/bin/env python3
"""Verify the selectable sparse linear solver backends (README M6 / A11).

Runs chflow-solver-backend on the tridiagonal system whose exact solution is
x[i] = 1, and checks that the native backend converges to that solution. When
Eigen is compiled in (the app reports `eigen eigen`), the checker additionally
requires the Eigen backend to converge and agree with the native solution to
1e-6 in inf-norm. The native path is what every other checker relies on for
bit-identity; this checker only adds a cross-backend physics agreement check.
"""

import re
import subprocess
import sys


NUM = r"(-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)"


def backend_line(text, label):
    """Return the single output line carrying `solver-backend <label> ...`."""
    lines = [ln for ln in text.splitlines()
             if re.search(r"solver-backend " + re.escape(label) + r" ", ln)]
    if len(lines) != 1:
        raise AssertionError("expected one '%s' line, got %d" % (label, len(lines)))
    return lines[0]


def field(line, key):
    m = re.search(r"\b" + re.escape(key) + r" " + NUM, line)
    if not m:
        raise AssertionError("missing '%s' in line: %s" % (key, line))
    raw = m.group(1)
    try:
        return int(raw)
    except ValueError:
        return float(raw)


def main():
    if len(sys.argv) < 3:
        print("usage: check_solver_backend.py <chflow-solver-backend> <case.dict>",
              file=sys.stderr)
        return 2
    solver, case = sys.argv[1], sys.argv[2]
    proc = subprocess.run([solver, case], capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise AssertionError("solver-backend exited %d" % proc.returncode)

    out = proc.stdout
    # Backend name tokens are read directly from the banner.
    m = re.search(r"solver-backend n \d+ tol \S+ eigen (\S+) petsc (\S+)", out)
    if not m:
        # Older app build without the petsc banner field.
        m = re.search(r"solver-backend n \d+ tol \S+ eigen (\S+)", out)
        if not m:
            raise AssertionError("missing backend banner in output:\n" + out)
        eigen_name, petsc_name = m.group(1), "none"
    else:
        eigen_name, petsc_name = m.group(1), m.group(2)

    native = backend_line(out, "native")
    nconv = field(native, "converged")
    nerr = field(native, "err")
    nres = field(native, "residual")
    niter = field(native, "iterations")
    print("native: converged=%d iterations=%d err=%.3e residual=%.3e"
          % (nconv, niter, nerr, nres))

    assert nconv == 1, "native backend did not converge"
    assert nerr < 1e-6, "native backend error %.3e exceeds 1e-6" % nerr
    assert nres < 1e-6, "native backend residual %.3e exceeds 1e-6" % nres

    if eigen_name == "eigen":
        eline = backend_line(out, "eigen")
        econv = field(eline, "converged")
        eerr = field(eline, "err")
        eres = field(eline, "residual")
        eiter = field(eline, "iterations")
        cross = field(backend_line(out, "crossBackendMaxDiff"), "crossBackendMaxDiff")
        print("eigen:  converged=%d iterations=%d err=%.3e residual=%.3e"
              % (econv, eiter, eerr, eres))
        print("cross-backend max diff: %.3e" % cross)
        assert econv == 1, "Eigen backend did not converge"
        assert eerr < 1e-6, "Eigen backend error %.3e exceeds 1e-6" % eerr
        assert eres < 1e-6, "Eigen backend residual %.3e exceeds 1e-6" % eres
        assert cross < 1e-6, "cross-backend diff %.3e exceeds 1e-6" % cross
    else:
        print("eigen backend not compiled in (native-only build); cross-check skipped")

    if petsc_name == "petsc":
        pline = backend_line(out, "petsc")
        pconv = field(pline, "converged")
        perr = field(pline, "err")
        pres = field(pline, "residual")
        piter = field(pline, "iterations")
        pcross = field(backend_line(out, "petscCrossMaxDiff"), "petscCrossMaxDiff")
        print("petsc:  converged=%d iterations=%d err=%.3e residual=%.3e"
              % (pconv, piter, perr, pres))
        print("petsc cross-backend max diff: %.3e" % pcross)
        assert pconv == 1, "PETSc backend did not converge"
        assert perr < 1e-6, "PETSc backend error %.3e exceeds 1e-6" % perr
        assert pres < 1e-6, "PETSc backend residual %.3e exceeds 1e-6" % pres
        assert pcross < 1e-6, "PETSc cross-backend diff %.3e exceeds 1e-6" % pcross
    else:
        print("petsc backend not compiled in (no PETSc library); cross-check skipped")

    print("PASS: sparse solver backend(s) converge to the analytic solution")
    return 0


if __name__ == "__main__":
    sys.exit(main())
