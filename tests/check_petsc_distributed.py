#!/usr/bin/env python3
"""Cross-check PETSc's own MPI machinery against the homegrown distributed CG
(README M6 / A11, §5.133).

Both apps assemble the SAME global graph-Laplacian CSR from the case and solve
b = A*xExact for a known polynomial, but through genuinely independent
distributed paths:

  chflow-mpi-linear-solver      homegrown decomposition (coordinate-median
                                block partition) + diagonal-preconditioned CG
                                with Sendrecv halo exchange + Allreduce dots
  chflow-3d-petsc-distributed   PETSc itself (MatCreateAIJ + KSPCG/PCJACOBI)
                                under PETSc's default contiguous-block
                                ownership -- deliberately a different partition

Each run dumps the same six-column format (cell, xDist/xPetsc, xSerial, xExact,
centre xyz). This checker runs BOTH apps under `mpiexec -np N` for every N in
`mpi ranks`, then verifies from the logs AND the per-cell dumps:

  1. every rank converged, with a sane iteration count (both apps)
  2. distributed and serial solutions == the analytic solution, re-derived
     here from `exactCoeffs` and the dump centres (both apps)
  3. rank-count invariance of each app's distributed solution
  4. CROSS: the PETSc distributed solution == the homegrown distributed
     solution on the same cells, and the two serial reference columns agree
     too -- two independent decompositions of the same system converge to the
     same answer

The dumps of the two apps share a file name (`<base>_<n>.dat`), so every run
executes in a private temporary working directory and the checker reads each
dump before it is overwritten.

usage: check_petsc_distributed.py <mpiexec> <petsc_app> <mpi_solver> <case.dict>
"""
import os
import re
import subprocess
import sys
import tempfile

# ------------------------------------------------------------------ case file

def strip_comments(text):
    out = []
    for line in text.splitlines():
        i = line.find('//')
        if i >= 0:
            line = line[:i]
        out.append(line)
    return "\n".join(out)


def get_block(text, name):
    m = re.search(r'\b' + name + r'\s*\{', text)
    if m is None:
        return ''
    open_pos = text.index('{', m.start())
    depth = 0
    for i in range(open_pos, len(text)):
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
            if depth == 0:
                return text[open_pos + 1:i]
    raise SystemExit("case: unbalanced block %s" % name)


def list_value(body, key):
    m = re.search(r'\b' + key + r'\s*\(([^)]*)\)', body)
    return None if m is None else [float(v) for v in m.group(1).replace(';', ' ').split()]


def scalar_value(body, key):
    m = re.search(r'\b' + key + r'\s*([-+0-9.eE]+)\s*;', body)
    return None if m is None else float(m.group(1))


def int_value(body, key):
    m = re.search(r'\b' + key + r'\s*(\d+)\s*;', body)
    return None if m is None else int(m.group(1))


def word_value(body, key):
    m = re.search(r'\b' + key + r'\s*([A-Za-z0-9_./+-]+)\s*;', body)
    return None if m is None else m.group(1)


# ------------------------------------------------------------------- log lines

RE_PETSC_START = (r'\[chflow\] petscDist start ranks (\d+) cells (\d+) nnz (\d+) '
                  r'shift (\S+) tol (\S+) maxIter (\d+)')
RE_PETSC_SOLVE = (r'\[chflow\] petscDist solve ranks (\d+) iters (\d+) serIters (\d+) '
                  r'converged (\d+) reason (\d+) resNorm (\S+) bnorm (\S+)')
RE_PETSC_VERIFY = (r'\[chflow\] petscDist verify ranks (\d+) analyticDiff (\S+) '
                   r'petscSerialDiff (\S+) serAnalyticDiff (\S+)')

RE_MPI_START = (r'\[chflow\] mpiLinearSolver start ranks (\d+) cells (\d+) nnz (\d+) '
                r'shift (\S+) tol (\S+) maxIter (\d+)')
RE_MPI_SOLVE = (r'\[chflow\] mpiLinearSolver solve ranks (\d+) iters (\d+) serIters (\d+) '
                r'converged (\d+) convRanks (\d+)/(\d+) resRel (\S+) matvecDev (\S+)')
RE_MPI_VERIFY = (r'\[chflow\] mpiLinearSolver verify ranks (\d+) maxDiff (\S+) '
                 r'analyticDiff (\S+) serAnalyticDiff (\S+) bnorm (\S+)')


def read_dump(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            rows.append(tuple(float(v) for v in parts))
    return rows


def poly(c, x, y, z):
    return (c[0] + c[1] * x + c[2] * y + c[3] * z
            + c[4] * x * x + c[5] * y * y + c[6] * z * z
            + c[7] * x * y + c[8] * y * z + c[9] * z * x)


def maxdiff(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def run_and_dump(mpiexec, app, case_path, n, base, workdir):
    """Run one app at -np n inside workdir; return (log, dump rows)."""
    proc = subprocess.run([mpiexec, '-np', str(n), app, case_path],
                          cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("mpiexec -np %d %s exited %d\n%s%s"
                         % (n, os.path.basename(app), proc.returncode,
                            proc.stdout, proc.stderr))
    path = os.path.join(workdir, '%s_%d.dat' % (base, n))
    if not os.path.exists(path):
        raise SystemExit("the -np %d %s run did not write %s"
                         % (n, os.path.basename(app), path))
    rows = read_dump(path)
    os.remove(path)
    return proc.stdout + proc.stderr, rows


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    mpiexec = os.path.abspath(sys.argv[1])
    petsc_app = os.path.abspath(sys.argv[2])
    mpi_solver = os.path.abspath(sys.argv[3])
    case_path = os.path.abspath(sys.argv[4])
    text = strip_comments(open(case_path).read())

    block = get_block(text, 'block')
    mpi_block = get_block(text, 'mpi')
    ls = get_block(text, 'linearSolver')
    verify = get_block(text, 'verify')
    out_block = get_block(text, 'output')

    cells = list_value(block, 'cells')
    ranks = [int(v) for v in (list_value(mpi_block, 'ranks') or [])]
    coeffs = list_value(ls, 'exactCoeffs')
    tol = scalar_value(ls, 'tol')
    max_iter = int_value(ls, 'maxIter')
    shift = scalar_value(ls, 'shift')
    dump = word_value(out_block, 'dump') or 'solution.dat'

    tol_analytic = scalar_value(verify, 'analyticDiffTol')
    tol_dist = scalar_value(verify, 'maxDiffTol')
    if (not cells or not ranks or coeffs is None or len(coeffs) != 10
            or tol is None or max_iter is None or shift is None
            or tol_analytic is None or tol_dist is None):
        raise SystemExit("case is missing block/mpi/linearSolver/verify entries")

    n_cells = int(round(cells[0] * cells[1] * cells[2]))
    base = dump if '.' not in dump else dump[:dump.rfind('.')]

    with tempfile.TemporaryDirectory(prefix='cf_petsc_xcheck_') as workdir:
        # Each app runs in the same private workdir, one rank count at a time;
        # the dump is read (and removed) before the next run can overwrite it.
        petsc_runs = {}   # n -> (log, rows)
        mpi_runs = {}     # n -> (log, rows)
        for n in ranks:
            petsc_runs[n] = run_and_dump(mpiexec, petsc_app, case_path,
                                         n, base, workdir)
            mpi_runs[n] = run_and_dump(mpiexec, mpi_solver, case_path,
                                       n, base, workdir)

    failures = []

    def check(ok, what):
        if not ok:
            failures.append(what)

    # Per-app, per-rank checks; returns the distributed column for the
    # rank-invariance and cross-app comparisons below.
    def audit(app_name, re_start, re_solve, re_verify, runs):
        base_dist, base_ser = None, None
        for n in ranks:
            out, rows = runs[n]
            st = re.search(re_start, out)
            sv = re.search(re_solve, out)
            vr = re.search(re_verify, out)
            if st is None or sv is None or vr is None:
                raise SystemExit("ranks=%d %s: missing log lines\n%s"
                                 % (n, app_name, out[-3000:]))

            check(int(st.group(1)) == n,
                  "%s ranks=%d: start line reports %s ranks" % (app_name, n, st.group(1)))
            check(int(st.group(2)) == n_cells,
                  "%s ranks=%d: start cells %s, the block says %d"
                  % (app_name, n, st.group(2), n_cells))
            check(int(st.group(3)) > n_cells,
                  "%s ranks=%d: nnz %s is implausibly small" % (app_name, n, st.group(3)))
            check(float(st.group(4)) == shift and float(st.group(5)) == tol
                  and int(st.group(6)) == max_iter,
                  "%s ranks=%d: logged solver parameters differ from the case"
                  % (app_name, n))

            iters, ser_iters = int(sv.group(2)), int(sv.group(3))
            conv = int(sv.group(4))
            check(conv == 1,
                  "%s ranks=%d: the distributed solve did not converge" % (app_name, n))
            check(1 <= iters <= max_iter,
                  "%s ranks=%d: iteration count %d is not in [1,%d]"
                  % (app_name, n, iters, max_iter))
            check(abs(iters - ser_iters) <= 2,
                  "%s ranks=%d: distributed iters %d vs serial iters %d"
                  % (app_name, n, iters, ser_iters))

            if app_name == 'petsc':
                check(int(vr.group(1)) == n,
                      "%s ranks=%d: verify line ranks mismatch" % (app_name, n))
                check(float(sv.group(6)) > 0, "%s ranks=%d: bnorm is zero" % (app_name, n))
            else:
                conv_ranks, conv_of = int(sv.group(5)), int(sv.group(6))
                check(conv_ranks == n == conv_of,
                      "%s ranks=%d: only %d/%d ranks converged" % (app_name, n, conv_ranks, conv_of))
                check(float(vr.group(5)) > 0, "%s ranks=%d: bnorm is zero" % (app_name, n))

            # ---- independently re-derived from the per-cell dump ----
            check(len(rows) == n_cells,
                  "%s ranks=%d: the dump has %d rows, expected %d"
                  % (app_name, n, len(rows), n_cells))
            check([int(r[0]) for r in rows] == list(range(n_cells)),
                  "%s ranks=%d: the dump is not in global cell order 0..%d"
                  % (app_name, n, n_cells - 1))
            x_dist = [r[1] for r in rows]
            x_ser = [r[2] for r in rows]
            x_exact_dump = [r[3] for r in rows]
            centres = [(r[4], r[5], r[6]) for r in rows]
            x_exact = [poly(coeffs, *p) for p in centres]

            d_da = maxdiff(x_dist, x_exact)
            d_sa = maxdiff(x_ser, x_exact)
            d_ds = maxdiff(x_dist, x_ser)
            d_da_dump = maxdiff(x_exact_dump, x_exact)
            check(d_da <= tol_analytic,
                  "%s ranks=%d: distributed vs analytic solution max diff %g (tol %g)"
                  % (app_name, n, d_da, tol_analytic))
            check(d_sa <= tol_analytic,
                  "%s ranks=%d: serial vs analytic solution max diff %g (tol %g)"
                  % (app_name, n, d_sa, tol_analytic))
            check(d_ds <= tol_dist,
                  "%s ranks=%d: distributed vs serial column max diff %g (tol %g)"
                  % (app_name, n, d_ds, tol_dist))
            check(d_da_dump <= 1e-12,
                  "%s ranks=%d: the dumped analytic column is not the case "
                  "polynomial (off %g)" % (app_name, n, d_da_dump))

            if base_dist is None:
                base_dist, base_ser = x_dist, x_ser
            else:
                d_rank = maxdiff(x_dist, base_dist)
                check(d_rank <= tol_dist,
                      "%s ranks=%d: the solution differs from the %d-rank run "
                      "by %g (tol %g)" % (app_name, n, ranks[0], d_rank, tol_dist))
                d_ser = maxdiff(x_ser, base_ser)
                check(d_ser <= tol_dist,
                      "%s ranks=%d: the serial column differs from the "
                      "%d-rank run by %g" % (app_name, n, ranks[0], d_ser))
        return base_dist, base_ser

    petsc_dist, _ = audit('petsc', RE_PETSC_START, RE_PETSC_SOLVE,
                          RE_PETSC_VERIFY, petsc_runs)
    homegrown_dist, _ = audit('mpi', RE_MPI_START, RE_MPI_SOLVE,
                              RE_MPI_VERIFY, mpi_runs)

    # ---- CROSS: two independent distributed paths on the same system ----
    for n in ranks:
        _, petsc_rows = petsc_runs[n]
        _, mpi_rows = mpi_runs[n]
        d_cross = maxdiff([r[1] for r in petsc_rows], [r[1] for r in mpi_rows])
        check(d_cross <= tol_dist,
              "ranks=%d: the PETSc and homegrown distributed solutions differ "
              "by %g (tol %g)" % (n, d_cross, tol_dist))
        d_ser = maxdiff([r[2] for r in petsc_rows], [r[2] for r in mpi_rows])
        check(d_ser <= tol_dist,
              "ranks=%d: the two apps' serial reference columns differ by %g "
              "(tol %g)" % (n, d_ser, tol_dist))
    check(petsc_dist is not None and homegrown_dist is not None, "both apps ran")

    if failures:
        raise SystemExit("PETSc distributed cross-check FAILED:\n  "
                         + "\n  ".join(failures))
    print("petsc distributed cross-check passed on %s for ranks %s "
          "(PETSc KSPCG vs homegrown CSR CG, cells=%d)"
          % (os.path.basename(case_path), ranks, n_cells))


if __name__ == '__main__':
    main()
