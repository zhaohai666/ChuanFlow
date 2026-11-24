#!/usr/bin/env python3
"""Verify the real-MPI distributed sparse linear solve (A11 / §5.122).

Criterion makeup: 3 independent physical criteria + mirror diagnostics
(regression safety net). The mirror no longer carries the verdict alone.

Independent physical criteria (built from the case + the per-cell dump only):

  I1. Independent residual. The checker assembles the graph Laplacian itself
      from the uniform grid (adjacency re-derived from the dumped cell
      centres by spacing, diagonal = face count + shift, off-diagonal -1),
      builds b = A * xExact from the case polynomial, and holds both the
      distributed and the serial solution to ||A x - b|| <~ 0 -- i.e. the
      solve really is the SPD-shifted Laplacian system, independent of the
      app's own matrix bookkeeping.

  I2. Partitioning never degrades the answer. The distributed-vs-analytic
      error and the distributed-vs-serial gap stay bounded by their first
      rank count's values (and the case tolerances) as np grows.

  I3. Rank-count invariance of the gathered solution (kept from the earlier
      round): the solution does not depend on how the rows were cut up.

Mirror diagnostics (regression safety net):

  M1. distributed solution == serial solveCG on the same global CSR (maxDiff)
  M2. distributed and serial solutions == the analytic solution
  M3. every rank converged, with a sane iteration count
  M4. the distributed matvec IS the serial CSR product (matvecDev)

The analytic field is re-derived here from the case's `exactCoeffs` and the
per-cell centres in the dump, so the app's analytic column is not trusted.

usage: check_mpi_linear_solver.py <mpiexec> <solver> <case.dict>
"""
import os
import re
import subprocess
import sys

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

RE_START = (r'\[chflow\] mpiLinearSolver start ranks (\d+) cells (\d+) nnz (\d+) '
            r'shift (\S+) tol (\S+) maxIter (\d+)')
RE_SOLVE = (r'\[chflow\] mpiLinearSolver solve ranks (\d+) iters (\d+) serIters (\d+) '
            r'converged (\d+) convRanks (\d+)/(\d+) resRel (\S+) matvecDev (\S+)')
RE_VERIFY = (r'\[chflow\] mpiLinearSolver verify ranks (\d+) maxDiff (\S+) '
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


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    mpiexec = os.path.abspath(sys.argv[1])
    solver = os.path.abspath(sys.argv[2])
    case_path = os.path.abspath(sys.argv[3])
    case_dir = os.path.dirname(case_path)
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

    tol_diff = scalar_value(verify, 'maxDiffTol')
    tol_matvec = scalar_value(verify, 'matvecDevTol')
    tol_analytic = scalar_value(verify, 'analyticDiffTol')
    if (not cells or not ranks or coeffs is None or len(coeffs) != 10
            or tol is None or max_iter is None or shift is None
            or tol_diff is None or tol_matvec is None or tol_analytic is None):
        raise SystemExit("case is missing block/mpi/linearSolver/verify entries")

    n_cells = int(round(cells[0] * cells[1] * cells[2]))
    base = dump if '.' not in dump else dump[:dump.rfind('.')]

    # The checker runs with cwd = case_dir, so mpiexec/solver/case must be
    # absolute paths (the caller passes the three absolute paths explicitly).
    outs, dumps = {}, {}
    for n in ranks:
        proc = subprocess.run([mpiexec, '-np', str(n), solver, case_path],
                              cwd=case_dir, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit("mpiexec -np %d solver exited %d\n%s%s"
                             % (n, proc.returncode, proc.stdout, proc.stderr))
        outs[n] = proc.stdout + proc.stderr
        path = os.path.join(case_dir, '%s_%d.dat' % (base, n))
        if not os.path.exists(path):
            raise SystemExit("the -np %d run did not write %s" % (n, path))
        dumps[n] = read_dump(path)
        if len(dumps[n]) != n_cells:
            raise SystemExit("the -np %d dump has %d rows, expected %d"
                             % (n, len(dumps[n]), n_cells))

    failures = []

    def check(ok, what):
        if not ok:
            failures.append(what)

    # ---- I1: the graph Laplacian, assembled independently of the app ----
    bmin = list_value(block, 'min')
    bmax = list_value(block, 'max')
    spacing = [(bmax[a] - bmin[a]) / cells[a] for a in range(3)]
    centres0 = [(r[4], r[5], r[6]) for r in dumps[ranks[0]]]
    deg = [0] * n_cells
    for i in range(n_cells):
        for j in range(i + 1, n_cells):
            off = 0
            axis_hit = False
            for a in range(3):
                dv = abs(centres0[j][a] - centres0[i][a])
                if dv > 1e-9:
                    if abs(dv - spacing[a]) > 1e-9:
                        off = 2
                        break
                    off += 1
                    axis_hit = True
            if off == 1 and axis_hit:
                deg[i] += 1
                deg[j] += 1
    maxdeg = max(deg)
    check(maxdeg > 0, "I1: the re-derived grid adjacency has no faces")

    A = [[0.0] * n_cells for _ in range(n_cells)]
    for i in range(n_cells):
        A[i][i] = deg[i] + shift
    for i in range(n_cells):
        for j in range(i + 1, n_cells):
            dx = [abs(centres0[j][a] - centres0[i][a]) for a in range(3)]
            near = [abs(dx[a] - spacing[a]) <= 1e-9 for a in range(3)]
            zero = [dx[a] <= 1e-9 for a in range(3)]
            if sum(near) == 1 and all(zero[a] or near[a] for a in range(3)):
                A[i][j] = A[j][i] = -1.0
    x_e0 = [poly(coeffs, *p) for p in centres0]
    b_ind = [sum(A[i][j] * x_e0[j] for j in range(n_cells)) for i in range(n_cells)]
    res_bound = (2 * maxdeg + shift) * tol_analytic * 10.0

    def matvec(x):
        return [sum(A[i][j] * x[j] for j in range(n_cells)) for i in range(n_cells)]

    base_dist = None
    base_err = None
    base_gap = None
    for n in ranks:
        out = outs[n]
        st = re.search(RE_START, out)
        sv = re.search(RE_SOLVE, out)
        vr = re.search(RE_VERIFY, out)
        if st is None or sv is None or vr is None:
            raise SystemExit("ranks=%d: missing log lines\n%s" % (n, out[-3000:]))

        check(int(st.group(1)) == n,
              "ranks=%d: start line reports %s ranks" % (n, st.group(1)))
        check(int(st.group(2)) == n_cells,
              "ranks=%d: start cells %s, the block says %d" % (n, st.group(2), n_cells))
        check(int(st.group(3)) > n_cells,
              "ranks=%d: nnz %s is implausibly small" % (n, st.group(3)))
        check(float(st.group(4)) == shift and float(st.group(5)) == tol
              and int(st.group(6)) == max_iter,
              "ranks=%d: logged solver parameters differ from the case" % n)

        iters, ser_iters = int(sv.group(2)), int(sv.group(3))
        conv, conv_ranks, conv_of = int(sv.group(4)), int(sv.group(5)), int(sv.group(6))
        res_rel, matvec_dev = float(sv.group(7)), float(sv.group(8))
        check(conv == 1, "ranks=%d: the distributed solve did not converge" % n)
        check(conv_ranks == n == conv_of,
              "ranks=%d: only %d/%d ranks converged" % (n, conv_ranks, conv_of))
        check(1 <= iters <= max_iter,
              "ranks=%d: iteration count %d is not in [1,%d]" % (n, iters, max_iter))
        check(abs(iters - ser_iters) <= 2,
              "ranks=%d: distributed iters %d vs serial iters %d" % (n, iters, ser_iters))
        check(res_rel <= tol * 10,
              "ranks=%d: relative residual %g exceeds tol %g" % (n, res_rel, tol))
        check(matvec_dev <= tol_matvec,
              "ranks=%d: distributed matvec differs from the serial CSR product by "
              "%g (tol %g)" % (n, matvec_dev, tol_matvec))
        check(int(vr.group(1)) == n, "ranks=%d: verify line ranks mismatch" % n)
        check(float(vr.group(4)) > 0, "ranks=%d: bnorm is zero" % n)

        # ---- independently re-derived from the per-cell dump ----
        rows = dumps[n]
        check([int(r[0]) for r in rows] == list(range(n_cells)),
              "ranks=%d: the dump is not in global cell order 0..%d" % (n, n_cells - 1))
        x_dist = [r[1] for r in rows]
        x_ser = [r[2] for r in rows]
        x_exact_dump = [r[3] for r in rows]
        centres = [(r[4], r[5], r[6]) for r in rows]
        x_exact = [poly(coeffs, *p) for p in centres]

        d_ds = maxdiff(x_dist, x_ser)
        d_da = maxdiff(x_dist, x_exact)
        d_sa = maxdiff(x_ser, x_exact)
        d_da_dump = maxdiff(x_exact_dump, x_exact)
        check(d_ds <= tol_diff,
              "ranks=%d: distributed vs serial solution max diff %g (tol %g)"
              % (n, d_ds, tol_diff))
        check(d_da <= tol_analytic,
              "ranks=%d: distributed vs analytic solution max diff %g (tol %g)"
              % (n, d_da, tol_analytic))
        check(d_sa <= tol_analytic,
              "ranks=%d: serial vs analytic solution max diff %g (tol %g)"
              % (n, d_sa, tol_analytic))
        check(d_da_dump <= 1e-12,
              "ranks=%d: the dumped analytic column is not the case polynomial "
              "(off %g)" % (n, d_da_dump))

        # ---- I1: independent residual against the re-assembled system ----
        r_dist = matvec(x_dist)
        r_ser = matvec(x_ser)
        rd = max(abs(r_dist[i] - b_ind[i]) for i in range(n_cells))
        rs = max(abs(r_ser[i] - b_ind[i]) for i in range(n_cells))
        check(rd <= res_bound,
              "ranks=%d: distributed residual of the independently assembled "
              "Laplacian system is %g (bound %g)" % (n, rd, res_bound))
        check(rs <= res_bound,
              "ranks=%d: serial residual of the independently assembled "
              "Laplacian system is %g (bound %g)" % (n, rs, res_bound))

        if base_dist is None:
            base_dist = x_dist
            base_err, base_gap = d_da, d_ds
        else:
            d_rank = maxdiff(x_dist, base_dist)
            check(d_rank <= tol_diff,
                  "ranks=%d: the solution differs from the %d-rank run by %g (tol %g)"
                  % (n, ranks[0], d_rank, tol_diff))
            # ---- I2: partitioning never degrades the answer ----
            check(d_da <= max(base_err, tol_analytic),
                  "ranks=%d: analytic error %g exceeds the %d-rank error %g"
                  % (n, d_da, ranks[0], base_err))
            check(d_ds <= max(base_gap, tol_diff),
                  "ranks=%d: distributed-vs-serial gap %g exceeds the "
                  "%d-rank gap %g" % (n, d_ds, ranks[0], base_gap))

    if failures:
        raise SystemExit("mpi linear-solver check FAILED:\n  " + "\n  ".join(failures))
    print("mpi linear-solver check passed on %s for ranks %s"
          % (os.path.basename(case_path), ranks))


if __name__ == '__main__':
    main()
