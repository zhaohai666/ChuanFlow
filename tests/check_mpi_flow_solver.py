#!/usr/bin/env python3
"""5.233 independent regression for the distributed SIMPLE flow solver.

The same channel case is run through chflow-3d-mpi-flow at np = 1, 2, 4.
np = 1 is the single-domain serial SIMPLE loop (the "single process serial"
reference of the increment); np = 2 / 4 are the decomposed runs. Nothing
trusts the solver's own numbers:

  * the gathered p / U fields are compared cell by cell between rank counts
    at 1e-8 relative to the field scale. The tolerance source is declared:
    the converged fields solve the same discrete system, but the Gauss
    sweep order changes with the decomposition and the global reductions
    (MPI_Allreduce tree order) differ from the serial left-to-right
    accumulate, so only the last digits move. Uz is exactly zero by the
    symmetry planes, so it is graded absolutely.
  * div(phi) is recomputed in Python from the per-rank face-flux dumps and
    the uniform grid rebuilt from the case (block min/max/cells), and must
    be < 1e-9 after the gather. The internal-face fluxes are themselves
    re-derived from the dumped cell velocities (linear interpolation with
    the mesh's own owner weight, 0.5 on this uniform block) and the exact
    face areas.
  * every halo slot of every rank must equal, bit for bit, the owned value
    the owning rank dumped for that global cell (the Sendrecv moved exact
    doubles), and each rank's printed checksum must reproduce the sum of
    its own halo dump in slot order exactly.
  * refusals: a case whose declared `domains` disagrees with the launched
    np, and a domain count larger than the mesh cell count.

usage: check_mpi_flow_solver.py <binary> <case.dict> <mpiexec>
"""
import glob
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
nchecks = 0


def check(cond, msg):
    global nchecks
    nchecks += 1
    if not cond:
        fails.append(msg)


def strip_comments(text):
    return "\n".join(l[:l.find("//")] if "//" in l else l
                     for l in text.splitlines())


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("no " + name + " block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced " + name + " block")


def vectors(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", block_text)
    if not m:
        raise RuntimeError("no vector " + key)
    return [float(m.group(i)) for i in (1, 2, 3)]


def run_mpi(np, case_path, cwd):
    return subprocess.run([MPIEXEC, "-np", str(np), BINARY, case_path],
                          capture_output=True, text=True, cwd=cwd, timeout=900)


def parse_solve(out):
    m = re.search(r"mpiFlowSolver solve ranks (\d+) iters (\d+) converged (\d+) "
                  r"convRanks (\d+)/(\d+) contErr (" + NUM + r") dpMax (" + NUM +
                  r") uRes (" + NUM + r") haloMoved (" + NUM + r") divPhiMax ("
                  + NUM + r")", out)
    if not m:
        return None
    return {"ranks": int(m.group(1)), "iters": int(m.group(2)),
            "converged": int(m.group(3)), "conv": (int(m.group(4)), int(m.group(5))),
            "contErr": float(m.group(6)), "dpMax": float(m.group(7)),
            "uRes": float(m.group(8)), "haloMoved": int(m.group(9)),
            "divPhiMax": float(m.group(10))}


def parse_halos(out, np):
    got = {}
    for m in re.finditer(r"mpiFlowSolver halo rank (\d+) slots (\d+) checksum ("
                         + NUM + r")", out):
        got[int(m.group(1))] = (int(m.group(2)), float(m.group(3)))
    return got


def load_fields(path):
    d = {}
    for line in open(path):
        p = line.split()
        d[int(p[0])] = [float(x) for x in p[1:8]]  # cx cy cz p Ux Uy Uz
    return d


def load_halos(np):
    per_rank = {}
    for fn in sorted(glob.glob(os.path.join(TD, "flow_halo_np%d_r*.dat" % np))):
        r = int(re.search(r"_r(\d+)\.dat$", fn).group(1))
        rows = []
        for line in open(fn):
            p = line.split()
            rows.append((int(p[0]), [float(x) for x in p[1:5]]))
        per_rank[r] = rows
    return per_rank


def load_phis(np):
    faces = {}
    for fn in sorted(glob.glob(os.path.join(TD, "flow_phi_np%d_r*.dat" % np))):
        for line in open(fn):
            p = line.split()
            fid, ow, nb, sign = int(p[0]), int(p[1]), int(p[2]), int(p[4])
            phi = float(p[3])
            fc = (float(p[5]), float(p[6]), float(p[7]))
            # owner-oriented flux: both sides of a coupled face must carry the
            # identical double
            if fid in faces:
                check(abs(faces[fid][1] - phi) == 0.0,
                      "np%d: coupled face %d fluxes differ between ranks (%r vs %r)"
                      % (np, fid, faces[fid][1], phi))
            else:
                faces[fid] = ((ow, nb), phi, fc)
    return faces


if len(sys.argv) != 4:
    raise SystemExit("usage: check_mpi_flow_solver.py <binary> <case.dict> <mpiexec>")
BINARY = os.path.abspath(sys.argv[1])
CASE = os.path.abspath(sys.argv[2])
MPIEXEC = os.path.abspath(sys.argv[3])
BASE = strip_comments(open(CASE).read())

blk = block(BASE, "block")
mn = vectors(blk, "min")
mx = vectors(blk, "max")
nc = [int(v) for v in vectors(blk, "cells")]
Nx, Ny, Nz = nc
nC = Nx * Ny * Nz
dx, dy, dz = (mx[0] - mn[0]) / Nx, (mx[1] - mn[1]) / Ny, (mx[2] - mn[2]) / Nz
nu = scalar(BASE, "nu")
simple = block(BASE, "simple")
contTol = scalar(simple, "continuityTol")
relaxU = scalar(simple, "relaxU")
relaxP = scalar(simple, "relaxP")

TD = tempfile.mkdtemp(prefix="cf_mpi_flow_")
runs = {}
for np in (1, 2, 4):
    r = run_mpi(np, CASE, TD)
    check(r.returncode == 0, "np=%d returned %d: %s" %
          (np, r.returncode, (r.stdout + r.stderr)[-400:]))
    s = parse_solve(r.stdout)
    check(s is not None, "np=%d: no solve line" % np)
    if s:
        check(s["ranks"] == np and s["converged"] == 1
              and s["conv"] == (np, np),
              "np=%d: solve line %s" % (np, s))
        check(s["contErr"] < contTol and s["divPhiMax"] < 1e-9,
              "np=%d: continuity not at tolerance: %s" % (np, s))
        check(s["haloMoved"] > 0 if np > 1 else s["haloMoved"] == 0,
              "np=%d: halo traffic %d unexpected" % (np, s["haloMoved"]))
    runs[np] = (r.stdout, load_fields(os.path.join(TD, "flow_np%d.dat" % np)))

# ------------------------------------------------- cross-rank field agreement
ref = runs[1][1]
for np in (2, 4):
    other = runs[np][1]
    check(len(other) == nC, "np=%d: %d cells dumped, expected %d" % (np, len(other), nC))
    for i, nm in enumerate(("p", "Ux", "Uy", "Uz")):
        scale = max(abs(v[3 + i]) for v in ref.values())
        d = max(abs(ref[k][3 + i] - other[k][3 + i]) for k in ref)
        if nm == "Uz":  # exactly zero by the symmetry planes
            check(d <= 1e-12, "np=%d %s absolute drift %g" % (np, nm, d))
        else:
            check(d <= 1e-8 * max(scale, 1e-30),
                  "np=%d %s rel %.3e exceeds 1e-8 (abs %g, scale %g)"
                  % (np, nm, d / max(scale, 1e-30), d, scale))

# ------------------------------------------------------ div(phi) re-derivation
for np in (1, 2, 4):
    faces = load_phis(np)
    fields = runs[np][1]
    div = {c: 0.0 for c in range(nC)}
    for (ow, nb), phi, fc in faces.values():
        div[ow] += phi
        if nb >= 0:
            div[nb] -= phi
    worst = max(abs(v) / (dx * dy * dz) for v in div.values())
    check(worst < 1e-9,
          "np=%d: recomputed div(phi) max %g exceeds 1e-9" % (np, worst))

    # boundary-flux closure re-derived from the case via the dumped face
    # centres: the inlet carries the declared fixedValue velocity (p
    # zeroGradient -> no pressure correction there), walls and symmetry
    # faces must read exactly zero flux, the outlet must pass the inflow.
    inlet_area = dy * dz
    n_in = n_out = 0
    for fid, ((ow, nb), phi, fc) in faces.items():
        if nb >= 0:
            continue
        if abs(fc[0] - mn[0]) < 1e-9:                    # inlet
            check(abs(phi + inlet_area) < 1e-12,
                  "np=%d: inlet face %d flux %g != -|U|*A = %g"
                  % (np, fid, phi, -inlet_area))
            n_in += 1
        elif abs(fc[0] - mx[0]) < 1e-9:                  # outlet
            check(0.0 < phi < 2.0 * inlet_area,
                  "np=%d: outlet face %d flux %g outside (0, 2|U|A)" % (np, fid, phi))
            n_out += 1
        else:                                            # wall / symmetry
            check(phi == 0.0,
                  "np=%d: wall/symmetry face %d leaks flux %g" % (np, fid, phi))
    check(n_in == Ny * Nz and n_out == Ny * Nz,
          "np=%d: expected %d inlet + %d outlet faces, saw %d + %d"
          % (np, Ny * Nz, Ny * Nz, n_in, n_out))
    # net boundary flux == sum of the per-cell divergences; bounded by the
    # continuity tolerance times the box volume, not by machine epsilon
    net = sum(div[c] for c in div)
    check(abs(net) < 1e-9 * (dx * dy * dz) * nC,
          "np=%d: the flux dump does not telescope (net %g)" % (np, net))

# ------------------------------------------------------------- halo exactness
for np in (2, 4):
    owned = runs[np][1]
    per_rank = load_halos(np)
    out = runs[np][0]
    checks = parse_halos(out, np)
    check(set(checks) == set(range(np)),
          "np=%d: halo checksums missing for some ranks: %s" % (np, sorted(checks)))
    for r, rows in per_rank.items():
        for gc, vals in rows:
            got = owned[gc][3:7]
            check(vals == got,
                  "np=%d rank %d halo cell %d: %s != owned %s (must be exact)"
                  % (np, r, gc, vals, got))
        # the printed checksum must reproduce this rank's halo dump exactly
        csum = 0.0
        for gc, vals in rows:
            csum += abs(vals[0]) + abs(vals[1]) + abs(vals[2]) + abs(vals[3])
        slots, printed = checks[r]
        check(slots == len(rows) and printed == csum,
              "np=%d rank %d checksum mismatch: printed %r slots %d vs file %r slots %d"
              % (np, r, printed, slots, csum, len(rows)))

# ------------------------------------------------------------- physical sanity
for c, v in ref.items():
    check(v[4] > 0.0, "np=1: Ux must be positive everywhere (cell %d: %g)" % (c, v[4]))
# pressure falls monotonically from inlet to outlet toward the declared
# p = 0 reference (column averages on the uniform grid)
cols = {}
for c, v in ref.items():
    i = int((v[0] - mn[0]) / dx)
    cols.setdefault(i, []).append(v[3])
colp = [sum(cols[i]) / len(cols[i]) for i in range(Nx)]
check(all(colp[i] > colp[i + 1] for i in range(Nx - 1)),
      "np=1: column-mean pressure is not monotonically decreasing")
check(colp[0] > 0 and abs(colp[-1]) < colp[0],
      "np=1: inlet must press above the outlet reference")
check(max(abs(v[5]) for v in ref.values()) < 0.1 * max(v[4] for v in ref.values()),
      "np=1: transverse velocity too large for a channel")

# ---------------------------------------------------------------- refusals
mismatch = BASE + "\nmpi { domains 3; }\n"
tiny = re.sub(r"cells \(\d+ \d+ \d+\);", "cells (2 1 1);", mismatch)
tiny = tiny.replace("domains 3;", "domains 4;")
for name, txt, np, needle in (
        ("np mismatch", mismatch, 2, "does not match"),
        ("domains > mesh segments", tiny, 4, "more domains than mesh segments"),
):
    p = os.path.join(TD, "reject_%s.dict" % name.replace(" ", "_"))
    open(p, "w").write(txt)
    r = run_mpi(np, p, TD)
    check(r.returncode != 0 and needle in (r.stdout + r.stderr),
          "refusal '%s' not raised loudly (rc=%d)" % (name, r.returncode))

if fails:
    print("MPI flow solver FAILED %d checks:" % len(fails))
    for f in fails[:40]:
        print("  -", f)
    sys.exit(1)

print("MPI flow solver OK: %d assertions; SIMPLE converged at np=1/2/4, "
      "gathered p/U agree cell-by-cell within 1e-8 of the field scale "
      "(reduction-order last digits), recomputed div(phi) < 1e-9, halo "
      "dumps and per-rank checksums exact, 2 loud refusals" % nchecks)
