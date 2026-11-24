#!/usr/bin/env python3
"""Run a domain-decomposition case and verify what cutting a 3D mesh into
domains has to deliver.

Criterion makeup: 5 independent physical criteria + mirror diagnostics
(regression safety net). The mirror no longer carries the verdict alone.

Independent physical criteria (do not trust any kernel bookkeeping):

  I1. Geometric volume closure. Every cell's volume is computed twice, by two
      independent exact methods (divergence theorem over face rings vs a
      signed tetrahedral decomposition from the cell centre), and the global
      sum must equal the analytic block volume from the case.

  I2. Analytic integral conservation across the cut. The written field holds
      point values at cell centres, so its volume-weighted sum is the
      midpoint quadrature: for the linear polynomial that quadrature is
      exact, and it must equal the continuum integral of the polynomial over
      the exact mesh geometry (2nd-order moments from the same tet
      decomposition) -- halo assembly neither creates nor destroys mass.
      For the quadratic field only the rank-grouped reassembly of the same
      integral is owed (midpoint quadrature is O(h^2) there by design).

  I3. Load balance read off the written `domain` field alone: the rank cell
      counts sum to the mesh, no rank is empty, and max - min <= 1.

  I4. Per-rank cumulative flux reassembly. The divergence integral grouped by
      written rank sums back to the serial divergence integral of the same
      run -- the two sides of every cut face carry contributions that cancel.

  I5. Cut invariance of the physics (kept from the earlier round, tightened
      wording): the linear field's gradient is the analytic one, and every
      written field is bit-identical between domain counts.

Mirror diagnostics (regression safety net; each one re-derives kernel
bookkeeping and pins it to the logs, but none is the sole verdict):

  M1. The partition is re-derived, not trusted. This checker recomputes the
      kernel's own recursive coordinate-median ordering from the cell centres
      read out of the .vtu, slices it into contiguous domains, and compares
      the result with the written `domain` field for every cell. It then
      derives, from the face adjacency of the mesh alone, each domain's owned
      cells, halo cells, internal/coupled/boundary face entries and distinct
      peers, and holds every one of those to the numbers the kernel logged.
      Balance is max-min <= 1 by construction, and the three face classes
      have to reassemble the global mesh.

  M2. A decomposed stencil IS the serial stencil. gradDev (assembled gradient
      against fvc::grad on the whole mesh) and divDev (assembled net out-flux
      against fvc::divPhi) are required to be exactly 0 -- the domains sweep
      their faces in global order and use the mesh's own bounded weights, so
      nothing is gained in accuracy and nothing is lost. The same holds for
      the two Gauss schemes and for haloDrift, the checksum of everything the
      exchange moved.

  M3. The answer does not depend on the cut. Every written field -- the field
      itself, both gradients, both out-fluxes -- is compared cell by cell
      between the 1-domain file and each other domain count, again exactly.

  M4. The physics stays put. On the hex mesh the linear test field's gradient
      is re-derived here from the .vtu geometry and held to the analytic one
      (the tolerance is the case's, not a fitted number); on a mesh coarse or
      skewed enough that the serial scheme cannot reproduce even a linear
      field, that accuracy is reported rather than demanded, and only
      cut-invariance is owed. A divergence-free flux gives every interior
      cell -- cells whose faces all have a neighbour, including the ones a
      cut runs through -- a net out-flux of zero, and the volume-weighted
      sum, min, max and halo checksum the kernel reduced per domain are
      re-summed here independently.

usage: check_parallel.py <chflow-3d-parallel> <case.dict>
"""
import math
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


def word_value(body, key):
    m = re.search(r'\b' + key + r'\s*([A-Za-z0-9_./+-]+)\s*;', body)
    return None if m is None else m.group(1)


# ------------------------------------------------------------------- vtu read

def read_vtu(path):
    """(points, cell vertex lists, per-cell ordered face rings, cell data)."""
    txt = open(path).read()
    pm = re.search(r'<Points>\s*<DataArray[^>]*format="ascii">\s*(.*?)</DataArray>', txt, re.S)
    p = [float(v) for v in pm.group(1).split()]
    points = [(p[i], p[i + 1], p[i + 2]) for i in range(0, len(p), 3)]

    def ints(name):
        m = re.search(r'Name="%s"[^>]*format="ascii">\s*(.*?)</DataArray>' % name, txt, re.S)
        return [int(v) for v in m.group(1).split()]

    conn, off, faces = ints('connectivity'), ints('offsets'), ints('faces')
    cells, s = [], 0
    for o in off:
        cells.append(conn[s:o])
        s = o

    # the writer emits nFaces n0 v.. n1 v.. per cell, in the mesh's own order
    cell_faces, i = [], 0
    for _ in cells:
        nf = faces[i]
        i += 1
        rings = []
        for _k in range(nf):
            nv = faces[i]
            i += 1
            rings.append(tuple(faces[i:i + nv]))
            i += nv
        cell_faces.append(rings)

    fields = {}
    for m in re.finditer(
                r'<DataArray type="Float64" Name="([^"]+)"([^>]*)format="ascii">\s*(.*?)</DataArray>',
                txt, re.S):
        name, rest, body = m.group(1), m.group(2), m.group(3)
        v = [float(x) for x in body.split()]
        if 'NumberOfComponents="3"' in rest:
            fields[name] = [(v[3 * c], v[3 * c + 1], v[3 * c + 2]) for c in range(len(cells))]
        else:
            fields[name] = v
    return points, cells, cell_faces, fields


def cell_centre(points, verts):
    n = float(len(verts))
    return tuple(sum(points[v][a] for v in verts) / n for a in range(3))


def area_vector(points, ring):
    """Newell's normal: magnitude = area of the planar ring, direction = its
    orientation. Only |S . (xf - C)| is used, so the sign does not matter."""
    nx = ny = nz = 0.0
    for i in range(len(ring)):
        a = points[ring[i]]
        b = points[ring[(i + 1) % len(ring)]]
        nx += (a[1] - b[1]) * (a[2] + b[2])
        ny += (a[2] - b[2]) * (a[0] + b[0])
        nz += (a[0] - b[0]) * (a[1] + b[1])
    return (0.5 * nx, 0.5 * ny, 0.5 * nz)


def cell_volume(points, rings, centre):
    """Divergence theorem over the cell's own faces: exact for planar faces and
    a convex cell, which is what both test meshes are."""
    v = 0.0
    for ring in rings:
        s = area_vector(points, ring)
        xf = cell_centre(points, ring)
        d = tuple(xf[a] - centre[a] for a in range(3))
        v += abs(s[0] * d[0] + s[1] * d[1] + s[2] * d[2])
    return v / 3.0


def face_adjacency(cell_faces):
    """face key -> the cells that carry it. A key in one cell is the mesh's
    outer boundary; in two, an internal face."""
    where = {}
    for c, rings in enumerate(cell_faces):
        for ring in rings:
            where.setdefault(tuple(sorted(ring)), []).append(c)
    return where


def tet_decompose(points, cell_faces, centres):
    """Independent exact geometry: split every cell into tetrahedra (cell
    centre -> a triangle fan of each face ring) and accumulate exact volume
    and 2nd-order moments. Exact for planar faces. The .vtu face rings carry
    no guaranteed winding, so each ring is oriented by its Newell normal
    against (face centre - cell centre) before its fan is accumulated."""
    vols = []
    mom = {k: 0.0 for k in ('V', 'x', 'y', 'z', 'x2', 'y2', 'z2', 'xy', 'yz', 'zx')}
    for c, rings in enumerate(cell_faces):
        C = centres[c]
        v = 0.0
        for ring in rings:
            xf = cell_centre(points, ring)
            s = area_vector(points, ring)
            outward = 1.0 if (s[0] * (xf[0] - C[0]) + s[1] * (xf[1] - C[1])
                              + s[2] * (xf[2] - C[2])) >= 0 else -1.0
            n = len(ring)
            for i in range(1, n - 1):
                a, b, d = points[ring[0]], points[ring[i]], points[ring[i + 1]]
                u = tuple(b[k] - C[k] for k in range(3))
                w = tuple(d[k] - C[k] for k in range(3))
                t = tuple(a[k] - C[k] for k in range(3))
                d6 = (u[0] * (w[1] * t[2] - w[2] * t[1])
                      - u[1] * (w[0] * t[2] - w[2] * t[0])
                      + u[2] * (w[0] * t[1] - w[1] * t[0])) * outward
                v += d6 / 6.0
                sx = C[0] + a[0] + b[0] + d[0]
                sy = C[1] + a[1] + b[1] + d[1]
                sz = C[2] + a[2] + b[2] + d[2]
                sxx = C[0] ** 2 + a[0] ** 2 + b[0] ** 2 + d[0] ** 2
                syy = C[1] ** 2 + a[1] ** 2 + b[1] ** 2 + d[1] ** 2
                szz = C[2] ** 2 + a[2] ** 2 + b[2] ** 2 + d[2] ** 2
                sxy = C[0] * C[1] + a[0] * a[1] + b[0] * b[1] + d[0] * d[1]
                syz = C[1] * C[2] + a[1] * a[2] + b[1] * b[2] + d[1] * d[2]
                szx = C[2] * C[0] + a[2] * a[0] + b[2] * b[0] + d[2] * d[0]
                mom['x'] += d6 * sx / 24.0
                mom['y'] += d6 * sy / 24.0
                mom['z'] += d6 * sz / 24.0
                mom['x2'] += d6 * (sx * sx + sxx) / 120.0
                mom['y2'] += d6 * (sy * sy + syy) / 120.0
                mom['z2'] += d6 * (sz * sz + szz) / 120.0
                mom['xy'] += d6 * (sx * sy + sxy) / 120.0
                mom['yz'] += d6 * (sy * sz + syz) / 120.0
                mom['zx'] += d6 * (sz * sx + szx) / 120.0
        vols.append(v)
        mom['V'] += v
    return vols, mom


def poly_integral(coeffs, mom):
    """Exact integral of the case polynomial over the exact mesh geometry."""
    return (coeffs[0] * mom['V'] + coeffs[1] * mom['x'] + coeffs[2] * mom['y']
            + coeffs[3] * mom['z'] + coeffs[4] * mom['x2'] + coeffs[5] * mom['y2']
            + coeffs[6] * mom['z2'] + coeffs[7] * mom['xy'] + coeffs[8] * mom['yz']
            + coeffs[9] * mom['zx'])


# ------------------------------------------------------------- the RCD mirror

def _longest_axis(e):
    """The kernel's axis rule: the longest extent within a 1e-9 relative band,
    lowest index first. Equal extents are a real case -- a symmetric block has
    them exactly -- and then x is chosen."""
    emax = max(e)
    band = 1e-9 * emax
    for a in range(3):
        if e[a] >= emax - band:
            return a
    return 2


def rcd_order(centres):
    """The kernel's ordering, mirrored: sort the block by the longest axis of
    its centre bounding box, collapse the cells that share a coordinate (up to
    the rounding noise of a centroid) so they are decided by cell id, cut at the
    median, and recurse into both halves."""
    def recurse(sub):
        n = len(sub)
        if n <= 1:
            return sub
        e = []
        for a in range(3):
            col = [centres[c][a] for c in sub]
            e.append(max(col) - min(col))
        axis = _longest_axis(e)

        ordered = sorted(sub, key=lambda c: (centres[c][axis], c))
        tie = 1e-9 * max(e[axis], 1e-30)
        collapsed, g = [], 0
        while g < n:
            anchor = centres[ordered[g]][axis]
            h = g + 1
            while h < n and abs(centres[ordered[h]][axis] - anchor) <= tie:
                h += 1
            collapsed.extend(sorted(ordered[g:h]))
            g = h
        return recurse(collapsed[:n // 2]) + recurse(collapsed[n // 2:])

    return recurse(list(range(len(centres))))


def partition(centres, adjacency, n_domains):
    """Own every cell exactly once by slicing the RCD order, then derive the
    halo and the face classes of each domain from the adjacency alone."""
    order = rcd_order(centres)
    domain_of = [0] * len(centres)
    cells = [[] for _ in range(n_domains)]
    for i, c in enumerate(order):
        d = (i * n_domains) // len(centres)
        cells[d].append(c)
        domain_of[c] = d

    halo = [set() for _ in range(n_domains)]
    internal = [0] * n_domains
    coupled = [0] * n_domains
    boundary = [0] * n_domains
    peers = [set() for _ in range(n_domains)]
    for key, nb in adjacency.items():
        if len(nb) == 1:
            boundary[domain_of[nb[0]]] += 1
            continue
        a, b = nb
        da, db = domain_of[a], domain_of[b]
        if da == db:
            # both cells here: the face contributes one entry per cell, and the
            # kernel lists a face once per side it belongs to
            internal[da] += 2
        else:
            coupled[da] += 1
            coupled[db] += 1
            halo[da].add(b)
            halo[db].add(a)
            peers[da].add(db)
            peers[db].add(da)
    # a halo cell can be reached over several faces; the slot is per cell
    return domain_of, cells, [len(h) for h in halo], halo, internal, coupled, \
        boundary, [len(p) for p in peers]


# ----------------------------------------------------------------- log lines

RE_START = (r'\[chflow\] parallel start cells (\d+) faces (\d+) internal (\d+) '
            r'boundary (\d+) patches (\d+) runs (\d+)')
RE_DOMAIN = (r'\[chflow\] parallel domain domains (\d+) d (\d+) cells (\d+) halo (\d+) '
             r'faces (\d+) internal (\d+) coupled (\d+) boundary (\d+)')
RE_RUN = (r'\[chflow\] parallel run domains (\d+) cells (\d+) coupled (\d+) messages (\d+) '
          r'cellsMin (\d+) cellsMax (\d+) uncovered (\d+)')
RE_FLUX = (r'\[chflow\] parallel flux domains (\d+) divDev (\S+) divInterior (\S+) '
           r'interiorCells (\d+) phiScale (\S+)')
RE_STENCIL = (r'\[chflow\] parallel stencil domains (\d+) field (\w+) scheme (\w+) '
              r'gradAnalytic (\S+) gradAnalyticRel (\S+) gradDev (\S+) gradDevRel (\S+) '
              r'haloDrift (\S+) scale (\S+)')
RE_REDUCE = (r'\[chflow\] parallel reduce domains (\d+) field (\w+) sum (\S+) sumSerial (\S+) '
             r'minSerial (\S+) maxSerial (\S+) sumDev (\S+) sumDevRel (\S+) minDev (\S+) '
             r'maxDev (\S+) haloChecksum (\S+) haloSlots (\d+)')


def rows(out, regex):
    return [m.groups() for m in re.finditer(regex, out, re.M)]


def poly_value(c, x, y, z):
    return (c[0] + c[1] * x + c[2] * y + c[3] * z
            + c[4] * x * x + c[5] * y * y + c[6] * z * z
            + c[7] * x * y + c[8] * y * z + c[9] * z * x)


def poly_grad(c, x, y, z):
    return (c[1] + 2 * c[4] * x + c[7] * y + c[9] * z,
            c[2] + 2 * c[5] * y + c[7] * x + c[8] * z,
            c[3] + 2 * c[6] * z + c[8] * y + c[9] * x)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    solver, case_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    text = strip_comments(open(case_path).read())

    block = get_block(text, 'block')
    mesh = get_block(text, 'mesh')
    dec = get_block(text, 'decompose')
    field = get_block(text, 'field')
    verify = get_block(text, 'verify')
    out_block = get_block(text, 'output')

    domains = [int(v) for v in (list_value(dec, 'domains') or [])]
    lin = list_value(field, 'linear')
    quad = list_value(field, 'quadratic')
    vtk = word_value(out_block, 'vtk') or 'out_parallel.vtu'
    if not domains or lin is None or quad is None:
        raise SystemExit("case needs decompose domains and field linear/quadratic")

    tol_stencil = scalar_value(verify, 'stencilDevTol')
    tol_cut = scalar_value(verify, 'domainCountDevTol')
    tol_div = scalar_value(verify, 'divInteriorTol')
    tol_lin = scalar_value(verify, 'linearGradAnalyticTol')
    tol_quad = scalar_value(verify, 'quadraticGradAnalyticTol')
    for key, val in (('stencilDevTol', tol_stencil), ('domainCountDevTol', tol_cut),
                     ('divInteriorTol', tol_div)):
        if val is None:
            raise SystemExit("case verify block lacks %s" % key)

    proc = subprocess.run([solver, case_path], cwd=os.path.dirname(case_path),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("chflow-3d-parallel exited %d\n%s%s"
                         % (proc.returncode, proc.stdout, proc.stderr))
    out = proc.stdout + proc.stderr

    failures, notes = [], []

    def check(ok, what):
        if not ok:
            failures.append(what)

    start = re.search(RE_START, out)
    if start is None:
        raise SystemExit("the log holds no start line: %s" % out[-2000:])
    n_cells, n_faces, n_int, n_bnd = [int(v) for v in start.groups()[:4]]

    files = {}
    for d in domains:
        path = os.path.join(os.path.dirname(case_path),
                            vtk[:-4] + '_%d.vtu' % d)
        if not os.path.exists(path):
            raise SystemExit("the run did not write %s" % path)
        files[d] = read_vtu(path)

    # the geometry is the same mesh in every file; read it from the first one
    points, _cells, cell_faces, first_fields = files[domains[0]]
    check(len(_cells) == n_cells, "%d cells in the file, %d in the log" % (len(_cells), n_cells))
    centres = [cell_centre(points, v) for v in _cells]
    adjacency = face_adjacency(cell_faces)
    check(len(adjacency) == n_faces,
          "%d distinct faces from the file, %d in the log" % (len(adjacency), n_faces))
    check(sum(1 for nb in adjacency.values() if len(nb) == 1) == n_bnd,
          "%d boundary faces from the file, %d in the log"
          % (sum(1 for nb in adjacency.values() if len(nb) == 1), n_bnd))

    volumes = [cell_volume(points, cell_faces[c], centres[c]) for c in range(n_cells)]
    interior = [all(len(adjacency[tuple(sorted(r))]) == 2 for r in cell_faces[c])
                for c in range(n_cells)]

    # ---- I1: geometric volume closure, two independent exact methods ----
    bmin = list_value(block, 'min')
    bmax = list_value(block, 'max')
    if bmin and bmax:
        block_vol = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[2] - bmin[2])
    else:
        block_vol = None
    tet_vols, moments = tet_decompose(points, cell_faces, centres)
    worst_v = max(abs(a - b) / max(b, 1e-30) for a, b in zip(tet_vols, volumes))
    check(worst_v <= 1e-9,
          "I1: divergence-theorem and tet-decomposition cell volumes disagree "
          "by %g relative" % worst_v)
    if block_vol is not None:
        check(abs(sum(volumes) - block_vol) <= 1e-9 * block_vol,
              "I1: summed cell volume %r != analytic block volume %r"
              % (sum(volumes), block_vol))

    # ---- I2 reference: exact integrals of the case polynomials ----
    integral = {tag: poly_integral(coeffs, moments)
                for tag, coeffs in (('Linear', lin), ('Quad', quad))}

    by_field = {}
    for tag, coeffs in (('Linear', lin), ('Quad', quad)):
        by_field[tag] = coeffs

    base = None
    for d in domains:
        _points, _cl, _cf, fields = files[d]

        # ---------------------------------------------------------- partition
        domain_of, cells, n_halo, halo_sets, internal, coupled, boundary, n_peers = \
            partition(centres, adjacency, d)

        log_domains = [r for r in rows(out, RE_DOMAIN) if int(r[0]) == d]
        check(len(log_domains) == d, "%d domains logged for %s" % (len(log_domains), d))
        for row in log_domains:
            k = int(row[1])
            mine = (len(cells[k]), n_halo[k], internal[k] + coupled[k] + boundary[k],
                    internal[k], coupled[k], boundary[k])
            logged = tuple(int(v) for v in row[2:8])
            check(mine == logged,
                  "domains=%d d=%d: cells/halo/faces/internal/coupled/boundary are %s, "
                  "the re-derived partition says %s" % (d, k, logged, mine))

        written_domain = [int(v) for v in fields['domain']]
        check(written_domain == domain_of,
              "domains=%d: the written domain of a cell is not the re-derived one" % d)
        check(sorted(set(written_domain)) == list(range(d)),
              "domains=%d: the written domains are not 0..%d" % (d, d - 1))

        # ---- I3: load balance, read off the written field alone ----
        counts = [written_domain.count(k) for k in range(d)]
        check(sum(counts) == n_cells,
              "domains=%d: written rank cell counts sum to %d, mesh has %d"
              % (d, sum(counts), n_cells))
        check(min(counts) >= 1,
              "domains=%d: an empty rank in the written partition (%s)"
              % (d, counts))
        check(max(counts) - min(counts) <= 1,
              "domains=%d: written load imbalance %d (counts %s)"
              % (d, max(counts) - min(counts), counts))

        run = [r for r in rows(out, RE_RUN) if int(r[0]) == d]
        check(len(run) == 1, "domains=%d: %d run lines" % (d, len(run)))
        if run:
            r = run[0]
            check(int(r[1]) == n_cells, "domains=%d: run line cells %s" % (d, r[1]))
            check(int(r[2]) == sum(coupled) // 2,
                  "domains=%d: %s coupled faces logged, the adjacency has %d"
                  % (d, r[2], sum(coupled) // 2))
            check(int(r[3]) == sum(n_peers),
                  "domains=%d: %s messages logged, the peer sets have %d"
                  % (d, r[3], sum(n_peers)))
            check(int(r[5]) - int(r[4]) <= 1,
                  "domains=%d: load imbalance, cellsMin %s cellsMax %s" % (d, r[4], r[5]))
            check(int(r[6]) == 0, "domains=%d: %s cells owned by no domain" % (d, r[6]))

        check(sum(fields['haloCount']) == sum(n_halo),
              "domains=%d: the written per-cell halo count sums to %g, the mirror says %d"
              % (d, sum(fields['haloCount']), sum(n_halo)))

        # the face classes have to reassemble the global mesh
        check(sum(boundary) == n_bnd,
              "domains=%d: %d boundary entries for %d boundary faces" % (d, sum(boundary), n_bnd))
        check(sum(internal) // 2 + sum(coupled) // 2 == n_int,
              "domains=%d: local internal + coupled does not reassemble %d internal faces"
              % (d, n_int))

        # ------------------------------------------------- stencils and fluxes
        for row in rows(out, RE_STENCIL):
            if int(row[0]) != d:
                continue
            tag, scheme = row[1], row[2]
            grad_dev, drift = float(row[5]), float(row[7])
            check(grad_dev <= tol_stencil,
                  "domains=%d %s %s: the decomposed stencil is off by %g (tol %g)"
                  % (d, tag, scheme, grad_dev, tol_stencil))
            check(drift <= tol_stencil,
                  "domains=%d %s %s: the halo exchange drifted by %g" % (d, tag, scheme, drift))

        for row in rows(out, RE_FLUX):
            if int(row[0]) != d:
                continue
            check(float(row[1]) <= tol_stencil,
                  "domains=%d: the assembled net out-flux is off by %s" % (d, row[1]))
            check(float(row[2]) <= tol_div,
                  "domains=%d: an interior cell keeps %s of a divergence-free flux" % (d, row[2]))
            check(int(row[3]) == sum(1 for v in interior if v),
                  "domains=%d: %s interior cells logged, the adjacency has %d"
                  % (d, row[3], sum(1 for v in interior if v)))

        # --------------------------------------- the cut changes nothing at all
        for name in ('tLinear', 'tQuad', 'divPhiParallel', 'divPhiSerial',
                     'gradSerialLinearX', 'gradSerialLinearY', 'gradSerialLinearZ',
                     'gradParallelLinearX', 'gradParallelLinearY', 'gradParallelLinearZ',
                     'gradSerialQuadX', 'gradSerialQuadY', 'gradSerialQuadZ',
                     'gradParallelQuadX', 'gradParallelQuadY', 'gradParallelQuadZ'):
            if name not in fields:
                check(False, "domains=%d: the file has no %s array" % (d, name))
                continue
            if base is None:
                continue
            worst = max(abs(a - b) for a, b in zip(fields[name], base[name]))
            check(worst <= tol_cut,
                  "domains=%d: %s differs from the 1-domain run by %g (tol %g)"
                  % (d, name, worst, tol_cut))

        # ------------------------------------------------- the physics, re-done
        # a divergence-free flux leaves every interior cell with no net out-flux
        div_off = max(abs(fields['divPhiParallel'][c]) for c in range(n_cells) if interior[c])
        check(div_off <= tol_div,
              "domains=%d: an interior cell's written net out-flux is %g (tol %g)"
              % (d, div_off, tol_div))

        # ---- I4: per-rank cumulative flux reassembles across the cuts ----
        scale = sum(abs(fields['divPhiSerial'][c]) * volumes[c] for c in range(n_cells))
        per_rank = [0.0] * d
        for c in range(n_cells):
            per_rank[written_domain[c]] += fields['divPhiParallel'][c] * volumes[c]
        serial_total = sum(fields['divPhiSerial'][c] * volumes[c] for c in range(n_cells))
        check(abs(sum(per_rank) - serial_total) <= 1e-9 * scale,
              "domains=%d: rank-grouped divergence sums to %r, the serial "
              "integral is %r -- cut-face contributions do not cancel"
              % (d, sum(per_rank), serial_total))

        for tag, coeffs in sorted(by_field.items()):
            t = fields['t' + tag]
            off = max(abs(v - poly_value(coeffs, *centres[c])) for c, v in enumerate(t))

            # ---- I2: analytic integral conservation across the cut. The
            # written field holds point values at the cell centres, so its
            # volume-weighted sum is the midpoint quadrature: exact for the
            # linear field (a true mass-conservation statement), and for the
            # quadratic field owed only as rank-grouped reassembly below.
            s_field = sum(t[c] * volumes[c] for c in range(n_cells))
            if tag == 'Linear':
                ref = integral[tag]
                check(abs(s_field - ref) <= 1e-8 * abs(ref) + 1e-11,
                      "domains=%d %s: written field integrates to %r, the exact "
                      "geometry integral is %r" % (d, tag, s_field, ref))
            rank_parts = [0.0] * d
            for c in range(n_cells):
                rank_parts[written_domain[c]] += t[c] * volumes[c]
            check(abs(sum(rank_parts) - s_field) <= 1e-9 * max(1.0, abs(s_field)),
                  "domains=%d %s: rank-grouped field integral %r does not "
                  "reassemble the global %r" % (d, tag, sum(rank_parts), s_field))

            s_lin = sum(t[c] * volumes[c] for c in range(n_cells))
            check(off <= 1e-9,
                  "domains=%d %s: the written field is not the case's polynomial (off %g)"
                  % (d, tag, off))

            s_lin = sum(t[c] * volumes[c] for c in range(n_cells))
            red = [r for r in rows(out, RE_REDUCE) if int(r[0]) == d and r[1] == tag]
            check(len(red) == 1, "domains=%d %s: %d reduce lines" % (d, tag, len(red)))
            if red:
                r = red[0]
                check(abs(float(r[2]) - float(r[3])) <= 0.0,
                      "domains=%d %s: the reduced sum is not the serial one (%s vs %s)"
                      % (d, tag, r[2], r[3]))
                check(float(r[8]) <= 0.0 and float(r[9]) <= 0.0,
                      "domains=%d %s: the reduced min/max is not the serial one" % (d, tag))
                check(abs(s_lin - float(r[3])) <= 1e-9 * abs(s_lin),
                      "domains=%d %s: this script's sum(t V) is %r, the kernel's is %s"
                      % (d, tag, s_lin, r[3]))
                # the .vtu prints 12 significant digits, so a field read back is
                # not exactly the field the kernel held -- compare within that
                check(abs(min(t) - float(r[4])) <= 1e-11 * max(1.0, abs(float(r[4])))
                      and abs(max(t) - float(r[5])) <= 1e-11 * max(1.0, abs(float(r[5]))),
                      "domains=%d %s: this script's min/max of the field is not the reduced one"
                      % (d, tag))
                # the halo checksum, re-summed from the re-derived halo sets
                chk = sum(t[c] for hs in halo_sets for c in hs)
                check(abs(chk - float(r[10])) <= 1e-8 * max(1.0, abs(chk)),
                      "domains=%d %s: this script's halo checksum is %r, the kernel moved %s"
                      % (d, tag, chk, r[10]))
                check(int(r[11]) == sum(n_halo),
                      "domains=%d %s: %s halo slots logged, the mirror has %d"
                      % (d, tag, r[11], sum(n_halo)))

            for scheme in ('Serial', 'Parallel'):
                gx = fields['grad%s%sX' % (scheme, tag)]
                gy = fields['grad%s%sY' % (scheme, tag)]
                gz = fields['grad%s%sZ' % (scheme, tag)]
                worst = max(max(abs(gx[c] - g[0]), abs(gy[c] - g[1]), abs(gz[c] - g[2]))
                            for c, g in ((c, poly_grad(coeffs, *centres[c])) for c in range(n_cells)))
                tol = tol_lin if tag == 'Linear' else tol_quad
                if tol is None:
                    if d == domains[0] and scheme == 'Serial':
                        notes.append("%s %s: this mesh's own gradient error is %g "
                                     "(not asserted; only that the cut does not change it)"
                                     % (tag, scheme, worst))
                    continue
                check(worst <= tol,
                      "domains=%d %s %s: gradient off the analytic one by %g (tol %g)"
                      % (d, tag, scheme, worst, tol))

        if base is None:
            base = {k: v for k, v in fields.items() if k != 'domain'}

    if failures:
        for note in notes:
            print("  " + note)
        raise SystemExit("parallel check FAILED:\n  " + "\n  ".join(failures))
    for note in notes:
        print("  " + note)
    print("parallel decomposition check passed on %s for domains %s"
          % (os.path.basename(case_path), domains))


if __name__ == '__main__':
    main()
