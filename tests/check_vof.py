#!/usr/bin/env python3
"""Run a prescribed-divergence-free-velocity VOF case and verify what the phase
fraction and the interface geometry have to deliver.

The two cases are the two things the increment claims:

  vofTranslation - a plane interface carried by a uniform velocity. The
      transport's contract is conservation and boundedness, and here both have a
      closed form: the phase volume grows by exactly the volume the inlet feeds,
      U * area * t, and alpha never leaves [0,1] because upwind convection of a
      divergence-free flux is an M-matrix. The plane is also the case where the
      curvature operator has to read zero, whatever the mesh.

  vofRotation - a droplet carried around a solid-body rotation. The exact
      solution is the initial field turned through omega * t, so the phase
      centroid has to sit at angle omega * t on its starting radius, the phase
      volume has to be conserved in a closed box, and the interface has to still
      read kappa = -1/R with nhat radial. This is the case that measures the
      geometry, which a plane cannot (its normal never varies).

Both velocities are linear in the coordinates, which is why sum_f U(x_f).Sf is
exactly zero per cell and the checks below are absolute rather than grid-rate.

usage: check_vof.py <chflow-vof> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile


def case_constants(path):
    """The numbers the checks are built from, read out of the case itself so it
    cannot drift away from them."""
    txt = open(path).read()

    def triple(key):
        return [float(v) for v in re.search(r'%s\s*\(([^)]*)\)' % key, txt).group(1).split()]

    lo, hi = triple('min'), triple('max')
    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]
    td = re.search(r'time\s*\{([^}]*)\}', txt, re.S).group(1)
    vel = re.search(r'velocity\s*\{([^}]*)\}', txt, re.S).group(1)
    al = re.search(r'alpha\s*\{(.*?)\n\}', txt, re.S).group(1)

    c = dict(
        lo=lo, hi=hi, cells=cells,
        dx=[(hi[i] - lo[i]) / cells[i] for i in range(3)],
        vol=[hi[i] - lo[i] for i in range(3)],
        dt=float(re.search(r'dt\s+([0-9.eE+-]+);', td).group(1)),
        end_time=float(re.search(r'endTime\s+([0-9.eE+-]+);', td).group(1)),
        vtu=re.search(r'vtk\s+([A-Za-z0-9_.-]+);', txt).group(1),
    )
    init = re.search(r'vtkInitial\s+([A-Za-z0-9_.-]+);', txt)
    c['vtu_initial'] = init.group(1) if init else None
    c['steps'] = int(round(c['end_time'] / c['dt']))

    if 'uniform' in vel:
        c['kind'] = 'translation'
        c['U'] = [float(v) for v in
                  re.search(r'U\s*\(([^)]*)\)', vel).group(1).split()]
        c['area'] = c['vol'][1] * c['vol'][2]      # the inlet is the xmin face
    elif 'rotation' in vel:
        c['kind'] = 'rotation'
        c['omega'] = float(re.search(r'omega\s+([0-9.eE+-]+);', vel).group(1))
        c['axis'] = [float(v) for v in
                     re.search(r'center\s*\(([^)]*)\)', vel).group(1).split()]
    else:
        raise SystemExit("the velocity block must be uniform or rotation")

    if 'plane' in al:
        c['shape'] = 'plane'
        c['point'] = [float(v) for v in
                      re.search(r'point\s*\(([^)]*)\)', al).group(1).split()]
        c['normal'] = [float(v) for v in
                       re.search(r'normal\s*\(([^)]*)\)', al).group(1).split()]
    elif 'circle' in al:
        c['shape'] = 'circle'
        c['center'] = [float(v) for v in
                       re.search(r'center\s*\(([^)]*)\)', al).group(1).split()]
        c['radius'] = float(re.search(r'radius\s+([0-9.eE+-]+);', al).group(1))
        c['smooth'] = float(re.search(r'smooth\s+([0-9.eE+-]+);', al).group(1))
        if c['kind'] == 'rotation':
            c['r0'] = math.hypot(c['center'][0] - c['axis'][0],
                                 c['center'][1] - c['axis'][1])
    else:
        raise SystemExit("the alpha block must be a plane or a circle")
    return c


# ---------------------------------------------------------------- vtu reading
# The geometry is checked on what the solver wrote, not on what it printed, so
# the cell centres it measured are the mesh's own.

def _arrays(txt):
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', txt, re.S):
        name = re.search(r'Name="([^"]+)"', m.group(1))
        if not name:
            continue
        out[name.group(1)] = [float(v) for v in m.group(2).split()]
    return out


def read_vtu(path):
    txt = open(path).read()
    head = re.search(r'<DataArray type="Float64" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
                     txt, re.S)
    p = [float(v) for v in head.group(1).split()]
    points = list(zip(p[0::3], p[1::3], p[2::3]))
    a = _arrays(txt)
    conn, off = [int(v) for v in a['connectivity']], [int(v) for v in a['offsets']]
    cent, s = [], 0
    for o in off:
        cl = [points[v] for v in conn[s:o]]
        s = o
        cent.append(tuple(sum(q[i] for q in cl) / len(cl) for i in range(3)))
    nf = len(cent)
    fields = {}
    for name, v in a.items():
        if name in ('connectivity', 'offsets', 'types') or len(v) % nf:
            continue
        if len(v) == nf:
            fields[name] = list(v)
        else:
            k = len(v) // nf
            fields[name] = [tuple(v[i * k:(i + 1) * k]) for i in range(nf)]
    return cent, fields


# ---------------------------------------------------------------- log reading

def log_numbers(out):
    """Every report line as a dict, plus the two lines that are not per-step."""
    keys, rows = None, []
    for m in re.finditer(r'^\[chflow\] vof t (.*?)(?: solver \d+ residual [0-9.eE+-]+)?$',
                         out, re.M):
        v = m.group(1).split()
        if keys is None:
            keys = ['t'] + v[1::2]
        rows.append(dict(zip(keys, [float(x) for x in v[0::2]])))
    start = re.search(r'\[chflow\] vof start: mass ([0-9.eE+-]+) dt ([0-9.eE+-]+) '
                      r'endTime ([0-9.eE+-]+) kappa ([0-9.eE+-]+) iface (\d+)', out)
    div = re.search(r'\[chflow\] vof velocity: (\w+) max \|div\(phi\)\| ([0-9.eE+-]+)', out)
    if not rows or start is None or div is None:
        raise SystemExit("the log does not hold a start line and the report lines "
                         "the checker needs")
    return dict(rows=rows, mass0=float(start.group(1)), kappa0=float(start.group(4)),
                iface0=int(start.group(5)), kind=div.group(1),
                max_div=float(div.group(2)))


def failures_exit(case, failures, notes):
    for n in notes:
        print("  " + n)
    if failures:
        raise SystemExit("VOF check FAILED:\n  " + "\n  ".join(failures))
    print("VOF %s check passed" % case['kind'])


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    solver, case_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    c = case_constants(case_path)
    work = tempfile.mkdtemp(prefix='chflow_vof_')
    try:
        shutil.copy(os.path.abspath(case_path), os.path.join(work, 'case.dict'))
        proc = subprocess.run([solver, 'case.dict'], cwd=work,
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit("chflow-vof exited %d\n%s" % (proc.returncode, proc.stderr))
        out = proc.stdout + proc.stderr
        log = log_numbers(out)
        cent, f = read_vtu(os.path.join(work, c['vtu']))
        alpha, nh, kap = f['alpha'], f['nHat'], f['kappa']
        init = None
        if c['vtu_initial']:
            init = read_vtu(os.path.join(work, c['vtu_initial']))
        band = [i for i, v in enumerate(nh) if math.hypot(v[0], v[1], v[2]) > 0.5]
        if not band:
            raise SystemExit("no cell carries an interface in the result")
        failures, notes = [], []
        check_common(c, log, cent, alpha, nh, kap, band, failures, notes)
        if c['kind'] == 'translation':
            check_translation(c, log, cent, alpha, nh, kap, band, failures, notes)
        else:
            check_rotation(c, log, cent, alpha, nh, kap, band, failures, notes)
        if init:
            check_initial_geometry(c, init, failures, notes)
        failures_exit(c, failures, notes)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def check_common(c, log, cent, alpha, nh, kap, band, failures, notes):
    # the prescribed field has to be discretely divergence-free, or nothing
    # below means anything: conservation is a property of the transport on such
    # a flux, not of the transport alone
    exact = c['kind'] == 'translation'
    limit = 1e-18
    notes.append("max |div(phi)| = %g for the %s velocity (%s)"
                 % (log['max_div'], log['kind'],
                    "the Gauss identity, so exactly zero is the expectation"
                    if exact else "the same identity, at round-off"))
    if log['max_div'] > limit:
        failures.append("the prescribed velocity is not discretely divergence-free "
                        "(max |div(phi)| = %g)" % log['max_div'])

    # boundedness: the M-matrix claim, with no clipping anywhere in the kernel
    lo, hi = min(alpha), max(alpha)
    notes.append("alpha over %d cells: [%.6g, %.6g]" % (len(alpha), lo, hi))
    if lo < -1e-12 or hi > 1.0 + 1e-12:
        failures.append("alpha left [0,1]: [%.6g, %.6g]" % (lo, hi))

    # the transport itself has to be solved and every report has to be finite
    for r in log['rows']:
        for k, v in r.items():
            if k != 't' and not math.isfinite(v):
                failures.append("the report at t=%g has a non-finite %s" % (r['t'], k))
    if len(log['rows']) < 2:
        failures.append("too few reports to follow the interface")


def crossing_x50(cent, alpha, ref_y, ref_z):
    """Where alpha crosses 0.5 along x on the row of cells nearest (ref_y,
    ref_z) -- the position of a plane interface."""
    best, row_ref = None, None
    for i in range(len(cent)):
        d = abs(cent[i][1] - ref_y) + abs(cent[i][2] - ref_z)
        if best is None or d < best:
            best, row_ref = d, cent[i]
    row = [i for i in range(len(cent))
           if abs(cent[i][1] - row_ref[1]) < 1e-9 and abs(cent[i][2] - row_ref[2]) < 1e-9]
    row.sort(key=lambda i: cent[i][0])
    for a, b in zip(row, row[1:]):
        a0, a1 = alpha[a], alpha[b]
        if (a0 - 0.5) * (a1 - 0.5) <= 0.0 and a0 != a1:
            x0, x1 = cent[a][0], cent[b][0]
            return x0 + (x1 - x0) * (0.5 - a0) / (a1 - a0)
    return float('nan')


def check_translation(c, log, cent, alpha, nh, kap, band, failures, notes):
    u = c['U'][0]
    if abs(c['U'][1]) or abs(c['U'][2]) or abs(c['normal'][1]) or abs(c['normal'][2]):
        raise SystemExit("this checker follows an interface along x, so the "
                         "velocity and the plane normal have to be +x")

    # conservation, in closed form: phase 1 enters at the inlet area times U and
    # the only way the domain can hold less is the phase the outlet carries away.
    # So the deficit has a bound, mass0 + U*A*t - mass(t) <= U*A*t*alpha(outlet),
    # and it has a sign. Both are what the M-matrix transport owes.
    area, m0 = c['area'], log['mass0']
    last = log['rows'][-1]
    fed = u * area * last['t']
    a_out = max(alpha[i] for i in range(len(cent))
                if abs(cent[i][0] - (c['lo'][0] + c['vol'][0] - c['dx'][0] / 2.0)) < 1e-9)
    deficit = m0 + fed - last['mass']
    bound = fed * a_out + 1e-12
    notes.append("phase volume %.9g after t=%g; the inlet fed %.9g, so the balance "
                 "is short by %.3g -- the outlet's own alpha is %.3g and that allows "
                 "%.3g" % (last['mass'], last['t'], fed, deficit, a_out, bound))
    if deficit < -1e-12:
        failures.append("the domain holds more phase than the inlet fed (%.3g over)"
                        % (-deficit))
    if deficit > bound:
        failures.append("the transport lost %.3g of phase volume, over the %.3g the "
                        "outlet can account for" % (deficit, bound))
    if deficit / m0 > 1e-5:
        failures.append("the balance error is %.3g relative, not the leakage the "
                        "front implies" % (deficit / m0))

    # the interface has to move with the flow. Its absolute position carries the
    # bias of reading a 0.5 crossing off cell centres, so the speed over the run
    # is the number that means advection.
    x_row = crossing_x50(cent, alpha, last['cy'], c['lo'][2] + c['dx'][2] / 2.0)
    rows = sorted(log['rows'], key=lambda r: r['t'])
    if not math.isfinite(x_row) or abs(x_row - last['x50']) > 1e-9:
        failures.append("the 0.5 crossing read from the result (%s) is not the one "
                        "the run reported (%s)" % (x_row, last['x50']))
    speed = (last['x50'] - rows[0]['x50']) / (last['t'] - rows[0]['t'])
    bias = last['x50'] - (c['point'][0] + u * last['t'])
    notes.append("0.5 crossing %.6f -> %.6f over t %.4f -> %.4f: speed %.6f for U=%g "
                 "(%.3f%%), position %.4f dx behind x0+U*t"
                 % (rows[0]['x50'], last['x50'], rows[0]['t'], last['t'], speed, u,
                    100.0 * (speed / u - 1.0), -bias / c['dx'][0]))
    if abs(speed / u - 1.0) > 0.01:
        failures.append("the interface moves at %.6f, not U=%g" % (speed, u))
    if abs(bias) > 0.4 * c['dx'][0]:
        failures.append("the 0.5 crossing is %.4f m from x0+U*t, over 0.4 dx" % bias)

    # the field is one-dimensional by construction: with no flux across the
    # y-faces nothing can make alpha depend on y, so the only spread the rows can
    # have is the linear solver leaving each cell at its own last bit
    by = {}
    for i in range(len(cent)):
        by.setdefault((cent[i][0], cent[i][2]), []).append(alpha[i])
    worst = max(max(v) - min(v) for v in by.values())
    notes.append("alpha across the %d rows of y differs by at most %g (%d columns "
                 "compared; the ulp of alpha=1 is %.3g)"
                 % (c['cells'][1], worst, len(by), 2.0 ** -52))
    if worst > 1e-12:
        failures.append("the plane interface is not y-invariant (spread %g)" % worst)

    # monotone across the interface, in every row: the phase is on the negative
    # side of the normal, so alpha can only fall as x grows
    rows = {}
    for i in range(len(cent)):
        rows.setdefault(cent[i][1], []).append(i)
    rises = 0
    for v in rows.values():
        v.sort(key=lambda i: cent[i][0])
        rises += sum(1 for a, b in zip(v, v[1:]) if alpha[b] > alpha[a] + 1e-14)
    notes.append("across x: %d rows of %d cells, %d rises in alpha (the interface "
                 "falls towards +x)" % (len(rows), c['cells'][0], rises))
    if rises:
        failures.append("alpha is not monotone across the interface (%d rises)" % rises)

    # a plane has no curvature and its normal has no direction to get wrong
    kmax = max(abs(v) for v in kap)
    dev = max(abs(nh[i][1]) + abs(nh[i][2]) for i in band)
    nmax = max(nh[i][0] for i in band)
    notes.append("kappa: max |kappa| = %g over the mesh; nhat is off the x axis by at "
                 "most %g and its largest x component is %.6f" % (kmax, dev, nmax))
    if kmax > 1e-8:
        failures.append("a plane interface reads |kappa| up to %g" % kmax)
    if dev > 1e-12:
        failures.append("the plane's normal is off the x axis by %.3g" % dev)
    if nmax > -0.999:
        failures.append("the plane's normal does not point towards the phase it "
                        "was initialised on (nhat_x up to %.6f)" % nmax)


def check_rotation(c, log, cent, alpha, nh, kap, band, failures, notes):
    om, t = c['omega'], log['rows'][-1]['t']

    # closed box, droplet nowhere near a wall: the phase volume can only change
    # by what the boundary exchanges
    rel = abs(log['rows'][-1]['mass'] / log['mass0'] - 1.0)
    notes.append("phase volume %.9g -> %.9g over t=%g: %.3e relative"
                 % (log['mass0'], log['rows'][-1]['mass'], t, rel))
    if rel > 1e-8:
        failures.append("the phase volume is not conserved in the closed box (%.3e)" % rel)

    # rigid rotation: the centroid goes round the axis at the radius the initial
    # field put it on
    ax, r0 = c['axis'], c['r0']
    ang = math.atan2(log['rows'][-1]['cy'] - ax[1], log['rows'][-1]['cx'] - ax[0])
    rad = math.hypot(log['rows'][-1]['cx'] - ax[0], log['rows'][-1]['cy'] - ax[1])
    # wrapped, so a case that turns past pi still reads the error it has
    turn = (ang - om * t + math.pi) % (2.0 * math.pi) - math.pi
    notes.append("centroid after t=%g: %.6f rad turned past omega*t=%.6f (%.3f%% of "
                 "it), radius %.6f for %.6f (%.3f%%)"
                 % (t, turn, om * t, 100.0 * turn / (om * t), rad, r0,
                    100.0 * (rad / r0 - 1.0)))
    if abs(turn / (om * t)) > 2e-3:
        failures.append("the droplet is %.5f rad off the angle omega*t=%.5g it should "
                        "have turned" % (turn, om * t))
    if abs(rad / r0 - 1.0) > 5e-3:
        failures.append("the droplet slid %.3f%% in radius: it is not on its circle "
                        "any more" % (100.0 * (rad / r0 - 1.0)))

    # the geometry that survives advection: the band's mean curvature is still
    # the curvature of the droplet. The cell by cell numbers are asserted against
    # the phase's own centroid and with a wider band than the initial field's,
    # because upwind diffusion is anisotropic -- |U| = omega r stretches the band
    # along the flow and tilts its normals -- which is the transport's doing, not
    # the operators'. The analytic initial field is measured in
    # check_initial_geometry, where the band's own geometry is all that shows.
    kmean = sum(kap[i] for i in band) / len(band)
    notes.append("%d interface cells after advection: band mean kappa %.4f for "
                 "-1/R=%.4f (%d cells at t=0 read %.4f)"
                 % (len(band), kmean, -1.0 / c['radius'], log['iface0'], log['kappa0']))
    if abs(kmean / (-1.0 / c['radius']) - 1.0) > 0.08:
        failures.append("the band's mean curvature after advection is %.4f, not "
                        "-1/R=%.4f" % (kmean, -1.0 / c['radius']))

    # the droplet is still a droplet: its normals still point at its own centre.
    # Measuring them against the centre the case started from would charge the
    # 0.4 rad the whole field has turned to the interface instead of to the flow.
    mass = sum(alpha)
    cx = sum(a * p[0] for a, p in zip(alpha, cent)) / mass
    cy = sum(a * p[1] for a, p in zip(alpha, cent)) / mass
    dev, kErr = [], []
    for i in band:
        rx, ry = cent[i][0] - cx, cent[i][1] - cy
        r = math.hypot(rx, ry)
        if r < 1e-9:
            continue
        n = math.hypot(nh[i][0], nh[i][1])
        # nhat points towards increasing alpha, so into the droplet: -rhat
        dev.append(math.degrees(math.acos(max(-1.0, min(1.0,
                    -(nh[i][0] * rx + nh[i][1] * ry) / (r * n))))))
        kErr.append(abs(kap[i] * r + 1.0))
    notes.append("after advection, off the direction towards the phase centroid: "
                 "%.2f deg mean, %.2f deg max over %d band cells; cell by cell "
                 "|kappa*r + 1| = %.3f mean, %.3f max"
                 % (sum(dev) / len(dev), max(dev), len(dev),
                    sum(kErr) / len(kErr), max(kErr)))
    if sum(dev) / len(dev) > 8.0 or max(dev) > 30.0:
        failures.append("the advected interface is no longer round: its normals are "
                        "%.2f deg mean / %.2f deg max off the radius"
                        % (sum(dev) / len(dev), max(dev)))


def check_initial_geometry(c, init, failures, notes):
    """The cell by cell geometry of the initial field, where the answer is
    analytic: a plane has no curvature and a circle of radius R has -1/r at the
    cell a distance r from its centre."""
    cent, f = init
    nh, kap = f['nHat'], f['kappa']
    band = [i for i, v in enumerate(nh) if math.hypot(v[0], v[1], v[2]) > 0.5]
    if not band:
        failures.append("the initial field carries no interface for the geometry "
                        "operators to measure")
        return

    if c['shape'] == 'plane':
        kmax = max(abs(v) for v in kap)
        notes.append("at t=0, %d band cells: max |kappa| = %g over the mesh"
                     % (len(band), kmax))
        if kmax > 1e-8:
            failures.append("the initial plane interface reads |kappa| up to %g" % kmax)
        return

    cx0, cy0 = c['center'][0], c['center'][1]
    ang_err, k_err = [], []
    for i in band:
        rx, ry = cent[i][0] - cx0, cent[i][1] - cy0
        r = math.hypot(rx, ry)
        if r < 1e-9:
            continue
        d = -(nh[i][0] * rx + nh[i][1] * ry) / r
        ang_err.append(math.degrees(math.acos(max(-1.0, min(1.0, d)))))
        k_err.append(abs(kap[i] * r + 1.0))
    mean_a, max_a = sum(ang_err) / len(ang_err), max(ang_err)
    mean_k, max_k = sum(k_err) / len(k_err), max(k_err)
    notes.append("at t=0, %d band cells: nhat off the radial direction by %.2f deg "
                 "mean, %.2f deg max; cell by cell |kappa*r + 1| = %.3f mean, %.3f "
                 "max (initial profile %.2f dx thick)"
                 % (len(band), mean_a, max_a, mean_k, max_k, c['smooth']))
    if mean_a > 2.5 or max_a > 6.0:
        failures.append("the initial interface normal is not radial (%.2f deg mean, "
                        "%.2f max)" % (mean_a, max_a))
    if mean_k > 0.2:
        failures.append("the initial kappa is %.3f off -1/r per cell on average" % mean_k)
    if max_k > 1.0:
        failures.append("the initial kappa reaches %.3f off -1/r in a cell" % max_k)


if __name__ == '__main__':
    main()
