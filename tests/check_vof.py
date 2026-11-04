#!/usr/bin/env python3
"""Run a VOF case and verify what the phase fraction, the interface geometry and
-- in the momentum mode -- the pressure-velocity response have to deliver.

The transport cases (chflow-vof with a velocity{} block):

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

The momentum case (chflow-vof with a surfaceTension{} block, no velocity):

  vofDroplet - a droplet held still while the coupling answers the surface force
      its interface implies. The number the discrete force is NOT given is the
      Laplace jump p_inside - p_outside = sigma/(rho R) for a 2D droplet; the
      solver has to find it by balancing -grad(p) against -(sigma/rho) kappa
      grad(alpha), which is what a collocated pressure-velocity coupling is
      tested by. The current that survives that balance is the spurious one, so
      the case reports its amplitude against the capillary-inertial and free-
      acceleration scales.

  vofOscillating - the two-way loop: the solved flux carries alpha and the
      surface force is rebuilt from the shape it leaves, so an ellipse released
      off the circle oscillates in its n = 2 mode and decays. The period is the
      number the kernel is not given, and the case exists to decide WHICH law:
      one-fluid momentum means the ambient carries half the inertia,
      w^2 = 3 sigma/(rho R^3) on R = sqrt(a b) -- not the free-drop
      6 sigma/(rho R^3), which would be 35% off.

The geometric case (chflow-vof with a reconstruction{} block):

  vofPlic - the same droplet, but filled by the disk's exact overlap area, so the
      interface is one cell thick and the field's slope cannot measure it. The
      reconstruction cuts a plane in every band cell, and the numbers it owes have
      closed forms: the plane reads its cell's own fraction back, the planes cover
      exactly the cells whose alpha is strictly between 0 and 1, the area across a
      one-cell band is the circle's extruded perimeter, and the volume of phase 1
      is pi R^2 dz to start with and conserved throughout.

  vofZalesak - the same droplet and the same reconstruction, but alpha is carried
      by the planes instead of by a flux of alpha (`reconstruction{ transport
      geometric }`), which is Zalesak's rigid-rotation test. The exact answer is
      known cell by cell at every turn: a disk turned through omega t is still the
      disk, moved, so the field the transport leaves can be compared with the
      overlap area of that disk and each cell. What the geometric path owes and the
      upwind path cannot: the phase volume to round-off, alpha inside [0,1] with
      nothing clipped, the interface one cell thick after a full turn rather than
      smeared over the front's own width, and the shape error of the disk reduced
      to a fraction of R.

  vofFit - the same droplet, the same turn and the same planes, but the plane of every
      band cell is the least-squares plane of the neighbourhood around it
      (`reconstruction{ normal fit }`) rather than the one the cell's own staircase
      reads. Everything owed above is owed here too, because it is the same cut; what
      the fit adds to the log is the distance the planes of one neighbourhood still
      disagree by and the rounds it took to agree, and the checker follows both.

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
    cannot drift away from them. The prose is stripped first: a case may describe
    a key in its header without that being what it asks the solver to do."""
    txt = "\n".join(re.sub(r'//.*', '', line) for line in open(path))

    def triple(key):
        return [float(v) for v in re.search(r'%s\s*\(([^)]*)\)' % key, txt).group(1).split()]

    lo, hi = triple('min'), triple('max')
    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]
    td = re.search(r'time\s*\{([^}]*)\}', txt, re.S).group(1)
    al = re.search(r'alpha\s*\{(.*?)\n\}', txt, re.S).group(1)

    c = dict(
        lo=lo, hi=hi, cells=cells,
        dx=[(hi[i] - lo[i]) / cells[i] for i in range(3)],
        vol=[hi[i] - lo[i] for i in range(3)],
        dt=float(re.search(r'dt\s+([0-9.eE+-]+);', td).group(1)),
        end_time=float(re.search(r'endTime\s+([0-9.eE+-]+);', td).group(1)),
        vtu=re.search(r'vtk\s+([A-Za-z0-9_.-]+);', txt).group(1),
        momentum='surfaceTension' in txt,
    )
    init = re.search(r'vtkInitial\s+([A-Za-z0-9_.-]+);', txt)
    c['vtu_initial'] = init.group(1) if init else None
    c['steps'] = int(round(c['end_time'] / c['dt']))

    if c['momentum']:
        # The checker needs the patch that carries U, so the velocity block's
        # own U() is not what is read here.
        st = re.search(r'surfaceTension\s*\{([^}]*)\}', txt, re.S).group(1)
        adv = re.search(r'advect\s+(\w+)\s*;', st)
        c['advect'] = bool(adv and adv.group(1) == 'true')
        c['kind'] = 'oscillating' if c['advect'] else 'droplet'
        c['sigma'] = float(re.search(r'sigma\s+([0-9.eE+-]+);', st).group(1))
        c['rho'] = float(re.search(r'rho\s+([0-9.eE+-]+);', st).group(1))
        c['nu'] = float(re.search(r'nu\s+([0-9.eE+-]+);', st).group(1))
        corr = re.search(r'correctors\s+([0-9.eE+-]+);', st)
        c['correctors'] = int(corr.group(1)) if corr else 2
    else:
        vel = re.search(r'velocity\s*\{([^}]*)\}', txt, re.S).group(1)
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
        sm = re.search(r'smooth\s+([0-9.eE+-]+);', al)
        c['smooth'] = float(sm.group(1)) if sm else 0.0
        fl = re.search(r'fill\s+(\w+)\s*;', al)
        c['fill'] = fl.group(1) if fl else 'cell'
        if c['kind'] == 'rotation':
            c['r0'] = math.hypot(c['center'][0] - c['axis'][0],
                                 c['center'][1] - c['axis'][1])
    elif 'ellipse' in al:
        c['shape'] = 'ellipse'
        c['center'] = [float(v) for v in
                       re.search(r'center\s*\(([^)]*)\)', al).group(1).split()]
        c['a'] = float(re.search(r'radiusA\s+([0-9.eE+-]+);', al).group(1))
        c['b'] = float(re.search(r'radiusB\s+([0-9.eE+-]+);', al).group(1))
        c['smooth'] = float(re.search(r'smooth\s+([0-9.eE+-]+);', al).group(1))
    else:
        raise SystemExit("the alpha block must be a plane, a circle or an ellipse")
    # A reconstruction{} block makes the case the geometric one: what is verified
    # is the planes cut into the band, not the smeared band's slope.
    rc = re.search(r'reconstruction\s*\{([^}]*)\}', txt, re.S)
    if rc:
        if c['momentum'] or c['shape'] != 'circle' or c.get('fill') != 'exact':
            raise SystemExit("a reconstruction{} case has to drive the transport "
                             "with a circle filled exact")
        if c['smooth'] > 0:
            raise SystemExit("an exactly filled interface has no smoothing to smear")
        nm = re.search(r'normal\s+(\w+)\s*;', rc.group(1))
        c['normal_choice'] = nm.group(1) if nm else 'mixed'
        bt = re.search(r'bandTolerance\s+([0-9.eE+-]+);', rc.group(1))
        c['band_tolerance'] = float(bt.group(1)) if bt else 1e-6
        tp = re.search(r'transport\s+(\w+)\s*;', rc.group(1))
        c['transport'] = tp.group(1) if tp else 'alpha'
        if c['transport'] == 'geometric':
            # the closed rotation is the case with an analytic answer at every
            # turn; on a plane front the geometric path would only be re-testing
            # what the reconstruction already proves
            if c['kind'] != 'rotation':
                raise SystemExit("transport geometric is verified on the rotation, "
                                 "whose exact field the checker can write down")
        elif c['transport'] != 'alpha':
            raise SystemExit("reconstruction transport must be alpha or geometric")
        c['kind'] = 'plic'
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
    """Every report line as a dict, plus the lines that are not per-step."""
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
    csf = re.search(r'\[chflow\] vof csf: sigma/rho ([0-9.eE+-]+) nu ([0-9.eE+-]+) '
                    r'advect (true|false) kappa ([0-9.eE+-]+) max \|f\| ([0-9.eE+-]+)', out)
    plic = []
    for m in re.finditer(r'^\[chflow\] plic (start|t [0-9.eE+-]+): band (\d+) cells '
                         r'area ([0-9.eE+-]+) readBack ([0-9.eE+-]+) kappa ([0-9.eE+-]+)'
                         r'(?: scatter ([0-9.eE+-]+) rounds (\d+))?',
                         out, re.M):
        plic.append(dict(stamp=m.group(1), band=int(m.group(2)), area=float(m.group(3)),
                         read_back=float(m.group(4)), kappa=float(m.group(5)),
                         scatter=None if m.group(6) is None else float(m.group(6)),
                         rounds=None if m.group(7) is None else int(m.group(7))))
    geo = []
    for m in re.finditer(r'^\[chflow\] geo t ([0-9.eE+-]+) flux (\d+) limited (\d+) '
                         r'rounds (\d+) minFactor ([0-9.eE+-]+) maxCourant ([0-9.eE+-]+) '
                         r'clamp ([0-9.eE+-]+) volErr ([0-9.eE+-]+) min ([0-9.eE+-]+) '
                         r'max ([0-9.eE+-]+)', out, re.M):
        geo.append(dict(t=float(m.group(1)), fluxes=int(m.group(2)),
                        limited=int(m.group(3)), rounds=int(m.group(4)),
                        min_factor=float(m.group(5)), courant=float(m.group(6)),
                        clamp=float(m.group(7)), vol_err=float(m.group(8)),
                        alpha_min=float(m.group(9)), alpha_max=float(m.group(10))))
    if not rows or start is None:
        raise SystemExit("the log does not hold a start line and the report lines "
                         "the checker needs")
    if (div is None) == (csf is None):
        raise SystemExit("the log has neither exactly one of the velocity line "
                         "(transport) and the csf line (momentum)")
    return dict(rows=rows, mass0=float(start.group(1)), kappa0=float(start.group(4)),
                iface0=int(start.group(5)), kind=div.group(1) if div else 'momentum',
                plic=plic, geo=geo,
                max_div=float(div.group(2)) if div else None,
                sigma_over_rho=float(csf.group(1)) if csf else None,
                nu=float(csf.group(2)) if csf else None,
                advect=(csf.group(3) == 'true') if csf else None,
                f_max=float(csf.group(5)) if csf else None)


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
        init_alpha = init[1]['alpha'] if init else None
        check_common(c, log, cent, alpha, nh, kap, band, failures, notes)
        if c['kind'] == 'oscillating':
            check_oscillating(c, log, cent, f, alpha, nh, kap, band, init,
                              failures, notes)
        elif c['kind'] == 'plic':
            check_plic(c, log, cent, f, init, failures, notes)
        elif c['kind'] == 'droplet':
            check_droplet(c, log, cent, f, alpha, nh, kap, band, init_alpha,
                          failures, notes)
        elif c['kind'] == 'translation':
            check_translation(c, log, cent, alpha, nh, kap, band, failures, notes)
        else:
            check_rotation(c, log, cent, alpha, nh, kap, band, failures, notes)
        if init and c['shape'] in ('plane', 'circle') and c['kind'] != 'plic':
            # the ellipse has no cell-by-cell closed form (its curvature varies
            # around the interface); its initial field is checked by
            # check_oscillating instead. A plic case is excluded because these
            # thresholds are the smoothed band's contract: on a one-cell front
            # the gradient normal is grid-aliased by construction, which is the
            # limitation the reconstruction removes and check_plic measures.
            check_initial_geometry(c, init, failures, notes)
        failures_exit(c, failures, notes)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def check_common(c, log, cent, alpha, nh, kap, band, failures, notes):
    # the prescribed field has to be discretely divergence-free, or nothing
    # below means anything: conservation is a property of the transport on such
    # a flux, not of the transport alone. The momentum mode has no prescribed
    # field -- its phi is the one the coupling solves -- and the continuity of
    # that is asserted where the coupling is checked.
    if log['max_div'] is not None:
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


def check_droplet(c, log, cent, f, alpha, nh, kap, band, init_alpha, failures, notes):
    """The momentum mode: a frozen interface and the pressure-velocity response
    to the surface force it implies."""
    U, p, csf = f['U'], f['p'], f['csf']
    dv = c['dx'][0] * c['dx'][1] * c['dx'][2]        # a block mesh, so uniform
    last = log['rows'][-1]
    r = c['radius']
    # the kernel's momentum is in acceleration form, so its p is p/rho and the
    # jump a 2D droplet of radius R has to carry is sigma/(rho R)
    target = c['sigma'] / (c['rho'] * r)

    # the jump, recomputed from the field rather than taken from the report
    ins = [p[i] for i in range(len(p)) if alpha[i] > 0.99]
    outs = [p[i] for i in range(len(p)) if alpha[i] < 0.01]
    dp = sum(ins) / len(ins) - sum(outs) / len(outs)
    notes.append("Laplace jump %.9g from the field (%.9g reported) for sigma/(rho R) "
                 "= %.9g: %.3f%% ; the far field sits at %.3g, its anchor"
                 % (dp, last['dp'], target, 100.0 * (dp / target - 1.0),
                    sum(outs) / len(outs)))
    if abs(dp / target - 1.0) > 0.03:
        failures.append("the pressure jump is %.6f, not sigma/(rho R) = %.6f"
                        % (dp, target))
    if dp <= 0:
        failures.append("the pressure jump has the wrong sign: the surface force "
                        "expands the droplet instead of squeezing it")
    if abs(sum(outs) / len(outs)) > 1e-3 * target:
        failures.append("the phase outside the band is %.3g off its own anchor, so "
                        "the reported jump is measuring the wrong thing"
                        % (sum(outs) / len(outs)))

    # one-way coupling: nothing moves the interface, so the phase volume and
    # alpha itself are exactly what the case started with
    if init_alpha is not None:
        frozen = max(abs(a - b) for a, b in zip(alpha, init_alpha))
        notes.append("alpha unchanged by the momentum mode: %.3e over %d cells"
                     % (frozen, len(alpha)))
        if frozen > 1e-11:
            failures.append("the momentum mode moved the phase fraction by %.3e: it "
                            "is meant to be frozen" % frozen)
    mass = sum(alpha) * dv
    notes.append("phase volume %.9g from the result for %.9g reported: alpha is "
                 "frozen, so the two have to agree to the writer's precision "
                 "(%.3e relative)" % (mass, log['mass0'], abs(mass / log['mass0'] - 1.0)))
    if abs(mass / log['mass0'] - 1.0) > 1e-10:
        failures.append("the frozen interface changed its volume by %.3e"
                        % (mass / log['mass0'] - 1.0))

    # the force is band-localised and internal: it can only bend the interface,
    # never translate the droplet, so its volume sum has to vanish
    net = [sum(dv * csf[i][k] for i in range(len(csf))) for k in range(3)]
    fsum = sum(dv * math.hypot(*csf[i]) for i in range(len(csf)))
    notes.append("surface force: |f| over the mesh %.6g, net %s = %.3e of it"
                 % (fsum, ['%.3e' % v for v in net],
                    math.hypot(*net) / fsum if fsum else float('nan')))
    if fsum == 0 or math.hypot(*net) / fsum > 1e-10:
        failures.append("the discrete surface force does not cancel as a pair: it "
                        "would accelerate the droplet as a whole")

    # the same statement one step later: the box is closed, so the momentum the
    # force puts in has nowhere to go but a circulation. The net of rho*U over
    # the mesh is the integral the surface force's own net would drive.
    mom = [sum(dv * U[i][k] for i in range(len(U))) for k in range(3)]
    ms = sum(dv * math.hypot(*U[i]) for i in range(len(U)))
    notes.append("momentum in the box: |integral U dV| %s = %.3e of integral |U| dV "
                 "%.6g" % (['%.3e' % v for v in mom],
                           math.hypot(*mom) / ms if ms else float('nan'), ms))
    if ms == 0 or math.hypot(*mom) / ms > 1e-3:
        failures.append("the closed box acquired a net translation (%.3e of its own "
                        "|U|): the surface force is not internal"
                        % (math.hypot(*mom) / ms if ms else float('nan')))

    # away from the band there is no force, so the momentum balance there is
    # grad(p) = 0 and each phase's own bulk is at one pressure: the ramp belongs
    # to the band and nowhere else. Measured well inside each bulk, away from the
    # profile's own tails.
    worst = 0.0
    for sel in ([v for v, a in zip(p, alpha) if a > 0.999],
                [v for v, a in zip(p, alpha) if a < 0.001]):
        worst = max(worst, (max(sel) - min(sel)) / abs(target))
    notes.append("pressure in the two bulks: flat to %.3f%% of the jump (no force "
                 "there, so grad(p) = 0)" % (100.0 * worst))
    if worst > 0.02:
        failures.append("the pressure ramps by %.3f%% of the jump inside a bulk "
                        "phase, where the force is zero" % (100.0 * worst))

    # and the response is mirror-symmetric about the droplet's own axes: the
    # force cancels in pairs, so what the coupling solves has nothing to lean on.
    # This is the identity of the unit test read out of the SOLVED fields.
    idx = {(round(q[0], 9), round(q[1], 9)): i for i, q in enumerate(cent)}
    umax_f = max(math.hypot(*v) for v in U)
    dp_asym = max(abs(p[i] - p[idx[(round(1.0 - q[0], 9), round(q[1], 9))]])
                  for i, q in enumerate(cent)
                  if (round(1.0 - q[0], 9), round(q[1], 9)) in idx)
    ux_asym = max(abs(U[i][0] + U[idx[(round(1.0 - q[0], 9), round(q[1], 9))]][0])
                  for i, q in enumerate(cent)
                  if (round(1.0 - q[0], 9), round(q[1], 9)) in idx)
    notes.append("mirror about x = 0.5: |dp| up to %.3e (%.4f%% of the jump), "
                 "|Ux + Ux'| up to %.3e (%.3f%% of max |U|)"
                 % (dp_asym, 100.0 * dp_asym / target, ux_asym,
                    100.0 * ux_asym / umax_f))
    if dp_asym > 1e-2 * target or ux_asym > 1e-2 * umax_f:
        failures.append("the response is not mirror-symmetric (dp %.3e, Ux %.3e): "
                        "the surface force left a net push" % (dp_asym, ux_asym))

    # the continuity of the SOLVED flux, which is what closes the balance
    cont = max(row['cont'] for row in log['rows'])
    notes.append("max |div(phi)| of the solved flux over the run: %.3e with up to "
                 "%d correctors (%d configured)"
                 % (cont, max(row['corr'] for row in log['rows']), c['correctors']))
    if cont > 1e-8:
        failures.append("the coupling leaves div(phi) at %.3e, so the pressure does "
                        "not close the force" % cont)

    # the spurious current: what the mismatched cancellation of f and grad(p)
    # leaves behind. Measured against the capillary-inertial velocity
    # sqrt(sigma/(rho R)), the speed the same surface tension would reach in an
    # inviscid drop of this size.
    u_ci = math.sqrt(target)
    umax = max(row['umax'] for row in log['rows'])
    notes.append("spurious current: max |U| %.6g = %.2f%% of sqrt(sigma/(rho R)) = "
                 "%.4f over %d steps (nu dt/dx^2 = %.3f)"
                 % (umax, 100.0 * umax / u_ci, u_ci, c['steps'],
                    c['nu'] * c['dt'] / c['dx'][0] ** 2))
    if not all(math.isfinite(row['umax']) for row in log['rows']):
        failures.append("the velocity diverged")
    if umax > 0.1 * u_ci:
        failures.append("the surface force drives |U| to %.2f%% of the capillary "
                        "velocity: the pressure is not balancing it"
                        % (100.0 * umax / u_ci))
    # and it has to still be building, not running away: the increments shrink
    g = [log['rows'][i]['umax'] - log['rows'][i - 1]['umax']
         for i in range(1, len(log['rows']))]
    if len(g) > 2 and max(g[1:]) > g[0]:
        failures.append("the spurious current grows faster as it goes (%s): that is "
                        "an instability, not a gyre building" % ['%.3g' % v for v in g])
    notes.append("|U| history %s (increments %s)"
                 % (['%.4g' % row['umax'] for row in log['rows']],
                    ['%.3g' % v for v in g]))

    # the interface geometry the force was built from, on the result the case
    # wrote rather than on the report line
    kmean = sum(kap[i] for i in band) / len(band)
    notes.append("%d interface cells carry the force: band mean kappa %.4f for "
                 "-1/R = %.4f (the same field the run started with)"
                 % (len(band), kmean, -1.0 / r))
    if abs(kmean / (-1.0 / r) - 1.0) > 0.08:
        failures.append("the band reads kappa %.4f, not -1/R = %.4f" % (kmean, -1.0 / r))


def check_oscillating(c, log, cent, f, alpha, nh, kap, band, init, failures, notes):
    """The two-way momentum mode: the solved flux carries alpha and the surface
    force is rebuilt from the shape it leaves, so an ellipse released off the
    circle oscillates in its n = 2 mode and decays. The period is the number
    the kernel is not given, and the case exists to decide WHICH law: one-fluid
    momentum means the ambient carries half the inertia, w^2 = 3 sigma/(rho
    R^3) on R = sqrt(a b), not the free-drop 6 sigma/(rho R^3) -- 35% apart."""
    if not log['advect']:
        failures.append("the oscillating case ran with advect false: the "
                        "interface would be frozen and nothing below means anything")
        return

    dv = c['dx'][0] * c['dx'][1] * c['dx'][2]
    r_eq = math.sqrt(c['a'] * c['b'])
    t_lin = 2.0 * math.pi * math.sqrt(c['rho'] * r_eq ** 3 / (3.0 * c['sigma']))
    target = c['sigma'] / (c['rho'] * r_eq)
    u_ci = math.sqrt(target)
    rows = log['rows']

    if abs(rows[-1]['t'] - c['end_time']) > c['dt']:
        failures.append("the run stopped at t=%g, not the endTime %g the case "
                        "asked for" % (rows[-1]['t'], c['end_time']))

    # mass: the box is closed and every wall carries alpha = 0, so the boundary
    # exchange is zero and the two-way loop conserves the phase volume to
    # round-off -- the flux telescopes over the interior faces, whatever it is.
    merr = max(abs(r['massErr']) for r in rows)
    notes.append("phase volume drift over the run: %.3e relative (closed box)"
                 % merr)
    if merr > 1e-11:
        failures.append("the two-way loop lost %.3e of the phase volume: the "
                        "advected transport is not conservative on the solved "
                        "flux" % merr)

    # nothing in an internal surface force can translate the drop, so the
    # centroid stays pinned; what moves it is the linear solver's sweep
    # asymmetry, which accumulates in x
    drift = max(max(abs(r['cx'] - c['center'][0]), abs(r['cy'] - c['center'][1]))
                for r in rows)
    notes.append("centroid drift over the run: %.3e (BiCGStab sweep asymmetry)"
                 % drift)
    if drift > 5e-4:
        failures.append("the droplet's centroid moved %.3e: the surface force is "
                        "internal and cannot translate it" % drift)

    # the crossing clock: the reported aspect - 1 linearly interpolated to zero.
    # cross[0] is the first compression, cross[1] the return, cross[2] the
    # second; consecutive same-direction crossings are one period apart.
    s = [(r['t'], r['aspect'] - 1.0) for r in rows]
    cross = [t0 + (t1 - t0) * s0 / (s0 - s1)
             for (t0, s0), (t1, s1) in zip(s, s[1:]) if (s0 > 0) != (s1 > 0)]
    if len(cross) < 3:
        failures.append("only %d crossings of aspect = 1: the shape is not "
                        "oscillating about the circle (a wrong-sign surface "
                        "force grows away from it instead)" % len(cross))
        return
    periods = [y - x for x, y in zip(cross, cross[2:])]
    t_meas = sum(periods) / len(periods)
    t_free = 2.0 * math.pi * math.sqrt(c['rho'] * r_eq ** 3 / (6.0 * c['sigma']))
    notes.append("aspect = 1 at t = %s; period %.4f s against the density-matched "
                 "2D Lamb period 2 pi sqrt(rho R^3/(3 sigma)) = %.4f s on "
                 "R = sqrt(a b) = %.5f (%.1f%%; the free-drop law would say "
                 "%.4f s)"
                 % (['%.3f' % t for t in cross], t_meas, t_lin, r_eq,
                    100.0 * (t_meas / t_lin - 1.0), t_free))
    if abs(t_meas / t_lin - 1.0) > 0.12:
        failures.append("the measured period %.4f s is %.1f%% off the "
                        "density-matched analytic %.4f s" % (t_meas,
                        100.0 * (t_meas / t_lin - 1.0), t_lin))
    # the first quarter period carries the full release amplitude, where the
    # finite-amplitude softening bites; it runs long but not arbitrarily so
    if not 0.22 < cross[0] / t_lin < 0.38:
        failures.append("the first crossing at t=%.3f is %.0f%% of T_lin/4: the "
                        "release phase does not read as the n = 2 mode"
                        % (cross[0], 400.0 * cross[0] / t_lin))

    # the swing is real and decaying on both axes
    amin = min(r['aspect'] for r in rows)
    amp0 = max(r['aspect'] for r in rows if r['t'] < cross[0]) - 1.0
    amp_in1 = 1.0 - min(r['aspect'] for r in rows if cross[0] < r['t'] < cross[1])
    amp_out2 = max(r['aspect'] for r in rows if cross[1] < r['t'] < cross[2]) - 1.0
    amp_in2 = 1.0 - min(r['aspect'] for r in rows if r['t'] > cross[2])
    notes.append("aspect range [%.4f, %.4f]; amplitudes: release %.3f, first "
                 "compression %.3f, second extension %.3f, second compression "
                 "%.3f" % (amin, max(r['aspect'] for r in rows), amp0, amp_in1,
                           amp_out2, amp_in2))
    if amin > 0.90:
        failures.append("the shape never swings past aspect %.3f: that is not an "
                        "n = 2 oscillation" % amin)
    if amp0 < 0.2:
        failures.append("the case is released from a/b = %.3f but the largest "
                        "reported amplitude is %.3f" % (c['a'] / c['b'], amp0))
    if amp_out2 > 0.75 * amp0 or amp_out2 < 0.02:
        failures.append("the second extension's amplitude %.3f against the "
                        "release's %.3f: the oscillation is either not decaying "
                        "or already dead" % (amp_out2, amp0))
    if amp_in2 >= amp_in1 or amp_in2 < 0.02:
        failures.append("the second compression %.3f against the first's %.3f: "
                        "the decay is not monotone across a period (or it is "
                        "dead)" % (amp_in2, amp_in1))

    # the jump breathes with the shape and the flow; its run-long mean sits on
    # the Laplace value of the volume-equivalent radius
    dps = [r['dp'] for r in rows]
    dmean = sum(dps) / len(dps)
    notes.append("dp = pIn - pOut over the run: [%.4f, %.4f], mean %.4f, against "
                 "sigma/(rho R) = %.4f"
                 % (min(dps), max(dps), dmean, target))
    if not 0.92 < dmean / target < 1.02:
        failures.append("the mean pressure jump %.4f is %.1f%% off the Laplace "
                        "value %.4f on R = sqrt(a b)"
                        % (dmean, 100.0 * (dmean / target - 1.0), target))

    # the release reaches the capillary speed (plus the band-edge force spike
    # the case header documents); past the ring-down the flow decays with the
    # shape amplitude
    umax_all = max(r['umax'] for r in rows)
    late = max(r['umax'] for r in rows if r['t'] > 0.25)
    notes.append("max |U|: %.3f in the release against the capillary-inertial "
                 "speed %.3f, %.3f past t = 0.25" % (umax_all, u_ci, late))
    if not all(math.isfinite(r['umax']) for r in rows):
        failures.append("the velocity diverged")
    if late > 0.5 * u_ci:
        failures.append("past the startup the flow still reaches %.3f = %.0f%% "
                        "of the capillary speed: it is not decaying with the "
                        "shape" % (late, 100.0 * late / u_ci))

    # the continuity of the solved flux, which carries alpha
    cont = max(r['cont'] for r in rows)
    notes.append("max |div(phi)| of the solved flux over the run: %.3e with up "
                 "to %d correctors (%d configured)"
                 % (cont, max(r['corr'] for r in rows), c['correctors']))
    if cont > 1e-8:
        failures.append("the coupling leaves div(phi) at %.3e" % cont)

    # the result on disk has to carry the same statement, read independently of
    # the report line
    mass = sum(alpha) * dv
    if abs(mass / log['mass0'] - 1.0) > 1e-10:
        failures.append("the written phase volume is %.3e off the run's own"
                        % (mass / log['mass0'] - 1.0))
    m = sum(alpha)
    cx = sum(v * p[0] for v, p in zip(alpha, cent)) / m
    cy = sum(v * p[1] for v, p in zip(alpha, cent)) / m
    ixx = sum(v * (p[0] - cx) ** 2 for v, p in zip(alpha, cent)) / m
    iyy = sum(v * (p[1] - cy) ** 2 for v, p in zip(alpha, cent)) / m
    aspect = math.sqrt(ixx / iyy)
    notes.append("from the written field: volume %.9g, aspect %.6f against the "
                 "reported %.6f" % (mass, aspect, rows[-1]['aspect']))
    if abs(aspect - rows[-1]['aspect']) > 1e-9:
        failures.append("the written field's aspect %.9f is not the reported "
                        "%.9f" % (aspect, rows[-1]['aspect']))

    # the normals are unit and point into the phase, on the deformed shape as
    # on the round one; the band's mean curvature is the ellipse-scale one
    bad = [i for i in band if abs(math.hypot(nh[i][0], nh[i][1], nh[i][2]) - 1.0)
           > 1e-6]
    outward = max(nh[i][0] * (cent[i][0] - cx) + nh[i][1] * (cent[i][1] - cy)
                  for i in band)
    kmean = sum(kap[i] for i in band) / len(band)
    notes.append("%d band cells at the end: %d non-unit normals, max outward "
                 "component %.3f, band mean kappa %.3f"
                 % (len(band), len(bad), outward, kmean))
    if bad:
        failures.append("%d band cells carry a non-unit normal" % len(bad))
    if outward > 0.05:
        failures.append("a band cell's normal points out of the phase (%.3f): "
                        "nhat has to read towards increasing alpha" % outward)
    if not -9.0 < kmean < -2.5:
        failures.append("the band mean kappa %.3f is not a closed deformed "
                        "interface's curvature in this geometry" % kmean)

    # the initial field the run started from: the released ellipse, its volume
    # the one the mass statement above is anchored to
    if init:
        icent, ifl = init
        ia, inh, ikap = ifl['alpha'], ifl['nHat'], ifl['kappa']
        iband = [i for i, v in enumerate(inh)
                 if math.hypot(v[0], v[1], v[2]) > 0.5]
        m0 = sum(ia)
        icx = sum(v * p[0] for v, p in zip(ia, icent)) / m0
        icy = sum(v * p[1] for v, p in zip(ia, icent)) / m0
        ixx = sum(v * (p[0] - icx) ** 2 for v, p in zip(ia, icent)) / m0
        iyy = sum(v * (p[1] - icy) ** 2 for v, p in zip(ia, icent)) / m0
        a0 = math.sqrt(ixx / iyy)
        outward0 = max(inh[i][0] * (icent[i][0] - icx)
                       + inh[i][1] * (icent[i][1] - icy) for i in iband)
        ikmean = sum(ikap[i] for i in iband) / len(iband)
        h = ((c['a'] - c['b']) / (c['a'] + c['b'])) ** 2
        k_true = -2.0 * math.pi / (math.pi * (c['a'] + c['b'])
                                   * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h))))
        notes.append("at t=0: volume %.9g for %.9g reported, aspect %.4f for the "
                     "sharp a/b = %.4f, %d band cells, max outward normal "
                     "component %.3f, band mean kappa %.3f for the arc-length "
                     "mean %.3f (the threshold band under-resolves the tips)"
                     % (m0 * dv, log['mass0'], a0, c['a'] / c['b'], len(iband),
                        outward0, ikmean, k_true))
        if abs(m0 * dv / log['mass0'] - 1.0) > 1e-10:
            failures.append("the written initial volume is %.3e off the run's "
                            "own" % (m0 * dv / log['mass0'] - 1.0))
        if not 1.2 < a0 < 1.5:
            failures.append("the initial aspect %.3f is not the released "
                            "ellipse's" % a0)
        if outward0 > 0.05:
            failures.append("an initial band normal points out of the phase "
                            "(%.3f)" % outward0)
        if not 0.6 < ikmean / k_true < 1.05:
            failures.append("the initial band mean kappa %.3f is %.0f%% of the "
                            "ellipse's arc-length mean %.3f"
                            % (ikmean, 100.0 * ikmean / k_true, k_true))


def check_plic(c, log, cent, f, init, failures, notes):
    """The reconstruction against the closed forms a cell-sharp disk gives: the
    band it covers, the fraction each plane reads back, the area a one-cell band
    spans and the volume the exact filling holds. The planes are re-read from what
    the solver wrote, so the geometry is checked on the mesh's own coordinates."""
    if not init:
        raise SystemExit("a plic case has to write vtkInitial: the field the "
                         "reconstruction is measured against is the analytic one")
    if not log['plic']:
        failures.append("the case asks for a reconstruction but the log holds no "
                        "plic line for the checker to follow")
        return
    n_cells = c['cells'][0] * c['cells'][1] * c['cells'][2]
    if len(cent) != n_cells:
        raise SystemExit("the mesh the case wrote is not the block it describes")
    v_cell = (c['vol'][0] * c['vol'][1] * c['vol'][2]) / n_cells
    R, dz = c['radius'], c['vol'][2]
    exact_area, perimeter = math.pi * R * R * dz, 2.0 * math.pi * R * dz
    tol = c['band_tolerance']
    # the file holds 12 significant digits, so anything re-derived from it is
    # compared with a little more slack than the solver's own report
    wt = 5e-10

    def plic_band(fields):
        return [i for i, v in enumerate(fields['nHatPlic'])
                if math.hypot(v[0], v[1], v[2]) > 0.5]

    def alpha_band(values):
        return [i for i, a in enumerate(values) if tol < a < 1.0 - tol]

    # ---- the exact filling: alpha is the disk's own overlap area ----
    icent, ifields = init
    ia = ifields['alpha']
    mass = sum(ia) * v_cell
    band0 = alpha_band(ia)
    notes.append("exact filling over %d cells: phase volume %.12g for pi R^2 dz = "
                 "%.12g (%.3g relative), %d band cells with %g < alpha < %g"
                 % (len(ia), mass, exact_area, abs(mass / exact_area - 1.0),
                    len(band0), tol, 1.0 - tol))
    if abs(mass / exact_area - 1.0) > 1e-8:
        failures.append("the exact filling holds %.6g of phase, not the disk's %.6g"
                        % (mass, exact_area))
    if abs(log['mass0'] / exact_area - 1.0) > 1e-8:
        failures.append("the solver's own phase volume %.6g is not the disk's %.6g"
                        % (log['mass0'], exact_area))
    # every cell the band does not name is full or empty to the last digits, which
    # is what makes this interface one cell thick and the gradient path blind on it
    band_set = set(band0)
    sharp = max(min(ia[i], 1.0 - ia[i]) for i in range(len(ia)) if i not in band_set)
    notes.append("outside the band the closest any cell gets to half-filled is %.3g"
                 % sharp)
    if sharp > 1e-10:
        failures.append("cells outside the band are not full or empty (closest %.3g), "
                        "so the interface is not one cell thick" % sharp)

    # ---- the band the planes cover is the band alpha names ----
    gplic0, gband0 = log['plic'][0], plic_band(ifields)
    notes.append("at t=0 the reconstruction reports band %d, the planes it wrote are "
                 "on %d cells and alpha names %d"
                 % (gplic0['band'], len(gband0), len(band0)))
    if set(gband0) != set(band0):
        failures.append("the planes are not on exactly the cells whose alpha is in "
                        "(%g, %g): %d vs %d" % (tol, 1.0 - tol, len(gband0), len(band0)))
    if gplic0['band'] != len(band0):
        failures.append("the reported band %d is not the %d cells alpha puts on the "
                        "interface" % (gplic0['band'], len(band0)))
    # the gradient path the smeared-band operators use is over-inclusive here: it
    # marks cells that hold no interface at all, which is why the reconstruction
    # takes its band from alpha instead
    ggrad0 = [i for i, v in enumerate(ifields['nHat']) if math.hypot(v[0], v[1], v[2]) > 0.5]
    missed = set(band0) - set(ggrad0)
    extra = set(ggrad0) - set(band0)
    notes.append("the gradient threshold marks %d cells at t=0: %d of them hold no "
                 "interface, and it misses %d of the cells the reconstruction covers"
                 % (len(ggrad0), len(extra), len(missed)))
    if missed:
        failures.append("the gradient path misses %d band cells the reconstruction "
                        "covers" % len(missed))
    if not extra:
        failures.append("the gradient threshold marks no cell outside the band, so "
                        "this case no longer measures what PLIC adds")

    # ---- the planes themselves ----
    nrm, off, area, ctr = (ifields['nHatPlic'], ifields['offsetPlic'],
                           ifields['areaPlic'], ifields['centrePlic'])
    unit = max(abs(math.hypot(*nrm[i]) - 1.0) for i in gband0)
    planar = max(abs(nrm[i][0] * ctr[i][0] + nrm[i][1] * ctr[i][1]
                     + nrm[i][2] * ctr[i][2] - off[i]) for i in gband0)
    if min(area[i] for i in gband0) <= 0.0:
        failures.append("a band cell was cut with no interface area")
    if max(abs(area[i]) for i in range(len(area)) if i not in set(gband0)) > 0.0:
        failures.append("a cell outside the band carries interface area")
    a_sum = sum(area[i] for i in gband0)
    notes.append("%d planes: |n| off 1 by %.3g, the cut face lies on its plane to "
                 "%.3g, area summed %.9g vs the reported %.9g"
                 % (len(gband0), unit, planar, a_sum, gplic0['area']))
    if unit > wt or planar > wt:
        failures.append("the planes are not the unit planes they report (|n|-1 up to "
                        "%.3g, offset off by %.3g)" % (unit, planar))
    if abs(a_sum / gplic0['area'] - 1.0) > 1e-10:
        failures.append("the written areas sum to %.9g but the log reported %.9g"
                        % (a_sum, gplic0['area']))
    if not (0.85 < a_sum / perimeter < 1.05):
        failures.append("a one-cell band spans %.9g of interface, %.1f%% of the "
                        "extruded perimeter %.9g" % (a_sum, 100.0 * a_sum / perimeter,
                                                     perimeter))

    # ---- the direction and curvature the planes imply ----
    kap_p = ifields['kappaPlic']
    ang = []
    for i in gband0:
        p = ctr[i]
        rx, ry = p[0] - c['center'][0], p[1] - c['center'][1]
        r = math.hypot(rx, ry)
        if r < 1e-9:
            continue
        d = -(nrm[i][0] * rx + nrm[i][1] * ry) / r
        ang.append(math.degrees(math.acos(max(-1.0, min(1.0, d)))))
    mean_a, max_a = sum(ang) / len(ang), max(ang)
    kmean = sum(kap_p[i] for i in gband0) / len(gband0)
    notes.append("normals at the cut faces, against the analytic radial: %.2f deg mean, "
                 "%.2f deg max -- sharper than the same normals read at cell centres, "
                 "because a cut face lies on the interface and a cell centre does not; "
                 "band mean kappaPlic %.4f for -1/R = %.4f (%.3f of it), log reported %.4f"
                 % (mean_a, max_a, kmean, -1.0 / R, kmean * R, gplic0['kappa']))
    if abs(kmean * R + 1.0) > 0.45 or kmean >= 0.0:
        failures.append("the reconstructed band's curvature is %.4f, not the -1/R = %.4f "
                        "a one-cell stencil can resolve" % (kmean, -1.0 / R))
    if abs(kmean - gplic0['kappa']) > 1e-9:
        failures.append("the written kappaPlic means %.6f but the log reported %.6f"
                        % (kmean, gplic0['kappa']))
    if mean_a > 12.0 or max_a > 20.0:
        failures.append("the reconstructed normal is off the radial direction by %.2f "
                        "deg mean / %.2f max" % (mean_a, max_a))

    # ---- what the fit says its neighbourhood was worth ----
    # The scatter is the distance the planes of one neighbourhood still disagree by
    # after the tilt, so it is the part of the interface a single plane cannot hold:
    # over the disk's own arc it should stay well under a cell, and the rounds should
    # stop inside the budget the library gives them.
    if c['normal_choice'] == 'fit':
        dx = min(c['vol'][0] / c['cells'][0], c['vol'][1] / c['cells'][1])
        sc = [l['scatter'] / dx for l in log['plic'] if l['scatter'] is not None]
        rs = [l['rounds'] for l in log['plic'] if l['rounds'] is not None]
        if len(sc) != len(log['plic']) or len(rs) != len(log['plic']):
            failures.append("a fit case reported its neighbourhood on %d of its %d "
                            "reconstructions" % (len(sc), len(log['plic'])))
        else:
            notes.append("the planes of one neighbourhood disagree by %.4g to %.4g of a cell "
                         "over the %d reports, and the fit stopped within %d rounds of the "
                         "four it is given"
                         % (min(sc), max(sc), len(sc), max(rs)))
            if max(sc) > 0.06:
                failures.append("the fitted neighbourhoods disagree by %.4g of a cell, more "
                                "than the sagitta of this disk over the cell" % max(sc))
            if not rs or max(rs) > 4 or min(rs) < 1:
                failures.append("the fit took %d rounds where the library budgets four"
                                % (max(rs) if rs else 0))

    # ---- through the run: the planes follow the field, not the analytic disk ----
    for line in log['plic']:
        if line['read_back'] > 1e-12:
            failures.append("at %s a plane reads its cell's fraction back with error "
                            "%.3g" % (line['stamp'], line['read_back']))
    for r in log['rows']:
        if abs(r['massErr']) > 1e-12:
            failures.append("the closed rotation lost %.3g of the phase volume by t=%g"
                            % (abs(r['massErr']), r['t']))
    fb, fband = plic_band(f), alpha_band(f['alpha'])
    area_sum = sum(f['areaPlic'][i] for i in fb)
    if c['transport'] == 'geometric':
        notes.append("after %d geometric steps the band is %d cells (was %d) and spans "
                     "%.9g of interface for the perimeter %.9g: the front stayed one "
                     "cell thick, so the area is the circle's own and not the width of "
                     "a smeared band"
                     % (c['steps'], len(fb), len(band0), area_sum, perimeter))
    else:
        notes.append("after %d steps the band is %d cells (was %d) and spans %.9g of "
                     "interface: a band thicker than one cell counts the cut face of "
                     "every cell in the smearing, so this is the front's thickness, "
                     "not the circle's perimeter"
                     % (c['steps'], len(fb), len(band0), area_sum))
    if set(fb) != set(fband):
        failures.append("the last reconstruction written is not the one for the last "
                        "field written: %d planes vs %d band cells" % (len(fb), len(fband)))
    if log['plic'][-1]['band'] != len(fb):
        failures.append("the last report says band %d, the file holds %d planes"
                        % (log['plic'][-1]['band'], len(fb)))
    if abs(area_sum / log['plic'][-1]['area'] - 1.0) > 1e-10:
        failures.append("the last report's area and the last file's do not agree")
    if len(fb) < len(band0):
        failures.append("advecting a sharp front narrowed its band from %d cells to %d, "
                        "which no conservative transport does" % (len(band0), len(fb)))
    notes.append("read-back stayed under %.3g over %d reports, so every plane keeps "
                 "the fraction of the cell it was cut in"
                 % (max(l['read_back'] for l in log['plic']), len(log['plic'])))

    if c['transport'] == 'geometric':
        check_geometric(c, log, cent, f, init, v_cell, perimeter, len(band0),
                        failures, notes)


def disk_cell_alpha(cent, dx, dy, centre, radius):
    """The exact alpha of a disk of `radius` centred at `centre`, cell by cell, as
    the area it shares with each axis-aligned cell divided by that cell's area. A
    disk turned through any angle is the same disk with its centre moved along its
    radius, so this is what a rigid rotation leaves - and the checker writes it
    from the chord's antiderivative rather than from what the solver filled with."""
    out = []
    for p in cent:
        x1, x2 = p[0] - centre[0] - 0.5 * dx, p[0] - centre[0] + 0.5 * dx
        y1, y2 = p[1] - centre[1] - 0.5 * dy, p[1] - centre[1] + 0.5 * dy
        a = _disk_rectangle(x1, x2, y1, y2, radius) / (dx * dy)
        out.append(min(1.0, max(0.0, a)))
    return out


def _circle_quadrant(x, y, radius):
    """The signed area of the rectangle (0,0)-(x,y) inside the disk centred at the
    origin: signed in each coordinate so that the four differences of it over a
    rectangle's corners give that rectangle's own overlap with the disk."""
    s = (1.0 if x >= 0.0 else -1.0) * (1.0 if y >= 0.0 else -1.0)
    x, y = abs(x), abs(y)
    r2 = radius * radius
    if x * x + y * y <= r2:
        return s * x * y

    def anti(t):  # integral of sqrt(R^2 - x^2)
        return 0.5 * (t * math.sqrt(max(0.0, r2 - t * t)) + r2 * math.asin(t / radius))

    a = min(x, radius)
    b = min(a, math.sqrt(max(0.0, r2 - y * y)))
    return s * (b * y + anti(a) - anti(b))


def _disk_rectangle(x1, x2, y1, y2, radius):
    return (_circle_quadrant(x2, y2, radius) - _circle_quadrant(x1, y2, radius)
            - _circle_quadrant(x2, y1, radius) + _circle_quadrant(x1, y1, radius))


def rotated_centre(c, angle):
    """Where the disk's centre sits after the case's own rotation carried it."""
    ax = c['axis']
    phi = math.atan2(c['center'][1] - ax[1], c['center'][0] - ax[0])
    return (ax[0] + c['r0'] * math.cos(phi + angle), ax[1] + c['r0'] * math.sin(phi + angle))


def check_geometric(c, log, cent, f, init, v_cell, perimeter, n_band0,
                    failures, notes):
    """What the planes, and not a flux of alpha, owe: the phase volume is the sum
    of what each face carried across it, so it is conserved by construction; alpha
    stays in [0,1] because the limiter is what keeps it there and nothing clips it
    afterwards; and the field after a turn is the analytic disk moved, cell by
    cell, which is the only shape a rigid rotation can leave."""
    if not log['geo']:
        failures.append("the case asks for geometric transport but the log holds no "
                        "geo line, so the checker cannot follow the fluxes")
        return

    # ---- conservation, boundedness and the limiter over every report ----
    worst_err = max(abs(g['vol_err']) for g in log['geo']) / log['mass0']
    worst_clamp = max(g['clamp'] for g in log['geo'])
    notes.append("%d reported geometric steps: the fluxes the faces of one cell "
                 "disagree about are at most %.3g of the phase volume and what the "
                 "bounds had to clip at most %.3g of it, so the phase moved face to "
                 "face and was not reconciled afterwards"
                 % (len(log['geo']), worst_err, worst_clamp))
    if worst_err > 1e-13:
        failures.append("a step's face fluxes disagree by %.3g of the phase volume"
                        % worst_err)
    if worst_clamp > 1e-15:
        failures.append("a step had to clip %.3g of the phase back into its cells, "
                        "so the limiter did not keep alpha in [0,1] on its own"
                        % worst_clamp)
    for g in log['geo']:
        if g['alpha_min'] < -1e-15 or g['alpha_max'] > 1.0 + 1e-15:
            failures.append("at t=%g the transported field left [0,1]: [%.6g, %.6g]"
                            % (g['t'], g['alpha_min'], g['alpha_max']))
    lo, hi = min(f['alpha']), max(f['alpha'])
    if lo < -1e-15 or hi > 1.0 + 1e-15:
        failures.append("the field written after %d steps is outside [0,1]: "
                        "[%g, %g]" % (c['steps'], lo, hi))
    rounds = max(g['rounds'] for g in log['geo'])
    limited = max(g['limited'] for g in log['geo'])
    smallest = min(g['min_factor'] for g in log['geo'])
    courant = max(g['courant'] for g in log['geo'])
    notes.append("the two-limit pass needed up to %d rounds to settle and restricted "
                 "up to %d of the %d faces carrying flux; the smallest factor any "
                 "cell gave was %.3g at a Courant number of %.3g, which is where two "
                 "outflow prisms of one cell overlap at a corner"
                 % (rounds, limited, max(g['fluxes'] for g in log['geo']),
                    smallest, courant))
    if rounds > 100:
        failures.append("the two-limit pass had not settled after %d rounds" % rounds)
    if smallest < 0.25:
        failures.append("the limiter took a face flux down to %.3g of what the swept "
                        "prism asked for, which is the geometry failing rather than "
                        "the corner overlap it is there for" % smallest)

    # ---- the closed form the checker writes for itself ----
    exact0 = disk_cell_alpha(cent, c['dx'][0], c['dx'][1],
                              (c['center'][0], c['center'][1]), c['radius'])
    d0 = max(abs(e - a) for e, a in zip(exact0, init[1]['alpha']))
    notes.append("the checker's own disk/rectangle integral reproduces the field the "
                 "solver filled to %.3g per cell, so the reference below is the "
                 "geometry and not a copy of the solver's arithmetic" % d0)
    if d0 > 1e-10:
        raise SystemExit("the checker's closed form and the solver's exact filling "
                         "disagree by %.3g, so neither can be the reference" % d0)

    turn = c['omega'] * log['rows'][-1]['t']
    exact = disk_cell_alpha(cent, c['dx'][0], c['dx'][1], rotated_centre(c, turn),
                            c['radius'])
    err = [abs(a - e) for a, e in zip(f['alpha'], exact)]
    sym_diff = sum(err) * v_cell
    disk = math.pi * c['radius'] ** 2 * c['vol'][2]
    max_err = max(err)
    fb = [i for i, a in enumerate(f['alpha'])
          if c['band_tolerance'] < a < 1.0 - c['band_tolerance']]
    mass = sum(f['alpha'])
    cx = sum(a * p[0] for a, p in zip(f['alpha'], cent)) / mass
    cy = sum(a * p[1] for a, p in zip(f['alpha'], cent)) / mass
    rc = rotated_centre(c, turn)
    off = math.hypot(cx - rc[0], cy - rc[1])
    dx = c['dx'][0]
    # a cell off by a whole half is one the front is on the wrong side of, which on
    # a one-cell interface is the only per-cell error with a meaning: a cell can sit
    # at 0.9 against an analytic 0.3 because the front crossed it a third of a cell
    # early, and that is the shape moving, not the field being wrong
    half = sum(1 for e in err if e > 0.5)
    notes.append("after t=%g the field is %.4f%% of the disk's area away from the disk "
                 "turned through %.4f rad - the boundary displaced by %.3f dx, %.2f%% "
                 "of R - with %d of the band's %d cells on the wrong side of their "
                 "front (worst %.3f of alpha) and its centroid %.3f dx off the "
                 "analytic centre"
                 % (log['rows'][-1]['t'], 100.0 * sym_diff / disk, turn,
                    sym_diff / perimeter / dx, 100.0 * (sym_diff / perimeter) / c['radius'],
                    half, len(fb), max_err, off / dx))
    if sym_diff / disk > 0.05:
        failures.append("the rotated field differs from the analytic disk by %.2f%% of "
                        "the disk's area: the planes deformed the interface more than "
                        "a rigid turn can" % (100.0 * sym_diff / disk))
    if half > 0.1 * len(fb):
        failures.append("%d of the %d interface cells put the front on the wrong side "
                        "of itself - %.0f%% of the band" % (half, len(fb),
                                                            100.0 * half / len(fb)))
    if off > 0.05 * c['radius']:
        failures.append("the droplet's centroid is %.4g off the centre its rotation "
                        "puts it at (%.2f%% of R)" % (off, 100.0 * off / c['radius']))

    # ---- what the turn did to the band and the area it spans ----
    area_sum = sum(f['areaPlic'][i] for i in fb)
    notes.append("the band after the turn is %d cells against the %d it started with, "
                 "and spans %.9g of interface for the %.9g the disk's edge extrudes: "
                 "a front carried by its planes widens only where the disk's own "
                 "curvature needs another cell to hold its fraction"
                 % (len(fb), n_band0, area_sum, perimeter))
    if len(fb) > 1.5 * n_band0:
        failures.append("the band grew from %d cells to %d, which is the transport "
                        "thickening the front instead of moving it" % (n_band0, len(fb)))
    if abs(area_sum / perimeter - 1.0) > 0.15:
        failures.append("after the turn the planes span %.9g of interface, %.1f%% off "
                        "the perimeter %.9g" % (area_sum, 100.0 * (area_sum / perimeter - 1.0),
                                                perimeter))




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
