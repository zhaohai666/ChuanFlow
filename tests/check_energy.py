#!/usr/bin/env python3
"""Run a 3D case with an energy block and verify the thermal field against the
closed forms each case is built to have.

Every mode compares against a number the solver is NOT given, and every one of
them is exact enough to be asserted tightly rather than absorbed into a slack
tolerance.

  fixed    - sealed box at rest, one isothermal face at each end. With no flow
      the equation is -div(alpha grad T) = 0 and its answer is the linear profile
      between the faces, so the field is compared point by point and the two wall
      rates have to be k A dT/dx.

  flux     - the same box with a prescribed heat flux on one face and a fixed
      temperature on the other. Fourier's law fixes the slope without fixing the
      level, so this is the mode that pins the SIGN convention of `heatFlux`:
      positive heats the fluid, the face that receives it reports +q A and the
      isothermal one -q A.

  source   - uniform volumetric generation between two isothermal faces. Two
      identities here are properties of the discrete problem rather than of the
      mesh: the field is mirror symmetric, and each face carries away exactly
      half the generated heat, S*Vol/2, to the linear-solver tolerance. The
      parabola is checked loosely on purpose: a cell-centred scheme imposes the
      wall value half a cell away from the first node, and the reported deviation
      is that curvature rather than a mystery.

  rest     - a uniformly heated fluid in a box closed on every face, under
      gravity. The Boussinesq force is then uniform, which is exactly what the
      pressure can absorb, so nothing has to move and the pressure has to be the
      hydrostatic of the force, dp/dz = beta (T - tRef) |g|. The residual this
      mode measures (~3e-4 m/s) is the closed box's own pressure-velocity
      coupling, not the energy equation's: a case with no thermal field driven by
      the equivalent body force reaches the same state. The assertions are
      therefore scaled to what a coupling that could NOT answer the force would
      produce - the Stokes velocity f H^2/nu, four orders above the residual that
      is actually left - and to the sign of the pressure gradient, which is what a
      misplaced buoyancy term gets wrong.

  decay    - a cosine mode in a closed box at rest, advanced with the transient
      solver. cos(m pi (i + 1/2)/N) is an exact eigenvector of the cell-centred
      Neumann Laplacian, so implicit Euler on a zero-flux field has the closed
      amplitude 1 / prod_n (1 + dt alpha lambda_d). The checker replays the
      solver's own step sequence and compares every reported point to 1e-8, then
      compares that factor with the continuous exp(-alpha k^2 t), which measures
      the diffusion operator's discretisation error instead of hiding it. The
      discrete mean is an invariant of an adiabatic box, so it is asserted.

  channel  - force-driven plane Poiseuille with the same heat flux entering
      through both walls, which couples convection to conduction and to the flux
      that closes continuity:
        * the mixing cup gains exactly the heat the surfaces feed,
          rho cp Phi (Tb_out - Tb_in) = sum of the wall rates, asserted to 1e-6
          relative with no analytic input at all;
        * the wall-to-bulk offset of the fully developed field is the Nusselt
          number of the geometry, Tw - Tb = 17 q h / (35 k), Nu_Dh = 140/17;
        * the shape of the profile across the gap, independently of q, k and the
          flow rate, is (Tb - T_centre)/(Tw - Tb) = 39/136 - a conductive field
          would give 0 there, so it is the convection-diffusion balance itself.
      Both are read at an interior station, where the axial conduction of the
      entrance cannot reach.

usage: check_energy.py <solver-binary> <case.dict> <mode>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile


# ------------------------------------------------------------------- case input
def case_constants(path):
    """Every number the checks are built from, read out of the case itself so a
    case edit cannot silently move the target."""
    txt = open(path).read()

    def num(key, default=None):
        m = re.search(r'(?<!\w)%s\s+(-?[0-9.eE+-]+)\s*;' % key, txt)
        return float(m.group(1)) if m else default

    def triple(key, source=None):
        m = re.search(r'(?<!\w)%s\s*\(([^)]*)\)' % key, source or txt)
        return [float(v) for v in m.group(1).split()] if m else None

    energy = re.search(r'energy\s*\{(.*?)\n\}', txt, re.S).group(1)
    block = re.search(r'block\s*\{(.*?)\n\}', txt, re.S).group(1)
    patches = dict(re.findall(r'(\w+)\s*\{([^}]*)\}',
                              re.search(r'patches\s*\{(.*?)\n\}', txt, re.S).group(1)))

    def e(key, default=None):
        m = re.search(r'(?<!\w)%s\s+(-?[0-9.eE+-]+)' % key, energy)
        return float(m.group(1)) if m else default

    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', block).group(1).split()]
    return {
        'nx': cells[0], 'ny': cells[1], 'nz': cells[2],
        'lo': triple('min', block), 'hi': triple('max', block),
        'nu': num('nu'),
        'vtu': re.search(r'vtk\s+([A-Za-z0-9_.-]+);', txt).group(1),
        'rho': e('rho', 1000.0), 'cp': e('cp', 4180.0), 'kappa': e('kappa', 0.6),
        'source': e('volumetricSource', 0.0),
        'beta': e('expansion', 0.0), 'tRef': e('referenceT', 0.0),
        'gravity': triple('gravity', energy) or [0.0, 0.0, 0.0],
        'initialT': e('initialT', 0.0), 'mode': int(e('initialMode', 0.0)),
        'amp': e('initialAmplitude', 0.0),
        'axis': re.search(r'initialAxis\s+(\w+)', energy).group(1)
                if re.search(r'initialAxis\s+(\w+)', energy) else 'z',
        'dt': num('dt'), 'endTime': num('endTime'), 'force': triple('force', txt),
        'patches': patches,
    }


def patch_num(c, name, key):
    m = re.search(r'(?<!\w)%s\s+(-?[0-9.eE+-]+)\s*;' % key, c['patches'].get(name, ''))
    return float(m.group(1)) if m else None


# ------------------------------------------------------------------ vtu output
def read_vtu(path):
    """Cell centres and cell data as the solver wrote them: the geometry the
    checks measure is the mesh's own, not a re-derivation of it."""
    txt = open(path).read()
    arrays = re.findall(r'<DataArray([^>]*)>(.*?)</DataArray>', txt, re.S)

    def raw(name):
        for head, body in arrays:
            m = re.search(r'Name="([^"]+)"', head)
            if m and m.group(1) == name:
                return body.split()
        return None

    pts = [float(v) for v in arrays[0][1].split()]
    points = list(zip(pts[0::3], pts[1::3], pts[2::3]))
    conn, offs = [int(v) for v in raw('connectivity')], [int(v) for v in raw('offsets')]
    cent, s = [], 0
    for o in offs:
        cl = [points[v] for v in conn[s:o]]
        s = o
        cent.append(tuple(sum(q[i] for q in cl) / len(cl) for i in range(3)))
    nf = len(cent)
    fields = {}
    for head, body in arrays:
        m = re.search(r'Name="([^"]+)"', head)
        if not m:
            continue
        v = [float(x) for x in body.split()]
        if len(v) == nf:
            fields[m.group(1)] = v
        elif len(v) == 3 * nf:
            fields[m.group(1)] = [tuple(v[3 * i:3 * i + 3]) for i in range(nf)]
    return cent, fields


# ------------------------------------------------------------------- log output
def log_values(out):
    """The numbers the solver prints about the thermal field. They are written at
    setprecision(12), which is what lets the checks below be relative to 1e-6."""
    d = {}
    for k in ('in', 'source', 'enthalpyOut', 'balance'):
        d[k] = float(re.search(r'%s (\S+)' % k, out).group(1))
    for m in re.finditer(r'heatPatch (\S+) Q (\S+)', out):
        d['Q:' + m.group(1)] = float(m.group(2))
    d['history'] = [(float(a), float(b), float(c)) for a, b, c in
                    re.findall(r'energy t (\S+) Tmin (\S+) Tmax (\S+)', out)]
    m = re.search(r'energy: Tmin (\S+) Tmax (\S+)', out)
    d['Tmin'], d['Tmax'] = float(m.group(1)), float(m.group(2))
    s = re.search(r'piso done: steps (\d+)', out)
    d['steps'] = int(s.group(1)) if s else None
    return d


def expect_scale(value, target, rel, failures, what):
    if abs(value - target) > rel * max(abs(target), 1e-30):
        failures.append('%s: got %.9g, expected %.9g (+-%.3g%%)'
                        % (what, value, target, 100 * rel))


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    solver, case_path, mode = sys.argv[1], sys.argv[2], sys.argv[3]
    c = case_constants(case_path)
    work = tempfile.mkdtemp(prefix='chflow_energy_')
    try:
        shutil.copy(os.path.abspath(case_path), os.path.join(work, 'case.dict'))
        proc = subprocess.run([solver, 'case.dict'], cwd=work,
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit('%s exited %d\n%s\n%s'
                             % (solver, proc.returncode, proc.stdout[-3000:], proc.stderr))
        out = proc.stdout + proc.stderr
        if 'NOT converged' in out:
            raise SystemExit('the case did not converge:\n' + out[-2000:])
        log = log_values(out)
        cent, f = read_vtu(os.path.join(work, c['vtu']))
        failures, notes = [], []
        check_common(c, log, f, failures, notes)
        {'fixed': check_fixed, 'flux': check_flux, 'source': check_source,
         'rest': check_rest, 'decay': check_decay, 'channel': check_channel}[mode](
            c, log, cent, f, failures, notes)
        for n in notes:
            print('  ' + n)
        if failures:
            raise SystemExit('energy %s check FAILED:\n  %s' % (mode, '\n  '.join(failures)))
        print('energy %s check passed' % mode)
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------- all modes
def check_common(c, log, f, failures, notes):
    """Whatever any thermal case owes: a finite field and a closed heat balance.
    The balance is the discrete statement that the assembled operators telescope -
    an internal face cancels against the neighbour's row - so what survives is
    exactly what the patches feed and the sources generate. It holds to the
    linear-solver tolerance, not to a tolerance chosen here."""
    if any(not math.isfinite(t) for t in f['T']):
        failures.append('temperature is not finite')
    scale = max(abs(log['in']), abs(log['source']), abs(log['enthalpyOut']), 1e-12)
    notes.append('balance %.3e W against a %.6g W scale' % (log['balance'], scale))
    if abs(log['balance']) > 1e-6 * scale:
        failures.append('heat balance did not close: %.3e W against a %.3g W scale'
                        % (log['balance'], scale))


# ------------------------------------------------------------------- conduction
def check_fixed(c, log, cent, f, failures, notes):
    L = c['hi'][0] - c['lo'][0]
    hot, cold = patch_num(c, 'hot', 'T'), patch_num(c, 'cold', 'T')
    area = (c['hi'][1] - c['lo'][1]) * (c['hi'][2] - c['lo'][2])
    line = lambda x: hot + (cold - hot) * x / L
    err = max(abs(t - line(p[0])) for p, t in zip(cent, f['T']))
    q = c['kappa'] * area * (hot - cold) / L
    notes.append('profile off the linear solution by %.3e K; wall rate %.9g W '
                 'against k A dT/dx = %.9g' % (err, log['Q:hot'], q))
    if err > 1e-7:
        failures.append('profile off the linear solution by %.3e K' % err)
    expect_scale(log['Q:hot'], q, 1e-8, failures, 'hot face heat rate')
    expect_scale(log['Q:cold'], -q, 1e-8, failures, 'cold face heat rate')
    # the extrema of a linear profile on cell centres sit half a cell in
    inset = (cold - hot) / (2.0 * c['nx'])
    expect_scale(log['Tmax'], hot + inset, 1e-9, failures, 'hottest cell')
    expect_scale(log['Tmin'], cold - inset, 1e-9, failures, 'coldest cell')


def check_flux(c, log, cent, f, failures, notes):
    L = c['hi'][0] - c['lo'][0]
    q, cold = patch_num(c, 'feed', 'heatFlux'), patch_num(c, 'cold', 'T')
    k = c['kappa']
    area = (c['hi'][1] - c['lo'][1]) * (c['hi'][2] - c['lo'][2])
    line = lambda x: cold + q * (L - x) / k
    err = max(abs(t - line(p[0])) for p, t in zip(cent, f['T']))
    notes.append('profile off the Fourier solution by %.3e K; hottest face reads '
                 '%.9g K, Fourier gives %.9g K' % (err, log['Tmax'], cold + q * L / k))
    if err > 1e-7:
        failures.append('profile off the Fourier solution by %.3e K' % err)
    expect_scale(log['Q:feed'], q * area, 1e-8, failures, 'heated face heat rate')
    expect_scale(log['Q:cold'], -q * area, 1e-8, failures, 'isothermal face heat rate')
    # a fixedGradient face has to report the flux it was prescribed, with the
    # sign that says the fluid gained it: that is the convention of the key
    if log['Q:feed'] <= 0:
        failures.append('a positive heatFlux cooled the fluid: sign convention broken')


def check_source(c, log, cent, f, failures, notes):
    L = c['hi'][0] - c['lo'][0]
    wall, k, s = patch_num(c, 'left', 'T'), c['kappa'], c['source']
    vol = L * (c['hi'][1] - c['lo'][1]) * (c['hi'][2] - c['lo'][2])
    ts = sorted(zip([p[0] for p in cent], f['T']))
    sym = max(abs(a[1] - b[1]) for a, b in zip(ts, ts[::-1]))
    generated = s * vol
    peak = max(t for _, t in ts)
    rise = peak - wall
    continuous = s * L * L / (8.0 * k)
    para = lambda x: wall + s / (2.0 * k) * x * (L - x)
    shape = max(abs(t - para(x)) for x, t in ts)
    notes.append('mirror symmetry %.3e K over a %.4g K rise' % (sym, rise))
    notes.append('peak %.9g K: the rise against the continuous S L^2/(8k) = %.9g K'
                 % (peak, wall + continuous))
    notes.append('profile %.3e K off the continuous parabola (%.2f%% of the rise): '
                 'the half-cell the scheme imposes the wall value from'
                 % (shape, 100.0 * shape / rise))
    if sym > 1e-9:
        failures.append('the discrete field is not mirror symmetric: %.3e K' % sym)
    expect_scale(log['source'], generated, 1e-9, failures, 'reported source power')
    for name in ('left', 'right'):
        expect_scale(log['Q:' + name], -generated / 2.0, 1e-8,
                     failures, '%s face heat rate, half the generation' % name)
    expect_scale(log['Q:left'], log['Q:right'], 1e-9, failures, 'the two face rates')
    expect_scale(rise, continuous, 0.005, failures, 'centre rise above the wall')
    if shape > 0.01 * rise:
        failures.append('profile shape %.3e K off the parabola, %.2f%% of the rise'
                        % (shape, 100.0 * shape / rise))


# --------------------------------------------------------------------- buoyancy
def check_rest(c, log, cent, f, failures, notes):
    g, H = c['gravity'][2], c['hi'][2] - c['lo'][2]
    fz = -c['beta'] * (c['initialT'] - c['tRef']) * g
    unbalanced = abs(fz) * H * H / c['nu']
    uMax = max(math.sqrt(sum(v * v for v in u)) for u in f['U'])
    pts = sorted((p[2], pr) for p, pr in zip(cent, f['p']))
    zs, ps = [a for a, _ in pts], [b for _, b in pts]
    n = len(zs)
    sz, sp = sum(zs), sum(ps)
    slope = (n * sum(z * q for z, q in zip(zs, ps)) - sz * sp) / \
        (n * sum(z * z for z in zs) - sz * sz)
    intercept = (sp - slope * sz) / n
    resid = max(abs(q - (slope * z + intercept)) for z, q in zip(zs, ps))
    notes.append('|U| max %.3e m/s; the same force with nothing to answer it would '
                 'reach %.3e' % (uMax, unbalanced))
    notes.append('dp/dz %.6f against f_z %.6f, %.3e off a straight line over '
                 'f H = %.4g' % (slope, fz, resid, abs(fz) * H))
    if uMax > 0.05 * unbalanced:
        failures.append('the box circulates at %.3e m/s: the pressure does not '
                        'answer the buoyancy (unanswered scale %.3e)' % (uMax, unbalanced))
    expect_scale(slope, fz, 0.05, failures, 'hydrostatic pressure gradient')
    if slope * fz <= 0:
        failures.append('pressure gradient %.6g does not oppose buoyancy %.6g: '
                        'the force enters with the wrong sign' % (slope, fz))
    if resid > 0.05 * abs(fz) * H:
        failures.append('pressure is not the linear hydrostatic: %.3e off' % resid)
    if log['Tmax'] - log['Tmin'] > 1e-5:
        failures.append('a uniform offset from the reference temperature evolved: '
                        'spread %.3e K' % (log['Tmax'] - log['Tmin']))
    expect_scale(0.5 * (log['Tmin'] + log['Tmax']), c['initialT'], 1e-8,
                 failures, 'temperature level held')


# -------------------------------------------------------------------- transient
def check_decay(c, log, cent, f, failures, notes):
    if c['axis'] != 'x':
        raise SystemExit('check_decay reads the mode along x; the case puts it on ' + c['axis'])
    alpha = c['kappa'] / (c['rho'] * c['cp'])
    n, L, m = c['nx'], c['hi'][0] - c['lo'][0], c['mode']
    dx = L / n
    lam = 4.0 * math.sin(m * math.pi / (2.0 * n)) ** 2 / dx ** 2
    k2 = (m * math.pi / L) ** 2
    notes.append('discrete eigenvalue %.6f against continuous k^2 = %.6f (alpha %.4g)'
                 % (lam, k2, alpha))

    def decayed(t):
        """The solver's own stepping loop: dt until the remaining time is shorter."""
        tn, factor = 0.0, 1.0
        while tn < t - 1e-12 * max(1.0, t):
            h = min(c['dt'], t - tn)
            tn += h
            factor /= (1.0 + h * alpha * lam)
        return factor

    total = decayed(c['endTime'])
    steps = 0
    tn = 0.0
    while tn < c['endTime'] - 1e-12 * max(1.0, c['endTime']):
        tn += min(c['dt'], c['endTime'] - tn)
        steps += 1
    continuous = math.exp(-alpha * k2 * c['endTime'])
    notes.append('%d steps (solver: %s), discrete factor %.9f, continuous %.9f: '
                 'the mode decays %.3f%% slower than the exact PDE'
                 % (steps, log['steps'], total, continuous, 100.0 * (total / continuous - 1.0)))
    if log['steps'] is not None and steps != log['steps']:
        failures.append('replayed %d steps, the solver took %d' % (steps, log['steps']))

    peak0 = c['amp'] * math.cos(m * math.pi / (2.0 * n))
    worst = 0.0
    for ti, tmin, tmax in log['history']:
        a = peak0 * decayed(ti)
        worst = max(worst, abs(tmax - a), abs(tmin + a))
    worst /= peak0
    notes.append('amplitude history off the discrete eigenmode by %.3e (relative, '
                 'over %d reports)' % (worst, len(log['history'])))
    if worst > 1e-8:
        failures.append('decay off the discrete eigenmode by %.3e' % worst)
    expect_scale(log['Tmax'], peak0 * total, 1e-8, failures, 'final amplitude')

    mean = sum(f['T']) / len(f['T'])
    notes.append('discrete mean %.3e (initialT %.4g)' % (mean, c['initialT']))
    if abs(mean - c['initialT']) > 1e-12:
        failures.append('the adiabatic box lost or gained enthalpy: mean %.9g' % mean)
    if abs(log['enthalpyOut']) > 1e-12:
        failures.append('a closed box reports an enthalpy flux of %.3e W' % log['enthalpyOut'])


# -------------------------------------------------------------------- convection
def check_channel(c, log, cent, f, failures, notes):
    T, U = f['T'], f['U']
    xs = sorted(set(round(p[0], 9) for p in cent))
    ys = sorted(set(round(p[1], 9) for p in cent))
    dy = ys[1] - ys[0]
    depth = c['hi'][2] - c['lo'][2]
    H = c['hi'][1] - c['lo'][1]
    h = H / 2.0
    idx = {}
    for i, p in enumerate(cent):
        ix = min(range(len(xs)), key=lambda j: abs(xs[j] - p[0]))
        iy = min(range(len(ys)), key=lambda j: abs(ys[j] - p[1]))
        idx[(ix, iy)] = i
    column = lambda ix: [idx[(ix, iy)] for iy in range(len(ys))]

    phi = [sum(U[i][0] for i in column(ix)) * dy * depth for ix in range(len(xs))]
    if max(phi) <= 0:
        failures.append('no through-flow: the channel is not driven')
        return
    notes.append('flow rate %.6g m3/s, its spread along the channel %.3e'
                 % (max(phi), (max(phi) - min(phi)) / max(phi)))

    tb = []
    for ix in range(len(xs)):
        cells = column(ix)
        num = sum(U[i][0] * T[i] for i in cells)
        tb.append(num / sum(U[i][0] for i in cells))

    # 1. the discrete convective identity, with no analytic input: the enthalpy
    #    the mixing cup gains from the prescribed inlet temperature to the last
    #    column is the heat every surface fed it.
    tin = patch_num(c, 'left', 'T')
    fed = sum(v for k, v in log.items() if k.startswith('Q:'))
    gained = c['rho'] * c['cp'] * phi[0] * (tb[-1] - tin)
    notes.append('mixing cup gained %.9g W, the patches fed %.9g W' % (gained, fed))
    expect_scale(gained, fed, 1e-6, failures, 'enthalpy carried against heat fed')

    # 2 and 3. the fully developed shape at an interior station
    q = patch_num(c, 'bottom', 'heatFlux')
    if abs(patch_num(c, 'top', 'heatFlux') - q) > 1e-12:
        failures.append('the two walls do not carry the same flux')
    k = c['kappa']
    ix = len(xs) // 2
    cells = column(ix)
    # a fixedGradient face value is T_cell + (q/k) dy/2, which is the wall
    # temperature the scheme actually imposes
    tw = [T[idx[(ix, 0)]] + (q / k) * (dy / 2.0),
          T[idx[(ix, len(ys) - 1)]] + (q / k) * (dy / 2.0)]
    offset = 0.5 * (tw[0] + tw[1]) - tb[ix]
    asym = abs(tw[0] - tb[ix] - (tw[1] - tb[ix]))
    centre = T[idx[(ix, len(ys) // 2)]] if len(ys) % 2 else \
        0.5 * (T[idx[(ix, len(ys) // 2 - 1)]] + T[idx[(ix, len(ys) // 2)]])
    shape = (tb[ix] - centre) / offset
    umax = max(U[i][0] for i in cells)
    notes.append('wall-to-bulk offset %.9g K against 17 q h /(35 k) = %.9g; '
                 'Nu_Dh %.4f against 140/17 = %.4f'
                 % (offset, 17.0 * q * h / (35.0 * k), q * 4.0 * h / (k * offset),
                    140.0 / 17.0))
    notes.append('asymmetry between the two walls %.3e K; shape ratio %.6f against '
                 '39/136 = %.6f' % (asym, shape, 39.0 / 136.0))
    expect_scale(offset, 17.0 * q * h / (35.0 * k), 0.02,
                 failures, 'wall-to-bulk temperature difference')
    expect_scale(shape, 39.0 / 136.0, 0.02, failures, 'cross-sectional shape ratio')
    if asym > 1e-9:
        failures.append('the two identical walls give different offsets: %.3e K' % asym)
    # the field the temperature rides on
    expect_scale(umax, c['force'][0] * H * H / (8.0 * c['nu']), 0.02,
                 failures, 'centreline velocity')


if __name__ == '__main__':
    main()
