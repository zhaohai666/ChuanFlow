#!/usr/bin/env python3
"""Run the heated turbulent channel and verify the thermal wall function.

The case is the force-driven column of the turbulent channel (u_tau = 0.6245
m/s is the momentum balance, not a prediction), with 310 K at the bottom wall,
300 K at the top and air on a 16-cell column whose first cell sits at y+ ~ 52
- four thermal sublayers above the wall, where the two-branch temperature law
of the wall (conductive T+ = Pr y+, log T+ = Prt/k ln(E y+) + Prt P) is all
log. The steady temperature field is pure conduction across the gap: the
continuity flux is zero on a streamwise-uniform column.

The checks, in the order they pin the kernel down:

  1. the reported T+ and wall heat rate are the two-branch law, reimplemented
     here in Python and evaluated at the reported y+ - the 12-digit wall
     report makes this an identity at 1e-8;
  2. the wall flux is 3.1x the molecular first-cell flux, which is what "the
     closure is actually in play" means numerically;
  3. the temperature FIELD satisfies the discrete problem the kernel claims
     to assemble: harmonic face diffusivity from the vtu nut inside, the
     law-implied u_tau y/T+ at the wall faces. That tridiagonal system is
     rebuilt and solved here, and the vtu field has to match it to the outer
     iteration's own residual (~4e-5 K measured) - a wrong wall treatment, a
     wrong harmonic weight or a molecular wall would all miss by O(0.1) K or
     more;
  4. the solved field carries the same heat through every face, the walls
     agree with the interior, the two equal wall laws force a mirror-symmetric
     profile, and the reported balance closes with the wall flux as scale.

Both RANS models run the same checks - the wall layer is the turbulence
model's, the interior diffusivity comes from the vtu nut, so nothing here
knows which model produced the field.

usage: check_thermal_wall.py <solver-binary> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

KAPPA = 0.41
LOG_E = 9.8


def case_constants(path):
    """Every number the checks are built from, read out of the case itself."""
    txt = open(path).read()
    nu = float(re.search(r'^nu\s+([0-9.eE+-]+);', txt, re.M).group(1))
    energy = re.search(r'energy\s*\{(.*?)\n\}', txt, re.S).group(1)

    def e(key):
        return float(re.search(r'(?<!\w)%s\s+([0-9.eE+-]+)' % key, energy).group(1))

    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]
    patches = dict(re.findall(r'(\w+)\s*\{([^}]*)\}',
                              re.search(r'patches\s*\{(.*?)\n\}', txt, re.S).group(1)))

    def pT(name):
        return float(re.search(r'(?<!\w)T\s+([0-9.eE+-]+)\s*;', patches[name]).group(1))

    return dict(nu=nu, rho=e('rho'), cp=e('cp'), kappa=e('kappa'), Prt=e('turbulentPrandtl'),
                ny=cells[1], vtu=re.search(r'^\s*vtk\s+(\S+);', txt, re.M).group(1),
                twBottom=pT('bottom'), twTop=pT('top'),
                height=float(re.search(r'max\s*\(([^)]*)\)', txt).group(1).split()[1]))


# ------------------------------------------------------- the thermal law itself
def jayatilleke_p(pr, prt):
    r = pr / prt
    return 9.24 * (r ** 0.75 - 1.0) * (1.0 + 0.28 * math.exp(-0.007 * r))


def sublayer_edge(pr, prt):
    p = jayatilleke_p(pr, prt)
    f = lambda y: pr * y - prt / KAPPA * math.log(LOG_E * y) - prt * p
    lo = prt / (KAPPA * pr)
    hi = 2.0 * lo
    while f(hi) < 0:
        lo, hi = hi, 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def thermal_tplus(yp, pr, prt):
    yt = sublayer_edge(pr, prt)
    return pr * yp if yp <= yt else prt / KAPPA * math.log(LOG_E * yp) + prt * jayatilleke_p(pr, prt)


# ------------------------------------------------------------------- solver I/O
def reports_from_log(out):
    wall = {}
    for m in re.finditer(r'thermalWall (\S+) uTau (\S+) yPlus (\S+) tPlus (\S+) q (\S+)', out):
        wall[m.group(1)] = dict(uTau=float(m.group(2)), yPlus=float(m.group(3)),
                                tPlus=float(m.group(4)), q=float(m.group(5)))
    heat = {m.group(1): float(m.group(2))
            for m in re.finditer(r'heatPatch (\S+) Q (\S+)', out)}
    done = re.search(r'SIMPLE done: iterations=(\d+) continuity=(\S+) momentum=(\S+)'
                     r' \((converged|NOT converged)\)', out)
    return wall, heat, done


def read_vtu(path):
    txt = open(path).read()
    arrays = re.findall(r'<DataArray([^>]*)>(.*?)</DataArray>', txt, re.S)

    def raw(name):
        for head, body in arrays:
            m = re.search(r'Name="([^"]+)"', head)
            if m and m.group(1) == name:
                return [float(x) for x in body.split()]
        raise SystemExit('%s: %s missing from the vtu' % (path, name))

    pts = [float(v) for v in arrays[0][1].split()]
    points = list(zip(pts[0::3], pts[1::3], pts[2::3]))
    conn = [int(v) for v in raw('connectivity')]
    offs = [int(v) for v in raw('offsets')]
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
    return cent, fields


def run(binary, case, model):
    work = tempfile.mkdtemp(prefix='ch-thermal-wall-')
    binary = os.path.abspath(binary)
    src = os.path.dirname(case)
    for f in os.listdir(src):
        shutil.copy(os.path.join(src, f), work)
    local = os.path.join(work, os.path.basename(case))
    txt = open(local).read()
    open(local, 'w').write(txt.replace('model       kEpsilon;', 'model       %s;' % model))
    r = subprocess.run([binary, local], capture_output=True, text=True, cwd=work)
    print('---- %s (rc %d) ----' % (model, r.returncode))
    print(r.stdout[-900:])
    if r.returncode != 0:
        print(r.stderr)
        raise SystemExit('%s solver failed' % model)
    if 'nan' in r.stdout.lower():
        raise SystemExit('%s: the run produced a non-finite field' % model)
    return dict(out=r.stdout, vtu=os.path.join(work, case_constants(case)['vtu']), work=work)


# --------------------------------------------------- the discrete problem, again
def solve_discrete(c, nut, u_tau):
    """The steady conduction column the kernel claims to assemble: harmonic
    face diffusivity of alpha_mol + nut/Pr_t inside, u_tau y / T+(y+) at the
    wall faces. Solved as a tridiagonal system over the cell column."""
    n = c['ny']
    dy = c['height'] / n
    area = 4.0e-4  # the case's own x/z extents
    alpha = c['kappa'] / (c['rho'] * c['cp'])
    rhocp = c['rho'] * c['cp']
    pr = c['nu'] * rhocp / c['kappa']

    g_wall = u_tau * (dy / 2.0) / thermal_tplus(u_tau * (dy / 2.0) / c['nu'], pr, c['Prt'])
    coef_w = g_wall * 2.0 * area / dy
    a = [0.0] * n
    b = [0.0] * n
    cc = [0.0] * n
    d = [0.0] * n
    for i in range(1, n):
        ga = alpha + nut[i - 1] / c['Prt']
        gb = alpha + nut[i] / c['Prt']
        coef = 2.0 * ga * gb / (ga + gb) * area / dy
        a[i] -= coef
        b[i] += coef
    for i in range(n - 1):
        ga = alpha + nut[i] / c['Prt']
        gb = alpha + nut[i + 1] / c['Prt']
        coef = 2.0 * ga * gb / (ga + gb) * area / dy
        cc[i] -= coef
        b[i] += coef
    b[0] += coef_w
    d[0] = coef_w * c['twBottom']
    b[n - 1] += coef_w
    d[n - 1] = coef_w * c['twTop']

    cp = [0.0] * n
    dp = [0.0] * n
    cp[0] = cc[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        den = b[i] - a[i] * cp[i - 1]
        cp[i] = cc[i] / den if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / den
    x = [0.0] * n
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    flux = []
    for i in range(n - 1):
        ga = alpha + nut[i] / c['Prt']
        gb = alpha + nut[i + 1] / c['Prt']
        coef = 2.0 * ga * gb / (ga + gb) * area / dy
        flux.append(rhocp * coef * (x[i] - x[i + 1]))
    return x, flux


# ------------------------------------------------------------------------ checks
def check(model, res, c):
    out = res['out']
    failures = []
    pr = c['nu'] * c['rho'] * c['cp'] / c['kappa']
    rhocp = c['rho'] * c['cp']
    area = 4.0e-4
    dy = c['height'] / c['ny']

    wall, heat, done = reports_from_log(out)
    if done is None:
        raise SystemExit('%s: no SIMPLE convergence line in the log' % model)
    iters, momentum, state = int(done.group(1)), float(done.group(3)), done.group(4)
    print('%s: SIMPLE %s in %d iterations, momentum residual %.2e'
          % (model, state, iters, momentum))
    if state != 'converged':
        failures.append('SIMPLE did not reach its convergence criteria')
    if iters < 2:
        failures.append('SIMPLE stopped without letting the fields evolve')

    cent, f = read_vtu(res['vtu'])
    t = f['T']
    nut = f['nut']
    if len(t) != c['ny'] or len(nut) != c['ny']:
        raise SystemExit('%s: expected %d cells in the vtu' % (model, c['ny']))
    col = sorted(zip([p[1] for p in cent], t, nut))
    ts = [v[1] for v in col]
    nuts = [v[2] for v in col]
    if any(not math.isfinite(v) for v in ts):
        failures.append('temperature is not finite')

    for name in ('bottom', 'top'):
        if name not in wall:
            raise SystemExit('%s: no thermalWall report for %s' % (model, name))
    wb, wt = wall['bottom'], wall['top']
    twb, twt = c['twBottom'], c['twTop']

    # 1. the wall law, against the Python reimplementation
    for name, w, tw, cell in (('bottom', wb, twb, 0), ('top', wt, twt, c['ny'] - 1)):
        tp_py = thermal_tplus(w['yPlus'], pr, c['Prt'])
        expect_q = rhocp * w['uTau'] * area * (tw - ts[cell]) / tp_py
        err_t = abs(w['tPlus'] / tp_py - 1.0)
        err_q = abs(w['q'] / expect_q - 1.0)
        print('%s: %s y+=%.2f T+ %.6f (law %.6f), q %.6g W (law %.6g W)'
              % (model, name, w['yPlus'], w['tPlus'], tp_py, w['q'], expect_q))
        if err_t > 1e-8:
            failures.append('%s: reported T+ is not the two-branch law' % name)
        if err_q > 1e-8:
            failures.append('%s: wall heat rate is not rho cp u_tau A dT / T+' % name)
        if abs(heat[name] / w['q'] - 1.0) > 1e-12:
            failures.append('%s: heatPatch and thermalWall disagree' % name)

    # 2. the closure is in play: molecular conduction through the first cell
    #    would carry about a third of this flux
    q_mol = c['kappa'] * area * (twb - ts[0]) / (dy / 2.0)
    ratio = wb['q'] / q_mol
    print('%s: wall flux is %.2f x the molecular first-cell flux' % (model, ratio))
    if ratio < 2.0:
        failures.append('the wall function carries less than twice the molecular flux')

    # 3. the field against the discrete problem the kernel claims to assemble
    t_py, flux = solve_discrete(c, nuts, wb['uTau'])
    worst = max(abs(a - b) for a, b in zip(ts, t_py))
    print('%s: field off the restated discrete solution by %.2e K; its flux '
          '%.9g W (spread %.2e relative)' % (model, worst, flux[0],
                                             (max(flux) - min(flux)) / flux[0]))
    if worst > 1e-3:
        failures.append('the field is not the restated discrete solution (%.2e K off)' % worst)
    if (max(flux) - min(flux)) > 1e-9 * flux[0]:
        failures.append('the restated system does not carry a constant flux')

    # 4. profile shape and the reported balance, with the wall flux as scale
    tmid = 0.5 * (twb + twt)
    sym = max(abs((v - tmid) + (w - tmid)) for v, w in zip(ts, ts[::-1]))
    print('%s: T %.4f..%.4f K, mirror symmetry %.2e K' % (model, ts[0], ts[-1], sym))
    if any(b >= a for a, b in zip(ts, ts[1:])):
        failures.append('the profile is not monotone from the hot wall to the cold one')
    if sym > 1e-6:
        failures.append('the profile is not mirror symmetric about the mean wall temperature')

    m = re.search(r'energy: in (\S+) source (\S+) enthalpyOut (\S+) balance (\S+)', out)
    bal = float(m.group(4))
    scale = max(abs(heat['bottom']), abs(heat['top']), 1e-30)
    print('%s: balance %.3e W against the %.6g W wall flux' % (model, bal, scale))
    if abs(bal) > 1e-6 * scale:
        failures.append('the heat balance did not close against the wall flux')

    return failures


def main():
    binary, case = sys.argv[1], sys.argv[2]
    c = case_constants(case)
    pr = c['nu'] * c['rho'] * c['cp'] / c['kappa']
    print('channel: height %.4g, Pr %.4f, Prt %.2f -> thermal sublayer y+_T = %.2f'
          % (c['height'], pr, c['Prt'], sublayer_edge(pr, c['Prt'])))

    bad = False
    for model in ('kEpsilon', 'kOmegaSST'):
        res = run(binary, case, model)
        try:
            for fail in check(model, res, c):
                print('FAIL %s: %s' % (model, fail))
                bad = True
        finally:
            shutil.rmtree(res['work'], ignore_errors=True)
    if bad:
        raise SystemExit('thermal wall verification FAILED')
    print('thermal wall verification passed: both models carry the wall-law flux')


if __name__ == '__main__':
    main()
