#!/usr/bin/env python3
"""Verify that a coupled edge's rheology reaches the 3D detailed region.

The case is the shipped steady coupled network with ONE change: pipe X (the
edge handed to the 3D block) declares a generalized-Newtonian fluid with the
same three keys a wall-law pipe uses. The 1D row must not evaluate it -- its
loss is the partner's statement R*q -- and the 3D block must assemble its
momentum diffusion with the per-face apparent viscosity of THAT law.

Five arms are generated from one case dict in a temporary directory:

    plain    no declaration                  (the constant-nu equation)
    newton   n=1, K=mu, tau0=0              (declared, still Newtonian)
    thin     n=0.5, K=1.0
    thick    n=1.5, K=0.00913
    hb       n=0.8, K=0.3, tau0=0.5        (Herschel-Bulkley)

Families of statements:

  1. bit collapse        - plain and newton produce byte-identical solver logs
                           (every rate, drop, sweep count, interface history),
                           although only newton goes through the per-face path
  2. constitutive law    - the printed wall viscosity/shear/stress are rebuilt
                           from K*gdot^(n-1) + tau0/gdot by the checker
  3. discrete momentum   - |probeSlope|*H/2 equals the wall shear stress
  4. closed form         - the block's drop is the plane-channel integral of
                           the SAME law at the rate the interface delivered
  5. interface identity  - each pass's residual is |R*q - dpBlock|/dpBlock,
                           every seed converges to one operating point
  6. no double charge    - the 1D report still uses rho*v*d/mu for Re and the
                           parallel branches still share their drop
  7. transient path      - the same declaration on a transient coupled edge
                           reaches the same fluid and the same operating point
  8. rejection boundary  - a missing K, negative n / tau0, and RANS + rheology
                           are still refused with a named reason

usage: check_coupled_rheology.py <steady-app> <steady-case.dict>
                                 <transient-app> <transient-case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

RHO, MU = 800.0, 0.1
H, W, L, D = 0.02, 0.02, 0.5, 0.02   # duct height, width, length, 1D bore
Q_IN = 2.4e-4
SHEAR_FLOOR, REG, RELAX = 20.0, 100.0, 0.15
TOL_INTERFACE = 1e-4

RHEO_LINE = re.compile(r'flowIndex [^;]+; consistency [^;]+;(?: yieldStress [^;]+;)?')
GN_BLOCK = ("    generalizedNewtonian {\n"
            "        shearFloor    20.0;\n"
            "        regularization 100.0;\n"
            "        relaxation    0.15;\n"
            "    }\n")
TURB_BLOCK = ("    turbulence {\n"
              "        model kOmegaSST;\n"
              "        intensity 0.05;\n"
              "        lengthScale 1.4e-3;\n"
              "        relaxation 0.7;\n"
              "    }\n")

ARMS = [
    ('plain',  '',                                                       1.0, MU,      0.0),
    ('newton', 'flowIndex 1.0; consistency 0.1;',                       1.0, MU,      0.0),
    ('thin',   'flowIndex 0.5; consistency 1.0;',                       0.5, 1.0,     0.0),
    ('thick',  'flowIndex 1.5; consistency 0.00913;',                   1.5, 0.00913, 0.0),
    ('hb',     'flowIndex 0.8; consistency 0.3; yieldStress 0.5;',      0.8, 0.3,     0.5),
]

CHECKS = 0
FAILURES = []


def check(cond, label, detail=''):
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILURES.append('%s %s' % (label, detail))
        print('FAIL', label, detail)


def close(a, b, rel, label, detail=''):
    scale = max(abs(a), abs(b), 1e-30)
    check(abs(a - b) <= rel * scale, label,
          '%s: %.12g vs %.12g (rel %.4g)' % (detail, a, b, abs(a - b) / scale))


def value(line, key):
    m = re.search(r'(?:^|\s)' + re.escape(key) + r'\s+([-+0-9.eE]+)', line)
    if m is None:
        raise SystemExit("log line has no '%s': %s" % (key, line))
    return float(m.group(1))


def tagged(out, prefix):
    return [l for l in out.splitlines() if l.startswith('[chflow] ' + prefix)]


def run(app, case_text, root, name, expect_ok=True, timeout=180):
    work = os.path.join(root, name)
    os.makedirs(work, exist_ok=True)
    path = os.path.join(work, 'case.dict')
    with open(path, 'w') as f:
        f.write(case_text)
    try:
        r = subprocess.run([app, path], capture_output=True, text=True,
                           cwd=work, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise SystemExit('%s: app timed out after %s s' % (name, timeout))
    out = r.stdout + r.stderr
    if expect_ok and r.returncode != 0:
        print(out[-3000:])
        raise SystemExit('%s: app failed (rc=%d)' % (name, r.returncode))
    return out, r.returncode


def arm_case(base, clause):
    return RHEO_LINE.sub(clause, base, count=1)


def canonical(out):
    """Solver log without the lines only a declared fluid produces."""
    return [l for l in out.splitlines() if 'block fluid' not in l]


def apparent_mu(n, K, tau0, gdot):
    g = max(gdot, SHEAR_FLOOR)
    mu = K * g ** (n - 1.0)
    if tau0 > 0:
        mu += tau0 * (-math.expm1(-REG * g)) / g
    return mu


def closed_dp(n, K, tau0, q):
    """Plane-channel developed drop of tau = tau0 + K*gdot^n, inverted in the
    pressure gradient. The plug region contributes its own integral; tau0=0
    collapses to 2*W*n/(2n+1)*(G/K)^(1/n)*(H/2)^((2n+1)/n) = q."""
    h = H / 2.0

    def q_of(G):
        T = G * h - tau0
        if T <= 0:
            return 0.0
        a = K ** (-1.0 / n) / (G * G)
        first = n / (2.0 * n + 1.0) * T ** ((2.0 * n + 1.0) / n)
        plug = n / (n + 1.0) * tau0 * T ** ((n + 1.0) / n)
        return 2.0 * W * a * (first + plug)

    lo, hi = 1e-9, 1e9
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if q_of(mid) < q:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi) * L


def stats(out, name):
    coupled = tagged(out, 'coupled seed ')
    blocks = tagged(out, 'block seed ')
    fluids = tagged(out, 'block fluid seed ')
    interfaces = tagged(out, 'interface seed ')
    networks = tagged(out, 'network seed ')
    finals = tagged(out, 'pipe final ')
    check(len(coupled) == 3, name + ': three seeds', str(len(coupled)))
    check(len(blocks) == 3 and len(networks) == 3,
          name + ': one block/network report per seed')
    check(all('(converged)' in l for l in coupled),
          name + ': every seed converged')
    check('NOT converged' not in out, name + ': no stalled interface')

    qs = [value(l, 'q') for l in coupled]
    dps = [value(l, 'dp') for l in coupled]
    spread = (max(qs) - min(qs)) / max(qs)
    check(spread < 1e-4, name + ': seed-independent rate', 'spread %.3g' % spread)

    # 5. interface identity, on every pass of every seed
    for l in interfaces:
        R, q, dp = value(l, 'R'), value(l, 'q'), value(l, 'dpBlock')
        res = value(l, 'residual')
        close(abs(R * q - dp) / dp, res, 2e-8,
              name + ': residual is the cross-solver gap', l[:70])
    for s, net in enumerate(networks):
        it = [l for l in interfaces if l.startswith('[chflow] interface seed %d ' % s)]
        close(value(net, 'dpNetwork'), value(it[-1], 'R') * value(it[-1], 'q'),
              2e-9, name + ': network charge == R*q')
        check(value(it[-1], 'residual') <= TOL_INTERFACE,
              name + ': final residual under tolerance')

    # block conservation and its reported drop belong to the same solve
    for s, l in enumerate(blocks):
        close(value(l, 'fluxIn'), value(l, 'fluxOut'), 1e-8,
              name + ': block conservation')
        close(value(l, 'dpPa'), dps[s], 1e-12,
              name + ': reported drop == interface drop')
        check(value(l, 'blockConverged') == 1.0,
              name + ': block solve converged')

    # 6. the 1D row did not read the fluid: Re is still rho*v*d/mu, and the
    # parallel branches still share exactly one drop
    px = [l for l in finals if 'pipe final X ' in l]
    by = [l for l in finals if 'pipe final By ' in l]
    lead = [l for l in finals if 'pipe final A ' in l or 'pipe final B ' in l]
    check(len(px) == 1 and len(by) == 1 and len(lead) == 2,
          name + ': final network report present')
    if px and by and len(lead) == 2:
        vx = value(px[0], 'v')
        close(value(px[0], 'Re'), RHO * vx * D / MU, 2e-9,
              name + ': 1D Re still uses the material mu')
        close(value(px[0], 'dp'), value(by[0], 'dp'), 1e-12,
              name + ': parallel branches share the drop')
        for l in lead:
            close(value(l, 'Q'), Q_IN, 1e-9, name + ': lead carries the injection')

    # 2/3/4 on the last seed's fluid report
    slope = abs(value(blocks[-1], 'probeSlope'))
    q_last, dp_last = qs[-1], dps[-1]
    fluid = None
    if fluids:
        f = fluids[-1]
        n_f, K_f, tau0_f = value(f, 'n'), value(f, 'K'), value(f, 'tau0')
        nu_w, gdot_w, tau_w = (value(f, 'nuWall'), value(f, 'gdotWall'),
                               value(f, 'tauWall'))
        fluid = (n_f, K_f, tau0_f, nu_w, gdot_w, tau_w)
        mu_pred = apparent_mu(n_f, K_f, tau0_f, gdot_w)
        close(nu_w, mu_pred / RHO, 0.06, name + ': wall viscosity is the law',
              'nuWall')
        close(tau_w, RHO * nu_w * gdot_w, 0.08, name + ': tau_w = rho nu_w gdot_w')
        close(tau_w, K_f * gdot_w ** n_f + tau0_f, 0.10,
              name + ': tau_w = K gdot^n + tau0')
        close(slope * H / 2.0, tau_w, 0.09,
              name + ': developed gradient carries the wall shear')
    return dict(q=q_last, dp=dp_last, slope=slope, fluid=fluid,
                qs=qs, dps=dps, interfaces=interfaces, fluids=fluids)


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    steady_app, steady_case, transient_app, transient_case = sys.argv[1:]
    base = open(steady_case).read()
    tbase = open(transient_case).read()
    check(RHEO_LINE.search(base) is not None, 'steady case declares a fluid')
    check('coupledResistance' in tbase, 'transient case has a coupled edge')

    with tempfile.TemporaryDirectory(prefix='cf_rheology_') as root:
        outputs, results = {}, {}
        for name, clause, n, K, tau0 in ARMS:
            out, _ = run(steady_app, arm_case(base, clause), root, name)
            outputs[name] = out
            results[name] = stats(out, name)

            # 4. developed closed form of the declared law
            q, dp = results[name]['q'], results[name]['dp']
            dp_c = closed_dp(n, K, tau0, q)
            ratio = dp / dp_c
            print('%-8s q=%.10g dp=%12.6f closed=%12.6f ratio=%.5f' %
                  (name, q, dp, dp_c, ratio))
            check(0.95 < ratio < 1.13, name + ': drop follows its own closed form',
                  'ratio %.5f' % ratio)

        # 1. a Newtonian declaration is a change of code path, not of a number
        a, b = canonical(outputs['plain']), canonical(outputs['newton'])
        check(a == b, 'declared-but-Newtonian is bit-identical',
              'first difference: ' + str(next(
                  ((x, y) for x, y in zip(a, b) if x != y), ('length', (len(a), len(b))))))
        check(not any('block fluid' in l for l in outputs['plain'].splitlines()),
              'plain arm has no viscosity model')
        check(any('block fluid seed' in l for l in outputs['newton'].splitlines()),
              'newton arm goes through the per-face viscosity path')
        f = results['newton']['fluid']
        check(f is not None, 'newton arm reports its wall viscosity')
        if f:
            close(f[3], MU / RHO, 1e-13, 'newton wall viscosity == mu/rho')
            close(f[5], MU * f[4], 1e-12, 'newton tau_w == mu*gdot_w')

        # the three genuinely non-Newtonian arms must not be the Newtonian one
        plain_dp = results['plain']['dp']
        for name in ('thin', 'thick', 'hb'):
            d = abs(results[name]['dp'] - plain_dp) / plain_dp
            check(d > 0.01, name + ': the fluid moves the operating point',
                  'relative dp change %.4g' % d)
        check(results['hb']['dp'] > results['plain']['dp'],
              'yield stress raises the block drop')

        # 7. the same declaration on a transient coupled edge
        tcase = tbase.replace('coupledResistance 3.75e6; }',
                              'coupledResistance 3.75e6; '
                              'flowIndex 0.5; consistency 1.0; }')
        tcase = tcase.replace('detail {\n', 'detail {\n' + GN_BLOCK)
        tout, _ = run(transient_app, tcase, root, 'transient', timeout=240)
        tf = tagged(tout, 'transient block fluid ')
        check(len(tf) >= 1, 'transient block reports the fluid')
        if tf:
            check(abs(value(tf[0], 'n') - 0.5) < 1e-12 and
                  abs(value(tf[0], 'K') - 1.0) < 1e-12,
                  'transient block received n=0.5 K=1.0')
            thin_f = results['thin']['fluid']
            if thin_f:
                close(value(tf[0], 'nuWall'), thin_f[3], 0.06,
                      'fluid wall viscosity matches the steady arm')
        steps = tagged(tout, 'transient step ')
        check(len(steps) > 0, 'transient produced time levels')
        for l in steps:
            close(value(l, 'dpNetwork'), value(l, 'R') * value(l, 'q'), 2e-9,
                  'transient momentum row charges R*q', l[:70])
        last = tagged(tout, 'transient coupled seed ')
        check(len(last) >= 1, 'transient interface settled')
        if last:
            close(value(last[-1], 'q'), results['thin']['q'], 0.02,
                  'transient settles on the steady operating point')
            close(value(last[-1], 'dp'), results['thin']['dp'], 0.03,
                  'transient settles on the steady drop')

        # 8. the guards that must still fire
        bad = [
            # The 3D viscosity model is built before the 1D network is loaded,
            # so these two are refused by the constitutive law itself; the 1D
            # parser carries the same checks for an uncoupled pipe.
            ('missing_K', arm_case(base, 'flowIndex 0.5;'),
             'needs a positive consistency K'),
            ('negative_n', arm_case(base, 'flowIndex -0.5; consistency 1.0;'),
             'negative'),
            ('negative_tau0',
             arm_case(base, 'flowIndex 0.8; consistency 0.3; yieldStress -1.0;'),
             'negative flow index or yield stress'),
            ('rans_plus_fluid',
             arm_case(base, ARMS[2][1]).replace('detail {\n',
                                                'detail {\n' + TURB_BLOCK),
             'cannot both supply'),
        ]
        for name, text, expected in bad:
            out, rc = run(steady_app, text, root, name, expect_ok=False)
            check(rc != 0, name + ': refused', 'rc=%d' % rc)
            check(expected in out, name + ': names the reason', expected)

    print('\n%d checks, %d failures' % (CHECKS, len(FAILURES)))
    if FAILURES:
        for f in FAILURES[:20]:
            print('  -', f)
        raise SystemExit(1)
    print('OK: coupled-edge rheology drives the 3D wall law without touching the 1D row')


if __name__ == '__main__':
    main()
