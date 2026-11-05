#!/usr/bin/env python3
"""Run the Sod shock tube and verify what a density-based Euler increment has
to deliver.

sodShockTube (chflow-3d-compressible): the canonical Riemann problem, on the
polyhedral kernel, with every boundary an impermeable slip wall. A 1 m pipe of
one-cell square section is split at x = 0.5 into rho 1 / p 1 and rho 0.125 /
p 0.1, both at rest; the diaphragm bursts into a left rarefaction, a contact
and a right shock, and at t = 0.2 none of the fronts has reached a wall. Four
things are owed, and each is checked a different way:

  1. The discrete solution is THE discrete solution. The checker re-runs the
     kernel's own first-order Godunov + HLLC recursion cell by cell in Python
     -- primitives from the conserved totals, per-face HLLC flux with the same
     wave-speed and star-state formulas, wall flux p*Sf, the same scatter and
     update order -- and holds the final written fields to 1e-9. Nothing about
     the physics below can mean anything unless this holds.

  2. Conservation in the form a closed wall-bounded domain actually has it.
     Mass and total energy telescope over the interior faces and come back
     exactly (1e-12 relative); momentum does not, because the walls take the
     pressure force -- the run accumulates the impulse integral
     impulse(t) = integral dt * sum_boundary (p Sf)_x and the checker verifies
     mom(t) + impulse(t) = mom(0) to 1e-12 for every report. Until the fronts
     reach the walls the wall pressures are the initial ones, so the impulse
     has the closed form dt*(pR - pL)*A per step and mom grows linearly at
     (pL - pR)*A -- both are asserted.

  3. The exact Riemann solution (Toro's star state by Newton on the pressure
     function, then sampling) at t = 0.2. The solver is NOT given p*, u* or
     either star density. First-order Godunov keeps uniform plateaus uniform but
     diffuses every front, and the diffusion is set by the scheme's numerical
     viscosity |lambda| dx / 2, so each front's profile is a few sigma_lambda
     wide with sigma_lambda = sqrt(|lambda| dx t / 2). The windows are therefore
     exempted by that physical width and not by a cell count, which is what
     makes the same tolerances mean the same thing on any mesh: far left 1%,
     the star plateaus 3%, the rarefaction compared pointwise against the
     self-similar profile rather than exempted. The fronts themselves are read
     as the profile's 50%-density crossing, interpolated between cell centres,
     and have to sit within half a sigma of the exact position. The exact
     solver is itself anchored to Toro's published Sod star state before it is
     used.

  4. Convergence, over the whole field: the L1 error against the exact sample
     is reported for the case's mesh and for the checker's own refined mesh
     (cells doubled in x, dt halved so the CFL number is the same), and the
     order log2(coarse/fine) is held to [0.4, 1.5] -- the band a first-order
     scheme on a solution with three discontinuities lives in, where the
     smearing of the fronts contributes sqrt(dx) rather than dx.

  5. Positivity: p and rho stay positive in every report (the invariant a
     conservative Godunov scheme is allowed to promise at these CFL numbers),
     and the transverse velocity the 1D setup cannot generate stays exactly
     zero in the written field.

usage: check_euler.py <chflow-3d-compressible> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

GAMMA_TORO_SOD = 1.4


def case_constants(path):
    """The numbers the checks are built from, read out of the case itself so it
    cannot drift away from them. The prose is stripped first: a case may
    describe a key in its header without that being what it asks the solver to
    do."""
    txt = "\n".join(re.sub(r'//.*', '', line) for line in open(path))

    def triple(key):
        return [float(v) for v in re.search(r'%s\s*\(([^)]*)\)' % key, txt).group(1).split()]

    lo, hi = triple('min'), triple('max')
    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]
    td = re.search(r'time\s*\{([^}]*)\}', txt, re.S).group(1)
    ini = re.search(r'initial\s*\{(.*?)\n\}', txt, re.S).group(1)
    left = re.search(r'left\s*\{([^}]*)\}', ini, re.S).group(1)
    right = re.search(r'right\s*\{([^}]*)\}', ini, re.S).group(1)

    c = dict(
        lo=lo, hi=hi, cells=cells,
        dx=[(hi[i] - lo[i]) / cells[i] for i in range(3)],
        vol=[hi[i] - lo[i] for i in range(3)],
        gamma=float(re.search(r'gamma\s+([0-9.eE+-]+);', txt).group(1)),
        diaphragm=float(re.search(r'diaphragm\s+([0-9.eE+-]+);', ini).group(1)),
        rhoL=float(re.search(r'rho\s+([0-9.eE+-]+);', left).group(1)),
        uL=float(re.search(r'u\s+([0-9.eE+-]+);', left).group(1)),
        pL=float(re.search(r'p\s+([0-9.eE+-]+);', left).group(1)),
        rhoR=float(re.search(r'rho\s+([0-9.eE+-]+);', right).group(1)),
        uR=float(re.search(r'u\s+([0-9.eE+-]+);', right).group(1)),
        pR=float(re.search(r'p\s+([0-9.eE+-]+);', right).group(1)),
        dt=float(re.search(r'dt\s+([0-9.eE+-]+);', td).group(1)),
        end_time=float(re.search(r'endTime\s+([0-9.eE+-]+);', td).group(1)),
        vtu=re.search(r'vtk\s+([A-Za-z0-9_.-]+);', txt).group(1),
        text=txt,
    )
    c['steps'] = int(round(c['end_time'] / c['dt']))
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
    """Every report line as a dict, plus the start line's totals."""
    keys, rows = None, []
    for m in re.finditer(r'^\[chflow\] euler t (.*)$', out, re.M):
        v = m.group(1).split()
        if keys is None:
            keys = ['t'] + v[1::2]
        rows.append(dict(zip(keys, [float(x) for x in v[0::2]])))
    start = re.search(r'\[chflow\] euler start: cells (\d+) mass ([0-9.eE+-]+) '
                      r'mom ([0-9.eE+-]+) energy ([0-9.eE+-]+) dt ([0-9.eE+-]+) '
                      r'endTime ([0-9.eE+-]+) gamma ([0-9.eE+-]+)', out)
    if not rows or start is None:
        raise SystemExit("the log does not hold a start line and the report lines "
                         "the checker needs")
    return dict(rows=rows, cells=int(start.group(1)), mass0=float(start.group(2)),
                mom0=float(start.group(3)), energy0=float(start.group(4)),
                dt=float(start.group(5)), end_time=float(start.group(6)),
                gamma=float(start.group(7)))


# ------------------------------------------------------- exact Riemann solver
# Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics", ch. 4: the
# star pressure by Newton on f(p) = fL(p) + fR(p) + uR - uL, then the sampled
# solution. The checker needs it only at one time, along x.

def _f_k(p, rho, u, p0, gamma):
    """Toro (4.6)/(4.7): one side's pressure function."""
    if p > p0:
        a = 2.0 / ((gamma + 1.0) * rho)
        b = (gamma - 1.0) / (gamma + 1.0) * p0
        return (p - p0) * math.sqrt(a / (p + b))
    c0 = math.sqrt(gamma * p0 / rho)
    return (2.0 * c0 / (gamma - 1.0)) * ((p / p0) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)


def _df_k(p, rho, u, p0, gamma):
    if p > p0:
        a = 2.0 / ((gamma + 1.0) * rho)
        b = (gamma - 1.0) / (gamma + 1.0) * p0
        s = math.sqrt(a / (p + b))
        return s + (p - p0) * (-0.5 * a / (p + b) ** 1.5)
    c0 = math.sqrt(gamma * p0 / rho)
    return (1.0 / (rho * c0)) * (p / p0) ** (-(gamma + 1.0) / (2.0 * gamma))


def riemann_star(gamma, rhoL, uL, pL, rhoR, uR, pR):
    """The star pressure and velocity; Newton from Toro's initial guess, with a
    bisection fallback on the (monotone) pressure function."""
    def f(p):
        return _f_k(p, rhoL, uL, pL, gamma) + _f_k(p, rhoR, uR, pR, gamma) + (uR - uL)

    p = 0.5 * (pL + pR)
    for _ in range(200):
        fp = f(p)
        if abs(fp) < 1e-14:
            break
        d = _df_k(p, rhoL, uL, pL, gamma) + _df_k(p, rhoR, uR, pR, gamma)
        step = fp / d
        if not (1e-10 < p - step < 1e4):
            break
        p -= step
    lo, hi = 1e-10, 1e4
    for _ in range(300):  # bisection polish / fallback: f is increasing in p
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    u = 0.5 * (uL + uR) + 0.5 * (_f_k(p, rhoR, uR, pR, gamma)
                                 - _f_k(p, rhoL, uL, pL, gamma))
    return p, u


def star_density(p_star, p0, rho0, gamma):
    """The density the star state carries on each side of the contact."""
    if p_star > p0:
        num = p_star / p0 + (gamma - 1.0) / (gamma + 1.0)
        den = (gamma - 1.0) / (gamma + 1.0) * p_star / p0 + 1.0
        return rho0 * num / den
    return rho0 * (p_star / p0) ** (1.0 / gamma)


def sample_wave_speeds(gamma, rhoL, uL, pL, rhoR, uR, pR, p_star):
    """The wave speeds the exact solution fronts move at, head to tail."""
    cL = math.sqrt(gamma * pL / rhoL)
    cR = math.sqrt(gamma * pR / rhoR)
    cl_star = cL * (p_star / pL) ** ((gamma - 1.0) / (2.0 * gamma))
    cr_star = cR * (p_star / pR) ** ((gamma - 1.0) / (2.0 * gamma))
    return dict(
        fan_head=uL - cL,                       # rarefaction head (leftward)
        fan_tail=None,                          # set by the caller with u_star
        contact=None,
        shock=uR + cR * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p_star / pR
                                  + (gamma - 1.0) / (2.0 * gamma)),
        cl_star=cl_star, cr_star=cr_star,
    )


# ------------------------------------------------- the discrete HLLC mirror
# The kernel's own flux, per unit area through a face whose normal points from
# the left (owner) to the right (neighbour) state. The checker re-runs the full
# recursion with this, so the formulas have to be the kernel's, not the book's.

def two_rarefaction_pressure(pL, uL, cL, pR, uR, cR, gamma):
    """The star pressure the two-rarefaction Riemann solver reports - the
    estimate the kernel's wave speeds are built from, mirrored here so this
    stays a check of the discrete scheme and not of a different one."""
    k = (gamma - 1.0) / (2.0 * gamma)
    return ((0.5 * (gamma - 1.0) * (uL - uR) + cL + cR)
            / (cL * pL ** -k + cR * pR ** -k)) ** (1.0 / k)


def hllc_1d(rhoL, uL, pL, rhoE_L, rhoR, uR, pR, rhoE_R, gamma):
    """Per-unit-area HLLC flux (mass, momentum, energy) with n = +x."""
    cL = math.sqrt(gamma * pL / rhoL)
    cR = math.sqrt(gamma * pR / rhoR)
    pStar = two_rarefaction_pressure(pL, uL, cL, pR, uR, cR, gamma)

    def q(pK):
        if pStar <= pK:
            return 1.0
        return math.sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pK)
                         + (gamma - 1.0) / (2.0 * gamma))

    sL = uL - cL * q(pL)
    sR = uR + cR * q(pR)
    sM = (pR - pL + rhoL * uL * (sL - uL) - rhoR * uR * (sR - uR)) \
        / (rhoL * (sL - uL) - rhoR * (sR - uR))
    fL = (rhoL * uL, rhoL * uL * uL + pL, (rhoE_L + pL) * uL)
    fR = (rhoR * uR, rhoR * uR * uR + pR, (rhoE_R + pR) * uR)
    if sL >= 0.0:
        return fL
    if sR <= 0.0:
        return fR
    rhoLs = rhoL * (sL - uL) / (sL - sM)
    rhoRs = rhoR * (sR - uR) / (sR - sM)
    eLs = rhoLs * (rhoE_L / rhoL + (sM - uL) * (sM + pL / (rhoL * (sL - uL))))
    eRs = rhoRs * (rhoE_R / rhoR + (sM - uR) * (sM + pR / (rhoR * (sR - uR))))
    if sM >= 0.0:
        us, fs, ss, uo = (rhoLs, rhoLs * sM, eLs), fL, sL, (rhoL, rhoL * uL, rhoE_L)
    else:
        us, fs, ss, uo = (rhoRs, rhoRs * sM, eRs), fR, sR, (rhoR, rhoR * uR, rhoE_R)
    return (fs[0] + ss * (us[0] - uo[0]),
            fs[1] + ss * (us[1] - uo[1]),
            fs[2] + ss * (us[2] - uo[2]))


def replicate(c):
    """The kernel's own stepping, cell by cell: conserved totals, primitives,
    per-face HLLC, wall flux p*Sf, the same scatter and update order."""
    n = c['cells'][0]
    dx, dy, dz = c['dx']
    a = dy * dz
    vol = dx * a
    gamma, dt = c['gamma'], c['dt']

    mass = [0.0] * n
    mom = [0.0] * n
    energy = [0.0] * n
    for i in range(n):
        xc = c['lo'][0] + (i + 0.5) * dx
        left = xc < c['diaphragm']
        rho = c['rhoL'] if left else c['rhoR']
        p = c['pL'] if left else c['pR']
        u = c['uL'] if left else c['uR']
        mass[i] = rho * vol
        mom[i] = rho * u * vol
        energy[i] = rho * (p / ((gamma - 1.0) * rho) + 0.5 * u * u) * vol

    impulse = 0.0
    t = 0.0
    rho, vel, pre = [0.0] * n, [0.0] * n, [0.0] * n
    for _ in range(c['steps']):
        for i in range(n):
            rho[i] = mass[i] / vol
            vel[i] = mom[i] / mass[i]
            e = energy[i] / mass[i] - 0.5 * vel[i] * vel[i]
            pre[i] = (gamma - 1.0) * rho[i] * e
        cfl = max(abs(vel[i]) + math.sqrt(gamma * pre[i] / rho[i])
                  for i in range(n)) * dt / dx
        assert cfl <= 1.0 + 1e-12, "the replicated recursion broke its own CFL"

        dm = [0.0] * n
        dmom = [0.0] * n
        dE = [0.0] * n
        for f in range(n - 1):
            phi = hllc_1d(rho[f], vel[f], pre[f], energy[f] / vol,
                          rho[f + 1], vel[f + 1], pre[f + 1], energy[f + 1] / vol,
                          gamma)
            dm[f] -= phi[0] * a
            dmom[f] -= phi[1] * a
            dE[f] -= phi[2] * a
            dm[f + 1] += phi[0] * a
            dmom[f + 1] += phi[1] * a
            dE[f + 1] += phi[2] * a
        # walls: the xmin face (owner 0, Sf = -A x) and the xmax face
        # (owner n-1, Sf = +A x); the transverse walls contribute nothing
        dmom[0] -= pre[0] * (-a)
        impulse += dt * (pre[0] * (-a))
        dmom[n - 1] -= pre[n - 1] * a
        impulse += dt * (pre[n - 1] * a)
        for i in range(n):
            # the totals are extensive: dU/dt = -div Phi, so dt * flux
            mass[i] += dt * dm[i]
            mom[i] += dt * dmom[i]
            energy[i] += dt * dE[i]
        t += dt

    for i in range(n):
        rho[i] = mass[i] / vol
        vel[i] = mom[i] / mass[i]
        e = energy[i] / mass[i] - 0.5 * vel[i] * vel[i]
        pre[i] = (gamma - 1.0) * rho[i] * e
    return dict(x=[c['lo'][0] + (i + 0.5) * dx for i in range(n)],
                rho=rho, u=vel, p=pre, impulse=impulse, t=t)


def failures_exit(case, failures, notes):
    for note in notes:
        print("  " + note)
    if failures:
        raise SystemExit("Euler check FAILED:\n  " + "\n  ".join(failures))
    print("Euler %s check passed" % case)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    solver, case_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    c = case_constants(case_path)
    work = tempfile.mkdtemp(prefix='chflow_euler_')
    try:
        shutil.copy(case_path, os.path.join(work, 'case.dict'))
        proc = subprocess.run([solver, 'case.dict'], cwd=work,
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit("chflow-3d-compressible exited %d\n%s"
                             % (proc.returncode, proc.stderr))
        out = proc.stdout + proc.stderr
        log = log_numbers(out)
        cent, f = read_vtu(os.path.join(work, c['vtu']))
        failures, notes = [], []

        gamma, dx, n = c['gamma'], c['dx'][0], c['cells'][0]
        area = c['dx'][1] * c['dx'][2]
        vol = dx * area

        # ---- the mesh the case describes is the mesh the result carries ----
        if len(cent) != n or log['cells'] != n:
            failures.append("the result holds %d cells for the block's %d"
                            % (len(cent), n))
            failures_exit("sod", failures, notes)
        v0 = vol
        mass_ic = sum((c['rhoL'] if cent[i][0] < c['diaphragm'] else c['rhoR']) * v0
                      for i in range(n))
        energy_ic = sum((c['rhoL'] if cent[i][0] < c['diaphragm'] else c['rhoR'])
                        * (c['pL'] if cent[i][0] < c['diaphragm'] else c['pR'])
                        / ((gamma - 1.0)
                           * (c['rhoL'] if cent[i][0] < c['diaphragm'] else c['rhoR']))
                        * v0 for i in range(n))

        # ---- the exact star state, anchored to Toro's published Sod values ----
        ps, us = riemann_star(gamma, c['rhoL'], c['uL'], c['pL'],
                              c['rhoR'], c['uR'], c['pR'])
        rhoLs = star_density(ps, c['pL'], c['rhoL'], gamma)
        rhoRs = star_density(ps, c['pR'], c['rhoR'], gamma)
        notes.append("exact star state at the diaphragm: p* = %.5f, u* = %.5f, "
                     "rho*L = %.5f, rho*R = %.5f (Toro's Sod values: 0.30313, "
                     "0.92745, 0.42632, 0.26557)" % (ps, us, rhoLs, rhoRs))
        for name, got, want in (('p*', ps, 0.30313), ('u*', us, 0.92745),
                                ('rho*L', rhoLs, 0.42632), ('rho*R', rhoRs, 0.26557)):
            if abs(got - want) > 5e-4:
                failures.append("the checker's exact solver reads %s = %.5f, not "
                                "Toro's %.5f: it cannot anchor the plateau checks"
                                % (name, got, want))

        # ---- 1: the discrete solution is THE discrete solution ----
        rep = replicate(c)
        worst = dict(rho=0.0, u=0.0, p=0.0)
        for i in range(n):
            worst['rho'] = max(worst['rho'], abs(f['rho'][i] - rep['rho'][i]))
            worst['u'] = max(worst['u'], abs(f['U'][i][0] - rep['u'][i]))
            worst['p'] = max(worst['p'], abs(f['p'][i] - rep['p'][i]))
        notes.append("replicated recursion vs the written field: max |drho| %.3g, "
                     "|du| %.3g, |dp| %.3g over %d cells (the contract is 1e-9)"
                     % (worst['rho'], worst['u'], worst['p'], n))
        for k, v in worst.items():
            if v > 1e-9:
                failures.append("the written %s is not the replicated recursion's "
                                "(max error %.3g): the kernel is not running the "
                                "scheme its checker mirrors" % (k, v))
        trans = max(max(abs(f['U'][i][1]), abs(f['U'][i][2])) for i in range(n))
        notes.append("transverse velocity in the written field: %.3g (the 1D setup "
                     "cannot generate any)" % trans)
        if trans > 1e-12:
            failures.append("the transverse velocity is %.3g, not zero: the wall "
                            "fluxes leak momentum across the section" % trans)

        # ---- 2: conservation in the form the walls leave it ----
        m0, e0, p0 = log['mass0'], log['energy0'], log['mom0']
        if abs(m0 / mass_ic - 1.0) > 1e-12 or abs(e0 / energy_ic - 1.0) > 1e-12:
            failures.append("the reported totals (%.9g, %.9g) are not the initial "
                            "condition's (%.9g, %.9g)" % (m0, e0, mass_ic, energy_ic))
        for r in log['rows']:
            for k, v in r.items():
                if k != 't' and not math.isfinite(v):
                    failures.append("the report at t=%g has a non-finite %s"
                                    % (r['t'], k))
            if abs(r['mass'] / m0 - 1.0) > 1e-12:
                failures.append("mass drifted %.3g relative by t=%g (closed walls)"
                                % (abs(r['mass'] / m0 - 1.0), r['t']))
            if abs(r['energy'] / e0 - 1.0) > 1e-12:
                failures.append("total energy drifted %.3g relative by t=%g "
                                "(closed walls, no flux through them)"
                                % (abs(r['energy'] / e0 - 1.0), r['t']))
            resid = abs(r['mom'] + r['impulse'] - p0)
            if resid > 1e-12:
                failures.append("the momentum balance mom + impulse = mom(0) is "
                                "off by %.3g at t=%g" % (resid, r['t']))
            if r['pmin'] <= 0.0 or r['pmax'] <= 0.0:
                failures.append("pressure lost positivity at t=%g: [%.6g, %.6g]"
                                % (r['t'], r['pmin'], r['pmax']))
        notes.append("conservation over %d reports: mass and energy exact, "
                     "mom + impulse - mom(0) <= %.3g, p in [%.4g, %.4g]"
                     % (len(log['rows']),
                        max(abs(r['mom'] + r['impulse'] - p0) for r in log['rows']),
                        min(r['pmin'] for r in log['rows']),
                        max(r['pmax'] for r in log['rows'])))

        # while no front has touched a wall, both wall pressures are the
        # initial ones and the impulse has a closed form, per step and in total
        cL = math.sqrt(gamma * c['pL'] / c['rhoL'])
        t_fan_to_wall = (c['diaphragm'] - c['lo'][0]) / abs(c['uL'] - cL) \
            if c['uL'] < cL else float('inf')
        ws = sample_wave_speeds(gamma, c['rhoL'], c['uL'], c['pL'],
                                c['rhoR'], c['uR'], c['pR'], ps)
        t_shock_to_wall = (c['vol'][0] - c['diaphragm']) / ws['shock']
        last = log['rows'][-1]
        notes.append("fronts reach the walls at t = %.3f (fan head, x=0) and "
                     "t = %.3f (shock, x=1): the run ends at t = %g inside both"
                     % (t_fan_to_wall, t_shock_to_wall, last['t']))
        if last['t'] < min(t_fan_to_wall, t_shock_to_wall):
            want = c['dt'] * (c['pR'] - c['pL']) * area * c['steps']
            if abs(last['impulse'] / want - 1.0) > 1e-9:
                failures.append("the accumulated impulse %.6g is not the closed form "
                                "dt*(pR - pL)*A*steps = %.6g" % (last['impulse'], want))
            if abs(last['mom'] / (-want) - 1.0) > 1e-9:
                failures.append("the momentum %.6g is not growing at the closed-form "
                                "rate (pL - pR)*A per unit time (%.6g expected)"
                                % (last['mom'], -want))
            notes.append("impulse %.6g against the closed form %.6g; momentum "
                         "%.6g against (pL - pR)*A*t = %.6g"
                         % (last['impulse'], want, last['mom'], -want))

        # ---- 3: the exact Riemann solution, away from the fronts' own width ----
        t = last['t']
        x_head = c['diaphragm'] + (c['uL'] - cL) * t
        x_tail = c['diaphragm'] + (us - ws['cl_star']) * t
        x_contact = c['diaphragm'] + us * t
        x_shock = c['diaphragm'] + ws['shock'] * t
        xs = [cent[i][0] for i in range(n)]
        uf = [f['U'][i][0] for i in range(n)]
        rho_, p_ = f['rho'], f['p']

        # sigma: how far a front of speed lambda has diffused by t. The scheme's
        # numerical viscosity is |lambda| dx / 2, so the profile is a few sigma
        # wide; exempting that physical width, rather than a cell count, is what
        # lets one tolerance serve any mesh.
        def sigma(lam):
            return math.sqrt(abs(lam) * dx / 2.0 * t)

        s_head, s_tail, s_contact, s_shock = (sigma(c['uL'] - cL),
                                              sigma(us - ws['cl_star']),
                                              sigma(us), sigma(ws['shock']))

        def exact_at(xi):
            """The exact sample at x = xi: inside the fan this is the
            self-similar profile, not a plateau constant, so the fan is compared
            pointwise instead of being exempted."""
            if xi <= x_head:
                return (c['rhoL'], c['uL'], c['pL'])
            if xi <= x_tail:
                s = (xi - c['diaphragm']) / t
                fan_u = 2.0 / (gamma + 1.0) * (cL + (gamma - 1.0) / 2.0 * c['uL'] + s)
                fan_c = 2.0 / (gamma + 1.0) * (cL - (gamma - 1.0) / 2.0 * (s - c['uL']))
                r = c['rhoL'] * (fan_c / cL) ** (2.0 / (gamma - 1.0))
                return (r, fan_u, c['pL'] * (r / c['rhoL']) ** gamma)
            if xi <= x_contact:
                return (rhoLs, us, ps)
            if xi <= x_shock:
                return (rhoRs, us, ps)
            return (c['rhoR'], c['uR'], c['pR'])

        regions = dict(
            far_left=(0.0, x_head - 3.0 * s_head),
            fan=(x_head - 3.0 * s_head, x_tail + 2.0 * s_tail),
            star_left=(x_tail + 2.0 * s_tail, x_contact - 4.0 * s_contact),
            star_right=(x_contact + 4.0 * s_contact, x_shock - 2.0 * s_shock),
            far_right=(x_shock + 2.0 * s_shock, c['vol'][0]),
        )

        def window(lo, hi):
            return [i for i in range(n) if lo <= xs[i] <= hi]

        notes.append("fronts at t=%g: fan head %.4f, fan tail %.4f, contact %.4f, "
                     "shock %.4f; their diffusive widths sigma %.4f, %.4f, %.4f, %.4f"
                     % (t, x_head, x_tail, x_contact, x_shock,
                        s_head, s_tail, s_contact, s_shock))

        # The star plateaus' density is the tightest of these: the contact's
        # profile is 4 sigma wide, the shock's 2, and the plateau between them is
        # only (x_shock - x_contact) wide -- but its pressure and velocity are
        # uniform across the whole run and read to well under 1%.
        TOL = dict(far_left=(0.01, 0.01, 0.01 * abs(us)),
                   fan=(0.08, 0.10, 0.10 * abs(us)),
                   star_left=(0.03, 0.03, 0.02 * abs(us)),
                   star_right=(0.03, 0.03, 0.02 * abs(us)),
                   far_right=(0.03, 0.03, 0.02 * abs(us)))
        for region in ('far_left', 'fan', 'star_left', 'star_right', 'far_right'):
            idx = window(*regions[region])
            if not idx:
                failures.append("the %s window [%g, %g] is empty on %d cells: the "
                                "fronts are not where the exact solution puts them"
                                % (region, *regions[region], n))
                continue
            tr, tp, tu = TOL[region]
            er = max(abs(rho_[i] - exact_at(xs[i])[0]) / exact_at(xs[i])[0] for i in idx)
            ep = max(abs(p_[i] - exact_at(xs[i])[2]) / exact_at(xs[i])[2] for i in idx)
            eu = max(abs(uf[i] - exact_at(xs[i])[1]) for i in idx)
            notes.append("%s over %d cells: max |drho| %.2f%% (allowed %.0f%%), "
                         "|dp| %.2f%% (allowed %.0f%%), |du| %.5f (allowed %.5f)"
                         % (region.replace('_', ' '), len(idx), 100 * er, 100 * tr,
                            100 * ep, 100 * tp, eu, tu))
            if er > tr:
                failures.append("%s: density %.2f%% off the exact sample "
                                "(allowed %.0f%%)" % (region, 100 * er, 100 * tr))
            if ep > tp:
                failures.append("%s: pressure %.2f%% off the exact sample "
                                "(allowed %.0f%%)" % (region, 100 * ep, 100 * tp))
            if eu > tu:
                failures.append("%s: velocity %.5f off the exact sample "
                                "(allowed %.5f)" % (region, eu, tu))

        def crossing(level, lo, hi):
            """Where the density profile passes level, linearly between the two
            cell centres that bracket it: the front's own position, which the
            smeared profile still carries, unlike the steepest single cell."""
            for i in range(n - 1):
                if lo <= xs[i] and xs[i + 1] <= hi:
                    a, b = rho_[i], rho_[i + 1]
                    if (a - level) * (b - level) <= 0.0 and a != b:
                        return xs[i] + (level - a) / (b - a) * dx
            return None

        for name, level, lo, hi, want, s in (
                ('contact', 0.5 * (rhoLs + rhoRs), x_tail, x_shock, x_contact, s_contact),
                ('shock', 0.5 * (rhoRs + c['rhoR']), x_contact, c['vol'][0],
                 x_shock, s_shock)):
            found = crossing(level, lo, hi)
            if found is None:
                failures.append("the density profile never passes the %s's "
                                "mid-jump level %.5f between %g and %g: the front "
                                "was not captured at all" % (name, level, lo, hi))
                continue
            off = abs(found - want)
            notes.append("%s: the profile crosses its mid-jump density at %.4f, "
                         "%.3f sigma (%.2f dx) from the exact position %.4f"
                         % (name, found, off / s, off / dx, want))
            if off > 0.5 * s:
                failures.append("the %s sits %.4f from the exact position %.4f, "
                                "over half its own diffusive width %.4f"
                                % (name, off, want, s))

        # the velocity plateau's max is u* itself: cheap and sharp
        umax = max(uf)
        notes.append("max u in the field %.5f against u* = %.5f (%.2f%%)"
                     % (umax, us, 100.0 * (umax / us - 1.0)))
        if abs(umax / us - 1.0) > 0.02:
            failures.append("the field's max velocity %.5f is not the star "
                            "velocity %.5f" % (umax, us))

        # ---- 4: the whole-field error and its rate with the mesh ----
        def l1(centres, rho_f, u_f, p_f, h):
            out = [0.0, 0.0, 0.0]
            for i, xi in enumerate(centres):
                ex = exact_at(xi)
                out[0] += abs(rho_f[i] - ex[0])
                out[1] += abs(u_f[i] - ex[1])
                out[2] += abs(p_f[i] - ex[2])
            return [h * v for v in out]

        coarse = l1(xs, rho_, uf, p_, dx)
        fine_txt = re.sub(r'cells\s*\(([^)]*)\)',
                          lambda m: 'cells (%d %s)'
                          % (2 * c['cells'][0], m.group(1).split()[1] + ' '
                             + m.group(1).split()[2]), c['text'])
        fine_txt = re.sub(r'dt\s+([0-9.eE+-]+);',
                          'dt %g;' % (c['dt'] / 2.0), fine_txt, count=1)
        fine_txt = re.sub(r'vtk\s+[A-Za-z0-9_.-]+;', 'vtk sodFine.vtu;',
                          fine_txt, count=1)
        with open(os.path.join(work, 'caseFine.dict'), 'w') as fh:
            fh.write(fine_txt)
        fine_proc = subprocess.run([solver, 'caseFine.dict'], cwd=work,
                                   capture_output=True, text=True)
        if fine_proc.returncode != 0:
            raise SystemExit("the refined run exited %d\n%s"
                             % (fine_proc.returncode, fine_proc.stderr))
        fcent, ff = read_vtu(os.path.join(work, 'sodFine.vtu'))
        fxs = [cc[0] for cc in fcent]
        fine = l1(fxs, ff['rho'], [v[0] for v in ff['U']], ff['p'], dx / 2.0)
        for k, name in enumerate(('rho', 'u', 'p')):
            order = math.log(coarse[k] / fine[k]) / math.log(2.0)
            notes.append("L1 %s: %.6f on %d cells, %.6f on %d (order %.2f)"
                         % (name, coarse[k], n, fine[k], 2 * n, order))
            if not 0.4 <= order <= 1.5:
                failures.append("the L1 %s error refines at order %.2f, outside "
                                "[0.4, 1.5]: the scheme is not converging to the "
                                "exact solution" % (name, order))

        failures_exit('sod', failures, notes)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
