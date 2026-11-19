#!/usr/bin/env python3
"""Standalone verification of the MUSCL face-reconstruction operator.

The checker never trusts the app's numbers. It reads the case dict, rebuilds
the analytic field and the limiter in Python, and for every printed face it
recomputes the cell values from the printed cell centres (the only geometry it
accepts from the app), re-derives the limited slope and the face value, and
compares against what the app printed. Four profile modes are driven by the
case itself:

  linear    face values must be exact to 1e-12 (constant gradient, no limiting)
  quadratic independent re-derivation must match; no claimed exactness
  step      every face value must stay in [min(U,D), max(U,D)] (TVD), while the
            *unlimited* upwind extrapolation must overshoot somewhere
  sine      re-derivation on a coarse mesh, then L1 refinement order >= 1.7
            for muscl and in [0.8, 1.3] for upwind; a constant field must be
            reproduced exactly

usage: check_reconstruction.py <exe> <case.dict> [fine_case.dict]
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

TOL_MATCH = 1e-10       # app vs independent re-derivation
TOL_EXACT = 1e-12       # linear field / constant field exactness
TOL_BOUND = 1e-12       # TVD bound slack


# ------------------------------------------------------------------ case parsing
def strip_comments(txt):
    return "\n".join(re.sub(r'//.*', '', line) for line in txt.splitlines())


def triple(txt, key):
    m = re.search(r'%s\s*\(([^)]*)\)' % key, txt)
    if not m:
        raise SystemExit("case is missing '%s (...)'" % key)
    return [float(v) for v in m.group(1).split()]


def block(txt, key):
    m = re.search(r'\b%s\s*\{(.*?)\n\}' % key, txt, re.S)
    if not m:
        raise SystemExit("case is missing '%s {}'" % key)
    return m.group(1)


def word_of(txt, key, default=None):
    m = re.search(r'\b%s\s+([A-Za-z_][A-Za-z0-9_]*)\s*;' % key, txt)
    if m:
        return m.group(1)
    if default is not None:
        return default
    raise SystemExit("case is missing '%s <word>;'" % key)


def scalar_of(txt, key, default=None):
    m = re.search(r'\b%s\s+([0-9eE+.\-]+)\s*;' % key, txt)
    if m:
        return float(m.group(1))
    if default is not None:
        return default
    raise SystemExit("case is missing '%s <scalar>;'" % key)


def parse_case(path):
    txt = strip_comments(open(path).read())
    lo, hi = triple(txt, 'min'), triple(txt, 'max')
    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]

    vm = re.search(r'\bvelocity\s*\(([^)]*)\)', txt)
    if vm:
        vel = [float(v) for v in vm.group(1).split()]
        vel_kind = 'uniform'
    else:
        vb = block(txt, 'velocity')
        vel_kind = word_of(vb, 'type')
        if vel_kind == 'uniform':
            vel = triple(vb, 'U')
        elif vel_kind == 'rotation':
            vel = [word_of(vb, 'omega'), triple(vb, 'center')]
        else:
            raise SystemExit("unknown velocity type '%s'" % vel_kind)

    fb = block(txt, 'field')
    profile = word_of(fb, 'profile')
    field = dict(profile=profile, lo=lo, hi=hi)
    if profile == 'constant':
        field['value'] = scalar_of(fb, 'value')
    elif profile == 'linear':
        field['gradient'] = triple(fb, 'gradient')
        field['offset'] = scalar_of(fb, 'offset', 0.0)
    elif profile == 'quadratic':
        field['coefficient'] = scalar_of(fb, 'coefficient')
        field['offset'] = scalar_of(fb, 'offset', 0.0)
        field['center'] = triple(fb, 'center') if 'center' in fb else [0.0, 0.0, 0.0]
    elif profile == 'step':
        field['axis'] = 'xyz'.index(word_of(fb, 'axis', 'x'))
        field['threshold'] = scalar_of(fb, 'threshold')
        field['low'] = scalar_of(fb, 'low', 0.0)
        field['high'] = scalar_of(fb, 'high', 1.0)
    elif profile == 'sine':
        field['axis'] = 'xyz'.index(word_of(fb, 'axis', 'x'))
        field['amplitude'] = scalar_of(fb, 'amplitude', 1.0)
        field['waves'] = scalar_of(fb, 'waves')
        field['offset'] = scalar_of(fb, 'offset', 0.0)
    else:
        raise SystemExit("unknown field profile '%s'" % profile)

    sb_txt = txt if re.search(r'\bschemes\s*\{', txt) else ''
    mode = word_of(sb_txt, 'reconstruction', 'none') if sb_txt else 'none'
    limiter = word_of(sb_txt, 'limiter', 'minmod') if sb_txt else 'minmod'
    if mode not in ('none', 'muscl'):
        raise SystemExit("unknown reconstruction mode '%s'" % mode)
    if limiter not in ('minmod', 'vanLeer', 'MC', 'mc', 'superbee'):
        raise SystemExit("unknown limiter '%s'" % limiter)
    return dict(lo=lo, hi=hi, cells=cells, vel=vel, vel_kind=vel_kind,
                field=field, mode=mode, limiter=limiter, text=txt)


# ------------------------------------------------------------- analytic mirror
def analytic(f, x):
    p = f['profile']
    if p == 'constant':
        return f['value']
    if p == 'linear':
        return sum(g * xi for g, xi in zip(f['gradient'], x)) + f['offset']
    if p == 'quadratic':
        r = [x[i] - f['center'][i] for i in range(3)]
        return f['coefficient'] * sum(v * v for v in r) + f['offset']
    if p == 'step':
        return f['low'] if x[f['axis']] < f['threshold'] else f['high']
    if p == 'sine':
        a = f['axis']
        L = f['hi'][a] - f['lo'][a]
        return f['offset'] + f['amplitude'] * math.sin(
            2.0 * math.pi * f['waves'] * (x[a] - f['lo'][a]) / L)
    raise SystemExit("analytic: unknown profile " + p)


def limiter_phi(r, limiter):
    """Audit's ratio-form limiter phi(r). phi(1)==1 for every limiter."""
    if r <= 0.0:
        return 0.0
    if limiter == 'minmod':
        return min(r, 1.0)
    if limiter == 'vanLeer':
        return 2.0 * r / (1.0 + r)
    if limiter in ('MC', 'mc'):
        return min(2.0 * r, 0.5 * (1.0 + r), 2.0)
    if limiter == 'superbee':
        return max(min(2.0 * r, 1.0), min(r, 2.0))
    raise SystemExit("limiter_phi: unknown limiter " + limiter)


def limited_slope(s1, s2, limiter):
    """gP-like upwind slope s1, face slope s2; returns phi(s2/s1)*s1."""
    if s1 == 0.0:
        return 0.0
    return limiter_phi(s2 / s1, limiter) * s1


def parse_vec(s):
    return tuple(float(v) for v in s.split())


FACE_RE = re.compile(
    r'^\[chflow\] recon f=(\d+) fc=\(([^)]*)\) U=\(([^)]*)\) D=\(([^)]*)\) '
    r'F=\(([^)]*)\) ghost=([01]) exact=(\S+) recon=(\S+) err=(\S+)\s*$')


def run_app(exe, case_path, cwd):
    proc = subprocess.run([exe, case_path], cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit("app exited %d\n%s\n%s" % (proc.returncode, proc.stdout, proc.stderr))
    faces = []
    summary = None
    for line in (proc.stdout + proc.stderr).splitlines():
        m = FACE_RE.match(line.strip())
        if m:
            faces.append(dict(
                f=int(m.group(1)), fc=parse_vec(m.group(2)), U=parse_vec(m.group(3)),
                D=parse_vec(m.group(4)), F=parse_vec(m.group(5)),
                ghost=(m.group(6) == '1'), exact=float(m.group(7)),
                recon=float(m.group(8)), err=float(m.group(9))))
        if 'reconstruction summary:' in line:
            summary = dict(
                mode=re.search(r'\bmode=(\S+)', line).group(1),
                limiter=re.search(r'\blimiter=(\S+)', line).group(1),
                faces=int(re.search(r'faces=(\d+)', line).group(1)),
                L1=float(re.search(r'L1=(\S+)', line).group(1)),
                Linf=float(re.search(r'Linf=(\S+)', line).group(1)))
    if not faces or summary is None:
        raise SystemExit("no face lines / summary in the app output:\n" + proc.stdout + proc.stderr)
    return faces, summary


def recompute(faces, f, mode, limiter):
    """Independently rebuild every face value from the printed geometry.
    Returns (abs_errors, unbounded_overshoot, max_match_err)."""
    errs, match_err, over = [], 0.0, 0.0
    for e in faces:
        Tu = analytic(f, e['U'])
        Td = analytic(f, e['D'])
        if e['ghost']:
            Tf = 2.0 * Tu - Td
        else:
            Tf = analytic(f, e['F'])
        exact = analytic(f, e['fc'])

        du = math.dist(e['fc'], e['U'])
        dd = math.dist(e['D'], e['U'])
        r = du / dd

        # unlimited (backward-difference) extrapolation, the control probe
        unlimited = Tu + (Tu - Tf) * r
        lo, hi = min(Tu, Td), max(Tu, Td)
        if unlimited < lo - TOL_BOUND:
            over = max(over, lo - unlimited)
        if unlimited > hi + TOL_BOUND:
            over = max(over, unlimited - hi)

        if mode == 'muscl':
            L = limited_slope(Tu - Tf, Td - Tu, limiter)
            recon = Tu + L * r
            recon = min(max(recon, lo), hi)
        else:
            recon = Tu

        match_err = max(match_err, abs(recon - e['recon']))
        if abs(e['exact'] - exact) > TOL_MATCH:
            raise SystemExit("exact value at face %d disagrees: app %.17g, checker %.17g"
                             % (e['f'], e['exact'], exact))
        if abs(e['err'] - abs(exact - recon)) > TOL_MATCH:
            raise SystemExit("error at face %d disagrees: app %.17g, checker %.17g"
                             % (e['f'], e['err'], abs(exact - recon)))
        errs.append(abs(exact - recon))
    return errs, over, match_err


def l1(errs):
    return sum(errs) / len(errs) if errs else 0.0


def run_case(exe, case_path, work):
    os.makedirs(work, exist_ok=True)
    c = parse_case(case_path)
    faces, summary = run_app(exe, os.path.abspath(case_path), work)
    return c, faces, summary


def check_common(c, faces, summary, label):
    notes = []
    errs, over, match_err = recompute(faces, c['field'], c['mode'], c['limiter'])
    notes.append("%s: %d faces, mode=%s, limiter=%s; re-derivation max |dT| %.3g"
                 % (label, len(faces), c['mode'], c['limiter'], match_err))
    if match_err > TOL_MATCH:
        raise SystemExit("%s: app reconstruction disagrees with the independent "
                         "re-derivation by %.3g" % (label, match_err))
    expect_n = (c['cells'][0] - 1) * c['cells'][1] * c['cells'][2] \
        + c['cells'][0] * (c['cells'][1] - 1) * c['cells'][2] \
        + c['cells'][0] * c['cells'][1] * (c['cells'][2] - 1)
    if len(faces) != expect_n or summary['faces'] != expect_n:
        raise SystemExit("%s: expected %d internal faces, app printed %d / summary %d"
                         % (label, expect_n, len(faces), summary['faces']))
    if abs(l1(errs) - summary['L1']) > TOL_MATCH or \
       abs(max(errs) - summary['Linf']) > TOL_MATCH:
        raise SystemExit("%s: printed L1/Linf is not the recomputed mean/max error"
                         % label)
    return notes, errs, over


def rewrite_case(txt, work, name, subs):
    out = txt
    for old, new in subs:
        if old not in out:
            raise SystemExit("cannot rewrite '%s' in the case for %s" % (old, name))
        out = out.replace(old, new, 1)
    path = os.path.join(work, name)
    with open(path, 'w') as fh:
        fh.write(out)
    return path


def check_linear(exe, case, work):
    c, faces, summary = run_case(exe, case, work)
    notes, errs, _ = check_common(c, faces, summary, 'linear')
    worst = max(errs)
    notes.append("linear: max |exact - reconstructed| = %.3e (must be <= %.0e)"
                 % (worst, TOL_EXACT))
    if worst > TOL_EXACT:
        raise SystemExit("linear MUSCL is not exact: worst error %.3e" % worst)
    return notes


def check_quadratic(exe, case, work):
    c, faces, summary = run_case(exe, case, work)
    notes, errs, _ = check_common(c, faces, summary, 'quadratic')
    notes.append("quadratic: re-derivation matched; max face error %.3e" % max(errs))
    return notes


def check_step(exe, case, work):
    c, faces, summary = run_case(exe, case, work)
    notes, errs, over = check_common(c, faces, summary, 'step')
    # TVD: every reconstructed face value is within the two cells it joins.
    worst_violation = 0.0
    for e in faces:
        Tu = analytic(c['field'], e['U'])
        Td = analytic(c['field'], e['D'])
        lo, hi = min(Tu, Td), max(Tu, Td)
        if e['recon'] < lo - TOL_BOUND or e['recon'] > hi + TOL_BOUND:
            worst_violation = max(worst_violation,
                                  max(lo - e['recon'], e['recon'] - hi))
    notes.append("step: largest bound violation %.3e; unlimited probe overshoot %.4g"
                 % (worst_violation, over))
    if worst_violation > TOL_BOUND:
        raise SystemExit("step reconstruction overshot by %.3e (not TVD)" % worst_violation)
    if over <= 1e-3:
        raise SystemExit("the unlimited control probe did not overshoot (%.3g): "
                         "the step cannot demonstrate the limiter's value" % over)
    return notes


def check_sine(exe, case, fine_case, work):
    cc, cf, cs = run_case(exe, case, os.path.join(work, 'coarse'))
    coarse_notes, coarse, _ = check_common(cc, cf, cs, 'sine coarse')
    notes = coarse_notes
    notes.append("sine coarse: max face error %.3e" % max(coarse))

    fine_c, fine_f, fine_s = run_case(exe, fine_case, os.path.join(work, 'fine'))
    fine_notes, fine, _ = check_common(fine_c, fine_f, fine_s, 'sine fine')
    notes.extend(fine_notes)

    if len(fine) <= len(coarse):
        raise SystemExit("the fine sine mesh is not finer than the coarse one")
    order = math.log(l1(coarse) / l1(fine)) / math.log(2.0)
    notes.append("sine muscl: L1 %.6e -> %.6e, order %.3f (want >= 1.7)"
                 % (l1(coarse), l1(fine), order))
    if not order >= 1.7:
        raise SystemExit("muscl L1 order %.3f is below 1.7" % order)

    # Upwind (reconstruction none) on both meshes, generated from each case.
    coarse_txt = strip_comments(open(case).read())
    fine_txt = strip_comments(open(fine_case).read())
    up_coarse_path = rewrite_case(coarse_txt, work, 'upwindCoarse.dict',
                                  [('reconstruction muscl;', 'reconstruction none;')])
    up_fine_path = rewrite_case(fine_txt, work, 'upwindFine.dict',
                                [('reconstruction muscl;', 'reconstruction none;')])
    uc, uf, us = run_case(exe, up_coarse_path, os.path.join(work, 'upCoarse'))
    vfc, vff, vfs = run_case(exe, up_fine_path, os.path.join(work, 'upFine'))
    _, up_c_errs, _ = check_common(uc, uf, us, 'sine upwind coarse')
    _, up_f_errs, _ = check_common(vfc, vff, vfs, 'sine upwind fine')
    up_order = math.log(l1(up_c_errs) / l1(up_f_errs)) / math.log(2.0)
    notes.append("sine upwind: L1 %.6e -> %.6e, order %.3f (want [0.8,1.3])"
                 % (l1(up_c_errs), l1(up_f_errs), up_order))
    if not 0.8 <= up_order <= 1.3:
        raise SystemExit("upwind L1 order %.3f is outside [0.8,1.3]" % up_order)

    # Constant field: replace the sine block and assert exactness everywhere.
    const_txt = re.sub(r'profile\s+sine\s*;', 'profile constant; value 3.25;', coarse_txt)
    if const_txt == coarse_txt:
        raise SystemExit("could not rewrite the sine profile to constant")
    const_path = os.path.join(work, 'constant.dict')
    with open(const_path, 'w') as fh:
        fh.write(const_txt)
    cc, cf, cs = run_case(exe, const_path, os.path.join(work, 'constant'))
    _, c_errs, _ = check_common(cc, cf, cs, 'sine constant')
    notes.append("constant field: max face error %.3e (must be <= %.0e)"
                 % (max(c_errs), TOL_EXACT))
    if max(c_errs) > TOL_EXACT:
        raise SystemExit("constant field is not reproduced exactly: %.3e" % max(c_errs))
    return notes


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    exe = os.path.abspath(sys.argv[1])
    case = os.path.abspath(sys.argv[2])
    fine_case = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else \
        os.path.join(os.path.dirname(os.path.dirname(case)), 'reconstructionSineFine',
                     'case.dict')
    c = parse_case(case)
    profile = c['field']['profile']
    work = tempfile.mkdtemp(prefix='chflow_recon_')
    try:
        if profile == 'linear':
            notes = check_linear(exe, case, work)
        elif profile == 'quadratic':
            notes = check_quadratic(exe, case, work)
        elif profile == 'step':
            notes = check_step(exe, case, work)
        elif profile == 'sine':
            notes = check_sine(exe, case, fine_case, work)
        else:
            raise SystemExit("no checker mode for profile '%s'" % profile)
        for n in notes:
            print("  " + n)
        print("reconstruction %s check passed" % profile)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
