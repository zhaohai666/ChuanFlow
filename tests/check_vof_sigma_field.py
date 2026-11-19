#!/usr/bin/env python3
"""Verify 3D VOF per-cell Weinaug-Katz surface tension (D13/D10, section 5.120).

usage: check_vof_sigma_field.py <exe> <case.dict>

The case declares a `sigmaFromP{}` block: an analytic absolute thermodynamic
pressure p_cell on cell centres, a feed composition, SI parachors and T. The
solver flashes every cell at its own pressure and feeds the resulting kinematic
sigma into a cell-wise continuum surface force. This checker re-derives all of
that from scratch, in pure Python (stdlib only):

  flash    a verbatim Python transcription of the project's already-verified
           PR flash (src/pvt/src/{components,pengRobinson,flash}.cpp): Wilson
           K-values, Rachford-Rice, successive substitution on fugacity
           equality, and the two-phase/single-phase verdict.
  formula  the Weinaug-Katz parachor equation, the same transcription that
           tests/check_pvt_surface_tension.py pins to 1e-12 against the
           `chemicals` worked examples:
               sigma = [ sum_i P_i (x_i/Vml - y_i/Vmg) ]^4
  field    every written cell's sigma must equal WK(p_cell), and sigma really
           has to vary from cell to cell; sigmaOverRho = sigma / rho.
  force    csf[u] = gradAlpha[u] * (-sigmaOverRho[u] * kappa[u]) cell by cell,
           using the gradient the solver itself wrote to the .vtu.
  gate     a scalar twin (sigmaFromP removed, scalar sigma declared instead)
           has no per-cell fields in its .vtu and its force uses exactly that
           scalar; removing sigmaFromP without declaring sigma is rejected.
  reject   unknown component, missing table parachor, malformed composition,
           non-positive temperature / p0, wrong parachor count and an unknown
           pressure-field type all have to be loud errors.

No platform library, no numpy, no dependency on chflow-pvt being built.
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case_path = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])

# ---------------------------------------------------------------- PR + WK ---
# Constants, component data and algorithm transcribed from the verified PVT
# implementation; see src/pvt/src/{components,pengRobinson,flash}.cpp.
R = 8.314462618
SQRT2 = 1.4142135623730951
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846
CRIT_V_OVER_B = 0.30740 / 0.07779607390388846

# name -> (molar mass [g/mol], Tc [K], Pc [bar], acentric)
COMPONENTS = {
    "methane": (16.043, 190.564, 45.992, 0.0115),
    "ethane": (30.070, 305.322, 48.722, 0.0995),
    "propane": (44.097, 369.834, 42.477, 0.1523),
    "iButane": (58.123, 407.817, 36.480, 0.1835),
    "nButane": (58.123, 425.125, 37.962, 0.2002),
    "nPentane": (72.150, 469.700, 33.700, 0.2515),
    "nHexane": (86.178, 507.600, 30.250, 0.3013),
    "carbonDioxide": (44.010, 304.128, 73.773, 0.2239),
    "hydrogenSulphide": (34.081, 373.400, 89.628, 0.0900),
    "nitrogen": (28.014, 126.192, 33.958, 0.0372),
    "water": (18.01528, 647.096, 220.64, 0.3443),
}
# Weinaug & Katz (1943) cgs parachors (Table 18.1 as transcribed in
# src/pvt/src/componentData.cpp).
PARACHORS_CGS = {
    "carbonDioxide": 78.0, "nitrogen": 41.0, "methane": 77.0,
    "ethane": 108.0, "propane": 150.3, "iButane": 181.5,
    "nButane": 189.9, "iPentane": 225.0, "nPentane": 231.5,
    "nHexane": 271.0, "nHeptane": 312.5, "nOctane": 351.5,
}
CGS_TO_SI = 10.0 ** -6.75


def component(name):
    mw, tc, pc, w = COMPONENTS[name]
    return (mw * 1.0e-3, tc, pc * 1.0e5, w)


def wilson_k(comp, T, p):
    _, tc, pc, w = comp
    return pc / p * math.exp(5.373 * (1.0 + w) * (1.0 - tc / T))


def cubic_real_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p * p * p / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        s = math.sqrt(disc)
        return [math.cbrt(-q / 2.0 + s) + math.cbrt(-q / 2.0 - s) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p * p * p / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    roots = [m * math.cos((theta + 2.0 * k * math.pi) / 3.0) + shift for k in range(3)]
    return sorted(roots)


def mixture_parameters(comps, x, T, p):
    n = len(comps)
    a_i, b_i = [], []
    a_ij = [[0.0] * n for _ in range(n)]
    RT = R * T
    for i, comp in enumerate(comps):
        _, tc, pc, w = comp
        kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w
        alpha = (1.0 + kappa * (1.0 - math.sqrt(T / tc))) ** 2
        rtc = R * tc
        a_i.append(OMEGA_A * rtc * rtc / pc * alpha)
        b_i.append(OMEGA_B * rtc / pc)
    a_mix = b_mix = 0.0
    for i in range(n):
        b_mix += x[i] * b_i[i]
        for j in range(n):
            aij = math.sqrt(a_i[i] * a_i[j])  # kij = 0
            a_ij[i][j] = aij
            a_mix += x[i] * x[j] * aij
    A = a_mix * p / (RT * RT)
    B = b_mix * p / RT
    return A, B, a_mix, b_mix, a_i, b_i, a_ij


def compressibility_roots(A, B):
    raw = cubic_real_roots(-(1.0 - B), A - 3.0 * B * B - 2.0 * B,
                           -(A * B - B * B - B * B * B))
    return [z for z in raw if z > B * (1.0 + 1e-12)]


def peng_robinson_state(comps, x, T, p, pref):
    A, B, a_mix, b_mix, a_i, b_i, a_ij = mixture_parameters(comps, x, T, p)
    roots = compressibility_roots(A, B)
    if not roots:
        raise RuntimeError("PR cubic has no physical root")
    n = len(x)
    RT = R * T
    log_term_base = 2.0 * SQRT2 * B
    mix_m = sum(x[i] * comps[i][0] for i in range(n))

    def build(Z):
        vm = Z * RT / p
        log_term = math.log((Z + (1.0 + SQRT2) * B) / (Z + (1.0 - SQRT2) * B))
        ln_phi, gibbs = [], 0.0
        for i in range(n):
            sum_a = sum(x[j] * a_ij[i][j] for j in range(n))
            bi_b = b_i[i] / b_mix
            lp = (bi_b * (Z - 1.0) - math.log(Z - B)
                  - A / log_term_base * (2.0 * sum_a / a_mix - bi_b) * log_term)
            ln_phi.append(lp)
            gibbs += x[i] * lp
        return dict(Z=Z, Vm=vm, lnPhi=ln_phi, gibbs=gibbs, VoverB=vm / b_mix)

    distinct = []
    for z in roots:
        if z not in distinct:
            distinct.append(z)
    if pref == "Liquid":
        pick = min(distinct)
    elif pref == "Vapor":
        pick = max(distinct)
    else:
        pick = distinct[0]
        for z in distinct:
            if build(z)["gibbs"] < build(pick)["gibbs"]:
                pick = z
    return build(pick)


def normalise(v):
    s = sum(v)
    if s <= 0.0:
        raise RuntimeError("zero composition")
    return [q / s for q in v]


def rr_residue(z, K, beta):
    return sum(zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0))
               for zi, ki in zip(z, K))


def solve_rr(z, K):
    lo, hi, beta = 0.0, 1.0, 0.5
    for _ in range(100):
        f = rr_residue(z, K, beta)
        if f > 0.0:
            lo = beta
        else:
            hi = beta
        df = 0.0
        for zi, ki in zip(z, K):
            d = ki - 1.0
            w = 1.0 + beta * d
            df -= zi * d * d / (w * w)
        nxt = beta - f / df if df != 0.0 else 0.5 * (lo + hi)
        if not (lo < nxt < hi):
            nxt = 0.5 * (lo + hi)
        if abs(nxt - beta) < 1e-15:
            beta = nxt
            break
        beta = nxt
    return beta


def flash(comps, feed, T, p):
    n = len(comps)
    z = normalise(list(feed))
    K = [wilson_k(c, T, p) for c in comps]
    for _ in range(1000):
        f0, f1 = rr_residue(z, K, 0.0), rr_residue(z, K, 1.0)
        x, y = [0.0] * n, [0.0] * n
        if f0 <= 0.0:
            x = list(z)
            ys = [K[i] * z[i] for i in range(n)]
            s = sum(ys)
            if s <= 0.0:
                raise RuntimeError("zero vapour composition")
            y = [q / s for q in ys]
        elif f1 >= 0.0:
            y = list(z)
            xs = [z[i] / K[i] for i in range(n)]
            s = sum(xs)
            if s <= 0.0:
                raise RuntimeError("zero liquid composition")
            x = [q / s for q in xs]
        else:
            beta = solve_rr(z, K)
            x = [z[i] / (1.0 + beta * (K[i] - 1.0)) for i in range(n)]
            y = [K[i] * x[i] for i in range(n)]
            x, y = normalise(x), normalise(y)
        sl = peng_robinson_state(comps, x, T, p, "Liquid")
        sv = peng_robinson_state(comps, y, T, p, "Vapor")
        change = 0.0
        for i in range(n):
            new_k = math.exp(sl["lnPhi"][i] - sv["lnPhi"][i])
            change = max(change, abs(math.log(new_k / K[i])))
            K[i] = new_k
        if change < 1.0e-12:
            break

    f0, f1 = rr_residue(z, K, 0.0), rr_residue(z, K, 1.0)
    spread = max(abs(k - 1.0) for k in K)
    two_phase = False
    beta_final = 0.0
    if f0 > 0.0 and f1 < 0.0 and spread > 1e-6:
        beta = solve_rr(z, K)
        if 0.0 < beta < 1.0:
            xt = [z[i] / (1.0 + beta * (K[i] - 1.0)) for i in range(n)]
            yt = [K[i] * xt[i] for i in range(n)]
            xt, yt = normalise(xt), normalise(yt)
            if max(abs(yt[i] - xt[i]) for i in range(n)) > 1e-8:
                two_phase = True
                beta_final = beta

    if not two_phase:
        one = peng_robinson_state(comps, z, T, p, "LowestGibbs")
        return True, list(z), list(z), one["Vm"], one["Vm"]
    x = [z[i] / (1.0 + beta_final * (K[i] - 1.0)) for i in range(n)]
    y = [K[i] * x[i] for i in range(n)]
    x, y = normalise(x), normalise(y)
    vl = peng_robinson_state(comps, x, T, p, "Liquid")["Vm"]
    vg = peng_robinson_state(comps, y, T, p, "Vapor")["Vm"]
    return False, x, y, vl, vg


def weinaug_katz(Ps, Vml, Vmg, xs, ys):
    """The same formula check_pvt_surface_tension.py pins to 1e-12."""
    rho_l, rho_v = 1.0 / Vml, 1.0 / Vmg
    tot = sum(p * (x * rho_l - y * rho_v) for p, x, y in zip(Ps, xs, ys))
    return tot ** 4


def wk_sigma(T, p, names, z, parachors):
    comps = [component(n) for n in names]
    single, x, y, vl, vg = flash(comps, z, T, p)
    if single or not (vl > 0.0 and vg > 0.0):
        return 0.0
    return weinaug_katz(parachors, vl, vg, x, y)


# ------------------------------------------------------------ case / vtu ---
def read_file(path):
    with open(path) as fh:
        return fh.read()


def extract_block(text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", text)
    if not m:
        raise SystemExit("case has no %s{} block" % key)
    depth, i = 1, m.end()
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    if depth:
        raise SystemExit("unbalanced %s{} block" % key)
    return text[m.end():i - 1]


def remove_block(text, key):
    """Delete `key { ... }`, matching braces, so nested blocks go with it."""
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", text)
    if not m:
        raise SystemExit("case has no %s{} block" % key)
    depth, i = 1, m.end()
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    return text[:m.start()] + text[i:]


def num(pattern, text):
    m = re.search(pattern, text)
    if not m:
        raise SystemExit("cannot read number matching %r" % pattern)
    return float(m.group(1))


def parse_case(path):
    txt = read_file(path)
    sfp = extract_block(txt, "sigmaFromP")
    T = num(r"temperature\s+([0-9.eE+-]+)\s*;", sfp)
    comp = extract_block(sfp, "composition")
    pairs = re.findall(r"([A-Za-z][A-Za-z0-9]*)\s+([0-9.eE+-]+)\s*;", comp)
    names = [n for n, _ in pairs]
    z = [float(v) for _, v in pairs]
    pn = re.search(r"parachors\s+([^;]+);", sfp).group(1).strip()
    if pn == "fromTable":
        parachors = [PARACHORS_CGS[n] * CGS_TO_SI for n in names]
    else:
        parachors = [float(v) for v in re.findall(r"[0-9.eE+-]+", pn)]
    pf = extract_block(sfp, "pressureField")
    p0 = num(r"p0\s+([0-9.eE+-]+)\s*;", pf)
    gradient = [float(v) for v in re.search(r"gradient\s*\(([^)]*)\)", pf).group(1).split()]
    ref = [float(v) for v in re.search(r"reference\s*\(([^)]*)\)", pf).group(1).split()]
    st = extract_block(txt, "surfaceTension")
    rho = num(r"\brho\s+([0-9.eE+-]+)\s*;", st)
    vtu = re.search(r"\bvtk\s+([A-Za-z0-9_.-]+)\s*;", extract_block(txt, "output")).group(1)
    return dict(txt=txt, T=T, names=names, z=z, parachors=parachors,
                p0=p0, gradient=gradient, reference=ref, rho=rho, vtu=vtu,
                sfp=sfp)


def read_vtu(path):
    txt = read_file(path)
    arrays = {}
    for m in re.finditer(r"<DataArray([^>]*)>(.*?)</DataArray>", txt, re.S):
        name = re.search(r'Name="([^"]+)"', m.group(1))
        if name:
            arrays[name.group(1)] = [float(v) for v in m.group(2).split()]
    head = re.search(r'<DataArray type="Float64" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
                     txt, re.S)
    pts = [float(v) for v in head.group(1).split()]
    points = list(zip(pts[0::3], pts[1::3], pts[2::3]))
    conn = [int(v) for v in arrays["connectivity"]]
    off = [int(v) for v in arrays["offsets"]]
    centres, s = [], 0
    for o in off:
        cl = [points[v] for v in conn[s:o]]
        s = o
        centres.append(tuple(sum(q[i] for q in cl) / len(cl) for i in range(3)))
    nf = len(centres)
    fields = {}
    for name, v in arrays.items():
        if name in ("connectivity", "offsets", "types") or len(v) % nf:
            continue
        if len(v) == nf:
            fields[name] = list(v)
        else:
            k = len(v) // nf
            fields[name] = [tuple(v[i * k:(i + 1) * k]) for i in range(nf)]
    return centres, fields


def run(exe, work, text, name="case.dict"):
    path = os.path.join(work, name)
    with open(path, "w") as fh:
        fh.write(text)
    proc = subprocess.run([exe, name], cwd=work, capture_output=True, text=True)
    return proc, path


# ------------------------------------------------------------------ checks ---
fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(("%s %s" % (name, detail)).strip())


def close(name, got, want, rtol=1e-9, atol=1e-14):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, "got %.12g want %.12g" % (got, want))


def main():
    cfg = parse_case(case_path)
    work = tempfile.mkdtemp(prefix="chflow_vof_sigma_")
    try:
        proc, _ = run(exe, work, cfg["txt"])
        if proc.returncode != 0:
            raise SystemExit("chflow-vof exited %d\n%s\n%s" %
                             (proc.returncode, proc.stdout, proc.stderr))
        out = proc.stdout + proc.stderr
        centres, f = read_vtu(os.path.join(work, cfg["vtu"]))
        n = len(centres)
        needed = ("alpha", "kappa", "csf", "sigma", "sigmaOverRho", "pThermo", "gradAlpha")
        for key in needed:
            check("field %s written" % key, key in f, "missing")
        if any(k not in f for k in needed):
            raise SystemExit("the .vtu is missing per-cell sigma fields")

        sigma, sOverRho = f["sigma"], f["sigmaOverRho"]
        pThermo, kappa = f["pThermo"], f["kappa"]
        csf, galpha = f["csf"], f["gradAlpha"]
        check("one sigma per cell", len(sigma) == n, "%d vs %d" % (len(sigma), n))

        # p_cell on the cell centre, recomputed independently from the mesh.
        p_ref, sig_ref = [], []
        for c in centres:
            pc = cfg["p0"] + sum(g * (c[i] - cfg["reference"][i])
                                 for i, g in enumerate(cfg["gradient"]))
            p_ref.append(pc)
            sig_ref.append(wk_sigma(cfg["T"], pc, cfg["names"],
                                    cfg["z"], cfg["parachors"]))
        for i in range(n):
            close("pThermo cell %d" % i, pThermo[i], p_ref[i], 2e-10, 1e-6)
            close("sigma WK cell %d" % i, sigma[i], sig_ref[i], 1e-8, 1e-15)
            close("sigma/rho cell %d" % i, sOverRho[i],
                  sig_ref[i] / cfg["rho"], 1e-8, 1e-16)

        smin, smax = min(sigma), max(sigma)
        check("sigma varies between cells", smax - smin > 1e-6 * smax,
              "min %.6g max %.6g" % (smin, smax))
        check("sigma all two-phase and positive", smin > 0.0, "min %.6g" % smin)
        check("sigmaFromP reported", "vof sigmaFromP:" in out, "no log line")

        fmax = max(math.sqrt(sum(v * v for v in q)) for q in csf)
        atol_f = 1e-9 * fmax + 1e-14
        worst = 0.0
        for i in range(n):
            for k in range(3):
                want = galpha[i][k] * (-sOverRho[i] * kappa[i])
                worst = max(worst, abs(csf[i][k] - want))
                check("csf cell %d comp %d" % (i, k),
                      abs(csf[i][k] - want) <= atol_f,
                      "got %.9g want %.9g" % (csf[i][k], want))

        # ---- scalar twin: the untouched legacy path must still be the one
        # a case without sigmaFromP takes. Same mesh and alpha, so the field
        # case's own gradAlpha is the reference gradient for the twin's force.
        sigma_scalar = wk_sigma(cfg["T"], cfg["p0"], cfg["names"],
                                cfg["z"], cfg["parachors"])
        twin = remove_block(cfg["txt"], "sigmaFromP")
        twin = re.sub(r"(?m)^(\s*)nu\s+([0-9.eE+-]+)\s*;",
                      r"\1sigma %.17g;\n\1nu \2;" % sigma_scalar, twin, count=1)
        check("twin dropped sigmaFromP", re.search(r"sigmaFromP\s*\{", twin) is None,
              "still there")
        check("twin declares scalar sigma", re.search(r"\bsigma\s+[0-9.eE+-]+\s*;", twin))
        tp, _ = run(exe, work, twin, "twin.dict")
        if tp.returncode != 0:
            check("scalar twin runs", False, tp.stdout + tp.stderr)
        else:
            # the twin writes the same output filename, overwriting the field
            # case's file on disk; the field case was already read into memory.
            _, tf = read_vtu(os.path.join(work, cfg["vtu"]))
            for key in ("sigma", "sigmaOverRho", "pThermo", "gradAlpha"):
                check("twin has no %s field" % key, key not in tf, "leaked into legacy path")
            tsigma = sigma_scalar / cfg["rho"]
            reported = re.search(r"vof csf: sigma/rho ([0-9.eE+-]+)",
                                 tp.stdout + tp.stderr)
            check("twin reports scalar sigma/rho", reported is not None, "no csf line")
            if reported:
                close("twin csf report", float(reported.group(1)), tsigma, 1e-11, 1e-16)
            worst_twin = 0.0
            for i in range(n):
                for k in range(3):
                    want = galpha[i][k] * (-tsigma * tf["kappa"][i])
                    worst_twin = max(worst_twin, abs(tf["csf"][i][k] - want))
            check("twin force uses scalar sigma", worst_twin <= atol_f,
                  "worst %.6g" % worst_twin)
            diff = max(math.sqrt(sum((f["csf"][i][k] - tf["csf"][i][k]) ** 2
                                     for k in range(3))) for i in range(n))
            check("per-cell force differs from scalar twin", diff > 1e-8 * fmax,
                  "max diff %.6g of %.6g" % (diff, fmax))

        # ---- rejections. Each is a loud error at case-loading time.
        def reject(label, text, needle):
            p, _ = run(exe, work, text, "reject.dict")
            combined = p.stdout + p.stderr
            check("reject " + label, p.returncode != 0 and needle in combined,
                  "rc=%d msg=%s" % (p.returncode, combined.strip()[:160]))

        bad = cfg["txt"].replace("methane 0.40;", "bogus 0.40;", 1)
        reject("unknown component", bad, "component 'bogus'")
        bad = cfg["txt"].replace("parachors fromTable;",
                                 "parachors (1e-5 1e-5);", 1)
        reject("wrong parachor count", bad, "one SI value per component")
        bad = cfg["txt"].replace("temperature 300.0;", "temperature -1;", 1)
        reject("negative temperature", bad, "temperature must be positive")
        bad = cfg["txt"].replace("p0       3.5e6;", "p0       -1;", 1)
        reject("negative p0", bad, "p0 must be positive")
        bad = cfg["txt"].replace("methane 0.40;", "methane -0.40;", 1)
        reject("negative fraction", bad, "must not be negative")
        bad = cfg["txt"].replace("type linear;", "type quadratic;", 1)
        reject("unknown pressure field", bad, "type must be linear")
        bad = cfg["txt"].replace("parachors fromTable;",
                                 "parachors fromTable;", 1).replace("methane 0.40;",
                                 "water 0.40;", 1)
        reject("no table parachor", bad, "no Weinaug-Katz")
        bad = remove_block(cfg["txt"], "sigmaFromP")
        reject("no sigma at all", bad, "missing required key 'sigma'")

        if fails:
            print("\n".join("FAIL: " + m for m in fails[:40]))
            if len(fails) > 40:
                print("... and %d more" % (len(fails) - 40))
            raise SystemExit("%d vof-sigma-field check(s) failed" % len(fails))
        print("vof-sigma-field: %d cells, sigma %.6g..%.6g N/m, "
              "WK(p_cell) matched, CSF cell-wise matched (worst %.3g), "
              "scalar twin clean, %d rejections - all checks passed" %
              (n, smin, smax, worst, 8))
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
