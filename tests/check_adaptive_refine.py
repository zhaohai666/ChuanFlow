#!/usr/bin/env python3
"""Independent acceptance test for the F13 adaptive refinement operator.

The checker never trusts the C++ app. It parses the case dictionaries, rebuilds
the dyadic slab refinement in pure Python, evaluates the analytic scalar field
and its exact box integral, and compares every reported number. It asserts:

  * volume conservation: sum(cell volumes) == box volume to 1e-9
  * scalar integral conservation: sum(phi_c V_c) == exact analytic integral
  * the uniform-field twin is not refined at all
  * the step-field refinement is localised around the interface
  * minVolume > 0 and minVolume <= V <= maxVolume for every cell
  * invalid thresholds / depths / field declarations are rejected loudly

usage: check_adaptive_refine.py <exe> <case.dict> <caseUniform.dict> <caseReject.dict>
"""
import math
import os
import re
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
TOL = 1e-9

fails = []
notes = []
nchecks = 0


def check(cond, msg):
    global nchecks
    nchecks += 1
    if not cond:
        fails.append(msg)


def close(a, b, rel=1e-9, abs_=1e-12):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def strip_comments(text):
    out = []
    i = 0
    while i < len(text):
        if text[i:i + 2] == "//" or text[i] == "#":
            j = text.find("\n", i)
            if j < 0:
                break
            i = j
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def extract_block(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("case is missing '%s {}'" % key)
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced braces in '%s {}'" % key)


def scalar_of(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if m:
        return float(m.group(1))
    if default is not None:
        return default
    raise RuntimeError("missing scalar '%s'" % key)


def word_of(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", txt)
    if m:
        return m.group(1)
    if default is not None:
        return default
    raise RuntimeError("missing word '%s'" % key)


def triple_of(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", txt)
    if not m:
        raise RuntimeError("missing list '%s (...)'" % key)
    v = [float(x) for x in m.group(1).split()]
    if len(v) != 3:
        raise RuntimeError("'%s' needs 3 components" % key)
    return v


def parse_case(path):
    txt = strip_comments(open(path, encoding="utf-8").read())
    blk = extract_block(txt, "block")
    adapt = extract_block(txt, "adapt")
    lo = triple_of(blk, "min")
    hi = triple_of(blk, "max")
    cells = [int(v) for v in triple_of(blk, "cells")]
    if any(not (n >= 1) for n in cells):
        raise RuntimeError("block cells must be positive")

    f = dict(kind=word_of(adapt, "field"))
    if f["kind"] == "linear":
        f["gradient"] = triple_of(adapt, "gradient")
        f["offset"] = scalar_of(adapt, "offset", 0.0)
    elif f["kind"] == "step":
        f["axis"] = "xyz".index(word_of(adapt, "axis", "x"))
        f["face"] = scalar_of(adapt, "interface")
        f["low"] = scalar_of(adapt, "low", 0.0)
        f["high"] = scalar_of(adapt, "high", 1.0)
    else:
        raise RuntimeError("unknown field '%s'" % f["kind"])

    return dict(lo=lo, hi=hi, cells=cells, field=f,
                threshold=scalar_of(adapt, "gradientThreshold"),
                maxDepth=int(scalar_of(adapt, "maxDepth")),
                dump=(word_of(adapt, "dumpCells", "no") == "yes"), text=txt)


# --------------------------------------------------------------- field mirror
def field_at(f, x):
    if f["kind"] == "linear":
        return sum(g * xi for g, xi in zip(f["gradient"], x)) + f["offset"]
    return f["low"] if x[f["axis"]] < f["face"] else f["high"]


def volume_of(lo, hi):
    v = 1.0
    for a in range(3):
        v *= hi[a] - lo[a]
    return v


def box_integral(f, lo, hi):
    vol = volume_of(lo, hi)
    if f["kind"] == "linear":
        centre = [0.5 * (lo[a] + hi[a]) for a in range(3)]
        return vol * (sum(g * c for g, c in zip(f["gradient"], centre)) + f["offset"])
    a = f["axis"]
    t0, t1 = lo[a], hi[a]
    clamped = min(max(f["face"], t0), t1)
    len_low = max(0.0, clamped - t0)
    other = 1.0
    for b in range(3):
        if b != a:
            other *= hi[b] - lo[b]
    return f["low"] * len_low * other + f["high"] * (vol - len_low * other)


def cell_average(f, lo, hi):
    return box_integral(f, lo, hi) / volume_of(lo, hi)


def dir_grad(f, lo, hi, a):
    edge = hi[a] - lo[a]
    p = [0.5 * (lo[b] + hi[b]) for b in range(3)]
    p[a] = lo[a]
    v0 = field_at(f, p)
    p[a] = hi[a]
    v1 = field_at(f, p)
    return abs(v1 - v0) / edge


# ------------------------------------------------------ independent refinement
def split_bin(lo, hi, depth):
    edges = [lo, hi]
    for _ in range(depth):
        nxt = []
        for i in range(len(edges) - 1):
            nxt.append(edges[i])
            nxt.append(0.5 * (edges[i] + edges[i + 1]))
        nxt.append(edges[-1])
        edges = nxt
    return edges


def build_edges(lo, hi, ncell, depth):
    edges = []
    for i in range(ncell):
        a = lo + (hi - lo) * i / ncell
        b = lo + (hi - lo) * (i + 1) / ncell
        sub = split_bin(a, b, depth[i])
        for k in range(0 if i == 0 else 1, len(sub)):
            edges.append(sub[k])
    return edges


def base_bin(mid, lo, hi, n):
    b = int(math.floor((mid - lo) / (hi - lo) * n))
    return max(0, min(n - 1, b))


def refine_depths(case):
    lo, hi, ncell = case["lo"], case["hi"], case["cells"]
    f, thr, maxd = case["field"], case["threshold"], case["maxDepth"]
    depth = [[0] * n for n in ncell]
    eps = 1e-12 * max(1.0, thr)
    for _ in range(maxd):
        edges = [build_edges(lo[a], hi[a], ncell[a], depth[a]) for a in range(3)]
        bump = [[0] * n for n in ncell]
        for ix in range(len(edges[0]) - 1):
            for iy in range(len(edges[1]) - 1):
                for iz in range(len(edges[2]) - 1):
                    clo = [edges[0][ix], edges[1][iy], edges[2][iz]]
                    chi = [edges[0][ix + 1], edges[1][iy + 1], edges[2][iz + 1]]
                    g = [dir_grad(f, clo, chi, a) for a in range(3)]
                    if math.sqrt(sum(x * x for x in g)) <= thr:
                        continue
                    bins = [base_bin(0.5 * (clo[a] + chi[a]), lo[a], hi[a], ncell[a])
                            for a in range(3)]
                    for a in range(3):
                        if ncell[a] <= 1 or g[a] <= eps:
                            continue
                        b = bins[a]
                        if depth[a][b] < maxd:
                            bump[a][b] = 1
        changed = False
        for a in range(3):
            for b in range(ncell[a]):
                if bump[a][b]:
                    depth[a][b] += 1
                    changed = True
        if not changed:
            break
    return depth


def final_cells(case, depth):
    lo, hi, ncell = case["lo"], case["hi"], case["cells"]
    edges = [build_edges(lo[a], hi[a], ncell[a], depth[a]) for a in range(3)]
    cells = []
    for ix in range(len(edges[0]) - 1):
        for iy in range(len(edges[1]) - 1):
            for iz in range(len(edges[2]) - 1):
                clo = [edges[0][ix], edges[1][iy], edges[2][iz]]
                chi = [edges[0][ix + 1], edges[1][iy + 1], edges[2][iz + 1]]
                bins = [base_bin(0.5 * (clo[a] + chi[a]), lo[a], hi[a], ncell[a])
                        for a in range(3)]
                lvl = max(depth[a][bins[a]] for a in range(3))
                cells.append(dict(lo=clo, hi=chi, vol=volume_of(clo, chi),
                                  phi=cell_average(case["field"], clo, chi),
                                  level=lvl))
    return cells


def expectation(case):
    depth = refine_depths(case)
    cells = final_cells(case, depth)
    exact = box_integral(case["field"], case["lo"], case["hi"])
    vol = sum(c["vol"] for c in cells)
    level = max(max(d) for d in depth)
    return dict(cells=cells, n=len(cells), level=level, volume=vol,
                integral=sum(c["phi"] * c["vol"] for c in cells),
                exact=exact,
                minvol=min(c["vol"] for c in cells),
                maxvol=max(c["vol"] for c in cells))


# ------------------------------------------------------------------ app output
ADAPT_RE = re.compile(
    r"ADAPT level=(\d+)\s+cells=(\d+)\s+volumeSum=(" + NUM + r")\s+"
    r"integral=(" + NUM + r")\s+integralError=(" + NUM + r")\s+"
    r"minVolume=(" + NUM + r")\s+maxVolume=(" + NUM + r")")

CELL_RE = re.compile(
    r"ADAPT_CELL c=(\d+)\s+level=(\d+)\s+"
    r"lo=\((" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\)\s+"
    r"hi=\((" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\)\s+"
    r"phi=(" + NUM + r")\s+vol=(" + NUM + r")")

TOPO_RE = re.compile(r"adaptiveRefine: nPoints (\d+)\s+nCells (\d+)\s+"
                     r"nFaces (\d+)\s+nInternalFaces (\d+)\s+"
                     r"nBoundaryFaces (\d+)\s+nPatches (\d+)")


def run(exe, case):
    import subprocess
    proc = subprocess.run([exe, case], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=120)
    return proc.returncode, proc.stdout


def parse_report(text):
    m = ADAPT_RE.search(text)
    if not m:
        raise RuntimeError("app output has no ADAPT report line:\n" + text)
    t = TOPO_RE.search(text)
    topo = None
    if t:
        topo = dict(nPoints=int(t.group(1)), nCells=int(t.group(2)),
                    nFaces=int(t.group(3)), internal=int(t.group(4)),
                    boundary=int(t.group(5)), patches=int(t.group(6)))
    cells = []
    for cm in CELL_RE.finditer(text):
        cells.append(dict(c=int(cm.group(1)), level=int(cm.group(2)),
                          lo=[float(cm.group(i)) for i in (3, 4, 5)],
                          hi=[float(cm.group(i)) for i in (6, 7, 8)],
                          phi=float(cm.group(9)), vol=float(cm.group(10))))
    return dict(level=int(m.group(1)), cells=int(m.group(2)),
                volumeSum=float(m.group(3)), integral=float(m.group(4)),
                integralError=float(m.group(5)), minVolume=float(m.group(6)),
                maxVolume=float(m.group(7)), topo=topo, dumped=cells)


def run_ok(exe, case, tag):
    rc, out = run(exe, case)
    check(rc == 0, tag + ": app returned %d\n%s" % (rc, out))
    if rc != 0:
        return None, out
    try:
        return parse_report(out), out
    except RuntimeError as e:
        check(False, tag + ": " + str(e))
        return None, out


def canon_cell(c):
    return (tuple(round(v, 10) for v in c["lo"]),
            tuple(round(v, 10) for v in c["hi"]),
            c["level"],
            round(c["vol"], 12),
            round(c["phi"], 9))


def check_summary(rep, exp, case, tag):
    check(rep["cells"] == exp["n"],
          "%s: cells should be %d, got %d" % (tag, exp["n"], rep["cells"]))
    check(rep["level"] == exp["level"],
          "%s: level should be %d, got %d" % (tag, exp["level"], rep["level"]))
    boxvol = volume_of(case["lo"], case["hi"])
    check(abs(rep["volumeSum"] - boxvol) <= 1e-9,
          "%s: volumeSum %.12g should equal box volume %.12g" % (tag, rep["volumeSum"], boxvol))
    check(abs(rep["volumeSum"] - exp["volume"]) <= 1e-9,
          "%s: volumeSum %.12g differs from independent %.12g" % (tag, rep["volumeSum"], exp["volume"]))
    check(close(rep["integral"], exp["exact"], rel=1e-9, abs_=1e-11),
          "%s: integral %.12g should equal exact %.12g" % (tag, rep["integral"], exp["exact"]))
    check(rep["integralError"] <= 1e-9,
          "%s: integralError %.3g exceeds 1e-9" % (tag, rep["integralError"]))
    check(rep["minVolume"] > 0,
          "%s: minVolume must be positive, got %.12g" % (tag, rep["minVolume"]))
    check(rep["minVolume"] <= rep["maxVolume"],
          "%s: minVolume must not exceed maxVolume" % tag)
    check(abs(rep["minVolume"] - exp["minvol"]) <= 1e-9,
          "%s: minVolume %.12g differs from independent %.12g" % (tag, rep["minVolume"], exp["minvol"]))
    check(abs(rep["maxVolume"] - exp["maxvol"]) <= 1e-9,
          "%s: maxVolume %.12g differs from independent %.12g" % (tag, rep["maxVolume"], exp["maxvol"]))
    mean = rep["volumeSum"] / rep["cells"]
    # All six report numbers are printed at 12 significant digits; the mean
    # formed from volumeSum/cells can sit a fraction of a last digit outside
    # the printed min/max even when the underlying cell volumes are identical.
    tol = 1e-12 * max(1.0, abs(mean))
    check(rep["minVolume"] - tol <= mean <= rep["maxVolume"] + tol,
          "%s: mean cell volume must lie between min and max" % tag)
    if rep["topo"]:
        t = rep["topo"]
        check(t["nCells"] == rep["cells"],
              "%s: polyMesh nCells %d != report cells %d" % (tag, t["nCells"], rep["cells"]))
        check(t["nFaces"] == t["internal"] + t["boundary"],
              "%s: nFaces != internal + boundary" % tag)
        check(t["patches"] == 6,
              "%s: expected six boundary patches, got %d" % (tag, t["patches"]))


def rewrite(text, subs):
    out = text
    for old, new in subs:
        if old not in out:
            raise RuntimeError("cannot rewrite '%s'" % old)
        out = out.replace(old, new, 1)
    return out


def write_case(td, name, text):
    path = os.path.join(td, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    exe = os.path.abspath(sys.argv[1])
    maincase, uniform, reject = [os.path.abspath(p) for p in sys.argv[2:5]]

    # ---- 1) step field: full independent re-derivation -------------------
    case = parse_case(maincase)
    exp = expectation(case)
    rep, out = run_ok(exe, maincase, "step")
    if rep:
        check_summary(rep, exp, case, "step")
        boxvol = volume_of(case["lo"], case["hi"])
        check(abs(boxvol - 0.1) <= 1e-12, "closed form: box volume should be 0.1")
        # analytic integral: high*(1 - interface)*Ly*Lz = 1*0.45*1*0.1
        check(abs(exp["exact"] - 0.045) <= 1e-12,
              "closed form: exact step integral should be 0.045, got %.12g" % exp["exact"])
        check(exp["n"] > case["cells"][0] * case["cells"][1] * case["cells"][2],
              "step: refinement must increase the cell count")
        check(rep["level"] == case["maxDepth"],
              "step: refinement must reach maxDepth")

        # per-cell verification + localisation
        got = rep["dumped"]
        check(len(got) == rep["cells"],
              "step: dumped %d cells, report says %d" % (len(got), rep["cells"]))
        if len(got) == len(exp["cells"]):
            want_set = sorted(canon_cell(c) for c in exp["cells"])
            got_set = sorted(canon_cell(c) for c in got)
            check(got_set == want_set,
                  "step: per-cell geometry/phi disagrees with the independent refinement")
        axis = case["field"]["axis"]
        refined = [c for c in got if c["level"] > 0]
        check(len(refined) > 0, "step: no cell was refined")
        check(len(refined) < len(got), "step: refinement was not localised (all cells refined)")
        # only the base bin holding the interface may contain refined cells
        a, b = 0.5, 0.625
        slop = 1e-9
        check(all(c["lo"][axis] >= a - slop and c["hi"][axis] <= b + slop for c in refined),
              "step: refined cells escaped the interface bin [0.5, 0.625]")
        check(all(c["level"] == 0 for c in got
                  if c["lo"][axis] > b + slop or c["hi"][axis] < a - slop),
              "step: cells away from the interface must stay at level 0")

        # divergence-theorem volume from the per-cell polyhedra
        polyvol = sum(c["vol"] for c in got)
        check(abs(polyvol - boxvol) <= 1e-9,
              "step: per-cell volume sum %.12g != box volume %.12g" % (polyvol, boxvol))

    # ---- 2) uniform twin: must not refine -------------------------------
    ucase = parse_case(uniform)
    uexp = expectation(ucase)
    urep, _ = run_ok(exe, uniform, "uniform")
    if urep:
        check_summary(urep, uexp, ucase, "uniform")
        check(urep["level"] == 0, "uniform: level must be 0, got %d" % urep["level"])
        check(urep["cells"] == ucase["cells"][0] * ucase["cells"][1] * ucase["cells"][2],
              "uniform: cell count must equal the base grid")
        check(abs(urep["integral"] - 3.5 * 0.1) <= 1e-9,
              "uniform: integral should be 3.5*0.1 = 0.35")

    # ---- 3) generated linear cases: depth cap and threshold gate --------
    with tempfile.TemporaryDirectory(prefix="cf_adapt_check_") as td:
        base = ucase["text"]
        lin = rewrite(base, [("gradient (0 0 0);", "gradient (1 0 0);"),
                             ("gradientThreshold 3;", "gradientThreshold 0.5;")])
        lin_path = write_case(td, "linear.dict", lin)
        lcase = parse_case(lin_path)
        lexp = expectation(lcase)
        lrep, _ = run_ok(exe, lin_path, "linearX")
        if lrep:
            check_summary(lrep, lexp, lcase, "linearX")
            check(lrep["level"] == lcase["maxDepth"],
                  "linearX: all cells must refine up to maxDepth")
            check(lrep["cells"] == 64 * 8 * 1,
                  "linearX: expected 64 x-bins x 8 y-bins = 512 cells, got %d" % lrep["cells"])

        gate = rewrite(base, [("gradient (0 0 0);", "gradient (1 0 0);"),
                              ("gradientThreshold 3;", "gradientThreshold 10;")])
        gate_path = write_case(td, "gate.dict", gate)
        gcase = parse_case(gate_path)
        gexp = expectation(gcase)
        greport, _ = run_ok(exe, gate_path, "linearGate")
        if greport:
            check_summary(greport, gexp, gcase, "linearGate")
            check(greport["level"] == 0,
                  "linearGate: gradient below threshold must not refine")

        multi = rewrite(base, [("gradient (0 0 0);", "gradient (1 1 0);"),
                               ("gradientThreshold 3;", "gradientThreshold 0.5;")])
        multi_path = write_case(td, "multi.dict", multi)
        mcase = parse_case(multi_path)
        mexp = expectation(mcase)
        mreport, _ = run_ok(exe, multi_path, "linearXY")
        if mreport:
            check_summary(mreport, mexp, mcase, "linearXY")
            check(mreport["level"] == mcase["maxDepth"],
                  "linearXY: refinement must reach maxDepth")

        # ---- 4) loud rejections -----------------------------------------
        probes = [
            (reject, "negative threshold", "gradientThreshold"),
            (write_case(td, "badDepth.dict", rewrite(base, [
                ("maxDepth 3;", "maxDepth -1;")])), "negative depth", "maxDepth"),
            (write_case(td, "fracDepth.dict", rewrite(base, [
                ("maxDepth 3;", "maxDepth 1.5;")]), ), "fractional depth", "maxDepth"),
            (write_case(td, "bigDepth.dict", rewrite(base, [
                ("maxDepth 3;", "maxDepth 99;")]), ), "depth above cap", "maxDepth"),
            (write_case(td, "badField.dict", rewrite(base, [
                ("field linear;", "field quadratic;")])), "unknown field", "field"),
            (write_case(td, "badAxis.dict", rewrite(base, [
                ("field linear;", "field step;\n    axis w;\n    interface 0.5;")])),
             "unknown axis", "axis"),
        ]
        for path, label, token in probes:
            rc, out = run(exe, path)
            check(rc != 0, "reject %s: app must fail (rc=%d)" % (label, rc))
            check(token.lower() in out.lower(),
                  "reject %s: message must mention '%s'\n%s" % (label, token, out))

    notes.append("step: %d cells, level %d, exact integral %.12g"
                 % (exp["n"], exp["level"], exp["exact"]))
    for n in notes:
        print("  " + n)
    if fails:
        print("check_adaptive_refine: FAIL (%d of %d assertions)" % (len(fails), nchecks))
        for m in fails:
            print("  - " + m)
        return 1
    print("check_adaptive_refine: PASS (%d assertions)" % nchecks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
