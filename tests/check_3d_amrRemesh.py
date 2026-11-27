#!/usr/bin/env python3
"""5.237 (H6): independent acceptance test for the opt-in AMR remeshing loop.

The checker mirrors the whole loop in pure Python from the case dictionaries
(seed depths from adapt{}, explicit conservative diffusion, per-axis slab
refine / coarsen indicator, conservative box-intersection field transfer and
the uniform conservation shift) and compares every AMR_CELL / AMR_CYCLE /
AMR_DONE number. It asserts:

  1) refine arm (cycles 1): bins 3,4,5 deepen to 1 (88 cells); each parent's
     child volumes sum to the parent volume (1e-14), the two children of one
     parent share the parent seed phi plus one identical round-off shift
     (|shift| < 1e-12);
  2) coarsen arm (cycles 2): after 400 diffusion steps the merged parent phi
     equals the volume-weighted average of its two (diffused) children
     (1e-10 against the mirror) and the count returns to 64;
  3) loop arm (cycles 3): cycle 3 does nothing, finalCells == initialCells,
     massError <= 1e-10 and the mass stays at the exact step integral 0.045;
  4) passive twin: without amr{} no AMR line is printed and the output is
     byte-identical to the amr-stripped case;
  5) loud refusals: coarsenThreshold >= refineThreshold, a negative
     refineThreshold and cycles < 1.

usage: check_3d_amrRemesh.py <exe> <case> <caseRefine> <caseCoarsen> <caseTwin> <caseReject>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
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


def block_span(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", txt)
    if not m:
        return None
    start = m.start()
    depth = 0
    for j in range(m.end() - 1, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return start, j + 1
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
    f = dict(kind=word_of(adapt, "field"))
    if f["kind"] == "step":
        f["axis"] = "xyz".index(word_of(adapt, "axis", "x"))
        f["face"] = scalar_of(adapt, "interface")
        f["low"] = scalar_of(adapt, "low", 0.0)
        f["high"] = scalar_of(adapt, "high", 1.0)
    elif f["kind"] == "linear":
        f["gradient"] = triple_of(adapt, "gradient")
        f["offset"] = scalar_of(adapt, "offset", 0.0)
    else:
        raise RuntimeError("unknown field '%s'" % f["kind"])
    amr = dict(has=False)
    span = block_span(txt, "amr")
    if span:
        a = txt[span[0]:span[1]]
        amr = dict(has=True,
                   refine=scalar_of(a, "refineThreshold"),
                   coarsen=scalar_of(a, "coarsenThreshold"),
                   cycles=int(scalar_of(a, "cycles")),
                   every=int(scalar_of(a, "every", 1)),
                   diffusion=scalar_of(a, "diffusion", 0.0),
                   dt=scalar_of(a, "dt", 0.0))
    return dict(lo=lo, hi=hi, cells=cells, field=f,
                threshold=scalar_of(adapt, "gradientThreshold"),
                maxDepth=int(scalar_of(adapt, "maxDepth")),
                amr=amr, text=txt)


# ------------------------------------------------------------- field mirror
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


# ------------------------------------------------- seed refinement (F13 law)
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
                        if depth[a][bins[a]] < maxd:
                            bump[a][bins[a]] = 1
        changed = False
        for a in range(3):
            for b in range(ncell[a]):
                if bump[a][b]:
                    depth[a][b] += 1
                    changed = True
        if not changed:
            break
    return depth


# ------------------------------------------------------- AMR loop (mirror)
def cell_box(e, c):
    lo = [e[a][c[a]] for a in range(3)]
    hi = [e[a][c[a] + 1] for a in range(3)]
    return lo, hi


def run_mirror(case):
    """Re-derive the whole remesh loop from the case, without the C++ app."""
    lo, hi, ncell = case["lo"], case["hi"], case["cells"]
    amr = case["amr"]
    depth = refine_depths(case)
    e = [build_edges(lo[a], hi[a], ncell[a], depth[a]) for a in range(3)]
    nc = [len(e[a]) - 1 for a in range(3)]
    phi = {}
    for c in [(ix, iy, iz) for ix in range(nc[0]) for iy in range(nc[1])
              for iz in range(nc[2])]:
        clo, chi = cell_box(e, c)
        phi[c] = cell_average(case["field"], clo, chi)
    box_vol = volume_of(lo, hi)
    mass_target = sum(p * volume_of(*cell_box(e, c)) for c, p in phi.items())

    def neighbours(c, a):
        for d_ in (-1, 1):
            j = c[a] + d_
            if 0 <= j < nc[a]:
                nb = list(c)
                nb[a] = j
                yield tuple(nb), j

    def diffuse_once():
        new = {}
        for c, p in phi.items():
            flux = 0.0
            clo, chi = cell_box(e, c)
            for a in range(3):
                hc = chi[a] - clo[a]
                area = 1.0
                for b in range(3):
                    if b != a:
                        area *= chi[b] - clo[b]
                for nb, j in neighbours(c, a):
                    nlo, nhi = cell_box(e, nb)
                    d = 0.5 * (hc + (nhi[a] - nlo[a]))
                    flux += amr["diffusion"] * area * (phi[nb] - p) / d
            new[c] = p + amr["dt"] * flux / volume_of(clo, chi)
        return new

    def indicator():
        g = {}
        for c, p in phi.items():
            clo, chi = cell_box(e, c)
            ga = []
            for a in range(3):
                hc = chi[a] - clo[a]
                gmax = 0.0
                for nb, j in neighbours(c, a):
                    nlo, nhi = cell_box(e, nb)
                    d = 0.5 * (hc + (nhi[a] - nlo[a]))
                    gmax = max(gmax, abs(phi[nb] - p) / d)
                ga.append(gmax)
            g[c] = ga
        return g

    cycles = []
    for cyc in range(1, amr["cycles"] + 1):
        if cyc > 1 and amr["diffusion"] > 0:
            for _ in range(amr["every"]):
                phi = diffuse_once()
        # the diffused state of the CURRENT cells, i.e. what the transfer of
        # this cycle consumes
        pre_diff = dict(phi)
        pre_e = [list(x) for x in e]
        g = indicator()
        n_ref = n_coa = 0
        for a in range(3):
            bin_max = [0.0] * ncell[a]
            for c in phi:
                clo, chi = cell_box(e, c)
                b = base_bin(0.5 * (clo[a] + chi[a]), lo[a], hi[a], ncell[a])
                bin_max[b] = max(bin_max[b], g[c][a])
            for b in range(ncell[a]):
                if bin_max[b] > amr["refine"] and depth[a][b] < case["maxDepth"]:
                    depth[a][b] += 1
                    n_ref += 1
                elif bin_max[b] < amr["coarsen"] and depth[a][b] > 0:
                    depth[a][b] -= 1
                    n_coa += 1
        old_e, old_nc, old_phi = e, nc, phi
        e = [build_edges(lo[a], hi[a], ncell[a], depth[a]) for a in range(3)]
        nc = [len(e[a]) - 1 for a in range(3)]
        phi = {}
        for c in [(ix, iy, iz) for ix in range(nc[0]) for iy in range(nc[1])
                  for iz in range(nc[2])]:
            clo, chi = cell_box(e, c)
            vnew = volume_of(clo, chi)
            acc = 0.0
            for oc in old_phi:
                olo, ohi = cell_box(old_e, oc)
                vint = 1.0
                for a in range(3):
                    ov = min(chi[a], ohi[a]) - max(clo[a], olo[a])
                    if ov <= 0:
                        vint = 0.0
                        break
                    vint *= ov
                if vint > 0:
                    acc += old_phi[oc] * vint
            phi[c] = acc / vnew
        mass = sum(p * volume_of(*cell_box(e, c)) for c, p in phi.items())
        shift = (mass_target - mass) / box_vol
        for c in phi:
            phi[c] += shift
        cells = []
        for c, p in phi.items():
            clo, chi = cell_box(e, c)
            lvl = max(depth[a][base_bin(0.5 * (clo[a] + chi[a]), lo[a], hi[a], ncell[a])]
                      for a in range(3))
            cells.append(dict(lo=clo, hi=chi, level=lvl, phi=p,
                              vol=volume_of(clo, chi)))
        cycles.append(dict(cycle=cyc, cells=cells, n=len(cells),
                           refined=n_ref, coarsened=n_coa,
                           mass=sum(cd["phi"] * cd["vol"] for cd in cells),
                           pre_diff=pre_diff, pre_e=pre_e))
    return dict(cycles=cycles, mass_target=mass_target, box_vol=box_vol,
                seed=dict(phi=phi))


# ------------------------------------------------------------- app parsing
CYCLE_RE = re.compile(r"AMR_CYCLE cycle=(\d+)\s+cells=(\d+)\s+refined=(\d+)\s+"
                      r"coarsened=(\d+)\s+mass=(" + NUM + r")")
DONE_RE = re.compile(r"AMR_DONE initialCells=(\d+)\s+finalCells=(\d+)\s+"
                     r"massInitial=(" + NUM + r")\s+massFinal=(" + NUM + r")\s+"
                     r"massError=(" + NUM + r")")
CELL_RE = re.compile(r"AMR_CELL cycle=(\d+)\s+c=(\d+)\s+level=(\d+)\s+"
                     r"lo=\((" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\)\s+"
                     r"hi=\((" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\)\s+"
                     r"phi=(" + NUM + r")\s+vol=(" + NUM + r")")


def run(exe, case):
    proc = subprocess.run([exe, case], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=300)
    return proc.returncode, proc.stdout


def run_text(exe, text, tag):
    with tempfile.TemporaryDirectory(prefix="cf_amr_") as td:
        p = os.path.join(td, tag + ".dict")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        return run(exe, p)


def parse_amr(out):
    cycles = {}
    for m in CYCLE_RE.finditer(out):
        cycles[int(m.group(1))] = dict(cells=int(m.group(2)),
                                       refined=int(m.group(3)),
                                       coarsened=int(m.group(4)),
                                       mass=float(m.group(5)))
    dumped = {}
    for m in CELL_RE.finditer(out):
        k = int(m.group(1))
        dumped.setdefault(k, []).append(dict(
            level=int(m.group(3)),
            lo=[float(m.group(i)) for i in (4, 5, 6)],
            hi=[float(m.group(i)) for i in (7, 8, 9)],
            phi=float(m.group(10)), vol=float(m.group(11))))
    d = DONE_RE.search(out)
    done = None
    if d:
        done = dict(initial=int(d.group(1)), final=int(d.group(2)),
                    m0=float(d.group(3)), m1=float(d.group(4)),
                    err=float(d.group(5)))
    return cycles, dumped, done


def canon(c):
    return (tuple(round(v, 10) for v in c["lo"]),
            tuple(round(v, 10) for v in c["hi"]),
            c["level"], round(c["vol"], 12), round(c["phi"], 8))


def compare_mirror(case, path, tag):
    """Full per-cycle comparison of the app against the Python mirror."""
    rc, out = run(EXE, path)
    check(rc == 0, "%s: app rc=%d\n%s" % (tag, rc, out[-300:]))
    if rc != 0:
        return None, None
    cycles, dumped, done = parse_amr(out)
    mir = run_mirror(case)
    check(len(cycles) == len(mir["cycles"]),
          "%s: app has %d cycles, mirror %d" % (tag, len(cycles), len(mir["cycles"])))
    for mc in mir["cycles"]:
        k = mc["cycle"]
        if k not in cycles or k not in dumped:
            check(False, "%s: cycle %d missing from app" % (tag, k))
            continue
        ap = cycles[k]
        check(ap["cells"] == mc["n"],
              "%s c%d: app cells %d != mirror %d" % (tag, k, ap["cells"], mc["n"]))
        check(ap["refined"] == mc["refined"],
              "%s c%d: app refined %d != mirror %d" % (tag, k, ap["refined"], mc["refined"]))
        check(ap["coarsened"] == mc["coarsened"],
              "%s c%d: app coarsened %d != mirror %d" % (tag, k, ap["coarsened"], mc["coarsened"]))
        check(abs(ap["mass"] - mc["mass"]) <= 1e-9,
              "%s c%d: app mass %.12g != mirror %.12g" % (tag, k, ap["mass"], mc["mass"]))
        want = sorted(canon(c) for c in mc["cells"])
        got = sorted(canon(c) for c in dumped[k])
        check(want == got, "%s c%d: per-cell geometry/phi set differs from the mirror"
              % (tag, k))
    return mir, (cycles, dumped, done)


def main():
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    global EXE
    EXE = os.path.abspath(sys.argv[1])
    paths = [os.path.abspath(p) for p in sys.argv[2:7]]
    maincase, refinecase, coarsencase, twincase, rejectcase = paths

    # ---- 3) loop arm: full mirror + qualitative expectations --------------
    case = parse_case(maincase)
    mir, got = compare_mirror(case, maincase, "loop")
    if got:
        cycles, dumped, done = got
        check(done is not None, "loop: no AMR_DONE line")
        if done:
            check(done["initial"] == 64 and done["final"] == 64,
                  "loop: initial/final cells %d/%d must both be 64"
                  % (done["initial"], done["final"]))
            check(done["err"] <= 1e-10, "loop: massError %.3g > 1e-10" % done["err"])
            check(abs(done["m0"] - 0.045) <= 1e-12 and abs(done["m1"] - 0.045) <= 1e-12,
                  "loop: masses %.12g/%.12g != exact step integral 0.045"
                  % (done["m0"], done["m1"]))
        if 1 in cycles:
            check(cycles[1]["cells"] == 88 and cycles[1]["refined"] == 3,
                  "loop: cycle1 must refine bins 3,4,5 to 88 cells, got %d/%d"
                  % (cycles[1]["cells"], cycles[1]["refined"]))
        if 2 in cycles:
            check(cycles[2]["cells"] == 64 and cycles[2]["coarsened"] == 3,
                  "loop: cycle2 must coarsen back to 64 cells, got %d/%d"
                  % (cycles[2]["cells"], cycles[2]["coarsened"]))
        if 3 in cycles:
            check(cycles[3]["refined"] == 0 and cycles[3]["coarsened"] == 0
                  and cycles[3]["cells"] == 64,
                  "loop: cycle3 must be inert")
        for k, ap in cycles.items():
            check(abs(ap["mass"] - 0.045) <= 1e-12,
                  "loop c%d: mass %.12g != 0.045" % (k, ap["mass"]))
        notes.append("loop: %s" % ", ".join(
            "c%d=%d cells (+%d/-%d)" % (k, v["cells"], v["refined"], v["coarsened"])
            for k, v in sorted(cycles.items())))

    # ---- 1) refine arm: child volumes and single-parent copy --------------
    rcase = parse_case(refinecase)
    rmir, rgot = compare_mirror(rcase, refinecase, "refine")
    if rgot:
        _, rdumped, rdone = rgot
        check(rdone["final"] == 88, "refine: final cells %d != 88" % rdone["final"])
        seed = rcase  # 64 seed cells from the closed-form step averages
        sdepth = refine_depths(seed)
        se = [build_edges(seed["lo"][a], seed["hi"][a], seed["cells"][a], sdepth[a])
              for a in range(3)]
        snc = [len(x) - 1 for x in se]
        parents = {}
        for c in rdumped.get(1, []):
            key = tuple(base_bin(0.5 * (c["lo"][a] + c["hi"][a]), seed["lo"][a],
                                 seed["hi"][a], snc[a]) for a in range(3))
            parents.setdefault(key, []).append(c)
        n_parents = snc[0] * snc[1] * snc[2]
        check(len(parents) == n_parents,
              "refine: %d parent groups != %d seed cells" % (len(parents), n_parents))
        two_child = 0
        for key, kids in parents.items():
            clo = [se[a][key[a]] for a in range(3)]
            chi = [se[a][key[a] + 1] for a in range(3)]
            pv = volume_of(clo, chi)
            check(abs(sum(k["vol"] for k in kids) - pv) <= 1e-14,
                  "refine: child volumes of parent %s sum to %.14g != %.14g"
                  % (key, sum(k["vol"] for k in kids), pv))
            pphi = cell_average(seed["field"], clo, chi)
            if len(kids) == 2:
                two_child += 1
                check(abs(kids[0]["phi"] - kids[1]["phi"]) <= 1e-12,
                      "refine: children of %s differ: %.12g vs %.12g"
                      % (key, kids[0]["phi"], kids[1]["phi"]))
                shift = kids[0]["phi"] - pphi
                check(abs(shift) < 1e-12,
                      "refine: child phi %.12g != parent seed %.12g (shift %.3g)"
                      % (kids[0]["phi"], pphi, shift))
        check(two_child == 24, "refine: expected 3 bins x 8 y = 24 split parents, got %d"
              % two_child)

    # ---- 2) coarsen arm: volume-weighted merge ----------------------------
    ccase = parse_case(coarsencase)
    cmir, cgot = compare_mirror(ccase, coarsencase, "coarsen")
    if cgot and cmir:
        _, cdumped, cdone = cgot
        check(cdone["final"] == 64, "coarsen: final cells %d != 64" % cdone["final"])
        check(cdumped.get(2) and len(cdumped[2]) == 64, "coarsen: cycle2 must have 64 cells")
        mc1 = cmir["cycles"][0]
        mc2 = cmir["cycles"][1]
        # the children of the merge, AFTER the diffusion that precedes cycle 2
        pdiff, pe = mc2["pre_diff"], mc2["pre_e"]
        n_checked = 0
        for c in cdumped.get(2, []):
            kids = []
            for oc, op in pdiff.items():
                olo = [pe[a][oc[a]] for a in range(3)]
                ohi = [pe[a][oc[a] + 1] for a in range(3)]
                if all(olo[a] >= c["lo"][a] - 1e-12 and ohi[a] <= c["hi"][a] + 1e-12
                       for a in range(3)):
                    kids.append((op, volume_of(olo, ohi)))
            if len(kids) == 2:
                n_checked += 1
                vsum = sum(v for _, v in kids)
                acc = sum(p * v for p, v in kids)
                check(abs(vsum - c["vol"]) <= 1e-14,
                      "coarsen: merged children volumes %.14g != parent %.14g"
                      % (vsum, c["vol"]))
                check(abs(c["phi"] - acc / vsum) <= 1e-10,
                      "coarsen: parent phi %.12g != weighted child mean %.12g"
                      % (c["phi"], acc / vsum))
        check(n_checked == 24, "coarsen: checked %d merged parents, expected 24"
              % n_checked)

    # ---- 4) passive twin ---------------------------------------------------
    trc, tout = run(EXE, twincase)
    check(trc == 0, "twin: app rc=%d" % trc)
    check("AMR" not in tout, "twin: AMR lines printed without an amr{} block")
    # strip comments first: the case header mentions amr{}, and a raw span
    # search would remove the wrong block
    src = case["text"]
    s = block_span(src, "amr")
    stripped_raw = src[:s[0]].rstrip() + "\n" + src[s[1]:].lstrip()
    src2, out2 = run_text(EXE, stripped_raw, "stripped")
    check(src2 == 0 and out2 == tout,
          "twin: amr-stripped case is not byte-identical to caseTwin.dict")
    check("AMR" not in out2, "twin: stripped case printed AMR lines")

    # ---- 5) loud refusals ---------------------------------------------------
    rc, out = run(EXE, rejectcase)
    check(rc != 0 and "coarsenThreshold must be below refineThreshold" in out,
          "reject coarsen>=refine not raised loudly (rc=%d)" % rc)
    probes = [
        (re.sub(r"refineThreshold(\s+)" + NUM + r";", r"refineThreshold\g<1>-1;",
                src, count=1), "refineThreshold must be a finite positive"),
        (re.sub(r"cycles(\s+)" + NUM + r";", r"cycles\g<1>0;", src, count=1),
         "cycles must be an integer"),
        (re.sub(r"every(\s+)" + NUM + r";", r"every\g<1>0;", src, count=1),
         "every must be an integer"),
    ]
    for txt, needle in probes:
        check(txt != src, "reject probe did not rewrite: " + needle)
        rr, ro = run_text(EXE, txt, "rej")
        check(rr != 0 and needle in ro,
              "reject '%s' not raised loudly (rc=%d)" % (needle, rr))

    for n in notes:
        print("  " + n)
    if fails:
        print("check_3d_amrRemesh: FAIL (%d of %d assertions)" % (len(fails), nchecks))
        for m in fails[:40]:
            print("  - " + m)
        return 1
    print("check_3d_amrRemesh: PASS (%d assertions: full per-cycle mirror of "
          "diffusion/indicator/transfer/shift, refine+coarsen arms, byte-"
          "identical passive twin, 4 loud rejections)" % nchecks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
