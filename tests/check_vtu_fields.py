#!/usr/bin/env python3
"""VTU field-set and value-alignment verification (batch 2, item 3).

Two groups:

  A. chflow-piso-foam VTU补场: a field is written only when the solver really
     built it.
       * laminar / isothermal case  -> vtu holds U, p only (no T, k, nut,
         epsilon, omega, gamma, desLength);
       * energy-only case           -> U, p, T; still no turbulence fields;
       * k-epsilon + energy case   -> U, p, nut, k, epsilon, T (no omega);
       * explicit `turbulence { model laminar; }` -> still no k/nut placeholder.
     Where the solver prints a range or a maximum (Tmin/Tmax, uMax), the vtu
     array's own min/max must be the same number.

  B. chflow-1d-steady VTU 落盘: one VTK_VERTEX per pipe at its centre, with
     x / p / T / alpha / regimeIndex. Rows align to the stdout `pipe <name> ...`
     lines in declaration order. A case that declares no `output { vtk ...; }`
     writes no file at all.

  Regime encoding (mirrors steadyNetworkApp.cpp::regimeName):
     -1 none / single-phase (no two-phase law priced)
      0 segregated   1 transition   2 intermittent   3 distributed

usage: check_vtu_fields.py <chflow-piso-foam> <chflow-1d-steady> [cases-root]
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

REGIME_TO_INDEX = {"segregated": 0, "transition": 1, "intermittent": 2,
                   "distributed": 3}

# Turbulence / energy field names that must never appear in a case that did not
# build the corresponding model.
TURB_FIELDS = {"k", "nut", "epsilon", "omega", "gamma", "desLength"}
ENERGY_FIELDS = {"T"}


# ---------------------------------------------------------------------------
# minimal ASCII vtu parser (no numpy / vtk dependency)
# ---------------------------------------------------------------------------
def read_vtu(path):
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
    conn = [int(v) for v in raw("connectivity")]
    offs = [int(v) for v in raw("offsets")]
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
        vals = [float(x) for x in body.split()]
        if len(vals) == nf:
            fields[m.group(1)] = vals
        elif len(vals) == 3 * nf:
            fields[m.group(1)] = [tuple(vals[3 * i:3 * i + 3])
                                  for i in range(nf)]
    return cent, fields


def run_case(binary, case_path, work, inject=None):
    """Copy the case into `work`, run the solver there, return (stdout, vtu or
    None). The vtu path is read back from the local dictionary so a relative
    landing dir matches what the solver was given."""
    local = os.path.join(work, os.path.basename(case_path))
    shutil.copy(case_path, local)
    if inject is not None:
        content = open(local).read()
        open(local, "w").write(inject(content))
    r = subprocess.run([binary, local], capture_output=True, text=True, cwd=work)
    if r.returncode != 0:
        raise SystemExit(f"solver failed rc={r.returncode}\n{r.stdout}\n{r.stderr}")
    m = re.search(r'^\s*vtk\s+(\S+);', open(local).read(), re.M)
    vtu = os.path.join(work, m.group(1)) if m else None
    return r.stdout, vtu


def parse_pipe_decls(txt):
    """name -> (from, to, length) for every pipe block in the dictionary.
    Comments are stripped first: a `//` line above a pipe can otherwise make
    the brace regex swallow the pipe's own opening brace."""
    clean = "\n".join(re.sub(r"//.*", "", ln) for ln in txt.splitlines())
    decl = {}
    for m in re.finditer(
            r'(\w+)\s*\{\s*from\s+(\w+);\s*to\s+(\w+);[^}]*?'
            r'length\s+([0-9.]+)', clean):
        decl[m.group(1)] = (m.group(2), m.group(3), float(m.group(4)))
    return decl


# ---------------------------------------------------------------------------
# group A: 3D pisoFoamApp
# ---------------------------------------------------------------------------
def check_3d(piso, cases_root, fails):
    def expect_fields(tag, vtu, required, forbidden):
        _, f = read_vtu(vtu)
        names = set(f.keys())
        for req in required:
            if req not in names:
                fails.append(f"[A {tag}] required field '{req}' missing; "
                             f"vtu fields = {sorted(names)}")
        for bad in forbidden:
            if bad in names:
                fails.append(f"[A {tag}] field '{bad}' must not be written "
                             f"(field set = {sorted(names)})")
        return f

    # A1: laminar, isothermal, no turbulence / energy block at all.
    work = tempfile.mkdtemp(prefix="cf_vtu_a1_")
    out, vtu = run_case(piso, os.path.join(cases_root, "startedChannel",
                                           "case.dict"), work)
    if vtu is None:
        fails.append("[A A1] startedChannel declared no vtk output")
    else:
        f = expect_fields("A1 laminar/isothermal", vtu, {"U", "p"},
                          TURB_FIELDS | ENERGY_FIELDS)
        umax = float(re.findall(r"uMax (\S+)", out)[-1])
        umax_vtu = max(mag(u) for u in f["U"])
        if not math.isclose(umax, umax_vtu, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A1] uMax stdout {umax} != vtu max |U| "
                         f"{umax_vtu}")

    # A2: energy only, no turbulence.
    work = tempfile.mkdtemp(prefix="cf_vtu_a2_")
    out, vtu = run_case(piso, os.path.join(cases_root, "energyDecay",
                                           "case.dict"), work)
    if vtu is None:
        fails.append("[A A2] energyDecay declared no vtk output")
    else:
        f = expect_fields("A2 energy-only", vtu, {"U", "p", "T"}, TURB_FIELDS)
        t = list(re.finditer(r"energy: Tmin (\S+) Tmax (\S+)", out))[-1]
        tmin, tmax = float(t.group(1)), float(t.group(2))
        Tv = f["T"]
        if not math.isclose(min(Tv), tmin, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A2] Tmin stdout {tmin} != vtu {min(Tv)}")
        if not math.isclose(max(Tv), tmax, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A2] Tmax stdout {tmax} != vtu {max(Tv)}")

    # A3: k-epsilon turbulence plus energy.
    work = tempfile.mkdtemp(prefix="cf_vtu_a3_")
    out, vtu = run_case(piso, os.path.join(cases_root, "turbulentChannel",
                                           "case.dict"), work)
    if vtu is None:
        fails.append("[A A3] turbulentChannel declared no vtk output")
    else:
        f = expect_fields("A3 kEpsilon+energy", vtu,
                           {"U", "p", "nut", "k", "epsilon", "T"}, {"omega"})
        for scal in ("k", "nut", "epsilon", "T", "p"):
            bad = [v for v in f[scal] if not math.isfinite(v)]
            if bad:
                fails.append(f"[A A3] field '{scal}' has {len(bad)} non-finite "
                             f"entries")
        if any(v < -1e-12 for v in f["nut"]) or any(v < -1e-12 for v in f["k"]):
            fails.append("[A A3] nut/k went negative")
        t = list(re.finditer(r"energy: Tmin (\S+) Tmax (\S+)", out))[-1]
        tmin, tmax = float(t.group(1)), float(t.group(2))
        Tv = f["T"]
        if not math.isclose(min(Tv), tmin, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A3] Tmin stdout {tmin} != vtu {min(Tv)}")
        if not math.isclose(max(Tv), tmax, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A3] Tmax stdout {tmax} != vtu {max(Tv)}")
        umax = float(re.findall(r"uMax (\S+)", out)[-1])
        umax_vtu = max(mag(u) for u in f["U"])
        if not math.isclose(umax, umax_vtu, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[A A3] uMax stdout {umax} != vtu max |U| "
                         f"{umax_vtu}")

    # A4: explicit laminar model block: the model object exists but carries no
    # transported field, so no zero-placeholder nut/k may be written.
    def inject_laminar(txt):
        return txt.replace("nu 0.01;",
                           "nu 0.01;\n\nturbulence\n{\n    model laminar;\n}")
    work = tempfile.mkdtemp(prefix="cf_vtu_a4_")
    out, vtu = run_case(piso, os.path.join(cases_root, "startedChannel",
                                           "case.dict"), work, inject_laminar)
    if vtu is None:
        fails.append("[A A4] laminar-injected case declared no vtk output")
    else:
        expect_fields("A4 explicit laminar", vtu, {"U", "p"},
                      TURB_FIELDS | ENERGY_FIELDS)


def mag(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


# ---------------------------------------------------------------------------
# group B: 1D steady network
# ---------------------------------------------------------------------------
def check_1d(steady, cases_root, fails):
    # B1: two-phase network that declares the vtk dump.
    case = os.path.join(cases_root, "steadyNetworkVtu", "case.dict")
    work = tempfile.mkdtemp(prefix="cf_vtu_b1_")
    out, vtu = run_case(steady, case, work)
    if vtu is None:
        fails.append("[B B1] steadyNetworkVtu declared no vtk output")
        return
    if not os.path.exists(vtu):
        fails.append(f"[B B1] vtu {vtu} was not written")
        return

    cent, f = read_vtu(vtu)
    n = len(cent)
    for key in ("x", "p", "T", "alpha", "regimeIndex"):
        if key not in f:
            fails.append(f"[B B1] required field '{key}' missing; "
                         f"vtu fields = {sorted(f.keys())}")
        elif len(f[key]) != n:
            fails.append(f"[B B1] field '{key}' has {len(f[key])} entries, "
                         f"expected {n} (pipe centres)")

    nodes, pipes = {}, []
    for line in out.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node":
            nodes[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            d = {"name": tok[1]}
            it = iter(tok[2:])
            for a, b in zip(it, it):
                try:
                    d[a] = float(b)
                except ValueError:
                    d[a] = b
            pipes.append(d)

    if n != len(pipes):
        fails.append(f"[B B1] vtu rows {n} != stdout pipe lines {len(pipes)}")
        return

    decl = parse_pipe_decls(open(case).read())
    feed_t = float(re.search(r"temperature\s+([0-9.]+)",
                             open(case).read()).group(1))

    arc = 0.0
    for i, row in enumerate(pipes):
        name = row["name"]
        if name not in decl:
            fails.append(f"[B B1] pipe {name} not parsed from dictionary")
            continue
        frm, to, length = decl[name]
        x_exp = arc + 0.5 * length
        arc += length

        if not math.isclose(f["x"][i], x_exp, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[B B1] {name}: x vtu {f['x'][i]} != centre "
                         f"{x_exp}")

        # centre pressure = mean of the two solved node pressures; the pipe
        # drop itself (dpTotal = p(from)-p(to)) is the stdout `dp` token, not
        # duplicated in the vtu.
        p_exp = 0.5 * (nodes[frm] + nodes[to])
        if not math.isclose(f["p"][i], p_exp, rel_tol=1e-9, abs_tol=1e-6):
            fails.append(f"[B B1] {name}: centre p vtu {f['p'][i]} != mean "
                         f"node pressure {p_exp}")

        if "HL" in row:
            if not math.isclose(f["alpha"][i], row["HL"], rel_tol=1e-9,
                                abs_tol=1e-12):
                fails.append(f"[B B1] {name}: alpha vtu {f['alpha'][i]} != "
                             f"stdout HL {row['HL']}")
            if row["regime"] not in REGIME_TO_INDEX:
                fails.append(f"[B B1] {name}: unknown regime "
                             f"'{row['regime']}'")
            else:
                idx_exp = REGIME_TO_INDEX[row["regime"]]
                if int(f["regimeIndex"][i]) != idx_exp:
                    fails.append(f"[B B1] {name}: regimeIndex vtu "
                                 f"{int(f['regimeIndex'][i])} != {idx_exp} "
                                 f"({row['regime']})")
        else:
            if f["alpha"][i] != 0.0:
                fails.append(f"[B B1] {name}: single-phase pipe alpha "
                             f"{f['alpha'][i]} != 0")
            if int(f["regimeIndex"][i]) != -1:
                fails.append(f"[B B1] {name}: single-phase pipe regimeIndex "
                             f"{int(f['regimeIndex'][i])} != -1")

        if not math.isclose(f["T"][i], feed_t, rel_tol=1e-9, abs_tol=1e-12):
            fails.append(f"[B B1] {name}: T vtu {f['T'][i]} != feed "
                         f"temperature {feed_t}")

    # B2: no output declaration -> no file, stdout untouched.
    work2 = tempfile.mkdtemp(prefix="cf_vtu_b2_")
    out2, vtu2 = run_case(steady, os.path.join(cases_root, "steadyNetwork",
                                                "case.dict"), work2)
    if vtu2 is not None:
        fails.append("[B B2] undeclared case resolved a vtk path")
    stray = [f for f in os.listdir(work2) if f.endswith(".vtu")]
    if stray:
        fails.append(f"[B B2] undeclared case wrote vtu file(s): {stray}")
    if "steady solved" not in out2:
        fails.append("[B B2] baseline steadyNetwork case did not solve")


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    piso, steady = sys.argv[1], sys.argv[2]
    cases_root = sys.argv[3] if len(sys.argv) > 3 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases")

    fails = []
    check_3d(piso, cases_root, fails)
    check_1d(steady, cases_root, fails)

    if fails:
        for msg in fails:
            print("FAIL:", msg)
        raise SystemExit("vtu field verification FAILED")

    print("vtu field verification OK: 3D field sets and stdout alignment "
          "match; 1D network vtu rows align to the pipe report; no vtu on the "
          "undeclared baseline case.")


if __name__ == "__main__":
    main()
