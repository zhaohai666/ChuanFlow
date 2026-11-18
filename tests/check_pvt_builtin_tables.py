#!/usr/bin/env python3
"""Increment D (5.110): opt-in built-in parachor table (Weinaug-Katz 1943).

Verifies tests/cases/pvtBuiltinTables/case.dict, which declares
`surfaceTension { parachors fromTable; }` instead of an explicit list:

1.  Independently transcribes the published Weinaug-Katz (1943) pure-substance
    parachors (Table 18.1 of Adewumi, Penn State; original: Weinaug & Katz,
    Ind. Eng. Chem. 35, 239 (1943)) and applies the SI conversion
    P_SI = P_cgs * 10^-6.75 exactly as SurfaceTension documents.
2.  Re-derives sigma from the reported PR phase state (Vml = ML/rhoL,
    Vmg = MV/rhoV, x, y) using the table parachors and compares to the kernel.
3.  Gated twin: a case with the same feed declaring the table-derived SI
    parachors explicitly must produce byte-identical output except for the
    `fromTable` provenance line.
4.  Refusals: `fromTable` with a component absent from the published table
    (water) is a loud error; a non-list non-`fromTable` parachors entry is
    rejected; kij stays declared-only (no built-in table is consulted even
    when binaryInteraction is absent).
"""
import re
import subprocess
import sys
import tempfile
import os

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

# Published Weinaug-Katz (1943) parachors, cgs [cm3/mol (dyn/cm)^0.25].
WK_CGS = {
    "carbonDioxide": 78.0,
    "nitrogen": 41.0,
    "methane": 77.0,
    "ethane": 108.0,
    "propane": 150.3,
    "iButane": 181.5,
    "nButane": 189.9,
    "iPentane": 225.0,
    "nPentane": 231.5,
    "nHexane": 271.0,
    "nHeptane": 312.5,
    "nOctane": 351.5,
}
K_CGS_TO_SI = 10.0 ** -6.75

FEED = {"methane": 0.40, "ethane": 0.30, "propane": 0.30}


def weinaug_katz(Ps, Vml, Vmg, xs, ys):
    rhoL = 1.0 / Vml
    rhoV = 1.0 / Vmg
    tot = sum(P * (x * rhoL - y * rhoV) for P, x, y in zip(Ps, xs, ys))
    return tot ** 4


fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    return proc


# ---- 1. kernel accepts fromTable and reports provenance ----
out = run(case)
check("fromTable case runs", out.returncode == 0, out.stderr[:200])
check("fromTable provenance line present",
      "parachors fromTable (Weinaug-Katz 1943)" in out.stdout)

# ---- 2. re-derive sigma from the table parachors and reported state ----
table_si = {k: v * K_CGS_TO_SI for k, v in WK_CGS.items()}
Ps = [table_si[k] for k in FEED]

# split stdout into per-flash segments: each begins at "pvt flash T ... p ..."
starts = [m.start() for m in re.finditer(r"pvt flash T ", out.stdout)]
starts.append(len(out.stdout))
segs = [out.stdout[starts[i]:starts[i + 1]] for i in range(len(starts) - 1)]
check("flash segments parsed", len(segs) >= 3, f"got {len(segs)}")

checked_two_phase = 0
for seg in segs[:3]:
    pm = re.match(r"pvt flash T \S+ p (\S+) beta (\S+) single (\S+)", seg)
    check("flash header parsed", pm is not None, seg[:60])
    if not pm:
        continue
    p, beta, single = float(pm.group(1)), float(pm.group(2)), int(pm.group(3))
    rhoL = float(re.search(r"pvt liquid rho (\S+)", seg).group(1))
    rhoV = float(re.search(r"pvt vapour rho (\S+)", seg).group(1))
    ML = float(re.search(r"pvt liquid rho \S+ M (\S+)", seg).group(1))
    MV = float(re.search(r"pvt vapour rho \S+ M (\S+)", seg).group(1))
    xs = [float(m[1]) for m in re.findall(r"pvt x (\S+) (\S+)", seg)]
    ys = [float(m[1]) for m in re.findall(r"pvt y (\S+) (\S+)", seg)]
    st = float(re.search(r"pvt surfaceTension (\S+)", seg).group(1))
    if single:
        check(f"single-phase sigma 0 at p={p:.4g}", st == 0.0, f"got {st}")
        continue
    checked_two_phase += 1
    Vml = ML / rhoL
    Vmg = MV / rhoV
    close(f"two-phase sigma from table parachors at p={p:.4g}",
          weinaug_katz(Ps, Vml, Vmg, xs, ys), st, 1e-9)

check("at least one two-phase point verified", checked_two_phase >= 1)

# ---- 3. gated twin: explicit list == fromTable, byte-identical modulo
#         provenance line ----
src = open(case).read()
si_list = " ".join(repr(table_si[k]) for k in FEED)
twin = re.sub(r"parachors fromTable;", f"parachors ( {si_list} );", src)
check("twin differs", twin != src)

with tempfile.TemporaryDirectory(prefix="cf_pvt_table_") as td:
    twin_path = os.path.join(td, "twin.dict")
    open(twin_path, "w").write(twin)
    tout = run(twin_path)
    check("twin runs", tout.returncode == 0, tout.stderr[:200])
    strip = lambda t: "\n".join(l for l in t.splitlines()
                                if "fromTable" not in l)
    check("twin output byte-identical modulo provenance line",
          strip(out.stdout) == strip(tout.stdout))

# ---- 4. refusals ----
with tempfile.TemporaryDirectory(prefix="cf_pvt_table_") as td:
    # 4a. fromTable + water (not in the published table) must fail loudly
    water = src.replace("propane   0.30;", "water     0.30;")
    r1 = run_td = run_path = os.path.join(td, "water.dict")
    open(run_path, "w").write(water)
    r1 = run(run_path)
    check("fromTable with table-absent component fails loudly",
          r1.returncode != 0 and "no Weinaug-Katz" in r1.stderr,
          f"rc={r1.returncode}")

    # 4b. parachors as a bare word (not fromTable, not a list) is rejected
    bad = src.replace("parachors fromTable;", "parachors bogus;")
    bad_path = os.path.join(td, "bad.dict")
    open(bad_path, "w").write(bad)
    r2 = run(bad_path)
    check("non-list non-fromTable parachors rejected",
          r2.returncode != 0 and "must be a list or 'fromTable'" in r2.stderr,
          f"rc={r2.returncode}")

    # 4c. without binaryInteraction the flash path never consults any built-in
    #     kij table: output must equal a case that declares kij 0 for all pairs
    kij0 = src.replace("parachors fromTable;",
                       "parachors ( " + si_list + " );") + """
binaryInteraction
{
    methane { ethane 0.0; propane 0.0; }
    ethane  { propane 0.0; }
}
"""
    kij_path = os.path.join(td, "kij0.dict")
    open(kij_path, "w").write(kij0)
    kout = run(kij_path)
    check("kij-0 twin runs", kout.returncode == 0, kout.stderr[:200])
    strip2 = lambda t: "\n".join(l for l in t.splitlines()
                                 if "fromTable" not in l)
    check("no built-in kij is consulted (explicit kij=0 identical)",
          strip2(kout.stdout) == strip2(tout.stdout))

if fails:
    print("FAIL")
    for f in fails:
        print(" -", f)
    sys.exit(1)
print("PASS: built-in parachor table (Weinaug-Katz 1943) verified")
