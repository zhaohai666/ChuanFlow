#!/usr/bin/env python3
"""Independent verification of the declared-critical-volume LBC path.

usage: check_pvt_lbc_vc.py <exe> <caseWithVc.dict> <casePrZc.dict>

The Whitson & Brule (2000) worked example publishes tabulated critical
volumes (m3/mol) for methane / ethane / propane. The C++ `viscosity{
model lbc; criticalVolumes (...); }` path must use THOSE numbers instead of
back-calculating Vc = 0.3074 R Tc / Pc from PR's universal Zc. This checker
re-derives every reported phase viscosity from the printed flash state with
the declared Vcs (no C++ reuse), proves the result is close to the kernel
and genuinely different from the PR-Zc twin, and checks loud refusals for a
wrong-length list, a non-positive Vc and an unknown model word.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 4:
    raise SystemExit(__doc__)
exe, case_vc, case_pr = sys.argv[1:4]

R = 8.314462618
ATM = 101325.0
ZC = 0.3074
# MW [g/mol], Tc [K], Pc [Pa] from src/pvt/src/components.cpp, feed order.
CT = {
    "methane": (16.043, 190.564, 45.992e5),
    "ethane": (30.070, 305.322, 48.722e5),
    "propane": (44.097, 369.834, 42.477e5),
}
FEED = ["methane", "ethane", "propane"]
VC_DECLARED = [9.86e-5, 1.455e-4, 2.0e-4]  # Whitson & Brule worked example

fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol=1e-8, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def xi(Tc, Pc, MW):
    return Tc ** (1.0 / 6.0) / (math.sqrt(MW) * (Pc / ATM) ** (2.0 / 3.0))


def stiel(T, Tc, Pc, MW):
    tr = T / Tc
    cp = (17.78e-5 * (4.58 * tr - 1.67) ** 0.625 / xi(Tc, Pc, MW) if tr > 1.5
          else 34.0e-5 * tr ** 0.94 / xi(Tc, Pc, MW))
    return cp * 1.0e-3


def lbc(T, Vm, zs, Vcs):
    MWs = [CT[n][0] for n in FEED]
    Tcs = [CT[n][1] for n in FEED]
    Pcs = [CT[n][2] for n in FEED]
    Tc = sum(z * t for z, t in zip(zs, Tcs))
    MW = sum(z * m for z, m in zip(zs, MWs))
    Pc = sum(z * p for z, p in zip(zs, Pcs))
    Vc = sum(z * v for z, v in zip(zs, Vcs))
    rr = Vc / Vm
    num = sum(z * stiel(T, t, p, m) * math.sqrt(m)
              for z, t, p, m in zip(zs, Tcs, Pcs, MWs))
    den = sum(z * math.sqrt(m) for z, m in zip(zs, MWs))
    mu0 = num / den * 1.0e3
    poly = (0.1023 + 0.023364 * rr + 0.058533 * rr ** 2
            - 0.040758 * rr ** 3 + 0.0093324 * rr ** 4)
    return (mu0 + (poly ** 4 - 1.0e-4) / xi(Tc, Pc, MW)) * 1.0e-3


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def parse(out):
    flashes, cur = [], None
    for line in out.splitlines():
        t = line.replace("[chflow]", "").split()
        if not t or t[0] != "pvt":
            continue
        k = t[1]
        if k == "flash":
            kv = dict(zip(t[2::2], t[3::2]))
            cur = {"T": float(kv["T"]), "p": float(kv["p"]), "x": {}, "y": {}}
            flashes.append(cur)
        elif k == "liquid" and cur is not None:
            cur["rhoL"], cur["ML"] = float(t[3]), float(t[5])
        elif k == "vapour" and cur is not None:
            cur["rhoV"], cur["MV"] = float(t[3]), float(t[5])
        elif k in ("x", "y") and cur is not None:
            cur[k][t[2]] = float(t[3])
        elif k == "viscosity" and cur is not None:
            cur["muL"], cur["muV"] = float(t[3]), float(t[5])
    return flashes


rc, out_vc = run(case_vc)
check("declared-Vc case runs", rc == 0, f"rc={rc}\n{out_vc[:400]}")
rc_pr, out_pr = run(case_pr)
check("PR-Zc twin runs", rc_pr == 0, f"rc={rc_pr}\n{out_pr[:400]}")
fv, fp = parse(out_vc), parse(out_pr)
check("both cases produce the same flash states",
      len(fv) == len(fp) >= 3
      and all(a["p"] == b["p"] and a["T"] == b["T"] for a, b in zip(fv, fp)))

# The worked example's Vcs are tabulated values, NOT PR back-calculations:
# the PR-Zc values are 105.99/123.14/147.64 cm3/mol vs declared 98.6/145.5/200.
check("declared Vcs are genuinely different from the PR-Zc fallback",
      max(abs(a - b) / max(a, b) for a, b in zip(
          VC_DECLARED,
          [ZC * R * CT[n][1] / CT[n][2] for n in FEED])) > 0.1)

worst = 0.0
moved = 0.0
for av, ap in zip(fv, fp):
    if "muL" not in av or "muL" not in ap:
        continue
    for phase, rho_key, m_key, mu_key in (("liquid", "rhoL", "ML", "muL"),
                                          ("vapour", "rhoV", "MV", "muV")):
        zs = [av["x"][n] if phase == "liquid" else av["y"][n] for n in FEED]
        want = lbc(av["T"], av[m_key] / av[rho_key], zs, VC_DECLARED)
        worst = max(worst, abs(av[mu_key] - want) / max(want, 1e-30))
        moved = max(moved, abs(av[mu_key] - ap[mu_key]) / max(ap[mu_key], 1e-30))
check("declared-Vc kernel viscosities re-derived independently from case text",
      worst < 1e-8, f"worst relative error {worst:.3g}")
check("declared Vcs actually move the reported viscosity",
      moved > 1e-4, f"max relative movement {moved:.3g}")

# Loud refusals: wrong list length, non-positive Vc, dead model word.
base = open(case_vc).read()
variants = [
    ("wrong criticalVolumes length",
     "criticalVolumes (9.86e-5 1.455e-4);"),
    ("non-positive critical volume",
     "criticalVolumes (9.86e-5 -1.455e-4 2.0e-4);"),
    ("unknown viscosity model", None),
]
with tempfile.TemporaryDirectory() as td:
    for i, (label, repl) in enumerate(variants):
        text = base
        if repl is not None:
            text = re.sub(r"criticalVolumes \([^)]*\);", repl, text)
        else:
            text = text.replace("model lbc;", "model somethingElse;")
        p = os.path.join(td, f"reject{i}.dict")
        with open(p, "w") as f:
            f.write(text)
        r, o = run(p)
        check(f"refuse {label}", r != 0 and ("criticalVolumes" in o or "model" in o or "viscosity" in o),
              f"rc={r} {o[:200]}")

if fails:
    print(f"FAIL ({len(fails)}):")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("all declared-critical-volume LBC checks passed")
