#!/usr/bin/env python3
"""Independent verification of the Griffith bubble / mist branches of the
classic vertical Hagedorn-Brown method (2026-11-25 fix).

usage: check_steady_hb_griffith.py <steady-exe> <steady-exe>

No case file is read: all three cases are generated from one black-oil
one-pipe dictionary, and every physical number is re-derived from the
kernel's printed `v` and `lambda` using the published equations, not from
the C++ holdup code. Checks:
  * bubble branch at low rate: HL = 1 - H_g with the Griffith & Wallis
    (1961) bubble-slip relation, and the L_B criterion that selects it;
  * mist branch at high gas/oil ratio: HL = lambda (Duns-Ros Region III
    no-slip) while the plain Whitson chart gives a larger holdup;
  * hbGriffith false / absent is byte-identical to the undeclared run;
  * loud refusals for hbGriffith with the wrong correlation / unknown word.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 2:
    raise SystemExit(__doc__)
STEADY = sys.argv[1]

fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol=1e-8, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


FT = 3.28083989501312


def run(body):
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "t.dict")
        with open(p, "w") as f:
            f.write(body)
        r = subprocess.run([STEADY, p], capture_output=True, text=True, cwd=tmp)
    return r


def one_pipe(decl, gor, rate):
    return f"""
fluid {{ blackoil {{ gasGravity 0.75; oilDensity 850; gor {gor}; waterCut 0;
    temperature 340; surfaceTension 0.025; }} }}
nodes {{ W {{}}; Outlet {{ pressure 4000000; }} }}
pipes {{ P {{ from W; to Outlet; length 800; diameter 0.1; roughness 4.6e-5;
    elevation 0; {decl} }} }}
sources {{ s {{ node W; rate {rate}; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


def pipe_line(out):
    lines = [l for l in out.splitlines() if l.startswith("[chflow] pipe P")]
    if not lines:
        return None
    t = lines[0].replace("[chflow]", "").split()
    return dict(zip(t[2::2], t[3::2]))


# ---- bubble branch ------------------------------------------------------
decl = "correlation hagedornBrown; hbGriffith true;"
r = run(one_pipe(decl, 80, 0.0001))
check("bubble-arm case runs", r.returncode == 0, r.stdout[-300:])
b = pipe_line(r.stdout)
check("bubble arm reports holdup", b is not None and "HL" in b)
if b:
    vm, lam, hl = float(b["v"]), float(b["lambda"]), float(b["HL"])
    vsl, vsg = lam * vm, (1.0 - lam) * vm
    vmFt, dFt, vsgFt = vm * FT, 0.1 * FT, vsg * FT
    lb = max(1.071 - 0.2218 * vmFt * vmFt / dFt, 0.13)
    gas_content = vsg / vm
    check("bubble branch is selected by the published L_B criterion",
          gas_content < lb, f"lambda_g {gas_content:.6g} < L_B {lb:.6g}")
    vs = 0.8  # ft/s, Griffith & Wallis bubble slip velocity
    rr = vmFt / vs
    disc = (1.0 + rr) ** 2 - 4.0 * vsgFt / vs
    hg = 0.5 * ((1.0 + rr) - math.sqrt(disc))
    want = min(max(1.0 - hg, lam), 1.0)
    close("bubble holdup re-derived from the Griffith slip relation", hl, want, 1e-6)
    check("bubble holdup exceeds the no-slip content (slip is visible)",
          hl > lam + 0.1, f"HL {hl:.6g} vs lambda {lam:.6g}")
    plain = pipe_line(run(one_pipe("correlation hagedornBrown;", 80, 0.0001)).stdout)
    check("Griffith bubble branch actually changes the Hagedorn-Brown holdup",
          plain is not None and abs(float(plain["HL"]) - hl) > 0.1,
          f"plain {plain['HL'] if plain else '?'} vs Griffith {hl:.6g}")

# ---- mist branch --------------------------------------------------------
r2 = run(one_pipe(decl, 4000, 0.003))
check("mist-arm case runs", r2.returncode == 0, r2.stdout[-300:])
m = pipe_line(r2.stdout)
check("mist arm reports holdup", m is not None and "HL" in m)
if m:
    close("mist branch uses the Duns-Ros Region III no-slip holdup",
          float(m["HL"]), float(m["lambda"]), 1e-9)
    plain2 = pipe_line(run(one_pipe("correlation hagedornBrown;", 4000, 0.003)).stdout)
    check("Griffith mist branch actually changes the Whitson chart holdup",
          plain2 is not None and float(plain2["HL"]) > float(m["lambda"]) + 1e-3,
          f"plain {plain2['HL'] if plain2 else '?'} vs mist {m['lambda']}")

# ---- twin and refusals --------------------------------------------------
base_stdout = run(one_pipe("correlation hagedornBrown;", 80, 0.0001)).stdout
for word in ("", " hbGriffith false;"):
    twin = run(one_pipe("correlation hagedornBrown;" + word, 80, 0.0001)).stdout
    check(f"hbGriffith{word or ' absent'} is byte-identical to the chart",
          twin == base_stdout)

for label, body, want in (
    ("hbGriffith on beggsBrill", one_pipe("correlation beggsBrill; hbGriffith true;", 80, 0.0001), "hbGriffith needs"),
    ("hbGriffith with no correlation", one_pipe("hbGriffith true;", 80, 0.0001), "hbGriffith needs"),
    ("unknown hbGriffith word", one_pipe("correlation hagedornBrown; hbGriffith maybe;", 80, 0.0001), "hbGriffith must be"),
):
    rr = run(body)
    check(f"refuse {label}", rr.returncode != 0 and want in rr.stdout + rr.stderr,
          f"rc={rr.returncode} {(rr.stdout + rr.stderr)[-160:]}")

if fails:
    print(f"FAIL ({len(fails)}):")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("all Hagedorn-Brown Griffith-branch checks passed")
