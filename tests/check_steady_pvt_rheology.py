#!/usr/bin/env python3
"""5.76 / A5 increment 2: a declared wall law on the PVT march.

The PVT table still supplies rho(p) and the Newtonian mu(p); a declared
flowIndex/consistency/yieldStress reprices that baseline on the segment's
wall rate gdot = 8|v|/D with the Rabinowitsch K' and Papanastasiou's
regularised yield share, BEFORE it enters the same Churchill call. Zero
declaration returns the table value bit for bit.

The shipped case uses a FLAT table (rho 1000, mu 0.01 everywhere) and two
fixed-pressure nodes, so there is no compressibility drift and no network
solve: v = Q/A is one number in all 24 segments and the laminar wall stress

    tau_w = dp D / (4 L) = 0.25 Pa

is purely geometric -- the SAME number for every pipe, Newtonian or not.
Each pipe is therefore an independent closed-form statement of its own wall
law inverted for the rate:

    tau_w = K' gdot^n + tau_0 (1 - exp(-100 gdot)),
    K' = K ((3n+1)/(4n))^n,  n = 0 -> mu_table gdot,
    gdot = 8 v / D,  Q = A v,  Re = rho v D / (tau_w/gdot).

The checker re-parses the table, the geometry and each pipe's declaration
from the case dictionary, inverts that law, re-marches the 24 segments with
the full Churchill blend, and asserts:

  * convergence, no choke, iterations = 0 (no free node),
  * reported Q / v / Re / f / dp against the closed form,
  * f Re == 64 and Re < 2040 (the laminar member is exact, not a fit),
  * the undeclared pipe against Hagen-Poiseuille pi d^4 dp / (128 mu L),
  * flowIndex 1 with K = mu is bit-identical to the undeclared pipe,
  * K' is the Rabinowitsch consistency, not the textbook K (a 1.25x rate),
  * the five laws are visibly distinct at that one wall stress,
  * four parse-time refusals, loud and naming their reason.

usage: check_steady_pvt_rheology.py <exe> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]
case = os.path.abspath(case)
text = re.sub(r"//[^\n]*", "", open(case).read())

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol, atol=0.0):
    return abs(got - want) <= rtol * abs(want) + atol


def block(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start + 1:j]
    raise SystemExit(f"unterminated {key}{{}} block")


def entries(src):
    out = {}
    for m in re.finditer(r"(\w+)\s*\{", src):
        start = m.end() - 1
        depth = 0
        for j in range(start, len(src)):
            if src[j] == "{":
                depth += 1
            elif src[j] == "}":
                depth -= 1
                if depth == 0:
                    out[m.group(1)] = src[start + 1:j]
                    break
    return out


def num(src, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", src)
    return float(m.group(1)) if m else default


def table(key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", block(text, "pvt"))
    return [float(x) for x in m.group(1).split()]


P_GRID, RHO_TAB, MU_TAB = table("pressure"), table("rho"), table("mu")
check("flat rho/mu table (no compressibility drift in this case)",
      len(P_GRID) == len(RHO_TAB) == len(MU_TAB) and len(P_GRID) >= 2
      and max(RHO_TAB) == min(RHO_TAB) and max(MU_TAB) == min(MU_TAB))
RHO, MU = RHO_TAB[0], MU_TAB[0]

NODE = {k: num(v, "pressure") for k, v in entries(block(text, "nodes")).items()}
PIPES = {}
for name, body in entries(block(text, "pipes")).items():
    PIPES[name] = {
        "from": re.search(r"\bfrom\s+(\w+)\s*;", body).group(1),
        "to": re.search(r"\bto\s+(\w+)\s*;", body).group(1),
        "L": num(body, "length"), "D": num(body, "diameter"),
        "eps": num(body, "roughness"), "z": num(body, "elevation"),
        "n": num(body, "flowIndex"), "K": num(body, "consistency"),
        "tau0": num(body, "yieldStress"),
    }


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n"
                         f"{proc.stdout}\n{proc.stderr}")
    nodes, pipes, raw, converged, iters = {}, {}, {}, False, None
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node":
            nodes[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            d, rest = {}, tok[2:]
            for a, b in zip(rest[0::2], rest[1::2]):
                d[a] = b
            for key in ("Q", "v", "Re", "f", "dp"):
                if key in d:
                    d[key] = float(d[key])
            pipes[tok[1]] = d
            raw[tok[1]] = line.strip()
        elif tok[0] == "steady":
            converged = "(converged)" in proc.stdout
            m = re.search(r"iterations=(\d+)", proc.stdout)
            iters = int(m.group(1)) if m else None
    return proc.stdout, nodes, pipes, raw, converged, iters


stdout, nodes, pipes, raw, converged, iters = run(case)
check("steady solve converged with no free node (iterations = 0)",
      converged and iters == 0, f"iterations={iters}")
check("no pipe choked", all("CHOKED" not in line for line in raw.values()))
check("the two fixed pressures are the case's", near(nodes["A"], NODE["A"], 1e-14)
      and near(nodes["B"], NODE["B"], 1e-14))

REG = 100.0  # Papanastasiou regularisation in wallRateApparentViscosityPaS
SHEAR_FLOOR = 1e-6


def wall_stress(p, gdot):
    """tau_w of the declared law at gdot; n == 0 uses the table's mu."""
    g = max(gdot, SHEAR_FLOOR)
    n, K, t0 = p["n"], p["K"], p["tau0"]
    if n > 0:
        kp = K * ((3.0 * n + 1.0) / (4.0 * n)) ** n
        tau = kp * g ** n
    else:
        tau = MU * g
    if t0 > 0:
        tau += t0 * (-math.expm1(-REG * g))
    return tau


def invert(p, tau_w):
    lo, hi = 0.0, 1.0
    while wall_stress(p, hi) < tau_w:
        hi *= 2.0
        if hi > 1e12:
            raise SystemExit("wall law has no finite inverse")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if wall_stress(p, mid) < tau_w:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def churchill(re, eps):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * eps
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


for name, p in PIPES.items():
    if name not in pipes:
        raise SystemExit(f"missing pipe {name} in solver output")
    r = pipes[name]
    D, L = p["D"], p["L"]
    area = math.pi * D * D / 4.0
    dp = NODE[p["from"]] - NODE[p["to"]]
    tau_w = dp * D / (4.0 * L)
    gdot = invert(p, tau_w)
    v = gdot * D / 8.0
    q = area * v
    mu_app = tau_w / gdot
    re = RHO * v * D / mu_app
    f = churchill(re, p["eps"] / D)

    check(f"{name}: dp is the imposed drop", near(r["dp"], dp, 1e-12))
    check(f"{name}: Q = the inverted wall law", near(r["Q"], q, 1e-8),
          f"got {r['Q']:.10g} want {q:.10g}")
    check(f"{name}: v = gdot D / 8", near(r["v"], v, 1e-8),
          f"got {r['v']:.10g} want {v:.10g}")
    check(f"{name}: Re on the apparent viscosity", near(r["Re"], re, 1e-8))
    check(f"{name}: f is the full Churchill blend", near(r["f"], f, 1e-8))
    check(f"{name}: laminar member f Re == 64", near(r["f"] * r["Re"], 64.0, 1e-8))
    check(f"{name}: Re < 2040 (no turbulent member)", re < 2040.0,
          f"Re={re:.4g}")

    # restate the kernel's 24-segment Euler recursion on the reported rate
    pcur = NODE[p["from"]]
    for _ in range(24):
        pcur -= churchill(re, p["eps"] / D) * (L / 24.0 / D) * RHO * v * v / 2.0
    check(f"{name}: re-march lands on the downstream node",
          near(pcur, NODE[p["to"]], 1e-12, 1e-6))

# ---- the undeclared pipe is exactly Hagen-Poiseuille -------------------
n0 = PIPES["Newt"]
dp0 = NODE[n0["from"]] - NODE[n0["to"]]
tau_w0 = dp0 * n0["D"] / (4.0 * n0["L"])
q_hp = math.pi * n0["D"] ** 4 * dp0 / (128.0 * MU * n0["L"])
check("undeclared PVT pipe is Hagen-Poiseuille pi d^4 dp / (128 mu L)",
      near(pipes["Newt"]["Q"], q_hp, 1e-9),
      f"{pipes['Newt']['Q']:.10g} vs {q_hp:.10g}")

# ---- declared-but-Newtonian is bit-identical on the printed line -------
def payload(line):
    return " ".join(line.split()[3:])  # drop "[chflow] pipe <name>"


check("flowIndex 1 with K = mu is bit-identical to the undeclared pipe",
      payload(raw["Newt"]) == payload(raw["PowN1"]),
      f"Newt [{payload(raw['Newt'])}] PowN1 [{payload(raw['PowN1'])}]")

# ---- K' is the Rabinowitsch consistency, not the textbook K -----------
thin = PIPES["Thin05"]
n, K = thin["n"], thin["K"]
kp = K * ((3.0 * n + 1.0) / (4.0 * n)) ** n
gdot_k = (tau_w0 / kp) ** (1.0 / n)
gdot_textbook = (tau_w0 / K) ** (1.0 / n)
check("Rabinowitsch K' = K ((3n+1)/(4n))^n is what the wall carries",
      near(pipes["Thin05"]["v"], gdot_k * n0["D"] / 8.0, 1e-8))
check("the textbook K would give a 1.25x larger rate (so it is NOT used)",
      near(gdot_textbook / gdot_k, 1.25, 1e-12)
      and abs(pipes["Thin05"]["v"] - gdot_textbook * n0["D"] / 8.0)
      > 0.1 * pipes["Thin05"]["v"],
      f"K'/K = {kp / K:.6g}, gdot ratio {gdot_textbook / gdot_k:.6g}")

# ---- the five laws are distinct at one fixed wall stress --------------
qs = [pipes[k]["Q"] for k in ("Newt", "Thin05", "Yield", "HB")]
check("Newtonian / power-law / Bingham / Herschel-Bulkley are distinct",
      all(qs[i] != qs[j] for i in range(4) for j in range(i + 1, 4)),
      f"Q = {[f'{x:.6g}' for x in qs]}")
check("ordering follows the wall laws at tau_w = 0.25 Pa",
      qs[3] < qs[1] < qs[2] < qs[0])

# ---- parse-time refusals: loud, not silent defaults -------------------
probes = [
    ("flowIndex without a consistency",
     text.replace("flowIndex 0.5; consistency 0.08944271909999159;",
                  "flowIndex 0.5;", 1), "positive consistency"),
    ("consistency without a flowIndex",
     text.replace("flowIndex 0.5; consistency 0.08944271909999159;",
                  "consistency 0.08944271909999159;", 1),
     "not a quantity the wall law reads"),
    ("negative flowIndex",
     text.replace("flowIndex 0.5;", "flowIndex -0.5;", 1),
     "negative"),
    ("negative yieldStress",
     text.replace("yieldStress 0.05; }", "yieldStress -0.05; }", 1),
     "negative yield"),
]
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, "probe.dict")
    for what, mutated, needle in probes:
        if mutated == text:
            raise SystemExit(f"probe '{what}' did not mutate the case")
        open(path, "w").write(mutated)
        rr = subprocess.run([exe, path], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {what}", ok, f"needle '{needle}'")

if fails:
    raise SystemExit(f"steadyPvtRheology FAILED: {fails}")
print("steadyPvtRheology OK: flat-table PVT march reprices mu on the wall rate;"
      " 5 laws x closed form + K' + bit-identity + 4 refusals")
