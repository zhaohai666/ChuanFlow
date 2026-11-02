#!/usr/bin/env python3
"""Analytic verification of the 1D steady network solver (chflow-1d-steady).

The shipped case is a laminar two-well gathering network, so the exact
Hagen-Poiseuille series/parallel solution (including the rho*g*dz gravity
heads) is known in closed form. This script runs the solver and compares
node pressures, pipe flows and the mass balance against that solution.

usage: check_steady_network.py <path-to-chflow-1d-steady> <path-to-case.dict>
"""
import math
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")

nodes, pipes = {}, {}
for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "node":
        nodes[tok[1]] = float(tok[3])
    elif tok[0] == "pipe":
        d = {}
        it = iter(tok[2:])
        for a, b in zip(it, it):
            d[a] = float(b)
        pipes[tok[1]] = d

# ---- case geometry (keep in sync with tests/cases/steadyNetwork/case.dict) ----
MU, RHO, D, G = 0.1, 850.0, 0.05, 9.80665
AREA = 0.25 * math.pi * D ** 2
P_OUT = 1.0e6
PIPES = {
    "P1": dict(frm="W1", to="M", L=1000.0, Q=0.002, dz=5.0),
    "P2": dict(frm="W2", to="M", L=1500.0, Q=0.001, dz=10.0),
    "P3": dict(frm="M", to="Outlet", L=800.0, Q=0.003, dz=0.0),
}


def resistance(L):
    return 128.0 * MU * L / (math.pi * D ** 4)


# exact node pressures by back-substitution from the terminal
exp_p = {"Outlet": P_OUT}
exp_p["M"] = P_OUT + resistance(800.0) * 0.003
exp_p["W1"] = exp_p["M"] + resistance(1000.0) * 0.002 + RHO * G * 5.0
exp_p["W2"] = exp_p["M"] + resistance(1500.0) * 0.001 + RHO * G * 10.0

fails = []


def check(name, got, want, rtol):
    if want == 0:
        ok = abs(got) < rtol
    else:
        ok = abs(got - want) / abs(want) < rtol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got={got:.6e} want={want:.6e}")
    if not ok:
        fails.append(name)


# ---- node pressures ----
for n, want in exp_p.items():
    if n not in nodes:
        fails.append(f"missing node {n}")
        continue
    check(f"p({n})", nodes[n], want, 1e-6)

# ---- pipe flows (mass balance) and friction state ----
for pname, spec in PIPES.items():
    if pname not in pipes:
        fails.append(f"missing pipe {pname}")
        continue
    f = pipes[pname]
    check(f"Q({pname})", f["Q"], spec["Q"], 1e-6)
    v_want = spec["Q"] / AREA
    check(f"v({pname})", f["v"], v_want, 1e-6)
    re_want = RHO * v_want * D / MU
    check(f"Re({pname})", f["Re"], re_want, 1e-4)
    if f["Re"] > 2300:
        fails.append(f"{pname} left the laminar regime")
    dp_want = (resistance(spec["L"]) * spec["Q"] + RHO * G * spec["dz"])
    check(f"dp({pname})", f["dp"], dp_want, 1e-6)

# ---- global mass balance: what enters must reach the terminal ----
check("Q(P1)+Q(P2)==Q(P3)", pipes["P1"]["Q"] + pipes["P2"]["Q"], pipes["P3"]["Q"], 1e-9)
check("injection sum", 0.003, pipes["P3"]["Q"], 1e-9)

if fails:
    raise SystemExit(f"steady network verification FAILED: {fails}")
print("steady network case matches the analytic Hagen-Poiseuille solution")
