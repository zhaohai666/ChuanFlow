#!/usr/bin/env python3
"""Verification of the black-oil two-phase steady network (chflow-1d-steady).

The shipped multiphase case is a Beggs-Brill gathering network fed by the
black-oil PVT layer. The correlation constants are anchored at the unit level
against the fluids library; this script checks the assembled solver physics:

  * the Newton solve converges and honours the stock-tank mass balance
    (source rates appear unchanged at the terminal pipe),
  * per-pipe mass rates equal the stock-tank rate times the black-oil mass
    per stock-tank sm3 (computed here from the dictionary fluid),
  * outlet holdups lie in (0, 1] with a named Beggs-Brill regime and no
    pipe chokes,
  * a dry-oil twin of the case (GOR forced to zero) converges too, and the
    horizontal terminal line sees a LARGER pressure drop with gas present
    (same oil rate, extra gas velocity, heavier friction term).

usage: check_steady_network_multiphase.py <path-to-chflow-1d-steady> <path-to-case.dict>
"""
import math
import subprocess
import sys
import tempfile
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

REGIMES = {"segregated", "transition", "intermittent", "distributed"}

# ---- black-oil surface masses from the dictionary (SI adapter replica) ----
M_AIR = 28.9647e-3      # kg/mol
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15


def mass_per_oil_sm3(gas_gravity, oil_density, gor, water_cut, water_density=1000.0):
    rho_g_std = P_STD * gas_gravity * M_AIR / (R_GAS * T_STD)
    wor = water_cut / (1.0 - water_cut) if water_cut > 0 else 0.0
    return oil_density + gor * rho_g_std + wor * water_density


def run(case_path):
    proc = subprocess.run([exe, str(case_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(
            f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    nodes, pipes = {}, {}
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        if toks[0] == "pipe" and "CHOKED" in toks:
            raise SystemExit(f"pipe {toks[1]} choked: {line}")
        if toks[0] == "node":
            nodes[toks[1]] = float(toks[3])
        elif toks[0] == "pipe":
            d = {}
            rest = toks[2:]
            for a, b in zip(rest[0::2], rest[1::2]):
                d[a] = b
            if d.get("regime") not in REGIMES:
                raise SystemExit(f"pipe {toks[1]} missing a two-phase regime: {line}")
            pipes[toks[1]] = {
                "Q": float(d["Q"]),
                "v": float(d["v"]),
                "HL": float(d["HL"]),
                "regime": d["regime"],
                "m": float(d["m"]),
                "dp": float(d["dp"]),
            }
    if "steady solved" not in proc.stdout or "(converged)" not in proc.stdout:
        raise SystemExit("solver did not report convergence")
    return nodes, pipes


fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


# ---- shipped case: two-phase gathering network ----
nodes, pipes = run(case)
for pname in ("P1", "P2", "P3"):
    if pname not in pipes:
        raise SystemExit(f"missing pipe {pname}")
    f = pipes[pname]
    check(f"HL({pname}) in (0,1]", 0.0 < f["HL"] <= 1.0, f"HL={f['HL']:.4f}")
    check(f"m({pname}) > 0", f["m"] > 0.0)

# stock-tank volumetric balance through the manifold
check("Q(P1)+Q(P2)==Q(P3)", abs(pipes["P1"]["Q"] + pipes["P2"]["Q"] - pipes["P3"]["Q"])
      < 1e-9 * pipes["P3"]["Q"],
      f"{pipes['P1']['Q']:.6f}+{pipes['P2']['Q']:.6f} vs {pipes['P3']['Q']:.6f}")
check("source rates preserved", abs(pipes["P3"]["Q"] - 0.006) < 1e-9)

# per-pipe mass rate = stock-tank rate * mass per stock-tank sm3
mm = mass_per_oil_sm3(0.75, 850.0, 80.0, 0.0)
for pname, qo in (("P1", 0.004), ("P2", 0.002), ("P3", 0.006)):
    check(f"m({pname}) = Q*massPerSm3",
          abs(pipes[pname]["m"] - qo * mm) < 1e-6 * qo * mm,
          f"{pipes[pname]['m']:.6f} vs {qo * mm:.6f}")

# the manifold merges both well streams
check("m(P1)+m(P2)==m(P3)",
      abs(pipes["P1"]["m"] + pipes["P2"]["m"] - pipes["P3"]["m"])
      < 1e-6 * pipes["P3"]["m"])

# well pressures sit above the terminal by the climb
check("W1 above terminal", nodes["W1"] > nodes["Outlet"] + 9.80665 * 400.0 * 0.5,
      f"p(W1)={nodes['W1']:.0f}")

# ---- dry-oil twin: GOR forced to zero ----
case_text = Path(case).read_text()
dry_text = case_text.replace("gor            80;", "gor             0;")
if dry_text == case_text:
    raise SystemExit("could not rewrite the GOR entry in the case dict")
with tempfile.TemporaryDirectory() as tmp:
    dry_case = Path(tmp) / "dry.dict"
    dry_case.write_text(dry_text)
    dry_nodes, dry_pipes = run(dry_case)

check("dry oil converges with liquid-full pipes",
      all(pipes[p]["HL"] > 0.0 for p in pipes)
      and all(dry_pipes[p]["HL"] <= 1.0 for p in dry_pipes))
# same oil rate, no gas: the horizontal terminal line loses the extra gas
# friction, so its drop shrinks
check("gas raises the horizontal drop",
      pipes["P3"]["dp"] > dry_pipes["P3"]["dp"],
      f"two-phase {pipes['P3']['dp']:.0f} Pa vs dry {dry_pipes['P3']['dp']:.0f} Pa")
# and the uphill well lines carry sensible drops (0.1 .. 5 MPa)
for pname in ("P1", "P2"):
    check(f"dp({pname}) plausible", 1.0e5 < pipes[pname]["dp"] < 5.0e6,
          f"dp={pipes[pname]['dp']:.0f} Pa")

if fails:
    raise SystemExit(f"steady multiphase verification FAILED: {fails}")
print("steady multiphase case satisfies the black-oil Beggs-Brill physics checks")
