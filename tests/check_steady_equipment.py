#!/usr/bin/env python3
"""Closed-form verification of the steady-network equipment units
(chflow-1d-steady): inline valve (loss coefficient K) and pump (quadratic
head curve) on laminar single-phase pipes, where the pipe law reduces to an
explicit quadratic and the checker never iterates.

usage: check_steady_equipment.py <exe> <case.dict> <mode>
  mode = valve    : both ends fixed, horizontal; q is the positive root of
            128 mu L q/(pi D^4) + K rho q^2/(2 A^2) = dp; dpValve is the
            quadratic term at the same root.
  mode = pump     : a closed loop (both end pressures equal); the operating
            point is where the pump curve rho g H(q) crosses the laminar
            system curve, and the pump lift equals the friction drop exactly.
  mode = network  : rate-fed valve pipe into a free junction, then a pump
            lifting over a hill to a fixed terminal; the free pressures
            follow the closed-form chain and the pump nearly balances the
            lift (a wrong pump-head sign would flip p_J1 by ~4 bar).
"""
import math
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

exe, case, mode = sys.argv[1], sys.argv[2], sys.argv[3]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")

nodes, pipes = {}, {}
for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "node" and len(tok) >= 4:
        nodes[tok[1]] = float(tok[3])
    elif tok[0] == "pipe":
        d = {}
        for j in range(2, len(tok) - 1, 2):
            d[tok[j]] = float(tok[j + 1])
        pipes[tok[1]] = d
if not pipes:
    raise SystemExit("no pipe lines parsed")

fails = []


def check(ok, msg):
    if not ok:
        fails.append(msg)


RHO, MU, G = 850.0, 0.1, 9.80665


def laminar_res(L, D):
    """d(q) = R q for a laminar pipe."""
    return 128.0 * MU * L / (math.pi * D ** 4)


def quadratic_root(a, b, c):
    """positive root of a q^2 + b q + c = 0 (a > 0)"""
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


if mode == "valve":
    K, DP = 500.0, 2.0e5
    L, D = 1000.0, 0.05
    A = math.pi * D * D / 4.0
    R = laminar_res(L, D)
    Kv = K * RHO / (2.0 * A * A)
    q = quadratic_root(Kv, R, -DP)
    p1 = nodes["W1"]
    check(abs(p1 - (1.0e5 + DP)) < 1e-6, f"W1 pressure {p1} != outlet + dp")
    qv = pipes["P1"]["Q"]
    check(abs(qv - q) < 1e-9 * q, f"rate {qv} != quadratic root {q}")
    dpv = pipes["P1"]["dpValve"]
    check(abs(dpv - Kv * q * q) < 1e-6, f"dpValve {dpv} != K rho q^2/2A^2 "
                                        f"{Kv * q * q}")
    dpv_share = dpv / DP
    check(0.02 < dpv_share < 0.04, f"valve share {dpv_share} not in the "
                                   f"visible band")
    print(f"valve OK: q {qv:.6e} m3/s (root {q:.6e}), dpValve {dpv:.2f} Pa "
          f"({100 * dpv_share:.1f}% of the drop), Re {pipes['P1']['Re']:.1f}")

elif mode == "pump":
    H0, H1, Q1 = 40.0, 25.0, 3.0e-4
    L, D = 1000.0, 0.05
    R = laminar_res(L, D)
    # rho g (h0 - (h0 - h1)(q/q1)^2) = R q  ->  quadratic
    a = RHO * G * (H0 - H1) / (Q1 * Q1)
    q = quadratic_root(a, R, -RHO * G * H0)
    qp = pipes["P1"]["Q"]
    check(abs(qp - q) < 1e-9 * q, f"rate {qp} != operating point {q}")
    H = H0 - (H0 - H1) * (q / Q1) ** 2
    dpp = pipes["P1"]["dpPump"]
    check(abs(dpp + RHO * G * H) < 1e-6, f"dpPump {dpp} != -rho g H(q) "
                                         f"{-RHO * G * H}")
    # closed loop: the pump lift is consumed by friction exactly
    fr = R * q
    check(abs(fr + dpp) < 1e-6 * abs(dpp),
          f"pump lift {dpp} != friction drop {fr}")
    check(abs(pipes["P1"]["dp"]) < 1e-6, "closed loop reports a net dp")
    check(60 < pipes["P1"]["Re"] < 80, "operating point left the laminar"
                                       " regime the checker assumes")
    print(f"pump OK: q {qp:.6e} m3/s, H {H:.4f} m, dpPump {dpp:.2f} Pa = "
          f"friction {-fr:.2f} Pa")

elif mode == "network":
    Q, H0, H1, Q1, KV, LIFT = 3.0e-4, 40.0, 25.0, 3.0e-4, 300.0, 10.0
    A = math.pi * 0.05 * 0.05 / 4.0
    R2 = laminar_res(500.0, 0.05)
    H = H0 - (H0 - H1) * (Q / Q1) ** 2
    fr1 = R2 * Q
    fr2 = R2 * Q
    dpv = KV * RHO * Q * Q / (2.0 * A * A)
    p_out = 1.0e5
    p_j1 = p_out + RHO * G * LIFT + fr2 - RHO * G * H
    p_w1 = p_j1 + fr1 + dpv
    check(abs(nodes["J1"] - p_j1) < 1e-6,
          f"J1 {nodes['J1']} != closed form {p_j1}")
    check(abs(nodes["W1"] - p_w1) < 1e-6,
          f"W1 {nodes['W1']} != closed form {p_w1}")
    qp1, qp2 = pipes["P1"]["Q"], pipes["P2"]["Q"]
    check(abs(qp1 - Q) < 1e-12 and abs(qp2 - Q) < 1e-12,
          f"rates {qp1}/{qp2} != the imposed {Q}")
    dpp = pipes["P2"]["dpPump"]
    check(abs(dpp + RHO * G * H) < 1e-6, f"dpPump {dpp} != -rho g H(q)")
    check(p_j1 < p_out, "junction sits above the terminal: the pump sign "
                        "convention looks flipped")
    print(f"network OK: p_W1 {nodes['W1']:.2f} Pa, p_J1 {nodes['J1']:.2f} Pa "
          f"(terminal {p_out:.0f}), dpPump {dpp:.2f} Pa vs lift "
          f"{RHO * G * LIFT:.0f} Pa")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"{mode}: {len(fails)} checks FAILED")
print(f"{mode}: all checks passed")
