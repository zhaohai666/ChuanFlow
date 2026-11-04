#!/usr/bin/env python3
"""Analytic verification of the two-fluid temperature equation
(chflow-1d-twofluid).

usage: check_two_fluid_thermal.py <exe> <case.dict> <mode>
  mode = relax      : sealed pipe, constant wall. Each section is an
            independent lumped capacity C_i dT/dt = k' dx (T_wall - T_i);
            implicit Euler has the exact discrete answer
            T^n = T_wall + (T_0 - T_wall) / (1 + dt k' dx / C_i)^n, applied
            per cell with the phase-weighted capacity (liquid-filled vs gas
            cells relax at visibly different rates).
  mode = exchanger  : steady single-phase flow, constant wall. The discrete
            steady state is the recurrence T_i - T_wall = (T_{i-1} -
            T_wall) * r with r = mcd / (mcd + k' dx), and the wall heat the
            sections absorb equals the enthalpy the stream loses,
            sum k' dx (T_i - T_wall) = mcd (T_in - T_out).
  mode = geothermal : sealed vertical well relaxing onto the formation
            profile T_wall(d) = surf + grad*d; the steady state is exact,
            T_i = T_wall(d_i) at the cell-centre depths, and the final wall
            heat rate vanishes.
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

# ---- parse the time series ----
times, frames = [], []
cur_t, cur_vm, cur_dp, cur_cells, cur_thermal = None, None, None, [], None


def flush():
    global cur_t, cur_vm, cur_dp, cur_cells, cur_thermal
    if cur_t is not None:
        times.append(cur_t)
        frames.append({"t": cur_t, "vm": cur_vm, "dp": cur_dp,
                       "cells": cur_cells, "thermal": cur_thermal})
        cur_t, cur_vm, cur_dp, cur_cells, cur_thermal = None, None, None, [], None


for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "t":
        flush()
        cur_t, cur_vm, cur_dp = float(tok[1]), float(tok[3]), float(tok[5])
    elif tok[0] == "cell":
        cur_cells.append({"alpha": float(tok[3]), "water": float(tok[5]),
                          "vsl": float(tok[7]), "vsg": float(tok[9]),
                          "T": float(tok[11])})
    elif tok[0] == "thermal:":
        cur_thermal = (float(tok[2]), float(tok[4]))
flush()
if not times:
    raise SystemExit("no time series parsed")

fails = []


def check(ok, msg):
    if not ok:
        fails.append(msg)


# ---- case parameters (mirrored from the dict) ----
RHO_L, RHO_G, C_L, C_G = 1000.0, 20.0, 4180.0, 2150.0
LENGTH, DIAMETER, SECTIONS = 10.0, 0.02, 20
K_WALL, DT, T0, T_WALL = 2.0, 2.0, 320.0, 280.0
N_STEPS = 30
if mode == "exchanger":
    LENGTH, DIAMETER, SECTIONS = 2.0, 0.02, 100
    K_WALL, DT, T0, T_WALL = 20.0, 0.8, 300.0, 280.0
    V_M, T_IN = 0.02, 330.0
elif mode == "geothermal":
    LENGTH, DIAMETER, SECTIONS = 200.0, 0.05, 100
    K_WALL, DT, T0 = 500.0, 10000.0, 320.0
    SURF, GRAD, DEPTH_IN = 288.0, 0.03, 300.0

DX = LENGTH / SECTIONS
AREA = math.pi * DIAMETER * DIAMETER / 4.0
KW = K_WALL * DX

last = frames[-1]
Ts = [c["T"] for c in last["cells"]]
if len(Ts) != SECTIONS:
    raise SystemExit(f"expected {SECTIONS} cell temperatures, got {len(Ts)}")

if mode == "relax":
    # liquid column occupies cells whose centre (i+0.5) dx <= 5
    for i, T in enumerate(Ts):
        alpha = 1.0 if (i + 0.5) * DX <= 5.0 else 0.0
        cap = AREA * DX * (alpha * RHO_L * C_L + (1 - alpha) * RHO_G * C_G)
        lam = 1.0 / (1.0 + DT * KW / cap)
        exact = T_WALL + (T0 - T_WALL) * lam ** N_STEPS
        check(abs(T - exact) < 1e-9,
              f"cell {i} (alpha={alpha}) T {T} != discrete closed form "
              f"{exact}")
    # the two capacity groups must be distinct: the gas cells sit on the wall
    check(abs(Ts[SECTIONS - 1] - T_WALL) < 0.1,
          f"gas cells did not reach the wall temperature: {Ts[-1]}")
    check(Ts[0] > T_WALL + 30.0,
          f"liquid cells barely relaxed: {Ts[0]}")
    print(f"relax OK: liquid T {Ts[0]:.6f}, gas T {Ts[-1]:.6f} "
          f"(wall {T_WALL}), all cells on the discrete closed form")

elif mode == "exchanger":
    mcd = AREA * V_M * RHO_L * C_L
    r = mcd / (mcd + KW)
    for i, T in enumerate(Ts):
        exact = T_WALL + (T_IN - T_WALL) * r ** (i + 1)
        check(abs(T - exact) < 1e-9,
              f"cell {i} T {T} != recurrence {exact}")
    q_wall = sum(KW * (T - T_WALL) for T in Ts)
    q_stream = mcd * (T_IN - Ts[-1])
    check(abs(q_wall - q_stream) < 1e-9 * max(abs(q_wall), 1.0),
          f"steady identity broken: wall {q_wall} vs stream {q_stream}")
    # the continuous exponential differs at O(dx): report it, assert nothing
    print(f"exchanger OK: T_out {Ts[-1]:.6f} K (recurrence "
          f"{T_WALL + (T_IN - T_WALL) * r ** SECTIONS:.6f}), "
          f"continuous Ramey {T_WALL + (T_IN - T_WALL) * math.exp(-K_WALL * LENGTH / mcd):.4f} K, "
          f"wall heat {q_wall:.4f} W = stream loss {q_stream:.4f} W")

elif mode == "geothermal":
    for i, T in enumerate(Ts):
        depth = DEPTH_IN - (i + 0.5) * DX  # angleDeg 90: s along == TVD drop
        exact = SURF + GRAD * depth
        check(abs(T - exact) < 1e-9,
              f"cell {i} T {T} != geothermal {exact} at depth {depth}")
    check(Ts[0] > Ts[-1], "profile not monotone down the depth axis")
    qw = last["thermal"][0] if last["thermal"] else None
    check(qw is not None and abs(qw) < 1e-9,
          f"final wall heat rate {qw} did not vanish")
    print(f"geothermal OK: T_top {Ts[-1]:.6f} K, T_bottom {Ts[0]:.6f} K, "
          f"qWall {qw:.3e} W")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"{mode}: {len(fails)} checks FAILED")
print(f"{mode}: all checks passed")
