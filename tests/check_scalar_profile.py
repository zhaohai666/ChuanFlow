#!/usr/bin/env python3
"""Run chflow-scalar-transport on a case and check the solution against
an analytic profile. Stdlib only; parses the ASCII .vtu cell data.

usage: check_scalar_profile.py <solver-binary> <case.dict> <mode> <tol>
  mode: advection  -> 1D steady convection-diffusion, both Dirichlet
        linearZ    -> pure diffusion, T = z / zmax
"""
import math
import re
import subprocess
import sys


def read_vtu(path):
    txt = open(path).read()
    pts = re.search(
        r'<DataArray type="Float64" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
        txt, re.S).group(1)
    coords = [float(v) for v in pts.split()]
    points = list(zip(coords[0::3], coords[1::3], coords[2::3]))

    conn = re.search(r'Name="connectivity" format="ascii">\s*(.*?)</DataArray>', txt, re.S).group(1)
    offs = re.search(r'Name="offsets" format="ascii">\s*(.*?)</DataArray>', txt, re.S).group(1)
    values = None
    for m in re.finditer(
        r'<DataArray type="Float64" Name="([^"]+)" format="ascii">\s*(.*?)</DataArray>', txt, re.S):
        if m.group(1) == "T":
            values = [float(v) for v in m.group(2).split()]
    if values is None:
        raise SystemExit("field T not found in " + path)

    ci = [int(v) for v in conn.split()]
    oi = [int(v) for v in offs.split()]
    centers = []
    prev = 0
    for o in oi:
        vids = ci[prev:o]
        prev = o
        n = len(vids)
        cx = sum(points[v][0] for v in vids) / n
        cy = sum(points[v][1] for v in vids) / n
        cz = sum(points[v][2] for v in vids) / n
        centers.append((cx, cy, cz))
    return centers, values


def analytic_advection(x, U, L, D):
    pe = U * L / D
    if abs(pe) < 1e-12:
        return 1.0 - x / L
    e1 = math.exp(pe * x / L)
    e2 = math.exp(pe)
    return (e1 - e2) / (1.0 - e2)


def main():
    binary, case, mode, tol = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    import os
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        raise SystemExit("solver failed")

    vtu = os.path.join(workdir,
                       "out_advection.vtu" if mode == "advection" else "out_poisson.vtu")
    centers, T = read_vtu(vtu)

    if mode == "advection":
        U, L, D = 0.1, 1.0, 0.01
        err = max(abs(t - analytic_advection(c[0], U, L, D)) for c, t in zip(centers, T))
    elif mode == "linearZ":
        zmax = max(c[2] for c in centers) + 0.05  # cell centre vs domain: use analytic slope
        # domain is 0..0.1 in z with T(0)=0, T(0.1)=1: exact linear field, cell-centre values
        err = max(abs(t - c[2] / 0.1) for c, t in zip(centers, T))
    else:
        raise SystemExit("unknown mode")

    print(f"max abs error vs analytic: {err:.4e} (tol {tol:.1e})")
    if err > tol:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
