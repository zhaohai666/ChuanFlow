#!/usr/bin/env python3
"""Run chflow-simple-foam on the channel case and verify the developed
Poiseuille profile, flux conservation and |v| ~ 0.

usage: check_poiseuille.py <solver-binary> <case.dict>
"""
import os
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
    uvals = re.search(
        r'<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
        txt, re.S).group(1)
    uv = [float(v) for v in uvals.split()]

    ci = [int(v) for v in conn.split()]
    oi = [int(v) for v in offs.split()]
    cells = []
    prev = 0
    k = 0
    for o in oi:
        vids = ci[prev:o]
        prev = o
        n = len(vids)
        cx = sum(points[v][0] for v in vids) / n
        cy = sum(points[v][1] for v in vids) / n
        cells.append((cx, cy, uv[3 * k], uv[3 * k + 1]))
        k += 1
    return cells


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    print(r.stdout[-2000:])
    if r.returncode != 0:
        print(r.stderr)
        raise SystemExit("solver failed")

    cells = read_vtu(os.path.join(workdir, "out_poiseuille.vtu"))

    H, HC, UM = 0.2, 0.1, 0.1  # channel height, centreline, mean velocity

    # vertical velocity must vanish in the developed region
    dev = [c for c in cells if c[0] > 3.0]
    if len(dev) < 10:
        raise SystemExit("not enough developed cells")
    vMax = max(abs(c[3]) for c in dev)

    def exact(y):
        return 1.5 * UM * (1.0 - ((y - HC) / (H / 2.0)) ** 2)

    profErr = max(abs(c[2] - exact(c[1])) for c in dev)
    uCentre = max(c[2] for c in dev)

    # flux conservation through one downstream column (centroids are spaced
    # by dx = 0.04, so this selects exactly the 20 cells at x = 3.5)
    col = [c for c in cells if abs(c[0] - 3.5) < 0.005]
    if len(col) != 20:
        raise SystemExit(f"expected a single cell column, got {len(col)}")
    dy = H / 20.0
    flux = sum(c[2] * dy for c in col)
    fluxErr = abs(flux / (UM * H) - 1.0)

    print(f"vMax={vMax:.3e} profileErr={profErr:.3e} centreline={uCentre:.4f} "
          f"(exact {1.5 * UM:.4f}) fluxErr={fluxErr:.3e}")

    ok = (vMax < 5e-3 and profErr < 0.012
          and abs(uCentre / (1.5 * UM) - 1.0) < 0.08 and fluxErr < 0.02)
    if not ok:
        raise SystemExit("Poiseuille verification FAILED")


if __name__ == "__main__":
    main()
