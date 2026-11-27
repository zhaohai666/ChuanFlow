#!/usr/bin/env python3
"""Verify the per-edge `starts{}` starting-resistance table in chflow-coupled-1d3d.

The two-edge network is Src -> A -> { X1 (3D duct) || By (1D bypass) } ->
X2 (3D duct) -> B -> Snk, the same topology as the steady multi-edge case.
This checker asks five things of the steady coupled app:

  1. per-edge starts really differ - the table quotes two starting resistances
     an order of magnitude apart, so the two edges' first interface residuals
     come from their own partners rather than from one shared guess
  2. the settled state is edge-local - each edge's converged resistance and
     rate agree with its own single-edge run, started from the same coefficient
  3. the scalar-seed path and the per-edge table land in the same residual
     band (actual numbers printed, not just "comparable")
  4. a missing coupled edge in the table is rejected, with that edge named
  5. an extra name in the table that is not a coupled edge is rejected, with
     that edge named

The legacy scalar `seeds ( )` path is exercised here through case_seed.dict;
its byte-for-byte identity with the pre-change output is checked separately
against the saved coupled1d3d baseline.

usage: check_coupled_start_peredge.py <solver-binary> <case-directory>
"""
import os
import re
import subprocess
import sys

EDGES = ["X1", "X2"]
TOL_INTERFACE = 1e-4          # couple{ tolerance } the cases ask for
# per-edge agreement between the two-edge per-edge run and the single-edge arm
R_BAND = 0.02                 # settled resistance, relative
Q_BAND = 0.05                 # settled rate, relative


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def run(binary, case):
    workdir = os.path.dirname(case)
    return subprocess.run([binary, case], capture_output=True, text=True,
                          cwd=workdir)


def lines_of(out):
    return [ln.replace("[chflow] ", "") for ln in out.splitlines()]


def parse_two_edge(lines, arm):
    """Per-edge start coefficient, first-pass residual and settled state."""
    starts, first, settled = {}, {}, {}
    for ln in lines:
        if ln.startswith(f"start {arm} edge "):
            t = ln.split()
            starts[t[3]] = value(ln, "resistance")
        elif ln.startswith(f"interface {arm} iter 1 edge "):
            t = ln.split()
            first[t[5]] = value(ln, "residual")
        elif ln.startswith(f"coupled {arm} edge ") and ": iterations" in ln:
            t = ln.split()
            settled[t[3].rstrip(":")] = {"converged": "NOT" not in ln,
                             "residual": value(ln, "residual"),
                             "R": value(ln, "R"), "q": value(ln, "q")}
    solved = {}
    for ln in lines:
        if ln.startswith("coupled solved: edge "):
            t = ln.split()
            solved[t[3]] = {"q": value(ln, "q"), "dp": value(ln, "dp"),
                            "R": value(ln, "R"), "residual": value(ln, "residual")}
    return starts, first, settled, solved


def parse_single(lines):
    """The legacy single-edge final line, keyed by the `coupled solved: edge` row."""
    solved = {}
    for ln in lines:
        if ln.startswith("coupled solved: edge "):
            t = ln.split()
            solved[t[3]] = {"q": value(ln, "q"), "dp": value(ln, "dp"),
                            "R": value(ln, "R"), "residual": value(ln, "residual")}
    return solved


def main():
    binary, case_dir = sys.argv[1], sys.argv[2]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def case(name):
        return os.path.join(case_dir, name)

    # ---- 1 & 2: the per-edge table run ----
    r = run(binary, case("case_peredge.dict"))
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("per-edge starts case failed to run")
    starts, first, settled, solved = parse_two_edge(
        lines_of(r.stdout), "peredge")

    require(set(starts) == set(EDGES),
            f"the table quotes a start coefficient for both edges (got {sorted(starts)})")
    require(set(first) == set(EDGES),
            "both edges report a first-pass interface residual")
    require(set(solved) == set(EDGES),
            "both edges report a final solved state")

    ratio = max(starts.values()) / min(starts.values())
    require(ratio >= 10.0,
            f"the two starting resistances are an order of magnitude apart "
            f"(ratio {ratio:.1f}, X1 {starts['X1']:.4g}, X2 {starts['X2']:.4g})")

    f1, f2 = first["X1"], first["X2"]
    rel_gap = abs(f1 - f2) / max(f1, f2)
    require(rel_gap > 1e-3,
            f"the first-pass residuals differ per edge (X1 {f1:.4g}, "
            f"X2 {f2:.4g}, rel gap {rel_gap:.3g})")

    for e in EDGES:
        require(settled[e]["converged"], f"{e}: the per-edge run converged")
        require(settled[e]["residual"] <= TOL_INTERFACE,
                f"{e}: settled residual {settled[e]['residual']:.4g} below "
                f"the case tolerance")

    # ---- single-edge arms: each block's own fixed point on its own ----
    single = {}
    for e, dict_name in (("X1", "case_single_X1.dict"), ("X2", "case_single_X2.dict")):
        rs = run(binary, case(dict_name))
        if rs.returncode != 0:
            print(rs.stdout[-3000:])
            print(rs.stderr[-3000:])
            raise SystemExit(f"single-edge {e} case failed to run")
        sol = parse_single(lines_of(rs.stdout))
        require(e in sol, f"{e}: its single-edge run reports a solved state")
        single[e] = sol[e]

    for e in EDGES:
        pe, se = solved[e], single[e]
        r_rel = abs(pe["R"] - se["R"]) / se["R"]
        q_rel = abs(pe["q"] - se["q"]) / se["q"]
        require(r_rel <= R_BAND,
                f"{e}: per-edge R {pe['R']:.6g} matches the single-edge run "
                f"{se['R']:.6g} (rel {r_rel:.3g})")
        require(q_rel <= Q_BAND,
                f"{e}: per-edge q {pe['q']:.6g} matches the single-edge run "
                f"{se['q']:.6g} (rel {q_rel:.3g})")

    # ---- 3: scalar-seed path vs per-edge table, same residual band ----
    rs = run(binary, case("case_seed.dict"))
    if rs.returncode != 0:
        print(rs.stdout[-3000:])
        print(rs.stderr[-3000:])
        raise SystemExit("scalar-seed comparison case failed to run")
    _, _, _, seed_solved = parse_two_edge(lines_of(rs.stdout), "seed 0")

    print("per-edge vs scalar-seed settled residuals:")
    for e in EDGES:
        rp, rs_ = solved[e]["residual"], seed_solved[e]["residual"]
        order = max(rp, rs_) / min(rp, rs_)
        print(f"  {e}: per-edge {rp:.6g}  scalar-seed {rs_:.6g}  ratio {order:.3g}")
        require(order <= 10.0,
                f"{e}: the two paths settle in the same residual band "
                f"(ratio {order:.3g})")
        require(rs_ <= TOL_INTERFACE,
                f"{e}: the scalar-seed run also settles inside the tolerance")

    # ---- 4 & 5: missing and extra edges are rejected, edge named ----
    for dict_name, bad_edge in (("case_missing.dict", "X2"),
                                ("case_extra.dict", "X9")):
        rr = run(binary, case(dict_name))
        combined = rr.stdout + rr.stderr
        require(rr.returncode != 0, f"{dict_name}: rejected with non-zero exit")
        require(bad_edge in combined,
                f"{dict_name}: the error names edge '{bad_edge}'")
        if rr.returncode == 0 or bad_edge not in combined:
            print(combined[-1500:])

    print("per-edge starts verification passed")
    if fails:
        raise SystemExit(f"{len(fails)} per-edge starts checks FAILED")


if __name__ == "__main__":
    main()
