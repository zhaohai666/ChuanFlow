#!/usr/bin/env python3
"""steadyAllBelowThreshold (B6): every pipe of the network sits below its own
yield threshold at once, so the steady network has a whole-network dead fixed
point -- and the solver finds it without entering Newton at all.

The law is the one 5.46 wrote down:

    dpStatic = 4 tau_0 L / d,      Q = 0 iff dp <= dpStatic
    Q = Q_Hagen-Poiseuille(dp) * g(xi),  xi = dpStatic / dp,
    g(xi) = 1 - (4/3) xi + (1/3) xi^4

with Q_HP = pi d^4 dp / (128 mu L). The kernel's free-node initial guess is the
mean of the fixed pressures (network.cpp:538-547), so with no external rate the
residual of every free node is exactly 0 = 0 at that guess: solveSteadyNetwork
evaluates once, sees relativeResidual == 0, never enters the Newton loop, and
reports iterations = 0. The HP cap on a dead branch (network.cpp:519) is what
keeps the Jacobian invertible for the arm where a source *does* force Newton
to run; it is a property of the matrix, not a zero row.

Six arms, all derived from one case file in a temporary directory; the two
checked-in dicts are the first two arms (5 pipes / 4 pipes):

    allDead     5 pipes, all thresholds above their offered dp
    topology    4 pipes, same fixed point on a different topology
    shifted     p(A) 3.0e4 -> 3.1e4: the free node follows the mean, not 25000
    knifeEdge   Stk1 tau_0 2.5: dp == dpStatic exactly, still Q = 0
    oneFlowing  Stk1 tau_0 2.0: that one pipe flows, yet iterations is still 0
    injected    a 1e-6 m3/s source at the free node: Newton runs 19 steps and
                raises E until one series pipe yields at Buckingham's exact Q

Nothing here remembers a printed number: thresholds, offered drops, the initial
guess, Buckingham's Q and the inverse dp are all re-derived from the dict and
the law. Zero kernel changes.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

BINARY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

MU = 0.1
RHO = 850.0
FIELDS = ("Q", "v", "Re", "f", "dp")

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def close(a, b, rel=1e-8):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-30)


# ---------------------------------------------------------------- dict -----

def strip_comments(src):
    return "\n".join(ln.split("//", 1)[0] for ln in src.splitlines())


def block(src, name):
    """Return the inside of a top-level `name { ... }` block, brace-counted so
    comments mentioning the word `pipes` cannot be mistaken for the block."""
    m = re.search(r"(?m)^\s*" + name + r"\s*\n?\s*\{", src)
    if not m:
        raise RuntimeError("no " + name + " block in case")
    i = src.index("{", m.start())
    depth = 0
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[i + 1:j]
    raise RuntimeError("unbalanced " + name + " block")


def entries(src, name):
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^{}]*)\}", strip_comments(block(src, name)), re.S):
        kv = {}
        for k, v in re.findall(r"(\w+)\s+([^;\s]+)\s*;", m.group(2)):
            try:
                kv[k] = float(v)
            except ValueError:
                kv[k] = v
        out[m.group(1)] = kv
    return out


def case_of(src):
    nodes = entries(src, "nodes")
    fixed = {k: v["pressure"] for k, v in nodes.items() if "pressure" in v}
    p_init = sum(fixed.values()) / len(fixed)
    pipes = entries(src, "pipes")
    for name, p in pipes.items():
        p["threshold"] = 4.0 * p.get("yieldStress", 0.0) * p["length"] / p["diameter"]
        p["offered"] = abs(fixed.get(p["from"], p_init) - fixed.get(p["to"], p_init))
    return {"nodes": nodes, "fixed": fixed, "p_init": p_init, "pipes": pipes}


# ------------------------------------------------------------- the law -----

def q_hp(dp, length, diam):
    return math.pi * diam ** 4 * dp / (128.0 * MU * length)


def g_xi(xi):
    return 1.0 - 4.0 * xi / 3.0 + xi ** 4 / 3.0


def buckingham_q(dp, length, diam, tau0):
    dps = 4.0 * tau0 * length / diam
    if dp <= dps:
        return 0.0
    return q_hp(dp, length, diam) * g_xi(dps / dp)


def invert_buckingham(q, length, diam, tau0):
    """Monotone bisection on dp; Q(dp) is strictly increasing above dpStatic."""
    dps = 4.0 * tau0 * length / diam
    lo, hi = dps, max(dps * 100.0, 1.0)
    while buckingham_q(hi, length, diam, tau0) < q:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if buckingham_q(mid, length, diam, tau0) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------- output -------

def run(case_path):
    return subprocess.run([BINARY, case_path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def parse(stdout):
    pipes, nodes, solved, net = {}, {}, {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if not tok:
            continue
        if tok[0] == "pipe" and len(tok) > 2:
            fields = {}
            i = 2
            while i + 1 < len(tok):
                key, val = tok[i], tok[i + 1]
                if key in FIELDS:
                    try:
                        fields[key] = float(val)
                    except ValueError:
                        break
                i += 2
            if set(FIELDS) <= set(fields):
                pipes[tok[1]] = fields
        elif len(tok) == 4 and tok[0] == "node" and tok[2] == "p":
            try:
                nodes[tok[1]] = float(tok[3])
            except ValueError:
                pass
        elif tok[0] == "steady" and "solved:" in tok:
            m = re.search(r"iterations=(\d+)\s+residual=([^\s]+)\s+\((\w+)", line)
            if m:
                solved = {"iterations": int(m.group(1)),
                          "residual": float(m.group(2)),
                          "state": m.group(3)}
        elif tok[0] == "network:":
            m = re.search(r"(\d+)\s+nodes,\s*(\d+)\s+pipes", line)
            if m:
                net = {"nodes": int(m.group(1)), "pipes": int(m.group(2))}
    return {"pipes": pipes, "nodes": nodes, "solved": solved, "net": net}


# --------------------------------------------------------------- arms -------

BASE = open(CASE).read()
PROBE = open(os.path.join(CASE_DIR, "probe.dict")).read()
STK1 = "Stk1 { from A; to B; length 20; diameter 0.02; yieldStress 3.0; }"
require(STK1 in BASE, "the case declares Stk1 with tau_0 = 3.0, the arm to edit")

SOURCES = ("sources\n{\n"
           "    feed { node E; rate 1.0e-6; }\n"
           "    take { node A; rate -1.0e-6; }\n"
           "}\n\n")
Q_INJECT = 1.0e-6

ARMS = {
    "allDead":   BASE,
    "topology":  PROBE,
    "shifted":   BASE.replace("A { pressure 3.0e4; }", "A { pressure 3.1e4; }"),
    "knifeEdge": BASE.replace(STK1, STK1.replace("yieldStress 3.0", "yieldStress 2.5")),
    "oneFlowing": BASE.replace(STK1, STK1.replace("yieldStress 3.0", "yieldStress 2.0")),
    "injected":  BASE.replace("pipes\n{", SOURCES + "pipes\n{", 1),
}

out, meta = {}, {}
with tempfile.TemporaryDirectory() as td:
    for name, src in ARMS.items():
        path = os.path.join(td, name + ".dict")
        open(path, "w").write(src)
        r = run(path)
        txt = r.stdout + r.stderr
        require(r.returncode == 0 and "nan" not in txt.lower() and "inf" not in txt.lower(),
                f"{name}: solves with finite output (rc={r.returncode})")
        out[name] = parse(r.stdout)
        meta[name] = case_of(src)
        require(len(out[name]["pipes"]) == len(meta[name]["pipes"]),
                f"{name}: all {len(meta[name]['pipes'])} declared pipes are reported")

# --- 1. the declaration really is "every edge below its own threshold" ------
for name in ("allDead", "topology", "shifted"):
    c = meta[name]
    require(all(p["threshold"] > p["offered"] for p in c["pipes"].values()),
            f"{name}: 4 tau_0 L/d > offered dp on every one of"
            f" {len(c['pipes'])} pipes -- computed from the dict, not remembered")

# --- 2. the whole-network dead fixed point ---------------------------------
for name in ("allDead", "topology", "shifted"):
    o, c = out[name], meta[name]
    s = o["solved"]
    require(s.get("iterations") == 0 and s.get("residual") == 0.0
            and s.get("state") == "converged",
            f"{name}: iterations=0 residual=0 converged -- Newton was never entered")
    require(all(p["Q"] == 0.0 and p["v"] == 0.0 and p["Re"] == 0.0 and p["f"] == 0.0
                for p in o["pipes"].values()),
            f"{name}: every pipe prints Q = v = Re = f = 0 exactly")
    require(all(close(p["dp"], c["pipes"][k]["offered"]) for k, p in o["pipes"].items()),
            f"{name}: each pipe holds exactly the offered drop, no pressure invented")
    require(close(o["nodes"]["E"], c["p_init"]),
            f"{name}: free node E stays at the mean of the fixed pressures"
            f" ({c['p_init']:.1f} Pa), its own equation is 0 = 0")

require(close(out["allDead"]["nodes"]["E"], 25000.0)
        and close(out["shifted"]["nodes"]["E"], 25500.0)
        and out["allDead"]["nodes"]["E"] != out["shifted"]["nodes"]["E"],
        "shifted: E moves 25000 -> 25500 with p(A) -- the fixed point is not a"
        " hardcoded pressure")
require(out["topology"]["net"]["pipes"] == 4 and out["allDead"]["net"]["pipes"] == 5,
        "topology: 4-pipe twin vs 5-pipe case -- the same dead fixed point on"
        " another topology")

# --- 3. the knife edge is a strict inequality ------------------------------
ke = meta["knifeEdge"]["pipes"]["Stk1"]
require(close(ke["threshold"], ke["offered"]) and close(ke["threshold"], 10000.0),
        "knifeEdge: Stk1 tau_0 = 2.5 makes dpStatic == dp == 10000 Pa exactly")
require(out["knifeEdge"]["pipes"]["Stk1"]["Q"] == 0.0
        and out["knifeEdge"]["solved"]["iterations"] == 0,
        "knifeEdge: Q is still exactly 0 at equality -- the law releases only"
        " for dp > dpStatic")

# --- 4. iterations == 0 does NOT mean the whole network is dead ------------
of = out["oneFlowing"]["pipes"]
q1 = of["Stk1"]["Q"]
stk1 = meta["oneFlowing"]["pipes"]["Stk1"]
q_closed = buckingham_q(stk1["offered"], stk1["length"], stk1["diameter"],
                        stk1["yieldStress"])
require(stk1["threshold"] < stk1["offered"],
        "oneFlowing: Stk1 tau_0 = 2.0 puts dpStatic 8000 Pa below the offered 10000 Pa")
require(q1 > 0 and close(q1, q_closed, 1e-9),
        f"oneFlowing: Stk1 Q = {q1:.12g} reproduces Buckingham's exact"
        f" Q = {q_closed:.12g} (1e-9)")
require(all(of[k]["Q"] == 0.0 for k in of if k != "Stk1"),
        "oneFlowing: the other four pipes stay exactly stuck")
require(out["oneFlowing"]["solved"]["iterations"] == 0
        and close(out["oneFlowing"]["nodes"]["E"], 25000.0),
        "oneFlowing: iterations is STILL 0 and E unchanged -- the zero is the"
        " free node's residual, not a 'network is dead' shortcut")
xi = stk1["threshold"] / stk1["offered"]
require(close(of["Stk1"]["f"] * of["Stk1"]["Re"] * g_xi(xi), 64.0, 1e-9),
        "oneFlowing: f Re g(xi) = 64 on the one flowing pipe -- the reported"
        " factor belongs to the solution, not to Churchill")

# --- 5. a source forces Newton to run, and the HP cap holds ----------------
inj, ci = out["injected"], meta["injected"]
require(inj["solved"]["iterations"] > 0 and inj["solved"]["state"] == "converged",
        f"injected: Newton runs {inj['solved']['iterations']} steps and converges"
        " -- the free-node row is invertible (HP cap on dead branches)")
require(inj["nodes"]["E"] > ci["p_init"],
        "injected: E rises above the initial guess to open a path out")
flowing = [k for k, p in inj["pipes"].items() if p["Q"] != 0.0]
require(flowing == ["Sr1b"],
        f"injected: exactly one pipe yields ({flowing}) -- the source leaves"
        " through the lower fixed-pressure side")
sr1b = ci["pipes"]["Sr1b"]
dp_closed = invert_buckingham(Q_INJECT, sr1b["length"], sr1b["diameter"],
                              sr1b["yieldStress"])
dp_printed = inj["pipes"]["Sr1b"]["dp"]
require(close(inj["pipes"]["Sr1b"]["Q"], Q_INJECT, 1e-9),
        "injected: the only flowing rate equals the imposed source 1e-6 m3/s")
require(close(dp_printed, dp_closed, 1e-9),
        f"injected: Sr1b dp = {dp_printed:.11g} is Buckingham's inverse"
        f" {dp_closed:.11g} at that rate, not a fitted number")
require(close(inj["nodes"]["E"], ci["fixed"]["B"] + dp_closed, 1e-9),
        f"injected: E = p(B) + dp = {ci['fixed']['B'] + dp_closed:.11g} --"
        " the node pressure is the closed form's own answer")
require(close(dp_printed, inj["nodes"]["E"] - ci["fixed"]["B"], 1e-9),
        "injected: the pipe's own drop closes with the two node pressures")
require(abs(inj["pipes"]["Sr1a"]["dp"]) < ci["pipes"]["Sr1a"]["threshold"]
        and inj["pipes"]["Sr1a"]["Q"] == 0.0,
        "injected: Sr1a is driven backwards by only 3978.9 Pa and stays stuck")

print()
print(f"{len(failures)} failure(s)")
if failures:
    sys.exit(1)
print("OK: the whole-network dead fixed point is 0 = 0 at the initial guess;"
      " a source reopens it through Buckingham, without changing the kernel")
