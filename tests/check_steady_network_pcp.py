#!/usr/bin/env python3
"""5.191 / G7: the positive-displacement machine as a row of the 1D network.

The increment claims four things, and every one is checkable from the case
text alone:

  1. the row's machine rise is the machine's own slip law read as a head,
       H(q) = (nV - q) / (C rho g),
     so the pipes line's `dpPcp` is exactly -(nV - q)/C at the rate it
     reports -- not a curve fitted to the rate;
  2. in a closed loop with equal end pressures the row's differential is
     zero, so the rate solves R q = (nV - q)/C against the laminar wall,
     R = 128 mu L / (pi D^4), i.e. q* = nV/(1 + R C), and the machine's
     rise is exactly the wall's drop, R q*;
  3. the law is linear in n, so doubling the speed doubles the rate and the
     rise exactly, and doubling the slip halves the denominator's excess --
     both arms are re-derived from the same closed form;
  4. a pipe without the machine carries no `dpPcp` token at all and, at
     equal end pressures, no rate -- the token is the machine's, and nothing
     else's.

Refusals (seven): a non-positive displacement, speed or slip, a pcp beside a
centrifugal pumpHead, beside a compressor, on a coupled edge, and on a
black-oil network.

usage: check_steady_network_pcp.py <steady-exe> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)

STEADY, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
CASE_DIR = os.path.dirname(CASE)

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def run(text_or_path):
    path = text_or_path
    tmp = None
    if "\n" in text_or_path:
        fd, tmp = tempfile.mkstemp(suffix=".dict", dir=CASE_DIR)
        with os.fdopen(fd, "w") as f:
            f.write(text_or_path)
        path = tmp
    p = subprocess.run([STEADY, path], capture_output=True, text=True,
                       cwd=CASE_DIR)
    if tmp:
        os.unlink(tmp)
    return p.returncode, p.stdout + p.stderr


def strip_comments(text):
    return re.sub(r"//[^\n]*", "", text)


CASE_TEXT = strip_comments(open(CASE).read())


def block_of(name, text=CASE_TEXT):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", text)
    if not m:
        return None
    i = m.end()
    depth = 1
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    return text[m.end():i - 1]


def scalar(body, key, default=None):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)\s*;", body)
    return float(m.group(1)) if m else default


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


FLUID = block_of("fluid")
RHO = scalar(FLUID, "rho")
MU = scalar(FLUID, "mu")
PIPES = block_of("pipes")
P1B = block_of("P1", PIPES)
P2B = block_of("P2", PIPES)
PCB = block_of("pcp", P1B)
L = scalar(P1B, "length")
D = scalar(P1B, "diameter")
DEV = scalar(PCB, "displacement")
SPEED = scalar(PCB, "speed")
SLIP = scalar(PCB, "slip")

NV = DEV * SPEED
R_WALL = 128.0 * MU * L / (math.pi * D ** 4)
Q_EXPECT = NV / (1.0 + R_WALL * SLIP)
DP_EXPECT = (NV - Q_EXPECT) / SLIP

# ---- the run ------------------------------------------------------------
RC, OUT = run(CASE)
check("case runs", RC == 0, f"rc={RC}")

PIPE = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] pipe (\S+) (.*)$", line)
    if not m:
        continue
    tok = dict(re.findall(r"(\w+) ([-+0-9.eE]+)", m.group(2)))
    PIPE[m.group(1)] = {k: float(v) for k, v in tok.items()}

P1 = PIPE["P1"]
P2 = PIPE["P2"] if "P2" in PIPE else {}

check("the machine row carries the closed-form rate",
      rel(P1["Q"], Q_EXPECT) < 1e-9,
      f"printed {P1['Q']:.12e} vs nV/(1+RC) = {Q_EXPECT:.12e} "
      f"(nV = {NV:.9e}, R = {R_WALL:.9e}, RC = {R_WALL * SLIP:.9e})")
check("the machine's rise is -(nV - q)/C",
      "dpPcp" in P1 and rel(-P1["dpPcp"], DP_EXPECT) < 1e-9,
      f"printed {-P1.get('dpPcp', float('nan')):.9e} vs {DP_EXPECT:.9e}")
check("the rise prices the wall's drop exactly",
      rel(-P1["dpPcp"], R_WALL * P1["Q"]) < 1e-9,
      f"-dpPcp = {-P1['dpPcp']:.9e}, R q = {R_WALL * P1['Q']:.9e}")
check("the row's differential stays the fixed nodes' (dp = 0)",
      abs(P1["dp"]) == 0.0, f"dp = {P1['dp']}")
check("dpPump is the same rise: one machine, one curve",
      "dpPump" in P1 and rel(P1["dpPump"], P1["dpPcp"]) < 1e-12,
      f"{P1.get('dpPump')} vs {P1.get('dpPcp')}")
check("the reported head is the slip law's",
      rel(-P1["dpPcp"] / (RHO * 9.81), DP_EXPECT / (RHO * 9.81)) < 1e-9,
      f"H = {-P1['dpPcp'] / (RHO * 9.81):.6f} m")

_AREA = math.pi * D * D / 4.0
_RE = RHO * abs(P1["Q"]) * D / (_AREA * MU)
check("the rate is the laminar wall's (Re < 2040 and recomputed)",
      P1["Re"] < 2040 and rel(P1["Re"], _RE) < 1e-9,
      f"printed Re {P1['Re']:.6f} vs {_RE:.6f}")
check("velocity is q over the bore",
      rel(P1["v"], P1["Q"] / _AREA) < 1e-9, f"v = {P1['v']:.9f}")

check("the machine-less pipe reports no dpPcp at all",
      "dpPcp" not in P2 and "dpPump" not in P2, f"{sorted(P2)}")
check("and carries no rate at equal end pressures", abs(P2["Q"]) == 0.0,
      f"Q = {P2.get('Q')}")

# ---- 3. the law is linear in n and 1/(1+RC) in C ------------------------
def with_pcp_key(text, key, value):
    """Set `key value;` inside the pcp block only (the case declares each
    once, but the scope is the machine's, so the rewrite is anchored there)."""
    def sub(m):
        return re.sub(r"\b" + key + r"\s+[-+0-9.eE]+",
                      f"{key} {value}", m.group(0), count=1)
    return re.sub(r"pcp \{[^}]*\}", sub, text, count=1)


def pipe_p1(out):
    for line in out.splitlines():
        m = re.match(r"\[chflow\] pipe P1 (.*)$", line)
        if m:
            return {k: float(v) for k, v in
                    dict(re.findall(r"(\w+) ([-+0-9.eE]+)", m.group(1))).items()}
    return None


T_SPEED = with_pcp_key(CASE_TEXT, "speed", repr(2.0 * SPEED))
RC2, OUT2 = run(T_SPEED)
R2 = pipe_p1(OUT2)
check("doubling the speed doubles the rate exactly (the law is linear in n)",
      R2 is not None and rel(R2["Q"], 2.0 * P1["Q"]) < 1e-9,
      f"{R2 and R2['Q']:.12e} vs {2.0 * P1['Q']:.12e}")
check("and doubles the rise exactly",
      R2 is not None and rel(-R2["dpPcp"], -2.0 * P1["dpPcp"]) < 1e-9,
      f"{-R2['dpPcp']:.9e} vs {-2.0 * P1['dpPcp']:.9e}")

T_SLIP = with_pcp_key(CASE_TEXT, "slip", repr(2.0 * SLIP))
RC3, OUT3 = run(T_SLIP)
Q_S = NV / (1.0 + R_WALL * 2.0 * SLIP)
R3 = pipe_p1(OUT3)
check("doubling the slip lands on nV/(1+2RC)",
      R3 is not None and rel(R3["Q"], Q_S) < 1e-9,
      f"{R3 and R3['Q']:.12e} vs {Q_S:.12e}")

# ---- the machine-less twin ----------------------------------------------
T_NOPCP = re.sub(r"pcp \{[^}]*\}", "", CASE_TEXT)
RC4, OUT4 = run(T_NOPCP)
TWIN = pipe_p1(OUT4)
check("without a pcp the same row carries neither rate nor dpPcp",
      TWIN is not None and abs(TWIN["Q"]) == 0.0 and "dpPcp" not in TWIN,
      f"{TWIN}")

# ---- the tabulated fill: the rise is in Pascals, so density cancels ------
# dp = rho g H = rho g (nV - q)/(C rho g) = (nV - q)/C: the same identity the
# constant-property arm asserts must survive a rho(p) march, where the
# density the head is referred to changes along the line.
PVT_FILL = """fluid
{
    pvt
    {
        pressure (1.0e5 5.0e6 9.0e6 1.3e7);
        rho (5.0 200.0 360.0 520.0);
        mu (1.0e-5 1.4e-5 1.8e-5 2.2e-5);
        temperature 300.0;
    }
}

nodes
{
    A { pressure 5.0e6; }
    B { pressure 5.0e6; }
}

pipes
{
    M { from A; to B; length 500; diameter 0.1; roughness 0.0;
        pcp { displacement 1.2e-4; speed 8.333333333333333; slip 1.0e-11; } }
}

solve { maxIter 100; tol 1e-12; }
"""
RC5, OUT5 = run(PVT_FILL)
M = pipe_p1(OUT5)
if M is None:
    for line in OUT5.splitlines():
        m = re.match(r"\[chflow\] pipe M (.*)$", line)
        if m:
            M = {k: float(v) for k, v in
                 dict(re.findall(r"(\w+) ([-+0-9.eE]+)", m.group(1))).items()}
NV_P = 1.2e-4 * 8.333333333333333
check("the tabulated fill runs a machine row", RC5 == 0 and M is not None,
      f"rc={RC5} {OUT5.strip().splitlines()[-1][:100] if OUT5 else ''}")
check("on the rho(p) march the rise is still (nV - q)/C",
      # nV - q is a difference of two 1e-3 numbers worth 2.4e-10, so ~7 digits
      # cancel: the identity is asserted where it lives, at 1e-6.
      M is not None and rel(-M["dpPcp"], (NV_P - M["Q"]) / 1.0e-11) < 1e-6,
      f"{-M['dpPcp']:.9e} vs {(NV_P - M['Q']) / 1.0e-11:.9e}")
check("the slip always leaks: the driven rate sits below nV",
      M is not None and 0 < M["Q"] < NV_P, f"q = {M['Q']:.12e} < nV = {NV_P:.9e}")
check("with no centrifugal curve on the row, dpPump is not a separate rise",
      M is not None and "dpPump" not in M, f"{sorted(M)}")

# the compressor keys must sit on the PIPE, not inside the pcp block
PVT_COMP = re.sub(
    r"pcp \{[^}]*\}",
    "pressureRatio 2.0; ratedRate 0.05; polytropicEfficiency 0.9;"
    " kappa 1.3; specificGasConstant 518.28; compressibility 0.95;"
    " pcp { displacement 1.2e-4; speed 8.333333333333333; slip 1.0e-11; }",
    PVT_FILL, count=1)

# ---- refusals -----------------------------------------------------------
def refuses(name, text, needle):
    rc, out = run(text)
    last = out.strip().splitlines()[-1][:140] if out else f"rc={rc}"
    check(name, rc != 0 and needle in out, last)


refuses("zero displacement refused",
        with_pcp_key(CASE_TEXT, "displacement", 0),
        "positive displacement")
refuses("zero speed refused", with_pcp_key(CASE_TEXT, "speed", 0),
        "positive speed")
refuses("zero slip refused", with_pcp_key(CASE_TEXT, "slip", 0),
        "positive slip")
refuses("pcp beside a centrifugal pumpHead refused",
        T_NOPCP.replace("roughness 0; }",
                        "roughness 0; pumpHead 20; pumpRate 1.0e-4; pcp {"
                        " displacement 1.2e-5; speed 8.33; slip 1.0e-12; } }"),
        "one row carries one machine curve")
refuses("pcp beside a compressor refused", PVT_COMP, "one row, one machine")
refuses("pcp on a coupled edge refused",
        CASE_TEXT.replace("pipes\n{",
                          "pipes\n{\n    C1 { from NodeA; to NodeB; length 10;"
                          " diameter 0.05; coupledResistance 1.0e6; pcp {"
                          " displacement 1.2e-5; speed 8.0; slip 1.0e-12; } }"),
        "coupled element")
refuses("pcp on a black-oil network refused",
        CASE_TEXT.replace("fluid\n{\n    rho 850; mu 0.1;\n}",
                          "fluid\n{\n    blackoil\n    {\n        oilDensity"
                          " 850;\n        gor 0;\n        waterCut 0;\n"
                          "        gasGravity 0.75;\n        temperature 340;\n"
                          "    }\n}"),
        "black-oil network")

if fails:
    print(f"\n{len(fails)} failure(s): {fails}")
    raise SystemExit(1)
print("\nall pcp network claims hold")
