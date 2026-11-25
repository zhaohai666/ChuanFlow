#!/usr/bin/env python3
"""5.192 / G8: a jet pump as a row of the 1D network (the two-stream slot).

The increment makes four claims, and each is re-derived here from the case
text alone, never from the report it is checking:

  1. THE ROW LAW IS THE ELEMENT'S OWN THROAT BALANCE. With both streams
     drawn at the row's upstream node, the rate the nodal balance solves is
     the COMBINED q_d = q_n + q_s, so the entrained share is q - q_n and
        H(q) = rho v_n^2/m^2 [ (m-1-r)^2/(m-1) - F (1+r)^2/2 ],
        r = (q - q_n)/q_n, m = A_t/A_n, v_n = q_n/A_n, F = f L_t/D_t,
     which is `jetPumpSolve`'s imposed-rate mode. In a loop with equal end
     pressures the row's differential is zero, so the balance is
        R q = H(q),  R = 128 mu L/(pi D^4),
     and the root is re-solved here by bisection -- the printer's q, rise,
     r, qs and supply are all compared against that root.
  2. THE TWO STREAMS ADD UP. q_s = q_d - q_n, q_d is the pipe's reported
     rate, m and v_n are the declared geometry, and `supply` is the motive
     pressure the nozzle needs at the row's own upstream node
     (p_s + rho v_n^2/2). A card that does not conserve those is not the
     machine this row claims to carry.
  3. THE ELEMENT SOLVER AGREES. The `jet` line also carries `cardRise`, the
     rise `jetPumpSolve` returns for the same declared suction rate, and
     `resid`, ITS momentum residual: both must match the row's rise to
     round-off. This is the cross-check against the independent card -- the
     row is not a second implementation of the algebra, it is the same one.
  4. THE ROW IS MONOTONE IN ITS WALL. Weakening the wall (L 10 -> 5 -> 0.01)
     raises the entrainment ratio and the rate monotonically, and the
     printed r always stays strictly under the free-jet ceiling r_free
     (bracket(r_free) = 0), which is the element's own envelope.

Two clamp arms: a viscous enough liquid leaves NO-SUCTION (r clamps to 0,
the stall head holds the wall, and qs is reported NEGATIVE because the row
carries less than the declared nozzle stream -- that is the reduction's own
edge and it says so rather than fabricating reverse entrainment). FREE-JET
(H <= 0) is the clamp's other guard: it is NOT reachable from a positive
pressure balance (the root always sits where the bracket is still positive),
so nothing here pretends to have tested it.

Refusals (nine): a non-positive nozzle area, a throat no larger than the
nozzle, a non-positive motive rate, negative throat length / friction,
beside a centrifugal pumpHead, beside a pcp, beside a compressor, on a
coupled edge, and on a black-oil network.

usage: check_steady_network_jet_pump.py <steady-exe> <case.dict>
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
G = 9.81

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def run(text_or_path):
    path, tmp = text_or_path, None
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
    i, depth = m.end(), 1
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
JB = block_of("J1", PIPES)
JET = block_of("jetPump", JB)
L = scalar(JB, "length")
D = scalar(JB, "diameter")

QN = scalar(JET, "motiveRate")
AN = scalar(JET, "nozzleArea")
AT = scalar(JET, "throatArea")
LT = scalar(JET, "throatLength")
FF = scalar(JET, "frictionFactor")

M = AT / AN
VN = QN / AN
DT = math.sqrt(4.0 * AT / math.pi)
FL = FF * LT / DT
R_WALL = 128.0 * MU * L / (math.pi * D ** 4)


def bracket(r):
    s = M - 1.0 - r
    return s * s / (M - 1.0) - 0.5 * FL * (1.0 + r) * (1.0 + r)


def head(r):
    """Metres of the combined stream -- what the row law subtracts."""
    return VN * VN / (M * M * G) * bracket(r)


def rise(r):
    """Pascals -- what the card and the pipe line report."""
    return RHO * VN * VN / (M * M) * bracket(r)


def solve_root(r_wall):
    """The row's own balance, R q = rise(r), by bisection on r."""
    r_free = None
    lo, hi = 0.0, M - 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bracket(mid) > 0.0:
            lo = mid
        else:
            hi = mid
            r_free = mid
    def resid(r):
        return r_wall * QN * (1.0 + r) - rise(r)
    lo, hi = 0.0, (r_free if r_free else M - 1.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if resid(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), (r_free if r_free else M - 1.0), resid


R_ROOT, R_FREE, RESID = solve_root(R_WALL)
Q_ROOT = QN * (1.0 + R_ROOT)

# ---- the run ------------------------------------------------------------
RC, OUT = run(CASE)
check("case runs", RC == 0, f"rc={RC}")

PIPE = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] pipe (\S+) (.*)$", line)
    if m:
        PIPE[m.group(1)] = {k: float(v) for k, v in
                            dict(re.findall(r"(\w+) ([-+0-9.eE]+)",
                                            m.group(2))).items()}
JETLINE = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] jet (\S+) (.*)$", line)
    if m:
        body = m.group(2)
        state = re.match(r"state (\S+) ", body).group(1)
        JETLINE[m.group(1)] = dict(
            state=state,
            **{k: float(v) for k, v in
               dict(re.findall(r"(\w+) ([-+0-9.eE]+)", body)).items()})

NODE_P = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] node (\S+) p (\S+)", line)
    if m:
        NODE_P[m.group(1)] = float(m.group(2))

J1, JL = PIPE["J1"], JETLINE.get("J1", {})
check("the machine row reports a jet line", bool(JL), f"{sorted(JL)}")

# ---- 1. the root --------------------------------------------------------
check("the rate is the row law's own root, re-solved",
      rel(J1["Q"], Q_ROOT) < 1e-9,
      f"printed {J1['Q']:.12e} vs R q = H(q) at {Q_ROOT:.12e} "
      f"(r = {R_ROOT:.9f}, R = {R_WALL:.9e})")
check("the entrainment ratio is (q_d - q_n)/q_n",
      rel(JL["r"], (J1["Q"] - QN) / QN) < 1e-9,
      f"{JL['r']:.9f} vs {(J1['Q'] - QN) / QN:.9f}")
check("the rise is the throat balance at that r",
      rel(JL["rise"], rise(R_ROOT)) < 1e-9,
      f"printed {JL['rise']:.9e} Pa vs {rise(R_ROOT):.9e} Pa")
check("the machine's rise pays the wall's drop exactly (dp = 0 loop)",
      rel(JL["rise"], R_WALL * J1["Q"]) < 1e-9,
      f"rise {JL['rise']:.9e} vs R q {R_WALL * J1['Q']:.9e}")
check("dpJet is the same rise with the row's sign",
      rel(-J1["dpJet"], JL["rise"]) < 1e-12, f"{J1['dpJet']} vs {JL['rise']}")
check("dpPump is the same rise: one machine on the row",
      rel(J1["dpPump"], J1["dpJet"]) < 1e-12,
      f"{J1['dpPump']} vs {J1['dpJet']}")
check("the row's differential stays the fixed nodes' (dp = 0)",
      abs(J1["dp"]) == 0.0, f"dp = {J1['dp']}")

# ---- 2. the two streams -------------------------------------------------
check("q_s + q_n is the row's rate (mass closes on one rate)",
      # the report prints ~12 significant digits, so the sum is asserted
      # against that print precision, not against the kernel's full word.
      rel(JL["qs"] + JL["qn"], J1["Q"]) < 1e-10,
      f"{JL['qs']} + {JL['qn']} vs {J1['Q']}")
check("q_s is the printed entrained stream",
      rel(JL["qs"], J1["Q"] - QN) < 1e-10, f"{JL['qs']:.12e}")
check("q_n and v_n are the declared nozzle stream",
      abs(JL["qn"] - QN) == 0.0 and rel(JL["vn"], VN) < 1e-12,
      f"qn {JL['qn']}, vn {JL['vn']} vs {VN}")
check("m is the declared area ratio", abs(JL["m"] - M) == 0.0,
      f"m = {JL['m']} vs {M}")
SUP = NODE_P["Suction"] + 0.5 * RHO * VN * VN
check("supply is the motive pressure the nozzle needs upstream",
      rel(JL["supply"], SUP) < 1e-12,
      f"printed {JL['supply']:.6f} vs p_s + rho v_n^2/2 = {SUP:.6f}")

# ---- 3. the element solver agrees --------------------------------------
check("the row is inside the element's envelope (PUMPING)",
      JL["state"] == "PUMPING" and 0.0 < JL["r"] < R_FREE,
      f"state {JL['state']}, r {JL['r']:.6f} < r_free {R_FREE:.6f}")
check("jetPumpSolve's card rise is the row's rise", 
      "cardRise" in JL and rel(JL["cardRise"], JL["rise"]) < 1e-9,
      f"card {JL.get('cardRise')} vs row {JL['rise']:.9e}")
check("and its momentum residual is zero at that rate",
      "resid" in JL and abs(JL["resid"]) <= 1e-6 * max(abs(JL["rise"]), 1.0),
      f"resid {JL.get('resid')}")
check("the efficiency and head ratio are the element's own numbers",
      "eff" in JL and "N" in JL
      and rel(JL["eff"], JL["r"] * JL["N"]) < 1e-9,
      f"eff {JL.get('eff')} vs r N = {JL['r'] * JL['N']:.6e}")

# ---- 4. monotone in the wall -------------------------------------------
def jet_of(text):
    rc, out = run(text)
    p = j = None
    for line in out.splitlines():
        m = re.match(r"\[chflow\] pipe J1 (.*)$", line)
        if m:
            p = {k: float(v) for k, v in
                 dict(re.findall(r"(\w+) ([-+0-9.eE]+)", m.group(1))).items()}
        m = re.match(r"\[chflow\] jet J1 (.*)$", line)
        if m:
            body = m.group(1)
            j = dict(state=re.match(r"state (\S+) ", body).group(1),
                     **{k: float(v) for k, v in
                        dict(re.findall(r"(\w+) ([-+0-9.eE]+)",
                                        body)).items()})
    return p, j


TEN = re.sub(r"length 5; diameter", "length 10; diameter", CASE_TEXT, count=1)
SHORT = re.sub(r"length 5; diameter", "length 0.01; diameter", CASE_TEXT,
               count=1)
P10, J10 = jet_of(TEN)
PSH, JSH = jet_of(SHORT)
check("a heavier wall entrains less: r falls 5 -> 10 m",
      P10 is not None and J10 is not None and J10["r"] < JL["r"],
      f"r(10 m) = {J10 and J10['r']:.6f} vs r(5 m) = {JL['r']:.6f}")
check("and the rate and the rise follow the root, not a clamp",
      P10 is not None and J10 is not None
      and rel(P10["Q"], QN * (1.0 + J10["r"])) < 1e-9
      and rel(J10["rise"], R_WALL * (10.0 / 5.0) * P10["Q"]) < 1e-9,
      f"q {P10 and P10['Q']:.9e}, rise {J10 and J10['rise']:.6e}")
check("a vanishing wall entrains more: r rises 5 -> 0.01 m",
      PSH is not None and JSH is not None and JSH["r"] > JL["r"],
      f"r(0.01 m) = {JSH and JSH['r']:.6f} vs r(5 m) = {JL['r']:.6f}")
check("and still stays strictly under the free-jet ceiling",
      PSH is not None and JSH is not None and JSH["r"] < R_FREE,
      f"r {JSH and JSH['r']:.9f} < r_free {R_FREE:.9f}")

# ---- the clamp arm: NO-SUCTION says so, it does not invent a rate -------
VISCOUS = CASE_TEXT.replace("mu 0.1", "mu 0.5", 1)
PV, JV = jet_of(VISCOUS)
R_V = 128.0 * 0.5 * L / (math.pi * D ** 4)
Q_V = rise(0.0) / R_V
check("a wall past the stall head clamps to NO-SUCTION",
      JV is not None and JV["state"] == "NO-SUCTION" and JV["r"] == 0.0,
      f"state {JV and JV['state']}, r {JV and JV['r']}")
check("the stall head holds that wall: q = H(0)/R",
      PV is not None and rel(PV["Q"], Q_V) < 1e-9,
      f"printed {PV and PV['Q']:.12e} vs {Q_V:.12e}")
check("and the entrained share is reported NEGATIVE, not fabricated",
      JV is not None and JV["qs"] < 0.0 and rel(JV["qs"], PV["Q"] - QN) < 1e-10,
      f"qs {JV and JV['qs']:.6e} (q_d {PV and PV['Q']:.6e} < q_n {QN})")
check("a clamped row still writes no machine number it cannot price",
      JV is not None and "cardRise" not in JV and "resid" not in JV,
      f"{sorted(JV or {})}")

# ---- the machine-less pipe ---------------------------------------------
check("the machine-less pipe reports no dpJet at all",
      "dpJet" not in PIPE.get("P2", {}) and "P2" not in JETLINE,
      f"{sorted(PIPE.get('P2', {}))}")
check("and carries no rate at equal end pressures",
      abs(PIPE["P2"]["Q"]) == 0.0, f"Q = {PIPE['P2']['Q']}")

# ---- refusals -----------------------------------------------------------
def refuses(name, text, needle):
    rc, out = run(text)
    last = out.strip().splitlines()[-1][:140] if out else f"rc={rc}"
    check(name, rc != 0 and needle in out, last)


def with_jet_key(text, key, value):
    def sub(m):
        return re.sub(r"\b" + key + r"\s+[-+0-9.eE]+", f"{key} {value}",
                      m.group(0), count=1)
    return re.sub(r"jetPump \{[^}]*\}", sub, text, count=1)


NOPUMP = re.sub(r"jetPump \{[^}]*\}", "", CASE_TEXT, count=1)
refuses("zero nozzle area refused",
        with_jet_key(CASE_TEXT, "nozzleArea", 0), "positive nozzleArea")
refuses("a throat no larger than the nozzle refused",
        with_jet_key(CASE_TEXT, "throatArea", AN), "throatArea > nozzleArea")
refuses("zero motive rate refused",
        with_jet_key(CASE_TEXT, "motiveRate", 0), "positive motiveRate")
refuses("negative throat length refused",
        with_jet_key(CASE_TEXT, "throatLength", -0.1),
        "length and friction factor are non-negative")
refuses("jet beside a centrifugal pumpHead refused",
        NOPUMP.replace("roughness 0; }",
                       "roughness 0; pumpHead 20; pumpRate 1.0e-4; jetPump {"
                       " motiveRate 8.0e-4; nozzleArea 1.0e-4; throatArea"
                       " 4.0e-4; throatLength 0.2; frictionFactor 0.02; } }"),
        "one row carries one machine curve")
refuses("jet beside a pcp refused",
        NOPUMP.replace("roughness 0; }",
                       "roughness 0; pcp { displacement 1.2e-5; speed 8.0;"
                       " slip 1.0e-12; } jetPump { motiveRate 8.0e-4;"
                       " nozzleArea 1.0e-4; throatArea 4.0e-4; throatLength"
                       " 0.2; frictionFactor 0.02; } }"),
        "one row carries one machine curve")
PVT_CASE = CASE_TEXT.replace("fluid\n{\n    rho 850; mu 0.1;\n}",
                              "fluid\n{\n    pvt\n    {\n        pressure"
                              " (1.0e5 5.0e6 9.0e6 1.3e7);\n        rho"
                              " (5.0 200.0 360.0 520.0);\n        mu (1.0e-5"
                              " 1.4e-5 1.8e-5 2.2e-5);\n        temperature"
                              " 300.0;\n    }\n}")
# the compressor keys sit on the PIPE (outside the jetPump block) and a
# compressor needs the tabulated fill, so this arm is read on a pvt case
refuses("jet beside a compressor refused",
        re.sub(r"jetPump \{[^}]*\}",
               "pressureRatio 2.0; ratedRate 0.05; polytropicEfficiency 0.9;"
               " kappa 1.3; specificGasConstant 518.28; compressibility 0.95;"
               " jetPump { motiveRate 8.0e-4; nozzleArea 1.0e-4; throatArea"
               " 4.0e-4; throatLength 0.2; frictionFactor 0.02; }",
               PVT_CASE, count=1),
        "one row carries one machine curve")
refuses("jet on a coupled edge refused",
        NOPUMP.replace("pipes\n{",
                       "pipes\n{\n    C1 { from NodeA; to NodeB; length 10;"
                       " diameter 0.05; coupledResistance 1.0e6; jetPump {"
                       " motiveRate 8.0e-4; nozzleArea 1.0e-4; throatArea"
                       " 4.0e-4; throatLength 0.2; frictionFactor 0.02; } }"),
        "coupled element")
refuses("jet on a black-oil network refused",
        CASE_TEXT.replace("fluid\n{\n    rho 850; mu 0.1;\n}",
                          "fluid\n{\n    blackoil\n    {\n        oilDensity"
                          " 850;\n        gor 0;\n        waterCut 0;\n"
                          "        gasGravity 0.75;\n        temperature 340;\n"
                          "    }\n}"),
        "black-oil network")

if fails:
    print(f"\n{len(fails)} failure(s): {fails}")
    raise SystemExit(1)
print("\nall jet-pump row claims hold")
