#!/usr/bin/env python3
"""5.234 (H3): independent verification of the QMOM pipe block.

The moment march is re-derived in Python from the case text alone: 24
segments, residence time dt = ds/|v| with v = massRate/(rho A), the
constant-kernel Smoluchowski source of 5.107

    dM0/dt = -beta M0^2 / 2,  dM1/dt = 0,
    dM2/dt = beta M1^2,       dM3/dt = 3 beta M1 M2,

and the Sauter-gated Brinkman viscosity mu (1 + 2.5 phi) for d32 = M3/M2 past
dThreshold. The whole per-segment friction march (Churchill 1977) is replayed
from the solved rate, so the reported dp is checked against a hydraulic
re-derivation that only knows the case text and the reported Q -- proving the
emulsion viscosity really charged the gradient.

Arms:
  * no-flow arm (A fixed at pC, sources removed): the moments do not evolve;
  * coagulation arm: M0 monotonically down, M1 conserved to 1e-10 (the
    5.107 invariant), outlet moments matching the replayed march;
  * d32 feedback arm: muEff == muBase (1 + 2.5 phi) above the threshold,
    muEff == muBase below it, and a higher phi raises muEff;
  * an undeclared twin run carries no qmom line at all;
  * rejections: negative M0, missing kernel, dThreshold <= 0, no phi
    statement, qmom on a non-PVT fill.

usage: check_network_qmom.py <steadyExe> <caseDict>
"""
import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

STEADY_EXE, CASE = sys.argv[1:3]
CASE = os.path.abspath(CASE)

fails = []
NCHECK = [0]


def check(name, cond, detail=""):
    NCHECK[0] += 1
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol=1e-10, atol=0.0):
    return abs(got - want) <= rtol * abs(want) + atol


def scalar(src, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)\s*;", src)
    if m:
        return float(m.group(1))
    if default is None:
        raise SystemExit("case: missing " + key)
    return default


def block_body(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        return None
    start = src.index("{", m.start())
    depth = 0
    for j in range(start, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start + 1:j]
    raise SystemExit("unbalanced " + key)


def numlist(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", src)
    return [float(x) for x in m.group(1).split()]


def churchill(re_, rr):
    re_ = max(re_, 1e-12)
    a1 = (7.0 / re_) ** 0.9 + 0.27 * rr
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re_) ** 16.0
    return 8.0 * ((8.0 / re_) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def table_lookup(pressures, values, p):
    # the kernel's clamped linear pvtLookup, re-derived
    if p <= pressures[0]:
        return values[0]
    if p >= pressures[-1]:
        return values[-1]
    for i in range(len(pressures) - 1):
        if pressures[i] <= p <= pressures[i + 1]:
            t = (p - pressures[i]) / (pressures[i + 1] - pressures[i])
            return values[i] + t * (values[i + 1] - values[i])
    raise AssertionError


def parse(out):
    nodes, pipes, qmom = {}, {}, None
    for line in out.splitlines():
        m = re.search(r"node (\w+) p ([-+0-9.eE]+)", line)
        if m:
            nodes[m.group(1)] = float(m.group(2))
        m = re.search(r"pipe (\w+) Q ([-+0-9.eE]+).*dp ([-+0-9.eE]+)", line)
        if m:
            pipes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        m = re.search(r"qmom (\w+) m0 ([-+0-9.eE]+) m1 ([-+0-9.eE]+) "
                      r"m2 ([-+0-9.eE]+) m3 ([-+0-9.eE]+) d32 ([-+0-9.eE]+) "
                      r"muBase ([-+0-9.eE]+) muEff ([-+0-9.eE]+)", line)
        if m:
            qmom = dict(zip(("m0", "m1", "m2", "m3", "d32", "muBase",
                             "muEff"), (float(x) for x in m.groups()[1:])))
    return nodes, pipes, qmom


def run(case, expect_fail=False):
    r = subprocess.run([STEADY_EXE, case], capture_output=True, text=True)
    if expect_fail:
        return r.returncode != 0, r.stderr
    return r.returncode, r.stdout + r.stderr


src = open(CASE).read()
pvt = block_body(src, "pvt")
pTab = numlist(pvt, "pressure")
rhoTab = numlist(pvt, "rho")
muTab = numlist(pvt, "mu")
pipes_body = block_body(src, "pipes")
p1 = block_body(pipes_body, "P1")
qm = block_body(p1, "qmom")
L, D = scalar(p1, "length"), scalar(p1, "diameter")
ROUGH = scalar(p1, "roughness", 0.0)
beta = scalar(qm, "beta")
m = [scalar(qm, "m0"), scalar(qm, "m1"), scalar(qm, "m2"), scalar(qm, "m3")]
thr = scalar(qm, "dThreshold")
phi = scalar(qm, "phi") if "phi " in qm else None
phiFromM0 = scalar(qm, "phiFromM0") if "phiFromM0" in qm else None
mC = re.search(r"C\s*\{\s*pressure ([-+0-9.eE]+)", src).group(1)
pC = float(mC)

rc, out = run(CASE)
check("solve converged", rc == 0 and "converged" in out)
nodes, pipes, qmom = parse(out)
check("qmom line present", qmom is not None)
pA = nodes["A"]  # the feed node is free: its pressure is a solved answer

# --- replay the march from the case text and the solved rate ---
q = pipes["P1"][0]
area = 0.25 * math.pi * D * D
nSeg = 24
ds = L / nSeg
rhoEntry = table_lookup(pTab, rhoTab, pA)
massRate = rhoEntry * q
mm = list(m)
p = pA
dp_sum = 0.0
mu_base = mu_eff = None
d32 = None
for i in range(nSeg):
    rho = table_lookup(pTab, rhoTab, p)
    v = massRate / (rho * area)
    if v != 0.0:
        dm = [-0.5 * beta * mm[0] * mm[0], 0.0,
              beta * mm[1] * mm[1], 3.0 * beta * mm[1] * mm[2]]
        dt = ds / abs(v)
        mm = [mm[k] + dt * dm[k] for k in range(4)]
    d32 = mm[3] / mm[2] if mm[2] > 0 else 0.0
    ph = phiFromM0 * mm[0] if phiFromM0 else phi
    mu0 = table_lookup(pTab, muTab, p)
    mu_base = mu0
    mu_eff = mu0 * (1.0 + 2.5 * ph) if d32 > thr else mu0
    re_ = rho * v * D / mu_eff
    f = churchill(re_, ROUGH / D)
    dp_i = (f * ds / D) * rho * v * v / 2.0
    dp_sum += dp_i
    p -= dp_i

check("outlet M0 matches the replayed march", near(qmom["m0"], mm[0]),
      f"{qmom['m0']!r} vs {mm[0]!r}")
check("outlet M1 matches", near(qmom["m1"], mm[1]))
check("outlet M2 matches", near(qmom["m2"], mm[2]))
check("outlet M3 matches", near(qmom["m3"], mm[3]))
check("d32 == M3/M2", near(qmom["d32"], mm[3] / mm[2]))

# --- 5.107 identities ---
check("M1 conserved to 1e-10", abs(qmom["m1"] - m[1]) <= 1e-10 * abs(m[1]),
      f"{qmom['m1']!r} vs {m[1]!r}")
check("M0 monotonically down", qmom["m0"] < m[0],
      f"{qmom['m0']!r} < {m[0]!r}")
check("M2 monotonically up", qmom["m2"] > m[2])
check("M3 monotonically up", qmom["m3"] > m[3])

# --- the Brinkman feedback charged the gradient ---
check("muBase is the clean table viscosity",
      near(qmom["muBase"], mu_base, 1e-9))
check("muEff == muBase (1 + 2.5 phi) above threshold",
      near(qmom["muEff"], mu_base * (1.0 + 2.5 * phi), 1e-9))
check("replayed dp matches the reported pipe dp",
      near(pipes["P1"][1], dp_sum, 1e-9),
      f"{pipes['P1'][1]!r} vs {dp_sum!r}")
check("node balance closes at B",
      near(pipes["P1"][0], pipes["P2"][0], 1e-9))

# --- no-flow arm: A fixed equal to C, sources dropped -> frozen moments ---
noflow = re.sub(r"(?m)^\s*sources\s*\{.*\}\s*$", "", src)
noflow = noflow.replace("A {};", f"A {{ pressure {mC}; }}")
with open("/tmp/_cf_qmom_noflow.dict", "w") as fh:
    fh.write(noflow)
rc2, out2 = run("/tmp/_cf_qmom_noflow.dict")
check("no-flow solve converged", rc2 == 0 and "converged" in out2)
_, pipes2, qmom2 = parse(out2)
check("no-flow rate is zero", abs(pipes2["P1"][0]) < 1e-12,
      f"Q={pipes2['P1'][0]!r}")
check("no-flow moments do not evolve",
      all(near(qmom2[k], m[i]) for i, k in
          enumerate(("m0", "m1", "m2", "m3"))),
      f"{qmom2['m0']!r} vs {m[0]!r}")

# --- feedback arms: threshold gate and phi slope ---
below = src.replace("dThreshold 5.0e-10;", "dThreshold 1.0e-3;")
with open("/tmp/_cf_qmom_below.dict", "w") as fh:
    fh.write(below)
rc3, out3 = run("/tmp/_cf_qmom_below.dict")
_, _, qmom3 = parse(out3)
check("below threshold muEff == muBase",
      rc3 == 0 and near(qmom3["muEff"], qmom3["muBase"], 1e-12),
      f"{qmom3['muEff']!r} vs {qmom3['muBase']!r}")

bigphi = src.replace("phi 0.1;", "phi 0.2;")
with open("/tmp/_cf_qmom_bigphi.dict", "w") as fh:
    fh.write(bigphi)
rc4, out4 = run("/tmp/_cf_qmom_bigphi.dict")
_, pipes4, qmom4 = parse(out4)
check("higher phi raises muEff",
      rc4 == 0 and qmom4["muEff"] > qmom["muEff"]
      and near(qmom4["muEff"], qmom4["muBase"] * 1.5, 1e-9))
check("higher phi raises the line's dp",
      pipes4["P1"][1] > pipes["P1"][1])

# --- phi priced from M0 ---
bym0 = src.replace("phi 0.1;", "phiFromM0 1.0e-10;")
with open("/tmp/_cf_qmom_m0.dict", "w") as fh:
    fh.write(bym0)
rc5, out5 = run("/tmp/_cf_qmom_m0.dict")
_, _, qmom5 = parse(out5)
check("phiFromM0 prices phi = coeff * M0(outlet)",
      rc5 == 0 and near(qmom5["muEff"],
                        qmom5["muBase"] * (1.0 + 2.5 * 1.0e-10 * qmom5["m0"]),
                        1e-9), f"{qmom5['muEff']!r}")

# --- undeclared twin: no qmom line anywhere ---
plain = re.sub(r"\n\s*qmom \{[^}]*\}", "", src)
with open("/tmp/_cf_qmom_plain.dict", "w") as fh:
    fh.write(plain)
rc6, out6 = run("/tmp/_cf_qmom_plain.dict")
check("undeclared twin runs", rc6 == 0 and "converged" in out6)
check("undeclared twin prints no qmom line", "qmom " not in out6)

# --- rejections ---
bad = src.replace("m0 1.0e10;", "m0 -1.0;")
with open("/tmp/_cf_qmom_negm0.dict", "w") as fh:
    fh.write(bad)
failed, err = run("/tmp/_cf_qmom_negm0.dict", True)
check("reject M0 < 0", failed and "m0" in err, err[-120:])

bad = re.sub(r"beta 1\.0e-18;\s*#", "#", src)
with open("/tmp/_cf_qmom_nokernel.dict", "w") as fh:
    fh.write(bad)
failed, err = run("/tmp/_cf_qmom_nokernel.dict", True)
check("reject missing kernel", failed and "beta" in err, err[-120:])

bad = src.replace("dThreshold 5.0e-10;", "dThreshold 0.0;")
with open("/tmp/_cf_qmom_thr0.dict", "w") as fh:
    fh.write(bad)
failed, err = run("/tmp/_cf_qmom_thr0.dict", True)
check("reject dThreshold <= 0", failed and "dThreshold" in err, err[-120:])

bad = src.replace("phi 0.1;", "")
with open("/tmp/_cf_qmom_nophi.dict", "w") as fh:
    fh.write(bad)
failed, err = run("/tmp/_cf_qmom_nophi.dict", True)
check("reject no phi statement", failed and "phi" in err, err[-120:])

bo_case = """
fluid { blackoil { gasGravity 0.75; oilDensity 850.0; gor 100.0;
    waterCut 0.1; temperature 320.0; } }
nodes { A { pressure 3.0e6; } C { pressure 1.0e6; } }
pipes { P1 { from A; to C; length 600.0; diameter 0.1;
    qmom { beta 1.0e-18; m0 1.0e10; m1 10.0; m2 1.0e-8; m3 1.0e-17;
           dThreshold 5.0e-10; phi 0.1; } } }
solve { maxIter 10; tol 1e-9; }
"""
with open("/tmp/_cf_qmom_bo.dict", "w") as fh:
    fh.write(bo_case)
failed, err = run("/tmp/_cf_qmom_bo.dict", True)
check("reject qmom on a black-oil network", failed and "qmom" in err,
      err[-140:])

print(f"\n{NCHECK[0] - len(fails)}/{NCHECK[0]} checks passed")
if fails:
    raise SystemExit("FAILED: " + ", ".join(fails))
