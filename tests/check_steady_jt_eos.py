#!/usr/bin/env python3
"""5.189: the Joule-Thomson coefficient taken from the fluid's own EOS.

Two claims, both re-derivable from the case text and the printed readout:

  1. the march's temperature drop is the integral of the coefficient it
     printed, over the pressure grid it printed:
       T_out = T_in + sum_k mu_k (p_{k+1} - p_k),
     so `thermal <pipe> pt` / `jtpt` lines are the march itself, not a
     decoration next to it;
  2. `jouleThomson eosSrk` prices mu from SRK's analytic dZ/dT:
       mu = R T^2 dZ/dT / (p M c_p),
     which is re-derived here from the published SRK constants and the
     cubic's implicit derivative, for every printed station.

The Peng-Robinson branch (`jouleThomson eos`) is NOT recomputed here: this
layer has no analytic PR derivative (the PR branch takes a central
difference of the state), so the checker only pins that it lands in a
documented band around the SRK value at the same state - a sanity band, not
an identity claim. The identity and the derivative chain are anchored in
tests/unit/srkTests.cpp instead.

Cross-checks: the declared-constant pipe (P3) still prices by the constant
and prints no `jt` line, and every mu is positive and smooth along the pipe
(a gas cooling on expansion).

Refusals: `jouleThomson eos` on a fluid without a composition, without a
positive heatCapacity, an unknown word, and the word on a plain PVT table.

usage: check_steady_jt_eos.py <steady-exe> <case.dict>
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


def run(path):
    p = subprocess.run([STEADY, path], capture_output=True, text=True,
                       cwd=os.path.dirname(path) or ".")
    return p.returncode, p.stdout + p.stderr


CASE_TEXT = re.sub(r"//[^\n]*", "", open(CASE).read())


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


def numbers(key, text):
    m = re.search(re.escape(key) + r"\s*\(([^)]*)\)", text)
    return [float(t) for t in m.group(1).split()] if m else []


# ---- the published constants the binary also uses (Components.cpp) ----
COMPONENTS = {
    "methane": dict(M=16.043e-3, Tc=190.564, Pc=45.992e5, w=0.0115),
    "ethane": dict(M=30.070e-3, Tc=305.322, Pc=48.722e5, w=0.0995),
}
R = 8.31446261815324
OMEGA_A, OMEGA_B = 0.42748023354034, 0.08664007007534

feed = {}
comp = block_of("composition")
for name in COMPONENTS:
    m = re.search(r"\b" + name + r"\s+([0-9.eE+-]+)", comp or "")
    if m:
        feed[name] = float(m.group(1))
if not feed:
    raise SystemExit("case has no composition: " + CASE)


def srk_mu(T, p, cp):
    """The SRK analytic Joule-Thomson coefficient, K/Pa."""
    x = {k: v / sum(feed.values()) for k, v in feed.items()}
    ai, bi, dai = {}, {}, {}
    for k, v in COMPONENTS.items():
        if k not in x:
            continue
        m = 0.480 + 1.574 * v["w"] - 0.176 * v["w"] ** 2
        tr = T / v["Tc"]
        alpha = (1.0 + m * (1.0 - math.sqrt(tr))) ** 2
        a0 = OMEGA_A * (R * v["Tc"]) ** 2 / v["Pc"]
        ai[k] = a0 * alpha
        bi[k] = OMEGA_B * R * v["Tc"] / v["Pc"]
        dai[k] = a0 * (-m * (1.0 + m * (1.0 - math.sqrt(tr)))
                       / math.sqrt(T * v["Tc"]))
    a = b = da = 0.0
    for i in x:
        b += x[i] * bi[i]
        for j in x:
            aij = math.sqrt(ai[i] * ai[j])       # no declared k_ij in the case
            daij = 0.5 * (dai[i] * ai[j] + ai[i] * dai[j]) / math.sqrt(
                ai[i] * ai[j])
            a += x[i] * x[j] * aij
            da += x[i] * x[j] * daij
    A = a * p / (R * T) ** 2
    B = b * p / (R * T)
    dA = da * p / (R * T) ** 2 - 2.0 * A / T
    dB = -B / T
    # Z^3 - Z^2 + (A - B - B^2) Z - A B = 0, single stable root here
    roots = cubic_roots(-1.0, A - B - B * B, -(A * B))
    Z = max(roots)
    FZ = 3.0 * Z * Z - 2.0 * Z + (A - B - B * B)
    dZ = -((dA - (1.0 + 2.0 * B) * dB) * Z - (A * dB + B * dA)) / FZ
    molar = sum(x[k] * COMPONENTS[k]["M"] for k in x)
    return R * T * T * dZ / (p * molar * cp)


def cubic_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p ** 3 / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        s = math.sqrt(disc)
        return [math.cbrt(-q / 2.0 + s) + math.cbrt(-q / 2.0 - s) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p ** 3 / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    th = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    return sorted(m * math.cos((th + 2 * k * math.pi) / 3.0) + shift
                  for k in range(3))


cp_of = {}
for name in ("P1", "P2", "P3"):
    body = block_of(name) or ""
    m = re.search(r"heatCapacity\s+([0-9.eE+-]+)", body)
    jt = re.search(r"jouleThomson\s+(\S+)\s*;", body)
    cp_of[name] = (float(m.group(1)) if m else None,
                   jt.group(1) if jt else None)

rc, out = run(CASE)
check("the case solves", rc == 0, f"(exit {rc})")

# ---- the printed readout, per pipe ----
profiles, jt_lines, jtpts, pipes = {}, {}, {}, {}
for line in out.splitlines():
    m = re.match(r"\[chflow\] thermal (\w+) pt (\d+) x (\S+) p (\S+) T (\S+)",
                 line)
    if m:
        name, k = m.group(1), int(m.group(2))
        profiles.setdefault(name, {})[k] = (float(m.group(4)),
                                            float(m.group(5)))
        continue
    m = re.match(r"\[chflow\] thermal (\w+) jt (eos\w*) mu (\S+)", line)
    if m:
        jt_lines[m.group(1)] = (m.group(2), float(m.group(3)))
        continue
    m = re.match(r"\[chflow\] thermal (\w+) jtpt (\d+) p (\S+) T (\S+) mu (\S+)",
                 line)
    if m:
        jtpts.setdefault(m.group(1), {})[int(m.group(2))] = (
            float(m.group(3)), float(m.group(4)), float(m.group(5)))
        continue
    m = re.match(r"\[chflow\] pipe (\w+) .* T (\S+)$", line)
    if m:
        pipes[m.group(1)] = float(m.group(2))

check("P1 declares eos (Peng-Robinson)", jt_lines.get("P1", ("", 0))[0] == "eos",
      str(jt_lines.get("P1")))
check("P2 declares eosSrk", jt_lines.get("P2", ("", 0))[0] == "eosSrk",
      str(jt_lines.get("P2")))
check("P3 (declared constant) prints no jt line", "P3" not in jt_lines)

TIN = 320.0
for name, cp, word in (("P1", cp_of["P1"][0], "eos"),
                       ("P2", cp_of["P2"][0], "eosSrk")):
    prof = profiles.get(name, {})
    pts = jtpts.get(name, {})
    check(f"{name} prints the 24-station mu curve", len(pts) == 24,
          f"(got {len(pts)})")
    if not (prof and pts and cp):
        continue
    # 1. the drop is the integral of the printed coefficient
    drops = 0.0
    for k in range(24):
        p_next = prof[k + 1][0]
        drops += pts[k][2] * (p_next - pts[k][0])
    lhs = pipes[name]
    rhs = TIN + drops
    check(f"{name}: T_out = T_in + sum mu dp", abs(lhs - rhs) < 1e-9,
          f"({lhs:.9f} vs {rhs:.9f})")
    # 2. mu positive, and smooth (a gas cooling on expansion)
    mus = [pts[k][2] for k in range(24)]
    check(f"{name}: mu > 0 everywhere", all(m > 0 for m in mus))
    check(f"{name}: mu is smooth along the pipe",
          max(abs(mus[k + 1] - mus[k]) / mus[k] for k in range(23)) < 5e-3)

# 3. the SRK branch re-derived here from the published constants
if cp_of["P2"][0]:
    worst = 0.0
    for k, (p, T, mu) in jtpts.get("P2", {}).items():
        ref = srk_mu(T, p, cp_of["P2"][0])
        worst = max(worst, abs(mu - ref) / abs(ref))
    check("P2: every station mu equals the SRK recomputation (1e-9 relat.)",
          worst < 1e-9, f"(worst {worst:.3e})")

# 4. the PR branch, sanity band against the same state's SRK value
if cp_of["P1"][0]:
    worst = 0.0
    for k, (p, T, mu) in jtpts.get("P1", {}).items():
        ref = srk_mu(T, p, cp_of["P1"][0])
        worst = max(worst, abs(mu - ref) / abs(ref))
    check("P1: PR mu tracks the SRK value at the same state", worst < 0.08,
          f"(worst {worst:.3f}; band, not an identity)")
    check("P1: PR and SRK are not the same number",
          worst > 1e-3, f"(worst {worst:.3f})")

# ---- refusals ----
probes = [
    ("no composition", "steadyJtEosNoComposition",
     "fluid\n{\n    pvt\n    {\n        pressure ( 1.0e6 3.0e6 );\n"
     "        rho ( 10 30 );\n        mu ( 1e-5 1.2e-5 );\n    }\n}\n"
     "nodes { W {} T { pressure 2.0e6; } }\n"
     "pipes { L { from W; to T; length 100; diameter 0.1;"
     " jouleThomson eos; heatCapacity 2400; } }\n"
     "sources { f { node W; rate 0.01; } }\nsolve { maxIter 100; tol 1e-10; }\n",
     "composition"),
    ("no heatCapacity", "steadyJtEosNoCp",
     "fluid\n{\n    pvt\n    {\n        pressure ( 1.0e6 3.0e6 );\n"
     "        temperature 320.0;\n        compositional\n        {\n"
     "            composition { methane 1.0; }\n            temperature 320.0;\n"
     "        }\n    }\n}\n"
     "nodes { W {} T { pressure 2.0e6; } }\n"
     "pipes { L { from W; to T; length 100; diameter 0.1;"
     " jouleThomson eos; } }\n"
     "sources { f { node W; rate 0.01; } }\nsolve { maxIter 100; tol 1e-10; }\n",
     "heatCapacity"),
    ("unknown word", "steadyJtEosBogus",
     "fluid\n{\n    pvt\n    {\n        pressure ( 1.0e6 3.0e6 );\n"
     "        temperature 320.0;\n        compositional\n        {\n"
     "            composition { methane 1.0; }\n            temperature 320.0;\n"
     "        }\n    }\n}\n"
     "nodes { W {} T { pressure 2.0e6; } }\n"
     "pipes { L { from W; to T; length 100; diameter 0.1;"
     " jouleThomson bogus; heatCapacity 2400; } }\n"
     "sources { f { node W; rate 0.01; } }\nsolve { maxIter 100; tol 1e-10; }\n",
     "unknown jouleThomson"),
    ("black-oil fluid", "steadyJtEosBlackOil",
     "fluid\n{\n    blackoil\n    {\n        gasGravity 0.75; oilDensity 850;"
     " gor 80; waterCut 0;\n        temperature 320;\n    }\n}\n"
     "nodes { W {} T { pressure 2.0e6; } }\n"
     "pipes { L { from W; to T; length 100; diameter 0.1;"
     " jouleThomson eos; heatCapacity 2400; } }\n"
     "sources { f { node W; rate 0.004; } }\nsolve { maxIter 100; tol 1e-10; }\n",
     "composition"),
]
for label, stem, text, needle in probes:
    path = os.path.join(tempfile.mkdtemp(), stem + ".dict")
    with open(path, "w") as fh:
        fh.write(text)
    rc, out2 = run(path)
    check(f"refused: {label}", rc != 0 and needle in out2,
          f"(exit {rc}, needle {'hit' if needle in out2 else 'missing'})")

print(f"\n{len(fails)} failed" if fails else "\nsteadyJtEos: all checks passed")
sys.exit(1 if fails else 0)
