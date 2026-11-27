#!/usr/bin/env python3
"""5.223 (H3): independent verification of the separator AS a network row.

The steady network solve is re-derived from the case text alone, with no
solver internals:

  * mass closure: liquidMassRate + gasMassRate == feedMassRate (1e-10), and
    the two edge volumes sum to the feed the source delivers (1e-9);
  * the dP statement: the solved upstream-downstream node pressure difference
    IS the declared dP (1e-9 relative) at BOTH outlets, which are equal by
    construction; the row's Jacobian presence is what produced it, and the
    node balance closing on the line demands proves it;
  * the twin case: the same feed flashed at the solved stage pressure through
    the independent chflow-separator app must agree on rs/bo/bg/gasLiberated
    (1e-9);
  * the liberated-gas identity gasLiberated == gor - rs, re-derived from the
    reported rs and the case's gor;
  * the full rejection set: efficiency outside [0,1], negative dP, a declared
    non-positive feedRate, an unknown gasNode and a converged row with no
    feed at all.

usage: check_network_separator_row.py <steadyExe> <caseDict> <separatorExe>
"""
import os
import re
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

STEADY_EXE, CASE, SEP_EXE = sys.argv[1:4]
CASE = os.path.abspath(CASE)

fails = []
NCHECK = [0]


def check(name, cond, detail=""):
    NCHECK[0] += 1
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol=1e-9, atol=1e-12):
    return abs(got - want) <= rtol * abs(want) + atol


def run(exe, case, expect_fail=False):
    r = subprocess.run([exe, case], capture_output=True, text=True)
    if expect_fail:
        return r.returncode != 0, r.stderr
    return r.returncode, r.stdout + r.stderr


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


def parse(out):
    nodes, pipes, sep = {}, {}, None
    for line in out.splitlines():
        m = re.search(r"node (\w+) p ([-+0-9.eE]+)", line)
        if m:
            nodes[m.group(1)] = float(m.group(2))
        m = re.search(r"pipe (\w+) Q ([-+0-9.eE]+).*dp ([-+0-9.eE]+)", line)
        if m:
            pipes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        if " separator SEP " in line:
            sep = {k: float(v) for k, v in
                   re.findall(r"(\w+) ([-+0-9.eE]+)", line)}
    return nodes, pipes, sep


src = open(CASE).read()
bo = block_body(src, "blackoil")
gor = scalar(bo, "gor")
feed_decl = scalar(block_body(src, "sources"), "rate")
sep_body = block_body(block_body(src, "pipes"), "separator")
eff = scalar(sep_body, "efficiency")
dp_decl = scalar(sep_body, "dP")
t_stage = scalar(sep_body, "temperature")

rc, out = run(STEADY_EXE, CASE)
check("solve converged", rc == 0 and "converged" in out, out[-200:] if rc else "")
nodes, pipes, sep = parse(out)

# --- mass closure (1e-10) ---
check("mass closure |liquid+gas-feed| <= 1e-10",
      abs(sep["liquid"] + sep["gas"] - sep["feed"]) <= 1e-10,
      f"closure={sep['closure']:.3e}")
check("reported closure token is the same residual",
      abs(sep["closure"] - (sep["liquid"] + sep["gas"] - sep["feed"])) <= 1e-16)

# --- volume: the two edges sum to the feed the source delivers ---
check("edges sum to the source feed",
      near(sep["qLiq"] + sep["qGas"], feed_decl, 1e-9),
      f"{sep['qLiq']+sep['qGas']!r} vs {feed_decl!r}")
check("row flowRate token equals the declared volume feed",
      near(pipes["SEP"][0], feed_decl, 1e-9))
# feed mass rate re-derived: oil+water stock-tank mass of the feed
rho_o = scalar(bo, "oilDensity")
rho_w = scalar(bo, "waterDensity", 1000.0)
wc = scalar(bo, "waterCut", 0.0)
# rhoGasStd re-derived from the published ideal-gas standard density
R_GAS, M_AIR, P_STD, T_STD = 8.314462618, 28.9647e-3, 101325.0, 288.15
rho_gas_std = P_STD * scalar(bo, "gasGravity") * M_AIR / (R_GAS * T_STD)
wor = wc / (1.0 - wc) if wc > 0 else 0.0
feed_mass = feed_decl * (rho_o + gor * rho_gas_std + wor * rho_w)
check("feedMassRate == feed*(rho_o + gor*rho_g + wor*rho_w)",
      near(sep["feed"], feed_mass, 1e-9), f"{sep['feed']!r} vs {feed_mass!r}")
check("entrained == (1-eff)/eff * liquid",
      near(sep["entrained"], (1.0 - eff) / eff * sep["liquid"], 1e-9))
lib_mass = sep["gas"] - sep["entrained"]
check("liberated mass == (gor-rs)*feed*rho_g",
      near(lib_mass, (gor - sep["rs"]) * feed_decl * rho_gas_std, 1e-9))

# --- dP enters the solve: both outlet node differences ARE dP ---
check("p[IN] - p[OUT] == dP",
      abs(nodes["IN"] - nodes["OUT"] - dp_decl) <= 1e-9 * dp_decl,
      f"got {nodes['IN']-nodes['OUT']!r}")
check("p[IN] - p[GAS] == dP",
      abs(nodes["IN"] - nodes["GAS"] - dp_decl) <= 1e-9 * dp_decl,
      f"got {nodes['IN']-nodes['GAS']!r}")
check("the two outlets are equal",
      abs(nodes["OUT"] - nodes["GAS"]) <= 1e-9 * nodes["OUT"])
check("row dp token equals dP",
      abs(pipes["SEP"][1] - dp_decl) <= 1e-9 * dp_decl)

# --- the lines carry the edges ---
check("LIQ carries the liquid edge", near(pipes["LIQ"][0], sep["qLiq"], 1e-9))
check("GASL carries the gas edge", near(pipes["GASL"][0], sep["qGas"], 1e-9))

# --- liberated-gas identity, re-derived ---
check("gasLiberated == gor - rs",
      near(sep["gasLiberated"], gor - sep["rs"], 1e-9),
      f"{sep['gasLiberated']!r} vs {gor - sep['rs']!r}")
check("qGasStd == (gor - rs) * feed",
      near(sep["qGasStd"], (gor - sep["rs"]) * feed_decl, 1e-9))

# --- stage pressure is the inlet node pressure ---
check("stageP == p[IN]", near(sep["stageP"], nodes["IN"], 1e-12, 1e-6))
check("outletP == stageP - dP",
      near(sep["outletP"], sep["stageP"] - dp_decl, 1e-12, 1e-6))

# --- twin against the independent chflow-separator app ---
twin = f"""separator {{
  pressure ({sep['stageP']!r});
  temperature {t_stage!r};
  fluid {{ blackoil {{ gasGravity {scalar(bo, 'gasGravity')};
      oilDensity {rho_o}; gor {gor}; waterCut {wc};
      waterDensity {rho_w}; temperature {t_stage}; }} }}
}}
"""
with open("/tmp/_cf_sep_twin.dict", "w") as fh:
    fh.write(twin)
rc2, out2 = run(SEP_EXE, "/tmp/_cf_sep_twin.dict")
check("twin separator app ran", rc2 == 0)
m = re.search(r"stage 0 .*rs ([-+0-9.eE]+) gasSm3 ([-+0-9.eE]+).*"
              r"bo ([-+0-9.eE]+) bg ([-+0-9.eE]+)", out2)
check("twin parsed", m is not None)
if m:
    t_rs, t_gl, t_bo, t_bg = (float(x) for x in m.groups())
    check("twin rs (1e-9)", near(sep["rs"], t_rs, 1e-9), f"{sep['rs']!r} vs {t_rs!r}")
    check("twin gasLiberated (1e-9)", near(sep["gasLiberated"], t_gl, 1e-9))
    check("twin bo (1e-9)", near(sep["bo"], t_bo, 1e-9))
    check("twin bg (1e-9)", near(sep["bg"], t_bg, 1e-9))

# --- rejections ---
bad_eff = src.replace("efficiency 0.9;", "efficiency 1.5;")
with open("/tmp/_cf_bad_eff.dict", "w") as fh:
    fh.write(bad_eff)
failed, err = run(STEADY_EXE, "/tmp/_cf_bad_eff.dict", True)
check("reject efficiency > 1", failed and "efficiency" in err, err[-120:])

bad_dp = src.replace("dP 1.0e4;", "dP -1.0;")
with open("/tmp/_cf_bad_dp.dict", "w") as fh:
    fh.write(bad_dp)
failed, err = run(STEADY_EXE, "/tmp/_cf_bad_dp.dict", True)
check("reject dP < 0", failed and "dP" in err, err[-120:])

bad_feed = src.replace("feedRate 2.0e-3;", "feedRate 0.0;")
with open("/tmp/_cf_bad_feed.dict", "w") as fh:
    fh.write(bad_feed)
failed, err = run(STEADY_EXE, "/tmp/_cf_bad_feed.dict", True)
check("reject declared non-positive feedRate", failed and "feedRate" in err,
      err[-120:])

bad_gas = src.replace("gasNode GAS;", "gasNode NOPE;")
with open("/tmp/_cf_bad_gas.dict", "w") as fh:
    fh.write(bad_gas)
failed, err = run(STEADY_EXE, "/tmp/_cf_bad_gas.dict", True)
check("reject unknown gasNode", failed and "NOPE" in err, err[-120:])

# no feed at all: a sourceless network whose every node rests at the same
# pressure delivers nothing to the stage -> the row converges with q_in = 0
nofeed = src.replace("sources { s1 { node IN; rate 2.0e-3; } }",
                     "sources { }")
nofeed = nofeed.replace("feedRate 2.0e-3;", "")
nofeed = nofeed.replace("dP 1.0e4;", "dP 0.0;")
nofeed = nofeed.replace("TANK { pressure 1.0e5; }",
                        "TANK { pressure 1.0e5; }\n    SRC { pressure 1.0e5; }")
nofeed = nofeed.replace("pipes {\n    SEP {",
                        "pipes {\n    FEED { from SRC; to IN; length 500.0;"
                        " diameter 0.05; roughness 0.0; }\n    SEP {")
with open("/tmp/_cf_nofeed.dict", "w") as fh:
    fh.write(nofeed)
failed, err = run(STEADY_EXE, "/tmp/_cf_nofeed.dict", True)
check("reject a row with no feed rate",
      failed and "no feed rate" in err, err[-160:])

print(f"\n{NCHECK[0] - len(fails)}/{NCHECK[0]} checks passed")
if fails:
    raise SystemExit("FAILED: " + ", ".join(fails))
