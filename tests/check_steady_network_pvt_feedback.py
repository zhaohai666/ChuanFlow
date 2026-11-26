#!/usr/bin/env python3
"""Verification of TEMPERATURE FEEDBACK on the pvt march (5.194, G6d).

The increment: the steady temperature march used to be PASSIVE -- it priced T
along the pressure profile and nothing read the result, so every station's
rho/mu came from the feed's single-temperature snapshot. With
`thermal{ temperatureFeedback true; }` the segment's properties are evaluated
at the T the march has reached, and the mass rate is booked at the entry state
(pUp, T_in). T now moves dp, and dp moves T.

Criterion makeup. The fluid is a light single-phase gas, so the checker can
re-derive every property from the published equations instead of trusting the
kernel's flash:

  P1. THE STATION STATE IS A REAL FLASH AT (p, T). Every printed
      `flashStation` rho/mu is re-derived in Python from the Peng-Robinson
      cubic (PR76 alpha, vdW mixing, the gas root) for density and the
      Lohrenz-Bray-Clark correlation on the flash's own molar volume for
      viscosity. On the FEEDBACK line the T used is the printed per-station
      T; on the FEEDBACK-OFF line it is the feed's 330 K. Both must agree to
      1e-9, and the two lines' rho values must visibly differ -- otherwise the
      switch would be a no-op that happens to parse.

  P2. THE RATE IS THE MARCH PRICED ON THOSE STATES. Marching the case's own
      24 segments (properties from the P1 relation at each printed (p, T),
      Churchill friction at the local Reynolds number, dz = 0 here because the
      line is level) and bisecting the rate that lands the outlet on the
      fixed 1 MPa reproduces the printed Q to 1e-9. This is the criterion that
      makes the COUPLING falsifiable, not just the declaration.

  P3. THE MASS IS BOOKED AT THE ENTRY STATE. `pvtfb` prints rhoIn and m;
      rhoIn must be the flash at (pUp, T_in) to 1e-9, m = rhoIn Q, and the
      outlet identity rho_last v_out A = m must close to 1e-12. A march that
      used the table density at the inlet while pricing the segments at the
      marched T would fail exactly here.

  P4. THE SWITCH IS HONOURED IN BOTH DIRECTIONS. The `attrib` line (no thermal
      block, no switch) and the FEEDBACK-OFF line (`false`) must be identical
      token for token -- an absent key is the quiet default -- and neither may
      print a `pvtfb` line.

  P5. THE GRAMMAR IS REFUSED. Four malformed declarations must fail loudly:
      a non-boolean word, feedback on a fluid without a live flash (a
      tabulated rho(p)/mu(p) is one temperature's snapshot), feedback on a
      black-oil network (its layer is isothermal by construction), and
      feedback with no heat declared at all (inert).

usage: check_steady_network_pvt_feedback.py <exe> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

EXE, CASE = sys.argv[1], sys.argv[2]
SEGMENTS = 24
G = 9.80665
R = 8.314462618
SQRT2 = 1.4142135623730951
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846
PR_ZC = 0.3074

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def run(path, cwd=None):
    p = subprocess.run([EXE, path], capture_output=True, text=True,
                       cwd=cwd or (os.path.dirname(path) or "."))
    return p.returncode, p.stdout + p.stderr


def strip_comments(text):
    return re.sub(r"//[^\n]*", "", text)


def block_of(name, text):
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


def numlist(body, key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", body)
    return [float(x) for x in m.group(1).split()] if m else []


# ------------------------------------------------------------ published EOS
CT = {  # (MW g/mol, Tc K, Pc Pa, acentric) -- the component table the kernel uses
    "methane": (16.043, 190.564, 45.992e5, 0.0115),
    "ethane": (30.070, 305.322, 48.722e5, 0.0995),
}


def cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def cubic_real_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p ** 3 / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        sq = math.sqrt(disc)
        return [cbrt(-q / 2.0 + sq) + cbrt(-q / 2.0 - sq) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p ** 3 / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    return sorted(m * math.cos((theta + 2.0 * k * math.pi) / 3.0) + shift
                  for k in range(3))


def pr_kappa(w):
    return 0.37464 + 1.54226 * w - 0.26992 * w * w


def pr_z(T, p, zs, names):
    """Gas root of the PR76 cubic with vdW mixing (kij = 0)."""
    aI, bI = [], []
    for n in names:
        _, Tc, Pc, w = CT[n]
        k = pr_kappa(w)
        alpha = (1.0 + k * (1.0 - math.sqrt(T / Tc))) ** 2
        aI.append(OMEGA_A * (R * Tc) ** 2 / Pc * alpha)
        bI.append(OMEGA_B * R * Tc / Pc)
    n = len(names)
    aMix = sum(zs[i] * zs[j] * math.sqrt(aI[i] * aI[j])
               for i in range(n) for j in range(n))
    bMix = sum(z * b for z, b in zip(zs, bI))
    A = aMix * p / (R * T) ** 2
    B = bMix * p / (R * T)
    roots = [x for x in cubic_real_roots(-(1.0 - B), A - 3.0 * B * B - 2.0 * B,
                                         -(A * B - B * B - B ** 3))
             if x > B * (1.0 + 1e-12)]
    if not roots:
        raise SystemExit("PR cubic has no physical root")
    return max(roots)


ATM = 101325.0


def xi(Tc, Pc, MW):
    return Tc ** (1.0 / 6.0) / (math.sqrt(MW) * (Pc / ATM) ** (2.0 / 3.0))


def stiel_thodos(T, Tc, Pc, MW):
    Tr = T / Tc
    x = xi(Tc, Pc, MW)
    cp = (17.78e-5 * (4.58 * Tr - 1.67) ** 0.625 / x if Tr > 1.5
          else 34.0e-5 * Tr ** 0.94 / x)
    return cp * 1.0e-3


def lbc_mu(T, Vm, zs, names):
    """Lohrenz-Bray-Clark on the PR molar volume, Vc from the PR critical Zc."""
    Tcs = [CT[n][1] for n in names]
    Pcs = [CT[n][2] for n in names]
    MWs = [CT[n][0] for n in names]
    Vcs = [PR_ZC * R * CT[n][1] / CT[n][2] for n in names]
    Tc = sum(z * t for z, t in zip(zs, Tcs))
    MW = sum(z * m for z, m in zip(zs, MWs))
    Pc = sum(z * p for z, p in zip(zs, Pcs))
    Vc = sum(z * v for z, v in zip(zs, Vcs))
    rr = Vc / Vm
    num = sum(z * stiel_thodos(T, t, p, m) * math.sqrt(m)
              for z, t, p, m in zip(zs, Tcs, Pcs, MWs))
    den = sum(z * math.sqrt(m) for z, m in zip(zs, MWs))
    mu0_cp = num / den * 1.0e3
    poly = (0.1023 + 0.023364 * rr + 0.058533 * rr ** 2
            - 0.040758 * rr ** 3 + 0.0093324 * rr ** 4)
    return (mu0_cp + (poly ** 4 - 1.0e-4) / xi(Tc, Pc, MW)) * 1.0e-3


def state(T, p, zs, names):
    Z = pr_z(T, p, zs, names)
    MW = sum(z * CT[n][0] for z, n in zip(zs, names)) * 1.0e-3  # kg/mol
    Vm = Z * R * T / p
    return MW / Vm, lbc_mu(T, Vm, zs, names)


def churchill(re, eps_rel):
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * eps_rel))) ** 16.0
    b = (37530.0 / re) ** 16.0
    return 8.0 * ((8.0 / re) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


# ---------------------------------------------------------------- the case
TEXT = strip_comments(open(CASE).read())
RC, OUT = run(CASE)
check("the case runs", RC == 0, f"rc={RC}")
if RC != 0:
    print(OUT)
    raise SystemExit(1)

PIPES_BODY = block_of("pipes", TEXT)
FLUID = block_of("pvt", block_of("fluid", TEXT))
COMP = block_of("compositional", FLUID)
FEED = {k: float(v) for k, v in
        re.findall(r"\b([A-Za-z]+)\s+([-+0-9.eE]+)\s*;", block_of("composition", COMP))}
NAMES = list(FEED)
NORM = sum(FEED.values())
ZS = [FEED[k] / NORM for k in NAMES]
T_IN = scalar(COMP, "temperature")
P_TABLE = numlist(FLUID, "pressure")

PIPE = {}
FLASH = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] pipe (\S+) Q (\S+) v (\S+)", line)
    if m:
        d = PIPE.setdefault(m.group(1), {"Q": float(m.group(2)),
                                         "v": float(m.group(3))})
        for key in ("Re", "f", "dp", "T"):
            mm = re.search(r" " + key + r" (\S+)", line)
            if mm:
                d[key] = float(mm.group(1))
        continue
    m = re.match(r"\[chflow\] flashStation (\S+) i (\d+) p (\S+) rho (\S+) "
                 r"mu (\S+) sigma (\S+)", line)
    if m:
        FLASH.setdefault(m.group(1), []).append(
            (int(m.group(2)), float(m.group(3)), float(m.group(4)),
             float(m.group(5))))
    m = re.match(r"\[chflow\] thermal (\S+) pt (\d+) x (\S+) p (\S+) T (\S+)", line)
    if m:
        FLASH.setdefault("_T_" + m.group(1), []).append(
            (int(m.group(2)), float(m.group(5))))

for name in ("FB", "NO"):
    FLASH["_T_" + name] = sorted(FLASH.get("_T_" + name, []))

PIPE_T = {k[3:]: [(i, t) for i, t in v] for k, v in FLASH.items()
          if k.startswith("_T_")}

PVTFB = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] pvtfb (\S+) Tin (\S+) rhoIn (\S+) m (\S+)", line)
    if m:
        PVTFB[m.group(1)] = dict(Tin=float(m.group(2)),
                                 rhoIn=float(m.group(3)),
                                 m=float(m.group(4)))

check("three lines report a rate", sorted(PIPE) == ["FB", "NO", "attrib"],
      f"{sorted(PIPE)}")
check("the feedback line prints its entry state, the off lines do not",
      sorted(PVTFB) == ["FB"], f"{sorted(PVTFB)}")

# ------------------------------------------------- P1: the station states
def station_states(name, t_of):
    """re-derive (rho, mu) at every printed station of a line"""
    out = []
    for idx, p, rho, mu in FLASH[name]:
        T = t_of(idx, p)
        rho_w, mu_w = state(T, p, ZS, NAMES)
        out.append(dict(i=idx, p=p, T=T, rho=rho, mu=mu,
                        rho_w=rho_w, mu_w=mu_w))
    return out


FB = station_states("FB", lambda i, p: PIPE_T["FB"][i][1])
NO = station_states("NO", lambda i, p: T_IN)
check("each line prints one station per march segment",
      len(FLASH["FB"]) == SEGMENTS and len(FLASH["NO"]) == SEGMENTS,
      f"{len(FLASH['FB'])} / {len(FLASH['NO'])}")
check("the thermal profile is one point per segment plus the exit",
      len(PIPE_T["FB"]) == SEGMENTS + 1
      and len(PIPE_T["NO"]) == SEGMENTS + 1,
      f"{len(PIPE_T['FB'])} / {len(PIPE_T['NO'])}")

for tag, rows, want_T in (("FB", FB, None), ("NO", NO, T_IN)):
    worst_rho = max(rel(r["rho"], r["rho_w"]) for r in rows)
    worst_mu = max(rel(r["mu"], r["mu_w"]) for r in rows)
    check(f"{tag}: every station rho is the PR gas root at its own (p, T)",
          worst_rho < 1e-9, f"worst rel {worst_rho:.3e}")
    check(f"{tag}: every station mu is LBC on that same state",
          worst_mu < 1e-9, f"worst rel {worst_mu:.3e}")

T_LAST = PIPE_T["FB"][-1][1]
check("the feedback line really cools before the last station",
      T_IN - T_LAST > 50.0, f"Tin {T_IN} -> Tlast {T_LAST:.3f} K")
drho = abs(FB[-1]["rho"] - NO[-1]["rho"]) / NO[-1]["rho"]
check("the cooled state is a different density, not a rounding",
      drho > 0.05, f"rho_last {FB[-1]['rho']:.6f} vs {NO[-1]['rho']:.6f}"
                   f" ({drho * 100:.2f}%)")

# ------------------------------------- P2: the march on those printed states
LEN = scalar(block_of("FB", PIPES_BODY), "length")
DIA = scalar(block_of("FB", PIPES_BODY), "diameter")
ROUGH = scalar(block_of("FB", PIPES_BODY), "roughness", 0.0)
P_UP = 3.0e6
P_LOW = 1.0e6
AREA = 0.25 * math.pi * DIA ** 2


P_FLOOR = max(0.02 * P_UP, 1000.0)


def march(q, T_of):
    """restate the kernel march, floor and early return included: below the
    floor the flash is meaningless, so the kernel stops there and the residual
    turns into a one-way signal instead of an extrapolated nonsense pressure"""
    rho_in, _ = state(T_of(0, P_UP), P_UP, ZS, NAMES)
    mass = rho_in * q
    p = P_UP
    ds = LEN / SEGMENTS
    for i in range(SEGMENTS):
        rho, mu = state(T_of(i, p), p, ZS, NAMES)
        v = mass / (rho * AREA)
        f = churchill(rho * v * DIA / mu, ROUGH / DIA)
        p -= f * ds / DIA * rho * v * v / 2.0
        if p < P_FLOOR:
            return P_FLOOR
    return p


def rate_for(T_of):
    def resid(q):
        return march(q, T_of) - P_LOW
    lo, hi = 1e-6, 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def t_fb(i, p):
    return PIPE_T["FB"][i][1]


def t_no(i, p):
    return T_IN


q_fb = rate_for(t_fb)
q_no = rate_for(t_no)
check("FB: the restated coupled march lands on the printed rate",
      rel(q_fb, PIPE["FB"]["Q"]) < 1e-9,
      f"march {q_fb:.12f} vs printed {PIPE['FB']['Q']:.12f} m3/s")
check("NO: the restated snapshot march lands on the printed rate",
      rel(q_no, PIPE["NO"]["Q"]) < 1e-9,
      f"march {q_no:.12f} vs printed {PIPE['NO']['Q']:.12f} m3/s")

dq = (PIPE["FB"]["Q"] - PIPE["NO"]["Q"]) / PIPE["NO"]["Q"]
check("the coupling moves the rate by a lot, not a rounding",
      dq > 0.05, f"Q(FB)/Q(NO) - 1 = {dq * 100:.2f}%")

# ------------------------------- P3: the mass is booked at the entry state
rho_in_w, mu_in_w = state(T_IN, P_UP, ZS, NAMES)
check("rhoIn is the flash at the entry state (pUp, T_in)",
      rel(PVTFB["FB"]["rhoIn"], rho_in_w) < 1e-9,
      f"printed {PVTFB['FB']['rhoIn']:.12f} vs PR {rho_in_w:.12f}")
check("Tin is the declared feed temperature",
      rel(PVTFB["FB"]["Tin"], T_IN) < 1e-15, f"{PVTFB['FB']['Tin']} K")
check("m = rhoIn Q",
      rel(PVTFB["FB"]["m"], PVTFB["FB"]["rhoIn"] * PIPE["FB"]["Q"]) < 1e-12,
      f"{PVTFB['FB']['m']:.12f} vs {PVTFB['FB']['rhoIn'] * PIPE['FB']['Q']:.12f}"
      " kg/s")
check("the outlet closes rho v A = m (mass conserved through the repricing)",
      rel(FB[-1]["rho"] * PIPE["FB"]["v"] * AREA, PVTFB["FB"]["m"]) < 1e-12,
      f"{FB[-1]['rho'] * PIPE['FB']['v'] * AREA:.12f} vs "
      f"{PVTFB['FB']['m']:.12f} kg/s")
check("the outlet Reynolds and friction are the last station's own",
      rel(PIPE["FB"]["Re"], FB[-1]["rho"] * PIPE["FB"]["v"] * DIA / FB[-1]["mu"])
      < 1e-12,
      f"{PIPE['FB']['Re']:.6f} vs "
      f"{FB[-1]['rho'] * PIPE['FB']['v'] * DIA / FB[-1]['mu']:.6f}")

# --------------------------------- P4: the switch in the other direction
for token in ("Q", "v", "Re", "f", "dp"):
    check(f"'false' and 'absent' agree on {token}",
          rel(PIPE["NO"][token], PIPE["attrib"][token]) < 1e-15,
          f"{PIPE['NO'][token]:.12g} vs {PIPE['attrib'][token]:.12g}")
check("the line without a thermal block prints no thermal profile nor entry state",
      "attrib" not in PIPE_T and "attrib" not in PVTFB,
      f"thermal {sorted(PIPE_T)} entry {sorted(PVTFB)}")
check("that line still flashes at the feed temperature, station for station",
      FLASH["attrib"] == FLASH["NO"],
      f"{len(FLASH['attrib'])} stations, rho/mu/p identical to NO:"
      f" {FLASH['attrib'] == FLASH['NO']}")

# ------------------------------------------- P5: the grammar is refused
def refuse(label, text, needle):
    d = tempfile.mkdtemp()
    try:
        path = os.path.join(d, "case.dict")
        open(path, "w").write(text)
        rc, out = run(path, cwd=d)
        ok = rc != 0 and needle in out
        last = "" if ok else " :: " + out.strip().splitlines()[-1]
        check(f"refused: {label}", ok, f"rc={rc}{last}")
    finally:
        shutil.rmtree(d)


def sub1(pattern, repl, text):
    new, n = re.subn(pattern, repl, text, count=1)
    if n != 1:
        raise AssertionError(f"pattern not found: {pattern}")
    return new


refuse("a non-boolean word",
       sub1(r"temperatureFeedback true;", "temperatureFeedback maybe;", TEXT),
       "must be true or false")
refuse("feedback on a fluid without a live flash",
       sub1(r"liveFlash true;", "", TEXT),
       "one temperature's snapshot")
refuse("feedback with no heat at all",
       sub1(r"wallHeatTransfer 200; wallTemperature 260; heatCapacity 2400;\n"
            r"        jouleThomson eosSrk;\n"
            r"        temperatureFeedback true;",
            "temperatureFeedback true;", TEXT),
       "declares thermal keys but no heat")

BO = """fluid
{
    blackoil
    {
        gasGravity 0.75; oilDensity 850; gor 80; waterCut 0;
        temperature 340; surfaceTension 0.025;
    }
}

nodes
{
    A { pressure 6.0e6; }
    B { pressure 4.0e6; }
}

pipes
{
    P
    {
        from A; to B;
        length 1000; diameter 0.1; roughness 4.6e-5;
        jouleThomson 2.0e-6; heatCapacity 2000;
        temperatureFeedback true;
    }
}

solve { maxIter 100; tol 1e-10; }
"""
refuse("feedback on a black-oil network", BO, "black-oil layer is isothermal")

print()
if fails:
    print(f"{len(fails)} criterion(s) failed")
    raise SystemExit(1)
print("all criteria passed")
