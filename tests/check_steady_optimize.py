#!/usr/bin/env python3
"""5.160 / P21: the steady diameter optimization and the gas-lift budget
allocation.

Both blocks are PURE APP-LAYER reports on top of the existing kernels: the
diameter search re-solves the network once per enumerated combination, the
allocation samples the nodal operating point per well per gas step (the
injection gas rides the black-oil GOR through a fixed point in q) and solves
a multiple-choice knapsack. An undeclared case must not gain a single byte.
This checker re-derives everything independently:

  * diameter: every combination is re-run through the solver itself, the
    feasible filter, the objective ranking, the tie-by-enumeration-order and
    the cost model are re-picked in Python and compared to the report lines;
    five objective/constraint variants each re-derive their own best;
  * gas lift: the whole nodal chain is restated in Python (the PVT, the
    Beggs-Brill march, the expanding-bracket + Illinois-secant inlet solve,
    the largest-root operating point, the gas/GOR fixed point), the DP is
    re-solved, and the exchange optimality of the reported allocation is
    verified against the re-derived curves;
  * byte-identity of the stripped twins and of the standing cases;
  * twenty-two refusals of the two declaration grammars.

usage: check_steady_optimize.py <steady-exe> <diameter.dict> <gasLift.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 4:
    raise SystemExit(__doc__)

STEADY = os.path.abspath(sys.argv[1])
DCASE = os.path.abspath(sys.argv[2])
GCASE = os.path.abspath(sys.argv[3])
HERE = os.path.dirname(os.path.abspath(__file__))
dtext = re.sub(r"//[^\n]*", "", open(DCASE).read())
gtext = re.sub(r"//[^\n]*", "", open(GCASE).read())

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def near(got, want, rtol):
    return rel(got, want) <= rtol


def block(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start + 1:j]
    raise SystemExit(f"unterminated {key}{{}} block")


def entries(src):
    out = {}
    for m in re.finditer(r"(\w+)\s*\{", src):
        start = m.end() - 1
        depth = 0
        for j in range(start, len(src)):
            if src[j] == "{":
                depth += 1
            elif src[j] == "}":
                depth -= 1
                if depth == 0:
                    out[m.group(1)] = src[start + 1:j]
                    break
    return out


def entries_top(src):
    """named blocks, outermost only: a child block's extent is skipped."""
    out = {}
    i = 0
    while True:
        m = re.compile(r"(\w+)\s*\{").search(src, i)
        if not m:
            return out
        start = m.end() - 1
        depth = 0
        j = start
        while j < len(src):
            if src[j] == "{":
                depth += 1
            elif src[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        out[m.group(1)] = src[start + 1:j]
        i = j + 1


def num(src, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", src)
    return float(m.group(1)) if m else default


def word(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(\w+)\s*;", src)
    return m.group(1) if m else None


def without_block(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block to strip")
    start = m.start()
    depth = 0
    j = m.end() - 1
    while j < len(src):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    return src[:start] + src[j + 1:]


def run_exe(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        proc = subprocess.run([STEADY, path], capture_output=True,
                              text=True, cwd=tmp)
    return proc


def run_full(body):
    proc = run_exe(body)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n"
                         f"{proc.stdout}\n{proc.stderr}")
    if "steady solved" not in proc.stdout or "(converged)" not in proc.stdout:
        raise SystemExit("solver did not report convergence")
    nodes, pipes, lines, dopt, glift = {}, {}, [], [], []
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        lines.append(" ".join(toks))
        if toks[0] == "node":
            nodes[toks[1]] = float(toks[3])
        elif toks[0] == "pipe":
            d = {}
            i = 2
            while i < len(toks):
                t = toks[i]
                try:
                    float(t)
                    i += 1
                    continue
                except ValueError:
                    pass
                if t in ("segregated", "transition", "intermittent",
                         "distributed", "FORCED", "CHOKED", "SONIC", "FULL",
                         "NO-LOAD", "STALLED", "CLOSED", "OPEN", "IN-SPREAD"):
                    i += 1
                    continue
                d[t] = toks[i + 1]
                i += 2
            pipes[toks[1]] = d
        elif toks[0] == "diameterOpt":
            dopt.append(toks)
        elif toks[0] == "gasLift":
            glift.append(toks)
    return dict(nodes=nodes, pipes=pipes, lines=lines, dopt=dopt,
                glift=glift, raw=proc.stdout)


def refuses(label, body, want, want2=None):
    proc = run_exe(body)
    out = (proc.stdout + proc.stderr).strip()
    ok_msg = (proc.returncode != 0 and want in out
              and (want2 is None or want2 in out))
    check(f"refused: {label}", ok_msg,
          out.splitlines()[-1][:130] if out else "no message")


# ===========================================================================
# A. diameterOptimization{}
# ===========================================================================
dopt = block(dtext, "diameterOptimization")
DPIPES = re.findall(r"\w+", re.search(r"pipes\s*\(([^)]*)\)", dopt).group(1))
DCANDS = [float(x) for x in
          re.search(r"candidates\s*\(([^)]*)\)", dopt).group(1).split()]
DCON = block(dopt, "constraints")
DCNODE = word(entries_top(DCON)["mLimit"], "node")
DCMAX = num(entries_top(DCON)["mLimit"], "maxPressure")
DREPORT = int(num(dopt, "reportTop", 3))
DCPM = num(dopt, "costPerMeter", 0.0)
DL = {p: num(entries(block(dtext, "pipes"))[p], "length") for p in DPIPES}

dstripped = without_block(dtext, "diameterOptimization")


def set_dia(body, name, d):
    m = re.search(re.escape(name) + r" \{[^}]*\}", body)
    blk = re.sub(r"diameter [0-9.eE+-]+;", f"diameter {d};", m.group(0))
    return body[:m.start()] + blk + body[m.end():]


# enumerate exactly like the app: pipe 0 is the slowest digit
combos = []
for idx in range(len(DCANDS) ** len(DPIPES)):
    dig, rest = [], idx
    for j in reversed(range(len(DPIPES))):
        dig.insert(0, rest % len(DCANDS))
        rest //= len(DCANDS)
    combos.append(dig)

table = []  # (order, diasc, pM, sumdp)
for order, dig in enumerate(combos):
    body = dstripped
    for j, name in enumerate(DPIPES):
        body = set_dia(body, name, DCANDS[dig[j]])
    r = run_full(body)
    pm = r["nodes"][DCNODE]
    sdp = sum(float(r["pipes"][p]["dp"]) for p in r["pipes"])
    table.append((order, dig, pm, sdp))

dcl = run_full(dtext)


def cost_of(dig, cpm=DCPM, tab=None):
    if tab is not None:
        return sum(DL[p] * tab[dig[j]] for j, p in enumerate(DPIPES))
    return sum(cpm * DL[p] * (DCANDS[dig[j]] / DCANDS[0]) ** 2
               for j, p in enumerate(DPIPES))


def rank(rows):
    return sorted(rows, key=lambda r: (r[0], r[1]))


feas = [(r, cost_of(r[1])) for r in table if r[2] <= DCMAX]
bestPy = rank([(r[3], r[0], r[1]) for r, _ in feas])[0]
nFeas, nComb = len(feas), len(table)

summary = next(t for t in dcl["dopt"] if t[1] == "combos")
check("diameter: the enumeration and feasibility counts",
      int(summary[2]) == nComb and int(summary[4]) == nFeas,
      f"{summary[2]}/{summary[4]} vs {nComb}/{nFeas}")
best = next(t for t in dcl["dopt"] if t[1] == "best")
wantBest = [x for j, p in enumerate(DPIPES)
            for x in (p, str(DCANDS[bestPy[2][j]]))]
check("diameter: the best combination is the re-picked tie winner",
      best[2:6] == wantBest, " ".join(best))
check("diameter: the best objective is the re-derived total drop",
      near(float(best[-3]), bestPy[0], 1e-9), f"{best[-3]}")
check("diameter: the tie resolves to the earlier enumeration",
      bestPy[1] == 1, f"order {bestPy[1]}")
check("diameter: the best cost line matches the cost model",
      near(float(best[-1]), cost_of(bestPy[2]), 1e-9), f"{best[-1]}")
runners = [t for t in dcl["dopt"] if t[1] == "runner"]
ranked = rank([(r[3], r[0], r[1]) for r, _ in feas])
check("diameter: the runner count is reportTop-1 over the feasible set",
      len(runners) == min(DREPORT, nFeas) - 1, str(len(runners)))
okRun = all(runners[k][2] == str(k + 2)
            and runners[k][3:3 + 2 * len(DPIPES)]
            == [x for j, p in enumerate(DPIPES)
                for x in (p, str(DCANDS[ranked[k + 1][2][j]]))]
            and near(float(runners[k][-3]), ranked[k + 1][0], 1e-9)
            for k in range(len(runners)))
check("diameter: every runner line matches the re-derived ranking", okRun)

# byte-identity
baseTwin = run_full(open(os.path.join(HERE, "cases", "steadyOptimize",
                                      "diameterBase.dict")).read())
kept = [l for l in dcl["lines"] if not l.startswith("diameterOpt")]
check("diameter: stripping the block leaves the base bytes exact",
      kept == baseTwin["lines"], f"{len(kept)} vs {len(baseTwin['lines'])}")
check("diameter: the stripped twin prints no diameterOpt token",
      "diameterOpt" not in baseTwin["raw"])
sn = run_full(open(os.path.join(HERE, "cases", "steadyNetwork",
                                "case.dict")).read())
check("the standing steadyNetwork case carries no optimization token",
      "diameterOpt" not in sn["raw"] and "gasLift" not in sn["raw"])

# --- variants: each re-derives its own best from the same subprocess table
def variant_best(keyfun):
    return rank([(keyfun(r), r[0], r[1]) for r in table])[0]


def dia_tokens(b):
    return [x for j, p in enumerate(DPIPES) for x in (p, str(DCANDS[b[2][j]]))]


v2 = dtext.replace(
    "constraints\n    {\n        mLimit { node M; maxPressure 1.69e6; }\n    }",
    "")
r2 = run_full(v2)
b2 = next(t for t in r2["dopt"] if t[1] == "best")
k2 = variant_best(lambda r: r[3])
check("variant no-constraints: the fattest pair minimizes the total drop",
      b2[2:6] == dia_tokens(k2) and near(float(b2[-3]), k2[0], 1e-9),
      " ".join(b2))

v3 = v2.replace("objective minTotalDp;",
               "objective nodePressure;\n objectiveNode M;") \
       .replace("objectiveSense min;", "objectiveSense max;")
r3 = run_full(v3)
b3 = next(t for t in r3["dopt"] if t[1] == "best")
k3 = variant_best(lambda r: -r[2])
check("variant nodePressure max: the fattest pair lifts the manifold",
      b3[2:6] == dia_tokens(k3) and near(float(b3[-3]), -k3[0], 1e-9),
      " ".join(b3))

v4 = v2.replace("objective minTotalDp;",
               "objective nodePressure;\n objectiveNode M;")
r4 = run_full(v4)
b4 = next(t for t in r4["dopt"] if t[1] == "best")
k4 = variant_best(lambda r: r[2])
check("variant nodePressure min: the leanest pair sinks the manifold",
      b4[2:6] == dia_tokens(k4) and near(float(b4[-3]), k4[0], 1e-9),
      " ".join(b4))

v5 = v2.replace("objective minTotalDp;", "objective minCost;")
r5 = run_full(v5)
b5 = next(t for t in r5["dopt"] if t[1] == "best")
k5 = variant_best(lambda r: cost_of(r[1]))
check("variant minCost per-meter: the smallest bores are the cheapest",
      b5[2:6] == dia_tokens(k5) and near(float(b5[-3]), k5[0], 1e-9)
      and near(float(b5[-1]), k5[0], 1e-9), " ".join(b5))

v6 = v5.replace("costPerMeter 500;", "costTable ( 1 4 9 );")
r6 = run_full(v6)
b6 = next(t for t in r6["dopt"] if t[1] == "best")
k6 = variant_best(lambda r: cost_of(r[1], tab=[1.0, 4.0, 9.0]))
check("variant minCost table: the per-candidate table prices the pair",
      b6[2:6] == dia_tokens(k6) and near(float(b6[-1]), k6[0], 1e-9),
      " ".join(b6))

# --- refusals
refuses("an unknown diameterOptimization key",
        dtext.replace("reportTop 3;", "reportTopX 3;"),
        "unknown key 'reportTopX'")
refuses("a pipe that is not in the network",
        dtext.replace("pipes ( P1 P2 );", "pipes ( P1 P9 );"),
        "not in the network")
refuses("descending candidates",
        dtext.replace("candidates ( 0.05 0.06 0.075 );",
                      "candidates ( 0.06 0.05 0.075 );"),
        "strictly ascending")
refuses("duplicate candidates",
        dtext.replace("candidates ( 0.05 0.06 0.075 );",
                      "candidates ( 0.05 0.05 0.075 );"),
        "strictly ascending")
refuses("a non-positive candidate",
        dtext.replace("candidates ( 0.05 0.06 0.075 );",
                      "candidates ( 0.05 -0.06 0.075 );"),
        "must be positive")
refuses("an empty candidate list",
        dtext.replace("candidates ( 0.05 0.06 0.075 );", "candidates ( );"),
        "candidates")
refuses("an enumeration above maxCombinations",
        dtext.replace("maxCombinations 64;", "maxCombinations 4;"),
        "exceeds maxCombinations")
refuses("an unknown objective word",
        dtext.replace("objective minTotalDp;", "objective minFriction;"),
        "unknown objective")
refuses("nodePressure without objectiveNode",
        dtext.replace("objective minTotalDp;", "objective nodePressure;"),
        "needs objectiveNode")
refuses("an unknown objectiveSense",
        dtext.replace("objectiveSense min;", "objectiveSense median;"),
        "objectiveSense must be min or max")
refuses("both cost models",
        dtext.replace("costPerMeter 500;",
                      "costPerMeter 500;\n costTable ( 1 4 9 );"),
        "cannot declare both")
refuses("minCost with no cost model",
        v2.replace("objective minTotalDp;", "objective minCost;")
        .replace("costPerMeter 500;\n    ", ""),
        "needs costPerMeter or costTable")
refuses("a costTable of the wrong length",
        v5.replace("costPerMeter 500;", "costTable ( 1 4 );"),
        "one entry per candidate")
refuses("an unknown constraint key",
        dtext.replace("maxPressure 1.69e6;", "maxPressur 1.69e6;"),
        "unknown key 'maxPressur'")
refuses("a constraint node outside the network",
        dtext.replace("mLimit { node M;", "mLimit { node M9;"),
        "not in the network")
refuses("a constraint with neither bound",
        dtext.replace("mLimit { node M; maxPressure 1.69e6; }",
                      "mLimit { node M; }"),
        "declares neither")
refuses("a constraint with minPressure above maxPressure",
        dtext.replace("mLimit { node M; maxPressure 1.69e6; }",
                      "mLimit { node M; minPressure 2.0e6; "
                      "maxPressure 1.69e6; }"),
        "minPressure above maxPressure")

# ===========================================================================
# B. gasLiftAllocation{}
# ===========================================================================
PSI = 6894.757293168361
SCF = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
G = 9.80665
RHO_WATER_STD = 999.014
LBFT3 = 16.0184634

fl = block(gtext, "blackoil")
GG = num(fl, "gasGravity", 0.65)
RHO_STD = num(fl, "oilDensity", 850.0)
GOR = num(fl, "gor", 0.0)
WC = num(fl, "waterCut", 0.0)
T_FLUID = num(fl, "temperature", 300.0)
SIG = num(fl, "surfaceTension", 0.025)
WATER_DENSITY = num(fl, "waterDensity", 1000.0)
WATER_MU = num(fl, "waterViscosity", 1e-3)

AB = (0.98, 0.0, 0.845, 1.065)
BB = (0.4846, 0.0, 0.5351, 0.5824)
CB = (0.0868, 0.0, 0.0173, 0.0609)


def boundaries(lam):
    return (316.0 * lam ** 0.302, 0.0009252 * lam ** -2.4684,
            0.1 * lam ** -1.4516, 0.5 * lam ** -6.738)


def regime_of(lam, nfr):
    l1, l2, l3, l4 = boundaries(lam)
    if (lam < 0.01 and nfr < l1) or (lam >= 0.01 and nfr < l2):
        return 0
    if lam >= 0.01 and l2 <= nfr <= l3:
        return 1
    if (0.01 <= lam < 0.4 and l3 < nfr <= l1) or (lam >= 0.4 and l3 < nfr <= l4):
        return 2
    return 3


def hl0_of(lam, nfr, r):
    return max(AB[r] * lam ** BB[r] * nfr ** -CB[r], lam)


def hl_horizontal(lam, nfr, r):
    if r == 1:
        l1, l2, l3, l4 = boundaries(lam)
        w = (l3 - nfr) / (l3 - l2)
        return w * hl0_of(lam, nfr, 0) + (1.0 - w) * lam
    return hl0_of(lam, nfr, r)


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg):
    hl = hl_horizontal(lam, nfr, r)
    if abs(angle_deg) > 1e-12 and lam > 1e-9 and r != 3:
        coeff = ((4.70, -0.3692, 0.1244, -0.5056) if angle_deg < 0.0
                 else ((0.011, -3.768, 3.539, -1.614) if r == 0
                       else (2.96, 0.305, -0.4473, 0.0978)))
        d, e, f, h = coeff
        nlv = vsl * math.sqrt(math.sqrt(rho_l / (G * SIG)))
        c = max((1.0 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** h), 0.0)
        s1 = math.sin(1.8 * math.radians(angle_deg))
        hl = hl * (1.0 + c * (s1 - s1 ** 3 / 3.0))
    return min(max(hl, 0.0), 1.0)


def colebrook(re, ed):
    if re < 2040.0:
        return 64.0 / max(re, 1e-12)
    ed = max(ed, 0.0)
    y = max(1.14 - 2.0 * math.log10(ed + 21.25 / re ** 0.9), 1.0)
    for _ in range(100):
        yn = -2.0 * math.log10(ed / 3.7 + 2.51 * y / re)
        err = abs(yn - y) / yn
        y = yn
        if err < 1e-14:
            break
    return 1.0 / (y * y)


class BlackOil:
    def __init__(self, gg, rho_std, gor, wc, temp, water_density, water_mu):
        self.gg, self.rho_std, self.gor, self.wc = gg, rho_std, gor, wc
        self.T, self.wd, self.wmu = temp, water_density, water_mu
        self.api = 141.5 / (rho_std / RHO_WATER_STD) - 131.5
        self.tf = temp * 1.8 - 459.67
        self.tr = temp * 1.8
        self.gor_scf = gor * SCF
        self.wor = wc / (1.0 - wc) if 0.0 < wc < 0.999 else 0.0
        rho_g_std = P_STD * gg * M_AIR / (R_GAS * T_STD)
        self.mass_per_sm3 = rho_std + gor * rho_g_std + self.wor * water_density

    def pb_psia(self):
        a = 0.00091 * self.tf - 0.0125 * self.api
        return max(18.2 * ((self.gor_scf / self.gg) ** 0.83 * 10 ** a - 1.4),
                   0.0)

    def rs_scf(self, p_psia):
        c1, c2, c3 = ((0.0362, 1.0937, 25.724) if self.api <= 30
                      else (0.0178, 1.187, 23.931))
        return (c1 * self.gg * max(p_psia, 0.0) ** c2
                * math.exp(c3 * self.api / (self.tf + 459.67)))

    def bo_of(self, p_psia, rs, pb):
        api, gg, tf = self.api, self.gg, self.tf
        c1, c2, c3 = ((4.670e-4, 1.100e-5, 1.337e-9) if api > 30
                      else (4.677e-4, 1.751e-5, -1.811e-8))
        tt = (tf - 60.0) * api / gg
        bob = 1 + c1 * rs + c2 * tt + c3 * rs * tt
        if p_psia <= pb:
            return bob
        co = ((-1433 + 5 * rs + 17.2 * tf - 1180 * gg + 12.61 * api)
              / (1e5 * p_psia))
        return bob * math.exp(co * (pb - p_psia))

    def z_dak(self, p_psia):
        gg, tr = self.gg, self.tr
        if gg < 0.75:
            tpc = 168 + 325 * gg - 12.5 * gg * gg
            ppc = 677 + 15 * gg - 37.5 * gg * gg
        else:
            tpc = 187 + 330 * gg - 71.5 * gg * gg
            ppc = 706 + 51 * gg - 11.1 * gg * gg
        A = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361,
             0.1844, 0.1056, 0.6134, 0.7210]
        pr, trr = max(p_psia, 1e-6) / ppc, tr / tpc
        t1 = A[0] + A[1] / trr + A[2] / trr ** 3 + A[3] / trr ** 4 + A[4] / trr ** 5
        t2 = A[5] + A[6] / trr + A[7] / trr ** 2

        def zofrho(rho):
            r2 = rho * rho
            return (1 + t1 * rho + t2 * r2
                    - A[8] * (A[6] / trr + A[7] / trr ** 2) * rho ** 5
                    + A[9] * (1 + A[10] * r2) * r2 / trr ** 3
                    * math.exp(-A[10] * r2))

        def res(z):
            return zofrho(0.27 * pr / (z * trr)) - z

        lo, hi = 0.05, 3.0
        if res(lo) <= 0:
            return lo
        if res(hi) >= 0:
            return hi
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if res(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-13 * hi:
                break
        return 0.5 * (lo + hi)

    def state(self, p_pa):
        ppsia = max(p_pa, 1.0) / PSI
        pb = self.pb_psia()
        if p_pa <= pb * PSI:
            rs = min(self.rs_scf(ppsia) / SCF, self.gor)
            b = self.bo_of(ppsia, rs * SCF, pb)
            fg = self.gor - rs
        else:
            rs, fg = self.gor, 0.0
            b = self.bo_of(ppsia, self.gor_scf, pb)
        z = self.z_dak(ppsia)
        bg = (P_STD / p_pa) * (self.T / T_STD) * z
        rhog = p_pa * self.gg * M_AIR / (z * R_GAS * self.T)
        rho_oil = ((350 * (self.rho_std / RHO_WATER_STD)
                    + 0.0764 * self.gg * rs * SCF) / (SCF * b)) * LBFT3
        x = 10 ** (3.0324 - 0.02023 * self.api) * self.tf ** -1.163
        mu_dead = 10 ** x - 1
        mu_oil = (10.715 * (rs * SCF + 100) ** -0.515) * mu_dead ** (
            5.44 * (rs * SCF + 150) ** -0.338) * 1e-3
        mw = self.gg * 28.9647
        k = (9.4 + 0.02 * mw) * self.tr ** 1.5 / (209 + 19 * mw + self.tr)
        xx = 3.5 + 986 / self.tr + 0.01 * mw
        yy = 2.4 - 0.2 * xx
        mu_gas = 1e-4 * k * math.exp(xx * (rhog / 1000) ** yy) * 1e-3
        v_t = b + self.wor
        fo = b / v_t
        rho_l = (rho_oil * b + self.wd * self.wor) / v_t
        mu_l = math.exp(fo * math.log(mu_oil) + (1.0 - fo) * math.log(self.wmu))
        return dict(bo=b, fg=fg, bg=bg, rhog=rhog, rhol=rho_l, mul=mu_l,
                    mug=mu_gas, wor=self.wor, rho_oil=rho_oil)

    def march_profile(self, q, pin, L, D, rough, sintheta=0.0, nseg=24):
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sintheta))))
        p = pin
        segs = []
        ok = True
        for _ in range(nseg):
            st = self.state(p)
            vsl = q * (st["bo"] + st["wor"]) / area
            vsg = q * st["fg"] * st["bg"] / area
            vm = vsl + vsg
            if vm <= 1e-9:
                # the kernel's static column: no-slip mixture density,
                # gravity-only gradient, regime 3
                vT = st["bo"] + st["wor"]
                rhoM = (st["rho_oil"] * st["bo"]
                        + WATER_DENSITY * st["wor"]) / vT
                pOld = p
                p -= G * sintheta * rhoM * ds
                if not p > p_floor:
                    p = p_floor
                    ok = False
                segs.append(dict(p=pOld, hl=st["bo"] / vT, lam=st["bo"] / vT,
                                  vm=0.0, rhom=rhoM, r=3, pNew=p))
                if not ok:
                    break
                continue
            lam = vsl / vm
            nfr = vm * vm / (G * D)
            r = regime_of(lam, nfr)
            hl = hl_priced(lam, nfr, r, vsl, st["rhol"], angle_deg)
            rhom = st["rhol"] * lam + st["rhog"] * (1 - lam)
            rhos = st["rhol"] * hl + st["rhog"] * (1 - hl)
            mum = st["mul"] * lam + st["mug"] * (1 - lam)
            re = rhom * vm * D / mum
            fns = colebrook(re, rough / D)
            if lam < 1e-9 or hl < 1e-9:
                sm = 0.0
            else:
                yy = lam / (hl * hl)
                ly = math.log(yy)
                if 1.0 < yy < 1.2:
                    sm = math.log(2.2 * yy - 1.2)
                else:
                    sm = (ly / (ly * (ly * (0.01853 * ly * ly - 0.8725)
                                      + 3.182) - 0.0523))
                sm = min(sm, 7.0)
            ftp = fns * math.exp(sm)
            ek = vsg * vm * rhos / max(p, 1000.0)
            ele = G * math.sin(math.radians(angle_deg)) * rhos
            dpds = (ele + ftp * rhom * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            pOld = p
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
            segs.append(dict(p=pOld, hl=hl, lam=lam, vm=vm, rhom=rhom,
                              r=r, pNew=p))
            if not ok:
                break
        return segs, p, ok


def inlet_for_rate(fluid, well, rate, pwh):
    """wellboreInletForRate restated: expanding bracket + Illinois secant."""
    L, D, rough, elev = well
    sinth = elev / L

    def g(p):
        _, pe, ok = fluid.march_profile(rate, p, L, D, rough, sinth)
        if ok:
            return pe - pwh
        return max(0.02 * p, 1000.0) - pwh

    pLo = pwh
    fLo = g(pLo)
    if fLo > 0.0:
        raise SystemExit("wellbore gains pressure toward the wellhead")
    if fLo > -1e-12 * pwh:
        return pLo
    pHi = pLo + max(1e4, 0.05 * pwh)
    fHi = g(pHi)
    expansions = 0
    while fHi < 0.0:
        pLo, fLo = pHi, fHi
        pHi *= 2.0
        fHi = g(pHi)
        expansions += 1
        if expansions > 100:
            raise SystemExit("the wellbore march cannot reach the wellhead")

    p = pLo + (pHi - pLo) * fLo / (fLo - fHi)
    for _ in range(100):
        fp = g(p)
        if abs(fp) < 1e-12 * pwh:
            break
        if fp > 0.0:
            pHi, fHi = p, fp * 0.5
        else:
            pLo, fLo = p, fp * 0.5
        pNew = pLo + (pHi - pLo) * fLo / (fLo - fHi)
        p = pNew if (pLo < pNew < pHi) else 0.5 * (pLo + pHi)
        if abs(pHi - pLo) < 1e-12 * pHi:
            break
    return p


def ipr_pw(ipr, q):
    form, pres = ipr["form"], ipr["pres"]
    if form == "linear":
        return pres - q / ipr["pi"]
    if form == "vogel":
        qq = min(max(q, 0.0), ipr["qmax"])
        x = (-0.2 + math.sqrt(3.24 - 3.2 * qq / ipr["qmax"])) / 1.6
        return pres * x
    if form == "jones":
        q2 = min(max(q, 0.0), ipr["c"] * (pres * pres) ** ipr["n"])
        return math.sqrt(max(pres * pres - (q2 / ipr["c"]) ** (1.0 / ipr["n"]),
                             0.0))
    raise SystemExit(f"ipr form {form} not restated")


def solve_nodal(fluid, ipr, well, pwh, tol=1e-9, maxIter=100):
    """solveNodalAnalysis restated: the largest IPR x VLP root, bracketed
    from the qMax end downward, Illinois secant with half-weighted updates."""
    pres = ipr["pres"]
    if ipr["form"] == "linear":
        qMax = ipr["pi"] * pres
    elif ipr["form"] == "vogel":
        qMax = ipr["qmax"]
    else:
        qMax = ipr["c"] * (pres * pres) ** ipr["n"]

    def f(q):
        return inlet_for_rate(fluid, well, q, pwh) - ipr_pw(ipr, q)

    qHi = qMax
    fHi = f(qHi)
    qLo = fLo = 0.0
    bracketed = False
    qProbe = qHi
    for _ in range(60):
        if bracketed:
            break
        qProbe *= 0.5
        fk = f(qProbe)
        if fk <= 0.0:
            qLo, fLo, bracketed = qProbe, fk, True
    if not bracketed and f(0.0) <= 0.0:
        qLo, fLo, bracketed = 0.0, f(0.0), True
    if not bracketed:
        return 0.0

    pScale = max(pres, pwh)
    q = qLo + (qHi - qLo) * fLo / (fLo - fHi)
    for _ in range(maxIter):
        fq = f(q)
        if abs(fq) < 1e-7 * pScale:
            break
        if fq > 0.0:
            qHi, fHi = q, fq * 0.5
        else:
            qLo, fLo = q, fq * 0.5
        qNew = qLo + (qHi - qLo) * fLo / (fLo - fHi)
        q = qNew if (qLo < qNew < qHi) else 0.5 * (qLo + qHi)
        if abs(qHi - qLo) < tol * max(qHi, 1e-30):
            break
    return q


def rate_at_gas(fluid_params, ipr, well, pwh, g, q0):
    """the gas/GOR fixed point: gor' = gor + g/q."""
    if g <= 0.0:
        return q0
    q = q0
    for _ in range(60):
        fl2 = BlackOil(GG, RHO_STD, GOR + g / q, WC, T_FLUID,
                       WATER_DENSITY, WATER_MU)
        r = solve_nodal(fl2, ipr, well, pwh)
        if abs(r - q) <= 1e-10 * r:
            return r
        q = r
    raise SystemExit("the gas/GOR fixed point did not converge")


gb = block(gtext, "gasLiftAllocation")
BUDGET = num(gb, "budget")
RES = int(num(gb, "resolution"))
GWELLS = entries_top(block(gb, "wells"))
STEP = BUDGET / RES

gcl = run_full(gtext)
gstripped = run_full(without_block(gtext, "gasLiftAllocation"))

# the base node pressures are the wellhead back-pressures
wells = []
for name, body in GWELLS.items():
    node = word(body, "node")
    iprb = block(body, "ipr")
    form = word(iprb, "form")
    ipr = dict(form=form, pres=num(iprb, "reservoirPressure"),
               pi=num(iprb, "productivityIndex"), qmax=num(iprb, "maxRate"),
               c=num(iprb, "deliverability"), n=num(iprb, "exponent", 1.0))
    wb = block(body, "well")
    well = (num(wb, "length"), num(wb, "diameter"), num(wb, "roughness"),
            num(wb, "elevation"))
    pwh = gstripped["nodes"][node]
    q0 = solve_nodal(BlackOil(GG, RHO_STD, GOR, WC, T_FLUID,
                              WATER_DENSITY, WATER_MU), ipr, well, pwh)
    curve = [q0]
    for k in range(1, RES + 1):
        curve.append(rate_at_gas(None, ipr, well, pwh, k * STEP, q0))
    wells.append(dict(name=name, q=curve))

NEG = float("-inf")
dp = [NEG] * (RES + 1)
dp[0] = 0.0
parent = [[0] * (RES + 1) for _ in wells]
for i, w in enumerate(wells):
    ndp = [NEG] * (RES + 1)
    for b in range(RES + 1):
        best, bk = NEG, 0
        for k in range(b + 1):
            if dp[b - k] == NEG:
                continue
            v = dp[b - k] + w["q"][k]
            if v > best:
                best, bk = v, k
        ndp[b], parent[i][b] = best, bk
    dp = ndp
alloc = [0] * len(wells)
b = RES
for i in range(len(wells) - 1, -1, -1):
    alloc[i] = parent[i][b]
    b -= alloc[i]
totalPy = sum(wells[i]["q"][alloc[i]] for i in range(len(wells)))
basePy = sum(w["q"][0] for w in wells)


def marg_py(q, k):
    if k == 0:
        return (q[1] - q[0]) / STEP
    if k == RES:
        return (q[RES] - q[RES - 1]) / STEP
    return (q[k + 1] - q[k - 1]) / (2.0 * STEP)


wl = next(t for t in gcl["glift"] if t[1] == "wells")
check("gasLift: the wells/budget/step header",
      int(wl[2]) == len(wells) and near(float(wl[4]), BUDGET, 1e-12)
      and near(float(wl[6]), STEP, 1e-12), " ".join(wl))
wlines = {t[2]: t for t in gcl["glift"] if t[1] == "well"}
check("gasLift: one line per declared well",
      set(wlines) == {w["name"] for w in wells})
okQ = okG = okGain = okMarg = True
gasSum = 0.0
for i, w in enumerate(wells):
    t = wlines[w["name"]]
    d = {t[j]: float(t[j + 1]) for j in range(3, len(t), 2)}
    gasSum += d["gas"]
    okQ &= near(d["q"], w["q"][alloc[i]], 1e-6)
    okG &= abs(d["gas"] - alloc[i] * STEP) <= 1e-12 * STEP
    okGain &= abs(d["gain"] - (w["q"][alloc[i]] - w["q"][0])) < 1e-9
    okMarg &= rel(d["marginal"], marg_py(w["q"], alloc[i])) <= 1e-3
check("gasLift: every reported rate is the re-derived nodal fixed point",
      okQ)
check("gasLift: every reported gas is the DP allocation on the grid", okG)
check("gasLift: the gains are q - q(0)", okGain)
check("gasLift: the marginals are the re-derived finite differences", okMarg)
check("gasLift: the budget is spent exactly",
      abs(gasSum - BUDGET) <= 1e-12 * BUDGET, f"{gasSum}")
tl = next(t for t in gcl["glift"] if t[1] == "total")
check("gasLift: the total rate matches the re-derived DP optimum",
      near(float(tl[3]), totalPy, 1e-6)
      and abs(float(tl[5]) - (totalPy - basePy)) < 1e-9, " ".join(tl))
ss = next(t for t in gcl["glift"] if t[1] == "slopeSpread")
margList = [marg_py(w["q"], alloc[i]) for i, w in enumerate(wells)]
check("gasLift: the slope spread is max - min marginal",
      rel(float(ss[2]), max(margList) - min(margList)) <= 1e-2, ss[2])
check("gasLift: the split is non-trivial (both wells get gas, neither gets all)",
      0 < alloc[0] < RES and 0 < alloc[1] < RES
      and alloc[0] + alloc[1] == RES, f"{alloc}")
okEx = True
for i in range(len(wells)):
    for j in range(len(wells)):
        if i == j or alloc[i] == 0 or alloc[j] == RES:
            continue
        moved = totalPy - wells[i]["q"][alloc[i]] + wells[i]["q"][alloc[i] - 1] \
                - wells[j]["q"][alloc[j]] + wells[j]["q"][alloc[j] + 1]
        okEx &= moved <= totalPy + 1e-9
check("gasLift: no single-step transfer improves the allocation", okEx)

# byte-identity for the allocation block
kept = [l for l in gcl["lines"] if not l.startswith("gasLift")]
check("gasLift: stripping the block leaves the base bytes exact",
      kept == gstripped["lines"], f"{len(kept)} vs {len(gstripped['lines'])}")
check("gasLift: the stripped twin prints no gasLift token",
      "gasLift" not in gstripped["raw"])

# --- refusals
refuses("a non-positive budget", gtext.replace("budget 0.12;", "budget 0;"),
        "budget must be positive")
refuses("a negative budget", gtext.replace("budget 0.12;", "budget -0.12;"),
        "budget must be positive")
refuses("a zero resolution", gtext.replace("resolution 4;", "resolution 0;"),
        "integer in [1, 64]")
refuses("a resolution above 64",
        gtext.replace("resolution 4;", "resolution 65;"),
        "integer in [1, 64]")
refuses("a non-integer resolution",
        gtext.replace("resolution 4;", "resolution 2.5;"),
        "integer in [1, 64]")
refuses("an empty wells block",
        re.sub(r"wells\s*\{.*?\n    \}", "wells { }", gtext, flags=re.S),
        "declares no wells")
refuses("a well node outside the network",
        gtext.replace("            node W1;", "            node W9;"),
        "not in the network")
refuses("an unknown gasLiftAllocation key",
        gtext.replace("resolution 4;", "resolutions 4;"),
        "unknown key 'resolutions'")
refuses("an unknown well key",
        gtext.replace("            node W1;",
                      "            node W1;\n            diameter 0.1;"),
        "unknown key 'diameter'")
refuses("a well without an ipr block",
        re.sub(r"ipr \{ form vogel; reservoirPressure 2.3e7; maxRate 0.003; \}\n",
               "", gtext),
        "missing required key 'ipr'")
refuses("gas lift on a constant-property network",
        """
fluid { rho 850; mu 0.001; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.1; roughness 4.6e-5;
    elevation 0; } }
sources { s { node W; rate 0.001; } }
solve { maxIter 100; tol 1e-10; }
gasLiftAllocation { budget 0.1; resolution 2;
    wells { W1 { node W;
        ipr { form vogel; reservoirPressure 2.0e7; maxRate 0.001; }
        well { length 1000; diameter 0.0762; roughness 5e-5; elevation 900; }
    } } }
""", "needs a black-oil")
refuses("a well that does not flow at zero gas",
        gtext.replace("reservoirPressure 2.3e7;", "reservoirPressure 1.0e7;"),
        "zero rate at zero gas")

if fails:
    raise SystemExit(f"steady optimization verification FAILED: {fails}")
print("the diameter search and the gas-lift allocation are pure app-layer "
      "reports: every combination, objective, tie order, cost, nodal fixed "
      "point, DP allocation and marginal slope reproduces the independent "
      "re-derivation, the undeclared byte stream is exact, and the two "
      "declaration grammars refuse every inert or impossible shape")
