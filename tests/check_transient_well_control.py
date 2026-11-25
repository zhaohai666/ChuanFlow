#!/usr/bin/env python3
"""5.185 (W-class) verification: the well-control reservoir inflow of the 1D
transient kernel.

usage: check_transient_well_control.py <exe> <case.dict> <tmpdir>

The checker never trusts a number the kernel prints for the reservoir. It reads
the declarations out of the case text, re-implements the two laws locally
    q = J * max(p_res - p_node, 0)                   [PI inflow]
    (V/K) dp/dt = q - Q_pipe                          [node continuity]
and cross-checks the marched trajectory against them:

  * every frame of the kick case: q against J (p_res - p_BH) from the printed
    node pressure; the gas share q*fg; the signed underbalance and its flag;
    the two static well-control numbers (rho_kill = p_res/(g depth) and
    SIDPP = p_res - rho_mud g depth);
  * every STEP: the implicit-Euler node row (V/K)(p_n - p_{n-1})/dt = q_n -
    Q_n, with q_n re-read from the printed pressure, and the cumulative volume
    increment dt*q_n (the account commitStep keeps);
  * after the valve closes the exact first-order geometric approach to p_res,
    (C/dt + J) p' = (C/dt) p + J p_res - Q', so the decay is a closed form,
    not a fitted curve;
  * a continuous UNDERBALANCED arm: the same laws without a shut-in, ending
    where the PI inflow feeds the pipe rate (q == Q to 1e-6);
  * an OVERBALANCED arm: q must be exactly 0 (never a negative "production"
    through the reservoir) and the whole run must be byte-identical to the
    same deck with reservoir{} deleted - the one-sided law charges nothing
    when it is inactive;
  * seven loud refusals (fixed node, sources overlap, two reservoirs on one
    node, gasFraction 1, zero depth, zero J, zero pressure).
"""

import math
import re
import subprocess
import sys
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

EXE, CASE, TMP = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
CASE_DIR = CASE.parent
G = 9.80665

checks = 0
failures = []


def ok(name, cond, detail=""):
    global checks
    checks += 1
    if cond:
        print("ok:  ", name)
    else:
        msg = name + (" | " + detail if detail else "")
        print("FAIL:", msg)
        failures.append(msg)


def close(name, got, want, rtol=1e-9, atol=0.0):
    scale = max(abs(got), abs(want), atol)
    ok(name, abs(got - want) <= rtol * scale + atol,
       f"got {got:.12g} want {want:.12g} diff {abs(got-want):.3g}")


# The app prints 12 significant digits, so a quantity of size X carries a
# quantisation of about X * 1e-12; a difference of two such values is the
# resolution floor of the report, not a physics disagreement.
QUANT_P = 1e-5          # Pa resolution at p ~ 1e7-1e8
QUANT_VOL = 2e-11       # m3 resolution at vol ~ 2 m3
QUANT_RHO = 2e-12       # relative resolution of a 12-digit number


def run(path):
    out = subprocess.run([EXE, str(path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


# ---- declaration parsing (from the case text, never from the report) --------
def scal(text, key, block=None):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", text)
    if not m:
        raise KeyError(key)
    return float(m.group(1))


def parse_case(path):
    text = Path(path).read_text()
    d = {}
    d["rho"] = scal(text, "rho")
    d["volume"] = scal(text, "volume")
    d["bulk"] = scal(text, "bulkModulus")
    nodes_blk = re.search(r"^nodes\s*\{(.*?)^\}", text, re.S | re.M).group(1)
    d["p0"] = scal(nodes_blk, "pressure")       # the free node's initial value
    d["dt"] = scal(text, "dt")
    rc = re.search(r"^reservoir\s*\{(.*?)^\}", text, re.S | re.M)
    d["res"] = rc.group(1) if rc else ""
    if d["res"]:
        body = d["res"]
        d["pres"] = float(re.search(r"pressure\s+([-+0-9.eE]+)", body).group(1))
        d["J"] = float(re.search(r"productivityIndex\s+([-+0-9.eE]+)",
                                 body).group(1))
        d["fg"] = float(re.search(r"gasFraction\s+([-+0-9.eE]+)",
                                  body).group(1))
        d["depth"] = float(re.search(r"depth\s+([-+0-9.eE]+)", body).group(1))
        d["node"] = re.search(r"node\s+(\w+)", body).group(1)
    return d


FRAME = re.compile(r"\[chflow\] t (\S+)")
NODE = re.compile(r"\[chflow\] node (\w+) p (\S+) alpha (\S+) hydrate (\S+)")
RES = re.compile(r"\[chflow\] reservoir (\w+) node (\w+) q (\S+) gasq (\S+) "
                 r"vol (\S+) gasvol (\S+) ub (\d+) dp (\S+) killmud (\S+) "
                 r"sidpp (\S+)")
PIPE = re.compile(r"\[chflow\] pipe (\w+) Q (\S+)")


def parse_frames(text, node):
    frames = []
    cur = None
    for line in text.splitlines():
        fm = FRAME.match(line)
        if fm:
            if cur is not None:
                frames.append(cur)
            cur = {"t": float(fm.group(1)), "p": {}, "res": None, "Q": {}}
            continue
        if cur is None:
            continue
        nm = NODE.match(line)
        if nm and nm.group(1) == node:
            cur["p"] = float(nm.group(2))
            cur["alpha"] = float(nm.group(3))
            continue
        rm = RES.match(line)
        if rm:
            cur["res"] = dict(name=rm.group(1), node=rm.group(2),
                              q=float(rm.group(3)), gasq=float(rm.group(4)),
                              vol=float(rm.group(5)), gasvol=float(rm.group(6)),
                              ub=int(rm.group(7)), dp=float(rm.group(8)),
                              killmud=float(rm.group(9)),
                              sidpp=float(rm.group(10)))
            continue
        pm = PIPE.match(line)
        if pm:
            cur["Q"][pm.group(1)] = float(pm.group(2))
    if cur is not None:
        frames.append(cur)
    return frames


# ---- the two independent laws ----------------------------------------------
def inflow(d, p):
    return d["J"] * max(d["pres"] - p, 0.0)


def continuity_residual(d, prev, cur, Q):
    C = d["volume"] / d["bulk"]
    return C * (cur["p"] - prev["p"]) / d["dt"] - (inflow(d, cur["p"]) - Q)


def main():
    # ---- 1. gas kick taken underbalanced, then shut in ----------------------
    d = parse_case(CASE)
    rc, out = run(CASE)
    ok("kick arm exits 0", rc == 0, f"rc={rc}")
    fr = parse_frames(out, d["node"])
    ok("kick arm has frames", len(fr) >= 3, f"n={len(fr)}")
    ok("kick arm has a reservoir line on every frame",
       all(f["res"] is not None for f in fr))

    def relerr(a, b, floor=0.0):
        scale = max(abs(a), abs(b), floor)
        return abs(a - b) / scale if scale > 0 else 0.0

    # Every law aggregated over the run: one judgement per closed form, with
    # the worst frame named (a 1200-frame trajectory needs no 1200 lines).
    wq = wg = wflag = wkill = wsid = wdp = 0.0
    for f in fr:
        r = f["res"]
        # q is priced from the printed dp, which is itself checked against
        # p_res - p; the pressure print resolution makes J*(p_res - p)
        # non-deterministic at 1e-12 when q is tiny, the printed dp is not.
        wq = max(wq, relerr(r["q"], d["J"] * max(r["dp"], 0.0), 1e-12))
        wg = max(wg, relerr(r["gasq"], d["fg"] * r["q"], 1e-12))
        wdp = max(wdp, abs(r["dp"] - (d["pres"] - f["p"])))
        wflag = max(wflag, 1.0 if (r["ub"] == 1) != (f["p"] < d["pres"]) else 0.0)
        wkill = max(wkill, relerr(r["killmud"],
                                  d["pres"] / (G * d["depth"])))
        wsid = max(wsid, relerr(r["sidpp"],
                                d["pres"] - d["rho"] * G * d["depth"]))
    ok("q == J dp on every frame (worst %.3g)" % wq, wq < 1e-9)
    ok("gasq == fg q on every frame (worst %.3g)" % wg, wg < 1e-9)
    ok("dp == p_res - p on every frame (worst %.3g Pa)" % wdp, wdp < 1e-2)
    ok("ub flag is p < p_res on every frame", wflag == 0.0)
    ok("kill mud == p_res/(g depth) on every frame (worst %.3g)" % wkill,
       wkill < 1e-10)
    ok("SIDPP == p_res - rho g depth on every frame (worst %.3g)" % wsid,
       wsid < 1e-10)

    # The implicit-Euler node row, every step (reportEvery 1 => one per step):
    #   (V/K)(p_n - p_{n-1})/dt = q_n - Q_n
    # and the cumulative volume is the same step's dt*q_n account.
    worst = 0.0
    worst_at = 0.0
    worst_abs = 0.0
    wvol = wgas = 0.0
    for prev, cur in zip(fr, fr[1:]):
        Q = cur["Q"].get("WELL", 0.0)
        res = continuity_residual(d, prev, cur, Q)
        scale = max(abs(inflow(d, cur["p"])), abs(Q), 1e-6)
        worst_abs = max(worst_abs, abs(res))
        # Above ~1e-3 m3/s a relative assertion carries the judgement; below
        # it the printed 12-digit pressure quantises at 1e-4 Pa (C/dt * 1e-4
        # Pa = 5e-11 m3/s of residual), so the honest floor is absolute.
        if abs(res) > 1e-9 + 1e-8 * scale:
            worst, worst_at = max(worst, abs(res) / scale), cur["t"]
        wvol = max(wvol, abs(cur["res"]["vol"] - prev["res"]["vol"]
                             - d["dt"] * inflow(d, cur["p"])))
        wgas = max(wgas, abs(cur["res"]["gasvol"] - prev["res"]["gasvol"]
                             - d["dt"] * d["fg"] * inflow(d, cur["p"])))
    ok("step continuity closes every step (largest residual %.3g, rel %.3g "
       "at t=%.2f)" % (worst_abs, worst, worst_at), worst < 1e-9,
       f"worst={worst:.3g}")
    ok("volume increment is dt q_n every step (worst %.3g m3)" % wvol,
       wvol < 2e-11)
    ok("gas volume increment is dt fg q_n every step (worst %.3g m3)" % wgas,
       wgas < 2e-11)

    # shut-in: the valve closes at t = 20 (declared), afterwards Q ~ 0 and the
    # pressure follows the exact geometric law
    C = d["volume"] / d["bulk"]
    tclose = scal(Path(CASE).read_text(), "start")
    post = [f for f in fr if f["t"] >= tclose + d["dt"] * 0.5]
    ok("shut-in leaves a sealed pipe", all(abs(f["Q"].get("WELL", 0.0)) < 1e-3
                                          for f in post),
       f"max Q = {max(abs(f['Q'].get('WELL',0.0)) for f in post):.3g}")
    worst = 0.0
    for prev, cur in zip(fr, fr[1:]):
        if cur["t"] < tclose + d["dt"] * 0.5:
            continue
        Q = cur["Q"].get("WELL", 0.0)
        want = ((C / d["dt"]) * prev["p"] + d["J"] * d["pres"] - Q) \
            / (C / d["dt"] + d["J"])
        worst = max(worst, abs(cur["p"] - want) / d["pres"])
    ok("shut-in follows the geometric approach to p_res", worst < 1e-9,
       f"worst rel {worst:.3g}")
    close("shut-in reaches p_res", post[-1]["p"], d["pres"], rtol=1e-4)
    ok("influx decays with the shut-in", post[-1]["res"]["q"] < 1e-3,
       f"q_end={post[-1]['res']['q']:.3g}")
    ok("cumulative volume is monotone and positive",
       all(b["res"]["vol"] >= a["res"]["vol"] for a, b in zip(fr, fr[1:])) and
       post[-1]["res"]["vol"] > 0.0)

    # ---- 2. continuous underbalanced drilling (no shut-in) -----------------
    ubcase = CASE_DIR / "caseUnderbalanced.dict"
    d2 = parse_case(ubcase)
    rc, out2 = run(ubcase)
    ok("underbalanced arm exits 0", rc == 0, f"rc={rc}")
    f2 = parse_frames(out2, d2["node"])
    wq = wg = wdp = 0.0
    wub = 0.0
    for f in f2:
        r = f["res"]
        wq = max(wq, relerr(r["q"], d2["J"] * max(r["dp"], 0.0), 1e-12))
        wg = max(wg, relerr(r["gasq"], d2["fg"] * r["q"], 1e-12))
        wdp = max(wdp, abs(r["dp"] - (d2["pres"] - f["p"])))
        wub = max(wub, 1.0 if r["ub"] != 1 else 0.0)
    ok("UB: q == J dp on every frame (worst %.3g)" % wq, wq < 1e-9)
    ok("UB: gasq == fg q on every frame (worst %.3g)" % wg, wg < 1e-9)
    ok("UB: dp == p_res - p on every frame (worst %.3g Pa)" % wdp,
       wdp < 1e-2)
    ok("UB: stays underbalanced on every frame", wub == 0.0)
    worst = 0.0
    for prev, cur in zip(f2, f2[1:]):
        Q = cur["Q"].get("WELL", 0.0)
        res = continuity_residual(d2, prev, cur, Q)
        scale = max(abs(inflow(d2, cur["p"])), abs(Q), 1e-6)
        rel = abs(res) / scale
        if (scale > 1e-4 and rel > worst) or (scale <= 1e-4
                                              and abs(res) > 1e-10
                                              and rel > worst):
            worst = rel
    ok("UB: step continuity closes", worst < 1e-9, f"worst={worst:.3g}")
    close("UB: steady inflow feeds the pipe rate",
          f2[-1]["res"]["q"], f2[-1]["Q"]["WELL"], rtol=1e-6)
    close("UB: kill mud closed form", f2[-1]["res"]["killmud"],
          d2["pres"] / (G * d2["depth"]), rtol=1e-11)
    close("UB: SIDPP closed form", f2[-1]["res"]["sidpp"],
          d2["pres"] - d2["rho"] * G * d2["depth"], rtol=1e-11)

    # ---- 3. overbalanced arm: the law must charge exactly nothing ----------
    obcase = CASE_DIR / "caseOverbalanced.dict"
    twin = CASE_DIR / "caseNoReservoir.dict"
    d3 = parse_case(obcase)
    rc, out3 = run(obcase)
    ok("overbalanced arm exits 0", rc == 0, f"rc={rc}")
    f3 = parse_frames(out3, d3["node"])
    ok("overbalanced: q is exactly zero on every frame",
       all(f["res"]["q"] == 0.0 for f in f3))
    ok("overbalanced: no volume booked",
       all(f["res"]["vol"] == 0.0 and f["res"]["gasvol"] == 0.0 for f in f3))
    ok("overbalanced: ub flag is 0", all(f["res"]["ub"] == 0 for f in f3))
    close("overbalanced: kill mud closed form", f3[-1]["res"]["killmud"],
          d3["pres"] / (G * d3["depth"]), rtol=1e-11)
    close("overbalanced: SIDPP closed form (negative = overbalanced)",
          f3[-1]["res"]["sidpp"],
          d3["pres"] - d3["rho"] * G * d3["depth"], rtol=1e-11)

    rc, outt = run(twin)
    ok("twin exits 0", rc == 0, f"rc={rc}")
    ok("twin prints no reservoir token", "reservoir" not in outt)
    stripped = "\n".join(l for l in out3.splitlines()
                         if "reservoir" not in l)
    ok("inactive reservoir is byte-identical to the undeclared twin",
       stripped.rstrip() == outt.rstrip())

    # ---- 4. loud refusals ---------------------------------------------------
    TMP.mkdir(parents=True, exist_ok=True)
    base = ubcase.read_text()

    def variant(tag, body):
        txt = re.sub(r"^reservoir\s*\{.*?^\}", body, base, flags=re.S | re.M)
        p = TMP / ("case_reject_" + tag + ".dict")
        p.write_text(txt)
        rc, out = run(p)
        ok("refuses " + tag, rc != 0 and "error" in out.lower(),
           out.splitlines()[-1] if out else "no output")

    variant("fixedNode", "reservoir\n{\n    RES { node SURF; pressure 4.2e7;"
                         " productivityIndex 2.0e-8; depth 3000; }\n}")
    variant("sourcesOverlap",
            "sources\n{\n    S1 { node BH; rate 0.01; fraction 1; }\n}"
            "\nreservoir\n{\n    RES { node BH; pressure 4.2e7;"
            " productivityIndex 2.0e-8; depth 3000; }\n}")
    variant("twoOnOneNode",
            "reservoir\n{\n    RES { node BH; pressure 4.2e7;"
            " productivityIndex 2.0e-8; depth 3000; }"
            "\n    RES2 { node BH; pressure 4.0e7;"
            " productivityIndex 1.0e-8; depth 3000; }\n}")
    variant("gasFractionOne", "reservoir\n{\n    RES { node BH; pressure 4.2e7;"
                              " productivityIndex 2.0e-8; gasFraction 1;"
                              " depth 3000; }\n}")
    variant("zeroDepth", "reservoir\n{\n    RES { node BH; pressure 4.2e7;"
                         " productivityIndex 2.0e-8; depth 0; }\n}")
    variant("zeroJ", "reservoir\n{\n    RES { node BH; pressure 4.2e7;"
                     " productivityIndex 0; depth 3000; }\n}")
    variant("zeroPressure", "reservoir\n{\n    RES { node BH; pressure 0;"
                            " productivityIndex 2.0e-8; depth 3000; }\n}")

    print()
    print("checks:", checks, "failures:", len(failures))
    if failures:
        for f in failures:
            print("FAILED:", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
