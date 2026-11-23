#!/usr/bin/env python3
"""E7a (§5.177) verification: transient wellbore equipment in the 1D kernel.

usage: check_transient_wellbore_equipment.py <exe> <case.dict> <case_dir>

The checker NEVER trusts a value the kernel prints for an element. It reads
every declaration from the case text, re-implements the three physical laws
locally, and cross-checks the marched trajectory two independent ways:

  * report lines: for each frame it recomputes SSSV K(t)/seal from the close/
    reopen ramps, the choke K from the declared orifice, and the ESP speed
    (exact exponential relaxation) and head (quadratic affinity curve) from
    the printed Q -- then compares to what the app printed;
  * momentum-row reversal: from printed Q(t), Q(t-dt), p(t) and the Churchill
    wall law it reverse-solves the local loss K the SSSV/choke row MUST have
    used, and the head pressure the ESP row MUST have used, and compares those
    to the locally computed declarations.

Physics anchors on this case:
  * SSSV: K ramps 0 -> 2000 over [2,3], the row SEALS (Q == 0) over [3,8)
    while the fed tank's pressure rises, reopen ramps back over [8,9];
  * choke: a single constant orifice K = (A_pipe/(Cd A_o))^2 for the whole run;
  * ESP: speed steps 1 -> 1.4 at t = 3 and relaxes with tau = 1.5 s, lifting
    the flow of a fixed-pressure line above its pre-command value.

A zero-declaration twin is checked to carry no element reports and to read no
extra loss or head (the exact pre-E7a path). Eleven loud rejections cover bad
targets, ramp ordering, negative K, valve/coupled overlap, double-claimed
pipes, malformed choke and bad ESP parameters. Broad byte-identity of the
undeclared path is also witnessed by the transient family regressions.
"""
import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

EXE, CASE, CASE_DIR = sys.argv[1], sys.argv[2], sys.argv[3]

RHO, MU, GRAV = 850.0, 0.1, 9.80665
L, D = 100.0, 0.05
AREA = 0.25 * math.pi * D * D
INERTIA = RHO * L / AREA
KQ_SCALE = RHO / (2.0 * AREA * AREA)
DT = 0.02

# The kernel accumulates t_ by repeated += dt; the schedule is evaluated at
# t_ + dt. Mirror that accumulation in double precision so a boundary at an
# exact multiple of dt lands on the same frame (7.999999999999917 < 8.0).
TACC = [0.0]
for _ in range(1000):
    TACC.append(TACC[-1] + DT)

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


def close(name, got, want, rtol, atol=0.0):
    scale = max(abs(got), abs(want), atol)
    ok(name, abs(got - want) <= rtol * scale,
       f"got {got:.12g} want {want:.12g} diff {abs(got-want):.3g}")


def rel(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def k_rel(a, b):
    """K comparison with an absolute floor: a ~1e-9 reverse-solved K next to
    the declared 0 is rounding noise, not a physics disagreement."""
    return abs(a - b) / max(1.0, abs(a), abs(b))


def T(i):
    """Physical time of step index i."""
    return i * DT


# ---- independent physics ---------------------------------------------------
def friction_pa(q):
    """Churchill (1977) all-regime drop, re-derived (smooth pipe)."""
    if q == 0.0:
        return 0.0
    v = q / AREA
    re = max(RHO * abs(v) * D / MU, 1e-12)
    a1 = (7.0 / re) ** 0.9
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    f = 8.0 * ((8.0 / re) ** 12 + (a + b) ** -1.5) ** (1.0 / 12.0)
    return f * L / D * RHO * v * v / 2.0


def valve_k(k_open, k_closed, cs, ce, os_, oe, t):
    """Continuous close/open ramps, held outside (independent of kernel)."""
    if cs < 0:
        return k_open
    if t < cs:
        return k_open
    if t < ce:
        return k_open + (t - cs) / (ce - cs) * (k_closed - k_open)
    if os_ < 0 or t < os_:
        return k_closed
    if t < oe:
        return k_closed + (t - os_) / (oe - os_) * (k_open - k_closed)
    return k_open


def sssv_sealed(cs, ce, os_, t):
    return cs >= 0 and t >= ce and (os_ < 0 or t < os_)


def esp_speed(s0, cmd, t_cmd, tau, t):
    if t < t_cmd:
        return s0
    if not (tau > 0):
        return cmd
    return cmd + (s0 - cmd) * math.exp(-(t - t_cmd) / tau)


def esp_head(stages, h0, hq, qr, q, n):
    if not (h0 > 0):
        return 0.0
    n2 = n * n
    h = stages * (n2 * h0 - (h0 - hq) * q * q / (qr * qr)) if qr > 0 else stages * n2 * h0
    return max(h, 0.0)


# ---- tiny dictionary reader (reads only what this case declares) -----------
def block_after(text, key):
    """Exact text inside the matching braces of a top-level `key { ... }`."""
    m = re.search(r"(?m)^[ \t]*" + re.escape(key) + r"\s*\{", text)
    if not m:
        return None
    i = m.end() - 1
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1:j]
    return None


def subblocks(text):
    """Immediate named child blocks (name, body) of a dictionary block."""
    out = []
    pat = re.compile(r"(?m)^[ \t]*([A-Za-z_][A-Za-z0-9_]*)\s*\{")
    pos = 0
    while True:
        m = pat.search(text, pos)
        if not m:
            break
        i = m.end() - 1
        depth = 0
        j = i
        for j in range(i, len(text)):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
        out.append((m.group(1), text[i + 1:j]))
        pos = j + 1
    return out


def scalar(body, key, default=None):
    m = re.search(r"(?<![A-Za-z0-9_])" + re.escape(key)
                  + r"[ \t]+([-+0-9.eE]+)[ \t]*;", body)
    return float(m.group(1)) if m else default


def word(body, key):
    m = re.search(r"(?<![A-Za-z0-9_])" + re.escape(key)
                  + r"[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*;", body)
    return m.group(1) if m else None


# ---- read the declarations: the case text is the single truth source -------
with open(CASE) as fh:
    source = fh.read()

s_blk, c_blk, e_blk = (block_after(source, k) for k in ("sssv", "choke", "esp"))
ok("all three wellbore blocks present in the case",
   s_blk is not None and c_blk is not None and e_blk is not None)

s_subs = subblocks(s_blk) if s_blk else []
c_subs = subblocks(c_blk) if c_blk else []
e_subs = subblocks(e_blk) if e_blk else []
ok("exactly one SSSV / choke / ESP declared",
   len(s_subs) == 1 and len(c_subs) == 1 and len(e_subs) == 1,
   f"{len(s_subs)}/{len(c_subs)}/{len(e_subs)}")

_, s_body = s_subs[0]
S_PIPE = word(s_body, "pipe")
S_OK = scalar(s_body, "openK", 0.0)
S_CK = scalar(s_body, "closedK", 1.0e12)
s_cl, s_op = block_after(s_body, "close"), block_after(s_body, "reopen")
S_CS, S_CE = scalar(s_cl, "start"), scalar(s_cl, "end")
S_OS, S_OE = scalar(s_op, "start"), scalar(s_op, "end")
ok("SSSV close/reopen windows read", None not in (S_PIPE, S_CS, S_CE, S_OS, S_OE),
   f"{S_PIPE} [{S_CS},{S_CE}] [{S_OS},{S_OE}]")

_, c_body = c_subs[0]
C_PIPE = word(c_body, "pipe")
C_D = scalar(c_body, "diameter")
C_CD = scalar(c_body, "dischargeCoeff")
C_AO = 0.25 * math.pi * C_D * C_D
C_K = (AREA / (C_CD * C_AO)) ** 2  # orifice -> the row's K convention
ok("choke orifice K re-derived", C_PIPE == "PB" and C_K > 0,
   f"{C_PIPE} K={C_K:.9g}")

_, e_body = e_subs[0]
E_PIPE = word(e_body, "pipe")
E_STAGES = scalar(e_body, "stages", 1.0)
E_H0 = scalar(e_body, "shutoffHead")
E_HQ = scalar(e_body, "headAtRate")
E_QR = scalar(e_body, "rate")
E_S0 = scalar(e_body, "speed0", 1.0)
E_CMD = scalar(e_body, "speedCmd", 1.0)
E_TCMD = scalar(e_body, "speedCmdTime", 0.0)
E_TAU = scalar(e_body, "speedTau", 0.0)
ok("ESP curve read", E_PIPE == "PC" and E_H0 > 0 and E_QR > 0,
   f"{E_PIPE} stages={E_STAGES} H0={E_H0} Hq={E_HQ} Qr={E_QR}")


def run(path):
    return subprocess.run([EXE, path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def walk_traj(text):
    """Parse printed frames into {step_index: {...}} (step = t/DT)."""
    traj = {}
    cur = None
    for ln in text.splitlines():
        m = re.search(r"\[chflow\] t\s+([-+0-9.eE]+)\s*$", ln)
        if m:
            cur = int(round(float(m.group(1)) / DT))
            traj[cur] = {}
            continue
        if cur is None:
            continue
        d = traj[cur]
        m = re.search(r"\[chflow\] node (\w+) p\s+([-+0-9.eE]+)", ln)
        if m:
            d.setdefault("p", {})[m.group(1)] = float(m.group(2)); continue
        m = re.search(r"\[chflow\] pipe (\w+) Q\s+([-+0-9.eE]+)", ln)
        if m:
            d.setdefault("q", {})[m.group(1)] = float(m.group(2)); continue
        m = re.search(r"\[chflow\] pipe (\w+) sssv K\s+([-+0-9.eE]+)\s+sealed\s+(\d+)", ln)
        if m:
            d["sssv"] = (float(m.group(2)), int(m.group(3))); continue
        m = re.search(r"\[chflow\] pipe (\w+) choke K\s+([-+0-9.eE]+)", ln)
        if m:
            d["choke"] = float(m.group(2)); continue
        m = re.search(r"\[chflow\] pipe (\w+) esp speed\s+([-+0-9.eE]+)\s+head\s+([-+0-9.eE]+)", ln)
        if m:
            d["esp"] = (float(m.group(2)), float(m.group(3)))
    return traj


# ---- main case: parse the printed trajectory -------------------------------
proc = run(CASE)
ok("main case runs", proc.returncode == 0, proc.stderr[-400:])
traj = walk_traj(proc.stdout)
times = sorted(traj)
ok("one frame per step plus the initial state",
   len(times) == 1001 and times[0] == 0 and times[-1] == 1000,
   f"{len(times)} frames, first {times[0]} last {times[-1]}")
ok("every frame printed all three element reports",
   all("sssv" in traj[i] and "choke" in traj[i] and "esp" in traj[i]
       for i in times))


def F(t):
    """Frame at the step nearest a physical time (for spot checks)."""
    return traj[int(round(t / DT))]


def q_of(i, pipe):
    return traj[i]["q"][pipe]


def dp_of(i, tank):
    return traj[i]["p"][tank] - traj[i]["p"]["OUT" + tank[-1]]


def reverse_k(i, tank, pipe):
    """Reverse-solve the local K from the implicit momentum row."""
    q, qp = q_of(i, pipe), q_of(i - 1, pipe)
    if abs(q) < 1e-6:
        return None
    return (dp_of(i, tank) - INERTIA / DT * (q - qp) - friction_pa(q)) \
        / (KQ_SCALE * abs(q) * q)


def reverse_head(i, tank, pipe):
    """Reverse-solve the ESP head pressure from the momentum row."""
    q, qp = q_of(i, pipe), q_of(i - 1, pipe)
    return INERTIA / DT * (q - qp) - dp_of(i, tank) + friction_pa(q)


# ---- SSSV: report lines, reversed K, and the seal/flow physics -------------
worst_sv_k, seal_mismatch = 0.0, 0
sv_k_got, sv_k_want = [], []
for i in times:
    kk, sealed = traj[i]["sssv"]
    tt = TACC[i]
    want_k = valve_k(S_OK, S_CK, S_CS, S_CE, S_OS, S_OE, tt)
    want_seal = sssv_sealed(S_CS, S_CE, S_OS, tt)
    worst_sv_k = max(worst_sv_k, k_rel(kk, want_k))
    if sealed != (1 if want_seal else 0):
        seal_mismatch += 1
    if i > 0 and not want_seal:
        got = reverse_k(i, "TANKA", S_PIPE)
        if got is not None:
            sv_k_got.append((tt, got))
            sv_k_want.append((tt, want_k))
ok("SSSV report K is the independent ramp at every frame",
   worst_sv_k < 1e-8, f"worst relative {worst_sv_k:.3g}")
ok("SSSV report seal is the independent latch at every frame",
   seal_mismatch == 0, f"{seal_mismatch} mismatches")
ok("enough unsealed frames to invert the SSSV momentum row", len(sv_k_got) > 400,
   str(len(sv_k_got)))
worst = max((k_rel(g, w), t, g, w) for (t, g), (_, w) in zip(sv_k_got, sv_k_want))
ok("SSSV K reverse-solved from Q(t) equals the ramp",
   worst[0] < 2e-3,
   f"worst at t={worst[1]:.2f}: got {worst[2]:.6g} want {worst[3]:.6g}")

sealed_steps = [i for i in times if sssv_sealed(S_CS, S_CE, S_OS, TACC[i])]
ok("the sealed window is non-empty and Q is exactly 0 in it",
   len(sealed_steps) > 100 and all(q_of(i, S_PIPE) == 0.0 for i in sealed_steps),
   f"{len(sealed_steps)} sealed frames")
p3, p77 = F(3.02)["p"]["TANKA"], F(7.9)["p"]["TANKA"]
ok("the fed tank pressure RISES through the sealed window",
   p77 > p3, f"p(3.02)={p3:.6g} p(7.9)={p77:.6g}")
ok("the close ramp ends in a seal jump (Q before > 0, Q at seal = 0)",
   q_of(149, S_PIPE) > 0 and q_of(150, S_PIPE) == 0.0,
   f"Q(2.98)={q_of(149,S_PIPE):.6g} Q(3.0)={q_of(150,S_PIPE):.6g}")
ok("reopen releases the seal (Q > 0 right after reopenStart)",
   q_of(400, S_PIPE) == 0 and q_of(401, S_PIPE) > 0
   and q_of(475, S_PIPE) > q_of(401, S_PIPE),
   f"Q(8)={q_of(400,S_PIPE):.6g} Q(8.02)={q_of(401,S_PIPE):.6g} "
   f"Q(9.5)={q_of(475,S_PIPE):.6g}")

# ---- choke: constant K in the reports and reversed from every row ----------
worst_ch_report = max(rel(traj[i]["choke"], C_K) for i in times)
ok("choke report K is the declared orifice K at every frame",
   worst_ch_report < 1e-9, f"worst relative {worst_ch_report:.3g}")
ch_k_got = []
for i in times[1:]:
    g = reverse_k(i, "TANKB", C_PIPE)
    if g is not None:
        ch_k_got.append((T(i), g))
ok("enough moving frames to invert the choke momentum row", len(ch_k_got) > 800,
   str(len(ch_k_got)))
worst = max((rel(g, C_K), t, g) for t, g in ch_k_got)
ok("choke K reverse-solved from Q(t) equals the orifice K",
   worst[0] < 1e-3, f"worst at t={worst[1]:.2f}: got {worst[2]:.6g} want {C_K:.6g}")
ok("the choke line blows down after spooling up (Q falls as the tank empties)",
   q_of(100, C_PIPE) > q_of(1000, C_PIPE) > 0,
   f"Q(2)={q_of(100,C_PIPE):.6g} Q(20)={q_of(1000,C_PIPE):.6g}")

# ---- ESP: speed/head reports, reversed head, and the lift physics ----------
worst_n, worst_h_report = 0.0, 0.0
esp_h_got, esp_h_want = [], []
for i in times:
    n_rep, h_rep = traj[i]["esp"]
    tt = TACC[i]
    n_want = esp_speed(E_S0, E_CMD, E_TCMD, E_TAU, tt)
    h_want = esp_head(E_STAGES, E_H0, E_HQ, E_QR, q_of(i, E_PIPE), n_want)
    worst_n = max(worst_n, rel(n_rep, n_want))
    worst_h_report = max(worst_h_report, rel(h_rep, h_want))
    if i > 0:
        # reverse_head is already in Pa; h_rep is in m of head.
        esp_h_got.append((tt, reverse_head(i, "TANKC", E_PIPE)))
        esp_h_want.append((tt, RHO * GRAV * h_rep))
ok("ESP report speed is the exact exponential relaxation at every frame",
   worst_n < 1e-9, f"worst relative {worst_n:.3g}")
ok("ESP report head is the affinity curve at every frame",
   worst_h_report < 1e-9, f"worst relative {worst_h_report:.3g}")
worst = max((rel(g, w), t, g, w) for (t, g), (_, w) in zip(esp_h_got, esp_h_want))
ok("ESP head reverse-solved from Q(t) equals the reported head",
   worst[0] < 2e-3,
   f"worst at t={worst[1]:.2f}: got {worst[2]:.6g} want {worst[3]:.6g}")
ok("ESP speed steps at t=3 then relaxes toward the command",
   abs(F(2.98)["esp"][0] - E_S0) < 1e-9
   and F(3.02)["esp"][0] > E_S0
   and F(20.0)["esp"][0] > F(3.02)["esp"][0],
   f"N(2.98)={F(2.98)['esp'][0]:.6g} N(3.02)={F(3.02)['esp'][0]:.6g} "
   f"N(20)={F(20.0)['esp'][0]:.6g}")
close("ESP speed reaches the commanded value by the end",
      F(20.0)["esp"][0], E_CMD, 1e-3)
ok("the lifted speed raises the flow above its pre-command value",
   q_of(1000, E_PIPE) > 1.2 * q_of(100, E_PIPE),
   f"Q(2)={q_of(100,E_PIPE):.6g} Q(20)={q_of(1000,E_PIPE):.6g}")

# ---- zero-declaration twin: the exact pre-E7a path -------------------------
twin_path = os.path.join(CASE_DIR, "case_noElements.dict")
tw = run(twin_path)
ok("zero-declaration twin runs", tw.returncode == 0, tw.stderr[-400:])
ok("twin prints no element report lines",
   all(tag not in tw.stdout for tag in (" sssv K ", " choke K ", " esp speed ")))
tw_traj = walk_traj(tw.stdout)


def twin_reverse_k(i, tank, pipe):
    q = tw_traj[i]["q"][pipe]
    qp = tw_traj[i - 1]["q"][pipe]
    dp = tw_traj[i]["p"][tank] - tw_traj[i]["p"]["OUT" + tank[-1]]
    if abs(q) < 1e-6:
        return 0.0
    return (dp - INERTIA / DT * (q - qp) - friction_pa(q)) / (KQ_SCALE * abs(q) * q)


for pipe, tank in (("PA", "TANKA"), ("PB", "TANKB")):
    worst = max(abs(twin_reverse_k(i, tank, pipe))
                for i in sorted(tw_traj)[1:])
    ok(f"twin: {pipe} reads no extra local loss", worst < 1e-3,
       f"worst |K| {worst:.3g}")
worst_twin_h = max(
    abs(INERTIA / DT * (tw_traj[i]["q"]["PC"] - tw_traj[i - 1]["q"]["PC"])
        - (tw_traj[i]["p"]["TANKC"] - tw_traj[i]["p"]["OUTC"])
        + friction_pa(tw_traj[i]["q"]["PC"]))
    for i in sorted(tw_traj)[1:])
ok("twin: PC reads no ESP head", worst_twin_h < 1e-2,
   f"worst |headPa| {worst_twin_h:.3g}")

# ---- loud rejections -------------------------------------------------------
rejects = [
    ("SSSV on unknown pipe", "case_rejectSssvBadPipe.dict", "targets unknown pipe"),
    ("SSSV close start > end", "case_rejectSssvCloseOrder.dict", "close start must be <= end"),
    ("SSSV reopen before close ends", "case_rejectSssvReopenEarly.dict", "reopen must not start before close ends"),
    ("SSSV negative K", "case_rejectSssvNegativeK.dict", "K must be >= 0"),
    ("SSSV on a base-valve pipe", "case_rejectSssvOnValve.dict", "one local-loss account"),
    ("choke with both k and diameter", "case_rejectChokeBoth.dict", "one resistance law"),
    ("choke dischargeCoeff out of range", "case_rejectChokeBadCd.dict", "dischargeCoeff must be in"),
    ("ESP with zero stages", "case_rejectEspBadStages.dict", "stages must be positive"),
    ("ESP without shutoff head", "case_rejectEspMissingHead.dict", "positive shutoffHead"),
    ("two elements on one pipe", "case_rejectTwoElements.dict", "two wellbore elements"),
    ("ESP with non-positive speed", "case_rejectEspBadSpeed.dict", "speeds must be positive"),
]
for label, fn, needle in rejects:
    r = run(os.path.join(CASE_DIR, fn))
    ok(f"rejection: {label}", r.returncode != 0 and needle in r.stderr,
       f"rc={r.returncode} stderr={r.stderr[-220:]}")

print(f"\n{checks} checks, {len(failures)} failure(s)")
if failures:
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("PASS: E7a transientWellboreEquipment -- SSSV seal, choke orifice and ESP "
      "affinity/lag re-derived from the marched trace; "
      f"{len(rejects)} rejections, zero-declaration twin clean")
