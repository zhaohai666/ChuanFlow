#!/usr/bin/env python3
"""transientBetaMomentum (A3): the 1D transient row learns beta, declarably.

beta = <v^2>/<v>^2 is the area-averaged momentum-flux factor, >= 1 by
Cauchy-Schwarz. A pipe whose end profiles differ owes its momentum balance

    dpConv = rho*(betaOut - betaIn)/A^2 * q|q|

-- the SAME account a coupled interface measures and installs as its
convective coefficient (5.43), now declarable on a pipe whose own geometry
(a section change at a node, an imposed profile) owns the imbalance. This
closes the cell named five times (5.34 boundary 1, 5.35/5.36/5.37 boundary 2)
and its last half (5.38 boundary 3: the term sits on the SAME row as the
static wall).

Families of judgement:
  1. declaring the identity is declaring nothing: betaIn 1.0 / betaOut 1.0 is
     bit-identical to a pipe with no declaration at all
  2. three geometrically identical pipes at three beta declarations land on
     three different rate traces, every state re-derived in closed form from
     the implicit row's quadratic -- billing lowers the rate, refunding
     raises it, and the direction is the physics, not a fit
  3. the account read back: dpConv == convective * q|q| in every reported
     state, with convective = rho*(betaOut - betaIn)/A^2 from the case's own
     numbers
  4. the term is even in the direction and odd in the flux: the same
     declaration driven backwards carries |q| bit-identical to the forward
     pipe and a dpConv of the opposite sign
  5. yield + beta on one row, moving: every state lands on the implicit root
     of Buckingham's wall AND the flux bill together (bisection in wall
     stress, 5.50's checker shape)
  6. the static branch and the account on the same row: a held pipe reports
     NO convective value worth naming -- Q = 0 kills the q|q| term, the same
     statement as 5.50's valveless static wall
  7. the refusals: a factor below 1 (no profile is fuller than flat), and a
     coupled pipe declaring beta (the interface owns that account)
"""

import math
import os
import re
import subprocess
import sys
import tempfile

TRANSIENT = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO, MU, L, D = 850.0, 0.1, 20.0, 0.04
DRIVE = 1.0e4
AREA = 0.25 * math.pi * D * D
RES = 128.0 * MU * L / (math.pi * D ** 4)
INERT = RHO * L / AREA

# name: (betaIn, betaOut, tau0); the app prints dpConv only when the pair
# differs, and static-hold only when tau0 > 0.
PIPES = {
    "Flat":      (1.0, 1.0, 0.0),
    "Identity":  (1.0, 1.0, 0.0),
    "Bump":      (1.0, 1.2, 0.0),
    "Recover":   (1.2, 1.0, 0.0),
    "BumpRev":   (1.0, 1.2, 0.0),
    "YieldMove": (1.0, 1.2, 1.0),
    "Held":      (1.0, 1.2, 6.0),
}
REVERSED = {"BumpRev"}
CONV = {n: RHO * (bo - bi) / (AREA * AREA) for n, (bi, bo, _) in PIPES.items()}

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def run(binary, path):
    return subprocess.run([binary, path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def q_of_tw(tw, tau0):
    """5.46's Buckingham bracket in wall-stress space (5.50's checker)."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    return (math.pi * (0.5 * D) ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def step_quadratic(q_prev, h, c, drive):
    """Implicit step of  M/h (q - q_prev) - drive + R q + c q|q| = 0, solved
    on the branch the drive points to. c > 0 bills (positive root of the
    upward quadratic), c < 0 refunds (small root; the discriminant is the
    assertion that the refund never outruns the row's own slope)."""
    b = RES + INERT / h
    q_inf_flat = abs(drive) / RES
    going = 1.0 if drive > 0 else -1.0
    # amplitude of the previous state in the direction the row is heading
    back = q_prev * going
    a = abs(drive) + INERT / h * back
    if c == 0.0:
        return going * (a / b)
    if c > 0.0:
        root = (-b + math.sqrt(b * b + 4.0 * c * a)) / (2.0 * c)
    else:
        disc = b * b - 4.0 * abs(c) * a
        assert disc > 0.0, f"refund outran the slope: disc {disc:.3e}"
        root = (b - math.sqrt(disc)) / (2.0 * abs(c))
    return going * root


def march(name, hs, drive):
    bi, bo, tau0 = PIPES[name]
    c = CONV[name]
    q = 0.0
    out = []
    for h in hs:
        if tau0 > 0:
            # bisection in wall stress: residual(tw) is monotone here because
            # INERT/h (~6.8e9) dwarfs 2*c*Q (~6.5e4)
            def residual(tw):
                qq = q_of_tw(tw, tau0)
                return (INERT / h * (qq - q) - drive + 4.0 * L / D * tw
                        + c * qq * abs(qq))
            lo, hi = tau0, max(tau0, 1.0)
            while residual(hi) < 0.0:
                hi = 2.0 * hi
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if residual(mid) > 0.0:
                    hi = mid
                else:
                    lo = mid
                if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
                    break
            q = q_of_tw(0.5 * (lo + hi), tau0)
        else:
            q = step_quadratic(q, h, c, drive)
        out.append(q)
    return out


def parse(stdout):
    times, q, dp_conv, hold = [], {}, {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            times.append(float(tok[1]))
        elif len(tok) >= 4 and tok[0] == "pipe":
            name, key = tok[1], tok[2]
            if tok[2] == "static" and tok[3] == "hold":
                hold.setdefault(name, []).append(int(tok[4]))
                continue
            try:
                val = float(tok[3])
            except ValueError:
                continue
            if key == "Q" and name in PIPES:
                q.setdefault(name, []).append(val)
            elif key == "dpConv":
                dp_conv.setdefault(name, []).append(val)
    return {"times": times, "q": q, "dpConv": dp_conv, "hold": hold}


base = run(TRANSIENT, CASE)
require(base.returncode == 0 and "nan" not in (base.stdout + base.stderr).lower(),
        f"the seven-arm case solves with finite output (rc {base.returncode})")
O = parse(base.stdout)

t_src = open(CASE).read()
t_clean = re.sub(r"//[^\n]*", "", t_src)
t_end = float(re.search(r"endTime ([0-9.]+)", t_clean).group(1))
t_dt = float(re.search(r"\bdt ([0-9.]+)", t_clean).group(1))
n_steps = int(math.ceil(t_end / t_dt - 1e-9))
hs = [min(t_dt, t_end - min((i - 1) * t_dt, t_end)) for i in
      range(1, n_steps + 1)]

# --- 1. the identity sentinel ------------------------------------------------
require(O["q"]["Flat"] == O["q"]["Identity"],
        f"betaIn 1.0 / betaOut 1.0 is bit-identical to declaring nothing, all"
        f" {len(O['q']['Flat'])} states of the two control pipes -- the"
        f" declaration adds a term, not a branch, and an empty term is the"
        f" old row exactly")

# --- 2. three identical geometries, three accounts, closed form --------------
worst = 0.0
for name in ("Bump", "Recover"):
    seq = [0.0] + march(name, hs, DRIVE)
    got = O["q"][name]
    ok = len(got) == len(seq)
    if ok:
        worst = max(worst, max(rel(a, b) for a, b in zip(got, seq)))
    require(ok and all(rel(a, b) < 5e-9 for a, b in zip(got, seq)),
            f"{name}: all {len(got)} states sit on the implicit root of the"
            f" row's own quadratic (worst rel {worst:.2e}) -- the account"
            f" enters the march as a term, not a story")
tail = n_steps - 100   # the developed part of the march
billed = all(b <= f for b, f in zip(O["q"]["Bump"], O["q"]["Flat"])) and all(
    b < f for b, f in zip(O["q"]["Bump"][tail:], O["q"]["Flat"][tail:]))
refunded = all(r >= f for r, f in zip(O["q"]["Recover"], O["q"]["Flat"])) and all(
    r > f for r, f in zip(O["q"]["Recover"][tail:], O["q"]["Flat"][tail:]))
sep_b = rel(O["q"]["Bump"][-1], O["q"]["Flat"][-1])
sep_r = rel(O["q"]["Recover"][-1], O["q"]["Flat"][-1])
require(billed and refunded,
        f"three geometrically identical pipes separate by the DECLARATION"
        f" alone, in the physics' direction (terminal separation {sep_b:.2e}"
        f" billed, {sep_r:.2e} refunded; early steps differ at second order"
        f" in the account, below the print noise, so the tail carries the"
        f" strict comparison): the billed pipe (betaOut > betaIn) runs"
        f" slower, the refunded one (betaIn > betaOut) runs faster than the"
        f" undeclared control")

# --- 3. the account read back -------------------------------------------------
worst_acc = 0.0
for name in ("Bump", "Recover", "YieldMove"):
    for got, qq in zip(O["dpConv"][name], O["q"][name]):
        worst_acc = max(worst_acc, rel(got, CONV[name] * qq * abs(qq)))
require(worst_acc < 1e-9,
        f"dpConv == rho*(betaOut-betaIn)/A^2 * q|q| to {worst_acc:.2e}"
        f" relative (print precision) in every reported state of the three billed pipes -- the"
        f" printed number is the term the row charges, not a separate"
        f" estimate of it")

# --- 4. even in the direction, odd in the flux -------------------------------
seq_rev = [0.0] + march("BumpRev", hs, -DRIVE)
got_rev = O["q"]["BumpRev"]
ok_rev = len(got_rev) == len(seq_rev) and all(
    rel(a, b) < 5e-9 for a, b in zip(got_rev, seq_rev))
require(ok_rev,
        f"BumpRev: all {len(got_rev)} states sit on the backward branch of"
        f" the same quadratic -- the negative root the row's drive asks for")
mag = all(abs(a) == abs(b) for a, b in zip(got_rev, O["q"]["Bump"]))
flip = all(a * b < 0 or (a == 0 and b == 0)
           for a, b in zip(O["dpConv"]["BumpRev"], O["dpConv"]["Bump"]))
require(mag and flip,
        "driven backwards, the SAME declaration carries |q| bit-identical to"
        " the forward pipe and a dpConv of the opposite sign -- the term is"
        " even in the direction, odd in the flux, which is what q|q| means")

# --- 5. yield + beta, one row, moving ----------------------------------------
seq_y = [0.0] + march("YieldMove", hs, DRIVE)
got_y = O["q"]["YieldMove"]
worst_y = max((rel(a, b) for a, b in zip(got_y, seq_y)), default=1.0)
require(len(got_y) == len(seq_y) and worst_y < 1e-8,
        f"YieldMove: all {len(got_y)} states land on the implicit root of"
        f" Buckingham's wall AND the flux bill together (bisection in wall"
        f" stress, worst rel {worst_y:.2e}) -- the account is on the same row"
        f" as the yield, not on a side ledger")

# --- 6. the static branch carries no flux to report ---------------------------
require(all(v == 0.0 for v in O["q"]["Held"]),
        f"Held: every one of its {len(O['q']['Held'])} states is exactly Q = 0"
        f" -- the 12000 Pa wall out-holds the 10000 Pa drive from step 0")
require(all(v == 0 for v in O["hold"]["Held"][:1])
        and all(v == 1 for v in O["hold"]["Held"][1:]),
        "and the latch says so: static hold 1 in every state a step has been"
        " solved for (t = 0 prints 0 by the app's convention -- no step, no"
        " latch yet)")
require(all(abs(v) < 1e-12 for v in O["dpConv"]["Held"]),
        "and it reports NO convective account: Q = 0 kills the q|q| term, the"
        " same statement as 5.50's valveless static wall -- beta and the"
        " static branch finally on one row")

# --- 7. the refusals ----------------------------------------------------------
PROBE_HEAD = """
solver { type transientNetwork; }
fluid { rho 850; mu 0.1; }
time { dt 0.002; endTime 0.01; reportEvery 5; }
nodes
{
    A { fixedPressure 3.0e4; }
    B { fixedPressure 2.0e4; }
}
pipes
{
    Probe { from A; to B; length 20; diameter 0.04; %s; }
}
"""


def variant(declare):
    """A one-pipe probe written from scratch: the refusal must be the ONLY
    thing a solvable file differs by, so the head is the case's own and the
    single pipe carries exactly the declaration under test."""
    return PROBE_HEAD % declare


for declare, why in (
        ("betaOut 0.9", "no profile is fuller than flat"),
        ("betaIn 0.5; betaOut 1.2", "no profile is fuller than flat"),
        ("betaOut -1", "no profile is fuller than flat")):
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "case_probe.dict")
        open(path, "w").write(variant(declare))
        r = run(TRANSIENT, path)
    require(r.returncode != 0 and "below 1" in (r.stdout + r.stderr),
            f"{declare} is refused with the naming message ({why}) -- beta is"
            f" an area average of v^2/vbar^2, and a refusal beats a silently"
            f" empty held set")
coupled = variant("coupledResistance 3.0e6; betaOut 1.2")
with tempfile.TemporaryDirectory() as td:
    path = os.path.join(td, "case_probe.dict")
    open(path, "w").write(coupled)
    r = run(TRANSIENT, path)
require(r.returncode != 0 and "coupled pipe" in (r.stdout + r.stderr),
        "a coupled pipe declaring betaIn/betaOut is refused: the momentum"
        " account of a partner region belongs to the interface, which"
        " measures it and installs the coefficient (5.43) -- a declaration"
        " here would be silently overwritten")

print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
