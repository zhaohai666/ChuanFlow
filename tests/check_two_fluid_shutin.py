#!/usr/bin/env python3
"""E2E shut-in defect fix (twoFluid + wallis closure, closure gravitySlip).

The minimal horizontal shut-in used to drive a spurious countercurrent: the
algebraic wallis drift is an absolute axial velocity, so at vM = 0 it kept
hauling liquid toward the closed inlet (EOS pressure unbounded) while the
thermal step - whose receive capacity flux exceeds its drain capacity flux at
the accumulating end - lost the max principle and ran the inlet cell to 383 K
against a 277 K wall. Fixed in two moves: the thermal step now ALWAYS
assembles the conservative receive-capacity form (5.132 default, max
principle exact, diag - |off| = cap/dt + kw), and the opt-in
`closure gravitySlip yes` projects the drift by the local inclination (zero
on horizontal) so a shut-in line cannot even start the spurious axial
countercurrent that used to feed the pressure runaway.

Judgements:
  fixed run (gravitySlip yes):
    1. every frame's every temperature inside
       [min(T_init, T_wall), max(T_init, T_wall, T_boundary)] - 1e-9;
    2. after the schedule zeroes the rates (t >= 60) every cell's vsl/vsg
       is exactly 0 (no spurious countercurrent on the horizontal line);
    3. the liquid inventory sum(alpha dx) is frozen after the shut-in step
       (nothing to transport, nothing moves);
    4. the inlet-cell pressure decreases monotonically after shut-in
       (relaxing toward the static 8 MPa outlet state);
    5. alpha stays in [0, 1] every frame (parsed, not just trusted).
  witness run (same dictionary without the key - unprojected drift):
    6. reproduces the transport bug: the post-shutdown inlet cell shows the
       spurious countercurrent (vsl < 0 with vsg > 0 at the closed rate
       boundary) and the liquid inventory keeps drifting (the fixed run
       freezes it exactly) - the two paths really differ and the key does
       real work;
    7. the witness temperature still stays inside the physical bound - the
       thermal half of the fix is now the DEFAULT assembly, so the runaway
       is confined to the transport/EOS layer (documents the 5.132
       migration: the historical witness at 383 K can no longer occur).
  refusal:
    8. gravitySlip beside closure model full is a loud cfFatal (nonzero
       exit): those closures already drive slip from buoyancy.
"""

import re
import subprocess
import sys
from pathlib import Path

APP = sys.argv[1]
CASE = Path(sys.argv[2])
TMP = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("/tmp")

fails = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))
    if not ok:
        fails.append(name)


def run_case(case_path):
    out = subprocess.run([APP, str(case_path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


def parse_frames(text):
    """frames: list of (t, {cell index -> field dict}); frames start at the
    per-frame 't <t> vm <..> dp <..>' line the app prints each report."""
    frames = []
    cur = None
    cur_t = None
    t_re = re.compile(r"\[chflow\] t (\S+) vm")
    c_re = re.compile(r"\[chflow\] cell (\d+) alpha (\S+) water (\S+) vsl (\S+) "
                      r"vsg (\S+) T (\S+) p (\S+)")
    for line in text.splitlines():
        tm = t_re.match(line)
        if tm:
            if cur is not None:
                frames.append((cur_t, cur))
            cur, cur_t = {}, float(tm.group(1))
            continue
        m = c_re.match(line)
        if m and cur is not None:
            cur[int(m.group(1))] = dict(alpha=float(m.group(2)),
                                        vsl=float(m.group(4)),
                                        vsg=float(m.group(5)),
                                        T=float(m.group(6)),
                                        p=float(m.group(7)))
    if cur is not None:
        frames.append((cur_t, cur))
    return frames


# ---- the fixed run -------------------------------------------------------
rc, text = run_case(CASE)
check("app exits zero", rc == 0, f"rc={rc}")

frames = parse_frames(text)
check("frames parsed", len(frames) >= 3, f"{len(frames)} frames")

T_INIT, T_WALL, T_INLET, T_OUTLET = 278.0, 277.0, 278.0, 277.0
T_LO = min(T_INIT, T_WALL, T_INLET, T_OUTLET) - 1e-9
T_HI = max(T_INIT, T_WALL, T_INLET, T_OUTLET) + 1e-9

worst_hi, worst_lo, worst_frame = -1e300, 1e300, -1
worst_alpha = 0.0
for k, (ft, fr) in enumerate(frames):
    for i, c in fr.items():
        if c["T"] > worst_hi:
            worst_hi, worst_frame = c["T"], ft
        if c["T"] < worst_lo:
            worst_lo = c["T"]
        worst_alpha = max(worst_alpha, c["alpha"], 1.0 - c["alpha"])
check("every temperature inside the physical bound (max principle)",
      worst_hi <= T_HI,
      f"max T {worst_hi:.6f} vs bound {T_HI} (frame t={worst_frame})")
check("no temperature below the physical floor", worst_lo >= T_LO,
      f"min T {worst_lo:.6f} vs floor {T_LO}")
check("alpha stays in [0, 1] every frame", worst_alpha <= 1.0 + 1e-12,
      f"worst excursion {worst_alpha:.3g}")

# shut-in frames: the schedule zeroes the rates at t = 60; frames land at
# t = 0, 10, 20, ..., 120 - everything with t > 60 is post-shutdown
post = [fr for ft, fr in frames if ft > 60.0]
check("post-shutdown frames present", len(post) >= 4, f"{len(post)} frames")
worst_v = max(max(abs(c["vsl"]), abs(c["vsg"])) for fr in post for c in fr.values())
check("post-shutdown velocities exactly zero (no spurious countercurrent)",
      worst_v == 0.0, f"worst |v| {worst_v:.3g}")

# liquid inventory frozen after shut-in (alpha transport has nothing to move)
inv = [(ft, sum(c["alpha"] for c in fr.values())) for ft, fr in frames]
inv_post = [v for ft, v in inv if ft > 60.0]
check("liquid inventory frozen after shut-in",
      max(inv_post) - min(inv_post) == 0.0,
      f"spread {max(inv_post) - min(inv_post):.3g}")

# inlet pressure relaxes monotonically toward the static state
p0 = [(ft, fr[0]["p"]) for ft, fr in frames]
post_p0 = [v for ft, v in p0 if ft > 60.0]
mono = all(post_p0[k] >= post_p0[k + 1] - 1e-6 for k in range(len(post_p0) - 1))
check("inlet pressure relaxes monotonically after shut-in", mono,
      f"p0 sequence {['%.4g' % v for v in post_p0]}")
pre_shut_p0 = [v for ft, v in p0 if ft == 60.0]
check("inlet pressure actually relaxes (drops)",
      post_p0[-1] < pre_shut_p0[0] - 1.0,
      f"{pre_shut_p0[0]:.4g} -> {post_p0[-1]:.4g}")

# ---- the witness run (unprojected drift, no gravitySlip) ------------------
# since 5.132 the thermal assembly is conservative by default, so the
# witness no longer shows the 383 K temperature overrun; the unprojected
# wallis drift alone still drives the spurious countercurrent and the
# EOS pressure runaway.
TMP.mkdir(parents=True, exist_ok=True)
witness = TMP / "twoFluidShutIn_witness.dict"
witness.write_text(CASE.read_text().replace(" gravitySlip yes;", ""))
rc2, text2 = run_case(witness)
check("witness run (unprojected drift) exits zero", rc2 == 0, f"rc={rc2}")
frames2 = parse_frames(text2)
check("witness frames parsed", len(frames2) >= 3, f"{len(frames2)} frames")
over = max(c["T"] for ft, fr in frames2 for c in fr.values())
check("witness temperature stays bounded (thermal fix is now the default)",
      over <= T_HI, f"witness max T {over:.4f} vs bound {T_HI}")
buggy = any(c["vsl"] < -1e-4 and c["vsg"] > 1e-4 for ft, c in
            [(ft, fr[0]) for ft, fr in frames2 if ft > 60.0])
check("witness shows the spurious post-shutdown countercurrent at the inlet",
      buggy)
inv2 = [(ft, sum(c["alpha"] for c in fr.values())) for ft, fr in frames2]
inv2_post = [v for ft, v in inv2 if ft > 60.0]
check("witness liquid inventory keeps drifting (not frozen)",
      max(inv2_post) - min(inv2_post) > 1e-9,
      f"spread {max(inv2_post) - min(inv2_post):.3g}")

# ---- refusal -------------------------------------------------------------
ref = TMP / "twoFluidShutIn_refusal.dict"
ref.write_text(CASE.read_text().replace(" model wallis;", " model full;"))
rc3, text3 = run_case(ref)
check("gravitySlip beside a dynamic/full closure is a loud refusal",
      rc3 != 0 and "gravitySlip" in text3,
      f"rc={rc3}")

print()
if fails:
    print(f"{len(fails)} FAILED: " + ", ".join(fails))
    sys.exit(1)
print("all judgements ok")
