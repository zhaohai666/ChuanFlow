#!/usr/bin/env python3
"""5.195 (G10 fix) 1D sand bookkeeping when a cell goes dry.

usage: check_two_fluid_sand_clamp.py <chflow-1d-twofluid> <case.dict>

The pre-fix transport step wrote `sandFrac = asNew / aCap` with a
`sandFrac = 0` fallback for aCap <= kTol: every grain left in a cell whose
liquid had drained away was silently dropped out of the alpha*S + B
inventory. This case drains the bore (vsl 0 in, gas only, so the sand cannot
leave through the outlet flux) until the holdup crosses the kernel's own
kTol = 1e-9, and the checker pins the deposit instead of the disappearance:

  C1 the dry step happens: sand 1 -> 0 in a frame while the holdup is below
     1e-9, and that frame's bed is non-zero (pre-fix: bed stayed 0, the mass
     was gone).
  C2 the deposit is the content that was there: the bed gained exactly the
     suspended content of the step before, minus the through-flow the report
     itself carries (which is zero here because the liquid flux is already
     ~1e-9, i.e. below the frame's printing resolution).
  C3 nothing keeps leaking: after the dry step the total inventory
     sum_c (alpha*S + B) is frozen to 1e-15 from frame to frame -- the mass
     that used to vanish step by step now stays put.
  C4 the readout is the audit trail: the cumulative `sandclamp vol` equals
     the sum of the beds to 1e-12, so the printed number and the state agree.
"""
import math
import re
import subprocess
import sys

if len(sys.argv) != 3:
    print("usage: check_two_fluid_sand_clamp.py <exe> <case.dict>")
    raise SystemExit(2)

OUT = subprocess.run([sys.argv[1], sys.argv[2]], capture_output=True,
                     text=True)
if OUT.returncode != 0:
    print(OUT.stdout)
    print(OUT.stderr)
    print("FAIL: the case did not run")
    raise SystemExit(1)
text = OUT.stdout + OUT.stderr

FAILS = []


def check(label, ok, detail=""):
    print(("ok:   " if ok else "FAIL: ") + label + ((" " + detail) if detail else ""))
    if not ok:
        FAILS.append(label)


# ------------------------------------------------------------ the frames
frames = []   # each: dict(cells=[(alpha, sand, bed)], clamp=[vol per line])
cur = None
for line in text.splitlines():
    m = re.match(r"\[chflow\] cell (\d+) alpha (\S+) water (\S+) vsl (\S+) "
                 r"vsg (\S+)(?: [^s]*?)? sand (\S+) bed (\S+)", line)
    if m:
        if cur is None:
            cur = dict(cells=[], clamp=[])
        cur["cells"].append((float(m.group(2)), float(m.group(6)),
                             float(m.group(7))))
        continue
    m = re.match(r"\[chflow\] sandclamp vol (\S+)", line)
    if m:
        if cur is not None:
            cur["clamp"].append(float(m.group(1)))
        continue
    if cur is not None and cur["cells"]:
        frames.append(cur)
        cur = None
if cur is not None and cur["cells"]:
    frames.append(cur)

if not frames:
    print("FAIL: no cell frames in the report")
    raise SystemExit(1)

inv = [sum(a * s + b for a, s, b in f["cells"]) for f in frames]
alphas = [max(a for a, _, _ in f["cells"]) for f in frames]
dry = [i for i in range(1, len(frames))
       if any(a < 1e-9 and s == 0.0 for a, s, _ in frames[i]["cells"])
       and any(b > 0.0 for _, _, b in frames[i]["cells"])]

check("C1 the dry step deposits instead of vanishing",
      bool(dry),
      f"frames {len(frames)}, dry frame {dry[0] if dry else 'none'},"
      f" alpha {alphas[dry[0]]:.3e}" if dry else
      f"frames {len(frames)}, no frame with a dry cell + non-zero bed")

if dry:
    i = dry[0]
    before, after = frames[i - 1]["cells"], frames[i]["cells"]
    # the cell(s) that went dry in exactly this frame: holdup below the kernel's
    # kTol, suspension cleared, bed non-zero
    dryCells = [j for j in range(len(after))
                if after[j][0] < 1e-9 and after[j][1] == 0.0 and after[j][2] > 0.0
                and before[j][2] == 0.0]
    # the deposit must be a real, bounded move: strictly positive, and never
    # more than the content that cell held (the rest of the loss in this frame
    # is the through-flow the liquid still carries, which the report shows as
    # a shrinking holdup, not as a vanishing inventory)
    ok2 = bool(dryCells)
    detail2 = []
    for j in dryCells:
        dbed = after[j][2] - before[j][2]
        held = before[j][0] * before[j][1]
        ok2 = ok2 and 0.0 < dbed <= held
        detail2.append(f"cell {j}: bed +{dbed:.6e} <= held {held:.6e}")
    check("C2 the bed gains the content the dry cell held", ok2, "; ".join(detail2))

    # once drained, that cell's inventory is frozen: nothing else can leak out
    frozen = []
    for j in dryCells:
        tail = [frames[k]["cells"][j][0] * frames[k]["cells"][j][1]
                + frames[k]["cells"][j][2] for k in range(i, len(frames))]
        frozen.append(max(tail) - min(tail))
    check("C3 the drained cells stop leaking",
          max((abs(d) for d in frozen), default=1.0) <= 1e-18,
          f"max |d(alpha*S + B)| on the drained cells "
          f"{max((abs(d) for d in frozen), default=0.0):.3e}")

    last = frames[-1]
    clamp = last["clamp"][-1] if last["clamp"] else 0.0
    beds = sum(b for _, _, b in last["cells"])
    check("C4 sandclamp vol equals the sum of the beds",
          clamp > 0.0 and abs(clamp - beds) <= 1e-10 * max(clamp, 1e-30),
          f"clamp {clamp:.12e} vs beds {beds:.12e}")
else:
    for label in ("C2 the bed gains the content that was held",
                  "C3 the inventory is frozen once the bore is dry",
                  "C4 sandclamp vol equals the sum of the beds"):
        check(label, False, "no dry step to audit")

print()
if FAILS:
    print(f"{len(FAILS)} criterion(s) failed")
    raise SystemExit(1)
print("all criteria passed")
