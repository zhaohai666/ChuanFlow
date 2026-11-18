#!/usr/bin/env python3
"""M6: independent verification of the parallel backend abstraction.

The app sweeps a write-conflict-free domain-local kernel through
cf::backend::parallelFor and prints its checksum, the serial-reference
checksum and the analytic closed form. This checker never imports solver
internals: it reads n / mode from the case, restates both closed forms with
exact rational arithmetic, and requires all three checksums to be bit-for-bit
identical to the restated value (the case size keeps every partial sum inside
the exact-integer range of a double).

It also verifies that the active backend and the serial reference produce the
same elementwise flag, then exercises three parse-time refusals (n = 0,
n < 0, unknown mode) and requires a loud, non-zero exit naming the reason.

usage: check_backend.py <chflow-3d-backend> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile
from fractions import Fraction

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]
case = os.path.abspath(case)
text = re.sub(r"//[^\n]*", "", open(case).read())

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def scalar(key):
    m = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)\s*;", text)
    if not m:
        raise SystemExit(f"case has no scalar '{key}'")
    return float(m.group(1))


def word(key):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_]\w*)\s*;", text)
    if not m:
        raise SystemExit(f"case has no word '{key}'")
    return m.group(1)


n = int(scalar("n"))
mode = word("mode")
if n <= 0 or mode not in ("square", "divide"):
    raise SystemExit(f"shipped case must have n>0 and a known mode, got n={n} mode={mode}")


def exact_checksum(n, mode):
    if mode == "square":
        return (Fraction((n - 1) * n * (2 * n - 1), 6), Fraction((n - 1) ** 2, 1))
    return (Fraction(n * (n + 1), 4), Fraction(n, 2))


exp_sum, exp_max = exact_checksum(n, mode)
want_sum, want_max = float(exp_sum), float(exp_max)
# The case is sized so the floating reduction is exact; document the bound.
check("case size keeps the reduction in the exact-double range",
      max(abs(want_sum), abs(want_max)) <= 2.0 ** 53,
      f"sum={want_sum}")


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    out = proc.stdout
    flags = {}
    sums = {}
    name = None
    seen_n = None
    seen_mode = None
    for line in out.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok or tok[0] != "backend":
            continue
        if tok[1] == "name":
            name = tok[2]
            seen_n = int(tok[4])
            seen_mode = tok[6]
        elif tok[1] == "checksum":
            sums[tok[2]] = (float(tok[4]), float(tok[6]))
        elif tok[1] == "elementwiseMatch":
            flags["elementwise"] = int(tok[2])
            flags["checksum"] = int(tok[4])
    return name, seen_n, seen_mode, sums, flags, out


name, seen_n, seen_mode, sums, flags, out = run(case)
check("app echoes the case n and mode", seen_n == n and seen_mode == mode,
      f"n={seen_n} mode={seen_mode}")
check("app reports a known backend", name in ("cpu-serial", "openmp", "kokkos"),
      f"name={name}")
check("all three checksums present", set(sums) == {"active", "serial", "analytic"},
      f"keys={sorted(sums)}")

for tag in ("active", "serial", "analytic"):
    if tag in sums:
        got_sum, got_max = sums[tag]
        check(f"{tag} sum == exact closed form (bitwise)", got_sum == want_sum,
              f"got {got_sum!r} want {want_sum!r}")
        check(f"{tag} max == exact closed form (bitwise)", got_max == want_max,
              f"got {got_max!r} want {want_max!r}")

check("active backend is bit-identical to the serial reference",
      sums.get("active") == sums.get("serial"))
check("app asserts elementwise and checksum agreement",
      flags.get("elementwise") == 1 and flags.get("checksum") == 1,
      f"flags={flags}")

# ---- parse-time refusals: loud, non-zero, naming the reason -------------
probes = [
    ("n = 0", re.sub(r"\bn\s+[-+0-9.eE]+\s*;", "n 0;", text, count=1), "n must be positive"),
    ("n < 0", re.sub(r"\bn\s+[-+0-9.eE]+\s*;", "n -5;", text, count=1), "n must be positive"),
    ("unknown mode", re.sub(r"\bmode\s+\w+\s*;", "mode bogus;", text, count=1), "unknown mode"),
]
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, "probe.dict")
    for what, mutated, needle in probes:
        if mutated == text:
            raise SystemExit(f"probe '{what}' did not mutate the case")
        open(path, "w").write(mutated)
        rr = subprocess.run([exe, path], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {what}", ok, f"needle '{needle}'")

if fails:
    raise SystemExit(f"backendCpu FAILED: {fails}")
print("backendCpu OK: parallel backend == serial == analytic closed form,"
      " with n<=0 / unknown-mode refusals")
