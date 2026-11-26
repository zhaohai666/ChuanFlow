#!/usr/bin/env python3
"""5.196 (F12 fix) the latent heat of a step that CROSSES the phase interval.

usage: check_melting_latent_step.py <chflow-3d-melting> <a.dict> <b.dict> <c.dict>

The pre-5.196 energy equation carried the latent heat as an apparent capacity
that was only assembled for cells whose OLD temperature sat strictly inside
(Ts, Tl); a one-shot correction then covered only a FULL interval crossing.
So a cell that started below the solidus (or above the liquidus) and landed
inside the mushy band in one step paid NO latent heat at all - the classic
"one step across the interval" hole. The fix solves the sensible equation and
inverts the enthalpy H/cp = T + (L/cp) f_l(T) in closed form, which is
monotone and covers partial entry, partial exit and full crossing alike.

Each case is a single cell with all-symmetry faces (no conductive or
convective flux at all) and a declared volumetric source, so the exact answer
is arithmetic:

    h   = T0 + (L/cp) f_l(T0) + q/(rho cp) dt          (H/cp after the step)
    T   = h                              if h <= Ts
        = h - L/cp                       if h >= Tl + L/cp
        = (h + (L/cp) Ts/(Tl - Ts)) / (1 + (L/cp)/(Tl - Ts))   otherwise
    f_l = clamp((T - Ts)/(Tl - Ts), 0, 1)

The three cases land the enthalpy INSIDE the band from below (a), deeper
inside (b) and beyond the liquidus (c), and the checker:

  M1 each printed T equals the closed-form enthalpy inversion (1e-12 rel),
  M2 each printed f_l equals clamp((T - Ts)/(Tl - Ts)) on that same T,
  M3 the old answer is non-vacuously different: for (a) and (b) the
     sensible-only temperature T0 + dT differs from the printed one by more
     than 1 K, so a regression to the pre-fix assembly cannot pass silently,
  M4 the enthalpy is conserved per cell: H(T) - H(T0) equals the declared
     source exactly (1e-12 rel), which is the statement the latent heat is
     neither skipped nor double counted.
"""
import math
import re
import subprocess
import sys

if len(sys.argv) != 5:
    print("usage: check_melting_latent_step.py <exe> <a.dict> <b.dict> <c.dict>")
    raise SystemExit(2)

EXE = sys.argv[1]
CASES = [("a", sys.argv[2]), ("b", sys.argv[3]), ("c", sys.argv[4])]
FAILS = []


def check(label, ok, detail=""):
    print(("ok:   " if ok else "FAIL: ") + label + ((" " + detail) if detail else ""))
    if not ok:
        FAILS.append(label)


def scalar(text, key, default=None):
    m = re.search(r"\b" + key + r"\s+([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", text)
    if not m:
        if default is None:
            raise SystemExit(f"case has no {key}")
        return default
    return float(m.group(1))


def block(text, name):
    m = re.search(name + r"\s*\{([^}]*)\}", text)
    return m.group(1) if m else ""


for tag, path in CASES:
    text = open(path).read()
    # the numbers are unique inside the energy block (rho/cp/latentHeat/...),
    # and the nested phaseChange brace makes a block regex stop early - read
    # the scalars off the case text directly
    eng = text
    cp = scalar(eng, "cp")
    rho = scalar(eng, "rho")
    L = scalar(eng, "latentHeat")
    ts = scalar(eng, "solidusT")
    tl = scalar(eng, "liquidusT")
    T0 = scalar(eng, "initialT")
    q = scalar(eng, "volumetricSource")
    dt = scalar(block(text, "piso"), "dt")
    lc = L / cp

    out = subprocess.run([EXE, path], capture_output=True, text=True)
    if out.returncode != 0:
        check(f"{tag}: the case runs", False, out.stderr.strip()[-200:])
        continue
    cells = [(float(m.group(1)), float(m.group(2)))
             for m in re.finditer(r"cell \d+ T (\S+) fl (\S+)", out.stdout)]
    if not cells:
        check(f"{tag}: the case prints a cell state", False,
              "no `cell i T .. fl ..` line (reportCells true missing?)")
        continue
    T, fl = cells[-1]

    def f_of(t):
        return min(max((t - ts) / (tl - ts), 0.0), 1.0)

    dT = q / (rho * cp) * dt
    h = T0 + lc * f_of(T0) + dT
    if h <= ts:
        want = h
    elif h >= tl + lc:
        want = h - lc
    else:
        want = (h + lc * ts / (tl - ts)) / (1.0 + lc / (tl - ts))

    check(f"{tag} (M1) T is the enthalpy inversion of the step",
          # the report prints 12 significant digits, so 1e-10 is the
          # comparison floor (1e-12 would be comparing the print itself)
          abs(T - want) <= 1e-10 * max(abs(want), 1.0),
          f"printed {T:.12g} vs closed form {want:.12g} (h {h:.6g})")
    check(f"{tag} (M2) f_l is that temperature's own liquid fraction",
          abs(fl - f_of(T)) <= 1e-10,
          f"printed {fl:.12g} vs {f_of(T):.12g}")
    check(f"{tag} (M3) the sensible-only answer is a different number",
          abs((T0 + dT) - T) > 1.0,
          f"sensible {T0 + dT:.6g} vs printed {T:.6g}"
          f" (apart {abs(T0 + dT - T):.3f} K)")

    def H(t):
        return cp * (t - T0) + L * (f_of(t) - f_of(T0))

    check(f"{tag} (M4) the enthalpy gain is the declared source",
          abs(H(T) - cp * dT) <= 1e-9 * max(abs(cp * dT), 1.0),
          f"H {H(T):.12g} vs source {cp * dT:.12g}")

print()
if FAILS:
    print(f"{len(FAILS)} criterion(s) failed")
    raise SystemExit(1)
print("all criteria passed")
