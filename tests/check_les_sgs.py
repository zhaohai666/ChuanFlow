#!/usr/bin/env python3
"""Independent regression checker for the chflow-3d-les SGS closures.

The truth here is analytic, never taken from the solver's own output.  The
case text is parsed independently (regex) to rebuild the mesh geometry and
the algebraic closures:

  * uniform block mesh: V = (Lx Ly Lz)/(Nx Ny Nz), filter width Delta = V^(1/3);
  * linearShear profile U = (rate*(y-y0), 0, 0) gives the constant velocity
    gradient dU_x/dy = rate.  With S_ij = (dU_i/dx_j + dU_j/dx_i)/2 and
    |S| = sqrt(2 S:S), a unidirectional shear has |S| = |rate|;
  * Smagorinsky: nut = (Cs*Delta)^2 |S|  (zero iff Cs = 0 or |S| = 0);
  * WALE: Sd = sym((grad U)^2) - tr/3 vanishes identically for a
    unidirectional shear, so nut = 0 while |S| = |rate| stays nonzero;
  * the volume-averaged kinetic energy of the seeded profile is also known:
    E0 = 0.5*rate^2*< (y-y0)^2 > = 0.5*rate^2*H^2*(4 - 1/Ny^2)/12 for cell
    centres on a uniform y-grid.

Assertions compare the reported t=0 nut/strain against these rebuilt values,
require WALE nut ~ 0 (and cross-run the Smagorinsky twin of the same case to
show a nonzero value is expected there), require Cs = 0 to switch SGS off
without changing the strain, check E is finite/non-negative/bounded and matches
the analytic E0, check the final filterWidth equals V^(1/3), and require loud
cf::Error (cfFatal) refusals for an unknown model and a negative Cs.

usage: check_les_sgs.py <binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

# Application defaults, mirrored from LesModel::Config / profile reading.
DEFAULT_MODEL = "smagorinsky"
DEFAULT_CS = 0.17
DEFAULT_CW = 0.325
DEFAULT_RATE = 1.0

passed = []
fails = []


def check(cond, msg):
    if cond:
        passed.append(msg)
    else:
        fails.append(msg)


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def block(txt, name):
    """Body of the FIRST named { } block, balancing nested braces."""
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("no " + name + " block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced " + name + " block")


def scalar(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        if default is not None:
            return default
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def word(txt, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_]\w*)\s*;", txt)
    if not m:
        if default is not None:
            return default
        raise RuntimeError("no word " + key)
    return m.group(1)


def vector(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", txt)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def run(binary, text, tag, td):
    path = os.path.join(td, tag + ".dict")
    with open(path, "w") as f:
        f.write(text)
    return subprocess.run([binary, path], capture_output=True, text=True)


REPORT = re.compile(
    r"^t\s+(" + NUM + r")\s+nutMin\s+(" + NUM + r")\s+nutMax\s+(" + NUM
    + r")\s+nutMean\s+(" + NUM + r")\s+strainMax\s+(" + NUM
    + r")\s+E\s+(" + NUM + r")\s*$", re.M)


def parse_rows(output):
    rows = []
    for m in REPORT.finditer(output):
        rows.append({"t": float(m.group(1)), "nutMin": float(m.group(2)),
                     "nutMax": float(m.group(3)), "nutMean": float(m.group(4)),
                     "strainMax": float(m.group(5)), "E": float(m.group(6))})
    return rows


def parse_filter_width(output):
    m = re.search(r"^filterWidth\s+(" + NUM + r")\s*$", output, re.M)
    return float(m.group(1)) if m else None


# ------------------------------------------------------------------ setup ---
if len(sys.argv) != 3:
    raise SystemExit("usage: check_les_sgs.py <binary> <case.dict>")
BINARY = sys.argv[1]
CASE = sys.argv[2]
BASE = open(CASE).read()
TAG = os.path.splitext(os.path.basename(CASE))[0]

blk = block(BASE, "block")
lo, hi, nc = vector(blk, "min"), vector(blk, "max"), vector(blk, "cells")
Nx, Ny, Nz = (int(nc[0]), int(nc[1]), int(nc[2]))
Lx, Ly, Lz = (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
nCells = Nx * Ny * Nz
V = Lx * Ly * Lz / nCells
DELTA = V ** (1.0 / 3.0)
nu = scalar(BASE, "nu")

profile = block(BASE, "profile") if re.search(r"\bprofile\s*\{", BASE) else ""
ptype = word(profile, "type", "")
rate = scalar(profile, "rate", DEFAULT_RATE)

les = block(BASE, "les") if re.search(r"\bles\s*\{", BASE) else ""
model = word(les, "model", DEFAULT_MODEL)
Cs = scalar(les, "Cs", DEFAULT_CS)
Cw = scalar(les, "Cw", DEFAULT_CW)

piso = block(BASE, "piso")
dt = scalar(piso, "dt")
endTime = scalar(piso, "endTime")

# Analytic strain for a unidirectional linear shear: |S| = |rate|.
strain_true = abs(rate)
nut_smag = Cs * Cs * DELTA * DELTA * strain_true


# ------------------------------------------------------------- refusals ----
def refusal(text, tag, needle):
    r = run(BINARY, text, tag, td)
    out = r.stdout + r.stderr
    check(r.returncode != 0 and "cf::Error" in out and needle in out,
          "%s: want loud cfFatal (%s), rc=%d out=%s"
          % (tag, needle, r.returncode, out.strip()[-180:]))


with tempfile.TemporaryDirectory(prefix="cf_les_") as td:
    invalid_model = model not in ("smagorinsky", "wale")
    invalid_Cs = not (Cs >= 0) or not math.isfinite(Cs)
    invalid_Cw = not (Cw >= 0) or not math.isfinite(Cw)

    if invalid_model or invalid_Cs or invalid_Cw:
        # The reject case itself: run it and require the validation branch
        # implied by the (first-failing) case values.
        if invalid_model:
            need = "les: unknown model '" + model + "'"
        elif invalid_Cs:
            need = "les: Cs must be finite and non-negative"
        else:
            need = "les: Cw must be finite and non-negative"
        refusal(BASE, TAG, need)
        if invalid_model and invalid_Cs:
            # caseReject carries both faults: repair the model so the
            # negative-Cs guard is reached on its own as well.
            fixed_model = re.sub(r"\bmodel\s+\w+\s*;", "model smagorinsky;",
                                 BASE, count=1)
            refusal(fixed_model, "bad_Cs", "les: Cs must be finite and non-negative")
    else:
        result = run(BINARY, BASE, TAG, td)
        out = result.stdout
        check(result.returncode == 0,
              "solver returned %d: %s" % (result.returncode,
                                          (result.stdout + result.stderr)[-300:]))
        rows = parse_rows(out)
        fw = parse_filter_width(out)
        check(len(rows) >= 2, "expected >=2 report lines, got %d" % len(rows))
        check(fw is not None, "no final filterWidth line")

        if rows and fw is not None:
            t0 = rows[0]
            check(close(t0["t"], 0.0, 0.0, 1e-12), "first report is not t=0")
            check(close(t0["strainMax"], strain_true, 1e-6, 1e-12),
                  "t=0 strainMax %g != analytic |S| %g" % (t0["strainMax"], strain_true))
            check(t0["strainMax"] > 0, "seeded shear has zero strain")

            if model == "smagorinsky":
                check(close(t0["nutMax"], nut_smag, 1e-6, 1e-15),
                      "t=0 nutMax %g != (Cs*Delta)^2|S| %g" % (t0["nutMax"], nut_smag))
                check(close(t0["nutMin"], nut_smag, 1e-6, 1e-15),
                      "t=0 nutMin %g != (Cs*Delta)^2|S| %g" % (t0["nutMin"], nut_smag))
                check(close(t0["nutMean"], nut_smag, 1e-6, 1e-15),
                      "t=0 nutMean %g != (Cs*Delta)^2|S| %g" % (t0["nutMean"], nut_smag))
                if Cs == 0:
                    check(t0["nutMax"] == 0.0 and t0["nutMin"] == 0.0
                          and t0["nutMean"] == 0.0,
                          "Cs=0 must give identically zero nut (got %g/%g/%g)"
                          % (t0["nutMin"], t0["nutMax"], t0["nutMean"]))
                    check(t0["strainMax"] > 0,
                          "Cs=0 must not erase the resolved strain")
                else:
                    check(t0["nutMax"] > 0, "Smagorinsky SGS is inactive")

            elif model == "wale":
                # Sd = 0 on pure shear: WALE nut is numerically zero while
                # Smagorinsky on the very same field is O(1e-3).
                tol = 1e-6 * max(nut_smag, 1e-12)
                check(0.0 <= t0["nutMax"] <= tol,
                      "WALE nutMax %g not ~0 on pure shear (tol %g)" % (t0["nutMax"], tol))
                check(0.0 <= t0["nutMin"] <= tol,
                      "WALE nutMin %g not ~0 on pure shear" % t0["nutMin"])
                check(all(0.0 <= r["nutMax"] <= tol for r in rows),
                      "WALE nut grows non-negligibly on the shear field")
                smag_twin = re.sub(r"\bmodel\s+wale\s*;", "model smagorinsky;",
                                   BASE, count=1)
                rr = run(BINARY, smag_twin, "smag_twin", td)
                twin_rows = parse_rows(rr.stdout)
                check(rr.returncode == 0 and len(twin_rows) >= 1,
                      "Smagorinsky twin of the WALE case did not run")
                if twin_rows:
                    tn = twin_rows[0]["nutMax"]
                    check(tn > 0 and close(tn, nut_smag, 1e-6, 1e-15),
                          "same shear: Smagorinsky nut %g != %g (WALE is ~0)"
                          % (tn, nut_smag))

            # filter width and reported cell count come from the mesh alone
            check(close(fw, DELTA, 1e-6, 1e-15),
                  "filterWidth %g != V^(1/3) %g" % (fw, DELTA))
            mcells = re.search(r"\bcells\s+(\d+)", out)
            check(mcells is not None and int(mcells.group(1)) == nCells,
                  "reported cell count != Nx*Ny*Nz=%d" % nCells)

            # analytic E0 of the seeded linear shear, at t=0
            if ptype == "linearShear":
                H = Ly
                E0 = 0.5 * rate * rate * H * H * (4.0 - 1.0 / (Ny * Ny)) / 12.0
                check(close(t0["E"], E0, 1e-9, 1e-15),
                      "t=0 E %g != analytic E0 %g" % (t0["E"], E0))

            # energy is finite, non-negative and bounded (viscous decay)
            check(all(math.isfinite(r["E"]) for r in rows), "E has non-finite values")
            check(all(r["E"] >= 0.0 for r in rows), "E has negative values")
            check(rows[-1]["E"] <= rows[0]["E"] * (1.0 + 1e-6) + 1e-15,
                  "E is not bounded/decaying: E0=%g final=%g"
                  % (rows[0]["E"], rows[-1]["E"]))
            check(all(rows[i + 1]["t"] >= rows[i]["t"] - 1e-12
                      for i in range(len(rows) - 1)), "time is not non-decreasing")
            check(close(rows[-1]["t"], endTime, 1e-6, 1e-12),
                  "last report t %g != endTime %g" % (rows[-1]["t"], endTime))

        # Both validation guards, reached by mutating this valid case in memory.
        bad_model = re.sub(r"\bmodel\s+\w+\s*;", "model bogus;", BASE, count=1)
        refusal(bad_model, "bad_model", "les: unknown model 'bogus'")
        bad_Cs = re.sub(r"\b(Cs\s+)" + NUM + r"\s*;", r"\g<1>-0.1;", BASE, count=1)
        refusal(bad_Cs, "bad_Cs", "les: Cs must be finite and non-negative")

if fails:
    print("les_sgs FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("les_sgs PASS %s: %d assertions" % (TAG, len(passed)))
