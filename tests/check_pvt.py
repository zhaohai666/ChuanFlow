#!/usr/bin/env python3
"""Analytic verification of the PVT layer (chflow-pvt).

usage: check_pvt.py <exe> <case.dict>

The app reports an isothermal TP flash at every pressure of the case, the
bubble point of the feed and a tabulated property. Checks are self-contained
thermodynamic identities, not stored numbers:

  flash      phase sums, Rachford-Rice material balance against the feed,
             beta in (0, 1) only between the envelope, liquid denser than
             vapour, light ends enriched in the vapour.
  envelope   beta falls monotonically with pressure, the sweep collapses to a
             single phase above the bubble point and to all vapour below the
             dew side, and the reported bubble point brackets it.
  bubble     the flash at the bubble point holds an infinitesimal vapour.
  table      a table node reproduces the directly flashed property, so the
             interpolated path and the EOS path stay consistent.
"""
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"pvt failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")

# ---- feed composition straight from the case file (exact literals) ----
with open(case) as fh:
    source = fh.read()
comp = re.search(r"composition\s*\{(.*?)\}", source, re.S)
if comp is None:
    raise SystemExit("case has no composition block")
feed = {}
for key, value in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s+([-+0-9.eE]+)\s*;", comp.group(1)):
    feed[key] = float(value)
total = sum(feed.values())
z = {k: v / total for k, v in feed.items()}

temperature = float(re.search(r"temperature\s+([-+0-9.eE]+)\s*;", source).group(1))

# ---- parse the app output ----
flashes, tables, bubble = [], {}, None
current = None

for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok or tok[0] != "pvt":
        continue
    kind = tok[1]
    if kind == "flash":
        kv = dict(zip(tok[2::2], tok[3::2]))
        current = {
            "T": float(kv["T"]), "p": float(kv["p"]), "beta": float(kv["beta"]),
            "single": kv["single"] == "1", "iterations": int(kv["iterations"]),
            "converged": kv["converged"] == "1", "residual": float(kv["residual"]),
            "x": {}, "y": {},
        }
        flashes.append(current)
    elif kind == "liquid" and current is not None:
        current["rhoL"] = float(tok[3])
        current["ML"] = float(tok[5])
    elif kind == "vapour" and current is not None:
        current["rhoV"] = float(tok[3])
        current["MV"] = float(tok[5])
    elif kind in ("x", "y") and current is not None:
        current[kind][tok[2]] = float(tok[3])
    elif kind == "bubble":
        kv = dict(zip(tok[2::2], tok[3::2]))
        bubble = {"T": float(kv["T"]), "p": float(kv["p"])}
        current = None
    elif kind == "table":
        prop = tok[2]
        kv = dict(zip(tok[3::2], tok[4::2]))
        tables[(prop, float(kv["T"]), float(kv["p"]))] = float(kv["value"])
        current = None

if len(flashes) < 3:
    raise SystemExit(f"expected several flashes in the output, got {len(flashes)}")

fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


# ---- per-flash identities ----
for f in flashes:
    label = f"p={f['p']:.6g}"
    check(f"flash {label} temperature", abs(f["T"] - temperature) < 1e-12,
          f"got {f['T']} want {temperature}")
    check(f"flash {label} converged", f["converged"], "successive substitution stalled")
    check(f"flash {label} residual", f["residual"] < 1e-9, f"residual {f['residual']:.3g}")
    check(f"flash {label} liquid composition complete",
          set(f["x"]) == set(z), f"{sorted(f['x'])} vs {sorted(z)}")
    check(f"flash {label} vapour composition complete",
          set(f["y"]) == set(z), f"{sorted(f['y'])} vs {sorted(z)}")
    close(f"flash {label} sum x", sum(f["x"].values()), 1.0, 1e-10)
    close(f"flash {label} sum y", sum(f["y"].values()), 1.0, 1e-10)
    check(f"flash {label} beta in range", 0.0 <= f["beta"] <= 1.0, f"beta {f['beta']}")
    for k in z:
        close(f"flash {label} material balance {k}",
              f["beta"] * f["y"][k] + (1.0 - f["beta"]) * f["x"][k], z[k], 1e-8, 1e-12)
    if not f["single"]:
        check(f"flash {label} liquid denser", f["rhoL"] > f["rhoV"],
              f"{f['rhoL']} vs {f['rhoV']}")
        check(f"flash {label} molar mass order", f["ML"] > f["MV"],
              f"{f['ML']} vs {f['MV']}")
    else:
        # a collapsed flash reports the feed itself in both slots
        close(f"flash {label} single-phase x equals feed",
              sum(abs(f["x"][k] - z[k]) for k in z), 0.0, 0.0, 1e-12)
        close(f"flash {label} single-phase y equals feed",
              sum(abs(f["y"][k] - z[k]) for k in z), 0.0, 0.0, 1e-12)

# ---- the light ends partition to the vapour ----
two_phase = [f for f in flashes if not f["single"]]
check("envelope has a two-phase branch", len(two_phase) >= 2,
      f"{len(two_phase)} two-phase flashes")
for f in two_phase:
    kvals = {c: f["y"][c] / f["x"][c] for c in z}
    most_volatile = max(kvals, key=lambda c: kvals[c])
    check(f"p={f['p']:.6g} methane is the most volatile",
          most_volatile == "methane", f"K = {kvals}")
    check(f"p={f['p']:.6g} vapour enriched in methane",
          f["y"]["methane"] > z["methane"] > f["x"]["methane"],
          f"y={f['y']['methane']:.6g} z={z['methane']:.6g} x={f['x']['methane']:.6g}")

# ---- monotone envelope ----
by_p = sorted(flashes, key=lambda f: f["p"])
betas = [f["beta"] for f in by_p]
for hi, lo in zip(betas, betas[1:]):
    check("vapour fraction falls with pressure", hi >= lo - 1e-12,
          f"{hi} then {lo}")

# ---- single-phase collapse agrees with the bubble point ----
if bubble is not None:
    close("bubble point temperature", bubble["T"], temperature, 1e-12)
    below = [f["p"] for f in flashes if not f["single"] and f["p"] < bubble["p"] * 0.99]
    above = [f["p"] for f in flashes if f["single"] and f["p"] > bubble["p"] * 1.01]
    check("bubble point sits above the two-phase sweep", bool(below),
          f"largest two-phase flash {max(below) if below else None}")
    check("bubble point sits below a single-phase flash", bool(above),
          f"smallest single-phase flash {min(above) if above else None}")
    if below:
        check("bubble point above the last two-phase pressure",
              bubble["p"] > max(below), f"{bubble['p']} vs {max(below)}")
    if above:
        check("bubble point below the first single-phase pressure",
              bubble["p"] < min(above), f"{bubble['p']} vs {min(above)}")
    at = [f for f in flashes if abs(f["p"] - bubble["p"]) <= 1e-6 * bubble["p"]]
    check("bubble point flash reported", bool(at), "no flash at the bubble pressure")
    if at:
        check("bubble point holds an incipient vapour", at[0]["beta"] < 1e-3,
              f"beta {at[0]['beta']:.3g}")

# ---- table nodes agree with the direct EOS evaluation ----
prop = re.search(r"property\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", source)
field = {"liquidDensity": "rhoL", "vapourDensity": "rhoV",
         "vapourFraction": "beta"}.get(prop.group(1) if prop else "")
if field:
    for f in flashes:
        key = (prop.group(1), f["T"], f["p"])
        if key in tables:
            close(f"table node at p={f['p']:.6g} matches the flash",
                  tables[key], f[field], 1e-9)
check("table produced values", bool(tables), "no pvt table lines parsed")

if fails:
    print("\n".join(f"FAIL: {m}" for m in fails))
    raise SystemExit(f"{len(fails)} pvt check(s) failed")

print(f"pvt: {len(flashes)} flashes, {len(tables)} table nodes, "
      f"bubble point {bubble['p'] if bubble else 'n/a'} Pa - all checks passed")
