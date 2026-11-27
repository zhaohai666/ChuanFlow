#!/usr/bin/env python3
"""5.226 (W-class) verification: distributed well control in the 1D transient
kernel -- a marched segmented annulus (drift-flux gas migration, pV=const
expansion, distributed Darcy-Weisbach friction) under a closing BOP ram.

usage: check_transient_well_control_distributed.py <exe> <case.dict> <tmpdir>

The checker never trusts a number the kernel prints for the column. It reads
the declarations out of the case text and re-implements the march locally:

  * migration arm (case.dict): the explicit scheme
        alphaG = clamp(vs / (p A L), 0, 1)          [pV = const expansion]
        vG     = vM + vdrift (1 - alphaG)           [graded drift-flux slip]
        flux   = vG vs / L                          [upwind, standard volume]
        vs    += dt (flux_in - flux_out)
    re-derived frame by frame and compared to every printed cell, plus:
    the gas inventory sum(alphaG A L p) conserved to 1e-10, and the inventory
    centroid rising at the declared drift speed;
  * BOP: kv == openKv before closeTime exactly, a linear ramp over the stroke,
    kv == closedKv exactly after closeTime + strokeTime, and the outflow law
    q == kv sqrt(|pTop - pSurface - sealPressure|) on every frame;
  * friction arm (caseFriction.dict): steady single-phase -- alphaG == 0, the
    per-cell dpf equals the closed-form Darcy-Weisbach number
    f rhoM |vM| vM L / (2 De) (Churchill on the declared roughness) from the
    printed vM, the total dp equals the sum over segments to 1e-9, and the
    reservoir PI law re-used against the column's bottom pressure holds;
  * passive twin: caseTwin.dict (the main deck with annulus{}/bop{} deleted)
    prints byte-identical output to the main deck minus its annulus/bop lines;
  * eight loud refusals (innerD >= outerD, negative length, negative
    roughness, closeTime < 0, strokeTime <= 0, name collision, bop without
    annulus, empty segments).
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


def run(path):
    out = subprocess.run([EXE, str(path)], capture_output=True, text=True)
    return out.returncode, out.stdout + out.stderr


# ---- the kernel's own closed forms, re-implemented --------------------------
def churchill(re_, rr):
    """Churchill (1977), the same expression the steady kernel reads."""
    re_ = max(re_, 1e-12)
    a1 = (7.0 / re_) ** 0.9 + 0.27 * rr
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re_) ** 16.0
    return 8.0 * ((8.0 / re_) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def bop_pos_at(close_time, stroke, t):
    if t < close_time:
        return 0.0
    if not stroke > 0 or t >= close_time + stroke:
        return 1.0
    return (t - close_time) / stroke


def bop_kv_at(open_kv, closed_kv, close_time, stroke, t):
    if t < close_time:
        return open_kv
    if not stroke > 0 or t >= close_time + stroke:
        return closed_kv
    return open_kv + (t - close_time) / stroke * (closed_kv - open_kv)


# ---- declaration parsing (from the case text, never from the report) --------
def block_inner(text, opener_idx):
    """Inner text of the {...} whose '{' sits at opener_idx."""
    depth = 0
    i = opener_idx
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[opener_idx + 1:i], i + 1
        i += 1
    raise ValueError("unbalanced braces")


def scal(text, key, default=None):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)", text)
    if m:
        return float(m.group(1))
    if default is not None:
        return default
    raise KeyError(key)


def parse_case(path):
    text = Path(path).read_text()
    d = {"text": text}
    fm = re.search(r"fluid\s*\{", text)
    fblk, _ = block_inner(text, fm.end() - 1)
    d["rho"] = scal(fblk, "rho")
    d["mu"] = scal(fblk, "mu")
    tm = re.search(r"time\s*\{", text)
    tblk, _ = block_inner(text, tm.end() - 1)
    d["dt"] = scal(tblk, "dt")
    d["endTime"] = scal(tblk, "endTime")
    d["reportEvery"] = int(scal(tblk, "reportEvery", 1))

    am = re.search(r"^annulus\s*\{", text, re.M)
    d["ann"] = None
    if am:
        ablk, _ = block_inner(text, am.end() - 1)
        nm = re.search(r"\w+\s*\{", ablk)
        ann_name = nm.group(0).split()[0]
        body, _ = block_inner(ablk, nm.end() - 1)
        a = {"name": ann_name}
        a["bottomPressure"] = scal(body, "bottomPressure", 0.0)
        a["drift"] = scal(body, "driftVelocity", 0.0)
        a["rhoG"] = scal(body, "gasDensity", 1.2)
        a["pSurf"] = scal(body, "surfacePressure", 1.0e5)
        a["T"] = scal(body, "temperature", 288.15)
        a["influx"] = scal(body, "influx", 0.0)
        a["fgIn"] = scal(body, "influxGasFraction", 0.0)
        a["fDecl"] = scal(body, "frictionFactor", 0.0)
        pm = re.search(r"reservoir\s+(\w+)\s*;", body)
        a["resName"] = pm.group(1) if pm else None
        if a["resName"]:
            rm = re.search(r"^reservoir\s*\{", text, re.M)
            rblk, _ = block_inner(text, rm.end() - 1)
            rbody, _ = block_inner(rblk, re.search(r"\w+\s*\{", rblk).end() - 1)
            a["pRes"] = scal(rbody, "pressure")
            a["J"] = scal(rbody, "productivityIndex")
            a["resFg"] = scal(rbody, "gasFraction", 0.0)
        segs = []
        sm = re.search(r"segments\s*\{", body)
        sblk, _ = block_inner(body, sm.end() - 1)
        i = 0
        while True:
            gm = re.search(r"\w+\s*\{", sblk[i:])
            if not gm:
                break
            gbody, end = block_inner(sblk, i + gm.end() - 1)
            segs.append({
                "innerD": scal(gbody, "innerD"),
                "outerD": scal(gbody, "outerD"),
                "length": scal(gbody, "length"),
                "roughness": scal(gbody, "roughness", 0.0),
                "a0": scal(gbody, "initialGasFraction", 0.0),
            })
            i = end
        a["segs"] = segs
        d["ann"] = a

    bm = re.search(r"^bop\s*\{", text, re.M)
    d["bop"] = None
    if bm:
        bblk, _ = block_inner(text, bm.end() - 1)
        nm = re.search(r"\w+\s*\{", bblk)
        b = {"name": nm.group(0).split()[0]}
        bbody, _ = block_inner(bblk, nm.end() - 1)
        b["annulus"] = re.search(r"annulus\s+(\w+)\s*;", bbody).group(1)
        b["closeTime"] = scal(bbody, "closeTime")
        b["stroke"] = scal(bbody, "strokeTime", 0.0)
        b["openKv"] = scal(bbody, "openKv", 0.0)
        b["closedKv"] = scal(bbody, "closedKv", 0.0)
        b["sealP"] = scal(bbody, "sealPressure", 0.0)
        d["bop"] = b
    return d


# ---- report parsing ----------------------------------------------------------
FRAME = re.compile(r"\[chflow\] t (\S+)")
ASEG = re.compile(r"\[chflow\] annulus (\w+) seg (\d+) z (\S+) alphaG (\S+) "
                  r"vM (\S+) p (\S+) T (\S+) dpf (\S+)")
BOPRE = re.compile(r"\[chflow\] bop (\w+) t (\S+) pos (\S+) kv (\S+) q (\S+) "
                   r"pTop (\S+)")


def parse_frames(text):
    frames = []
    cur = None
    for line in text.splitlines():
        fm = FRAME.match(line)
        if fm:
            if cur is not None:
                frames.append(cur)
            cur = {"t": float(fm.group(1)), "ann": [], "bop": None}
            continue
        if cur is None:
            continue
        am = ASEG.match(line)
        if am:
            cur["ann"].append(dict(name=am.group(1), seg=int(am.group(2)),
                                   z=float(am.group(3)),
                                   alphaG=float(am.group(4)),
                                   vM=float(am.group(5)), p=float(am.group(6)),
                                   T=float(am.group(7)),
                                   dpf=float(am.group(8))))
            continue
        bm = BOPRE.match(line)
        if bm:
            cur["bop"] = dict(name=bm.group(1), t=float(bm.group(2)),
                              pos=float(bm.group(3)), kv=float(bm.group(4)),
                              q=float(bm.group(5)), pTop=float(bm.group(6)))
    if cur is not None:
        frames.append(cur)
    return frames


# ---- the independent re-derivation of the column march ----------------------
class Column:
    """Mirrors TransientNetwork::annulusStep cell-for-cell, in the same order
    of operations, so the comparison is a scheme check, not a curve fit."""

    def __init__(self, d):
        a = d["ann"]
        self.a = a
        self.bop = d["bop"]
        self.rho = d["rho"]
        self.mu = d["mu"]
        self.n = len(a["segs"])
        self.A = [0.25 * math.pi * (s["outerD"] ** 2 - s["innerD"] ** 2)
                  for s in a["segs"]]
        self.De = [s["outerD"] - s["innerD"] for s in a["segs"]]
        self.L = [s["length"] for s in a["segs"]]
        self.rr = [s["roughness"] for s in a["segs"]]
        z = 0.0
        self.zc = []
        for s in a["segs"]:
            self.zc.append(z + 0.5 * s["length"])
            z += s["length"]
        # seed: declared holdups price the at-rest profile, inventory booked
        # as standard volume
        self.p = [0.0] * self.n
        pf = a["pSurf"] + (self.bop["sealP"] if self.bop else 0.0)
        for k in range(self.n - 1, -1, -1):
            a0 = a["segs"][k]["a0"]
            rhoM = a0 * a["rhoG"] + (1.0 - a0) * self.rho
            head = rhoM * G * self.L[k]
            self.p[k] = pf + 0.5 * head
            pf += head
        self.vs = [a["segs"][k]["a0"] * self.p[k] * self.A[k] * self.L[k]
                   for k in range(self.n)]
        self.pTop = a["pSurf"] + (self.bop["sealP"] if self.bop else 0.0)
        self.last = None

    def step(self, dt, t_new):
        a = self.a
        kv = None
        pos = 0.0
        if self.bop:
            pos = bop_pos_at(self.bop["closeTime"], self.bop["stroke"], t_new)
            kv = bop_kv_at(self.bop["openKv"], self.bop["closedKv"],
                           self.bop["closeTime"], self.bop["stroke"], t_new)
        q = (a["J"] * max(a["pRes"] - a["bottomPressure"], 0.0)
             if a["resName"] else a["influx"])
        if self.bop and kv <= 0.0:
            q = 0.0
        fg_in = a["resFg"] if a["resName"] else a["fgIn"]

        alpha, rhoM, vM, dpf = [], [], [], []
        for k in range(self.n):
            vc = self.A[k] * self.L[k]
            al = self.vs[k] / (self.p[k] * vc) if self.p[k] > 0 and vc > 0 else 0.0
            al = min(max(al, 0.0), 1.0)
            alpha.append(al)
            rhoM.append(al * a["rhoG"] + (1.0 - al) * self.rho)
            vM.append(q / self.A[k])
            if q != 0.0:
                f = a["fDecl"] if a["fDecl"] > 0 else churchill(
                    rhoM[k] * abs(vM[k]) * self.De[k] / self.mu,
                    self.rr[k] / self.De[k])
                dpf.append(f * rhoM[k] * abs(vM[k]) * vM[k] * self.L[k]
                           / (2.0 * self.De[k]))
            else:
                dpf.append(0.0)

        pf = a["pSurf"]
        if self.bop:
            pf += self.bop["sealP"]
            if kv > 0.0:
                pf += (q / kv) * (q / kv)
        p_top = pf
        pnew = [0.0] * self.n
        for k in range(self.n - 1, -1, -1):
            head = rhoM[k] * G * self.L[k] + dpf[k]
            pnew[k] = pf + 0.5 * head
            pf += head

        flux = [(vM[k] + a["drift"] * (1.0 - alpha[k])) * self.vs[k] / self.L[k]
                for k in range(self.n)]
        for k in range(self.n):
            inn = q * fg_in * pnew[0] if k == 0 else flux[k - 1]
            out = q * alpha[k] * pnew[k] if k + 1 == self.n else flux[k]
            self.vs[k] += dt * (inn - out)
            if self.vs[k] < 0.0:
                self.vs[k] = 0.0
        self.p = pnew
        self.pTop = p_top
        read = []
        for k in range(self.n):
            vc = self.A[k] * self.L[k]
            ag = (min(max(self.vs[k] / (self.p[k] * vc), 0.0), 1.0)
                  if self.p[k] > 0 and vc > 0 else 0.0)
            read.append(dict(seg=k, z=self.zc[k], alphaG=ag, vM=vM[k],
                             p=self.p[k], T=a["T"], dpf=dpf[k]))
        bop = None
        if self.bop:
            bop = dict(t=t_new, pos=pos, kv=kv, q=q, pTop=p_top)
        self.last = (read, bop)
        return read, bop

    def inventory(self):
        return sum(self.vs)


def march_frames(d):
    """The app's exact frame sequence: the constructor's t = 0 evaluation,
    then ceil(endTime/dt) steps with the last one clamped to the remainder."""
    col = Column(d)
    out = [col.step(0.0, 0.0)]
    dt, end = d["dt"], d["endTime"]
    n_steps = math.ceil(end / dt - 1e-9)
    t = 0.0
    for i in range(1, n_steps + 1):
        t_prev = min((i - 1) * dt, end)
        h = min(dt, end - t_prev)
        out.append(col.step(h, t + h))
        t += h
    return out, col


def worst_rel(pairs, key, floor=0.0):
    w = 0.0
    for got, want in pairs:
        scale = max(abs(got), abs(want), floor)
        w = max(w, abs(got - want) / scale if scale > 0 else 0.0)
    return w


def main():
    # ---- 1. migration arm: the scheme, the inventory, the ram ----------------
    d = parse_case(CASE)
    ok("case declares an annulus", d["ann"] is not None)
    ok("case declares a bop", d["bop"] is not None)
    rc, out = run(CASE)
    ok("migration arm exits 0", rc == 0, f"rc={rc}")
    fr = parse_frames(out)
    n = len(d["ann"]["segs"])
    ok("every frame prints one line per segment",
       all(len(f["ann"]) == n for f in fr), f"n={n}")
    ok("every frame prints one bop line", all(f["bop"] is not None for f in fr))

    want, col = march_frames(d)
    ok("frame count matches the re-derived march", len(fr) == len(want),
       f"got {len(fr)} want {len(want)}")

    wa = wp = wd = wv = 0.0
    wk = wq = wpt = wpos = 0.0
    for f, (read, bop) in zip(fr, want):
        for cell, exp in zip(f["ann"], read):
            wa = max(wa, worst_rel([(cell["alphaG"], exp["alphaG"])], "", 1e-9))
            wp = max(wp, worst_rel([(cell["p"], exp["p"])], "", 1.0))
            wd = max(wd, worst_rel([(cell["dpf"], exp["dpf"])], "", 1e-9))
            wv = max(wv, worst_rel([(cell["vM"], exp["vM"])], "", 1e-9))
        if f["bop"] and bop:
            wk = max(wk, worst_rel([(f["bop"]["kv"], bop["kv"])], "", 1e-12))
            wq = max(wq, worst_rel([(f["bop"]["q"], bop["q"])], "", 1e-12))
            wpt = max(wpt, worst_rel([(f["bop"]["pTop"], bop["pTop"])], "", 1.0))
            wpos = max(wpos, worst_rel([(f["bop"]["pos"], bop["pos"])], "", 1e-12))
    ok("alphaG matches the re-derived slip march (worst %.3g)" % wa, wa < 1e-9)
    ok("p matches the re-derived profile (worst %.3g)" % wp, wp < 1e-9)
    ok("dpf matches the re-derived friction (worst %.3g)" % wd, wd < 1e-9)
    ok("vM matches the re-derived string rate (worst %.3g)" % wv, wv < 1e-9)
    ok("bop kv matches the schedule (worst %.3g)" % wk, wk < 1e-9)
    ok("bop q matches the law (worst %.3g)" % wq, wq < 1e-9)
    ok("bop pTop matches the profile (worst %.3g)" % wpt, wpt < 1e-9)
    ok("bop pos matches the stroke (worst %.3g)" % wpos, wpos < 1e-9)

    # gas inventory: sum(alphaG * A * L * p) is the standard volume, conserved
    A = [0.25 * math.pi * (s["outerD"] ** 2 - s["innerD"] ** 2)
         for s in d["ann"]["segs"]]
    L = [s["length"] for s in d["ann"]["segs"]]

    def inventory(f):
        return sum(c["alphaG"] * A[c["seg"]] * L[c["seg"]] * c["p"]
                   for c in f["ann"])

    inv0 = inventory(fr[0])
    wiv = max(abs(inventory(f) - inv0) / inv0 for f in fr)
    ok("gas inventory conserved (worst %.3g)" % wiv, wiv < 1e-10)

    # the profile shifts upward at the declared drift speed: the inventory
    # centroid rises, bounded by the slip law's own limits
    def centroid(f):
        num = sum(c["alphaG"] * A[c["seg"]] * L[c["seg"]] * c["p"] * c["z"]
                  for c in f["ann"])
        return num / inventory(f)

    cz = [centroid(f) for f in fr]
    rise = cz[-1] - cz[0]
    drift = d["ann"]["drift"]
    a0 = d["ann"]["segs"][0]["a0"]
    ok("gas centroid rises monotonically",
       all(b >= a - 1e-12 for a, b in zip(cz, cz[1:])))
    ok("centroid shift matches the drift speed (%.3g m over %.3g s)" %
       (rise, d["endTime"]),
       0.3 * drift * (1.0 - a0) * d["endTime"] < rise < drift * d["endTime"])
    ok("gas reached the second segment", fr[-1]["ann"][1]["alphaG"] > 1e-6)

    # BOP schedule, asserted on the printed values against the case text
    ct, st = d["bop"]["closeTime"], d["bop"]["stroke"]
    ok("kv == openKv exactly before closeTime",
       all(f["bop"]["kv"] == d["bop"]["openKv"] for f in fr if f["t"] < ct))
    ok("kv == closedKv exactly after closeTime + strokeTime",
       all(f["bop"]["kv"] == d["bop"]["closedKv"]
           for f in fr if f["t"] >= ct + st))
    ramp = [f for f in fr if ct < f["t"] < ct + st]
    ok("kv ramps linearly during the stroke",
       all(abs(f["bop"]["kv"] - (d["bop"]["openKv"]
                                + (f["t"] - ct) / st
                                * (d["bop"]["closedKv"]
                                   - d["bop"]["openKv"]))) < 1e-12 for f in ramp)
       or not ramp)
    p_surf, seal = d["ann"]["pSurf"], d["bop"]["sealP"]
    wlaw = max(abs(f["bop"]["q"] - f["bop"]["kv"]
                   * math.sqrt(max(abs(f["bop"]["pTop"] - p_surf - seal), 0.0)))
               for f in fr)
    ok("q == kv sqrt(|dp|) on every frame (worst %.3g m3/s)" % wlaw, wlaw < 1e-9)

    # ---- 2. friction arm: steady single-phase, closed-form Darcy sum --------
    fcase = CASE_DIR / "caseFriction.dict"
    d2 = parse_case(fcase)
    rc, out2 = run(fcase)
    ok("friction arm exits 0", rc == 0, f"rc={rc}")
    fr2 = parse_frames(out2)
    ok("friction arm is single-phase (alphaG == 0 every frame)",
       all(c["alphaG"] == 0.0 for f in fr2 for c in f["ann"]))
    ok("friction arm is steady (dpf constant across frames)",
       all(abs(f["ann"][k]["dpf"] - fr2[1]["ann"][k]["dpf"])
           <= 1e-12 * (abs(fr2[1]["ann"][k]["dpf"]) + 1.0)
           for f in fr2 for k in range(len(fr2[1]["ann"]))))
    A2 = [0.25 * math.pi * (s["outerD"] ** 2 - s["innerD"] ** 2)
          for s in d2["ann"]["segs"]]
    De2 = [s["outerD"] - s["innerD"] for s in d2["ann"]["segs"]]
    L2 = [s["length"] for s in d2["ann"]["segs"]]
    rr2 = [s["roughness"] for s in d2["ann"]["segs"]]
    # the PI law re-used against the column bottom pressure
    q_pi = d2["ann"]["J"] * max(d2["ann"]["pRes"] - d2["ann"]["bottomPressure"],
                                0.0)
    ok("string rate is the reservoir PI law on the bottom pressure",
       all(f["bop"]["q"] == q_pi for f in fr2[1:]),
       f"q={fr2[-1]['bop']['q']} want {q_pi}")
    wsum = 0.0
    for f in fr2[1:]:
        tot_print = sum(c["dpf"] for c in f["ann"])
        tot_closed = 0.0
        for k, c in enumerate(f["ann"]):
            re_ = d2["rho"] * abs(c["vM"]) * De2[k] / d2["mu"]
            fr_ = (d2["ann"]["fDecl"] if d2["ann"]["fDecl"] > 0
                   else churchill(re_, rr2[k] / De2[k]))
            tot_closed += fr_ * d2["rho"] * abs(c["vM"]) * c["vM"] * L2[k] \
                / (2.0 * De2[k])
        wsum = max(wsum, abs(tot_print - tot_closed)
                   / max(abs(tot_closed), 1e-12))
    ok("total dp == closed-form Darcy-Weisbach sum (worst %.3g)" % wsum,
       wsum < 1e-9)
    ok("friction arm q == kv sqrt(dp) law holds",
       all(abs(f["bop"]["q"] - f["bop"]["kv"]
               * math.sqrt(max(f["bop"]["pTop"] - d2["ann"]["pSurf"]
                               - d2["bop"]["sealP"], 0.0))) < 1e-9
           for f in fr2))
    ok("friction arm kv == openKv (never closes in the run)",
       all(f["bop"]["kv"] == d2["bop"]["openKv"] for f in fr2))

    # ---- 3. passive twin: deleting the blocks is byte-inert -----------------
    twin = CASE_DIR / "caseTwin.dict"
    rc, out_t = run(twin)
    ok("twin exits 0", rc == 0, f"rc={rc}")
    ok("twin prints no annulus token", "annulus" not in out_t)
    ok("twin prints no bop token", " bop " not in out_t)
    stripped = "\n".join(l for l in out.splitlines()
                         if " annulus " not in l and " bop " not in l)
    ok("undeclared twin is byte-identical to the declared deck minus its lines",
       stripped.rstrip() == out_t.rstrip())

    # ---- 4. loud refusals ---------------------------------------------------
    TMP.mkdir(parents=True, exist_ok=True)
    base = CASE.read_text()

    def variant(tag, pattern, repl, needle):
        txt = re.sub(pattern, repl, base, count=1, flags=re.S | re.M)
        assert txt != base, "variant " + tag + " did not change the text"
        p = TMP / ("case_reject_" + tag + ".dict")
        p.write_text(txt)
        rc, out = run(p)
        ok("refuses " + tag, rc != 0 and needle in out,
           out.splitlines()[-1] if out else "no output")

    variant("innerGeOuter", r"innerD 0\.2; outerD 0\.4",
            "innerD 0.5; outerD 0.4", "innerD")
    variant("negativeLength", r"length 10; roughness",
            "length -10; roughness", "length")
    variant("negativeRoughness", r"outerD 0\.4; length 10; roughness 1\.0e-4",
            "outerD 0.4; length 10; roughness -1.0e-4", "roughness")
    variant("negativeCloseTime", r"closeTime 5;", "closeTime -5;", "closeTime")
    variant("zeroStroke", r"strokeTime 2", "strokeTime 0", "strokeTime")
    variant("nameCollision", r"annulus\s*\n\{\n    ANN", "annulus\n{\n    WELL",
            "collides")
    variant("bopWithoutAnnulus", r"(?ms)^annulus\s*\{.*?^\}\n\n", "", "without annulus")
    variant("emptySegments", r"segments\s*\n\s*\{\n.*?^\s*\}\n",
            "segments\n        {\n        }\n", "segments")

    print()
    print("checks:", checks, "failures:", len(failures))
    if failures:
        for f in failures:
            print("FAILED:", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
