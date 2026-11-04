#!/usr/bin/env python3
"""Run the wall-function channel with both RANS models and verify the friction
physics through either pressure-velocity coupling.

The case drives a streamwise-uniform column with a body force, so the wall
shear is not a prediction but the momentum balance: tau_w = u_tau^2 =
force * half-height. What the turbulence model has to deliver is the bulk
velocity that balance can sustain, and the shape of the near-wall profile.

The same assertions run on the transient (chflow-piso-foam, which reaches the
balance over 2000 steps) and on the steady loop (chflow-simple-foam, which
reaches it in a few dozen iterations); which one is in play is read off the
log. The profile always comes from the vtu cell centroids, so the two couplings
are compared on the same numbers.

The case also carries an energy block with turbulentDissipation on, so the
heat the turbulence releases through eps/cp is part of what is checked: the
discrete energy balance, the source against the volume integral of eps, and
the source against the pump power in the band the wall-function closure
leaves (assertion 8).

usage: check_turbulent_channel.py <solver-binary> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

KAPPA = 0.41
SPALDING_C = 5.5676
CMU = 0.09
BETA_STAR = 0.09


def case_constants(path):
    """The numbers the checks are built from, read out of the case itself so it
    cannot drift away from them."""
    txt = open(path).read()
    nu = float(re.search(r'^nu\s+([0-9.eE+-]+);', txt, re.M).group(1))
    maxv = [float(v) for v in
            re.search(r'max\s*\(([^)]*)\)', txt).group(1).split()]
    cells = [int(v) for v in re.search(r'cells\s*\(([^)]*)\)', txt).group(1).split()]
    force = float(re.search(r'force\s*\(([^)]*)\)', txt).group(1).split()[0])
    end_time = float(re.search(r'endTime\s+([0-9.eE+-]+);', txt).group(1))
    if cells[0] != 1 or cells[2] != 1:
        raise SystemExit("the channel has to be a single column in x and z")
    return dict(nu=nu, height=maxv[1], width=maxv[0], depth=maxv[2],
                force=force, ny=cells[1], end_time=end_time)


def spalding_u_plus(y_plus):
    """Invert Spalding's law of the wall, y+ = u+ + e^(-kC)[e^(ku+) - 1 - ku+
    - (ku+)^2/2 - (ku+)^3/6], by bisection on the monotone relation."""
    lo, hi = 0.0, 50.0
    a = math.exp(-KAPPA * SPALDING_C)

    def yp(up):
        ku = KAPPA * up
        return up + a * (math.exp(ku) - 1.0 - ku - ku * ku / 2.0 - ku ** 3 / 6.0)

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if yp(mid) < y_plus:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def pope_bulk_velocity(u_tau, d_h, nu):
    """Pope's smooth-wall friction law, 1/sqrt(f) = 2 log10(Re_b sqrt(f)) - 0.6.
    With u_tau = U_b sqrt(f/8) the product Re_b sqrt(f) is fixed by u_tau alone,
    so the law gives f - and then U_b - directly."""
    inv_sqrt_f = 2.0 * math.log10(d_h * u_tau * math.sqrt(8.0) / nu) - 0.6
    f = 1.0 / inv_sqrt_f ** 2
    return u_tau * math.sqrt(8.0 / f)


def reports_from_log(out):
    return [(float(m.group(1)), float(m.group(2))) for m in re.finditer(
        r'^\[chflow\] t ([0-9.eE+-]+) uMax ([0-9.eE+-]+) cont ([0-9.eE+-]+)', out, re.M)]


def profile_from_vtu(vtu_path):
    """Cell-centre (y, u) sorted across the column. The centroids come from the
    mesh vertices the vtu stores, so the profile is the same field whether the
    run was a transient or a steady loop -- and no log parsing is involved."""
    txt = open(vtu_path).read()
    pts = re.search(r'<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
                    r'\s*(.*?)</DataArray>', txt, re.S).group(1)
    c = [float(v) for v in pts.split()]
    points = list(zip(c[0::3], c[1::3], c[2::3]))
    conn = [int(v) for v in re.search(r'Name="connectivity" format="ascii">'
                                      r'\s*(.*?)</DataArray>', txt, re.S).group(1).split()]
    offs = [int(v) for v in re.search(r'Name="offsets" format="ascii">'
                                      r'\s*(.*?)</DataArray>', txt, re.S).group(1).split()]
    u = cell_vector(vtu_path, "U")
    prof, prev = [], 0
    for o in offs:
        vids = conn[prev:o]
        prev = o
        prof.append((sum(points[v][1] for v in vids) / len(vids), u[len(prof)]))
    prof.sort()
    return prof


def cell_scalars(vtu_path, name):
    txt = open(vtu_path).read()
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii">\s*(.*?)</DataArray>'
                  % name, txt, re.S)
    if not m:
        raise SystemExit("%s: field %s missing from the vtu" % (vtu_path, name))
    return [float(v) for v in m.group(1).split()]


def cell_vector(vtu_path, name):
    txt = open(vtu_path).read()
    m = re.search(r'<DataArray[^>]*Name="%s" NumberOfComponents="3"[^>]*format="ascii">'
                  r'\s*(.*?)</DataArray>' % name, txt, re.S)
    v = [float(x) for x in m.group(1).split()]
    return v[0::3]


def run(binary, case, model):
    """One model's run, in a copy of the case directory: the solver writes its
    vtu next to the case, and the repository copy stays clean."""
    work = tempfile.mkdtemp(prefix="ch-turbulent-channel-")
    binary = os.path.abspath(binary)
    src = os.path.dirname(case)
    for f in os.listdir(src):
        shutil.copy(os.path.join(src, f), work)
    local = os.path.join(work, os.path.basename(case))
    txt = open(local).read()
    open(local, "w").write(txt.replace("model       kEpsilon;", "model       %s;" % model))
    r = subprocess.run([binary, local], capture_output=True, text=True, cwd=work)
    print("---- %s (rc %d) ----" % (model, r.returncode))
    print(r.stdout[-1200:])
    if r.returncode != 0:
        print(r.stderr)
        raise SystemExit("%s solver failed" % model)
    if "nan" in r.stdout.lower():
        raise SystemExit("%s: the run produced a non-finite field" % model)
    vtu = os.path.join(work, os.path.basename(re.search(
        r'^\s*vtk\s+(\S+);', open(local).read(), re.M).group(1)))
    return dict(out=r.stdout, vtu=vtu, work=work)


def check(model, res, c):
    out = res["out"]
    half = c["height"] / 2.0
    u_tau_exact = math.sqrt(c["force"] * half)
    failures = []

    wall = re.search(r'turbulence: u_tau=([0-9.eE+-]+) y\+ mean=([0-9.eE+-]+)'
                     r' max=([0-9.eE+-]+)', out)
    if wall is None:
        raise SystemExit("%s: no wall report in the log" % model)
    u_tau, yp_mean, yp_max = (float(wall.group(i)) for i in (1, 2, 3))

    # 1. the integral balance: the wall function prescribes the shear, and the
    #    solved field has to agree with the driving force
    bal = abs(u_tau / u_tau_exact - 1.0)
    print("%s: u_tau=%.6f (balance %.6f, err %.1e)" % (model, u_tau, u_tau_exact, bal))
    if bal > 0.01:
        failures.append("wall shear does not balance the driving force")

    # 2. the mesh has to sit in the band the wall function is derived for
    print("%s: y+ = %.1f..%.1f" % (model, yp_mean, yp_max))
    if not (30.0 < yp_mean < 300.0):
        failures.append("mean y+ outside the wall-function band (30..300)")

    # 3. steady. Each coupling says it differently: the transient compares its
    #    last two reports, SIMPLE reports the fixed-point change of the velocity
    #    over the last outer iteration -- the only signal that means anything on
    #    a streamwise-uniform column, where every phi and every p' is zero.
    rep = reports_from_log(out)
    if rep:
        if len(rep) < 2:
            raise SystemExit("%s: too few reports to judge steadiness" % model)
        drift = abs(rep[-1][1] / rep[-2][1] - 1.0)
        print("%s: uMax %.6f at t=%.3g, drift over the last interval %.2e"
              % (model, rep[-1][1], rep[-1][0], drift))
        if drift > 0.01:
            failures.append("the profile is still moving at endTime")
    else:
        done = re.search(r"SIMPLE done: iterations=(\d+) continuity=([0-9.eE+-]+)"
                         r" momentum=([0-9.eE+-]+) \((converged|NOT converged)\)", out)
        if done is None:
            raise SystemExit("%s: neither a transient report nor a SIMPLE "
                             "convergence line in the log" % model)
        iters, res_u, state = int(done.group(1)), float(done.group(3)), done.group(4)
        print("%s: SIMPLE %s in %d iterations, velocity residual %.2e"
              % (model, state, iters, res_u))
        if state != "converged":
            failures.append("SIMPLE did not reach its convergence criteria")
        if res_u > 1e-4:
            failures.append("the velocity field is still moving")
        # one iteration means the stop fired on continuity alone, which is what
        # this case is identically zero at; the wall shear cannot be in balance
        # yet however small the number reads
        if iters < 2:
            failures.append("SIMPLE stopped without letting the momentum evolve")

    # 4. profile shape: monotone to the centreline, symmetric about it
    prof = profile_from_vtu(res["vtu"])
    if len(prof) != c["ny"]:
        raise SystemExit("%s: expected %d profile points, got %d"
                         % (model, c["ny"], len(prof)))
    n = len(prof) // 2
    lower = [u for _, u in prof[:n]]
    upper = [u for _, u in prof[n:]][::-1]
    monotone = all(b > a for a, b in zip(lower, lower[1:]))
    sym = max(abs(a - b) / b for a, b in zip(lower, upper))
    print("%s: Ub=%.4f Umax=%.4f symmetry=%.2e monotone=%s"
          % (model, sum(u for _, u in prof) / len(prof), max(lower), sym, monotone))
    if not monotone:
        failures.append("the profile is not monotone from the wall to the centreline")
    if sym > 1e-6:
        failures.append("the profile is not symmetric about the centreline")

    # 5. the law of the wall at the first cell centre, which is where the wall
    #    function closes the momentum equation
    y1 = prof[0][0]
    yp1 = u_tau * y1 / c["nu"]
    up1 = prof[0][1] / u_tau
    up_law = spalding_u_plus(yp1)
    err = abs(up1 / up_law - 1.0)
    print("%s: first cell y+=%.1f u+=%.3f (Spalding %.3f, err %.1e)"
          % (model, yp1, up1, up_law, err))
    if err > 0.02:
        failures.append("first-cell u+ does not follow the law of the wall")

    # 6. the friction law, the actual prediction of the case: the bulk velocity
    #    this shear can sustain. The mesh puts 8 cells between y+ = 52 and the
    #    centreline, so the band is the mesh's, not the model's.
    ub = sum(u for _, u in prof) / len(prof)
    expect = pope_bulk_velocity(u_tau, 2.0 * c["height"], c["nu"])
    err_f = ub / expect - 1.0
    print("%s: Ub/u_tau = %.2f, Pope: %.2f (err %+.1f%%)"
          % (model, ub / u_tau, expect / u_tau, 100.0 * err_f))
    if abs(err_f) > 0.2:
        failures.append("bulk velocity outside the friction-law band")

    # 7. the equilibrium wall cell, out of the vtu: the closure that keeps the
    #    near-wall transport off the wall function
    diss_name = "epsilon" if model == "kEpsilon" else "omega"
    k = cell_scalars(res["vtu"], "k")
    nut = cell_scalars(res["vtu"], "nut")
    diss = cell_scalars(res["vtu"], diss_name)
    sqrt_const = math.sqrt(CMU if model == "kEpsilon" else BETA_STAR)
    k_eq = u_tau * u_tau / sqrt_const
    nut_log = KAPPA * u_tau * y1
    # u_tau enters here through the wall report, whose printed digits bound how
    # closely the closure can be re-derived (3e-6 relative on u_tau, so 1e-4 on
    # k and 3e-4 on a cubic of it). A closure that was not applied at all leaves
    # the cell orders of magnitude away, not at this level.
    tol_pair = 5e-5
    if model == "kEpsilon":
        d_eq = u_tau ** 3 / (KAPPA * y1)
        ok_k = abs(k[0] / k_eq - 1.0) < tol_pair
        ok_d = abs(diss[0] / d_eq - 1.0) < tol_pair
        ok_nut = abs(nut[0] / nut_log - 1.0) < 1e-3
    else:
        d_eq = u_tau / (sqrt_const * KAPPA * y1)
        ok_k = abs(k[0] / k_eq - 1.0) < tol_pair
        ok_d = abs(diss[0] / d_eq - 1.0) < tol_pair
        # the SST shear branch may only reduce the log-layer viscosity
        ok_nut = 0.0 < nut[0] <= nut_log
    print("%s: wall cell k=%.4f (equilibrium %.4f), %s=%.1f (%.1f), nut=%.3e (log layer %.3e)"
          % (model, k[0], k_eq, diss_name, diss[0], d_eq, nut[0], nut_log))
    if not ok_k:
        failures.append("wall cell k is not the equilibrium value")
    if not ok_d:
        failures.append("wall cell dissipation is not the equilibrium value")
    if not ok_nut:
        failures.append("wall cell eddy viscosity is not the log-layer value")

    # 8. the heat the turbulence releases. turbulentDissipation puts eps/cp
    #    into the energy equation and the 300 K walls are the only exit, so
    #    the discrete balance closes with the wall flux carrying exactly the
    #    generated heat -- the residual that remains on the transient is the
    #    storage rate. The source itself is the volume integral of eps (the
    #    vtu stores it, through betaStar k omega for SST), which is what the
    #    k-equation actually transports. Against the pump power force*Ub*V
    #    the ratio only sits in a band: the wall-function closure holds the
    #    wall cells algebraically at the log-layer pair, which is not
    #    energy-conservative -- measured 0.76..1.18 over the four runs of
    #    this case. Without the eps source the ratio reads 0.004 (kEpsilon)
    #    to 0.03 (SST), far outside the band, so a missing or doubled source
    #    still shows.
    energy = re.search(r'energy: in (\S+) source (\S+) enthalpyOut (\S+)'
                       r' balance (\S+)', out)
    if energy is None:
        raise SystemExit("%s: no energy balance line in the log" % model)
    q_in, src, h_out, bal = (float(energy.group(i)) for i in (1, 2, 3, 4))
    heat = re.findall(r'heatPatch (\w+) Q (\S+)', out)
    if not heat:
        raise SystemExit("%s: no heatPatch report in the log" % model)
    eps_cell = cell_scalars(res["vtu"], "epsilon" if model == "kEpsilon"
                            else "omega")
    if model == "kOmegaSST":
        eps_cell = [BETA_STAR * kk * ww for kk, ww in zip(k, eps_cell)]
    v_cell = c["width"] * (c["height"] / c["ny"]) * c["depth"]
    src_vtu = v_cell * sum(eps_cell)
    pump = c["force"] * ub * c["width"] * c["height"] * c["depth"]
    ratio = src / pump
    print("%s: heat in %.6e over %d patches, source %.6e (integral of eps "
          "%.6e), balance %.1e; source/pump %.4f"
          % (model, q_in, len(heat), src, src_vtu, bal, ratio))
    if abs(bal) > 1e-4 * abs(src):
        failures.append("the energy balance does not close (storage left)")
    if abs(abs(q_in) / src - 1.0) > 1e-4:
        failures.append("the walls do not carry the generated heat out")
    if abs(h_out) > 1e-3 * abs(src):
        failures.append("enthalpy leaves through the ends of a uniform column")
    if abs(src_vtu / src - 1.0) > 1e-6:
        failures.append("the reported source is not the integral of eps")
    if not 0.5 < ratio < 2.0:
        failures.append("the heat source is not of the pump-power order")

    return ub, failures


def main():
    binary, case = sys.argv[1], sys.argv[2]
    c = case_constants(case)
    print("channel: height %.4g, nu %.3g, force %.4g -> u_tau %.6f, Re_tau %.1f"
          % (c["height"], c["nu"], c["force"], math.sqrt(c["force"] * c["height"] / 2.0),
             math.sqrt(c["force"] * c["height"] / 2.0) * c["height"] / 2.0 / c["nu"]))

    results = {}
    for model in ("kEpsilon", "kOmegaSST"):
        res = run(binary, case, model)
        try:
            results[model] = check(model, res, c)
        finally:
            shutil.rmtree(res["work"], ignore_errors=True)

    bad = False
    for model, (ub, failures) in results.items():
        for f in failures:
            print("FAIL %s: %s" % (model, f))
            bad = True
    if bad:
        raise SystemExit("turbulent channel verification FAILED")
    print("turbulent channel verification passed: %s"
          % ", ".join("%s Ub=%.3f" % (m, ub) for m, (ub, _) in results.items()))


if __name__ == "__main__":
    main()
