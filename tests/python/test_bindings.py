"""Python binding regression for the 1D kernels and the PVT layer
(M1 steady, M2 transient, shared src/pvt).

Run through CTest as `python.bindings` (PYTHONPATH is set by CMake to the
package dir and the chflow_core module dir). Compares against the exact
Hagen-Poiseuille closed form so the binding must round-trip every field.
"""
import math
import os
import unittest

import chflow


MU, RHO, D, G = 0.1, 850.0, 0.05, 9.80665
AREA = 0.25 * math.pi * D * D
CASE_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "cases")


def resistance(length, mu=MU, diameter=D):
    return 128.0 * mu * length / (math.pi * diameter ** 4)


class TestBindings(unittest.TestCase):

    def test_churchill_laminar_limit(self):
        self.assertLess(abs(chflow.churchill_friction(433.0, 0.0) - 64.0 / 433.0), 1e-9)

    def test_single_pipe_newton_matches_closed_form(self):
        net = chflow.OneDNetwork()
        net.fluid = chflow.FluidProperties(RHO, MU)
        net.add_node("W")
        net.add_node("Outlet")
        net.set_node_pressure("Outlet", 1.0e6)
        net.set_node_rate("W", 0.002)
        p = chflow.NetworkPipe()
        p.name, p.from_name, p.to_name = "T", "W", "Outlet"
        p.length, p.diameter, p.roughness, p.elevation = 1000.0, D, 0.0, 5.0
        net.add_pipe(p)

        res = chflow.solve_steady_network(net, chflow.SteadyControls())
        self.assertTrue(res.converged)
        self.assertAlmostEqual(res.node_pressure[0],
                               1.0e6 + resistance(1000.0) * 0.002 + RHO * G * 5.0,
                               delta=1e-6 * 2.3e6)
        self.assertAlmostEqual(res.pipe_flows[0].flow_rate, 0.002, delta=2e-9)
        self.assertAlmostEqual(res.pipe_flows[0].velocity * AREA,
                               res.pipe_flows[0].flow_rate, delta=1e-12)
        self.assertEqual(res.pipe_flows[0].name, "T")

    def test_error_translates_to_python_exception(self):
        net = chflow.OneDNetwork()
        net.fluid = chflow.FluidProperties(RHO, MU)
        net.add_node("A")
        with self.assertRaises(RuntimeError):
            net.node_index("missing")

    def test_case_dict_round_trip(self):
        net = chflow.load_network(os.path.join(CASE_DIR, "steadyNetwork", "case.dict"))
        ctl = chflow.load_steady_controls(os.path.join(CASE_DIR, "steadyNetwork", "case.dict"))
        self.assertEqual(net.n_nodes(), 4)
        self.assertEqual(len(net.pipes), 3)
        res = chflow.solve_steady_network(net, ctl)
        self.assertTrue(res.converged)
        pM_expected = 1.0e6 + resistance(800.0) * 0.003
        pM = res.node_pressure[net.node_index("M")]
        self.assertLess(abs(pM - pM_expected) / pM_expected, 1e-6)
        q3 = next(f for f in res.pipe_flows if f.name == "P3").flow_rate
        self.assertLess(abs(q3 - 0.003) / 0.003, 1e-6)


class TestTransientBindings(unittest.TestCase):

    BLOWDOWN = os.path.join(CASE_DIR, "transientBlowdown", "case.dict")

    def test_hydrate_kinetics_is_the_exact_exponential(self):
        vL, vH = chflow.hydrate_kinetics_step(1.0, 0.0, 0.5, 0.3)
        self.assertAlmostEqual(vL, math.exp(-0.15), delta=1e-15)
        self.assertAlmostEqual(vL + vH, 1.0, delta=1e-15)

    def test_blowdown_series_matches_the_rc_envelope(self):
        case = chflow.load_transient_case(self.BLOWDOWN)
        ctl = chflow.load_transient_controls(self.BLOWDOWN)
        self.assertEqual(len(case.nodes), 2)
        self.assertEqual(len(case.pipes), 1)
        self.assertEqual(ctl.dt, 0.05)
        series = chflow.run_transient(case, ctl.dt, ctl.end_time, ctl.report_every)
        self.assertEqual(len(series.samples), int(ctl.end_time / ctl.dt) + 1)

        tank = series.node_names.index("TANK")
        line = series.pipe_names.index("LINE")
        ps = [s.pressure[tank] for s in series.samples]
        self.assertTrue(all(b < a for a, b in zip(ps, ps[1:])))
        tau = resistance(1000.0, MU, 0.02) * 1.0 / 2.0e8
        envelope = 5.0e5 * math.exp(-series.samples[-1].t / tau)
        self.assertLess(abs(ps[-1] - envelope) / envelope, 0.02)
        self.assertEqual(len(series.samples[0].flow_rate), len(series.pipe_names))
        self.assertGreater(series.samples[-1].flow_rate[line], 0.0)
        # the tank drains, so its mobile-liquid fraction must fall
        self.assertLess(series.samples[-1].liquid_fraction[tank], 1.0)

    def test_steppable_transient_network_matches_batch_run(self):
        case = chflow.load_transient_case(self.BLOWDOWN)
        ctl = chflow.load_transient_controls(self.BLOWDOWN)
        series = chflow.run_transient(case, ctl.dt, ctl.end_time, ctl.report_every)

        sim = chflow.TransientNetwork(case)
        steps = int(round(ctl.end_time / ctl.dt))
        for _ in range(steps):
            sim.step(ctl.dt)
        self.assertAlmostEqual(sim.time(), ctl.end_time, delta=1e-12)
        self.assertEqual(sim.node_names(), series.node_names)
        self.assertEqual(sim.pipe_names(), series.pipe_names)
        last = series.samples[-1]
        for a, b in zip(sim.pressure(), last.pressure):
            self.assertTrue(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12))
        for a, b in zip(sim.pipe_flows(), last.flow_rate):
            self.assertTrue(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12))
        for a, b in zip(sim.liquid_fractions(), last.liquid_fraction):
            self.assertTrue(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12))
        self.assertEqual(sum(sim.hydrate_volumes()), 0.0)


class TestPvtBindings(unittest.TestCase):
    """README: src/pvt is the shared property layer - check the binding can
    reproduce the EOS invariants the C++ regression tests assert."""

    CASE = os.path.join(CASE_DIR, "pvtFlash", "case.dict")

    def feed(self):
        mix = chflow.Mixture()
        for name, frac in (("methane", 0.34), ("ethane", 0.18),
                           ("propane", 0.26), ("nPentane", 0.22)):
            mix.add(chflow.component(name), frac)
        return mix

    def test_unknown_component_raises(self):
        with self.assertRaises(RuntimeError):
            chflow.component("unobtainium")

    def test_critical_point_compressibility(self):
        c = chflow.component("methane")
        st = chflow.peng_robinson_state(self.single(c), c.tc, c.pc,
                                        chflow.RootPreference.lowest_gibbs)
        self.assertLess(abs(st.z - 0.30740) / 0.30740, 0.02)

    def single(self, component):
        mix = chflow.Mixture()
        mix.add(component, 1.0)
        return mix

    def test_pure_fluid_bubble_point_is_its_vapor_pressure(self):
        ethane = chflow.component("ethane")
        psat = chflow.pure_vapor_pressure(ethane, 250.0)
        self.assertGreater(psat, 1.0e5)
        self.assertLess(psat, ethane.pc)
        self.assertAlmostEqual(chflow.bubble_pressure(self.single(ethane), 250.0),
                               psat, delta=1e-5 * psat)
        self.assertEqual(chflow.pure_vapor_pressure(ethane, ethane.tc + 10.0), 0.0)

    def test_flash_satisfies_fugacity_equality(self):
        feed = self.feed()
        r = chflow.tp_flash(feed, 320.0, 4.0e6)
        self.assertTrue(r.converged)
        self.assertFalse(r.single_phase)
        self.assertLess(0.0, r.vapor_fraction)
        self.assertGreater(1.0, r.vapor_fraction)
        self.assertGreater(r.liquid_mass_density, r.vapour_mass_density)
        z = feed.fractions()
        for i, (x, y) in enumerate(zip(r.liquid_composition, r.vapour_composition)):
            self.assertAlmostEqual(r.vapor_fraction * y + (1.0 - r.vapor_fraction) * x,
                                   z[i], delta=1e-9)
            # f_i^L = f_i^V at equilibrium
            self.assertAlmostEqual(math.log(x) + r.liquid.ln_phi[i],
                                   math.log(y) + r.vapour.ln_phi[i], delta=1e-8)
        # beyond the bubble point the same feed collapses to one liquid phase
        self.assertTrue(chflow.tp_flash(feed, 320.0, 1.2e7).single_phase)

    def test_case_dict_and_bubble_point_bracket_the_envelope(self):
        case = chflow.load_flash_case(self.CASE)
        self.assertEqual(case.feed.size(), 4)
        self.assertEqual(len(case.pressures), 3)
        self.assertTrue(case.compute_bubble_point)
        pb = chflow.bubble_pressure(case.feed, case.temperature)
        self.assertLess(8.0e6, pb)
        self.assertGreater(9.0e6, pb)
        self.assertFalse(chflow.tp_flash(case.feed, case.temperature, 0.99 * pb).single_phase)
        self.assertTrue(chflow.tp_flash(case.feed, case.temperature, 1.01 * pb).single_phase)

    def test_table_matches_the_direct_evaluation(self):
        feed = self.feed()
        pGrid, tGrid = [2.0e6, 4.0e6, 6.0e6], [310.0, 320.0, 330.0]
        table = chflow.liquid_density_table(feed, pGrid, tGrid)
        direct = chflow.tp_flash(feed, 320.0, 4.0e6).liquid_mass_density
        self.assertAlmostEqual(table.sample(4.0e6, 320.0), direct, delta=1e-9 * direct)
        mid = chflow.tp_flash(feed, 320.0, 3.0e6).liquid_mass_density
        self.assertLess(abs(table.sample(3.0e6, 320.0) - mid) / mid, 0.02)
        self.assertEqual(table.size(), 9)


class TestNodalBindings(unittest.TestCase):
    """The IPR-VLP nodal kernel exported to Python: closed-form curve access,
    a case-dict round trip and parity with the app's golden anchor."""

    CASE = os.path.join(CASE_DIR, "nodalVogel", "case.dict")
    GOLDEN_Q = 0.00290785579167

    def test_ipr_curve_is_the_published_closed_form(self):
        ipr = chflow.IpR()
        ipr.form = chflow.IpRForm.vogel
        ipr.reservoir_pressure, ipr.max_rate = 2.2e7, 0.006
        self.assertAlmostEqual(ipr.rate(0.5 * 2.2e7), 0.7 * 0.006, delta=1e-12)
        self.assertAlmostEqual(ipr.p_wf(0.7 * 0.006), 1.1e7, delta=1e-3)
        self.assertAlmostEqual(ipr.rate_max(), 0.006, delta=1e-15)

    def test_case_round_trip_matches_the_app_golden_anchor(self):
        case = chflow.load_nodal_case(self.CASE)
        self.assertEqual(case.ipr.form, chflow.IpRForm.vogel)
        self.assertAlmostEqual(case.wellhead_pressure, 2.5e6, delta=1e-3)
        r = chflow.solve_nodal_analysis(case.ipr, case.fluid, case.well,
                                        case.wellhead_pressure)
        self.assertTrue(r.flows and r.converged)
        self.assertLess(abs(r.rate - self.GOLDEN_Q) / self.GOLDEN_Q, 1e-9)
        self.assertLess(r.bottom_hole_pressure, r.static_bottom_pressure)
        self.assertGreater(r.holdup, 0.0)
        self.assertLess(r.holdup, 1.0)

    def test_vlp_inversion_lands_on_the_wellhead(self):
        case = chflow.load_nodal_case(self.CASE)
        pwf = chflow.vlp_bottom_hole_pressure(case.fluid, case.well, 0.0,
                                              case.wellhead_pressure)
        self.assertGreater(pwf, 0.8 * 850.0 * 9.80665 * 2200.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
