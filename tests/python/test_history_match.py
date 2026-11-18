"""D9 increment 1 regression: Python black-box steady-network history matching.

Run through CTest as `python.history_match`.  The harness must recover pipe
diameters from structured solver readback alone; no C++ source is changed and
no CLI text report is parsed.  Parameters and observations are checked against
the independent Hagen-Poiseuille series-network closed form.
"""

import math
import unittest

import chflow
from chflow import history_match as hm


TRUE_DIAMETERS = (0.055, 0.045, 0.050)
START_DIAMETERS = (0.040, 0.060, 0.042)


def closed_form_pressure(diameters, problem):
    """Independent laminar solution for the three-pipe manifold (Pa)."""
    d1, d2, d3 = diameters
    rho, mu, g = problem.rho, problem.mu, 9.80665
    q1 = problem.source_rates["W1"]
    q2 = problem.source_rates["W2"]
    r = lambda length, d: 128.0 * mu * length / (math.pi * d ** 4)
    p_m = problem.outlet_pressure + r(800.0, d3) * (q1 + q2)
    return {
        "W1": p_m + r(1000.0, d1) * q1 + rho * g * 5.0,
        "W2": p_m + r(1500.0, d2) * q2 + rho * g * 10.0,
        "M": p_m,
        "Outlet": problem.outlet_pressure,
    }


class TestHistoryMatch(unittest.TestCase):

    def setUp(self):
        self.problem = hm.CalibrationProblem()

    def test_structured_observation_matches_laminar_closed_form(self):
        net, result = hm.evaluate(TRUE_DIAMETERS, self.problem)
        self.assertTrue(result.converged)
        obs = hm.observation_from_result(net, result)
        self.assertEqual(obs.node_names, hm.NODE_NAMES)
        self.assertEqual(obs.pipe_names, tuple(p[0] for p in hm.PIPE_LAYOUT))
        expected = closed_form_pressure(TRUE_DIAMETERS, self.problem)
        for name, pressure in zip(obs.node_names, obs.node_pressure):
            self.assertLess(abs(pressure - expected[name]) / expected[name], 1e-8)
        self.assertLess(
            abs(obs.pipe_flow_rate[0] - self.problem.source_rates["W1"])
            / self.problem.source_rates["W1"], 1e-8)
        self.assertLess(
            abs(obs.pipe_flow_rate[2]
                - self.problem.source_rates["W1"] - self.problem.source_rates["W2"])
            / 0.003, 1e-8)

    def test_objective_is_zero_at_truth_and_positive_when_perturbed(self):
        obs = hm.synthetic_observation(TRUE_DIAMETERS, self.problem)
        self.assertLess(hm.objective(TRUE_DIAMETERS, obs, self.problem), 1e-18)
        self.assertGreater(
            hm.objective(START_DIAMETERS, obs, self.problem), 1.0)

    def test_fit_recovers_true_diameters_from_offset_start(self):
        obs = hm.synthetic_observation(TRUE_DIAMETERS, self.problem)
        fit = hm.fit(obs, START_DIAMETERS, self.problem)
        self.assertTrue(fit.converged)
        self.assertLess(fit.cost, 1e-16)
        self.assertLess(fit.initial_cost, 1.0e30)
        for got, want in zip(fit.params, TRUE_DIAMETERS):
            self.assertLess(abs(got - want) / want, 1e-4)
        self.assertLess(fit.iterations, 400)
        self.assertGreater(fit.evaluations, 0)
        # best-so-far residual history never grows
        self.assertTrue(all(b <= a + 1e-15 for a, b in zip(
            fit.residual_history, fit.residual_history[1:])))

    def test_fit_is_deterministic(self):
        obs = hm.synthetic_observation(TRUE_DIAMETERS, self.problem)
        a = hm.fit(obs, START_DIAMETERS, self.problem)
        b = hm.fit(obs, START_DIAMETERS, self.problem)
        self.assertEqual(a.params, b.params)
        self.assertEqual(a.cost, b.cost)
        self.assertEqual(a.evaluations, b.evaluations)

    def test_scales_rescale_objective_not_fit(self):
        obs = hm.synthetic_observation(TRUE_DIAMETERS, self.problem)
        tight = hm.objective(START_DIAMETERS, obs, self.problem, scales=(1.0e5, 1.0e-4))
        loose = hm.objective(START_DIAMETERS, obs, self.problem, scales=(1.0e6, 1.0e-3))
        self.assertAlmostEqual(tight, 100.0 * loose, delta=1e-9 * max(1.0, tight))
        fit = hm.fit(obs, START_DIAMETERS, self.problem,
                     scales=(1.0e5, 1.0e-4))
        for got, want in zip(fit.params, TRUE_DIAMETERS):
            self.assertLess(abs(got - want) / want, 1e-4)

    def test_invalid_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            hm.build_network((0.05, 0.05), self.problem)
        with self.assertRaises(ValueError):
            hm.build_network((0.05, -0.05, 0.05), self.problem)
        with self.assertRaises(ValueError):
            hm.objective(TRUE_DIAMETERS, hm.synthetic_observation(
                TRUE_DIAMETERS, self.problem), self.problem, scales=(0.0, 1.0))


if __name__ == "__main__":
    unittest.main()
