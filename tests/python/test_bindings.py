"""Python binding regression for the 1D steady network kernel (M1).

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
        case = os.path.join(os.path.dirname(__file__), os.pardir,
                            "cases", "steadyNetwork", "case.dict")
        net = chflow.load_network(case)
        ctl = chflow.load_steady_controls(case)
        self.assertEqual(net.n_nodes(), 4)
        self.assertEqual(len(net.pipes), 3)
        res = chflow.solve_steady_network(net, ctl)
        self.assertTrue(res.converged)
        pM_expected = 1.0e6 + resistance(800.0) * 0.003
        pM = res.node_pressure[net.node_index("M")]
        self.assertLess(abs(pM - pM_expected) / pM_expected, 1e-6)
        q3 = next(f for f in res.pipe_flows if f.name == "P3").flow_rate
        self.assertLess(abs(q3 - 0.003) / 0.003, 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
