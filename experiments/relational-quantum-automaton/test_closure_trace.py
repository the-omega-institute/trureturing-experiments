"""Tests for the finite fixed-point closure trace."""
import unittest

from closure_trace import closure_trace


class RelationalQuantumAutomatonClosureTraceTests(unittest.TestCase):
    def test_trace_reaches_a_finite_fixed_point(self):
        result = closure_trace()
        self.assertTrue(result["fixed_point"])
        self.assertEqual(len(result["records"]), 6)
        self.assertEqual(result["rounds"][-1]["new_records"], 0)

    def test_rules_pass_in_dependency_order(self):
        result = closure_trace()
        rounds = result["rounds"]
        self.assertIn("R0_identity_norm", rounds[0]["passed"])
        self.assertNotIn("R1_binary_split", rounds[0]["passed"])
        self.assertIn("R1_binary_split", rounds[1]["passed"])
        self.assertIn("R2_oriented_ternary", rounds[2]["passed"])
        self.assertIn("R3_successor_proxy", rounds[3]["passed"])
        self.assertIn("R5_composition_proxy", rounds[3]["passed"])


if __name__ == "__main__":
    unittest.main()
