"""Tests for the finite relational quantum automaton certificate."""
import unittest

from verify import GATES, certificate, identity_checks, matrix_multiply


class RelationalQuantumAutomatonTests(unittest.TestCase):
    def test_gate_identities(self):
        self.assertTrue(all(identity_checks().values()))

    def test_projective_orbit_and_closure_checks(self):
        result = certificate()
        self.assertEqual(result["summary"]["projective_orbit_size"], 6)
        self.assertTrue(all(result["checks"].values()))

    def test_hadamard_is_its_own_inverse(self):
        square = matrix_multiply(GATES["H"], GATES["H"])
        self.assertAlmostEqual(square[0][0].real, 1.0)
        self.assertAlmostEqual(square[1][1].real, 1.0)
        self.assertAlmostEqual(square[0][1].real, 0.0)
        self.assertAlmostEqual(square[1][0].real, 0.0)


if __name__ == "__main__":
    unittest.main()
