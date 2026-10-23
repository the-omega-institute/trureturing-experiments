"""Independent finite checks of the exact Eisenstein arithmetic backend."""
import unittest
from eisenstein import Engine, UNITS, div_exact, mul, norm, primary, _modmul, _modpow


def power(x, k):
    # Repeated multiplication in the integral ring, not modular powering.
    y = (1, 0)
    for _ in range(k):
        y = mul(y, x)
    return y


def sub(x, y):
    return x[0]-y[0], x[1]-y[1]


class ArithmeticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = Engine(1000)

    def test_primary_units_and_complete_enumeration(self):
        self.assertEqual(tuple(power((1, 1), k) for k in range(6)), UNITS)
        for a in range(-12, 13):
            for b in range(-12, 13):
                x = a, b
                if norm(x) % 3:
                    associates = [mul(u, x) for u in UNITS]
                    admissible = [y for y in associates if y[0] % 3 == 1 and y[1] % 3 == 0]
                    self.assertEqual(admissible, [primary(x)])
        # An independent square enumeration is complete for norms <= 100.
        expected = {primary((a, b)) for a in range(-20, 21) for b in range(-20, 21)
                    if 0 < norm((a, b)) <= 100 and norm((a, b)) % 6 == 1}
        self.assertEqual(set(Engine(100).ideals), expected)
        self.assertEqual(self.engine.ideals[0], (1, 0))

    def test_factor_reconstruction_and_mobius(self):
        e = self.engine
        for x in e.ideals:
            y = (1, 0)
            for p, k in e.factors[x].items():
                y = mul(y, power(p, k))
                self.assertEqual(e.factors[p], {p: 1})
            self.assertEqual(y, x)
            self.assertEqual(e.mu(x), 0 if any(k > 1 for k in e.factors[x].values())
                             else (-1)**len(e.factors[x]))
        for p in [(-2, -3), (5, 0), (-5, -3)]:
            p = primary(p)
            self.assertEqual(e.mu(power(p, 2)), 0)
            self.assertEqual(e.liouville(power(p, 5)), -1)
            self.assertTrue(e.is_sixth_free(power(p, 5)))
            self.assertFalse(e.is_sixth_free(power(p, 6)))
        with self.assertRaises(ValueError):
            e.mu((2, 0))
        with self.assertRaises(ValueError):
            primary((1, -1))
        self.assertIsNone(div_exact((1, 0), (2, 0)))

    def test_all_small_quotient_fields_and_character_values(self):
        e = self.engine
        primes = [x for x in e.ideals if (norm(x) <= 49 or norm(x) == 121) and e.factors[x] == {x: 1}]
        self.assertEqual(sorted(norm(x) for x in primes), [7, 7, 13, 13, 19, 19, 25, 31, 31, 37, 37, 43, 43, 121])
        for prime in primes:
            p, root, _ = e._field(prime)
            elements = [(a, 0) for a in range(p)] if root is not None else [(a, b) for a in range(p) for b in range(p)]
            self.assertEqual(len(elements), norm(prime))
            counts = [0]*6
            for x in elements:
                k = e.prime_symbol_exponent(prime, x)
                self.assertEqual(k == -1, div_exact(x, prime) is not None)
                if k >= 0:
                    counts[k] += 1
                    direct = power(x, (norm(prime)-1) // 6)
                    self.assertIsNotNone(div_exact(sub(direct, UNITS[k]), prime))
                for y in elements:
                    product = _modmul(x, y, p)
                    self.assertIsNotNone(div_exact(sub(mul(x, y), product), prime))
                for k in range(8):
                    self.assertIsNotNone(div_exact(sub(power(x, k), _modpow(x, k, p)), prime))
            self.assertEqual(counts, [(norm(prime)-1) // 6]*6)
            # Check the split omega map against integral ideal divisibility.
            for x in [(a, b) for a in range(p) for b in range(p)]:
                residue = e._residue_raw(x, p, root)
                self.assertIsNotNone(div_exact(sub(x, residue), prime))

    def test_collision_multiplicativity_and_unit_weight(self):
        e = self.engine
        ideals = [x for x in e.ideals if norm(x) <= 50]
        for d in ideals:
            for k, unit in enumerate(UNITS):
                self.assertEqual(e.symbol_exponent(d, unit), k*e.unit_weight(d) % 6)
            for p in e.factors[d]:
                self.assertEqual(e.symbol_exponent(d, mul(p, (2, 1))), -1)
            for x in ideals:
                for y in ideals:
                    a, b = e.symbol_exponent(d, x), e.symbol_exponent(d, y)
                    expected = -1 if -1 in (a, b) else (a+b) % 6
                    self.assertEqual(e.symbol_exponent(d, mul(x, y)), expected)
                for other in ideals:
                    a, b = e.symbol_exponent(d, x), e.symbol_exponent(other, x)
                    expected = -1 if -1 in (a, b) else (a+b) % 6
                    self.assertEqual(e.symbol_exponent(mul(d, other), x), expected)
        with self.assertRaises(ValueError):
            e.prime_symbol_exponent((1, 0), (1, 0))


if __name__ == "__main__":
    unittest.main()
