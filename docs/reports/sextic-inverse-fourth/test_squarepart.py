#!/usr/bin/env python3
"""Independent ordered-four-factor/all-unit checks for the new split."""
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import unittest
from eisenstein import Engine, UNITS, mul, norm
from squarepart import run_scale


def literal_small_P(engine, U):
    """Visit original ordered factor slots and every physical unit row.

    No grouped inverse plan, unit orthogonality, pair plan or trace array
    from squarepart.py participates in this independent calculation.
    """
    D = math.sqrt(U)
    columns = [d for d in engine.ideals
               if D/2 < norm(d) < 1.5*D and engine.mu(d)]
    tuples = []
    for indices in itertools.product(range(len(columns)), repeat=4):
        occupancy = Counter()
        sign = 1
        for i in indices:
            occupancy.update(engine._factor(columns[i]).keys())
            sign *= engine.mu(columns[i])
        if not any(v % 2 for v in occupancy.values()):
            continue
        root_squarepart = math.prod(norm(p)**(v//2) for p, v in occupancy.items())
        if root_squarepart**4 < U:
            tuples.append((*indices, sign))
    # Build the six powers by exact ring multiplication (a+b*zeta)*zeta.
    powers, a, b = [], 1, 0
    for _ in range(6):
        powers.append((a, b))
        a, b = -b, a+b
    real, other = 0, 0
    physical_rows = 0
    for u in engine.ideals:
        if not U < norm(u) < 2*U or not engine.is_sixth_free(u):
            continue
        for unit in UNITS:
            row = mul(unit, u)
            exponents = [engine.symbol_exponent(d, row) for d in columns]
            physical_rows += 1
            for i, j, k, l, sign in tuples:
                x, y, z, w = [exponents[t] for t in (i, j, k, l)]
                if min(x, y, z, w) < 0:
                    continue
                a, b = powers[(x+y-z-w) % 6]
                real += sign*a
                other += sign*b
    if other != 0:
        raise ArithmeticError('ordered Hermitian sum has nonreal coefficient')
    return Fraction(real, U), physical_rows


class SmallSquarepartTests(unittest.TestCase):
    def test_positive_value_by_literal_original_factor_expansion(self):
        engine = Engine(4096)
        expected, rows = literal_small_P(engine, 2048)
        actual = run_scale(engine, 2048)
        q = actual['flat']['small_P']
        self.assertEqual(expected, Fraction(747, 64))
        self.assertEqual(Fraction(q['numerator'], q['denominator']), expected)
        self.assertEqual(actual['physical_rows'], rows)

    def test_fixed_bump_matches_existing_full_moment_baseline(self):
        source = Path(__file__).with_name('pilot.json')
        baseline = json.loads(source.read_text())['results']
        engine = Engine(8192)
        for U in [256, 1024, 4096]:
            actual = run_scale(engine, U)
            original = next(r for r in baseline if r['U'] == U)
            for new, old in [('fourth', 'mobius_fourth'), ('h2', 'h2'), ('e6', 'e6')]:
                self.assertAlmostEqual(actual['bump'][new], original[old], delta=2e-10)
            self.assertEqual(actual['physical_rows'], original['physical_rows'])

    def test_retained_exact_split_and_fixed_bump_sign_change(self):
        report = json.loads(Path(__file__).with_name('squarepart.json').read_text())
        values = {}
        for row in report['results']:
            def f(key):
                q = row['flat'][key]
                return Fraction(q['numerator'], q['denominator'])
            self.assertEqual(f('fourth'), f('h2')+f('small_P')+f('large_P'))
            values[row['U']] = row
        self.assertLess(values[128]['bump']['small_P'], 0)
        self.assertGreater(values[4096]['bump']['small_P'], 0)


if __name__ == '__main__':
    unittest.main()
