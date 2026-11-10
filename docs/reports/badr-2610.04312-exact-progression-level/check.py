#!/usr/bin/env python3
"""Exact, bounded checks for Badr 2610.04312v1; not a universal proof."""
from fractions import Fraction
from itertools import product
import sympy as s

X, Y, Z = s.symbols('X Y Z')

def degree(poly):
    return s.Poly(poly, Z, domain=s.QQ).degree()

def irreducible(poly):
    return s.Poly(poly, Z, domain=s.QQ).is_irreducible

def main():
    # Section 9's worked quintic, including its genus and projection degree.
    F = X**5 + Y*Z**4 + Y**4*Z
    assert F.subs({X: 0, Y: 1, Z: -1}) == 0
    assert Fraction((5-1)*(5-2), 2) == 6
    assert F.subs({X: 0, Y: 0, Z: 1}) == 0
    assert degree(F.subs({X: 1, Y: 1})) == 4
    zero_fiber = Z**4 + Z
    assert s.expand(Z*(Z**3+1)) == zero_fiber
    assert s.gcd(zero_fiber, s.diff(zero_fiber, Z)) == 1
    # The homogeneous partial derivatives have no common projective zero.
    partials = [s.diff(F, v) for v in (X, Y, Z)]
    for v in (X, Y, Z):
        basis = s.groebner(partials + [v-1], X, Y, Z, domain=s.QQ)
        assert list(basis) == [1]
    print('FIDELITY: Section 9 quintic: rational point, genus 6, degree 4, four simple zeros; smoothness checked exactly')

    comparisons = 0
    for d in range(4, 101):
        for D in (d-1, d):
            for e in range(2, 101):
                assert Fraction(D, e) <= Fraction(d, 2) < d-1
                comparisons += 1
    print(f'PB_NUMERICAL_PRECONDITION: {comparisons} exact comparisons, 4<=d<=100, 2<=e<=100, D in {{d-1,d}}')

    # Two smooth families: center on C, and center outside C respectively.
    # In the latter even-degree family C(R) is empty, hence C(Q) is empty.
    configs = 0
    checked_fibers = 0
    for d in range(4, 9):
        families = [('inner', X**d+Y*Z**(d-1)+Y**(d-1)*Z, d-1),
                    ('outer', X**d+Y**d+Z**d, d)]
        for name, form, D in families:
            assert (form.subs({X: 0, Y: 0, Z: 1}) == 0) == (name == 'inner')
            for coordinate in (X, Y, Z):
                equations = [s.diff(form, v) for v in (X, Y, Z)] + [coordinate-1]
                assert list(s.groebner(equations, X, Y, Z, domain=s.QQ)) == [1]
            largest_first_witness = 0
            for kind in ('AP', 'GP'):
                for a, b in product((-2, -1, 1, 2), repeat=2):
                    # GP ratios +/-1 are torsion, excluded by Definition 2.1.
                    if kind == 'GP' and abs(b) == 1:
                        continue
                    witness = None
                    for i in range(1, 17):
                        x = a+i*b if kind == 'AP' else a*b**i
                        fiber = form.subs({X: x, Y: 1})
                        assert degree(fiber) == D
                        checked_fibers += 1
                        if irreducible(fiber):
                            witness = i
                            # Every level n<D is obstructed by this full-degree fiber.
                            assert all(n < degree(fiber) for n in range(1, D))
                            break
                    assert witness is not None, (d, name, kind, a, b)
                    configs += 1
                    largest_first_witness = max(largest_first_witness, witness)
            print(f'SMALL_CASE: d={d}, center={name}, threshold={D}, largest first irreducible index={largest_first_witness}')
    print(f'SAMPLED_LOWER_BOUND: {configs} progression configurations; {checked_fibers} exact rational fiber factorizations')
    print('SCOPE: finite algebraic checks only; Hilbert irreducibility and the universal threshold are not proved by this program')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
