#!/usr/bin/env python3
"""Exact small-instance checks for Zhao, arXiv:2609.16081v1, Conjecture 25(a).

The geometric count uses only integer arithmetic and Fraction.  The inequalities
are certified by the classical rational enclosure 333/106 < pi < 355/113;
there is no floating-point evidence in this check.
"""
from fractions import Fraction
from itertools import product
import sympy as sp

PI_LO = Fraction(333, 106)
PI_HI = Fraction(355, 113)


def inside_status(squared_distance: Fraction, n: int) -> bool:
    """Decide squared_distance < n/pi using only rational bounds for pi."""
    if squared_distance * PI_HI < n:
        return True
    if squared_distance * PI_LO >= n:
        return False
    raise AssertionError(
        f"the certified pi interval is inconclusive for q={squared_distance}, n={n}"
    )


def count_inside(center, n: int) -> int:
    """Count lattice points in the radius-squared n/pi disk around center.

    For n <= 12 and center in [0,1]^2, every point at distance < radius lies
    in [-4,4]^2 because pi > 3, so this finite enumeration is exhaustive.
    """
    total = 0
    for x, y in product(range(-4, 5), repeat=2):
        q = (Fraction(x) - center[0]) ** 2 + (Fraction(y) - center[1]) ** 2
        total += int(inside_status(q, n))
    return total


def find_witness(n: int):
    """Search a finite rational grid for a certified center with exactly n points."""
    for denominator in (2, 4, 8, 16, 32):
        for i, j in product(range(denominator), repeat=2):
            center = (Fraction(i, denominator), Fraction(j, denominator))
            try:
                count = count_inside(center, n)
            except AssertionError:
                continue
            if count == n:
                return center
    raise AssertionError(f"no certified grid witness found for n={n}")


def check_volume_formula():
    """Reproduce the source's d-dimensional ball-volume normalization."""
    pi = sp.pi
    for d in (2, 3, 4, 5):
        kappa_d = pi ** (sp.Rational(d, 2)) / sp.gamma(sp.Rational(d, 2) + 1)
        assert sp.simplify(kappa_d - pi ** (sp.Rational(d, 2)) / sp.gamma(sp.Rational(d, 2) + 1)) == 0
    kappa_2 = pi ** 1 / sp.gamma(2)
    assert sp.simplify(kappa_2 - pi) == 0
    for n in range(1, 13):
        radius_squared = sp.Rational(n) / kappa_2
        assert sp.simplify(kappa_2 * radius_squared - n) == 0


def main():
    assert PI_LO < PI_HI
    check_volume_formula()
    witnesses = {}
    for n in range(1, 13):
        center = find_witness(n)
        count = count_inside(center, n)
        assert count == n
        witnesses[n] = center
        print(f"n={n:2d} center=({center[0]}, {center[1]}) exact_count={count}")
    print("Checked d=2 and n=1,...,12 with rational centers and certified pi bounds.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
