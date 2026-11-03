#!/usr/bin/env python3
"""Exact arithmetic checks for Angel--Ray--Spinka, Problem 6.2.

The checks use only Fraction and symbolic polynomial identities.  They verify
the d=2 recurrence on an exact local member of the known family, the induced
two-step operator, the inverse of h, and the defect-growth algebra used by the
proposed proof.
"""

from fractions import Fraction
from math import isqrt

import sympy as sp


def sqrt_fraction(x: Fraction) -> Fraction:
    if x < 0:
        raise ValueError("negative radicand")
    n = isqrt(x.numerator)
    d = isqrt(x.denominator)
    if n * n != x.numerator or d * d != x.denominator:
        raise ValueError(f"not a rational square: {x}")
    return Fraction(n, d)


def h(x: Fraction) -> Fraction:
    return 2 * sqrt_fraction(x) - x


def phi2(xm: Fraction, xp: Fraction) -> Fraction:
    return ((xm + xp) / 2) ** 2


def phi2_prime(xm: Fraction, xp: Fraction) -> Fraction:
    return 1 - (1 - (xm + xp) / 2) ** 2


def g2(x: Fraction, y: Fraction, z: Fraction) -> Fraction:
    return (1 - ((1 - (x + y) / 2) ** 2) / 2 - ((1 - (y + z) / 2) ** 2) / 2) ** 2


def check_known_family_local_identity() -> None:
    # Consecutive values u_i on one parity class.  All values are exact.
    u_prev = Fraction(4, 25)
    u = Fraction(16, 25)
    u_next = Fraction(24, 25)
    assert h(u_prev) == u
    assert h(u) == u_next
    assert phi2(u, u_next) == u
    assert g2(u_prev, u, u_next) == u

    # The complementary one-step identity used to pass from Phi to Phi'.
    assert phi2_prime(u_prev, u) + phi2(u_prev, u) == u_prev + u


def check_symbolic_h_identities() -> None:
    r = sp.symbols("r", nonnegative=True)
    t = r**2
    forward = sp.expand(((t + (2 * r - t)) / 2) ** 2 - t)
    assert forward == 0

    q = sp.symbols("q", nonnegative=True)
    s = (1 - q) ** 2
    inverse = sp.expand(2 * (1 - q) - s - (1 - q**2))
    assert inverse == 0

    # Exact monotonicity check on a small rational-square grid.
    grid = [Fraction(n * n, 100) for n in range(0, 11)]
    values = [h(x) for x in grid]
    assert values == sorted(values)
    assert all(Fraction(0) <= x <= 1 for x in values)


def check_defect_growth_under_perturbation() -> None:
    # Start from sqrt(u_i)=4/5 and perturb the midpoint s_i from 4/5 to 3/5.
    # The induced local u_{i+1}=2s_i-u_i is still exact rational data.
    u_i = Fraction(16, 25)
    s_i = Fraction(3, 5)
    a_i = sqrt_fraction(u_i) + s_i
    w_i = u_i - s_i * s_i
    w_prev = (2 - a_i) * w_i / a_i
    assert a_i == Fraction(7, 5)
    assert w_i == Fraction(7, 25)
    assert w_prev == Fraction(3, 25)
    assert (2 - a_i) * w_i == a_i * w_prev
    assert abs(w_i) >= 2 * abs(w_prev)


def main() -> None:
    check_known_family_local_identity()
    check_symbolic_h_identities()
    check_defect_growth_under_perturbation()
    print("known_family_local_identity: OK")
    print("symbolic_h_identities: OK")
    print("defect_growth_perturbation: OK (a_i=7/5, |w_i|/|w_prev|=7/3)")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
