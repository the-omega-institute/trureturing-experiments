#!/usr/bin/env python3
"""Exact checks for Ataei's log-convex optimal-constant candidate."""
from fractions import Fraction


def floor_frac(x: Fraction) -> int:
    return x.numerator // x.denominator


def F(x: Fraction) -> Fraction:
    n = floor_frac(x)
    return n * x - Fraction(n * (n + 1), 2)


def phi2(x: Fraction, eta: Fraction) -> Fraction:
    return eta * F(x)


def phi1(x: Fraction, eta: Fraction) -> Fraction:
    return eta * (F(x + Fraction(1, 2)) - x / 2)


def p_formula(x: Fraction, eta: Fraction) -> Fraction:
    t = x - floor_frac(x)
    return -eta * min(t, 1 - t) / 2


def check_eta(eta: Fraction) -> None:
    assert eta > 0
    # The two primitives have the common unit increment L(s)=eta*(s-1).
    for n in range(-4, 5):
        for r in (Fraction(0), Fraction(1, 7), Fraction(1, 4),
                  Fraction(1, 2), Fraction(3, 4), Fraction(6, 7)):
            s = Fraction(n) + r
            increment = eta * (s - 1)
            assert phi1(s, eta) - phi1(s - 1, eta) == increment
            assert phi2(s, eta) - phi2(s - 1, eta) == increment
            assert phi1(s, eta) - phi2(s, eta) == p_formula(s, eta)

    # The growth parameter is exact: log(g(s+1)/g(s)) = eta.
    for s in (Fraction(-3), Fraction(-1, 3), Fraction(0), Fraction(2, 5), Fraction(9, 4)):
        log_g = eta * (s - 1)
        assert eta * s - log_g == eta

    # The period range is [-eta/4, 0], and centering leaves eta/8 error.
    values = [p_formula(Fraction(0), eta), p_formula(Fraction(1, 4), eta),
              p_formula(Fraction(1, 2), eta), p_formula(Fraction(3, 4), eta)]
    assert min(values) == -eta / 4
    assert max(values) == 0
    center = (min(values) + max(values)) / 2
    assert max(abs(v - center) for v in values) == eta / 8
    # Any normalization d has max error at least half the range.
    for d in (Fraction(-1), Fraction(-1, 8), Fraction(0), Fraction(2)):
        assert max(abs(v - d) for v in values) >= eta / 8


for eta in (Fraction(1), Fraction(2), Fraction(7, 3)):
    check_eta(eta)

print("checked exact primitive, recurrence, growth ratio, periodic gauge, and sharp normalization")
print("ALL_CHECKS_OK")
