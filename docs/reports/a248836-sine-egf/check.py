#!/usr/bin/env python3
"""Exact checks for OEIS A248836 and Manyama's EGF conjecture."""

from fractions import Fraction
from itertools import product
from math import factorial


SOURCE = [1, 2, 6, 22, 96, 482, 2736, 17302, 120576, 917762]


def count_arrays(n: int) -> int:
    """Count p with 0 <= p(i) <= i and every value used at most twice."""
    total = 0
    for values in product(*(range(i + 1) for i in range(1, n + 1))):
        uses = {}
        for value in values:
            uses[value] = uses.get(value, 0) + 1
        if max(uses.values(), default=0) <= 2:
            total += 1
    return total


def sine_coefficients(n: int):
    coeff = [Fraction(0) for _ in range(n + 1)]
    for k in range((n - 1) // 2 + 1):
        degree = 2 * k + 1
        if degree <= n:
            coeff[degree] = Fraction((-1) ** k, factorial(degree))
    return coeff


def inverse_series(series):
    """Return coefficients of the inverse of a series with constant term 1."""
    n = len(series) - 1
    inv = [Fraction(0) for _ in range(n + 1)]
    inv[0] = 1
    for d in range(1, n + 1):
        inv[d] = -sum(series[j] * inv[d - j] for j in range(1, d + 1))
    return inv


def multiply(left, right):
    n = min(len(left), len(right)) - 1
    return [sum(left[j] * right[d - j] for j in range(d + 1)) for d in range(n + 1)]


def main() -> None:
    enumerated = [count_arrays(n) for n in range(len(SOURCE))]
    assert enumerated == SOURCE, (enumerated, SOURCE)
    print("enumeration n=0..9:", enumerated)

    n = len(SOURCE) - 1
    sine = sine_coefficients(n)
    base = [Fraction(int(d == 0)) - sine[d] for d in range(n + 1)]
    inverse = inverse_series(base)
    rhs = multiply(inverse, inverse)
    egf = [Fraction(a, factorial(d)) for d, a in enumerate(enumerated)]
    assert rhs == egf, (rhs, egf)
    print("EGF coefficients n=0..9:", [str(x) for x in rhs])
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
