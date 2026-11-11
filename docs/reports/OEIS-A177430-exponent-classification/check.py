#!/usr/bin/env python3
"""Exact checks for OEIS A177430/A177431 (no floating point)."""

from fractions import Fraction
from math import isqrt


def exp_coefficients(exponents, degree):
    """Coefficients of exp(sum exponents[n] x^n/n), as Fractions."""
    coeff = [Fraction(0) for _ in range(degree + 1)]
    coeff[0] = Fraction(1)
    for n in range(1, degree + 1):
        coeff[n] = sum(exponents[j] * coeff[n - j] for j in range(1, n + 1)) / n
    return coeff


def greedy_squares(degree):
    """Generate Hanna's least increasing square logarithmic derivative."""
    a = [0]
    coeff = [Fraction(1)]
    for n in range(1, degree + 1):
        previous = a[-1]
        root = isqrt(previous) + 1
        while True:
            candidate = root * root
            numerator = sum(a[j] * coeff[n - j] for j in range(1, n)) + candidate
            value = Fraction(numerator, n)
            if value.denominator == 1:
                a.append(candidate)
                coeff.append(value)
                break
            root += 1
    return a[1:], coeff


N = 50
a, base_coeff = greedy_squares(N)
assert a[:7] == [1, 9, 16, 25, 36, 144, 169]
assert base_coeff[:8] == [Fraction(x) for x in [1, 1, 5, 10, 24, 51, 122, 244]]
print("A177430 first 7 terms:", a[:7])
print("exp(A(x)) through degree 7:", [int(x) for x in base_coeff[:8]])

roots = [isqrt(x) for x in a]
assert roots[:5] == [1, 3, 4, 5, 6]
print("A177431 first 5 terms:", roots[:5])

integral_by_m = {}
first_failure_by_m = {}
for m in range(-20, 21):
    if m >= 0:
        powers = [Fraction(s**m) for s in roots]
    else:
        powers = [Fraction(1, s**(-m)) for s in roots]
    coeff = exp_coefficients([Fraction(0)] + powers, N)
    integral_by_m[m] = all(c.denominator == 1 for c in coeff)
    first_failure_by_m[m] = next((n for n, c in enumerate(coeff) if c.denominator != 1), None)

expected = {m for m in range(0, 21, 2)}
actual = {m for m, ok in integral_by_m.items() if ok}
assert actual == expected, (actual, expected)
assert all(first_failure_by_m[m] is not None for m in range(-20, 21) if m not in expected)
print("integral through degree 50 exactly for m:", sorted(actual))
print("first failures for m=-3,-1,1,3,5:", {m: first_failure_by_m[m] for m in [-3, -1, 1, 3, 5]})

# The four-coefficient obstruction is visible directly: g_4-g_2 = 5^m-3^m.
for m in range(-20, 21):
    if m >= 0:
        powers = [Fraction(s**m) for s in roots[:4]]
    else:
        powers = [Fraction(1, s**(-m)) for s in roots[:4]]
    first4 = exp_coefficients([Fraction(0)] + powers, 4)
    first4_integral = all(c.denominator == 1 for c in first4)
    assert first4_integral == (m in expected)
print("first four coefficients characterize the tested exponent range")
print("ALL_CHECKS_OK")
