#!/usr/bin/env python3
"""Exact finite checks of the OEIS A198631 Schulte weighted type-B transform.

This program checks finite cases; it does not prove the quantified conjecture.
Only integer and Fraction arithmetic supplies evidence.
"""
from fractions import Fraction
from math import comb, factorial, prod

LIMIT = 80


def type_b_rows(limit):
    rows = [[1]]
    for n in range(limit):
        previous = rows[-1]
        rows.append([
            (previous[k - 1] if k > 0 else 0)
            + (2 * k + 1) * (previous[k] if k <= n else 0)
            for k in range(n + 2)
        ])
    return rows


def logistic_egf(limit):
    """r_n = n![t^n] 2/(1+exp(-t)), by exact series inversion."""
    r = [Fraction(1)]
    for n in range(1, limit + 1):
        r.append(-sum(
            (Fraction((-1) ** j * comb(n, j)) * r[n - j]
             for j in range(1, n + 1)), Fraction(0)
        ) / 2)
    return r


def two_part(n):
    assert n >= 1
    return n & -n


def multiply(a, b, limit):
    return [sum((a[j] * b[n - j] for j in range(n + 1)), Fraction(0))
            for n in range(limit + 1)]


def main():
    rows = type_b_rows(LIMIT)
    r = logistic_egf(LIMIT)
    a = [(value / 2).numerator for value in r]
    # Fidelity: independently reproduce literal OEIS A039755 rows and A198631 terms.
    assert rows[:5] == [[1], [1, 1], [1, 4, 1], [1, 13, 9, 1], [1, 40, 58, 16, 1]]
    assert a[:13] == [1, 1, 0, -1, 0, 1, 0, -17, 0, 31, 0, -691, 0]
    assert r[:13] == list(map(Fraction, [1, '1/2', 0, '-1/4', 0, '1/2',
                                       0, '-17/8', 0, '31/2', 0, '-691/4', 0]))
    print('FIDELITY A039755 rows n=0..4:', rows[:5])
    print('FIDELITY A198631 numerators n=0..12:', a[:13])
    print('FIDELITY source r(n), n=0..12:', ', '.join(map(str, r[:13])))

    odd_factorials = [prod(range(1, 2 * k, 2)) for k in range(LIMIT + 1)]
    for n in range(LIMIT + 1):
        transformed = sum((Fraction((-1) ** k * odd_factorials[k] * rows[n][k], k + 1)
                           for k in range(n + 1)), Fraction(0))
        assert transformed == r[n], (n, transformed, r[n])
        assert Fraction(a[n], two_part(n + 1)) == r[n], (n, a[n], r[n])
        assert r[n].denominator == two_part(n + 1), (n, r[n])
        # A second, independent logistic recurrence from (1+exp(t))F=2exp(t).
        if n > 0:
            assert r[n] == 1 - sum((comb(n, j) * r[j] for j in range(n)), Fraction(0)) / 2
    print(f'TRANSFORM and literal numerator normalization: all n=0..{LIMIT}')
    print(f'EXACT reduced denominators and second logistic recurrence: all n=0..{LIMIT}')

    # Fidelity to the known type-B column EGF, using truncated ordinary series products.
    degree = 16
    exponential = [Fraction(1, factorial(n)) for n in range(degree + 1)]
    increment = [Fraction(0)] + [Fraction(2 ** n, factorial(n)) for n in range(1, degree + 1)]
    power = [Fraction(1)] + [Fraction(0)] * degree
    for k in range(degree + 1):
        column = [value / (2 ** k * factorial(k)) for value in multiply(exponential, power, degree)]
        for n in range(degree + 1):
            assert column[n] * factorial(n) == (rows[n][k] if k <= n else 0), (n, k)
        power = multiply(power, increment, degree)
    print(f'COLUMN EGF: all 0 <= n,k <= {degree}; exact coefficient products')
    print('SCOPE: finite corroboration only; no all-index settlement is certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
