#!/usr/bin/env python3
"""Exact finite checks for Browne's final FORMULA in OEIS A164660.

Source: https://oeis.org/search?q=id:A164660&fmt=text, revision 19.
Only Python integers and fractions.Fraction are used; no floating point.
The finite checks are evidence for the preregistration, not an all-n proof.
"""

from fractions import Fraction
from math import comb, prod

LIMIT = 300


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def chebyshev_recurrence(last):
    """Coefficient lists in increasing powers: T0=1, T1=X, Tn+2=2X Tn+1-Tn."""
    polys = [[1], [0, 1]]
    for n in range(last - 1):
        nxt = [0] + [2 * c for c in polys[-1]]
        for j, c in enumerate(polys[-2]):
            nxt[j] -= c
        polys.append(nxt)
    return polys[:last + 1]


def chebyshev_coefficients_closed(m):
    """Independent binomial coefficient formula, evaluated over exact rationals."""
    if m == 0:
        return [1]
    coefficients = [0] * (m + 1)
    for j in range(m // 2 + 1):
        c = Fraction((-1) ** j * m * comb(m - j, j) * 2 ** (m - 2 * j),
                     2 * (m - j))
        require(c.denominator == 1, ('nonintegral coefficient', m, j))
        coefficients[m - 2 * j] = c.numerator
    return coefficients


def integrated_row(coefficients):
    """Columns 1,...,m+1 of the zero-constant primitive in A164658/A164659."""
    return [Fraction(c, j + 1) for j, c in enumerate(coefficients)]


def radical(n):
    require(n >= 1, 'positive modulus required')
    remaining, result, p = n, 1, 2
    while p * p <= remaining:
        if remaining % p == 0:
            result *= p
            while remaining % p == 0:
                remaining //= p
        p += 1
    if remaining > 1:
        result *= remaining
    return result


def main():
    polys = chebyshev_recurrence(LIMIT + 1)
    for m, coefficients in enumerate(polys):
        require(coefficients == chebyshev_coefficients_closed(m),
                ('independent polynomial check', m))
    integrals = [sum(integrated_row(p), Fraction(0)) for p in polys]

    # Verbatim %e rational list of A164660 revision 19 (indices 0,...,10).
    published_rationals = [Fraction(s) for s in
                          '1 1/2 -1/3 -1/2 -1/15 1/6 -1/35 -1/6 -1/63 1/10 -1/99'.split()]
    require(integrals[:11] == published_rationals, 'source rational examples')
    # First five rows printed in %e of A164658 and A164659.
    published_rows = ['1', '0 1/2', '-1 0 2/3', '0 -3/2 0 1', '1 0 -8/3 0 8/5']
    for m, row in enumerate(published_rows):
        require(integrated_row(polys[m]) == [Fraction(s) for s in row.split()],
                ('source rational triangle', m))
    print('FIDELITY_OK: A164660 %e, m=0..10; A164658/A164659 %e, rows=0..4')
    print('POLYNOMIAL_CROSSCHECK_OK: recurrence versus binomial formula, m=0..301')

    indicators_checked = 0
    samples = []
    for n in range(1, LIMIT + 1):
        r = radical(n)
        require(n % r == 0, ('radical divisibility', n))
        increments, residues = [], []
        for k in range(1, n + 1):
            power = k ** n
            # Python // is exact floor division for these integer numerators.
            delta = power // n - (power - 1) // n
            divides = power % n == 0
            require(delta in (0, 1) and delta == int(divides),
                    ('floor indicator', n, k))
            require(divides == (k % r == 0), ('radical equivalence', n, k))
            increments.append(delta)
            residues.append(divides)
            indicators_checked += 1
        literal_product = prod(1 - 2 * delta for delta in increments)
        count = sum(residues)
        quotient = n // r
        require(count == quotient, ('nilpotent count', n))
        require(literal_product == (-1) ** quotient, ('product parity', n))
        require(integrals[n + 1].numerator == literal_product,
                ('Browne full identity', n))
        require((quotient % 2 == 0) == (n % 4 == 0), ('parity bridge', n))
        m = n + 1
        if m % 2 == 0:
            closed_integral = Fraction(-1, m * m - 1)
        elif m % 4 == 1:
            closed_integral = Fraction(1, m + 1)
        else:
            closed_integral = Fraction(-1, m - 1)
        require(integrals[m] == closed_integral, ('integral evaluation', m))
        if n in (1, 2, 4, 9, 12, 16, 300):
            samples.append((n, str(integrals[n + 1]), r, count, literal_product))
    print(f'IDENTITY_OK: both equalities for every n=1..{LIMIT}')
    print(f'INDICATORS_OK: {indicators_checked} pairs (n,k), 1<=k<=n<=300')
    print('INTEGRAL_AND_PARITY_BRIDGES_OK: m=2..301 and n=1..300')
    print('SAMPLES: n, I_(n+1), rad(n), count, product')
    for sample in samples:
        print(', '.join(map(str, sample)))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
