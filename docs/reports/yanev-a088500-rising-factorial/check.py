#!/usr/bin/env python3
"""Exact finite checks for the integral conjecture in OEIS A088500.

No numerical quadrature or floating point is used. The polynomial Q below
satisfies (Q-Q')=P, so -exp(-x)Q(x) is an exact antiderivative of exp(-x)P(x).
For each fixed polynomial the boundary at infinity is zero; the program checks
the algebra and its value Q(0), not that analytic boundary limit or all n.
"""
from fractions import Fraction
from math import comb, factorial

# OEIS revision 24, Jan 12 2026, %S/%T/%U (first 20 terms, offset 0).
SOURCE_TERMS = [
    1, 2, 10, 76, 772, 9808, 149552, 2660544, 54093696, 1237306560,
    31446049728, 879119219328, 26811313164672, 885830291432448,
    31518653868782592, 1201567079771092992, 48860409899753588736,
    2111033523652100407296, 96573225733453262979072,
    4663365352616114000510976,
]
LIMIT = 50


def egf_coefficients(limit):
    """Invert 1+2 log(1-z) using its exact ordinary coefficients."""
    denominator = [Fraction(1)] + [Fraction(-2, n) for n in range(1, limit + 1)]
    inv = [Fraction(1)]
    for n in range(1, limit + 1):
        inv.append(-sum(denominator[j] * inv[n-j] for j in range(1, n+1)))
    for n in range(limit + 1):
        assert sum(denominator[j] * inv[n-j] for j in range(n+1)) == (n == 0)
    values = [inv[n] * factorial(n) for n in range(limit + 1)]
    assert all(v.denominator == 1 for v in values)
    return [v.numerator for v in values]


def stirling_values(limit):
    """Unsigned first-kind recurrence, independent of series inversion."""
    row = [1]
    values = []
    for n in range(limit + 1):
        values.append(sum(c * factorial(k) * 2**k for k, c in enumerate(row)))
        new = [0] * (len(row) + 1)
        for k, c in enumerate(row):
            new[k] += n * c
            new[k+1] += c
        row = new
    return values


def rising_polynomial(n):
    """Coefficients of product_{j=0}^{n-1}(2*x+j), constant first."""
    p = [1]
    for j in range(n):
        new = [0] * (len(p) + 1)
        for k, c in enumerate(p):
            new[k] += j * c
            new[k+1] += 2 * c
        p = new
    return p


def primitive_at_zero(p):
    """Solve Q-Q'=P over integers and verify every coefficient."""
    q = [0] * len(p)
    for k in range(len(p)-1, -1, -1):
        q[k] = p[k] + ((k+1) * q[k+1] if k+1 < len(q) else 0)
    assert all(q[k] - ((k+1) * q[k+1] if k+1 < len(q) else 0) == p[k]
               for k in range(len(p)))
    return q[0]


def main():
    egf = egf_coefficients(LIMIT)
    stirling = stirling_values(LIMIT)
    primitive = [primitive_at_zero(rising_polynomial(n)) for n in range(LIMIT+1)]
    assert egf[:len(SOURCE_TERMS)] == SOURCE_TERMS
    assert egf == stirling == primitive
    # Also reproduce Gutkovskiy's recurrence in the same OEIS source.
    assert egf[0] == 1
    for n in range(1, LIMIT+1):
        assert egf[n] == 2 * sum(comb(n,k) * factorial(n-k-1) * egf[k]
                                 for k in range(n))
    assert rising_polynomial(0) == [1]
    assert primitive[0] == 1
    print('Source: OEIS A088500, revision 24, Jan 12 2026; offset 0')
    print('Fidelity: all 20 published initial terms agree exactly')
    print('First ten:', ','.join(map(str, egf[:10])))
    print(f'Exact EGF inversion, weighted Stirling sum, and primitive constants agree for n=0..{LIMIT}')
    print(f'Source recurrence agrees for n=1..{LIMIT}; n=0 empty product checked')
    print('Integral check uses exact polynomial antiderivatives; no numerical quadrature')
    print('These finite algebra checks do not prove the conjecture for all n or integrability')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
