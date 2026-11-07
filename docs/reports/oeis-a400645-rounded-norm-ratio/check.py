#!/usr/bin/env python3
"""Exact finite checks for the OEIS A400645 ratio-limit preregistration.

This program supplies finite evidence and algebraic identities, not a proof
of the infinite limiting assertion. No floating point is used.
"""

from fractions import Fraction as F
from math import isqrt

import sympy as sp


SOURCE_TERMS = [
    0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, 21, 28, 37,
    48, 63, 83, 109, 143, 188, 247, 324, 425, 558, 732, 961, 1261,
    1655, 2172, 2850, 3740, 4908, 6441, 8452, 11092, 14556, 19101,
    25066, 32894, 43166, 56646, 74336, 97549, 128012, 167988,
    220447, 289289,
]
N = 4096


def source_sequence(length):
    """Positive indices, with an unused dummy at zero."""
    a = [None, 0, 0, 0, 1]
    while len(a) <= length:
        n = len(a) - 1
        rad = a[n] ** 2 + a[n - 1] ** 2 + a[n - 2] ** 2 - a[n - 3] ** 2
        assert rad >= 0
        lower = isqrt(rad)
        a.append(lower + int(rad > lower * lower + lower))
    return a


def checked(condition, message):
    if not condition:
        raise AssertionError(message)


def sqrt_bracket(value, scale=10 ** 40):
    """Rational square-root enclosure certified only by integer squares."""
    value = F(value)
    checked(value >= 0, 'nonnegative rational square-root input')
    lower = F(isqrt(value.numerator * scale * scale // value.denominator), scale)
    upper = lower + F(1, scale)
    checked(lower ** 2 <= value < upper ** 2, 'rational square-root enclosure')
    return lower, upper


def main():
    a = source_sequence(N)
    checked(a[1:52] == SOURCE_TERMS, 'OEIS revision 21 data mismatch')
    print('fidelity: all 51 terms of OEIS A400645 revision 21 reproduced')
    print('positive-index prefix a_1..a_26:', a[1:27])

    for n in range(4, N):
        rad = a[n] ** 2 + a[n - 1] ** 2 + a[n - 2] ** 2 - a[n - 3] ** 2
        q = a[n + 1]
        checked(rad >= a[n] ** 2, ('signed radicand lower bound', n))
        checked(rad <= (a[n] + a[n - 1]) ** 2, ('radicand upper bound', n))
        checked(1 <= a[n] <= q <= a[n] + a[n - 1], ('growth', n))
        # These strict squared midpoint brackets independently certify rounding.
        checked((2 * q - 1) ** 2 < 4 * rad < (2 * q + 1) ** 2,
                ('independent nearest-integer bracket', n))
        error = q * q - rad
        checked(abs(error) <= 2 * a[n], ('signed squared error bound', n))
        b = lambda k: a[k] ** 2
        checked(b(n + 1) == b(n) + b(n - 1) + b(n - 2) - b(n - 3) + error,
                ('forced recurrence', n))
    print(f'finite recurrence, positivity, growth, midpoint and error checks: n=4..{N - 1}')

    growth_checks = 0
    for m in range(4, N + 1):
        for k in range(min(64, N - m) + 1):
            checked(a[m + k] * 3 ** k <= 2 * a[m] * 5 ** k,
                    ('relative Fibonacci envelope', m, k))
            growth_checks += 1
    checked(F(5, 3) ** 2 - F(5, 3) - 1 == F(1, 9), 'envelope constant')
    print(f'exact relative envelope checks: {growth_checks}; 4<=m<=4096, 0<=k<=min(64,4096-m)')

    y, x, z = sp.symbols('y x z')
    p = sp.Poly(y ** 4 - y ** 3 - y ** 2 - y + 1, y, domain=sp.QQ)
    octic = sp.Poly(x ** 8 - x ** 6 - x ** 4 - x ** 2 + 1, x, domain=sp.QQ)
    r = (1 + sp.sqrt(13) + sp.sqrt(2 * sp.sqrt(13) - 2)) / 4
    checked(sp.simplify(p.as_expr().subs(y, r)) == 0, 'radical root')
    checked(sp.expand(y ** 2 * ((y + 1 / y) ** 2 - (y + 1 / y) - 3)) == p.as_expr(),
            'reciprocal reduction')
    checked(p.gcd(p.diff()).degree() == 0, 'distinct roots')
    checked(p.count_roots(-sp.oo, sp.oo) == 2, 'two real roots')
    checked(p.count_roots(sp.Rational(17, 10), sp.Rational(7, 4)) == 1,
            'dominant real root interval')
    sqrt13_lo, sqrt13_hi = sqrt_bracket(F(13))
    inner_lo = sqrt_bracket(2 * sqrt13_lo - 2)[0]
    inner_hi = sqrt_bracket(2 * sqrt13_hi - 2)[1]
    r_lo = (1 + sqrt13_lo + inner_lo) / 4
    r_hi = (1 + sqrt13_hi + inner_hi) / 4
    checked(F(17, 10) < r_lo < r_hi < F(7, 4), 'rational radical-root enclosure')
    checked(-2 < (1 - sqrt13_hi) / 2 < (1 - sqrt13_lo) / 2 < 2,
            'unit-circle pair reduced root')
    checked(2 < (1 + sqrt13_lo) / 2, 'positive reciprocal real pair')
    checked(sp.expand(p.as_expr().subs(y, x ** 2)) == octic.as_expr(), 'octic substitution')
    print('exact algebra: radical, reciprocal reduction, distinct roots and octic substitution OK')
    print('root location: 17/10 < r < 7/4; other roots are 1/r and a unit-circle conjugate pair')
    print('p(17/10) =', p.eval(sp.Rational(17, 10)))

    companion = sp.Matrix([[1, 1, 1, -1], [1, 0, 0, 0],
                           [0, 1, 0, 0], [0, 0, 1, 0]])
    checked(companion.charpoly(y).as_expr() == p.as_expr(), 'companion characteristic')
    row = sp.Matrix([[1, y - 1, y ** 2 - y - 1, -1 / y]])
    for entry in row * companion - y * row:
        numerator = sp.cancel(entry).as_numer_denom()[0]
        checked(sp.rem(numerator, p.as_expr(), y) == 0, 'projection identity')
    print('left projection identity: all four numerator remainders modulo p are zero')

    seed = a[26]
    r_lower = F(17, 10)
    growth = F(5, 3)
    checked(r_lower - growth == F(1, 30), 'spectral gap lower threshold')
    checked(1 - 1 / r_lower == F(7, 17), 'projection lower coefficient')
    checked(r_lower ** 2 - r_lower - 1 == F(19, 100), 'second positive coefficient')
    certificate = (1 - 1 / r_lower) * seed ** 2 - 4 * seed / (r_lower - growth)
    checked(seed == 324 and certificate == F(73872, 17) and certificate > 0,
            'positive dominant-mode certificate')
    for j in range(N - 26):
        n = 26 + j
        error = a[n + 1] ** 2 - a[n] ** 2 - a[n - 1] ** 2 - a[n - 2] ** 2 + a[n - 3] ** 2
        checked(abs(error) * 3 ** j <= 4 * seed * 5 ** j, ('tail error envelope', j))
    print('positive-tail certificate: (7/17)*324^2 - 120*324 =', certificate)
    print('finite tail-envelope checks: j=0..4069')

    # Exact rational root isolation; largest real root is the rightmost interval.
    intervals = octic.intervals(eps=sp.Rational(1, 10 ** 40))
    ((lo, hi), multiplicity) = intervals[-1]
    checked(multiplicity == 1 and lo > 1, 'largest octic root isolation')
    lo, hi = F(int(lo.p), int(lo.q)), F(int(hi.p), int(hi.q))
    displayed = F(13122819079, 10 ** 10)
    half_unit = F(1, 2 * 10 ** 10)
    checked(displayed - half_unit < lo < hi < displayed + half_unit,
            'source approximate decimal, to nearest 10 decimal places')
    print('largest octic root exact isolating interval:', lo, hi)
    print('source decimal check: nearest 10-place rational is', displayed)
    tolerance = F(1, 10 ** 10)
    for n in range(128, N):
        ratio = F(a[n + 1], a[n])
        checked(hi - tolerance < ratio < lo + tolerance, ('finite ratio enclosure', n))
    print('finite ratio enclosure: |a_(n+1)/a_n - sqrt(r)| < 1/10^10 for n=128..4095')
    print('scope: finite checks and exact algebra only; no certification of the infinite limit')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
