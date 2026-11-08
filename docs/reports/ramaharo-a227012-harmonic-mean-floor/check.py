#!/usr/bin/env python3
"""Exact checks for Ramaharo's FORMULA conjecture in OEIS A227012.

Source: https://oeis.org/A227012, revision #45, final FORMULA (2018).
The finite checks are evidence, not a Lean proof of the unbounded claim.
Requires Python >= 3.12 and SymPy >= 1.14; all arithmetic is exact.
"""

from fractions import Fraction as F
import sys

import sympy as s

assert sys.version_info >= (3, 12)

# All three published data lines in OEIS revision #45, offset one.
SOURCE_DATA = [
    1, 4, 16, 43, 91, 166, 275, 422, 614, 857, 1158, 1521, 1953,
    2460, 3049, 3724, 4492, 5359, 6332, 7415, 8615, 9938, 11391,
    12978, 14706, 16581, 18610, 20797, 23149, 25672, 28373, 31256,
    34328, 37595, 41064, 44739, 48627, 52734, 57067, 61630,
    66430, 71473,
]
UPPER = [
    40, 1084, 13388, 99009, 484718, 1639762, 3880584, 6347425,
    6883430, 4501560, 1394840, 50080,
]
LOWER = [
    1120, 46928, 932176, 11660832, 102899870, 679866931, 3483079839,
    14134840232, 46005088504, 120792332859, 255897925487, 435032032522,
    586417085916, 613945294592, 482419274848, 268315711104,
    94408709120, 15840993280,
]


def reciprocal_sum(lo, hi):
    """Balanced exact summation on the closed positive-integer interval."""
    if lo > hi:
        return F(0)
    if hi - lo < 32:
        return sum((F(1, k) for k in range(lo, hi + 1)), F(0))
    mid = (lo + hi) // 2
    return reciprocal_sum(lo, mid) + reciprocal_sum(mid + 1, hi)


def floor(q):
    return q.numerator // q.denominator


def proposed(n):
    numerator = 2 * n**3 - 3 * n**2 + n + 2
    assert numerator % 2 == 0
    return numerator // 2 + max(0, n - 3) // 4


def check_symbolic():
    i = s.Symbol('i', integer=True)
    t = s.Symbol('t', positive=True, integer=True)
    formulas = {
        1: 0,
        2: (t**2 - 1) / 12,
        3: 0,
        4: (t**2 - 1) * (3 * t**2 - 7) / 240,
        5: 0,
        6: (t**2 - 1) * (3 * t**4 - 18 * t**2 + 31) / 1344,
    }
    for r, expected in formulas.items():
        actual = s.summation((i - (t - 1) / 2)**r, (i, 0, t - 1)) / t
        assert s.Poly(s.cancel(actual - expected), t).is_zero, r
    print('CENTRAL_MOMENTS: exact symbolic identities for degrees 1 through 6')

    v, c0, u = s.symbols('v c0 u')
    assert s.Poly(s.cancel(1 / (1 - v) - (1 + v + v**2 + v**3 / (1 - v))), v).is_zero
    assert s.Poly(s.cancel((1 / (c0 + u) + 1 / (c0 - u)) / 2 - c0 / (c0**2 - u**2)), c0, u).is_zero
    print('RECIPROCAL_IDENTITIES: symmetric pairing and exact geometric remainder')

    n, z = s.symbols('n z')
    p, q = (n - 1)**3, n**3
    N = q - p
    c, h = (p + q + 1) / 2, (N - 1) / 2
    P = n**3 - s.Rational(3, 2) * n**2 + s.Rational(3, 4) * n + s.Rational(3, 8)
    m2, m4, m6 = [s.sympify(formulas[r]).subs(t, N) for r in (2, 4, 6)]
    A = 1 + m2 / c**2 + m4 / c**4
    B = A + m6 / (c**6 * (1 - h**2 / c**2))
    d, e = 2 * n**2 - 3 * n + 3, n**2 - 3 * n + 3
    den_upper = 20 * n**3 * d**4
    den_lower = 560 * n**7 * e * d**4
    U = s.Poly(s.cancel(((P + s.Rational(1, 8)) * A - c) * den_upper), n)
    V = s.Poly(s.cancel((c - (P - s.Rational(1, 8)) * B) * den_lower), n)
    uc = s.Poly(U.as_expr().subs(n, z + 3).expand(), z).all_coeffs()
    vc = s.Poly(V.as_expr().subs(n, z + 3).expand(), z).all_coeffs()
    assert uc == UPPER
    assert vc == LOWER
    assert all(x > 0 for x in uc + vc)
    assert s.expand(d.subs(n, z + 3)) == 2 * z**2 + 9 * z + 12
    assert s.expand(e.subs(n, z + 3)) == z**2 + 3 * z + 3
    print('UPPER_CERTIFICATE: coefficients in descending powers of z=n-3:', list(uc))
    print('LOWER_CERTIFICATE: coefficients in descending powers of z=n-3:', list(vc))
    print('UNIFORM_CERTIFICATES: both rational identities exact; all shifted coefficients positive')

    for residue, numerator in enumerate([3, 5, 7, 1]):
        poly = s.Poly(P.subs(n, 4 * z + residue) - s.Rational(numerator, 8), z)
        assert all(x.q == 1 for x in poly.all_coeffs())
        frac = s.Rational(numerator, 8)
        assert frac - s.Rational(1, 8) >= 0
        assert frac + s.Rational(1, 8) <= 1
    I = (2 * n**3 - 3 * n**2 + n + 2) / 2
    assert s.expand(P - I - (2 * n - 5) / 8) == 0
    for residue in range(4):
        assert all(x.q == 1 for x in s.Poly(I.subs(n, 4 * z + residue), z).all_coeffs())
        assert (2 * residue - 5) // 8 == (residue - 3) // 4
    print('ROUNDING: period-four fractional pattern and floor translations checked exactly')


def check_finite():
    integral_means = []
    for n in range(1, 65):
        p, q = (n - 1)**3, n**3
        N = q - p
        H = reciprocal_sum(p + 1, q)
        M = F(N) / H
        assert N > 0 and H > 0
        assert floor(M) == proposed(n), n
        if n <= len(SOURCE_DATA):
            assert floor(M) == SOURCE_DATA[n - 1], n
        if M.denominator == 1:
            integral_means.append(n)
        if n == 1:
            assert M == 1
        if n == 2:
            assert H == F(481, 280) and M == F(1960, 481)
        if n >= 3:
            c, h = F(p + q + 1, 2), F(N - 1, 2)
            P = n**3 - F(3, 2) * n**2 + F(3, 4) * n + F(3, 8)
            m2 = F(N**2 - 1, 12)
            m4 = F((N**2 - 1) * (3 * N**2 - 7), 240)
            m6 = F((N**2 - 1) * (3 * N**4 - 18 * N**2 + 31), 1344)
            A = 1 + m2 / c**2 + m4 / c**4
            B = A + m6 / (c**6 * (1 - h**2 / c**2))
            assert 0 <= h < c
            assert A <= c * H / N <= B
            assert P - F(1, 8) < c / B <= M <= c / A < P + F(1, 8)
            assert floor(M) == floor(P)
            if n <= 12:
                for r, expected in [(2, m2), (4, m4), (6, m6)]:
                    assert sum(((F(i) - h)**r for i in range(N)), F(0)) / N == expected
        print(f'n={n}: floor(M)={floor(M)}; formula=OK; exact_bounds={"OK" if n >= 3 else "initial case"}')
    assert integral_means == [1]
    print(f'FIDELITY: all {len(SOURCE_DATA)} published OEIS data terms reproduced from reciprocal sums')
    print('INITIAL_CASES: M(1)=1; H(2)=481/280; M(2)=1960/481; floors=1,4')
    print('FINITE_CLAIM: n=1..64 checked with Fractions; only integral mean in this range is n=1')


if __name__ == '__main__':
    check_symbolic()
    check_finite()
    print('ALL_CHECKS_OK')
