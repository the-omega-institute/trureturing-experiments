#!/usr/bin/env python3
"""Exact checks for Zhai's Conjecture 1 (arXiv:2607.17474v1).

Only integers and Fraction are used as numerical evidence. Infinite C_h is
retained symbolically; its finite prefix and infinite remainder have rigorous
rational enclosures. No truncated series is reported as the full C_h.
"""

from collections import Counter
from fractions import Fraction as Q
from math import isqrt


def factors(n):
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def mu(n):
    ff = factors(n)
    return 0 if any(e > 1 for _, e in ff) else (-1) ** len(ff)


def tau(n):
    out = 1
    for _, e in factors(n):
        out *= e + 1
    return out


def a(x, m):
    assert x >= 1 and m >= 1
    return x // m - x // (m + 1)


def h(m):
    assert m >= 1
    t = 1024
    while t <= m:
        if m < 2 * t:
            return int(t * t // m == t * t // (m + 1))
        t *= t
    return 0


def source_fidelity():
    # The two exact identities immediately before Corollary 2, Section 1.2.
    for n in range(1, 1025):
        rhs = sum(tau(n // (d * d)) * mu(d)
                  for d in range(1, isqrt(n) + 1) if n % (d * d) == 0)
        assert 2 ** len(factors(n)) == rhs
        liouville = (-1) ** sum(e for _, e in factors(n))
        assert liouville == sum(mu(n // (d * d))
                                for d in range(1, isqrt(n) + 1)
                                if n % (d * d) == 0)
    # Equation (1.1), including noninteger x: S_1(x) = floor(x).
    for q in range(1, 8):
        for p in range(q, 101 * q):
            x = Q(p, q)
            nmax = x.numerator // x.denominator
            assert sum(1 for n in range(1, nmax + 1)
                       if x // n >= 1) == nmax
    print('FIDELITY: both source identities before Corollary 2 hold for 1<=n<=1024;'
          ' equation (1.1) for f=1 checked on rational x with denominators 1..7.')


def small_cases():
    for x in range(1, 129):
        counts = Counter(x // n for n in range(1, x + 1))
        for m in range(1, x + 2):
            assert counts[m] == a(x, m)
            e = a(x, m) - Q(x, m * (m + 1))
            fracdiff = Q(x % (m + 1), m + 1) - Q(x % m, m)
            assert e == fracdiff and abs(e) < 1
    for t in range(2, 257, 2):
        coefficients = [a(t * t, m) for m in range(t, 2 * t)]
        assert set(coefficients) <= {0, 1}
        assert coefficients.count(0) == t // 2
        assert sum(coefficients) == t // 2
    for n in range(1, 4097):
        assert 1 <= 1 + h(n) <= 2
    print('SMALL_CASES: quotient counts and fractional-part identity for X<=128;'
          ' hole counts for all even 2<=T<=256; fixed f takes values 1 or 2.')


def checked_blocks():
    # Each rational term lies in [floor(D/den)/D, (floor(D/den)+1)/D].
    # This is exact outward rounding by integer division, not floating point.
    den_scale = 1 << 100
    blocks = []
    t = 1024
    for k in range(2):
        x = t * t
        holes = []
        lattice_sum = 0
        count_coefficients = 0
        for m in range(t, 2 * t):
            coefficient = a(x, m)
            assert coefficient in (0, 1)
            count_coefficients += coefficient
            if coefficient == 0:
                holes.append(m)
                assert h(m) == 1
                assert 4 * x > m * (m + 1)
                lattice_sum += den_scale // (m * (m + 1))
        assert count_coefficients == t // 2
        assert len(holes) == t // 2
        lo = Q(lattice_sum, den_scale)
        hi = Q(lattice_sum + len(holes), den_scale)
        assert x * lo > Q(t, 8)
        # Endpoint of the next block, which begins at x, is not a hole.
        assert a(x * x, x) == 1 and h(x) == 0
        # At observation X=x, every subsequent supported m is strictly >x.
        assert a(x, x + 1) == 0
        blocks.append((t, holes, lo, hi))
        print(f'BLOCK {k}: T={t}; X={x}; holes={len(holes)};'
              f' current weighted mass > T/8={Q(t, 8)}; next endpoint h(X)=0')
        t = x
    return blocks, t


def infinite_error_enclosures(blocks, first_omitted_scale):
    # All unenumerated supported m are >=first_omitted_scale. Positivity and
    # sum_{m>=M} 1/(m(m+1))=1/M enclose the entire infinite remainder.
    c_lo = sum((b[2] for b in blocks), Q(0))
    c_hi = sum((b[3] for b in blocks), Q(0)) + Q(1, first_omitted_scale)
    assert c_lo > 0 and c_hi < 1
    all_enumerated_holes = [m for _, hh, _, _ in blocks for m in hh]
    print(f'C_h enclosure: [{c_lo}, {c_hi}];'
          f' infinite tail beyond the two blocks <= 1/{first_omitted_scale}')
    for k, (t, holes, lo, hi) in enumerate(blocks):
        x = t * t
        sh = sum(a(x, m) for m in all_enumerated_holes)
        # Later support is >x, so this integer is the exact full S_h(X).
        elo, ehi = Q(sh) - x * c_hi, Q(sh) - x * c_lo
        assert ehi < -Q(t, 16)
        # Cross-check earlier-block interference using exact Fractions.
        old = [m for tt, hh, _, _ in blocks if tt < t for m in hh]
        old_error = sum((a(x, m) - Q(x, m * (m + 1)) for m in old), Q(0))
        assert abs(old_error) <= len(old)
        if k:
            assert len(old) < 2 * isqrt(t)
            assert 2 * isqrt(t) <= Q(t, 16)
        print(f'ERROR {k}: exact S_h(X)={sh}; exact E_f(X)={sh}-{x}*C_h;'
              f' rational enclosure=[{elo}, {ehi}]; E_f(X)<-sqrt(X)/16={-Q(t, 16)}')
    # Exact comparison of exponents, rather than a decimal approximation.
    assert Q(2, 5) + Q(1, 20) == Q(9, 20)
    assert Q(1, 2) - Q(9, 20) == Q(1, 20)
    print('EXPONENTS: epsilon=1/20; 2/5+epsilon=9/20;'
          ' (1/2)-(9/20)=1/20>0. Infinite-scale implication is the prose argument.')


def main():
    source_fidelity()
    small_cases()
    blocks, first_omitted_scale = checked_blocks()
    infinite_error_enclosures(blocks, first_omitted_scale)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
