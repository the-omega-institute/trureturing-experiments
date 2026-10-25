#!/usr/bin/env python3
"""Exact finite checks for Song's OEIS A232624 coefficient-gap conjecture.

Integer coefficient lists are in increasing exponent order. Cyclotomic division
and reciprocal extraction are checked against an independent Chebyshev divisor
construction. No real or complex floating point arithmetic is used.
"""

from functools import lru_cache
from math import comb, prod


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b, scale=1):
    c = list(a) + [0] * max(0, len(b) - len(a))
    for i, x in enumerate(b):
        c[i] += scale * x
    return trim(c)


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def exact_div(a, b):
    r = trim(a)
    q = [0] * max(1, len(r) - len(b) + 1)
    while r != [0] and len(r) >= len(b):
        k = len(r) - len(b)
        assert r[-1] % b[-1] == 0
        c = r[-1] // b[-1]
        q[k] = c
        for j, x in enumerate(b):
            r[k + j] -= c * x
        r = trim(r)
    assert r == [0], (a, b, r)
    assert mul(q, b) == trim(a)
    return trim(q)


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def distinct_primes(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result.append(n)
    return result


def totient(n):
    r = n
    for p in distinct_primes(n):
        r = r // p * (p - 1)
    return r


@lru_cache(None)
def cyclotomic(n):
    f = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n)[:-1]:
        f = exact_div(f, cyclotomic(d))
    return tuple(f)


def reciprocal_lift(p, d):
    """Return z^d*p(z+z^(-1)) using integer binomial coefficients."""
    assert len(p) - 1 <= d
    r = [0] * (2 * d + 1)
    for j, c in enumerate(p):
        for k in range(j + 1):
            r[d + j - 2 * k] += c * comb(j, k)
    return trim(r)


@lru_cache(None)
def real_poly(n):
    if n == 1:
        return (-2, 1)
    if n == 2:
        return (2, 1)
    phi = list(cyclotomic(n))
    assert phi == phi[::-1] and (len(phi) - 1) % 2 == 0
    d = (len(phi) - 1) // 2
    r, p = phi[:], [0] * (d + 1)
    for j in range(d, -1, -1):
        c = r[d + j]
        p[j] = c
        for k in range(j + 1):
            r[d + j - 2 * k] -= c * comb(j, k)
    assert not any(r)
    assert p[-1] == 1 and d == totient(n) // 2
    assert reciprocal_lift(p, d) == phi
    return tuple(p)


@lru_cache(None)
def chebyshev(n):
    """t_n(x)=2*T_n(x/2), via t_0=2, t_1=x, t_n=x*t_(n-1)-t_(n-2)."""
    if n == 0:
        return (2,)
    if n == 1:
        return (0, 1)
    return tuple(add([0] + list(chebyshev(n - 1)), chebyshev(n - 2), -1))


@lru_cache(None)
def real_poly_by_chebyshev(n):
    if n == 1:
        return (-2, 1)
    if n % 2:
        q = add(chebyshev((n + 1) // 2), chebyshev((n - 1) // 2), -1)
    else:
        q = add(chebyshev(n // 2 + 1), chebyshev(n // 2 - 1), -1)
    for d in divisors(n)[:-1]:
        q = exact_div(q, real_poly_by_chebyshev(d))
    return tuple(q)


def coeff(p, j):
    return p[j] if j < len(p) else 0


def main():
    # Source EXAMPLE table, rows 1..20, not generated from the conjectured gap.
    source_rows = [
        [-2, 1], [2, 1], [1, 1], [0, 1], [-1, 1, 1], [-1, 1],
        [-1, -2, 1, 1], [-2, 0, 1], [1, -3, 0, 1], [-1, -1, 1],
        [1, 3, -3, -4, 1, 1], [-3, 0, 1], [-1, 3, 6, -4, -5, 1, 1],
        [1, -2, -1, 1], [1, 4, -4, -1, 1], [2, 0, -4, 0, 1],
        [1, -4, -10, 10, 15, -6, -7, 1, 1], [-1, -3, 0, 1],
        [1, 5, -10, -20, 15, 21, -7, -8, 1, 1], [5, 0, -5, 0, 1],
    ]
    for n, row in enumerate(source_rows, 1):
        assert list(real_poly(n)) == row
    print('FIDELITY: all 20 source coefficient-table rows reproduced exactly')
    examples = {3, 9, 15, 105, 225, 289, 441}
    count = 0
    for n in range(3, 452, 2):
        p, q = real_poly(n), real_poly(2 * n)
        assert p == real_poly_by_chebyshev(n)
        assert q == real_poly_by_chebyshev(2 * n)
        primes = distinct_primes(n)
        m, d, eps = n // prod(primes), totient(n) // 2, (-1) ** len(primes)
        delta = add(q, p, -1)
        assert all(c % 2 == 0 for c in delta)
        f = trim([c // 2 for c in delta])
        assert f != [0] and len(f) - 1 == d - m and f[-1] == eps
        assert d >= m > 0 and m % 2 == 1
        for k in range(1, m, 2):
            assert coeff(p, d - k) == coeff(q, d - k) == 0
        assert coeff(p, d - m) == -eps and coeff(q, d - m) == eps
        assert all(coeff(q, j) == (-1) ** (d - j) * coeff(p, j)
                   for j in range(d + 1))
        if n == 15:
            assert f == [0, -4, 0, 1]
            print('FIDELITY: n=15 half-difference is x^3 - 4*x')
        if n in examples:
            print(f'CASE n={n} d={d} n/rad(n)={m} degree(F)={len(f)-1} LC(F)={f[-1]}')
        count += 1
    assert count == 225
    print(f'EXACT_CHECK: {count} odd conductors 3..451; nonzero, degree, leading coefficient, all quoted coefficient consequences')
    print('CROSS_CHECK: both polynomials in every case equal the independent Chebyshev divisor construction')
    print('BOUNDARY: n=3 gives F=-1; prime powers and repeated-prime composite conductors included')
    print('SCOPE: finite exact checks only; no unbounded theorem or Lean proof established')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
