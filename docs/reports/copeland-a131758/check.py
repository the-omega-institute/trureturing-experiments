#!/usr/bin/env python3
"""Exact pre-probe checks for the OEIS A131758 polynomial-content conjecture.

Python 3 standard library only; no floating point, external service, or Lean proof.
Source: https://oeis.org/search?q=id:A131758&fmt=text (revision 48).
"""
from fractions import Fraction
from functools import reduce
from math import comb, factorial, gcd, isqrt


# Published COMMENTS, P(0,t) through P(7,t), in increasing powers of t.
SOURCE_ROWS = [
    [1], [0, 1], [-1, 1, 2], [4, -14, 10, 6],
    [-15, 83, -157, 89, 24], [56, -424, 1266, -1724, 826, 120],
    [-185, 1887, -8038, 17642, -19593, 8287, 720],
    [204, -4976, 36226, -126944, 239576, -234688, 90602, 5040],
]
# Published https://oeis.org/A131758/b131758.txt, indices 36 through 44.
BFILE_ROW_8 = [6209, -41417, 48645, 376259, -1724861, 3233925,
               -2977609, 1078849, 40320]


def laguerre(n, x):
    return sum((Fraction((-1)**j * comb(n, j) * x**j, factorial(j))
                for j in range(n + 1)), Fraction(0))


def source_row(n):
    row = []
    for m in range(n + 1):
        a = (-1)**n * factorial(n) * sum(
            ((-1)**k * comb(n + 1, k) * laguerre(n, m - k + 1)
             for k in range(m + 1)), Fraction(0))
        assert a.denominator == 1, (n, m, a)
        row.append(a.numerator)
    return row


def eulerian(n):
    """Eulerian polynomial of permutations, via the ascent insertion recurrence."""
    e = [1]
    for size in range(1, n + 1):
        old = e
        e = [(k + 1) * (old[k] if k < len(old) else 0)
             + (size - k) * (old[k - 1] if k > 0 else 0)
             for k in range(size)]
    return e


def inverse_eulerian_row(n):
    """P_n = sum_j binom(n,j) n!/j! E_j(t) (t-1)^(n-j)."""
    out = [0] * (n + 1)
    for j in range(n + 1):
        scale, rem = divmod(factorial(n), factorial(j))
        assert rem == 0
        scale *= comb(n, j)
        for u, coeff in enumerate(eulerian(j)):
            for v in range(n - j + 1):
                out[u + v] += scale * coeff * comb(n - j, v) * (-1)**(n - j - v)
    return out


def content(row):
    return reduce(gcd, map(abs, row), 0)


def small_divisors(n):
    assert n >= 1
    return [d for d in range(1, isqrt(n) + 1) if n % d == 0]


def factorize(n):
    assert n >= 1
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def divisors_from_factors(factors):
    ds = [1]
    for p, exponent in factors.items():
        ds = [d * p**k for d in ds for k in range(exponent + 1)]
    return ds


def q_at_one(n):
    value = 0
    for j in range(n + 1):
        quotient, rem = divmod(factorial(n), factorial(j))
        assert rem == 0
        value += (-1)**(n + j) * comb(n, j) * quotient
    return value


def main():
    for n, expected in enumerate(SOURCE_ROWS):
        assert source_row(n) == expected
    print('FIDELITY: all source COMMENTS rows P_0 through P_7 match exactly')
    for n in range(21):
        row = source_row(n)
        assert row == inverse_eulerian_row(n), n
        assert sum(row) == factorial(n)
        assert row[-1] == factorial(n)
    print('INDEPENDENT REPRESENTATIONS: Laguerre sum and inverse Eulerian transform agree for n=0..20')
    print('FIDELITY: row sums and leading coefficients equal n! for n=0..20')
    row = source_row(8)
    assert row == BFILE_ROW_8
    assert content(row) == 1
    bezout = -6417 * row[0] - 962 * row[1]
    assert bezout == 1
    assert small_divisors(8) == [1, 2]
    assert small_divisors(9) == [1, 3]
    assert content(row) != len(small_divisors(8))
    assert content(row) != len(small_divisors(9))
    print('WITNESS: n=8 row =', row)
    print('BFILE: all nine coefficients match published indices 36..44')
    print('BEZOUT: -6417*6209 - 962*(-41417) =', bezout)
    print('REFUTATION: c_8=1, D(8)=2, A060872(8)=16, D(9)=2')
    for n in range(2, 258, 2):
        assert q_at_one(n) % 2 == 1
    print('PARITY SANITY: q_n(1) is odd for every even n=2..256 (exact integer sums)')
    for s in range(3):
        r = 12 * s + 7
        n = 2**r
        ds = [2**j for j in range(r + 1)]
        d_n = sum(d * d <= n for d in ds)
        assert d_n == 6 * s + 4 and d_n % 2 == 0
        m, rem = divmod(n + 1, 3)
        assert rem == 0 and gcd(3, m) == 1
        assert m % 8 == 3 and isqrt(m)**2 != m
        factors = factorize(n + 1)
        shifted_divisors = divisors_from_factors(factors)
        d_shifted = sum(d * d <= n + 1 for d in shifted_divisors)
        assert factors[3] == 1 and len(shifted_divisors) % 4 == 0
        assert d_shifted % 2 == 0
        print(f'FAMILY SANITY: s={s}, r={r}, n={n}, D(n)={d_n}, D(n+1)={d_shifted}, factors(n+1)={factors}')
    for s in range(101):
        r = 12 * s + 7
        assert pow(2, r, 9) == 2
        assert ((2**r + 1) // 3) % 8 == 3
    print('FAMILY MODULAR SANITY: r=12s+7 residues verified for s=0..100')
    print('SCOPE: finite exact checks support the separately written infinite-family argument; no Lean proof claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
