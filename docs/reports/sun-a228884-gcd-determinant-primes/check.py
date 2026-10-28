#!/usr/bin/env python3
"""Exact, offline checks for Zhi-Wei Sun's OEIS A228884 clause (ii).

Source: https://oeis.org/A228884/internal, revision #17.
Only Python's standard library is required. No floating point is used.
The matrix determinant is computed from its entries, never defined by a spectrum.
"""

from fractions import Fraction
from math import gcd, prod


SOURCE_DATA = (
    1, 1, 3, 20, 128, 2304, 10800, 606528, 3932160, 141087744,
    1289945088, 210000000000, 335544320000, 222902511206400,
    804545281732608, 39137889484800000, 972777519512027136,
    608742554432415203328, 391804906912468697088,
    1455817098785971890290688, 968232702940866945220608,
)


def factor(n):
    assert n > 0
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def is_prime(n):
    return n >= 2 and factor(n) == {n: 1}


def phi(n):
    assert n > 0
    value = n
    for p in factor(n):
        value = value // p * (p - 1)
    return value


def matrix(n):
    # The source uses labels 1,...,n, with integer subtraction i-j.
    return [[gcd(i - j, n) for j in range(1, n + 1)]
            for i in range(1, n + 1)]


def bareiss(a):
    a = [row[:] for row in a]
    n = len(a)
    if n == 0:
        return 1
    previous, sign = 1, 1
    for k in range(n - 1):
        pivot_row = next((r for r in range(k, n) if a[r][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0, (n, k, i, j, numerator, previous)
                a[i][j] = quotient
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def rational_det(a):
    a = [[Fraction(x) for x in row] for row in a]
    value = Fraction(1)
    n = len(a)
    for j in range(n):
        row = next((r for r in range(j, n) if a[r][j]), None)
        if row is None:
            return 0
        if row != j:
            a[j], a[row] = a[row], a[j]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for k in range(j, n):
            a[j][k] /= pivot
        for i in range(j + 1, n):
            multiple = a[i][j]
            for k in range(j, n):
                a[i][k] -= multiple * a[j][k]
    assert value.denominator == 1
    return value.numerator


def spectrum(n):
    return [sum(d * phi(n // d) for d in range(1, n + 1)
                if gcd(j, n) % d == 0) for j in range(n)]


def local_det(p, a):
    # Frequency 0 followed by valuations 0,...,a-1.
    zero = p ** (a - 1) * ((a + 1) * p - a)
    return zero * prod(((s + 1) * (p - 1) * p ** (a - 1))
                       ** phi(p ** (a - s)) for s in range(a))


def factorized_det(n):
    return prod(local_det(p, a) ** (n // p ** a)
                for p, a in factor(n).items())


def main():
    determinants = {}
    composites = 0
    print('Source: OEIS A228884, revision #17; offset 0')
    print('Arithmetic: Python integers and fractions.Fraction only')
    for n in range(61):
        actual = bareiss(matrix(n))
        determinants[n] = actual
        if n <= 12:
            assert rational_det(matrix(n)) == actual, n
        if n < len(SOURCE_DATA):
            assert actual == SOURCE_DATA[n], (n, actual, SOURCE_DATA[n])
        if n:
            eigenvalues = spectrum(n)
            assert all(x > 0 for x in eigenvalues), n
            assert actual == prod(eigenvalues) == factorized_det(n), n
        if is_prime(n):
            # Source COMMENTS: a(p)=(p-1)^(p-1)*(2p-1).
            assert actual == (n - 1) ** (n - 1) * (2 * n - 1), n
        factors = factor(actual)
        if n >= 4 and not is_prime(n):
            composites += 1
            assert all(q < n for q in factors), (n, factors)
        print(f'n={n:2d} determinant={actual} prime_factors={factors}')
    pairs = 0
    for u in range(1, 61):
        for v in range(1, 60 // u + 1):
            if gcd(u, v) == 1:
                pairs += 1
                assert determinants[u * v] == determinants[u] ** v * determinants[v] ** u
                residues = [(i % u, i % v) for i in range(u * v)]
                assert len(set(residues)) == u * v
                for i in range(u * v):
                    for j in range(u * v):
                        assert gcd(i - j, u * v) == gcd(i - j, u) * gcd(i - j, v)
    local_cases = 0
    for p in range(2, 98):
        if is_prime(p):
            for a in range(2, 11):
                local_cases += 1
                b = (a + 1) * p - a
                assert 0 < b <= p ** a
                assert p * b - ((a + 2) * p - (a + 1)) == (a + 1) * (p - 1) ** 2
                assert all(q < p ** a for q in factor(b))
    assert determinants[2] == 3 and determinants[3] == 20
    print('Fidelity: all 21 source DATA terms (n=0..20) and every prime-index formula through 60 agree')
    print('Independent determinant cross-check: rational Gaussian elimination, n=0..12')
    print('Three determinant methods agree for every n=1..60; empty determinant n=0 equals 1')
    print(f'Composite prime-support assertion verified for all {composites} composite n in 4..60')
    print(f'CRT entry identities and determinant factorization: {pairs} ordered coprime pairs uv<=60')
    print(f'Local bound and recurrence: {local_cases} cases, prime p<=97 and exponent 2<=a<=10')
    print('Finite checks do not prove the unbounded conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
