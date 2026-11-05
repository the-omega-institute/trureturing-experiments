#!/usr/bin/env python3
"""Exact checks for Mikhail Kurkov's conjecture in OEIS A396563 revision 15.

Only Python standard-library integer and Fraction arithmetic is used. Finite
checks validate transcription and proposed bridges; they are not a proof.
"""
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import factorial, prod


@cache
def c(n, k):
    """Unsigned first-kind Stirling numbers, by their integer recurrence."""
    if n == 0:
        return int(k == 0)
    if k <= 0 or k > n:
        return 0
    return (n - 1) * c(n - 1, k) + c(n - 1, k - 1)


@cache
def S(n, k):
    """Second-kind Stirling numbers, independently by their recurrence."""
    if n == 0:
        return int(k == 0)
    if k <= 0 or k > n:
        return 0
    return k * S(n - 1, k) + S(n - 1, k - 1)


def det(matrix):
    """Gaussian elimination over Q; det of the empty matrix is 1."""
    a = [[Fraction(x) for x in row] for row in matrix]
    d = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            d = -d
        x = a[i][i]
        d *= x
        for j in range(i + 1, len(a)):
            q = a[j][i] / x
            for t in range(i + 1, len(a)):
                a[j][t] -= q * a[i][t]
            a[j][i] = Fraction(0)
    return d


def B(n, k, m):
    return det([[c(n + m + i, n + j) for j in range(1, k + 1)]
                for i in range(1, k + 1)])


def chains(bound, size):
    if size == 0:
        yield ()
    else:
        for j in range(bound + 1):
            for tail in chains(j, size - 1):
                yield (j,) + tail


def literal_sum(n, k, m):
    """The source's nested sum, before any change of variables."""
    total = Fraction(0)
    for js in chains(n, m):
        term = Fraction((-1) ** sum(js))
        term *= prod((js[q] - js[p] + p - q) ** 2
                     for p in range(m) for q in range(p))
        for t, j in enumerate(js, start=1):
            falling = prod(n + t - 1 - a for a in range(j))
            term *= Fraction(falling * (n - j + t) ** (n + k),
                             factorial(j + m - t))
        total += term
    return total / prod(factorial(n + i) for i in range(m))


def subset_sum(n, k, m):
    """Independent enumeration after z_t=n-j_t+t, with source sign."""
    N = n + m
    total = Fraction(0)
    for zs in combinations(range(1, N + 1), m):
        exponent = m * n + m * (m + 1) // 2 - sum(zs)
        assert exponent >= 0
        term = Fraction((-1) ** exponent)
        term *= prod((zs[p] - zs[q]) ** 2
                     for p in range(m) for q in range(p))
        term *= prod(Fraction(z ** (n + k),
                              factorial(z - 1) * factorial(N - z)) for z in zs)
        total += term
    return total


def symmetric_coefficients(N, degree):
    e = [1] + [0] * degree
    h = [1] + [0] * degree
    for z in range(1, N + 1):
        for r in range(degree, 0, -1):
            e[r] += z * e[r - 1]
        for r in range(1, degree + 1):
            h[r] += z * h[r - 1]
    return e, h


def coefficient(a, r):
    return a[r] if 0 <= r < len(a) else 0


def main():
    # Verbatim numeric data from the source's EXAMPLE array, rows n=0,...,5.
    rows = [
        [1, 2, 4, 8, 16, 32, 64],
        [1, 11, 85, 575, 3661, 22631, 137845],
        [1, 35, 725, 11935, 173901, 2360295, 30673525],
        [1, 85, 3850, 129360, 3665151, 93355395, 2216707900],
        [1, 175, 15190, 929334, 46149327, 2000578965, 79000290160],
        [1, 322, 48804, 5012406, 405247689, 27912905328, 1717078544896],
    ]
    print('Source: Mikhail Kurkov, OEIS A396563, revision 15; 2026-06-01 13:43:04')
    for n, row in enumerate(rows):
        actual = [B(n, k, 2) for k in range(len(row))]
        assert actual == row, (n, actual, row)
        print(f'Fidelity row n={n}: {row}')
    print('FIDELITY_OK: all 42 displayed array values reproduced')
    cases = 0
    chain_count = 0
    for n in range(7):
        for k in range(7):
            for m in range(1, 6):
                N = n + m
                lhs, rhs = B(n, k, m), literal_sum(n, k, m)
                assert lhs == rhs, ('claim', n, k, m, lhs, rhs)
                assert subset_sum(n, k, m) == rhs, ('subset', n, k, m)
                degree = max(m + k, 1)
                e, h = symmetric_coefficients(N, degree)
                elementary = det([[coefficient(e, m + i - j) for j in range(k)]
                                  for i in range(k)])
                complete = det([[coefficient(h, k + i - j) for j in range(m)]
                               for i in range(m)])
                moments = det([[S(n + k + 1 + i + j, N) for j in range(m)]
                              for i in range(m)])
                assert lhs == elementary == complete, ('duality', n, k, m)
                assert rhs == (-1) ** (m * (m - 1) // 2) * moments, (
                    'moments', n, k, m)
                # Independent finite-type bijection check, including n=0.
                source_nodes = {tuple(n - j + t for t, j in enumerate(js, 1))
                                for js in chains(n, m)}
                expected_nodes = set(combinations(range(1, N + 1), m))
                assert source_nodes == expected_nodes, ('bijection', n, k, m)
                chain_count += len(source_nodes)
                cases += 1
    print(f'IDENTITY_OK: {cases} triples, 0<=n,k<=6 and 1<=m<=5')
    print(f'BRIDGES_OK: subset normalization, elementary/complete determinants, '
          f'Stirling moments, and {chain_count} chain instances')
    print('BOUNDARIES_OK: all n=0 and all k=0 cases in the stated box')
    assert B(2, 2, 2) == literal_sum(2, 2, 2) == 725
    print('EXAMPLE_OK: B(2,2,2) = literal_sum(2,2,2) = 725')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
