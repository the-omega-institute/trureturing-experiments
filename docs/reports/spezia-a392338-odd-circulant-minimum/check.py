#!/usr/bin/env python3
"""Exact finite checks for Spezia's OEIS A392338 odd-order conjecture.

Source: https://oeis.org/A392338/internal, revision 22.
Standard library only. These finite checks do not prove the universal claim.
"""
from fractions import Fraction
from itertools import permutations
from math import factorial


def det(a):
    """Fraction-free Bareiss elimination, including signed row pivoting."""
    b = [row[:] for row in a]
    n, previous, sign = len(b), 1, 1
    if n == 0:
        return 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            sign = -sign
        value = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * b[i][j] - b[i][k] * b[k][j]
                assert numerator % previous == 0
                b[i][j] = numerator // previous
            b[i][k] = 0
        previous = value
    return sign * b[-1][-1]


def det_fraction(a):
    """Independent exact Gaussian elimination over rational numbers."""
    b = [[Fraction(x) for x in row] for row in a]
    n, answer = len(b), Fraction(1)
    for k in range(n):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            answer = -answer
        value = b[k][k]
        answer *= value
        for i in range(k + 1, n):
            scale = b[i][k] / value
            for j in range(k, n):
                b[i][j] -= scale * b[k][j]
    return answer


def circulant(first):
    n = len(first)
    return [[first[(j - i) % n] for j in range(n)] for i in range(n)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def fidelity():
    # Rotating a first row only permutes columns. Thus fixing its 0 in position
    # 0 enumerates every absolute determinant, with n rotations removed.
    published = {2: 1, 3: 9, 4: 48, 5: 50, 6: 2205, 7: 147, 8: 1120}
    for n, expected in published.items():
        values = [abs(det(circulant((0,) + p))) for p in permutations(range(1, n))]
        assert len(values) == factorial(n - 1)
        actual = min(x for x in values if x)
        assert actual == expected, (n, actual, expected)
        if n % 2:
            bound = n * n * ((n - 1) // 2)
            assert all(x % bound == 0 for x in values)
        print(f'FIDELITY n={n} normalized_permutations={len(values)} minimum={actual}')
    example = [[3, 0, 2, 1], [1, 3, 0, 2], [2, 1, 3, 0], [0, 2, 1, 3]]
    assert example == circulant(example[0])
    assert abs(det(example)) == 48
    assert det_fraction(example) == det(example)
    print('FIDELITY source_displayed_order_4_example abs_det=48')


def attaining_family():
    for n in range(3, 28, 2):
        m = (n - 1) // 2
        first = [j if j % 2 == 0 else 2 * m - j for j in range(n)]
        assert sorted(first) == list(range(n))
        c = circulant(first)
        p = [[int(j == (i + 1) % n) for j in range(n)] for i in range(n)]
        ip = [[int(i == j) + p[i][j] for j in range(n)] for i in range(n)]
        lap = [[int(i == j) - p[i][j] for j in range(n)] for i in range(n)]
        rhs = [[lap[i][j] + 4 * m for j in range(n)] for i in range(n)]
        assert multiply(c, multiply(ip, ip)) == rhs
        for k in range(n):
            expected = 4 * m + int(k == 0) - int(k == 1)
            assert first[k] + 2 * first[(k - 1) % n] + first[(k - 2) % n] == expected
        target = m * n * n
        assert det(c) == det_fraction(c) == target
        assert det(ip) == 2
        assert det(lap) == 0
        assert det(rhs) == 4 * m * n * n
        # Independent rank-one determinant checks, including negative t and 0.
        for t in (-2, 0, 1, 4 * m):
            assert det([[x + t for x in row] for row in lap]) == t * n * n
        if n <= 11:
            for i in range(n):
                for j in range(n):
                    minor = [[lap[r][s] for s in range(n) if s != j]
                             for r in range(n) if r != i]
                    assert (-1) ** (i + j) * det(minor) == 1
        assert all(sum(row) == n * m for row in c)
        assert all(sum(col) == n * m for col in zip(*c))
        print(f'FAMILY n={n} det={target} det_I_plus_P=2 identity=OK rational_crosscheck=OK')


def row_column_sum_checks():
    count = 0
    for n in (3, 5, 7, 9):
        for m in (1, 2, 3):
            for seed in range(12):
                a = [[m for _ in range(n)] for _ in range(n)]
                # Arbitrary integral 2x2 cycle moves preserve both margins.
                for k in range(1, 2 * n):
                    i, j = (k + seed) % n, (k * k + seed) % n
                    r, s = (i + 1) % n, (j + 1) % n
                    delta = (seed + 3 * k) % 11 - 5
                    a[i][j] += delta
                    a[r][s] += delta
                    a[i][s] -= delta
                    a[r][j] -= delta
                assert all(sum(row) == n * m for row in a)
                assert all(sum(col) == n * m for col in zip(*a))
                reduced = [[a[i][j] - a[-1][j] for j in range(n - 1)]
                           for i in range(n - 1)]
                assert all(sum(col) % n == 0 for col in zip(*reduced))
                assert det(a) == n * m * det(reduced)
                assert det(reduced) % n == 0
                assert det(a) % (m * n * n) == 0
                count += 1
    print(f'GENERAL_MARGIN_CHECKS matrices={count} determinant_factorization=OK divisibility=OK')


if __name__ == '__main__':
    fidelity()
    attaining_family()
    row_column_sum_checks()
    print('ALL_CHECKS_OK')
