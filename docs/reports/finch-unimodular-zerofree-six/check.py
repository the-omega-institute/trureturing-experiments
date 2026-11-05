#!/usr/bin/env python3
"""Exact checks for Finch arXiv:2605.03691v4, Conjecture 4.

Only Python's standard library is required. No floating point is used.
The finite family checks are evidence, not a proof for all dimensions.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def matmul(a, b):
    return [[sum(a[i][h] * b[h][j] for h in range(len(b)))
             for j in range(len(b))] for i in range(len(a))]


def det_bareiss(a):
    """Fraction-free elimination with exact-division assertions and row pivots."""
    b = [row[:] for row in a]
    n = len(b)
    if n == 0:
        return 1
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if b[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            b[k], b[pivot_row] = b[pivot_row], b[k]
            sign = -sign
        pivot = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = b[i][j] * pivot - b[i][k] * b[k][j]
                assert numerator % previous == 0
                b[i][j] = numerator // previous
            b[i][k] = 0
        previous = pivot
    return sign * b[-1][-1]


def det_permutations(a):
    """Independent Leibniz determinant; used on the order-six witness."""
    total = 0
    n = len(a)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = (-1) ** inversions
        for i in range(n):
            term *= a[i][p[i]]
        total += term
    return total


def inverse_gauss_jordan(a):
    """Independent rational inverse, with every arithmetic operation exact."""
    n = len(a)
    b = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if b[i][j])
        b[j], b[pivot] = b[pivot], b[j]
        divisor = b[j][j]
        b[j] = [x / divisor for x in b[j]]
        for i in range(n):
            if i != j:
                coefficient = b[i][j]
                b[i] = [x - coefficient * y for x, y in zip(b[i], b[j])]
    assert [row[:n] for row in b] == identity(n)
    return [row[n:] for row in b]


def max_abs(a):
    return max((abs(x) for row in a for x in row), default=0)


def zerofree(a):
    return all(x != 0 for row in a for x in row)


def structural_key(flat):
    # Finch orders all positive integers first, then the negative integers by magnitude.
    return tuple((x < 0, abs(x)) for x in flat)


def canonical_two(a):
    candidates = []
    for rp in permutations(range(2)):
        for cp in permutations(range(2)):
            for rs in product((-1, 1), repeat=2):
                for cs in product((-1, 1), repeat=2):
                    candidates.append(tuple(rs[i] * cs[j] * a[rp[i]][cp[j]]
                                            for i in range(2) for j in range(2)))
    return min(candidates, key=structural_key)


def fidelity_two(bound, expected):
    """Exhaust all zerofree 2 by 2 matrices with alpha = beta = bound."""
    alphabet = tuple(range(-bound, 0)) + tuple(range(1, bound + 1))
    orbits, raw_count = set(), 0
    for x in product(alphabet, repeat=4):
        a = [list(x[:2]), list(x[2:])]
        determinant = x[0] * x[3] - x[1] * x[2]
        if abs(determinant) != 1:
            continue
        b = [[determinant * x[3], -determinant * x[1]],
             [-determinant * x[2], determinant * x[0]]]
        if max_abs(a) != bound or max_abs(b) != bound:
            continue
        assert zerofree(a) and zerofree(b)
        assert matmul(a, b) == identity(2) == matmul(b, a)
        raw_count += 1
        orbits.add(canonical_two(a))
    assert orbits == expected, (bound, orbits, expected)
    print(f"FIDELITY n=2 alpha=beta={bound}: raw={raw_count}, orbits={len(orbits)}, "
          f"canonical={sorted(orbits, key=structural_key)}")


def family(signs, k):
    n = len(signs)
    a = [[k + (signs[i] if i == j else 0) for j in range(n)] for i in range(n)]
    b = [[-k * signs[i] * signs[j] + (signs[i] if i == j else 0)
          for j in range(n)] for i in range(n)]
    return a, b


def check_family(signs, k):
    assert sum(signs) == 0 and signs and k >= 2
    a, b = family(signs, k)
    n = len(signs)
    assert matmul(a, b) == identity(n) == matmul(b, a)
    determinant = (-1) ** signs.count(-1)
    assert det_bareiss(a) == determinant == det_bareiss(b)
    assert zerofree(a) and zerofree(b)
    assert max_abs(a) == max_abs(b) == k + 1


def check_polynomial_coefficients(signs):
    """Check coefficients of (D+tJ)(D-t ss^T) and the reverse over Z[t]."""
    n = len(signs)
    d = [[signs[i] if i == j else 0 for j in range(n)] for i in range(n)]
    j = [[1] * n for _ in range(n)]
    c = [[-signs[i] * signs[h] for h in range(n)] for i in range(n)]
    assert matmul(d, d) == identity(n)
    dc, jd, dj, cd = matmul(d, c), matmul(j, d), matmul(d, j), matmul(c, d)
    assert all(dc[i][h] + jd[i][h] == 0 for i in range(n) for h in range(n))
    assert all(dj[i][h] + cd[i][h] == 0 for i in range(n) for h in range(n))
    zero = [[0] * n for _ in range(n)]
    assert matmul(j, c) == zero == matmul(c, j)


def main():
    print("SOURCE arXiv:2605.03691v4 (12 June 2026), Sections 1 and 5")
    fidelity_two(2, {(1, 1, 1, 2)})
    fidelity_two(3, {(1, 2, 2, 3), (1, 2, 1, 3), (1, 1, 2, 3)})

    a, b = family([1, 1, 1, -1, -1, -1], 2)
    expected_a = [[3, 2, 2, 2, 2, 2], [2, 3, 2, 2, 2, 2], [2, 2, 3, 2, 2, 2],
                  [2, 2, 2, 1, 2, 2], [2, 2, 2, 2, 1, 2], [2, 2, 2, 2, 2, 1]]
    expected_b = [[-1, -2, -2, 2, 2, 2], [-2, -1, -2, 2, 2, 2], [-2, -2, -1, 2, 2, 2],
                  [2, 2, 2, -3, -2, -2], [2, 2, 2, -2, -3, -2], [2, 2, 2, -2, -2, -3]]
    assert a == expected_a and b == expected_b
    check_family([1, 1, 1, -1, -1, -1], 2)
    assert inverse_gauss_jordan(a) == b and inverse_gauss_jordan(b) == a
    assert det_permutations(a) == -1 == det_permutations(b)
    print("WITNESS A =", a)
    print("WITNESS B =", b)
    print("WITNESS AB=BA=I_6; det A=det B=-1 (Bareiss and all 720 permutations); "
          "exact Gauss-Jordan inverse=B; all entries nonzero; alpha=beta=3")

    count = 0
    for m in range(1, 9):
        signs = [1] * m + [-1] * m
        check_polynomial_coefficients(signs)
        for k in range(2, 11):
            check_family(signs, k)
            count += 1
    print(f"FAMILY canonical balanced signs: m=1..8, k=2..10, {count} pairs checked")
    print("POLYNOMIAL all coefficients of AB=BA=I in Z[t] checked separately for m=1..8")

    count = 0
    for m in range(1, 5):
        for plus in combinations(range(2 * m), m):
            signs = [1 if i in plus else -1 for i in range(2 * m)]
            for k in range(2, 5):
                check_family(signs, k)
                count += 1
    print(f"FAMILY all balanced sign placements: m=1..4, k=2..4, {count} pairs checked")

    a0, b0 = family([], 2)
    assert max_abs(a0) == max_abs(b0) == 0
    for m in range(1, 5):
        a1, b1 = family([1] * m + [-1] * m, 1)
        assert not zerofree(a1) and not zerofree(b1)
    au, bu = family([1, 1, 1, 1, -1, -1], 2)
    assert matmul(au, bu) != identity(6)
    print("BOUNDARIES empty dimension has maximum 0; k=1 has zero entries; unbalanced signs fail AB=I")
    print("SCOPE finite exact computation; no Lean proof or all-dimensions kernel theorem claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
