#!/usr/bin/env python3
"""Exact finite checks for Lang's OEIS A259056/A259057 conjectures.

Python standard library only. The matrices are filled from the source wedge,
then extended by its symmetry orbit. No determinant formula is used to fill
any entry. These finite checks are evidence, not an all-dimension proof.
"""
from fractions import Fraction
from math import factorial, prod

LIMIT = 32
SOURCE_DETS = [
    1, -3, -16, 60, 384, -1680, -12288, 60480, 491520, -2661120,
    -23592960, 138378240, 1321205760, -8302694400, -84557168640,
    564583219200, 6088116142080, -42908324659200, -487049291366400,
    3604299271372800,
]
SOURCE_QUOTIENTS = [
    -1, 20, -560, 20160, -887040, 46126080, -2767564800,
    188194406400, -14302774886400, 1201433090457600,
    -110531844322099200, 11053184432209920000, -1193743918678671360000,
    138474294566725877760000, -17170812526274008842240000,
]


def source_matrix(n):
    assert n >= 1
    b = [[None] * n for _ in range(n)]
    for i in range(1, (n + 1) // 2 + 1):
        for j in range(i, n + 2 - i):
            value = (i - 1) * n - (i - 1) ** 2 + j
            orbit = [(i, j), (j, i), (n + 1 - i, n + 1 - j),
                     (n + 1 - j, n + 1 - i)]
            for a, c in orbit:
                assert b[a - 1][c - 1] in (None, value)
                b[a - 1][c - 1] = value
    assert all(x is not None for row in b for x in row)
    assert {x for row in b for x in row} == set(range(1, (n + 1) ** 2 // 4 + 1))
    assert all(b[i][j] == b[j][i] == b[n - 1 - i][n - 1 - j]
               for i in range(n) for j in range(n))
    return b


def det_bareiss(matrix):
    """Fraction-free elimination with exact divisions and pivoting."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def det_rational(matrix):
    """Independent Gaussian elimination over fractions.Fraction."""
    a = [[Fraction(x) for x in row] for row in matrix]
    determinant = Fraction(1)
    for k in range(len(a)):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return Fraction(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            determinant = -determinant
        pivot = a[k][k]
        determinant *= pivot
        for i in range(k + 1, len(a)):
            ratio = a[i][k] / pivot
            for j in range(k, len(a)):
                a[i][j] -= ratio * a[k][j]
    return determinant


def mul(a, b):
    if not a:
        return []
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def reflection_basis(n):
    r = n // 2
    plus = r + n % 2
    p, q = [[0] * n for _ in range(n)], [[0] * n for _ in range(n)]
    for i in range(r):
        j = n - 1 - i
        p[i][i] = p[j][i] = 1
        p[i][plus + i], p[j][plus + i] = 1, -1
        q[i][i] = q[i][j] = Fraction(1, 2)
        q[plus + i][i], q[plus + i][j] = Fraction(1, 2), Fraction(-1, 2)
    if n % 2:
        p[r][r] = q[r][r] = 1
    assert mul(q, p) == identity(n) == mul(p, q)
    return p, q, plus


def min_kernel(f):
    assert f[0] == 0
    r = len(f) - 1
    k = [[f[min(i, j) + 1] for j in range(r)] for i in range(r)]
    delta = [f[i] - f[i - 1] for i in range(1, r + 1)]
    lower = [[int(j <= i) for j in range(r)] for i in range(r)]
    diagonal = [[delta[i] if i == j else 0 for j in range(r)] for i in range(r)]
    transpose = [list(row) for row in zip(*lower)]
    assert k == mul(mul(lower, diagonal), transpose)
    assert det_bareiss(lower) == 1
    assert det_bareiss(k) == prod(delta)
    return k


def check_dimension(n):
    b = source_matrix(n)
    determinant = det_bareiss(b)
    assert determinant == det_rational(b)
    if n <= len(SOURCE_DETS):
        assert determinant == SOURCE_DETS[n - 1]
    p, q, plus_size = reflection_basis(n)
    blocks = mul(mul(q, b), p)
    c = [row[:plus_size] for row in blocks[:plus_size]]
    d = [row[plus_size:] for row in blocks[plus_size:]]
    assert all(blocks[i][j] == 0 and blocks[j][i] == 0
               for i in range(plus_size) for j in range(plus_size, n))
    assert determinant == det_bareiss(c) * det_bareiss(d)
    reversed_d = [list(reversed(row)) for row in reversed(d)]
    if n % 2:
        m, r = (n + 1) // 2, (n - 1) // 2
        f = [0] + [(i - 1) * n - (i - 1) ** 2 + m for i in range(1, m + 1)]
        k = min_kernel(f)
        assert c == [[k[i][j] * (2 if j < r else 1) for j in range(m)]
                     for i in range(m)]
        assert reversed_d == [[-2 * min(i + 1, j + 1) for j in range(r)] for i in range(r)]
        min_kernel([0] + [-2 * i for i in range(1, r + 1)])
        assert det_bareiss(d) == (-2) ** r
        assert det_bareiss(c) == 4 ** r * factorial(m)
        assert determinant == (-1) ** (m - 1) * 8 ** (m - 1) * factorial(m)
    else:
        m = n // 2
        g = [0] + [2 * (i - 1) * n - 2 * (i - 1) ** 2 + n + 1
                   for i in range(1, m + 1)]
        assert c == min_kernel(g)
        assert reversed_d == [[-(2 * min(i + 1, j + 1) - 1)
                               for j in range(m)] for i in range(m)]
        min_kernel([0] + [-(2 * i - 1) for i in range(1, m + 1)])
        double_factorial = prod(range(1, 2 * m + 2, 2))
        assert det_bareiss(d) == (-1) ** m * 2 ** (m - 1)
        assert det_bareiss(c) == 2 ** (m - 1) * double_factorial
        assert determinant == (-1) ** m * 4 ** (m - 1) * double_factorial
        assert determinant % 3 == 0
        if m <= len(SOURCE_QUOTIENTS):
            assert determinant // 3 == SOURCE_QUOTIENTS[m - 1]
    return determinant


def main():
    assert source_matrix(1) == [[1]]
    assert source_matrix(2) == [[1, 2], [2, 1]]
    assert source_matrix(3) == [[1, 2, 3], [2, 4, 2], [3, 2, 1]]
    assert source_matrix(4) == [[1, 2, 3, 4], [2, 5, 6, 3], [3, 6, 5, 2], [4, 3, 2, 1]]
    assert min_kernel([0]) == []
    print('Arithmetic: integers and fractions.Fraction only; no floating point.')
    print('Fidelity: OEIS A259056 revision 7, A259057 revision 7.')
    print('Dimensions 1..32: symmetry-orbit construction, two determinant algorithms,')
    print('reflection basis and inverse, both block formulas, min-kernel factorizations.')
    determinants = [check_dimension(n) for n in range(1, LIMIT + 1)]
    print('Source matrices 1..4: matched.')
    print('Source determinant terms 1..20: matched.')
    print('Source quotient terms 1..15 (dimensions 2..30): matched.')
    print('Odd identity: m = 1..16; even identity and divisibility: m = 1..16.')
    print('Empty min-kernel: determinant = empty product = 1.')
    for n, d in enumerate(determinants, 1):
        print(f'N={n} det={d}')
    print('Finite checks only; no Lean proof or all-dimension certification.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
