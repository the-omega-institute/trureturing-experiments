#!/usr/bin/env python3
"""Exact bounded checks for OEIS A336186 and the proposed two-orbit route.

Only Python 3's standard library is needed. No floating point is used.
These checks support preregistration; they do not prove the infinite claim.
"""
from math import comb, isqrt

SOURCE_TERMS = [1, 17, 127, 1871, 13969, 205793, 1536463,
                22635359, 168996961, 2489683697]
T = ((55, 252), (12, 55))
U = ((55, -252), (-12, 55))
I = ((1, 0), (0, 1))
SEEDS = [(6, 1), (78, 17)]
SCAN_BOUND = 2_000_000
BASE_TABLE = [
    (1, 6, 0), (3, 14, 8), (5, 23, 11), (7, 32, 20),
    (9, 41, 35), (11, 50, 56), (13, 59, 83), (15, 68, 116),
    (17, 78, 0), (19, 87, 27), (21, 96, 60), (23, 105, 99),
    (25, 114, 144), (27, 123, 195), (29, 132, 252), (31, 142, 32),
    (33, 151, 83), (35, 160, 140), (37, 169, 203), (39, 178, 272),
    (41, 187, 347), (43, 197, 35), (45, 206, 104),
]


def board(i, j):
    """Same formula as the binding Lean conventions; no truncated subtraction."""
    return comb(i + j + 1, 2) + i + 1 if i + j else 1


def block_sum(m):
    return sum(board(i, j) for i in range(m) for j in range(m))


def quartic(m):
    numerator = m * m * (7 * m * m + 5)
    assert numerator % 12 == 0
    return numerator // 12


def square(n):
    r = isqrt(n)
    return r * r == n


def mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def transpose(a):
    return tuple(zip(*a))


def apply(a, pair):
    return tuple(sum(a[i][j] * pair[j] for j in range(2)) for i in range(2))


def norm(pair):
    x, y = pair
    return x * x - 21 * y * y


def main():
    displayed = [[1, 2, 4, 7, 11, 16], [3, 5, 8, 12, 17],
                 [6, 9, 13, 18], [10, 14, 19], [15, 20], [21]]
    for i, row in enumerate(displayed):
        assert [board(i, j) for j in range(len(row))] == row
    assert block_sum(1) == 1
    assert block_sum(17) == 48841 == 221**2
    assert [board(0, j) for j in range(17)] == [
        1, 2, 4, 7, 11, 16, 22, 29, 37, 46, 56, 67, 79, 92, 106, 121, 137]
    assert [board(i, 0) for i in range(17)] == [
        1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153]
    print('FIDELITY: displayed board and source example S(17)=48841=221^2')
    for m in range(65):
        assert block_sum(m) == quartic(m)
    print('BLOCK_SUM: direct double sums equal m^2*(7*m^2+5)/12 for 0<=m<=64')

    q = ((1, 0), (0, -21))
    assert mul(T, U) == mul(U, T) == I
    assert mul(mul(transpose(T), q), T) == q
    assert mul(mul(transpose(U), q), U) == q
    assert mul(T, T) == ((6049, 27720), (1320, 6049))
    assert mul(T, T) == tuple(tuple(110*T[i][j]-I[i][j]
                                   for j in range(2)) for i in range(2))
    assert apply(T, SEEDS[0]) == (582, 127)
    assert apply(T, SEEDS[1]) == (8574, 1871)
    assert all(norm(seed) == 15 for seed in SEEDS)
    print('MATRICES: inverses, norm preservation, T^2=110*T-I, seed successors')

    assert [y for y, _, _ in BASE_TABLE] == list(range(1, 47, 2))
    for y, r, d in BASE_TABLE:
        assert 21*y*y+15 == r*r+d and 0 <= d < 2*r+1
        assert isqrt(21*y*y+15) == r
    bases = [(isqrt(21*y*y+15), y) for y in range(1, 47)
             if square(21*y*y+15)]
    assert bases == SEEDS
    print('BASE_CASES: all 46 positive sizes checked; only (6,1) and (78,17)')
    print('BASE_TABLE:', BASE_TABLE)

    generated_pairs = []
    a, b = SEEDS
    for _ in range(100):
        assert a[0] < b[0] < apply(T, a)[0]
        assert a[1] < b[1] < apply(T, a)[1]
        generated_pairs.extend([a, b])
        a, b = apply(T, a), apply(T, b)
    terms = [y for _, y in generated_pairs]
    assert terms == sorted(set(terms))
    assert terms[:10] == SOURCE_TERMS
    for x, y in generated_pairs:
        assert norm((x, y)) == 15 and x % 6 == 0 and y % 2 == 1
        assert quartic(y) == (y*x//6)**2
        if y >= 47:
            X, Y = apply(U, (x, y))
            assert 0 < X and 0 < Y < y and norm((X, Y)) == 15
    print('SOURCE_TERMS:', terms[:10])
    print('ORBIT_CHECK: 200 terms strictly sorted, valid, and descend to the two seeds')

    scan = []
    for m in range(1, SCAN_BOUND+1):
        s2 = quartic(m)
        n2 = 21*m*m+15
        has_square, has_norm = square(s2), square(n2)
        assert has_square == has_norm
        if has_square:
            scan.append(m)
            pair = (isqrt(n2), m)
            if m >= 47:
                X, Y = apply(U, pair)
                assert 0 < X and 0 < Y < m
                assert norm((X, Y)) == 15
            while pair[1] >= 47:
                pair = apply(U, pair)
            assert pair in SEEDS
    assert scan == [m for m in terms if m <= SCAN_BOUND]
    print(f'EXHAUSTIVE_SCAN: 1<=m<={SCAN_BOUND}; valid sizes={scan}')

    coeff = [0] + terms
    for n in range(5, len(coeff)):
        assert coeff[n] == 110*coeff[n-2] - coeff[n-4]
    denominator = {0: 1, 2: -110, 4: 1}
    numerator = {1: 1, 2: 17, 3: 17, 4: 1}
    for n in range(len(coeff)):
        value = sum(d * coeff[n-k] for k, d in denominator.items() if k <= n)
        assert value == numerator.get(n, 0)
    print('RECURRENCE: all indices 5..200; GENERATING_FUNCTION: coefficients 0..200')
    print('SCOPE: exact finite checks; infinite completeness remains a proof obligation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
