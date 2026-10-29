#!/usr/bin/env python3
"""Exact finite checks for Hanna's full ternary coefficient conjecture, OEIS A078122.

The source array is constructed only from the shift-power recursion. Partition
coefficients are computed independently by ordinary unbounded coin-change DP.
Finite agreement checks transcription and a proposed proof route; it is no
proof of the universal conjecture. Python integers are used throughout.
"""

from itertools import product


SOURCE_ROWS = (
    (1,),
    (1, 1),
    (1, 3, 1),
    (1, 12, 9, 1),
    (1, 93, 117, 27, 1),
    (1, 1632, 3033, 1080, 81, 1),
    (1, 68457, 177507, 86373, 9801, 243, 1),
    (1, 7112055, 24975171, 15562314, 2371761, 88452, 729, 1),
    (1, 1879090014, 8786827629, 6734916423, 1291958181,
     64392813, 796797, 2187, 1),
)
SOURCE_CUBE = (
    (1, 0, 0, 0),
    (3, 1, 0, 0),
    (12, 9, 1, 0),
    (93, 117, 27, 1),
)


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def multiply(left, right):
    size = len(left)
    assert len(right) == size
    return [[sum(left[i][j] * right[j][k] for j in range(size))
             for k in range(size)] for i in range(size)]


def matrix_power(matrix, exponent):
    result = identity(len(matrix))
    for _ in range(exponent):
        result = multiply(result, matrix)
    return result


def source_matrix(base, last_row):
    """M^base(i,j)=M(i+1,j+1), first column and diagonal one."""
    rows = [[1]]
    for n in range(1, last_row + 1):
        prefix = [row + [0] * (n - len(row)) for row in rows]
        shifted = matrix_power(prefix, base)
        rows.append([1] + [shifted[n - 1][k - 1] for k in range(1, n)] + [1])
    size = last_row + 1
    return [row + [0] * (size - len(row)) for row in rows]


def partition_tables(base, max_d, limit):
    """P_d(N): counts of parts 1, base, ..., base**d, by coin change."""
    counts = [1] + [0] * limit
    result = []
    for d in range(max_d + 1):
        coin = base ** d
        for total in range(coin, limit + 1):
            counts[total] += counts[total - coin]
        result.append(counts.copy())
    return result


def tuple_count(base, d, total):
    """A second, direct multiplicity enumeration for tiny coefficient checks."""
    if d == 0:
        return 1
    coin = base ** d
    return sum(tuple_count(base, d - 1, total - copies * coin)
               for copies in range(total // coin + 1))


def main():
    ternary = source_matrix(3, 8)
    for n, expected in enumerate(SOURCE_ROWS):
        assert tuple(ternary[n][:n + 1]) == expected, (n, ternary[n], expected)
    first_four = [row[:4] for row in ternary[:4]]
    assert tuple(map(tuple, matrix_power(first_four, 3))) == SOURCE_CUBE
    print('FIDELITY: all 45 OEIS DATA entries through row 8 and the displayed 4x4 cube match')
    print('TERNARY_ROW_4:', tuple(ternary[4][:5]))
    print('TERNARY_ROW_5:', tuple(ternary[5][:6]))

    entries = row_sums = convolution = semigroup = recurrences = tiny_counts = 0
    for base in (2, 3, 4, 5):
        tables = partition_tables(base, 8 if base == 3 else 7,
                                  max(base ** 7, 10 * base ** 5, 3 ** 8 if base == 3 else 0))

        def f(d, t):
            assert d >= 0 and t >= 0
            if d == 0:
                return 1
            if t == 0:
                return 0
            target = base ** d * (t - 1)
            assert target < len(tables[d]), (base, d, t, target)
            return tables[d][target]

        matrix = source_matrix(base, 7)
        for n in range(8):
            for k in range(n + 1):
                degree = base ** n - base ** (n - k)
                assert matrix[n][k] == tables[n - k][degree], (base, n, k)
                entries += 1
            assert sum(matrix[n]) == tables[n][base ** n], (base, n)
            row_sums += 1

        for d, r, u in product(range(6), range(6), range(6)):
            rhs = sum(f(d - h, base ** h * r) * f(h, u) for h in range(d + 1))
            assert f(d, r + u) == rhs, (base, d, r, u)
            convolution += 1

        for d, t in product(range(1, 6), range(10)):
            assert f(d, t + 1) == f(d, t) + f(d - 1, base * t + 1)
            recurrences += 1

        def kernel(s):
            return [[f(n - k, s * base ** k) if k <= n else 0
                     for k in range(6)] for n in range(6)]

        assert kernel(0) == identity(6)
        for s, t in product(range(5), repeat=2):
            composed = multiply(kernel(s), kernel(t))
            expected = kernel(s + t)
            assert composed == expected, (base, s, t)
            semigroup += 21
        assert matrix_power(kernel(1), base) == kernel(base)

        for d, total in product(range(4), range(31)):
            assert tables[d][total] == tuple_count(base, d, total)
            tiny_counts += 1

        if base == 3:
            for n in range(9):
                for k in range(n + 1):
                    assert ternary[n][k] == tables[n - k][3 ** n - 3 ** (n - k)]

    print(f'COEFFICIENT_IDENTITY: {entries} entries, bases 2..5, rows 0..7')
    print(f'ROW_SUM_IDENTITY: {row_sums} rows, bases 2..5, rows 0..7')
    print(f'SCALED_CONVOLUTION: {convolution} identities, d,r,u=0..5, bases 2..5')
    print(f'FINITE_DIFFERENCE: {recurrences} identities, d=1..5, t=0..9, bases 2..5')
    print(f'SEMIGROUP: {semigroup} lower triangular entries, s,t=0..4, size 6, bases 2..5')
    print(f'INDEPENDENT_DP_FIDELITY: {tiny_counts} direct multiplicity counts, d=0..3, N=0..30')
    print('BOUNDARIES: F_0(0)=1, F_d(0)=0 for d>0, and N_0=I checked')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
