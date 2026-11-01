#!/usr/bin/env python3
"""Exact source-fidelity and finite identity checks for OEIS A395726.

Source: Mikhail Kurkov, approved revision 21, https://oeis.org/A395726/internal.
The published DATA and displayed rows below are source data, not generated oracles.
Python 3 standard library only; no floating-point arithmetic is used.
"""
from fractions import Fraction
from math import comb, factorial, isqrt

SOURCE_DATA = [
    1,0,1,0,1,1,0,1,2,1,0,1,4,3,1,0,1,8,7,4,1,0,1,16,15,13,5,1,
    0,1,32,31,40,19,6,1,0,1,64,63,121,65,25,7,1,0,1,128,127,364,211,
    90,35,8,1,0,1,256,255,1093,665,301,155,45,9,1,0,1,512,511,3280,
    2059,966,651,220,55,10,1,
]
SOURCE_ROWS = [
    [1,1,1,1,1,1,1,1,1,1],
    [0,1,2,3,4,5,6,7,8,9],
    [0,1,4,7,13,19,25,35,45,55],
    [0,1,8,15,40,65,90,155,220,285],
    [0,1,16,31,121,211,301,651,1001,1351],
    [0,1,32,63,364,665,966,2667,4368,6069],
    [0,1,64,127,1093,2059,3025,10795,18565,26335],
    [0,1,128,255,3280,6305,9330,43435,77540,111645],
    [0,1,256,511,9841,19171,28501,174251,320001,465751],
]
N_MAX, K_MAX = 20, 300


def f(k):
    """Literal source PARI formula on its positive-index domain."""
    return 0 if k == 0 else (isqrt(8 * k - 7) + 1) // 2


def triangular_block(k):
    """Independent implementation using cumulative block lengths."""
    m, endpoint = 0, 0
    while endpoint < k:
        m += 1
        endpoint += m
    return m


def position(k):
    return k - comb(f(k), 2)


def recurrence_table(n_max, k_max):
    a = [[1] * (k_max + 1)]
    for n in range(1, n_max + 1):
        row = [0] * (k_max + 1)
        for k in range(1, k_max + 1):
            m = triangular_block(k)
            assert 1 <= m <= k
            row[k] = row[k - m] + m * a[n - 1][k]
        a.append(row)
    return a


def source_list_column(n_max, k):
    """Direct translation of the source's PARI list_column, with zero-based lists."""
    v1 = [0] * f(k + 1)
    v1[0] = k
    v2 = [1] + [0] * n_max
    for i in range(1, len(v1)):
        v1[i] = v1[i - 1] - f(v1[i - 1])
    v1.reverse()
    for i in range(1, len(v1)):
        v3 = [1] + [0] * n_max
        weight = v1[i] - v1[i - 1]
        for n in range(1, n_max + 1):
            v3[n] = v2[n] + weight * v3[n - 1]
        v2 = v3
    return v2


def stirling_table(limit):
    s = [[0] * (limit + 1) for _ in range(limit + 1)]
    s[0][0] = 1
    for n in range(1, limit + 1):
        for h in range(1, n + 1):
            s[n][h] = h * s[n - 1][h] + s[n - 1][h - 1]
    return s


def conjecture_one(n, k):
    m, r = f(k), position(k)
    return sum((
        Fraction(comb(m - 1, j - 1) * (m - r - j) * (-1) ** (m + j + 1))
        * Fraction(j) ** (n + m - 2)
        for j in range(1, m + 1)
    ), Fraction(0)) / factorial(m - 1)


def column_weights(k):
    m, r = f(k + 1), position(k + 1)
    return [j + int(r > m - j) for j in range(m)]


def inverse_product_coefficients(weights, degree):
    """Multiply exact geometric series modulo x^(degree+1)."""
    coefficients = [1] + [0] * degree
    for weight in weights:
        new = [1] + [0] * degree
        for n in range(1, degree + 1):
            new[n] = coefficients[n] + weight * new[n - 1]
        coefficients = new
    return coefficients


def multiply_linear_factors(weights):
    coefficients = [1]
    for weight in weights:
        new = coefficients + [0]
        for i, coefficient in enumerate(coefficients):
            new[i + 1] -= weight * coefficient
        coefficients = new
    return coefficients


def main():
    a = recurrence_table(N_MAX, K_MAX)
    assert all(f(k) == triangular_block(k) for k in range(K_MAX + 2))
    assert all(1 <= position(k) <= f(k) for k in range(1, K_MAX + 2))
    print('block conventions: 302 indices checked by independent integer constructions')

    flattened = [a[d - k][k] for d in range(12) for k in range(d + 1)]
    assert flattened == SOURCE_DATA
    assert [row[:10] for row in a[:9]] == SOURCE_ROWS
    print('source fidelity: all 78 DATA terms and all 90 displayed array entries match')

    for k in range(K_MAX + 1):
        assert source_list_column(N_MAX, k) == [a[n][k] for n in range(N_MAX + 1)]
    print('source PARI column program: 6321 recurrence entries match (n=0..20, k=0..300)')

    s = stirling_table(N_MAX + f(K_MAX))
    for n in range(N_MAX + 1):
        for k in range(1, K_MAX + 1):
            m, r = f(k), position(k)
            assert a[n][k] == s[n + m - 1][m - 1] + r * s[n + m - 1][m]
            assert Fraction(a[n][k]) == conjecture_one(n, k), ('C1', n, k)
    print('triangular-block identification: 6300 pairs match (n=0..20, k=1..300)')
    print('Conjecture 1: 6300 pairs match, including n=0, m=1, j=1 with exponent -1')

    for k in range(K_MAX + 1):
        weights = column_weights(k)
        assert inverse_product_coefficients(weights, N_MAX) == [a[n][k] for n in range(N_MAX + 1)]
        denominator = multiply_linear_factors(weights)
        for n in range(N_MAX + 1):
            coefficient = sum(denominator[i] * a[n - i][k]
                              for i in range(min(n, len(denominator) - 1) + 1))
            assert coefficient == int(n == 0), ('C2 product', n, k)
    print('Conjecture 2: 6321 coefficients match; original denominator times column is 1 mod x^21')

    for n in range(1, N_MAX + 1):
        total = 0
        assert a[n][0] == total
        for k in range(1, K_MAX + 1):
            total += s[n + f(k) - 1][f(k)]
            assert a[n][k] == total, ('C3', n, k)
    print('Conjecture 3: 6020 pairs match (n=1..20, k=0..300)')

    # This non-target source line is printed without a conjecture label.
    # Its discrepancy is reported, never used as a premise or a target refutation.
    assert a[1][comb(2, 2)] == 1 and s[3][2] == 3
    print('separate unnumbered source line: n=1, k=2 gives recurrence 1 versus printed Stirling2 3')
    for k in [4, 6]:
        print(f'column {k}, n=0..8: {[a[n][k] for n in range(9)]}')
    print('finite exact agreement supports fidelity; it is not a proof of the unbounded claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
