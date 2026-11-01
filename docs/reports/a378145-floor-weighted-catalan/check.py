#!/usr/bin/env python3
"""Exact finite checks for Schulte's OEIS A378145 floor-weighted conjecture.

Source: https://oeis.org/A378145, revision 11.
This program checks finite cases; it does not prove the all-index conjecture.
Only Python integers and fractions.Fraction are used as numerical evidence.
"""
from fractions import Fraction
from math import comb

MAX_N = 80
SOURCE_DATA = [
    1, 1, 1, 1, 2, 1, 2, 4, 3, 1, 5, 10, 8, 4, 1,
    14, 28, 23, 13, 5, 1, 42, 84, 70, 42, 19, 6, 1,
    132, 264, 222, 138, 68, 26, 7, 1, 429, 858, 726, 462,
    240, 102, 34, 8, 1, 1430, 2860, 2431, 1573, 847,
    385, 145, 43, 9, 1, 4862, 9724, 8294, 5434, 3003,
    1430, 583, 198, 53, 10, 1,
]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def catalans_by_recurrence(bound):
    values = [1]
    for n in range(bound):
        values.append(sum(values[i] * values[n - i] for i in range(n + 1)))
    return values


def source_binomial(n, k):
    """The source's piecewise closed formula, including its diagonal branch."""
    require(0 <= k <= n, ('domain', n, k))
    if k == n:
        return 1
    denominator = (2 * n - k) * (2 * n - k - 1)
    require(denominator > 0, ('denominator', n, k, denominator))
    value = Fraction(comb(2 * n - k, n) *
                     (n * (3 * k + 1) - 2 * k * (k + 1)), denominator)
    require(value.denominator == 1, ('nonintegral', n, k, value))
    return value.numerator


def truncated_mul(left, right, bound):
    result = [0] * (bound + 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j in range(bound + 1 - i):
            result[i + j] += a * right[j]
    return result


def main():
    cats = catalans_by_recurrence(MAX_N + 1)
    for n, value in enumerate(cats):
        require(value * (n + 1) == comb(2 * n, n), ('Catalan formula', n))
    rows = [[source_binomial(n, k) for k in range(n + 1)]
            for n in range(MAX_N + 1)]
    prefix = [a for row in rows for a in row][:len(SOURCE_DATA)]
    require(prefix == SOURCE_DATA, 'OEIS DATA mismatch')
    print(f'FIDELITY: all {len(SOURCE_DATA)} published DATA terms match the closed formula')
    print(f'SOURCE_ROW_6: {rows[6]}')

    # Independently implement the defining Riordan coefficient formula.
    f = [0] + cats[:MAX_N]  # x*C(x), reduced modulo x^(MAX_N+1).
    g = f.copy()
    g[0] += 1
    power = [1] + [0] * MAX_N
    entry_count = 0
    outside_count = 0
    for k in range(MAX_N + 1):
        column = truncated_mul(g, power, MAX_N)
        for n in range(MAX_N + 1):
            if k <= n:
                require(column[n] == rows[n][k], ('Riordan coefficient', n, k))
                entry_count += 1
            else:
                require(column[n] == 0, ('above diagonal', n, k))
                outside_count += 1
        power = truncated_mul(power, f, MAX_N)
    print(f'FIDELITY: defining Riordan series matches {entry_count} triangle entries')
    print(f'FIDELITY: {outside_count} coefficients above the diagonal vanish')

    recurrence_count = 0
    for n in range(MAX_N + 1):
        for k in range(2, n + 1):
            require(rows[n][k] == rows[n][k - 1] - rows[n - 1][k - 2],
                    ('source recurrence', n, k))
            recurrence_count += 1
    print(f'FIDELITY: all {recurrence_count} source recurrence instances match')

    weighted = [sum((1 + k // 2) * value for k, value in enumerate(row))
                for row in rows]
    for n, value in enumerate(weighted):
        require(value == cats[n + 1], ('floor-weighted conjecture', n, value, cats[n + 1]))
    print(f'TARGET: all {len(weighted)} identities for 0 <= n <= {MAX_N} hold exactly')
    print(f'FIRST_11_WEIGHTED_SUMS: {weighted[:11]}')
    print(f'LAST_WEIGHTED_SUM_n_{MAX_N}: {weighted[-1]}')
    print('SCOPE: finite exact checks only; the universal statement remains a proof obligation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
