#!/usr/bin/env python3
"""Exact finite checks for Peter Bala's all-column conjecture, OEIS A337389.

Source: https://oeis.org/A337389, revision #40.
This checks finite truncations; it does not prove the unbounded conjecture.
Run with Python 3; only the standard library is required.
"""
from fractions import Fraction as Q
from math import comb

N = 30
PARAMETERS = list(map(Q, range(31))) + list(map(Q, (101, 1000000, -1, -4))) + [Q(1, 2)]
# 17 is already included in range(31): 36 distinct parameters, 1116 pairs.


def series(*coefficients):
    return [Q(coefficients[i]) if i < len(coefficients) else Q(0) for i in range(N + 1)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scale(c, a):
    return [c * x for x in a]


def mul(a, b):
    return [sum((a[j] * b[n - j] for j in range(n + 1)), Q(0)) for n in range(N + 1)]


def inverse(a):
    assert a[0] != 0
    b = [1 / a[0]]
    for n in range(1, N + 1):
        b.append(-sum((a[j] * b[n - j] for j in range(1, n + 1)), Q(0)) / a[0])
    assert mul(a, b) == series(1)
    return b


def sqrt_one(a):
    assert a[0] == 1
    b = [Q(1)]
    for n in range(1, N + 1):
        b.append((a[n] - sum((b[j] * b[n - j] for j in range(1, n)), Q(0))) / 2)
    assert mul(b, b) == a
    return b


def source_column(k):
    delta = series(1, -2 * (k + 4), (k - 4) ** 2)
    inner = sqrt_one(delta)
    radicand = scale(Q(1, 2), mul(add(series(1, k - 4), inner), inverse(delta)))
    return sqrt_one(radicand)


def diagonal_column(k):
    phi = mul(mul(series(1, k - 2, 1), series(1, 2, 1)), inverse(series(1, -2, 1)))
    power = series(1)
    diagonal = []
    for n in range(N + 1):
        diagonal.append(power[n])
        if n < N:
            power = mul(power, phi)
    return diagonal


def binomial_column(k):
    return [sum((k ** (n - j) * comb(2 * j, j) * comb(2 * n, 2 * j)
                 for j in range(n + 1)), Q(0)) for n in range(N + 1)]


def main():
    assert len(PARAMETERS) == len(set(PARAMETERS)) == 36
    published_rows = [
        [1, 1, 1, 1, 1, 1],
        [2, 3, 4, 5, 6, 7],
        [6, 19, 34, 51, 70, 91],
        [20, 141, 328, 587, 924, 1345],
        [70, 1107, 3334, 7123, 12870, 20995],
        [252, 8953, 34904, 89055, 184756, 337877],
    ]
    columns = {k: source_column(k) for k in PARAMETERS}
    actual_rows = [[columns[Q(k)][n] for k in range(6)] for n in range(6)]
    assert actual_rows == published_rows, 'source example fidelity'
    published_antidiagonals = [
        1, 1, 2, 1, 3, 6, 1, 4, 19, 20, 1, 5, 34, 141, 70, 1, 6, 51, 328, 1107,
        252, 1, 7, 70, 587, 3334, 8953, 924, 1, 8, 91, 924, 7123, 34904, 73789,
        3432, 1, 9, 114, 1345, 12870, 89055, 372436, 616227, 12870, 1, 10, 139,
        1856, 20995, 184756, 1135005, 4027216, 5196627, 48620,
    ]
    actual_flat = [columns[Q(s - n)][n] for s in range(10) for n in range(s + 1)]
    assert actual_flat == published_antidiagonals, 'source antidiagonal fidelity'
    print('Source: OEIS A337389 revision #40; Bala final FORMULA conjecture.')
    print('Arithmetic: fractions.Fraction and integers only; no floating point.')
    print('Fidelity: all 36 displayed entries and all 55 published prefix entries match.')
    print('Reproduced rows, n=0..5 and k=0..5:')
    for row in actual_rows:
        print(' '.join(str(x) for x in row))
    recurrence_checks = 0
    for k in PARAMETERS:
        source = columns[k]
        diagonal = diagonal_column(k)
        finite_sum = binomial_column(k)
        for n in range(N + 1):
            assert source[n] == diagonal[n] == finite_sum[n], (n, k, source[n], diagonal[n], finite_sum[n])
        assert source[0] == 1 and source[1] == k + 2
        for n in range(2, N + 1):
            lhs = n * (2 * n - 1) * (4 * n - 5) * source[n]
            rhs = ((4 * n - 3) * (4 * (k + 4) * n ** 2 - 6 * (k + 4) * n + k + 6) * source[n - 1]
                   - (k - 4) ** 2 * (n - 1) * (2 * n - 3) * (4 * n - 1) * source[n - 2])
            assert lhs == rhs, ('source recurrence', n, k)
            recurrence_checks += 1
    print('Parameters: ' + ', '.join(map(str, PARAMETERS)))
    print(f'Identity: {len(PARAMETERS) * (N + 1)} pairs, n=0..{N}, three independently specified expressions agree.')
    print(f'Additional source fidelity: {recurrence_checks} recurrence instances match exactly.')
    print('Finite evidence only: no unbounded proof or Lean settlement is claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
