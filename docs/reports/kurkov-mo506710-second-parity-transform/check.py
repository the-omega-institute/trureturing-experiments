#!/usr/bin/env python3
"""Exact finite checks for MathOverflow 506710's second new conjecture.

Python 3.12+, standard library only.  Finite checks are not a proof of the
unbounded claim.  T uses binary-block moment expansion, cross-checked against
the literal signed sum; P uses integer-valued Newton polynomials, cross-checked
against a separate direct recurrence table.
"""
from fractions import Fraction
from math import comb
import sys

assert sys.version_info >= (3, 12), sys.version
N, M = 30, 10


def evaluate_newton(coefficients, m):
    return sum(c * comb(m, i) for i, c in enumerate(coefficients) if i <= m)


def differences(values):
    result = []
    while values:
        result.append(values[0])
        values = [b - a for a, b in zip(values, values[1:])]
    return result


def source_p_rows(max_n):
    """P_n(m)=1+sum_{j=1}^m P_{n-1}(2j), in the binomial(m,i) basis.

    Discrete antidifferentiation uses
    sum_{j=1}^m binomial(j,i)=binomial(m+1,i+1)-[i=0].
    Thus the construction uses only the defining P recurrence.
    """
    rows = [[1]]
    for n in range(1, max_n + 1):
        d = differences([evaluate_newton(rows[-1], 2 * j) for j in range(n)])
        rows.append([1] + [d[i - 1] + (d[i] if i < n else 0)
                             for i in range(1, n + 1)])
    return rows


P_ROWS = source_p_rows(N)


def P(n, m):
    return evaluate_newton(P_ROWS[n], m)


# Independent direct implementation: each row is a prefix sum of the prior row.
D, W = 8, 20
previous = [1] * (W * 2 ** D + 1)
direct_p_count = 0
for n in range(D + 1):
    if n:
        current = [1]
        for m in range(1, W * 2 ** (D - n) + 1):
            current.append(current[-1] + previous[2 * m])
        previous = current
    for m in range(W + 1):
        assert P(n, m) == previous[m], ("direct P", n, m)
        direct_p_count += 1

# Additional recurrence samples include parameters beyond the transform grid.
p_recurrence_count = 0
for n in range(1, N + 1):
    assert P(n, 0) == 1
    for m in range(1, 81):
        assert P(n, m) == P(n, m - 1) + P(n - 1, 2 * m)
        p_recurrence_count += 1
assert all(P(0, m) == 1 for m in range(81))

# Source fidelity: OEIS A125790's first five displayed rows, first five columns.
expected_p = [
    [1, 1, 1, 1, 1],
    [1, 2, 3, 4, 5],
    [1, 4, 9, 16, 25],
    [1, 10, 35, 84, 165],
    [1, 36, 201, 656, 1625],
]
assert [[P(n, m) for m in range(5)] for n in range(5)] == expected_p

# T[n][k], including the shifted extra endpoint k=n+1.  Splitting off the
# highest binary bit gives the binomial expansion below (MO answer 504889).
T = [[0] * (N + 2) for _ in range(N + 1)]
for n in range(N + 1):
    T[n][0] = 1
for k in range(1, N + 2):
    for n in range(N + 1):
        T[n][k] = sum(comb(n, r) * 2 ** (k * (n - r)) * T[r][k - 1]
                      for r in range(n))

literal_t_count = 0
for k in range(11):
    literal = [0] * (N + 1)
    for j in range(2 ** k):
        sign = (-1) ** (j.bit_count() + k)
        power = 1
        for n in range(N + 1):
            literal[n] += sign * power
            power *= 2 * j + 1
    for n in range(N + 1):
        assert T[n][k] == literal[n], ("literal T", n, k)
        literal_t_count += 1
assert [[T[n][k] for k in range(n + 1)] for n in range(5)] == [
    [1], [1, 2], [1, 8, 16], [1, 26, 192, 384],
    [1, 80, 1696, 12288, 24576],
]
assert all(T[n][k] == 0 for n in range(N + 1)
           for k in range(n + 1, N + 2))

unsigned_count = signed_count = parity_count = 0
for n in range(N + 1):
    for m in range(M + 1):
        unsigned = sum((Fraction(T[n][k] * P(k, m), 2 ** (n * k))
                        for k in range(n + 1)), Fraction(0))
        signed = sum((Fraction(T[n][k] * P(k, m) * (-1) ** (n + k),
                               2 ** (n * k))
                      for k in range(n + 1)), Fraction(0))
        assert unsigned == (m + 2) ** n, ("unsigned", n, m, unsigned)
        assert signed == m ** n, ("signed", n, m, signed)
        unsigned_count += 1
        signed_count += 1
        if n:
            terms = [Fraction(T[n][k - 1] * P(k, m) *
                              (1 + (-1) ** (n + k)), 2 ** (n * k))
                     for k in range(1, n + 1)]
            assert sum(terms, Fraction(0)) == (m + 1) ** n, ("target", n, m)
            assert 1 + (-1) ** (n + n + 1) == 0  # omitted k=n+1 endpoint
            parity_count += 1

print("Arithmetic: Python integers and fractions.Fraction; no floating point")
print("Source fidelity: A125790 first 5 rows x 5 columns; A391910 rows 0..4")
print(f"P vs direct defining recurrence: {direct_p_count} pairs, n=0..8, m=0..20")
print(f"P recurrence checks: {p_recurrence_count} pairs, n=1..30, m=1..80")
print(f"T vs literal binary-sign sum: {literal_t_count} pairs, n=0..30, k=0..10")
print(f"Published unsigned transform: {unsigned_count} pairs, n=0..30, m=0..10")
print(f"Published signed transform: {signed_count} pairs, n=0..30, m=0..10")
print(f"Target shifted parity transform: {parity_count} pairs, n=1..30, m=0..10")
for n, m in [(3, 1), (4, 2), (30, 10)]:
    terms = [Fraction(T[n][k - 1] * P(k, m) *
                      (1 + (-1) ** (n + k)), 2 ** (n * k))
             for k in range(1, n + 1)]
    if n <= 4:
        print(f"Target n={n}, m={m}: terms=" + ", ".join(map(str, terms)) +
              f"; sum={sum(terms)}")
    else:
        print(f"Target n={n}, m={m}: sum={sum(terms)} = {m+1}^{n}")
print("Finite checks do not prove the unbounded assertion")
print("ALL_CHECKS_OK")
