#!/usr/bin/env python3
"""Exact checks for Oldekop Conjecture 1.8 and the paper's ML-degree tables.
All matrix ranks and determinants are SymPy exact integer computations.
"""
from itertools import combinations
from random import Random
from sympy import Matrix


def scaling(d, entry):
    return Matrix(d, d, lambda i, j: 0 if i == j else entry(i, j))


def all_ones(d):
    return scaling(d, lambda i, j: 1)


def rank_three(d):
    # The paper's rank-three construction uses squared differences.
    return scaling(d, lambda i, j: (i - j) ** 2)


def principal(M, S):
    return M.extract(S, S)


def parity_sum(M):
    d = M.rows
    total = 0
    for mask in range(1 << d):
        S = [i for i in range(d) if mask & (1 << i)]
        rank = principal(M, S).rank()
        total += (-1) ** (len(S) + 1) * (rank % 2)
    return total


def nonzero_principal_count(M, size):
    return sum(principal(M, list(S)).det() != 0 for S in combinations(range(M.rows), size))


def rank_three_submatrix_count(M, size):
    return sum(principal(M, list(S)).rank() == 3 for S in combinations(range(M.rows), size))


def fixed_rank_matrix(d, r):
    data = {
        (5, 4): [[0, -161, 320, 60, -66], [-161, 0, -69, -418, -12],
                 [320, -69, 0, 16, 198], [60, -418, 16, 0, -147],
                 [-66, -12, 198, -147, 0]],
        (6, 4): [[0, -161, 320, 60, -66, 30], [-161, 0, -69, -418, -12, -130],
                 [320, -69, 0, 16, 198, 130], [60, -418, 16, 0, -147, -72],
                 [-66, -12, 198, -147, 0, -9], [30, -130, 130, -72, -9, 0]],
        (6, 5): [[0, -675, -264, -350, -504, -168], [-675, 0, -144, -565, -74, -118],
                 [-264, -144, 0, 8, 43, 38], [-350, -565, 8, 0, -240, -376],
                 [-504, -74, 43, -240, 0, -129], [-168, -118, 38, -376, -129, 0]],
        (6, 6): [[0, -98, -196, 58, 78, 48], [-98, 0, 126, -11, 473, -100],
                 [-196, 126, 0, 159, 74, 136], [58, -11, 159, 0, 601, 156],
                 [78, 473, 74, 601, 0, 186], [48, -100, 136, 156, 186, 0]],
    }
    return Matrix(data[(d, r)])


def congruence(M, scales):
    D = Matrix.diag(*scales)
    return D * M * D


# Fidelity against the paper's generic and rank-three formulas.
generic = {3: 1, 4: 4, 5: 11, 6: 26, 7: 57}
rank3 = {3: 1, 4: 3, 5: 6, 6: 10, 7: 15}
for d, expected in generic.items():
    got = parity_sum(all_ones(d))
    assert got == expected, ("generic", d, got, expected)
for d, expected in rank3.items():
    got = parity_sum(rank_three(d))
    assert got == expected, ("rank-three", d, got, expected)

# The d=4 generic/non-generic values in Example 2 are 4 and 3.
assert parity_sum(all_ones(4)) == 4
assert parity_sum(rank_three(4)) == 3

# Check the table formulas in Section 6 on exact matrices in every rank.
for d, r, offset in [(5, 4, 1), (6, 4, 2), (6, 5, 3), (6, 6, 4)]:
    base = fixed_rank_matrix(d, r)
    for j in range(8):
        M = congruence(base, [j + offset + i for i in range(d)])
        k = nonzero_principal_count(M, 4)
        l = nonzero_principal_count(M, 5) if d == 6 else 0
        rr = rank_three_submatrix_count(M, 5) if d == 6 else 0
        if d == 5:
            expected = 5 + k
        elif r == 4:
            expected = 5 + k + rr
        elif r == 5:
            expected = 4 + k + l + rr
        else:
            expected = 5 + k + l + rr
        assert parity_sum(M) == expected, (d, r, k, l, rr, parity_sum(M), expected)

print("generic parity sums d=3..7:", [parity_sum(all_ones(d)) for d in range(3, 8)])
print("rank-three parity sums d=3..7:", [parity_sum(rank_three(d)) for d in range(3, 8)])
print("Section 6 exact table formulas: 32 rank-class instances checked (d=5,6).")
print("ALL_CHECKS_OK")
