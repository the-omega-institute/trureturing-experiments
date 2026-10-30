#!/usr/bin/env python3
"""Exact finite checks for Cuntz, arXiv:math/0611233v2, Section 7(5).

Standard library only. Enumeration concerns the larger class of square integer
matrices with first column one and nonzero mutually orthogonal rows. It is
exhaustive within the printed entry bound, up to row permutation. It is finite
experimental evidence, not a proof of the unrestricted assertion.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt


SOURCE_13 = (
    (1, 1, 1, 1, 2, 8, 8, 8, 10, 10, 20, 20, 20),
    (1, 1, 1, 1, 2, 8, 8, 8, 10, 10, 20, -20, -20),
    (1, 1, 1, 1, 2, 8, 8, 8, 10, 10, -20, 20, -20),
    (1, 1, 1, 1, 2, 8, 8, 8, 10, 10, -20, -20, 20),
    (1, 1, 1, 1, 2, 8, 8, 8, -10, -10, 0, 0, 0),
    (1, 1, 1, 1, 2, 3, -2, -2, 0, 0, 0, 0, 0),
    (1, 1, 1, 1, 2, -2, 3, -2, 0, 0, 0, 0, 0),
    (1, 1, 1, 1, 2, -2, -2, 3, 0, 0, 0, 0, 0),
    (1, 1, 1, 1, -2, 0, 0, 0, -2, 2, 0, 0, 0),
    (1, 1, 1, 1, -2, 0, 0, 0, 2, -2, 0, 0, 0),
    (1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
)


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def check_square(s):
    n = len(s)
    assert n and all(len(row) == n and row[0] == 1 for row in s)
    norms = [dot(row, row) for row in s]
    assert all(d > 0 for d in norms)
    for i in range(n):
        for j in range(n):
            assert dot(s[i], s[j]) == (norms[i] if i == j else 0)
    weights = [Q(1, d) for d in norms]
    assert sum(weights) == 1
    # Directly recompute every entry of s^T diag(1/d_i) s.
    for j in range(n):
        for k in range(n):
            assert sum(weights[i] * s[i][j] * s[i][k] for i in range(n)) == int(j == k)
    for k in range(n):
        assert sum(weights[i] * (s[i][k] ** 2 - 1) for i in range(n)) == 0
        zero_mass = sum(weights[i] for i in range(n) if s[i][k] == 0)
        excess = sum(weights[i] * (s[i][k] ** 2 - 1)
                     for i in range(n) if abs(s[i][k]) >= 2)
        assert zero_mass == excess
        if all(row[k] != 0 for row in s):
            assert all(abs(row[k]) == 1 for row in s)
    return norms


def fidelity():
    s = SOURCE_13
    n = len(s)
    norms = check_square(s)
    roots = [isqrt(d) for d in norms]
    assert all(r * r == d for r, d in zip(roots, norms))
    for i, j in product(range(n), repeat=2):
        assert roots[j] * s[i][j] == roots[i] * s[j][i]
    constants = [sum(Q(s[k][i] * s[k][j] * s[k][l], norms[k])
                     for k in range(n))
                 for i, j, l in product(range(n), repeat=3)]
    assert all(x.denominator == 1 for x in constants)
    w = roots[0]
    for i, j in product(range(n), repeat=2):
        assert (s[0][i] * s[i][j] - s[0][i] * s[0][j]) % w == 0
    # Source Section 4, 13-dimensional example: g_i*g_j need not be
    # divisible by mathbb S_ij = s_ij*g_i. Source indices i=j=6.
    i = j = 5
    scaled_entry = s[i][j] * s[0][i]
    degree_product = s[0][i] * s[0][j]
    assert scaled_entry == 24 and degree_product == 64
    assert degree_product % scaled_entry != 0
    print('FIDELITY: Section 4 dimension-13 example satisfies Fourier conditions (1)-(4)')
    print('FIDELITY: congruence conditions (5)-(6), w =', w)
    print('FIDELITY: row norms =', norms)
    print('FIDELITY: structure constants =', len(constants),
          '; negative =', sum(x < 0 for x in constants))
    print('FIDELITY: source indices (6,6): mathbb S_66 = 24 does not divide g_6*g_6 = 64')
    print('SOURCE_13: weighted Gram, zero-excess balance and zero-free-column tests pass')


def small_cases(n, bound):
    rows = [(1,) + tail for tail in product(range(-bound, bound + 1), repeat=n - 1)]
    later_orthogonal = [{j for j in range(i + 1, len(rows)) if dot(row, rows[j]) == 0}
                        for i, row in enumerate(rows)]
    count = 0
    def visit(chosen, candidates):
        nonlocal count
        if len(chosen) == n:
            check_square(tuple(rows[i] for i in chosen))
            count += 1
            return
        if len(candidates) < n - len(chosen):
            return
        for i in sorted(candidates):
            visit(chosen + (i,), candidates & later_orthogonal[i])
    visit((), set(range(len(rows))))
    print(f'EXHAUSTIVE: n={n}, entries in [-{bound},{bound}], first column one; '
          f'{count} matrices up to row permutation; all checks pass')
    return count


def boundaries():
    # These show why the unrestricted negative answer is not claimed for
    # rational entries, missing normalization, or rectangular matrices.
    rational = ((Q(1), Q(7, 5)), (Q(1), -Q(5, 7)))
    assert dot(*rational) == 0
    assert all(r[1] and abs(r[1]) != 1 for r in rational)
    unnormalized = ((1, 2), (2, -1))
    assert dot(*unnormalized) == 0
    assert all(r[0] != 0 for r in unnormalized) and unnormalized[1][0] == 2
    rectangular = ((1, 1, 2), (1, 1, -1))
    assert dot(*rectangular) == 0
    assert all(r[2] != 0 for r in rectangular) and rectangular[0][2] == 2
    print('BOUNDARIES: rational, unnormalized and rectangular counterexamples checked exactly')


def main():
    fidelity()
    total = sum(small_cases(n, 2) for n in range(1, 5))
    assert total > 0
    print('EXHAUSTIVE_TOTAL:', total)
    boundaries()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
