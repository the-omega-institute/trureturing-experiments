#!/usr/bin/env python3
"""Exact finite checks for the Shevelev--Moses OEIS A290540 conjecture.

Only Python integers and fractions are used. These checks do not establish
an asymptotic statement or a Lean theorem.
"""

from fractions import Fraction
from math import comb


def signed_residues(order, exponent):
    """Source K_(j+1)(exponent, order), for j=0,...,order-1."""
    return [sum((-1) ** a * comb(exponent, order * a + j)
                for a in range(exponent // order + 1)
                if order * a + j <= exponent)
            for j in range(order)]


def right_circulant(values):
    order = len(values)
    return [[values[(j - i) % order] for j in range(order)]
            for i in range(order)]


def bareiss(matrix):
    """Integer determinant, including singular matrices and row swaps."""
    a = [row[:] for row in matrix]
    order = len(a)
    previous = 1
    parity = 1
    for k in range(order - 1):
        pivot_row = next((i for i in range(k, order) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            parity = -parity
        pivot = a[k][k]
        for i in range(k + 1, order):
            for j in range(k + 1, order):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0, (order, k, i, j)
                a[i][j] = quotient
            a[i][k] = 0
        previous = pivot
    return parity * a[-1][-1]


def gaussian_fraction(matrix):
    """Independent rational Gaussian elimination for determinant cross-checks."""
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for k in range(len(a)):
        pivot_row = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            result = -result
        pivot = a[k][k]
        result *= pivot
        for i in range(k + 1, len(a)):
            ratio = a[i][k] / pivot
            for j in range(k, len(a)):
                a[i][j] -= ratio * a[k][j]
    assert result.denominator == 1
    return result.numerator


def determinant(order, exponent):
    k = signed_residues(order, exponent)
    return bareiss(right_circulant([(-1) ** j * x for j, x in enumerate(k)]))


def sign(value):
    return (value > 0) - (value < 0)


def phase_sign(order, exponent):
    """Exact sign of sin((exponent+1)*pi/order), without evaluating sine."""
    residue = (exponent + 1) % (2 * order)
    if residue in (0, order):
        return 0
    return 1 if residue < order else -1


def main():
    # The source recurrence is independent of the defining binomial sum.
    recurrence_count = 0
    for order in range(2, 21):
        state = [1] + [0] * (order - 1)
        for exponent in range(8 * order + 1):
            assert state == signed_residues(order, exponent), (order, exponent)
            recurrence_count += 1
            state = [state[0] - state[-1]] + [
                state[j] + state[j - 1] for j in range(1, order)]
    print(f"Source Definition 4 / equation (14): {recurrence_count} exact recurrence/binomial comparisons OK")

    # Fidelity: arXiv:1706.01454v4, Theorem 3(2).
    odd_count = 0
    for order in range(3, 14, 2):
        for exponent in range(1, 6 * order + 1):
            assert determinant(order, exponent) == 0, (order, exponent)
            odd_count += 1
    print(f"Source Theorem 3(2): {odd_count} odd-order determinants vanish (N=3,5,7,9,11,13; 1<=m<=6N)")

    # Independently transcribed data in OEIS A290540 revision #36.
    oeis = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -2276485387658524, -523547340003805770400,
            -39617190432735671861429500,
            -2896792542975174202888623380000,
            -95819032881785191861991031568287500,
            -1018409199709889673458815786392849200000]
    actual = [determinant(10, exponent) for exponent in range(len(oeis))]
    assert actual == oeis, (actual, oeis)
    print("OEIS A290540 revision #36: all 16 displayed terms reproduced exactly")

    even_count = 0
    initial_count = 0
    block_count = 0
    rational_count = 0
    for order in range(2, 21, 2):
        values = [determinant(order, exponent)
                  for exponent in range(8 * order + 1)]
        for exponent, value in enumerate(values):
            expected = phase_sign(order, exponent)
            if 1 <= exponent <= order - 2:
                assert value == 0
                assert expected == 1
                initial_count += 1
            else:
                assert sign(value) == expected, (order, exponent, sign(value), expected)
            even_count += 1
        for u in range(4):
            start = 2 * order * u
            if u > 0:
                assert all(values[start + r] > 0 for r in range(order - 1))
                block_count += 1
            assert values[start + order - 1] == 0
            assert all(values[start + order + r] < 0 for r in range(order - 1))
            assert values[start + 2 * order - 1] == 0
            block_count += 1
        for v in range(8):
            exponent = order * v + order - 1
            k = signed_residues(order, exponent)
            assert all(k[j] == (-1) ** v * k[order - 1 - j]
                       for j in range(order))
        for exponent in [0, 1, order - 1, order, order + 1, 2 * order]:
            k = signed_residues(order, exponent)
            matrix = right_circulant([(-1) ** j * x for j, x in enumerate(k)])
            assert gaussian_fraction(matrix) == values[exponent]
            rational_count += 1
    print(f"Even orders N=2,4,...,20; 0<=m<=8N: {even_count} exact determinant signs OK")
    print(f"Precisely {initial_count} initial zeros at 1<=m<=N-2 are exceptions to the phase sign rule")
    print(f"Strict positive/negative blocks with separating zeros: {block_count} complete blocks OK")
    print("Boundary reflection K_j(Nv+N-1)=(-1)^v K_(N-j+1)(Nv+N-1): all N=2,4,...,20 and 0<=v<=7 OK")
    print(f"Independent Fraction Gaussian determinant cross-check: {rational_count} cases OK")
    print("All evidence above is finite exact arithmetic; no infinite claim or Lean proof is asserted")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
