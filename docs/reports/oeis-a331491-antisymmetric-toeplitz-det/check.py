#!/usr/bin/env python3
"""Exact finite checks for the determinant conjecture in OEIS A331491."""

from fractions import Fraction


def source_matrix(order):
    """Zero-based translation of the displayed skew matrix and first row."""
    return [
        [0 if i == j else (j - i + 1) // 2 if i < j else -((i - j + 1) // 2)
         for j in range(order)]
        for i in range(order)
    ]


def pari_matrix(order):
    """Independent translation of the entry's one-based PARI tm(n)."""
    matrix = [[0] * order for _ in range(order)]
    for j in range(order):
        matrix[0][j] = (j + 1) // 2
    for i in range(1, order):
        matrix[i][0] = -((i + 1) // 2)
        for j in range(1, order):
            matrix[i][j] = matrix[i - 1][j - 1]
    return matrix


def identity(order):
    return [[int(i == j) for j in range(order)] for i in range(order)]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def multiply(left, right):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*right)]
            for row in left]


def congruence(matrix, basis):
    return multiply(multiply(transpose(basis), matrix), basis)


def determinant(matrix):
    """Bareiss elimination, with checked exact integer divisions and pivoting."""
    order = len(matrix)
    if not order:
        return 1
    work = [row[:] for row in matrix]
    denominator = 1
    sign = 1
    for k in range(order - 1):
        pivot_row = next((i for i in range(k, order) if work[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            work[k], work[pivot_row] = work[pivot_row], work[k]
            sign = -sign
        pivot = work[k][k]
        for i in range(k + 1, order):
            for j in range(k + 1, order):
                numerator = pivot * work[i][j] - work[i][k] * work[k][j]
                quotient, remainder = divmod(numerator, denominator)
                assert remainder == 0
                work[i][j] = quotient
            work[i][k] = 0
        denominator = pivot
    return sign * work[-1][-1]


def rational_determinant(matrix):
    """Independent ordinary Gaussian elimination over exact fractions."""
    work = [[Fraction(x) for x in row] for row in matrix]
    answer = Fraction(1)
    for k in range(len(work)):
        pivot_row = next((i for i in range(k, len(work)) if work[i][k]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != k:
            work[k], work[pivot_row] = work[pivot_row], work[k]
            answer = -answer
        pivot = work[k][k]
        answer *= pivot
        for i in range(k + 1, len(work)):
            factor = work[i][k] / pivot
            for j in range(k + 1, len(work)):
                work[i][j] -= factor * work[k][j]
            work[i][k] = Fraction(0)
    return answer


def permanent(matrix):
    """Subset recurrence sums signed integer products, without determinant signs."""
    order = len(matrix)
    partial = {0: 1}
    for i in range(order):
        updated = {}
        for mask, total in partial.items():
            for j in range(order):
                bit = 1 << j
                if not mask & bit:
                    updated[mask | bit] = updated.get(mask | bit, 0) + total * matrix[i][j]
        partial = updated
    return partial[(1 << order) - 1]


def constant_skew(order):
    return [[int(i < j) - int(j < i) for j in range(order)] for i in range(order)]


def check_reductions(n):
    order = 2 * n
    matrix = source_matrix(order)
    difference = identity(order)
    for i in range(1, order):
        difference[i][i - 1] = -1
    signs = [-1] + [(-1) ** (i + 1) for i in range(1, order)]
    reduced = multiply(multiply(difference, matrix), transpose(difference))
    reduced = [[-signs[i] * reduced[i][j] * signs[j] for j in range(order)]
               for i in range(order)]
    expected = constant_skew(order)
    for j in range(1, order):
        expected[0][j] = int(j % 2 == 1)
        expected[j][0] = -expected[0][j]
    assert reduced == expected
    assert determinant(difference) == 1
    assert determinant(reduced) == determinant(matrix) == n * n

    # Replace the final basis vector by e_2 - e_3 + ... + e_(2n).
    border_basis = identity(order)
    for i in range(1, order):
        border_basis[i][-1] = (-1) ** (i + 1)
    isolated = congruence(reduced, border_basis)
    assert determinant(border_basis) == 1
    assert isolated[0][-1] == n and isolated[-1][0] == -n
    assert all(isolated[i][-1] == isolated[-1][i] == 0 for i in range(1, order - 1))
    middle = [row[1:-1] for row in isolated[1:-1]]
    assert middle == constant_skew(order - 2)
    assert determinant(isolated) == n * n * determinant(middle)

    # The odd constant skew block annihilates its alternating vector.
    odd = constant_skew(order - 1)
    alternating = [(-1) ** i for i in range(order - 1)]
    assert all(sum(a * b for a, b in zip(row, alternating)) == 0 for row in odd)

    # Split J_(2n) into J_2 and J_(2n-2) by the stated integral basis change.
    basis = identity(order)
    for k in range(2, order):
        basis[0][k] = 1
        basis[1][k] = -1
    split = congruence(constant_skew(order), basis)
    assert determinant(basis) == 1
    assert [row[:2] for row in split[:2]] == constant_skew(2)
    assert [row[2:] for row in split[2:]] == constant_skew(order - 2)
    assert all(split[i][j] == split[j][i] == 0 for i in range(2) for j in range(2, order))
    assert determinant(constant_skew(order)) == 1


def main():
    displayed = [[0, 1, 1, 2], [-1, 0, 1, 1], [-1, -1, 0, 1], [-2, -1, -1, 0]]
    assert source_matrix(4) == displayed == pari_matrix(4)
    assert permanent(displayed) == 8
    print("FIDELITY: displayed M(4) and permanent a(2)=8 reproduced")

    expected_permanents = [1, -1, 8, -965, 301864]
    actual_permanents = [permanent(source_matrix(2 * n)) for n in range(5)]
    assert actual_permanents == expected_permanents
    print("FIDELITY: permanent data n=0..4:", actual_permanents)

    values = []
    for n in range(33):
        matrix = source_matrix(2 * n)
        if n:
            assert matrix == pari_matrix(2 * n)
        value = determinant(matrix)
        assert value == (n * n if n else 1)
        if n <= 8:
            assert rational_determinant(matrix) == value
        values.append(value)
    print("DETERMINANTS n=0..32:", values)
    print("CROSS_CHECK: fraction Gaussian elimination agrees for n=0..8")

    assert determinant(constant_skew(0)) == 1
    for n in range(1, 17):
        check_reductions(n)
    print("REDUCTIONS n=1..16: difference/sign, odd kernel, border isolation, J split all exact")
    print("BOUNDARY: empty determinant is 1; n=1 border uses J_0")
    print("Finite checks corroborate the outline; no all-size theorem or Lean proof is certified.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
