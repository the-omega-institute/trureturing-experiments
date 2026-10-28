#!/usr/bin/env python3
"""Exact small-order check for the Aztec-diamond local-correlation contrast.

The Kasteleyn determinants are evaluated over Gaussian integers by SymPy.  The
probability of a horizontal pair is the absolute determinant of the matrix with
its four occupied cells removed, divided by the absolute full determinant.
"""
from fractions import Fraction
from math import prod
import sympy as sp

I = sp.I
ORBIT_CENTERS = ((0, 0), (1, 1), (2, 0), (2, 2), (1, 0), (2, 1), (3, 0))
MULTIPLICITIES = (1, 4, 4, 4, 4, 8, 4)
EXPECTED_DELTA = {
    4: Fraction(-37, 576),
    5: Fraction(173, 2048),
    6: Fraction(-1191, 26624),
    7: Fraction(14445, 425984),
}
EXPECTED_MINOR_ABS = {
    4: (72, 120, 56, None, 160, 192, 64),
    5: (6400, 3840, 5888, 6144, 3072, 3200, 960),
    6: (262144, 163840, 235520, 225280, 327680, 270336, 356352),
    7: (33554432, 46137344, 32243712, 28835840, 28835840, 23920640, 29196288),
}
EXPECTED_FULL_ABS = {4: 1024, 5: 32768, 6: 2097152, 7: 268435456}


def aztec_cells(n):
    return [
        (i, j)
        for i in range(-n, n)
        for j in range(-n, n)
        if abs(2 * i + 1) + abs(2 * j + 1) <= 2 * n
    ]


def kasteleyn_matrix(n, removed=frozenset()):
    cells = [c for c in aztec_cells(n) if c not in removed]
    black = [c for c in cells if (c[0] + c[1]) % 2 == 0]
    white = [c for c in cells if (c[0] + c[1]) % 2 != 0]
    white_index = {c: k for k, c in enumerate(white)}
    matrix = sp.zeros(len(black), len(white))
    for row, (i, j) in enumerate(black):
        for di, dj, weight in ((0, 1, 1), (0, -1, -1), (1, 0, I), (-1, 0, -I)):
            w = (i + di, j + dj)
            if w in white_index:
                matrix[row, white_index[w]] = weight
    return matrix


def horizontal_pair_cells(x, y):
    # Cells occupy [j,j+1] x [i,i+1], so the square centered at (x,y)
    # has row/column cells (y-1,x-1), (y-1,x), (y,x-1), (y,x).
    return frozenset(((y - 1, x - 1), (y - 1, x), (y, x - 1), (y, x)))


def determinant_abs(matrix):
    value = sp.simplify(matrix.det(method="domain-ge"))
    assert value.is_real and value.is_integer, value
    return abs(int(value))


def contrast(n):
    full = determinant_abs(kasteleyn_matrix(n))
    values = []
    for center in ORBIT_CENTERS:
        x, y = center
        if abs(x) + abs(y) > n - 1:
            values.append(None)
        else:
            minor = determinant_abs(kasteleyn_matrix(n, horizontal_pair_cells(x, y)))
            values.append(Fraction(minor, full))
    black_indices = range(4) if values[3] is not None else range(3)
    black_count = sum(MULTIPLICITIES[k] for k in black_indices)
    black_mean = sum((MULTIPLICITIES[k] * values[k] for k in black_indices), Fraction()) / black_count
    white_count = sum(MULTIPLICITIES[4:])
    white_mean = sum((MULTIPLICITIES[k] * values[k] for k in range(4, 7)), Fraction()) / white_count
    return full, values, black_mean - white_mean


def main():
    for n in range(4, 8):
        full, values, delta = contrast(n)
        assert full == EXPECTED_FULL_ABS[n], (n, full)
        for actual, expected in zip(values, EXPECTED_MINOR_ABS[n]):
            assert actual is None if expected is None else actual * full == expected, (n, actual, expected)
        assert delta == EXPECTED_DELTA[n], (n, delta)
        assert ((-1) ** (n + 1)) * delta > 0, (n, delta)
        print(f"n={n}: Delta={delta}, sign={((-1) ** (n + 1)) * delta > 0}")
    assert Fraction(9, 128) == Fraction(EXPECTED_MINOR_ABS[4][0], EXPECTED_FULL_ABS[4])
    print("central P_H(4)=9/128")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
