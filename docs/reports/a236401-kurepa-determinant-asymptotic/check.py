#!/usr/bin/env python3
"""Exact checks for OEIS A236401 and the proposed determinant identity."""

from fractions import Fraction
from math import factorial

import sympy as sp


OEIS_VALUES = [
    15,
    -47,
    197,
    -1029,
    6439,
    -46927,
    390249,
    -3645737,
    37792331,
    -430400211,
    5341017373,
    -71724018781,
    1036207363983,
    -16024176975479,
    264083895859409,
    -4620276321889617,
    85520275455047059,
    -1669635965205539227,
    34287733935303686661,
]


def matrix(n: int) -> sp.Matrix:
    """The (n-4)-by-(n-4) matrix in Mestrovic's Definition 1, Eq. (4)."""
    m = n - 4
    rows = []
    for i in range(m):
        row = []
        for j in range(m):
            if i == m - 1:
                value = 1 if j == m - 2 else -4 if j == m - 1 else 0
            elif j == m - 1:
                value = 3 if i == 0 else 2
            elif j >= i or i == j + 2:
                value = 1
            elif i == j + 1:
                value = i + 2
            else:
                value = 0
            row.append(value)
        rows.append(row)
    return sp.Matrix(rows)


def kurepa(n: int) -> int:
    return int(matrix(n).det())


def derangements(limit: int) -> list[int]:
    d = [1, 0]
    for j in range(limit - 1):
        d.append((j + 1) * (d[j + 1] + d[j]))
    return d[: limit + 1]


def main() -> None:
    values = [kurepa(n) for n in range(7, 26)]
    assert values == OEIS_VALUES, (values, OEIS_VALUES)

    d = derangements(40)
    for j, value in enumerate(d):
        exact_partial_sum = sum((Fraction((-1) ** r, factorial(r)) for r in range(j + 1)), Fraction())
        assert value == factorial(j) * exact_partial_sum

    for n in range(7, 41):
        claimed = (-1) ** (n + 1) * (6 * factorial(n - 5) + d[n - 4] + d[n - 5])
        assert kurepa(n) == claimed, (n, kurepa(n), claimed)

    print("OEIS values n=7..25: exact determinant match")
    print("Derangement recurrence and finite alternating-sum formula: exact match through n=40")
    print("K_n = (-1)^(n+1) * (6*(n-5)! + d_(n-4) + d_(n-5)): exact match for 7<=n<=40")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
