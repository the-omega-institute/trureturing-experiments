#!/usr/bin/env python3
"""Exact checks for OEIS A183683 (five-column nonsingular binary strips)."""

from itertools import product


ROWS = list(product((0, 1), repeat=5))


def compatible(r, s):
    return all(r[j] * s[j + 1] - r[j + 1] * s[j] != 0 for j in range(4))


NEIGHBORS = {r: tuple(s for s in ROWS if compatible(r, s)) for r in ROWS}


def counts(length, edges=NEIGHBORS):
    """Number of labelled arrays with length+1 rows, by exact row DP."""
    values = {r: 1 for r in ROWS}
    for _ in range(length):
        values = {r: sum(values[s] for s in edges[r]) for r in ROWS}
    return sum(values.values())


def reduced_counts(length):
    q = ((0, 2, 1, 1), (2, 2, 0, 0), (1, 0, 0, 1), (4, 0, 4, 0))
    classes = (
        ("01010", "01011", "11010", "11011"),
        ("01101", "10110", "10111", "11101"),
        ("01110", "01111", "11110", "11111"),
        ("10101",),
    )
    cls = {tuple(map(int, row)): i for i, group in enumerate(classes) for row in group}
    v = [1, 1, 1, 1]
    for _ in range(length):
        v = [sum(q[i][j] * v[j] for j in range(4)) for i in range(4)]
    return 4 * v[0] + 4 * v[1] + 4 * v[2] + v[3], v, cls


def main():
    expected = [48, 208, 832, 3472, 14144, 58448, 239424, 986128, 4047680, 16650448]
    actual = [counts(n) for n in range(1, 11)]
    assert actual == expected, (actual, expected)

    # The 19 rows with a legal neighbor are exactly the four stated classes.
    nonisolated = {r for r in ROWS if NEIGHBORS[r]}
    class_rows = {
        tuple(map(int, row))
        for group in (("01010", "01011", "11010", "11011"),
                      ("01101", "10110", "10111", "11101"),
                      ("01110", "01111", "11110", "11111"),
                      ("10101",))
        for row in group
    }
    assert nonisolated == class_rows

    # Check the quotient matrix and the all-length recurrence through 200 terms.
    q_expected = ((0, 2, 1, 1), (2, 2, 0, 0), (1, 0, 0, 1), (4, 0, 4, 0))
    class_groups = (
        ("01010", "01011", "11010", "11011"),
        ("01101", "10110", "10111", "11101"),
        ("01110", "01111", "11110", "11111"),
        ("10101",),
    )
    for i, group in enumerate(class_groups):
        for row in group:
            row_tuple = tuple(map(int, row))
            got = tuple(sum(1 for s in NEIGHBORS[row_tuple]
                            if any(tuple(map(int, x)) == s for x in class_groups[j]))
                         for j in range(4))
            assert got == q_expected[i], (row, got)

    terms = [counts(n) for n in range(1, 201)]
    for n in range(4, 201):
        assert terms[n - 1] == 13 * terms[n - 3] + 16 * terms[n - 4]

    # Formal-series coefficient check: denominator times F has only the numerator.
    for n in range(1, 11):
        coefficient = terms[n - 1]
        if n >= 3:
            coefficient -= 13 * terms[n - 3]
        if n >= 4:
            coefficient -= 16 * terms[n - 4]
        expected_coefficient = {1: 48, 2: 208, 3: 208}.get(n, 0)
        assert coefficient == expected_coefficient, (n, coefficient)

    reduced = [reduced_counts(n)[0] for n in range(1, 11)]
    assert reduced == expected
    print("A183683 exact first terms:", ",".join(map(str, actual)))
    print("A183683 exact recurrence and quotient checks: 200 terms")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
