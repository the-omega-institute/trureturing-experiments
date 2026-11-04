#!/usr/bin/env python3
"""Exact checks for OEIS A206982 (Hardin; Barker's recurrence)."""

from itertools import product


ROWS = tuple(product((0, 1), repeat=3))
LEGAL_ROWS = tuple(row for row in ROWS if row not in ((0, 1, 0), (1, 0, 1)))


def alternating(x, y, z):
    return x != y and y != z


def compatible(a, b, c):
    """The five local tests for three successive rows of width three."""
    row_tests = (alternating(a[j], b[j], c[j]) for j in range(3))
    diagonals = (
        alternating(a[0], b[1], c[2]),
        alternating(a[2], b[1], c[0]),
    )
    return not any((*row_tests, *diagonals))


def literal_count(n):
    """Enumerate all legal row words, including every local direction."""
    total = 0
    for rows in product(LEGAL_ROWS, repeat=n):
        if all(compatible(rows[i], rows[i + 1], rows[i + 2])
               for i in range(max(0, n - 2))):
            total += 1
    return total


PAIRS = tuple(product(LEGAL_ROWS, repeat=2))
SUCCESSORS = {
    pair: tuple((pair[1], c) for c in LEGAL_ROWS if compatible(pair[0], pair[1], c))
    for pair in PAIRS
}


def transfer_count(n):
    if n == 1:
        return len(LEGAL_ROWS)
    if n == 2:
        return len(PAIRS)
    state = {pair: 1 for pair in PAIRS}
    for _ in range(n - 2):
        nxt = {pair: 0 for pair in PAIRS}
        for pair, weight in state.items():
            for successor in SUCCESSORS[pair]:
                nxt[successor] += weight
        state = nxt
    return sum(state.values())


def main():
    literal = [literal_count(n) for n in range(1, 8)]
    expected = [6, 36, 58, 158, 420, 1066, 2754]
    assert literal == expected, (literal, expected)

    values = [None] + [transfer_count(n) for n in range(1, 41)]
    assert values[1:8] == expected
    for n in range(8, 41):
        assert values[n] == (
            2 * values[n - 1] + values[n - 2] + 2 * values[n - 3]
            - values[n - 4] - 2 * values[n - 5]
        ), n

    # Coefficient check for Barker's rational generating function.  All
    # arithmetic here is over integers; no floating point is used.
    denominator = [1, -2, -1, -2, 1, 2]
    numerator = [0, 6, 24, -20, -6, -20, 0, 16]
    for n in range(1, 41):
        coefficient = sum(denominator[j] * (values[n - j] if n - j >= 1 else 0)
                         for j in range(min(n, 5) + 1))
        assert coefficient == (numerator[n] if n < len(numerator) else 0), n

    print("literal a(1..7) =", literal)
    print("transfer a(1..10) =", values[1:11])
    print("recurrence verified for n = 8..40")
    print("generating-function numerator verified through degree 40")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
