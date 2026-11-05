#!/usr/bin/env python3
"""Exact checks of the literal OEIS A300177 array count and the proposed route."""

from collections import Counter


SOURCE_TERMS = (
    4, 32, 255, 2033, 16208, 129217, 1030173, 8212978, 65477359,
    522013397, 4161713160, 33178950053, 264516722873, 2108839989442,
    16812570686523, 134036975068993, 1068599860225264,
    8519333271178937, 67919753770243365,
)
M = ((6, 6, 6), (2, 1, 1), (0, 1, 0))
IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def matrix_product(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def matrix_vector(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3))
                 for i in range(3))


def original_neighbour_masks(n):
    """Construct every cell's neighbours, directly from king distance."""
    cells = [(r, c) for r in range(n) for c in range(3)]
    return tuple(sum(1 << j for j, q in enumerate(cells)
                     if p != q and abs(p[0] - q[0]) <= 1
                     and abs(p[1] - q[1]) <= 1) for p in cells)


def literal_valid(bits, masks):
    """Each cell must have at most seven equal-valued neighbours."""
    for i, mask in enumerate(masks):
        matching_bits = bits if (bits >> i) & 1 else ~bits
        if (mask & matching_bits).bit_count() > 7:
            return False
    return True


def row_valid(rows):
    return all(not (rows[i] == rows[i + 1] == rows[i + 2]
                    and rows[i] in (0, 7)) for i in range(len(rows) - 2))


def state(rows):
    if rows[-1] not in (0, 7):
        return 0
    if len(rows) >= 2 and rows[-2] == rows[-1]:
        return 2
    return 1


def literal_enumeration(n):
    masks = original_neighbour_masks(n)
    assert all(mask.bit_count() <= 8 for mask in masks)
    counts = [0, 0, 0]
    # Bit zero is B(0,0), so even bit patterns impose the source's anchor.
    for bits in range(0, 1 << (3 * n), 2):
        rows = tuple((bits >> (3 * r)) & 7 for r in range(n))
        valid = literal_valid(bits, masks)
        assert valid == row_valid(rows), (n, bits)
        if valid:
            counts[state(rows)] += 1
    return tuple(counts)


def two_row_counts(limit):
    """Independent DP retaining all row symbols, without the three-state quotient."""
    counts = {1: 4}
    pairs = Counter((None, row) for row in (0, 2, 4, 6))
    for n in range(2, limit + 1):
        new = Counter()
        for (previous, last), multiplicity in pairs.items():
            for row in range(8):
                if not (previous == last == row and row in (0, 7)):
                    new[last, row] += multiplicity
        pairs = new
        counts[n] = sum(pairs.values())
    return counts


def check():
    square = matrix_product(M, M)
    cube = matrix_product(square, M)
    assert square == ((48, 48, 42), (14, 14, 13), (2, 1, 1))
    assert cube == ((384, 378, 336), (112, 111, 98), (14, 14, 13))
    rhs = tuple(tuple(7 * square[i][j] + 7 * M[i][j] + 6 * IDENTITY[i][j]
                      for j in range(3)) for i in range(3))
    assert cube == rhs
    print("Exact matrix identity: M^3 = 7 M^2 + 7 M + 6 I")

    counts = {0: 0}
    vector = (3, 1, 0)
    for n in range(1, 211):
        counts[n] = sum(vector)
        if n <= 6:
            direct = literal_enumeration(n)
            assert direct == vector, (n, direct, vector)
            print(f"Literal enumeration n={n}: states={direct}, a(n)={sum(direct)}")
        vector = matrix_vector(M, vector)

    assert tuple(counts[n] for n in range(1, 20)) == SOURCE_TERMS
    print("Fidelity: all 19 displayed OEIS A300177 terms reproduced")
    pair_counts = two_row_counts(210)
    assert all(counts[n] == pair_counts[n] for n in range(1, 211))
    print("Full two-row-symbol DP equals three-state counts for 1 <= n <= 210")
    assert all(counts[n] == 7 * counts[n - 1] + 7 * counts[n - 2]
               + 6 * counts[n - 3] for n in range(4, 211))
    print("Recurrence checked for 4 <= n <= 210")
    denominator = (1, -7, -7, -6)
    numerator = (0, 4, 4, 3)
    for n in range(211):
        coefficient = sum(denominator[k] * counts[n - k]
                          for k in range(min(n, 3) + 1))
        assert coefficient == (numerator[n] if n < 4 else 0)
    print("Exact generating-function coefficient identity checked through x^210")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    check()
