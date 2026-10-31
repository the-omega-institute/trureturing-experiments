#!/usr/bin/env python3
"""Exact finite checks of Wiseman's historical A083368 first-column conjecture.

No external packages, network calls, or floating-point operations are used.
This program tests a finite range; it does not prove the infinite conjecture.
"""

from itertools import groupby


MAX_BITS = 20
LIMIT = 1 << MAX_BITS


def profile(n):
    """Scan actual binary digits from the least significant end."""
    result = []
    run = 0
    while n:
        if n & 1:
            run += 1
        elif run:
            result.append(run)
            run = 0
        n >>= 1
    if run:
        result.append(run)
    return tuple(result)


def profile_by_indices(n):
    """Independent A048793-style binary-index split."""
    indices = [i + 1 for i in range(n.bit_length()) if (n >> i) & 1]
    return tuple(len(list(g)) for _, g in groupby(enumerate(indices), lambda t: t[1] - t[0]))


FIB = [0, 1]
for _ in range(MAX_BITS + 4):
    FIB.append(FIB[-1] + FIB[-2])


def fib_value(code):
    return sum(FIB[i + 2] for i in range(code.bit_length()) if (code >> i) & 1)


def minimum_representative(row):
    # row is little-endian in its run order. Exactly one zero separates runs.
    return int("0".join("1" * r for r in reversed(row)), 2)


def changed(u, v):
    # XOR pads missing high bits with zero. bit_length is a one-based position.
    return (u ^ v).bit_length()


def main():
    source_rows = [
        (), (1,), (1,), (2,), (1,), (1, 1), (2,), (3,), (1,), (1, 1),
        (1, 1), (2, 1), (2,), (1, 2), (3,), (4,), (1,), (1, 1),
        (1, 1), (2, 1), (1, 1), (1, 1, 1),
    ]
    assert [profile(n) for n in range(22)] == source_rows
    assert all(profile(n) == profile_by_indices(n) for n in range(4096))
    assert profile(53) == (1, 1, 2)
    print("FIDELITY A245563: all 22 printed rows n=0..21 reproduced")

    blocks = [[] for _ in range(MAX_BITS + 1)]
    blocks[0] = [0]
    first = {}
    for n in range(1, LIMIT):
        row = profile(n)
        if row not in first:
            first[row] = n
        if "00" not in format(n, "b"):
            blocks[n.bit_length()].append(n)

    lazy = [n for block in blocks for n in block]
    assert all(fib_value(code) == n for n, code in enumerate(lazy))
    assert lazy[27] == int("110111", 2)
    assert lazy[28] == int("111010", 2)
    assert changed(lazy[27], lazy[28]) == 4
    print("FIDELITY A083368: 27=110111, 28=111010 in lazy Fibonacci notation; a(28)=4")

    source_prefix = [
        1, 2, 1, 3, 2, 1, 4, 1, 3, 2, 1, 5, 2, 1, 4, 1, 3, 2, 1, 6,
        1, 3, 2, 1, 5, 2, 1, 4, 1, 3, 2, 1, 7, 2,
    ]
    a = [changed(lazy[n - 1], lazy[n]) for n in range(1, len(lazy))]
    assert a[:len(source_prefix)] == source_prefix
    print("FIDELITY A083368: all 34 terms on the source's first data line reproduced")

    # Dict insertion order implements first occurrence, never last occurrence.
    rows = list(first)
    codes = list(first.values())
    assert codes == sorted(codes)
    assert all(first[row] == minimum_representative(row) for row in rows)
    canonical = [n for block in blocks[1:] for n in block if n & 1]
    assert codes == canonical
    assert rows[15] == profile(53)  # A385817 row 16, after its empty row 0.
    assert [row[0] for row in rows] == a[:len(rows)]
    assert len(rows) == FIB[MAX_BITS + 2] - 1
    print(f"IDENTITY: every n=1..{len(rows)} agrees; scanned binary integers 0..{LIMIT - 1}")
    print(f"FIRST OCCURRENCES: {len(rows)}; every minimum is odd and avoids 00")
    print(f"LAZY REPRESENTATIONS: all values 0..{len(lazy) - 1} unique and consecutive")

    D = [[] for _ in blocks]
    T = [[] for _ in blocks]
    for m in range(1, MAX_BITS + 1):
        B = blocks[m]
        lo, hi = FIB[m + 2] - 1, FIB[m + 3] - 2
        assert [fib_value(n) for n in B] == list(range(lo, hi + 1))
        assert changed(blocks[m - 1][-1], B[0]) == m
        D[m] = [changed(u, v) for u, v in zip(B, B[1:])]
        T[m] = [profile(n)[0] for n in B if n & 1]
        assert T[m] == D[m - 1] + [m]
        if m >= 2:
            prefix10 = [(1 << (m - 1)) + n for n in blocks[m - 2]]
            prefix1 = [(1 << (m - 1)) + n for n in blocks[m - 1]]
            assert B == prefix10 + prefix1
            assert D[m] == D[m - 2] + [m - 1] + D[m - 1]
        if m >= 3:
            assert T[m] == T[m - 2] + T[m - 1][:-1] + [m]
    print(f"BLOCK CHECKS: ordered value intervals, boundaries, D and T recurrences, T_m=D_(m-1)||[m], m=1..{MAX_BITS}")
    print("TITLE CHECK: at n=3 the current NAME gives 2 but the historical change statistic gives 1")
    assert lazy[3].bit_length() == 2 and a[2] == 1
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
