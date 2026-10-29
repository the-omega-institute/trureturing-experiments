#!/usr/bin/env python3
"""Exact checks for Hoft's OEIS A237591 tail-band conjecture."""

from collections import Counter


def ceil_div(a: int, b: int) -> int:
    assert b > 0
    return -((-a) // b)


def count_consecutive(n: int, j: int) -> int:
    """Count positive first parts for j consecutive positive parts."""
    h = j * (j + 1) // 2
    return 0 if n < h else (n - h) // j + 1


def displayed_formula(n: int, k: int) -> int:
    a = ceil_div(2 * (n + 1) - k * (k + 1), 2 * k)
    kp = k + 1
    b = ceil_div(2 * (n + 1) - kp * (kp + 1), 2 * kp)
    return a - b


def direct_entry(n: int, k: int) -> int:
    return count_consecutive(n, k) - count_consecutive(n, k + 1)


def row(n: int) -> list[int]:
    values = []
    k = 1
    while k * (k + 1) <= 2 * n:
        values.append(displayed_formula(n, k))
        k += 1
    return values


def main() -> None:
    # Fidelity: these are the first rows printed in the OEIS entry.
    expected_rows = {
        1: [1],
        2: [2],
        3: [2, 1],
        4: [3, 1],
        5: [3, 2],
        6: [4, 1, 1],
        7: [4, 2, 1],
        8: [5, 2, 1],
        9: [5, 2, 2],
        10: [6, 2, 1, 1],
    }
    for n, expected in expected_rows.items():
        assert row(n) == expected, (n, row(n), expected)
    print("fidelity: OEIS rows 1-10 reproduced exactly")

    # The native count and the displayed ceiling formula agree independently.
    compared = 0
    for n in range(1, 151):
        k = 1
        while k * (k + 1) <= 2 * n:
            assert displayed_formula(n, k) == direct_entry(n, k), (n, k)
            compared += 1
            k += 1
    print(f"count bridge: {compared} entries agree")

    # Exact integer validation of the full claimed band for k <= 200.
    counts = Counter()
    total = 0
    for k in range(1, 201):
        K = k * (k + 1)
        for n in range((K + 1) // 2, K + 1):
            value = displayed_formula(n, k)
            assert K <= 2 * n and n <= K
            assert value in (1, 2), (n, k, value)
            counts[value] += 1
            total += 1
    assert total == 1_353_600
    assert counts == Counter({1: 1_015_150, 2: 338_450})
    print(f"tail band: {total} pairs, counts {dict(sorted(counts.items()))}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
