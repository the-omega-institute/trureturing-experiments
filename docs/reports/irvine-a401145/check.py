#!/usr/bin/env python3
"""Exact finite check for OEIS A401145 and its local-maximum conjecture."""

from __future__ import annotations

from math import prod

from sympy import isprime


def distinct_partitions(n: int) -> list[tuple[int, ...]]:
    """Enumerate each unordered distinct positive partition exactly once."""
    out: list[tuple[int, ...]] = []

    def extend(remaining: int, largest: int, chosen: tuple[int, ...]) -> None:
        if remaining == 0:
            out.append(chosen)
            return
        for part in range(min(remaining, largest), 0, -1):
            extend(remaining - part, part - 1, chosen + (part,))

    extend(n, n, ())
    return out


def a(n: int) -> int:
    return sum(isprime(n + prod(parts)) for parts in distinct_partitions(n))


expected_partitions = {
    7: ((7,), (6, 1), (5, 2), (4, 3), (4, 2, 1)),
    8: ((8,), (7, 1), (6, 2), (5, 3), (5, 2, 1), (4, 3, 1)),
    9: ((9,), (8, 1), (7, 2), (6, 3), (6, 2, 1), (5, 4), (5, 3, 1), (4, 3, 2)),
}

for n, expected in expected_partitions.items():
    actual = tuple(distinct_partitions(n))
    assert actual == expected, (n, actual, expected)

counts = {n: a(n) for n in (7, 8, 9)}
assert counts == {7: 3, 8: 1, 9: 3}, counts
assert counts[8] <= max(counts[7], counts[9])
assert isprime(8 + prod((5, 3)))

print("A401145 exact counts:", counts)
print("n=8 witness: partition (5, 3), 8 + 5*3 = 23 is prime")
print("Conjectured inequality at n=8:", counts[8], "> max(", counts[7], ",", counts[9], ")")
print("ALL_CHECKS_OK")
