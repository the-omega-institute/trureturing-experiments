#!/usr/bin/env python3
"""Exact finite checks for the OEIS A231147 subset-median interpretation."""
from __future__ import annotations

from itertools import combinations
from math import comb
from pathlib import Path


def subsets(n: int):
    """Yield every subset of {1,...,n} as a frozenset, excluding empty."""
    ground = tuple(range(1, n + 1))
    for size in range(1, n + 1):
        yield from (frozenset(s) for s in combinations(ground, size))


def median_twice(subset: frozenset[int]) -> int:
    values = sorted(subset)
    k = len(values)
    if k % 2:
        return 2 * values[k // 2]
    return values[k // 2 - 1] + values[k // 2]


def numerator_coefficients(n: int) -> list[int]:
    """Coefficients of sum_j x^(n-1-j)(1+x^2)^j, by exact integers."""
    coeffs = [0] * (2 * n - 1)
    for j in range(n):
        for q in range(j + 1):
            coeffs[n - 1 - j + 2 * q] += comb(j, q)
    return coeffs


def encode(n: int, subset: frozenset[int]) -> tuple[int, frozenset[int]]:
    values = sorted(subset)
    if len(values) % 2:
        a = b = values[len(values) // 2]
    else:
        a, b = values[len(values) // 2 - 1 : len(values) // 2 + 1]
    lower = {t for t in subset if t < a}
    upper = {t for t in subset if t > b}
    j = n - 1 - (b - a)
    q = {t for t in range(1, a) if t not in lower}
    q.update(a - 1 + t - b for t in upper)
    return j, frozenset(q)


def decode(n: int, j: int, q: frozenset[int]) -> frozenset[int]:
    a = len(q) + 1
    b = n + len(q) - j
    lower = {t for t in range(1, a) if t not in q}
    upper = {b + t - (a - 1) for t in q if t >= a}
    return frozenset(lower | {a, b} | upper)


def check_row(n: int) -> list[int]:
    counts = [0] * (2 * n - 1)
    for subset in subsets(n):
        r = median_twice(subset) - 2
        assert 0 <= r <= 2 * n - 2
        counts[r] += 1
    return counts


def check_bijection(n: int) -> int:
    seen: set[tuple[int, frozenset[int]]] = set()
    checked = 0
    for subset in subsets(n):
        j, q = encode(n, subset)
        assert 0 <= j < n
        assert q <= frozenset(range(1, j + 1))
        assert len(q) == sorted(subset)[(len(subset) - 1) // 2] - 1 if len(subset) % 2 else len(q) == sorted(subset)[len(subset) // 2 - 1] - 1
        assert decode(n, j, q) == subset
        assert (j, q) not in seen
        seen.add((j, q))
        checked += 1
    expected_pairs = sum(2**j for j in range(n))
    assert len(seen) == expected_pairs == (2**n - 1)
    for j in range(n):
        for mask in range(1 << j):
            q = frozenset(i + 1 for i in range(j) if mask & (1 << i))
            subset = decode(n, j, q)
            assert encode(n, subset) == (j, q)
    return checked


def main() -> None:
    expected_row_5 = [1, 1, 5, 4, 9, 4, 5, 1, 1]
    assert numerator_coefficients(5) == expected_row_5
    assert check_row(5) == expected_row_5
    print("n=5 coefficients:", expected_row_5)

    total_subsets = 0
    for n in range(1, 17):
        coeffs = numerator_coefficients(n)
        row = check_row(n)
        assert row == coeffs
        checked = check_bijection(n)
        total_subsets += checked
        print(f"n={n}: subsets={checked}, coefficients_match=True, bijection=True")
    assert total_subsets == 131054
    print("checked_nonempty_subsets=131054")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
