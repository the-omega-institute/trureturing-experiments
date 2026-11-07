#!/usr/bin/env python3
"""Exact checks for Michael Chu's A000384 balanced interval-cover conjecture."""

from itertools import combinations
from pathlib import Path


def hexagonal(n: int) -> int:
    return n * (2 * n - 1)


def construction(n: int):
    N = hexagonal(n)
    A = set(range(n)) | set(range(N - n + 1, N + 1))
    B = {j * n for j in range(2 * n)}
    return N, A, B


def sumset(A, B):
    return {a + b for a in A for b in B}


def bit_sum_mask(A, B):
    mask_b = sum(1 << b for b in B)
    mask = 0
    for a in A:
        mask |= mask_b << a
    return mask


def exhaustive_minimum(n: int) -> int:
    """Enumerate every equal-size pair for k < 2n and the first attaining k."""
    N = hexagonal(n)
    target = (1 << (2 * N + 1)) - 1
    for k in range(1, 2 * n + 1):
        subsets = list(combinations(range(N + 1), k))
        # Compute all pair sums exactly as bit masks.  This is exhaustive over A and B.
        for A in subsets:
            for B in subsets:
                if bit_sum_mask(A, B) == target:
                    return k
    raise AssertionError("no covering pair found")


def main() -> None:
    for n in range(1, 31):
        N, A, B = construction(n)
        assert A <= set(range(N + 1)) and B <= set(range(N + 1))
        assert len(A) == len(B) == 2 * n
        assert sumset(A, B) == set(range(2 * N + 1))
    print("construction n=1..30: OK")

    # Exhaustive equal-size search for n <= 3.  The largest search is
    # C(16,5)^2 = 19,079,424 pairs, checked with exact integer masks.
    for n in range(1, 4):
        found = exhaustive_minimum(n)
        assert found == 2 * n
        print(f"exhaustive minimum n={n}: k={found}")

    for n in range(1, 31):
        N = hexagonal(n)
        assert (2 * n - 1) ** 2 < 2 * N + 1
    print("counting lower bound n=1..30: OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
