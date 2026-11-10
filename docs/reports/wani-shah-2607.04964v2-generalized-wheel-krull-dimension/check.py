#!/usr/bin/env python3
"""Exact small-instance check for Wani--Shah generalized-wheel dimensions.

All quantities are integers.  A q-free face is a vertex subset whose induced
matching number is below q; the Stanley--Reisner dimension is its maximum size.
"""
from __future__ import annotations

from functools import lru_cache


def graph(n: int, m: int, independent_hubs: bool) -> tuple[int, tuple[int, ...]]:
    N = n + m
    adj = [0] * N
    def add(u: int, v: int) -> None:
        if u == v:
            return
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    for i in range(n):
        add(i, (i + 1) % n)
    for i in range(n):
        for j in range(n, N):
            add(i, j)
    if not independent_hubs:
        for i in range(n, N):
            for j in range(i + 1, N):
                add(i, j)
    return N, tuple(adj)


def matching_numbers(N: int, adj: tuple[int, ...]) -> list[int]:
    @lru_cache(maxsize=None)
    def mm(mask: int) -> int:
        if mask == 0:
            return 0
        low = mask & -mask
        u = low.bit_length() - 1
        best = mm(mask ^ low)
        nbrs = adj[u] & mask
        while nbrs:
            bit = nbrs & -nbrs
            best = max(best, 1 + mm(mask ^ low ^ bit))
            nbrs ^= bit
        return best
    return [mm(mask) for mask in range(1 << N)]


def exact_dimension(n: int, m: int, independent_hubs: bool, q: int) -> tuple[int, int]:
    N, adj = graph(n, m, independent_hubs)
    mm = matching_numbers(N, adj)
    nu = mm[(1 << N) - 1]
    dimension = max(mask.bit_count() for mask, value in enumerate(mm) if value < q)
    return nu, dimension


def predicted(n: int, m: int, independent_hubs: bool, q: int) -> int:
    h = n // 2
    hub_term = max(m, h) if independent_hubs else h
    return q - 1 + max(q, hub_term) if q <= (min(n, (n + m) // 2) if independent_hubs else (n + m) // 2) else n + m


def main() -> None:
    cases = 0
    for n in range(3, 9):
        for m in range(1, 5):
            for independent in (False, True):
                nu_plus, _ = exact_dimension(n, m, independent, 1)
                expected_nu = min(n, (n + m) // 2) if independent else (n + m) // 2
                assert nu_plus == expected_nu, (n, m, independent, nu_plus, expected_nu)
                for q in range(1, nu_plus + 2):
                    nu, dim = exact_dimension(n, m, independent, q)
                    assert nu == expected_nu
                    assert dim == predicted(n, m, independent, q), (
                        n, m, independent, q, dim, predicted(n, m, independent, q)
                    )
                    cases += 1
    # The source's q=1 dimensions are reproduced independently.
    for n in range(3, 9):
        for m in range(1, 5):
            _, plus_dim = exact_dimension(n, m, False, 1)
            _, minus_dim = exact_dimension(n, m, True, 1)
            assert plus_dim == n // 2
            assert minus_dim == max(n // 2, m)
    print(f"checked {cases} exact (n,m,q,family) cases; source q=1 formulas reproduced")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
