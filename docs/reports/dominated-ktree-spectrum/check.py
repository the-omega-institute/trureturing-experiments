#!/usr/bin/env python3
"""Exact finite checks of open-neighborhood dominated coloring of path powers.

Source: Sneha R, J. Geetha, K. Somasundaram, arXiv:2609.32439v1,
Section 2.3, path formula and Open Problems item 1.
Python standard library only; every evidence calculation uses integers.
"""
from functools import lru_cache


def path_power(k, n):
    return [sum(1 << j for j in range(n) if 0 < abs(i - j) <= k)
            for i in range(n)]


def exact_minimum(adj):
    """Enumerate ALL independent subsets with a common open-neighborhood witness.

    The subset recurrence computes both independence and the intersection of
    open neighborhoods. A partition DP then chooses the class containing the
    least remaining vertex; thus every unlabeled admissible partition is covered.
    """
    n = len(adj)
    full = (1 << n) - 1
    independent = [False] * (full + 1)
    independent[0] = True
    common = [full] * (full + 1)
    by_vertex = [[] for _ in range(n)]
    for mask in range(1, full + 1):
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        independent[mask] = independent[rest] and not (adj[v] & rest)
        common[mask] = common[rest] & adj[v]
        if independent[mask] and common[mask]:
            for w in range(n):
                if mask & (1 << w):
                    by_vertex[w].append(mask)
    choice = {}

    @lru_cache(None)
    def dp(mask):
        if mask == 0:
            return 0
        v = (mask & -mask).bit_length() - 1
        best = n + 1
        for cls in by_vertex[v]:
            if cls & mask == cls:
                cost = 1 + dp(mask ^ cls)
                if cost < best:
                    best, choice[mask] = cost, cls
        return best

    value = dp(full)
    assert value <= n, ('no admissible partition', adj)
    partition = []
    mask = full
    while mask:
        cls = choice[mask]
        partition.append(cls)
        mask ^= cls
    validate(adj, partition)
    assert len(partition) == value
    return value, partition


def validate(adj, partition):
    full = (1 << len(adj)) - 1
    used = 0
    for cls in partition:
        assert cls and not used & cls
        assert any(cls & neighborhood == cls for neighborhood in adj)
        assert all(not (cls & adj[v]) for v in range(len(adj)) if cls >> v & 1)
        used |= cls
    assert used == full


def partition_backtrack(adj):
    """Independent implementation: directly build restricted-growth partitions.

    No subset tables or DP are used. Each existing color is tried, followed by
    one new color; adjacency and common dominators are tested directly.
    """
    n = len(adj)
    best = n + 1
    classes = []

    def visit(v):
        nonlocal best
        if len(classes) >= best:
            return
        if v == n:
            best = len(classes)
            return
        for cls in classes:
            if all(not (adj[v] >> w & 1) for w in cls):
                proposed = cls + [v]
                if any(all(adj[d] >> w & 1 for w in proposed) for d in range(n)):
                    cls.append(v)
                    visit(v + 1)
                    cls.pop()
        if adj[v]:
            classes.append([v])
            visit(v + 1)
            classes.pop()

    visit(0)
    return best


def check_ktree(adj, k):
    """Actually delete last vertices, checking degrees and clique neighborhoods."""
    n = len(adj)
    remaining = set(range(n))
    for v in range(n - 1, k - 1, -1):
        neighbors = [w for w in remaining if adj[v] >> w & 1]
        assert len(neighbors) == k
        assert all(adj[a] >> b & 1 for a in neighbors for b in neighbors if a != b)
        remaining.remove(v)
    assert len(remaining) == k
    assert all(adj[a] >> b & 1 for a in remaining for b in remaining if a != b)
    reached = {0}
    while True:
        more = reached | {w for v in reached for w in range(n) if adj[v] >> w & 1}
        if more == reached:
            break
        reached = more
    assert len(reached) == n


def main():
    print('Exact arithmetic: Python integers only; open neighborhoods; used colors.')
    print('Fidelity: source Section 2.3 path formula, for n=2,...,16.')
    table = {}
    partitions = {}
    for k in range(1, 5):
        values = []
        for n in range(max(2, k + 1), 17):
            adj = path_power(k, n)
            assert all(not (adj[v] >> v & 1) for v in range(n))
            assert all((adj[v] >> w & 1) == (adj[w] >> v & 1)
                       for v in range(n) for w in range(n))
            check_ktree(adj, k)
            value, partition = exact_minimum(adj)
            table[k, n] = value
            partitions[k, n] = partition
            if n <= 8:
                assert partition_backtrack(adj) == value
            degree = max(x.bit_count() for x in adj)
            assert degree <= 2 * k
            assert n <= degree * value  # source maximum-degree lower bound
            assert n <= 2 * k * value  # proposed coarse capacity bound
            if k == 1:
                predicted = n // 2 + (0 if n % 4 == 0 else 1)
                assert value == predicted, (n, value, predicted)
            if n == k + 1:
                assert value == k + 1
            if (k, n - 1) in table:
                assert value <= table[k, n - 1] + 1
                validate(adj, partitions[k, n - 1] + [1 << (n - 1)])
            values.append(value)
        print('k=%d n=%d..16 values=%s' % (k, max(2, k + 1), values))
    print('Fidelity path formula: 15 cases passed.')
    print('Independent restricted-growth partition search: 22 cases passed.')
    for k in range(2, 5):
        top = max(v for (kk, n), v in table.items() if kk == k)
        witnesses = []
        for r in range(k + 1, top + 1):
            crossing = min(n for (kk, n), value in table.items() if kk == k and value >= r)
            assert table[k, crossing] == r
            assert crossing <= 2 * k * r
            witnesses.append((r - k, crossing))
        print('Spectrum k=%d, (t, first n)=%s' % (k, witnesses))
    assert [table[2, n] for n in range(3, 13)] == [3, 3, 3, 3, 4, 4, 5, 6, 6, 6]
    print('Family checks: 54 path powers, k=1..4, n<=16; all feasible colors in each tested prefix attained.')
    print('Scope: finite verification only; uniform k,t realization requires the mathematical proof.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
