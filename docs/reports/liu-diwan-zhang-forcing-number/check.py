#!/usr/bin/env python3
"""Exact definition-level checks for arXiv:2512.22761v1, Problem 5.5.

An n by n Boolean adjacency matrix represents a simple bipartite graph.
A permutation p represents the matching {(i, p[i])}. Forcing sets are
subsets of these edges, enumerated without using the proposed theorem.
Only Python integers and finite combinatorics are used.
"""
from itertools import permutations


def perfect_matchings(adj):
    n = len(adj)
    return [p for p in permutations(range(n))
            if all(adj[i][p[i]] for i in range(n))]


def is_forcing(p, chosen, matchings):
    return all(q == p or any((chosen >> i) & 1 and p[i] != q[i]
                            for i in range(len(p))) for q in matchings)


def forcing_data(matchings, n):
    masks = sorted(range(1 << n), key=lambda x: (x.bit_count(), x))
    answer = []
    for p in matchings:
        s = next(s for s in masks if is_forcing(p, s, matchings))
        answer.append((s.bit_count(), s))
    return answer


def interval_check(adj, matchings, data):
    n = len(adj)
    r, selected = min(data)
    p = matchings[data.index((r, selected))]
    residual = [i for i in range(n) if not (selected >> i) & 1]
    outgoing = {i: {j for j in residual if j != i and adj[i][p[j]]}
                for i in residual}
    order = []
    remaining = set(residual)
    while remaining:
        sources = [i for i in sorted(remaining)
                   if not any(i in outgoing[j] for j in remaining)]
        assert sources, ('forcing deletion has a directed cycle', adj, p, selected)
        v = sources[0]
        order.append(v)
        remaining.remove(v)
    m = len(order)
    left = [0] * n
    right = [m + 1] * n
    for rank, i in enumerate(order, 1):
        left[i] = rank
        right[p[i]] = rank
    assert all(not adj[i][j] or left[i] <= right[j]
               for i in range(n) for j in range(n))
    class_count = 0
    for q in matchings:
        intervals = [(left[i], right[q[i]], i) for i in range(n)]
        assert all(a <= b for a, b, _ in intervals)
        # Twice the coordinate samples every endpoint and every intervening gap.
        for t in range(2 * (m + 1) + 1):
            overlap = sum(2 * a <= t <= 2 * b for a, b, _ in intervals)
            expected = r + 1 if t % 2 == 0 and 0 < t < 2 * (m + 1) else r
            assert overlap == expected, ('overlap', intervals, t, overlap, expected)
        classes = [[] for _ in range(r + 1)]
        last_right = [-1] * (r + 1)
        for a, b, i in sorted(intervals):
            available = [c for c in range(r + 1) if last_right[c] < a]
            assert available, ('closed interval coloring', intervals)
            c = available[0]
            classes[c].append(i)
            last_right[c] = b
        assert sum(map(len, classes)) == n
        for group in classes:
            # Independent definition-level uniqueness check for the induced subgraph.
            columns = tuple(q[i] for i in group)
            induced_matchings = sum(all(adj[i][j] for i, j in zip(group, perm))
                                    for perm in permutations(columns))
            assert induced_matchings == 1
            for earlier, later in zip(group, group[1:]):
                assert right[q[earlier]] < left[later]
            chosen_complement = ((1 << n) - 1) ^ sum(1 << i for i in group)
            assert is_forcing(q, chosen_complement, matchings)
            class_count += bool(group)
    return len(matchings), class_count


def sharpness_graph(n, k):
    assert 2 <= k <= n - 1 and n % k == 1
    sizes = [n // k] * (k - 1) + [n // k + 1]
    starts = [0]
    for size in sizes:
        starts.append(starts[-1] + size)
    adj = [[False] * n for _ in range(n)]
    for block in range(k):
        for i in range(starts[block], starts[block + 1]):
            for j in range(starts[block], starts[block + 1]):
                adj[i][j] = True
        next_block = (block + 1) % k
        for t in range(min(sizes[block], sizes[next_block])):
            adj[starts[block] + t][starts[next_block] + t] = True
    return adj


def main():
    graph_total = matching_total = class_total = implication_total = 0
    for n in range(1, 5):
        per_graph = per_matching = per_class = 0
        off_diagonal = [(i, j) for i in range(n) for j in range(n) if i != j]
        for bits in range(1 << len(off_diagonal)):
            adj = [[i == j for j in range(n)] for i in range(n)]
            for t, (i, j) in enumerate(off_diagonal):
                adj[i][j] = bool((bits >> t) & 1)
            matchings = perfect_matchings(adj)
            assert matchings
            data = forcing_data(matchings, n)
            f = min(r for r, _ in data)
            F = max(r for r, _ in data)
            # Fidelity: Corollary 2.2, with all matchings computed literally.
            assert (F == n - 1) == all(map(all, adj))
            # Fidelity: Theorem 5.1.
            if F == n - 2:
                assert f >= (n + 1) // 2 - 1
            # Proposed uniform settlement, stronger than each fixed-k implication.
            assert (f + 1) * (n - F) >= n
            for k in range(1, n):
                if F >= n - k:
                    assert f >= (n + k - 1) // k - 1
                    implication_total += 1
            matched, classes = interval_check(adj, matchings, data)
            per_graph += 1
            per_matching += matched
            per_class += classes
        graph_total += per_graph
        matching_total += per_matching
        class_total += per_class
        print(f'n={n}: graphs={per_graph}, matching_forcing_numbers={per_matching}, '
              f'nonempty_interval_classes={per_class}')
    print(f'Fidelity Corollary 2.2 and Theorem 5.1: {graph_total} graphs passed')
    print(f'Product inequality: {graph_total} graphs passed; violations=0')
    print(f'Published fixed-k implication: {implication_total} applicable pairs passed')
    print(f'Interval construction: {matching_total} matchings, {class_total} nonempty classes passed')
    for n, k in [(3, 2), (4, 3), (5, 2), (7, 3)]:
        adj = sharpness_graph(n, k)
        matchings = perfect_matchings(adj)
        data = forcing_data(matchings, n)
        f, F = min(r for r, _ in data), max(r for r, _ in data)
        assert (f, F) == ((n + k - 1) // k - 1, n - k)
        print(f'Fidelity Remark 5.6: n={n}, k={k}, perfect_matchings={len(matchings)}, f={f}, F={F}')
    assert (graph_total, matching_total) == (4165, 16774)
    print('Finite sanity checks only; no universal conclusion follows from enumeration.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
