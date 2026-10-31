#!/usr/bin/env python3
"""Exact small-tree fidelity and consistency checks; no asymptotic proof."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
import sys
import networkx as nx

# Table 1 of Madaras--Paralic, arXiv:2609.09342v1, columns k=1,...,4.
SOURCE_ROWS = {
    1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (1, 0, 0, 0),
    4: (1, 1, 0, 0), 5: (3, 0, 0, 0), 6: (4, 2, 0, 0),
    7: (9, 1, 1, 0), 8: (18, 5, 0, 0), 9: (45, 2, 0, 0),
    10: (93, 10, 2, 1), 11: (225, 10, 0, 0), 12: (512, 39, 0, 0),
    13: (1244, 56, 1, 0),
}


def unique_path(adj, x, y):
    parent = {x: None}
    q = deque([x])
    while q:
        v = q.popleft()
        if v == y:
            break
        for w in adj[v]:
            if w not in parent:
                parent[w] = v
                q.append(w)
    path = [y]
    while path[-1] != x:
        path.append(parent[path[-1]])
    return path[::-1]


def deleted_components(adj, removed):
    unseen = set(adj) - {removed}
    result = []
    while unseen:
        stack = [unseen.pop()]
        component = set(stack)
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w in unseen:
                    unseen.remove(w)
                    component.add(w)
                    stack.append(w)
        result.append(component)
    return result


def check_tree(g):
    adj = {v: set(g[v]) for v in g}
    n = len(adj)
    N = n - 1
    # Literal shortest-path score: one shortest path for each unordered pair.
    score = dict.fromkeys(adj, 0)
    for x, y in combinations(adj, 2):
        for v in unique_path(adj, x, y)[1:-1]:
            score[v] += 1
    components = {v: deleted_components(adj, v) for v in adj}
    for v in adj:
        sizes = [len(c) for c in components[v]]
        assert sum(sizes) == N
        assert 2 * score[v] == N * N - sum(s * s for s in sizes)
    B = max(score.values())
    centers = {v for v in adj if score[v] == B}
    k = len(centers)
    if n == 1:
        assert score == {0: 0} and k == 1
        return k
    # Squaring is valid only when k-3 is positive.
    assert k <= 3 or (k - 3) ** 2 < 8 * N
    centroid = [v for v in adj if all(2 * len(c) <= n for c in components[v])]
    assert centroid
    assert 4 * B >= N * (N - 1)
    if k == 1:
        return k
    # Theorem 5: maxima lie on one path, with no requirement to fill its gaps.
    x, y = max(combinations(centers, 2), key=lambda xy: len(unique_path(adj, *xy)))
    path = unique_path(adj, x, y)
    assert centers <= set(path)
    selected = [v for v in path if v in centers]
    data = []
    for v in selected:
        left = next((c for c in components[v] if x in c), set())
        right = next((c for c in components[v] if y in c), set())
        off = [c for c in components[v] if c != left and c != right]
        a, b = len(left), len(right)
        r = sum(map(len, off))
        d = a - b
        assert N == a + b + r
        assert 4 * B == N * N + 2 * N * r - r * r - d * d - 2 * sum(len(c)**2 for c in off)
        assert d * d <= N * (2 * r + 1)
        data.append((a, b, r, d))
    for (a, b, r, d), (aa, bb, rr, dd) in zip(data, data[1:]):
        assert aa >= a + r + 1 and b >= bb + rr + 1
        delta = dd - d
        assert delta >= r + rr + 2 and delta >= 2
        assert Fraction(delta) >= 1 + Fraction(d * d + dd * dd, 2 * N)
        assert Fraction(delta) >= Fraction(d * dd, N)
    ds = [t[3] for t in data]
    # Verify the packing inequalities with all integer thresholds in [1,N].
    for h in range(1, N + 1):
        central = [d for d in ds if -h <= d <= h]
        assert len(central) <= h + 1
        for tail in ([d for d in ds if d > h], [-d for d in reversed(ds) if d < -h]):
            for xx, yy in zip(tail, tail[1:]):
                assert Fraction(1, xx) - Fraction(1, yy) >= Fraction(1, N)
            assert len(tail) < Fraction(N, h) + 1
        assert k < h + Fraction(2 * N, h) + 3
    return k


def main():
    assert nx.__version__ == '3.5', nx.__version__
    print('Exact arithmetic: integer path counts and Fraction; NetworkX 3.5 tree enumeration')
    example = nx.Graph([(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)])
    assert check_tree(example) == 3  # Source Example 1, seven-vertex binary tree.
    print('Fidelity Example 1: three maximizers, each of betweenness 9')
    # Explicit value check, independent of the assertion about their number.
    values = [sum(v in unique_path({u: set(example[u]) for u in example}, x, y)[1:-1]
                  for x, y in combinations(example, 2)) for v in [0, 1, 2]]
    assert values == [9, 9, 9]
    total = 0
    for n, expected in SOURCE_ROWS.items():
        trees = [nx.empty_graph(1)] if n == 1 else nx.nonisomorphic_trees(n)
        counts = Counter(check_tree(g) for g in trees)
        assert set(counts) <= {1, 2, 3, 4}
        row = tuple(counts[i] for i in range(1, 5))
        assert row == expected, (n, row, expected)
        checked = sum(row)
        total += checked
        print(f'n={n}: centers[1..4]={row}; trees={checked}; f(n)={max(counts)}')
    print(f'Checked all {total} nonisomorphic trees of orders 1 through 13')
    print('Cut identity, centroid bound, path containment, imbalance spacing, reciprocal packing: PASS')
    print('Proposed strict bound k < 2*sqrt(2*(n-1))+3: PASS (integer square comparison)')
    print('Finite checks do not certify the universal or asymptotic assertion')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
