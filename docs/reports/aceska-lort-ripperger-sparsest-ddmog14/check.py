#!/usr/bin/env python3
"""Exact checks for the Aceska--Lort--Ripperger sparse-14 conjecture.

The source Figure 6 is checked for fidelity (14 vertices, 24 edges, connected,
and zero DDM weights).  The nominated 23-edge graph is then checked against the
same literal oriented-graph, bijective-labeling, connectivity, and weight tests.
All arithmetic is integer arithmetic.
"""

from collections import deque


def check_graph(name, edges, expected_edges):
    vertices = set(range(1, 15))
    assert len(edges) == expected_edges, (name, len(edges), expected_edges)
    assert all(u in vertices and v in vertices and u != v for u, v in edges)
    assert len(set(edges)) == len(edges), f"{name}: repeated ordered edge"
    assert all((v, u) not in set(edges) for u, v in edges)

    undirected = {v: set() for v in vertices}
    for u, v in edges:
        undirected[u].add(v)
        undirected[v].add(u)
    seen = {1}
    queue = deque([1])
    while queue:
        u = queue.popleft()
        for v in undirected[u]:
            if v not in seen:
                seen.add(v)
                queue.append(v)
    assert seen == vertices, (name, "underlying graph is disconnected", seen)

    incoming = {v: [] for v in vertices}
    outgoing = {v: [] for v in vertices}
    for u, v in edges:
        outgoing[u].append(v)
        incoming[v].append(u)
    weights = {
        v: sum(incoming[v]) - sum(outgoing[v])
        for v in vertices
    }
    assert set(weights.values()) == {0}, (name, weights)
    print(f"{name}: vertices=14 edges={len(edges)} connected=True weights={weights}")


# Figure 6 (the source labels this construction Figure 6 in v3; it is
# Figure 3 in v1).  The first two vertices are the attached 2 K_1 nodes.
source_edges = [
    (3, 4), (4, 5), (5, 3),
    (6, 8), (8, 7), (7, 6),
    (9, 11), (11, 10), (10, 9),
    (12, 13), (13, 14), (14, 12),
    (3, 1), (2, 13), (14, 1),
    (12, 1), (2, 4), (5, 1),
    (7, 2), (1, 6), (1, 8),
    (10, 2), (1, 11), (1, 9),
]
check_graph("source Figure 6", source_edges, 24)

# Nominated counterexample: vertices are their labels 1,...,14.
counterexample_edges = [
    (1, 5), (1, 12), (2, 5), (2, 11), (3, 9),
    (4, 3), (4, 13), (5, 3), (6, 4), (6, 10),
    (7, 2), (7, 12), (8, 1), (8, 11), (9, 1),
    (9, 2), (10, 4), (10, 13), (11, 10), (12, 8),
    (13, 14), (14, 6), (14, 7),
]
check_graph("nominated 23-edge counterexample", counterexample_edges, 23)
assert 23 < 24
print("settlement: the nominated graph refutes the fixed-order 24-edge lower bound")
print("ALL_CHECKS_OK")
