#!/usr/bin/env python3
"""Exact checks for Antony--Lafond--Monti--Sinaimeri Problem 3.

All coordinates, labels, intervals, and distances are integers.  The script
checks the source's skew-grid labelling, the seven-vertex obstruction used for
the negative direction, and explicit caterpillar witnesses for every tested
thin strong grid.
"""
from __future__ import annotations

from collections import deque
from itertools import combinations, permutations
import json
from pathlib import Path

REPORT = Path(__file__).with_name("check-output.json")


def unordered(u, v):
    return frozenset((u, v))


def strong_adj(u: tuple[int, int], v: tuple[int, int]) -> bool:
    return u != v and max(abs(u[0] - v[0]), abs(u[1] - v[1])) == 1


def interval_labeling(vertices, adjacent, labels, interval):
    lo, hi = interval
    assert 0 <= lo <= hi
    for u, v in combinations(vertices, 2):
        actual = lo <= abs(labels[u] - labels[v]) <= hi
        assert actual == adjacent(u, v), (u, v, labels[u], labels[v], interval)


def check_source_skew_grid() -> dict:
    """Fidelity: verify the paper's [n+1,n+2] skew-grid labelling."""
    tested = []
    for n, m in [(1, 1), (2, 3), (3, 2), (4, 5), (5, 4)]:
        vertices = [
            (i, j)
            for i in range(1, n + 1)
            for j in range(n - i + 1, n - i + m + 1)
        ]

        def adjacent(u, v):
            return u != v and abs(u[0] - v[0]) + abs(u[1] - v[1]) == 1

        labels = {
            v: (v[1] - n + v[0] - 1) * (n + 1) + v[0] for v in vertices
        }
        interval_labeling(vertices, adjacent, labels, (n + 1, n + 2))
        tested.append({"n": n, "m": m, "vertices": len(vertices)})
    return {"source_skew_grid_cases": tested}


def obstruction_graph():
    # H is C4 a-b-c-d-a with ears x on ab, y on bc, z on cd.
    names = ("a", "b", "c", "d", "x", "y", "z")
    edge_names = (
        "ab", "bc", "cd", "da", "xa", "xb", "yb", "yc", "zc", "zd"
    )
    edges = {unordered(e[0], e[1]) for e in edge_names}

    def adjacent(u, v):
        return unordered(u, v) in edges

    return names, edges, adjacent


def one_sided_consecutive(order, adjacent) -> bool:
    """Necessary interval-difference condition from the source's Lemma 2."""
    for i, x in enumerate(order):
        for side in (order[:i], order[i + 1 :]):
            neighbor_indices = [j for j, y in enumerate(side) if adjacent(x, y)]
            if neighbor_indices and neighbor_indices[-1] - neighbor_indices[0] + 1 != len(neighbor_indices):
                return False
    return True


def check_obstruction() -> dict:
    names, edges, adjacent = obstruction_graph()
    coordinates = {
        "a": (2, 1),
        "b": (1, 0),
        "c": (0, 1),
        "d": (1, 2),
        "x": (2, 0),
        "y": (0, 0),
        "z": (0, 2),
    }
    assert len(edges) == 10
    assert all(strong_adj(coordinates[u], coordinates[v]) == adjacent(u, v)
               for u, v in combinations(names, 2))
    assert all(
        any(adjacent(u, w) != adjacent(v, w) for w in names if w not in (u, v))
        for u, v in combinations(names, 2)
    )

    orders = list(permutations(names))
    permissible = [o for o in orders if one_sided_consecutive(o, adjacent)]
    assert len(orders) == 5040
    assert not permissible

    # The source proof reduces the 16 permissible cyclic-core orders to five.
    core = [o for o in permutations(("a", "b", "c", "d"))
            if one_sided_consecutive(o, adjacent)]
    assert len(core) == 16
    automorphism = {"a": "d", "d": "a", "b": "c", "c": "b", "x": "z", "z": "x", "y": "y"}
    representatives = [("a", "b", "d", "c"), ("a", "c", "b", "d"),
                       ("a", "c", "d", "b"), ("a", "d", "b", "c"),
                       ("b", "d", "a", "c")]
    covered = set()
    for p in representatives:
        covered.update((p, p[::-1], tuple(automorphism[v] for v in p),
                        tuple(automorphism[v] for v in p[::-1])))
    assert covered == set(core)
    return {
        "obstruction_vertices": list(names),
        "embedding_coordinates": coordinates,
        "edges": sorted([sorted(e) for e in edges]),
        "nonedges": 11,
        "twin_pairs": [],
        "checked_orders": len(orders),
        "permissible_orders": len(permissible),
        "permissible_cycle_orders": len(core),
        "cycle_order_representatives": ["abdc", "acbd", "acdb", "adbc", "bdac"],
    }


def grid(m, n):
    vertices = [(i, j) for i in range(m) for j in range(n)]
    return vertices, strong_adj


def tree_distances(tree, start):
    distance = {start: 0}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in tree[u]:
            if v not in distance:
                distance[v] = distance[u] + 1
                queue.append(v)
    return distance


def check_thin_witness(m, n):
    # Orient the spine along the dimension with k columns; the other side has
    # one or two leaves per spine position.
    rows, k = (m, n) if m <= 2 else (n, m)
    assert rows in (1, 2) and rows * k >= 2
    graph_vertices = [(i, j) for i in range(rows) for j in range(k)]
    spine = [("s", j) for j in range(k)]
    tree = {u: set() for u in spine}
    leaf_of = {}
    for i, j in graph_vertices:
        leaf = ("v", i, j)
        leaf_of[(i, j)] = leaf
        tree[leaf] = set()
        tree[leaf].add(spine[j])
        tree[spine[j]].add(leaf)
    for j in range(k - 1):
        tree[spine[j]].add(spine[j + 1])
        tree[spine[j + 1]].add(spine[j])
    leaves = {u for u, nbrs in tree.items() if len(nbrs) == 1}
    assert leaves == set(leaf_of.values())

    def tree_adj(u, v):
        return u != v and v in tree[u]

    def graph_adj(u, v):
        return strong_adj(u, v)

    for u, v in combinations(graph_vertices, 2):
        du = tree_distances(tree, leaf_of[u])[leaf_of[v]]
        assert (2 <= du <= 3) == graph_adj(u, v)
    return {"m": m, "n": n, "vertices": len(graph_vertices), "interval": [2, 3]}


def check_all_thin_cases():
    cases = []
    for m in range(1, 7):
        for n in range(1, 7):
            if m * n >= 2 and min(m, n) <= 2:
                vertices, adjacent = grid(m, n)
                labels = {v: (v[0] + 1 if n <= 2 else v[1] + 1) for v in vertices}
                interval_labeling(vertices, adjacent, labels, (0, 1))
                cases.append(check_thin_witness(m, n))
    assert len(cases) == 19
    return cases


def check_large_embeddings():
    obstruction = {
        "a": (2, 1), "b": (1, 0), "c": (0, 1), "d": (1, 2),
        "x": (2, 0), "y": (0, 0), "z": (0, 2),
    }
    for m in range(3, 7):
        for n in range(3, 7):
            vertices = {(i, j) for i in range(m) for j in range(n)}
            assert set(obstruction.values()) <= vertices
            def host_adj(u, v):
                return u != v and max(abs(u[0] - v[0]), abs(u[1] - v[1])) == 1
            names, edges, adjacent = obstruction_graph()
            for u, v in combinations(names, 2):
                assert host_adj(obstruction[u], obstruction[v]) == adjacent(u, v)
    return {"tested_dimensions": 16}


def main():
    result = {}
    result.update(check_source_skew_grid())
    result["obstruction"] = check_obstruction()
    result["thin_witnesses"] = check_all_thin_cases()
    result["large_grid_embeddings"] = check_large_embeddings()
    REPORT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
