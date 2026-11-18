#!/usr/bin/env python3
"""Exact checks for Clow–van Bommel, arXiv:2308.00054, Conjecture 8.2."""

from collections import deque
from functools import cache
from itertools import combinations_with_replacement


def path_edges(order):
    return [(v, v + 1) for v in range(order - 1)]


def tree_distances(order, edges):
    """Check the tree condition and compute distances by BFS."""
    assert len(edges) == order - 1
    assert len({tuple(sorted(edge)) for edge in edges}) == len(edges)
    neighbors = [[] for _ in range(order)]
    for u, v in edges:
        assert u != v
        neighbors[u].append(v)
        neighbors[v].append(u)

    matrix = []
    for start in range(order):
        dist = [-1] * order
        dist[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in neighbors[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        assert min(dist) >= 0
        matrix.append(dist)
    return matrix


def pendant_paths(core_order, core_edges, length):
    blocks = [
        tuple(root * (length + 1) + i for i in range(length + 1))
        for root in range(core_order)
    ]
    edges = [(blocks[u][0], blocks[v][0]) for u, v in core_edges]
    for block in blocks:
        edges.extend(zip(block, block[1:]))
    return blocks, edges


def domination_number(dist, radius):
    order = len(dist)
    balls = [sum(1 << u for u, d in enumerate(row) if d <= radius) for row in dist]

    @cache
    def minimum(uncovered):
        if uncovered == 0:
            return 0
        # Branch over every possible dominator of one uncovered vertex.
        target = (uncovered & -uncovered).bit_length() - 1
        return 1 + min(
            minimum(uncovered & ~balls[guard])
            for guard, distance in enumerate(dist[target])
            if distance <= radius
        )

    return minimum((1 << order) - 1)


def movement_is_possible(dist, radius, old, new):
    """Match individual guard occurrences, including repeated vertices."""
    @cache
    def match(i, remaining):
        if not remaining:
            return True
        return any(
            dist[old[i]][new[j]] <= radius and match(i + 1, remaining ^ (1 << j))
            for j in range(len(new))
            if remaining & (1 << j)
        )

    return match(0, (1 << len(new)) - 1)


def eternal_number_by_game(dist, radius):
    """Compute the greatest set of positions surviving arbitrary attacks."""
    order = len(dist)
    for n in range(1, order + 1):
        states = list(combinations_with_replacement(range(order), n))
        occupied = [sum(1 << v for v in set(state)) for state in states]
        transitions = [[] for _ in states]
        for i, old in enumerate(states):
            for j in range(i, len(states)):
                if movement_is_possible(dist, radius, old, states[j]):
                    transitions[i].append(j)
                    if i != j:
                        transitions[j].append(i)
        full = (1 << order) - 1
        surviving = set(range(len(states)))
        while surviving:
            keep = set()
            for i in surviving:
                responses = 0
                for j in transitions[i]:
                    if j in surviving:
                        responses |= occupied[j]
                if responses == full:
                    keep.add(i)
            if keep == surviving:
                return n
            surviving = keep
    raise AssertionError("One guard per vertex must suffice.")


def check_block_strategy(dist, blocks, radius):
    """Check the proposed local response for every guard position and attack."""
    assert tuple(sorted(v for block in blocks for v in block)) == tuple(range(len(dist)))
    assert all(dist[u][u] == 0 for u in range(len(dist)))
    for block in blocks:
        for old_guard in block:
            for attacked in block:
                assert dist[old_guard][attacked] <= radius
                assert movement_is_possible(dist, radius, (old_guard,), (attacked,))


def main():
    # A published family from Theorem 7.2, solved by unrestricted game search.
    for k in (2, 3, 5):
        blocks, edges = pendant_paths(2, path_edges(2), k)
        order = len(blocks) * (k + 1)
        dist = tree_distances(order, edges)
        expected = (order + k) // (k + 1)
        actual = eternal_number_by_game(dist, k)
        assert actual == expected
        print(f"fidelity Theorem 7.2: k={k} n={order} eternal={actual}")

    # Check the arbitrary-core formula and the adjacent even case.
    for name, core_order, core_edges in (
        ("singleton", 1, []),
        ("four-vertex path", 4, path_edges(4)),
        ("four-vertex star", 4, [(0, v) for v in range(1, 4)]),
    ):
        core_dist = tree_distances(core_order, core_edges)
        for r in (1, 2, 3):
            for k in (2 * r, 2 * r + 1):
                blocks, edges = pendant_paths(core_order, core_edges, k)
                dist = tree_distances(core_order * (k + 1), edges)
                static = domination_number(dist, k)
                smaller = domination_number(dist, r)
                check_block_strategy(dist, blocks, k)
                extra = domination_number(core_dist, r) if k % 2 else 0
                assert static == core_order
                assert smaller == core_order + extra
                print(
                    f"core={name} r={r} k={k}: "
                    f"static={static} eternal={len(blocks)} smaller={smaller}"
                )

    for k, q in ((3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1)):
        core_order = q * k
        blocks, edges = pendant_paths(core_order, path_edges(core_order), k)
        order = core_order * (k + 1)
        dist = tree_distances(order, edges)
        assert max(sum(d == 1 for d in row) for row in dist) <= 3
        static = domination_number(dist, k)
        smaller = domination_number(dist, k // 2)
        check_block_strategy(dist, blocks, k)
        assert static == len(blocks) == q * k
        assert smaller == q * k + q
        print(f"refutation k={k} q={q} n={order}: ({static}, {len(blocks)}, {smaller})")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
