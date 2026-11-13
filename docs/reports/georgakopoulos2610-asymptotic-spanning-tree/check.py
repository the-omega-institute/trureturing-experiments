#!/usr/bin/env python3
"""Exact small checks for the proposed Cayley-graph refutation.

The script uses only integer arithmetic.  It checks reduced-word local structure
of Z^2 * Z and exhausts all spanning trees of the 2x2 and 3x3 square grids.
The finite grid calculation is a bounded sanity check for the one-ended-factor
lemma; the infinite conclusion is the mathematical argument in the draft.
"""

from __future__ import annotations

from collections import deque
from itertools import combinations, product
from typing import Iterable


Word = tuple[tuple[str, tuple[int, ...]], ...]
ID: Word = ()


def reduce_append(w: Word, letter: tuple[str, tuple[int, ...]]) -> Word:
    """Right multiply a reduced free-product word by one factor element."""
    tag, value = letter
    if not w:
        return ((tag, value),)
    last_tag, last_value = w[-1]
    if last_tag != tag:
        return w + ((tag, value),)
    if tag == "A":
        new_value = (last_value[0] + value[0], last_value[1] + value[1])
        if new_value == (0, 0):
            return w[:-1]
        return w[:-1] + (("A", new_value),)
    new_value = (last_value[0] + value[0],)
    if new_value == (0,):
        return w[:-1]
    return w[:-1] + (("B", new_value),)


GENERATORS = tuple(
    [("A", v) for v in ((1, 0), (-1, 0), (0, 1), (0, -1))]
    + [("B", (1,)), ("B", (-1,))]
)


def ball(radius: int) -> set[Word]:
    out = {ID}
    frontier = {ID}
    for _ in range(radius):
        nxt = {reduce_append(w, g) for w in frontier for g in GENERATORS}
        out |= nxt
        frontier = nxt
    return out


def word_length(w: Word) -> int:
    return sum(1 for _ in w)  # only used for the BFS radius bound


def first_factor(w: Word) -> str | None:
    return None if not w else w[0][0]


def check_free_product_structure() -> None:
    words = ball(5)
    assert ID in words and len(words) > 1
    # The six generators are inverse-closed and produce both factors.
    assert {(tag, value) for tag, value in GENERATORS} == {
        ("A", (1, 0)), ("A", (-1, 0)), ("A", (0, 1)), ("A", (0, -1)),
        ("B", (1,)), ("B", (-1,)),
    }
    # Away from the identity, right multiplication preserves the first factor.
    # This is the finite-radius check of the cut-vertex syllable argument.
    checked_edges = 0
    for w in words:
        if w == ID:
            continue
        for g in GENERATORS:
            z = reduce_append(w, g)
            if z != ID and z in words:
                assert first_factor(z) == first_factor(w)
                checked_edges += 1
    assert checked_edges > 0
    # The factor inclusions are isometric on this ball: projection to Z^2 gives
    # the lower bound and the A-generators give the matching path.
    a_vertices = [w for w in words if len(w) == 1 and w[0][0] == "A"]
    for u in a_vertices:
        ux, uy = u[0][1]
        for v in a_vertices:
            vx, vy = v[0][1]
            target = (abs(ux - vx) + abs(uy - vy))
            # BFS in the whole Cayley graph, restricted to the radius-5 ball.
            q: deque[Word] = deque([u])
            dist = {u: 0}
            while q:
                x = q.popleft()
                if x == v:
                    break
                for g in GENERATORS:
                    y = reduce_append(x, g)
                    if y in words and y not in dist:
                        dist[y] = dist[x] + 1
                        q.append(y)
            assert dist[v] == target


def grid_edges(n: int) -> list[tuple[int, int]]:
    def idx(x: int, y: int) -> int:
        return y * n + x

    out = []
    for y in range(n):
        for x in range(n):
            if x + 1 < n:
                out.append((idx(x, y), idx(x + 1, y)))
            if y + 1 < n:
                out.append((idx(x, y), idx(x, y + 1)))
    return out


def is_tree(nv: int, edges: Iterable[tuple[int, int]]) -> bool:
    parent = list(range(nv))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    used = 0
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        parent[ra] = rb
        used += 1
    return used == nv - 1 and len({find(i) for i in range(nv)}) == 1


def tree_distance(nv: int, edges: Iterable[tuple[int, int]], start: int) -> list[int]:
    adj = [[] for _ in range(nv)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    dist = [-1] * nv
    dist[start] = 0
    q: deque[int] = deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def grid_stretch(n: int) -> tuple[int, int, int]:
    edges = grid_edges(n)
    nv = n * n
    trees = 0
    best = None
    for chosen in combinations(edges, nv - 1):
        if not is_tree(nv, chosen):
            continue
        trees += 1
        distances = [tree_distance(nv, chosen, u) for u in range(nv)]
        stretch = max(distances[a][b] for a, b in edges)
        best = stretch if best is None else min(best, stretch)
    assert trees > 0 and best is not None
    return trees, best, len(edges)


def main() -> None:
    check_free_product_structure()
    print("free-product ball radius=5: syllable cut and Z^2 isometry checks passed")
    results = {}
    for n in (2, 3):
        results[n] = grid_stretch(n)
        trees, best, edge_count = results[n]
        print(f"grid {n}x{n}: edges={edge_count}, spanning_trees={trees}, minimum_edge_stretch={best}")
    # Every 3x3 spanning tree stretches at least one ambient edge by 3.
    assert results[3][1] >= 3
    print("bounded one-ended-factor obstruction check passed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
