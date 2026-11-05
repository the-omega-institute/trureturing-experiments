#!/usr/bin/env python3
"""Exact small-graph checks for Armfield--Mailler, arXiv:2610.18453v1, §1.2.

The symbolic checks use Fraction throughout.  The two simulations use only integer
random draws; their output is a sanity reading, not proof evidence.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from random import Random
from typing import Dict, Iterable, List, Sequence, Tuple

Vertex = int
Edge = Tuple[Vertex, Vertex]


def edge(u: Vertex, v: Vertex) -> Edge:
    return (u, v) if u < v else (v, u)


def connected(vertices: Sequence[Vertex], edges: Sequence[Edge]) -> bool:
    adj = {v: set() for v in vertices}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    seen = {vertices[0]}
    todo = [vertices[0]]
    while todo:
        u = todo.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                todo.append(v)
    return len(seen) == len(vertices)


def neighbours(vertices: Sequence[Vertex], edges: Sequence[Edge]) -> Dict[Vertex, List[Vertex]]:
    adj = {v: [] for v in vertices}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def solve_hitting(vertices: Sequence[Vertex], edges: Sequence[Edge], weights: Dict[Edge, int], start: Vertex, target: Vertex, avoid: Vertex) -> Fraction:
    """P_start(tau_target < tau_avoid), by exact rational Gaussian elimination."""
    adj = neighbours(vertices, edges)
    unknown = [v for v in vertices if v not in (target, avoid)]
    pos = {v: i for i, v in enumerate(unknown)}
    n = len(unknown)
    a = [[Fraction(0) for _ in range(n + 1)] for _ in range(n)]
    for v in unknown:
        i = pos[v]
        total = sum(weights[edge(v, u)] for u in adj[v])
        a[i][i] = Fraction(total)
        for u in adj[v]:
            w = weights[edge(v, u)]
            if u == target:
                a[i][-1] += Fraction(w)
            elif u != avoid:
                a[i][pos[u]] -= Fraction(w)
    for col in range(n):
        pivot = next(r for r in range(col, n) if a[r][col] != 0)
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [x / scale for x in a[col]]
        for r in range(n):
            if r == col or a[r][col] == 0:
                continue
            scale = a[r][col]
            a[r] = [a[r][j] - scale * a[col][j] for j in range(n + 1)]
    if start == target:
        return Fraction(1)
    if start == avoid:
        return Fraction(0)
    return a[pos[start]][-1]


def terminal_identity(vertices: Sequence[Vertex], edges: Sequence[Edge], weights: Dict[Edge, int], nest: Vertex, food: Vertex) -> None:
    adj = neighbours(vertices, edges)
    h = {v: solve_hitting(vertices, edges, weights, v, food, nest) for v in vertices}
    conductance = sum(Fraction(weights[edge(nest, u)]) * h[u] for u in adj[nest])
    assert conductance > 0
    common = sorted(set(adj[nest]).intersection(adj[food]))
    assert common, (vertices, edges)
    for i in common:
        a = Fraction(weights[edge(nest, i)])
        b = Fraction(weights[edge(i, food)])
        p_n = a * h[i] / conductance
        p_f = b * (1 - h[i]) / conductance
        assert p_n / a + p_f / b == 1 / conductance, (vertices, edges, i)
        assert Fraction(0) <= p_n <= 1 and Fraction(0) <= p_f <= 1


def enumerate_small_graphs() -> int:
    checked = 0
    for m in range(3, 6):
        vertices = tuple(range(m))
        nest, food = 0, 1
        candidates = [edge(u, v) for u, v in combinations(vertices, 2) if edge(u, v) != edge(nest, food)]
        for mask in range(1 << len(candidates)):
            edges = [candidates[j] for j in range(len(candidates)) if mask & (1 << j)]
            adj = neighbours(vertices, edges)
            if not connected(vertices, edges) or not set(adj[nest]).intersection(adj[food]):
                continue
            # Every tested graph has dist(nest, food)=2 because the direct edge is excluded.
            for variant in (0, 1):
                weights = {e: 1 + ((sum(e) + variant) % 3) for e in edges}
                terminal_identity(vertices, edges, weights, nest, food)
                checked += 1
    return checked


def loop_erase(path: Sequence[Vertex]) -> List[Vertex]:
    out: List[Vertex] = []
    where: Dict[Vertex, int] = {}
    for v in path:
        if v in where:
            out = out[: where[v] + 1]
            where = {x: i for i, x in enumerate(out)}
        else:
            where[v] = len(out)
            out.append(v)
    return out


def simulate(vertices: Sequence[Vertex], edges: Sequence[Edge], nest: Vertex, food: Vertex, ants: int, seed: int) -> Tuple[Dict[Edge, int], List[Tuple[int, List[Fraction]]]]:
    rng = Random(seed)
    adj = neighbours(vertices, edges)
    weights = {e: 1 for e in edges}
    terminal = [e for e in edges if nest in e or food in e]
    snapshots: List[Tuple[int, List[Fraction]]] = []
    for n in range(ants):
        path = [nest]
        current = nest
        while current != food:
            choices = [(u, weights[edge(current, u)]) for u in adj[current]]
            total = sum(w for _, w in choices)
            draw = rng.randrange(total)
            for u, w in choices:
                if draw < w:
                    current = u
                    break
                draw -= w
            path.append(current)
        erased = loop_erase(path)
        for u, v in zip(erased, erased[1:]):
            weights[edge(u, v)] += 1
        if n + 1 in (100, 1000, ants):
            vals = [Fraction(weights[e], n + 2) for e in terminal]
            snapshots.append((n + 1, vals))
    return weights, snapshots


def main() -> None:
    checked = enumerate_small_graphs()
    print(f"exact terminal-edge identity checks: {checked}")
    c4 = [edge(0, 1), edge(1, 2), edge(2, 3), edge(0, 3)]
    k23 = [edge(u, v) for u in (0, 1) for v in (2, 3, 4)]
    for name, vertices, edges, nest, food, seed in (
        ("C4 opposite", tuple(range(4)), c4, 0, 2, 2610),
        ("K2,3", tuple(range(5)), k23, 0, 1, 18453),
    ):
        weights, spread = simulate(vertices, edges, nest, food, 5000, seed)
        print(f"simulation {name}: terminal normalized vectors (n=100,1000,5000) = {spread}")
        assert all(all(x >= 0 for x in vals) for _, vals in spread)
        assert max(weights.values()) > 1
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
