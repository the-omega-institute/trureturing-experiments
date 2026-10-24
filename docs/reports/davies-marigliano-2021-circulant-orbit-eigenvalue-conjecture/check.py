#!/usr/bin/env python3
"""Exact finite checks for Davies–Marigliano Conjecture 4.4 (arXiv:2012.01905v5).

Python 3, SymPy and NetworkX. All mathematical evidence is integer arithmetic
or exact polynomial algebra. Automorphisms are exhaustively enumerated only
for cycles and the 30-vertex base example; other outputs are labeled formula
checks. The infinite-family argument is not established by finite enumeration.
"""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, product
from math import gcd, prod
import json

import networkx as nx
import sympy as sp

X = sp.Symbol("x")


def emit(label, value):
    print(label + " " + json.dumps(value, sort_keys=True))


def matrix(graph):
    vertices = list(graph)
    return sp.Matrix([[int(graph.has_edge(u, v)) for v in vertices] for u in vertices])


def spectral_formula(sizes):
    degree = sum(m - 1 for m in sizes)
    values = Counter()
    for mask in product((0, 1), repeat=len(sizes)):
        value = degree - sum(m for m, bit in zip(sizes, mask) if bit)
        values[value] += prod(m - 1 for m, bit in zip(sizes, mask) if bit)
    assert sum(values.values()) == prod(sizes)
    return values


def polynomial(values):
    return sp.Poly(prod((X - a) ** mult for a, mult in values.items()), X)


def complete_product(sizes):
    vertices = list(product(*(range(m) for m in sizes)))
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    for u, v in combinations_with_replacement(vertices, 2):
        if sum(a != b for a, b in zip(u, v)) == 1:
            graph.add_edge(u, v)
    return graph


def all_pair_orbits(graph):
    """Full uncoloured graph automorphism search, without product assumptions."""
    vertices = list(graph)
    rank = {u: i for i, u in enumerate(vertices)}
    pairs = list(combinations_with_replacement(range(len(vertices)), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    parent = list(range(len(pairs)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    count = 0
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph, graph)
    for automorphism in matcher.isomorphisms_iter():
        count += 1
        perm = [rank[automorphism[u]] for u in vertices]
        # Explicitly check every returned permutation against ALL adjacency bits.
        assert len(set(perm)) == len(vertices)
        assert all(graph.has_edge(u, v) == graph.has_edge(automorphism[u], automorphism[v])
                   for u, v in combinations_with_replacement(vertices, 2))
        for i, (a, b) in enumerate(pairs):
            image = tuple(sorted((perm[a], perm[b])))
            parent[root(i)] = root(index[image])
    classes = defaultdict(list)
    for i, pair in enumerate(pairs):
        classes[root(i)].append(pair)
    return count, list(classes.values()), vertices


def main():
    emit("environment", {"sympy": sp.__version__, "networkx": nx.__version__,
                          "arithmetic": "integers and exact polynomials"})
    # Source fidelity: Theorem 3.4(1) and Example 3.5, with diagonal pairs.
    for n in range(3, 11):
        graph = nx.cycle_graph(n)
        charpoly = matrix(graph).charpoly(X).as_poly()
        r = charpoly.sqf_part().degree()
        aut_count, orbits, _ = all_pair_orbits(graph)
        assert len(orbits) == r == n // 2 + 1
        assert aut_count == 2 * n
        if n == 6:
            assert charpoly == sp.Poly((X - 2) * (X + 2) * (X - 1)**2 * (X + 1)**2, X)
            assert sorted(map(len, orbits)) == [3, 6, 6, 6]
        emit("fidelity_cycle", {"n": n, "automorphisms": aut_count,
                               "unordered_pair_orbits": len(orbits), "distinct_eigenvalues": r,
                               "charpoly": str(charpoly.as_expr())})

    for b in (3, 5, 7):
        sizes = (2, b, b + 2)
        graph = complete_product(sizes)
        assert nx.is_connected(graph)
        assert len(graph) == 2 * b * (b + 2)
        assert set(dict(graph.degree()).values()) == {2 * b + 1}
        common_counts = defaultdict(set)
        for u, v in graph.edges:
            coordinate = next(i for i in range(3) if u[i] != v[i])
            common_counts[coordinate].add(len(set(graph[u]) & set(graph[v])))
        assert dict(common_counts) == {i: {m - 2} for i, m in enumerate(sizes)}
        # The simultaneous translation really visits every vertex once.
        zero = (0, 0, 0)
        visited = set()
        u = zero
        while u not in visited:
            visited.add(u)
            u = tuple((a + 1) % m for a, m in zip(u, sizes))
        assert u == zero and len(visited) == len(graph)
        shift = {u: tuple((a + 1) % m for a, m in zip(u, sizes)) for u in graph}
        assert all(graph.has_edge(u, v) == graph.has_edge(shift[u], shift[v])
                   for u, v in combinations_with_replacement(graph, 2))
        values = spectral_formula(sizes)
        expected = {2*b+1: 1, 2*b-1: 1, b+1: b-1, b-1: 2*b,
                    b-3: b+1, -1: b*b-1, -3: b*b-1}
        assert dict(values) == expected and len(values) == 7
        direct = matrix(graph).charpoly(X).as_poly()
        assert direct == polynomial(values)
        result = {"b": b, "vertices": len(graph), "degree": 2*b+1,
                  "translation_cycle_length": len(visited), "spectrum": sorted(values.items()),
                  "edge_common_neighbor_counts": {i: sorted(v) for i, v in common_counts.items()},
                  "direct_charpoly_matches_tensor_formula": True}
        if b == 3:
            aut_count, orbits, vertices = all_pair_orbits(graph)
            assert aut_count == 1440 and len(orbits) == 8
            support_sets = []
            for orbit in orbits:
                supports = {tuple(i for i in range(3) if vertices[a][i] != vertices[c][i])
                            for a, c in orbit}
                assert len(supports) == 1
                support_sets.append(next(iter(supports)))
            assert len(set(support_sets)) == 8
            assert sorted(map(len, orbits)) == [15, 30, 30, 30, 60, 60, 120, 120]
            result.update(full_automorphism_count=aut_count, full_unordered_pair_orbits=len(orbits),
                          orbit_sizes=sorted(map(len, orbits)),
                          factored_charpoly=str(sp.factor(direct.as_expr())))
        emit("direct_product_check", result)

    # These are finite arithmetic checks of the proposed family formula.
    for b in range(3, 102, 2):
        sizes = (2, b, b + 2)
        assert all(gcd(sizes[i], sizes[j]) == 1 for i in range(3) for j in range(i))
        sums = Counter(sum(m for m, bit in zip(sizes, mask) if bit)
                       for mask in product((0, 1), repeat=3))
        assert set(sums) == {0, 2, b, b+2, b+4, 2*b+2, 2*b+4}
        assert sorted(sums.values()) == [1, 1, 1, 1, 1, 1, 2]
        assert len(spectral_formula(sizes)) == 7
    emit("family_formula_check", {"odd_b_min": 3, "odd_b_max": 101, "cases": 50,
                                  "subset_sum_count_each": 7,
                                  "orbit_formula_each": 8,
                                  "orbit_formula_is_not_exhaustively_enumerated_here": True})
    sizes = [2, 3, 5]
    for t in range(5):
        if t:
            new = int(sp.nextprime(sum(sizes)))
            assert new > sum(sizes)
            sizes.append(new)
        assert len(set(sizes)) == len(sizes)
        assert all(gcd(sizes[i], sizes[j]) == 1 for i in range(len(sizes)) for j in range(i))
        r = len(spectral_formula(sizes))
        s_formula = 2 ** len(sizes)
        assert r == 7 * 2**t and s_formula - r == 2**t
        emit("extension_subset_sum_check", {"sizes": sizes.copy(), "spectral_values": r,
                                            "orbit_formula": s_formula, "gap_formula": s_formula-r})
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
