#!/usr/bin/env python3
"""Exact finite checks for Benedetti–Menara, arXiv:2611.08071v1, Q9.1.

Only Python integers, finite sets, and graph search are used. The infinite
family argument is mathematical; checking a finite range does not prove it.
"""

from collections import deque
from itertools import combinations, product
import json


def powerset(s):
    values = sorted(s)
    return {frozenset(t) for r in range(len(values) + 1)
            for t in combinations(values, r)}


def closure(facets):
    return set().union(*(powerset(f) for f in facets))


def maximal(sets):
    return {s for s in sets if not any(s < t for t in sets)}


def dual(facets, rank):
    facets = sorted(facets, key=lambda f: tuple(sorted(f)))
    return {f: {g for g in facets if f != g and len(f & g) + 1 == rank}
            for f in facets}


def components(graph, vertices):
    todo = set(vertices)
    result = []
    while todo:
        root = min(todo, key=lambda f: tuple(sorted(f)))
        seen, pending = {root}, [root]
        while pending:
            for v in graph[pending.pop()] & todo - seen:
                seen.add(v)
                pending.append(v)
        todo -= seen
        result.append(seen)
    return result


def normal_failures(facets, rank):
    graph = dual(facets, rank)
    return [(s, components(graph, {f for f in facets if s <= f}))
            for s in sorted(closure(facets), key=lambda t: (len(t), sorted(t)))
            if len(components(graph, {f for f in facets if s <= f})) != 1]


def diameter(graph, vertices=None):
    vertices = set(graph) if vertices is None else set(vertices)
    assert len(components(graph, vertices)) == 1
    answer = 0
    for root in vertices:
        dist, queue = {root: 0}, deque([root])
        while queue:
            u = queue.popleft()
            for v in graph[u] & vertices - dist.keys():
                dist[v] = dist[u] + 1
                queue.append(v)
        assert len(dist) == len(vertices)
        answer = max(answer, max(dist.values()))
    return answer


def one_skeleton(facets):
    vertices = set().union(*facets)
    edges = {s for s in closure(facets) if len(s) == 2}
    return {v: {w for w in vertices if frozenset((v, w)) in edges}
            for v in vertices}


def completion_facets(facets):
    """Bron–Kerbosch: enumerate every maximal clique, without a size cutoff."""
    graph, found = one_skeleton(facets), set()

    def visit(r, p, x):
        if not p and not x:
            found.add(frozenset(r))
        for v in sorted(p.copy()):
            visit(r | {v}, p & graph[v], x & graph[v])
            p.remove(v)
            x.add(v)

    visit(set(), set(graph), set())
    return found


def fvector(facets):
    faces = closure(facets)
    return [sum(len(s) == r for s in faces)
            for r in range(1, max(map(len, facets)) + 1)]


def minimal_nonfaces(vertices, faces):
    return {s for s in powerset(vertices) if s not in faces
            and all(s - {v} in faces for v in s)}


def gap(facets, rank, k):
    """Literal G_k. Enumerate only the forbidden cardinalities."""
    vertices, faces = set().union(*facets), closure(facets)
    for size in range(3, rank - k + 1):
        for values in combinations(sorted(vertices), size):
            s = frozenset(values)
            if s not in faces and all(s - {v} in faces for v in s):
                return False
    return True


def fidelity():
    """Reproduce Section 8.1, not the proposed new counterexample."""
    for d in range(3, 7):
        sphere = {frozenset(2 * i + signs[i] for i in range(d))
                  for signs in product((0, 1), repeat=d)}
        removed = frozenset(range(0, 2 * d, 2))
        delta = sphere - {removed}
        faces = closure(delta)
        assert maximal(faces) == delta
        assert {len(f) for f in delta} == {d}
        assert not normal_failures(delta, d)
        assert gap(delta, d, 1)
        assert completion_facets(delta) == sphere
        assert not normal_failures(sphere, d)
        expected_missing = {frozenset((2 * i, 2 * i + 1)) for i in range(d)}
        assert minimal_nonfaces(range(2 * d), faces) == expected_missing | {removed}
        graph = dual(delta, d)
        assert diameter(graph) == d
        local = max(diameter(graph, {f for f in delta if s <= f})
                    for s in faces if len(s) == d - 2)
        assert local == 2
        print(json.dumps(dict(fidelity="Section 8.1, cross-polytope minus facet",
                              d=d, vertices=2*d, facets=len(delta), normal=True,
                              G1=True, completion_facets=len(sphere),
                              dual_diameter=d, local_diameter_1=local)))


def family(n):
    a, z, w, b = range(n, n + 4)
    cycle = [(i, (i + 1) % n) for i in range(n)]
    ordered = [frozenset((z, w, i, j)) for i, j in cycle]
    ordered += [frozenset((a, z, i, i + 1)) for i in range(4)]
    wb_edges = [(0, 1)] + [(i, (i + 1) % n) for i in range(n - 1, 2, -1)]
    ordered += [frozenset((w, b, i, j)) for i, j in wb_edges]
    ordered += [frozenset((a, b, 0, 1))]
    return ordered, (a, z, w, b)


def check_family(n):
    ordered, (a, z, w, b) = family(n)
    facets, faces = set(ordered), closure(ordered)
    assert len(facets) == len(ordered) == 2 * n + 3
    assert maximal(faces) == facets
    assert {len(f) for f in facets} == {4}
    assert len(set().union(*facets)) == n + 4
    attachments = []
    for i, f in enumerate(ordered[1:], 1):
        generators = maximal({f & old for old in ordered[:i]})
        assert generators and {len(s) for s in generators} == {3}
        attachments.append([sorted(s) for s in sorted(generators, key=lambda s: sorted(s))])
    # Independently test every face-star, including the empty-face star.
    assert not normal_failures(facets, 4)
    assert gap(facets, 4, 2)
    assert not gap(facets, 4, 1)
    hat = completion_facets(facets)
    extra = frozenset((a, b, 3, 4))
    assert hat == facets | {extra}
    assert {len(f) for f in hat} == {4}
    hatfaces = closure(hat)
    assert hatfaces - faces == {extra, frozenset((a, b, 3)), frozenset((a, b, 4))}
    edge = frozenset((a, b))
    star = {f for f in hat if edge <= f}
    expected_star = {frozenset((a, b, 0, 1)), extra}
    assert star == expected_star
    assert len(components(dual(hat, 4), star)) == 2
    assert not normal_failures(facets, 4) and normal_failures(hat, 4)
    assert len(components(dual(hat, 4), hat)) == 1
    assert fvector(facets) == [n + 4, 4 * n + 8, 5 * n + 8, 2 * n + 3]
    assert fvector(hat) == [n + 4, 4 * n + 8, 5 * n + 10, 2 * n + 4]
    print(json.dumps(dict(N=n, d=4, k=2, delta_f=fvector(facets), hat_f=fvector(hat),
                          delta_pure=True, delta_normal=True, G2=True,
                          hat_pure=True, hat_normal=False, hat_dual_connected=True,
                          bad_face=sorted(edge), bad_star_facets=sorted(map(sorted, star)),
                          bad_star_components=2)))
    if n == 6:
        assert ordered == [frozenset(t) for t in
                           [(7,8,0,1),(7,8,1,2),(7,8,2,3),(7,8,3,4),(7,8,4,5),
                            (7,8,5,0),(6,7,0,1),(6,7,1,2),(6,7,2,3),(6,7,3,4),
                            (8,9,0,1),(8,9,5,0),(8,9,4,5),(8,9,3,4),(6,9,0,1)]]
        print(json.dumps(dict(N6_ridge_attachment_generators=attachments)))


def main():
    fidelity()
    for n in range(6, 25):
        check_family(n)
    print("FAMILY_RANGE_CHECKED: 6 <= N <= 24; arbitrary N requires the stated argument")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
