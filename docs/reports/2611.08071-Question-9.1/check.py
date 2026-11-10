#!/usr/bin/env python3
"""Exact finite checks for Benedetti--Menara, arXiv:2611.08071v1, Question 9.1.

Original implementation from the mathematical definitions; standard library only.
Faces include the empty face. Vertices are integers; auxiliary labels are printed.
This is a reproducible computation, not a Lean proof or an infinite-family certificate.
"""
from itertools import combinations
from collections import Counter, deque


def powerset(s):
    s = sorted(s)
    return {frozenset(t) for k in range(len(s) + 1) for t in combinations(s, k)}


def closure(facets):
    return set().union(*(powerset(f) for f in facets))


def facets_of(faces):
    # A face is maximal iff it has no one-vertex coface in a downward closed family.
    nonmaximal = {f - {v} for f in faces for v in f}
    return faces - nonmaximal


def components(vertices, edges):
    vertices = set(vertices)
    adjacency = {v: set() for v in vertices}
    for a, b in edges:
        assert a in vertices and b in vertices and a != b
        adjacency[a].add(b)
        adjacency[b].add(a)
    answer = []
    unseen = set(vertices)
    while unseen:
        start = min(unseen)
        todo, seen = [start], {start}
        while todo:
            for v in adjacency[todo.pop()] - seen:
                seen.add(v)
                todo.append(v)
        answer.append(frozenset(seen))
        unseen -= seen
    return sorted(answer, key=lambda s: (len(s), sorted(s)))


def dual_edges(facets, d):
    return [(a, b) for a, b in combinations(sorted(facets, key=sorted), 2)
            if len(a & b) + 1 == d]


def normal_bad_faces(faces, d):
    facets = facets_of(faces)
    edges = dual_edges(facets, d)
    bad = []
    for sigma in sorted(faces, key=lambda f: (len(f), sorted(f))):
        nodes = {f for f in facets if sigma <= f}
        induced = [(a, b) for a, b in edges if a in nodes and b in nodes]
        if len(components(nodes, induced)) != 1:
            bad.append(tuple(sorted(sigma)))
    return bad


def link(faces, sigma):
    return {f - sigma for f in faces if sigma <= f}


def link_components(faces, sigma):
    lk = link(faces, sigma)
    vertices = {next(iter(f)) for f in lk if len(f) == 1}
    edges = [tuple(sorted(f)) for f in lk if len(f) == 2]
    return components(vertices, edges)


def low_codimension_link_failures(faces, d):
    return [tuple(sorted(f)) for f in sorted(faces, key=lambda x: (len(x), sorted(x)))
            if len(f) + 2 <= d and len(link_components(faces, f)) != 1]


def diameter(facets, d):
    adjacency = {f: set() for f in facets}
    for a, b in dual_edges(facets, d):
        adjacency[a].add(b)
        adjacency[b].add(a)
    answer = 0
    for start in facets:
        distance = {start: 0}
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for w in adjacency[v] - distance.keys():
                distance[w] = distance[v] + 1
                queue.append(w)
        assert len(distance) == len(facets)
        answer = max(answer, max(distance.values()))
    return answer


def fvector(faces):
    return tuple(sum(len(f) == k for f in faces) for k in range(1, 5))


def flag_completion(faces):
    vertices = {next(iter(f)) for f in faces if len(f) == 1}
    adjacency = {v: set() for v in vertices}
    for f in faces:
        if len(f) == 2:
            a, b = sorted(f)
            adjacency[a].add(b)
            adjacency[b].add(a)
    cliques = {frozenset()}
    def visit(prefix, possible):
        for v in sorted(possible):
            new = prefix | {v}
            cliques.add(new)
            visit(new, {w for w in possible & adjacency[v] if w > v})
    visit(frozenset(), vertices)
    return cliques


def cyclic_boundary(n):
    cycle = [frozenset((i, (i + 1) % n)) for i in range(n)]
    return closure({e | f for e, f in combinations(cycle, 2) if not e & f})


def relative_subdivision(original, core):
    """Direct definition: tau in A and a strict chain of non-A original faces."""
    a = {f for f in original if f <= core}
    outside = sorted(original - a, key=lambda f: (len(f), sorted(f)))
    label = {f: 4 + i for i, f in enumerate(outside)}
    old_label = {v: i for i, v in enumerate(sorted(core))}
    faces = {frozenset(old_label[v] for v in f) for f in a}
    supersets = {f: [g for g in outside if f < g] for f in outside}
    def extend(chain, last, tau):
        faces.add(frozenset(old_label[v] for v in tau) |
                  frozenset(label[f] for f in chain))
        for g in supersets[last]:
            extend(chain + [g], g, tau)
    for f in outside:
        for tau in a:
            if tau <= f:
                extend([f], f, tau)
    return faces, label, old_label


def stellar_crosscheck(original, core, labels, old_labels):
    """Independent construction by descending stellar subdivisions, including singleton renaming."""
    facets = set(facets_of(original))
    outside = sorted(labels, key=lambda f: (-len(f), sorted(f)))
    # Shift the original vertices away from all auxiliary labels during the construction.
    shift = 1000000
    facets = {frozenset(shift + v for v in f) for f in facets}
    for f in outside:
        face = frozenset(shift + v for v in f)
        new = set()
        for g in facets:
            if face <= g:
                for v in face:
                    new.add((g - {v}) | {labels[f]})
            else:
                new.add(g)
        facets = new
    return closure({frozenset(old_labels[v - shift] if v >= shift else v for v in f)
                    for f in facets})


def fidelity():
    # Source Section 8.3: Gamma_7, its minimal nonfaces, normality and diameter 3 = n-d.
    raw = ['0256', '0356', '1245', '1246', '1256', '1345', '1346', '1356']
    facets = {frozenset(map(int, s)) for s in raw}
    faces = closure(facets)
    universe = frozenset(range(7))
    minimal_nonfaces = {s for s in powerset(universe) - faces
                        if all(s - {v} in faces for v in s)}
    assert minimal_nonfaces == {frozenset(map(int, s)) for s in ['01', '04', '23', '456']}
    assert {len(f) for f in facets_of(faces)} == {4}
    assert normal_bad_faces(faces, 4) == []
    assert low_codimension_link_failures(faces, 4) == []
    assert diameter(facets, 4) == 3
    completion = flag_completion(faces)
    assert max(map(len, completion)) == 5
    print('FIDELITY Gamma_7: listed facets and minimal nonfaces agree; normal; dual diameter=3=n-d; completion dimension=4')


def check(n):
    original = cyclic_boundary(n)
    core = frozenset((0, 1, 3, 5))
    a = {f for f in original if f <= core}
    assert {f for f in a if len(f) == 2} == {frozenset(e) for e in combinations(core, 2)}
    assert {f for f in a if len(f) == 3} == {frozenset((0, 1, 3)), frozenset((0, 1, 5))}
    assert core not in original
    faces, labels, old_labels = relative_subdivision(original, core)
    assert faces == closure(facets_of(faces))
    assert faces == stellar_crosscheck(original, core, labels, old_labels)
    assert {len(f) for f in facets_of(faces)} == {4}
    assert normal_bad_faces(faces, 4) == []
    assert low_codimension_link_failures(faces, 4) == []
    # Extra manifold check: every ridge belongs to exactly two tetrahedra.
    ridge_degrees = Counter(f - {v} for f in facets_of(faces) for v in f)
    assert set(ridge_degrees.values()) == {2}
    completion = flag_completion(faces)
    mapped_core = frozenset(old_labels.values())
    assert completion == faces | powerset(mapped_core)
    added = completion - faces
    assert added == {frozenset(old_labels[v] for v in f)
                     for f in [(0, 3, 5), (1, 3, 5), (0, 1, 3, 5)]}
    assert {len(f) for f in facets_of(completion)} == {4}
    for k in (faces, completion):
        fs = facets_of(k)
        assert len(components(fs, dual_edges(fs, 4))) == 1
    sigma = frozenset(old_labels[v] for v in (3, 5))
    old_link = link(faces, sigma)
    new_link = link(completion, sigma)
    old_vertices = {next(iter(f)) for f in old_link if len(f) == 1}
    old_edges = {f for f in old_link if len(f) == 2}
    assert len(old_vertices) == len(old_edges) == 6
    assert all(sum(v in e for e in old_edges) == 2 for v in old_vertices)
    assert len(link_components(faces, sigma)) == 1
    assert new_link == old_link | powerset(frozenset(old_labels[v] for v in (0, 1)))
    assert sorted(map(len, link_components(completion, sigma))) == [2, 6]
    bad_stars = normal_bad_faces(completion, 4)
    bad_links = low_codimension_link_failures(completion, 4)
    assert tuple(sorted(sigma)) in bad_stars
    assert tuple(sorted(sigma)) in bad_links
    # G_2 for d=4: there is no integer m with 3 <= m <= d-k=2.
    assert list(range(3, 4 - 2 + 1)) == []
    if n == 7:
        assert fvector(faces) == (62, 346, 568, 284)
        assert fvector(completion) == (62, 346, 570, 285)
    print(f'N={n}: K f={fvector(original)}; relative f={fvector(faces)}; completion f={fvector(completion)}')
    print(f'N={n}: relative normal at all {len(faces)} faces; all {sum(len(f)+2 <= 4 for f in faces)} codimension>=2 links connected; every ridge degree=2')
    print(f'N={n}: stellar crosscheck equal; all maximal cliques size=4; both dual graphs connected; added faces={sorted(map(sorted, added))}')
    print(f'N={n}: completion bad facet stars={bad_stars}; bad codimension>=2 links={bad_links}; edge 35 link components=[2,6]')
    reverse = {v: tuple(sorted(f)) for f, v in labels.items()}
    print(f'N={n}: old edge 35 link auxiliary labels={sorted((v, reverse[v]) for v in old_vertices)}')


if __name__ == '__main__':
    fidelity()
    for n in (7, 8):
        check(n)
    print('ALL_CHECKS_OK')
