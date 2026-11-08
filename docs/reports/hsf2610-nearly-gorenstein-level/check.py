#!/usr/bin/env python3
"""Exact checks for Hibi--Seyed Fakhari, arXiv:2610.24438v1, Conjecture 0.2.

No floating point or numerical convex hull is used. Products of edges determine
B and D. Discrete-polymatroid augmentation, rank submodularity, and the rank
inequalities are checked independently. Canonical monomials are strict integer
solutions of these inequalities. Normality of a discrete polymatroid supplies
the semigroup interpretation. Canonical generators have degree <= n+1; see the
README for the elementary convex-combination bound used to make the tests finite.
"""
from itertools import combinations, product
from functools import lru_cache
import json


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def le(a, b):
    return all(x <= y for x, y in zip(a, b))


def box(c):
    return product(*(range(k + 1) for k in c))


def discrete(D, n):
    assert (0,) * n in D
    for x in D:
        for i in range(n):
            if x[i]:
                z = list(x); z[i] -= 1
                assert tuple(z) in D, ('downward', x, i)
    for x in D:
        for y in D:
            if sum(x) < sum(y):
                assert any(x[i] < y[i] and
                           tuple(x[k] + (k == i) for k in range(n)) in D
                           for i in range(n)), ('augmentation', x, y)


def construction(n, edges, c):
    E = {tuple(int(k == i) + int(k == j) for k in range(n))
         for i, j in edges}
    B = {(0,) * n}; delta = 0
    while True:
        nxt = {add(b, e) for b in B for e in E if le(add(b, e), c)}
        if not nxt:
            break
        B = nxt; delta += 1
    D = set().union(*(set(box(b)) for b in B))
    return delta, B, D


def check_case(n, edges, c, label):
    delta, B, D = construction(n, edges, c)
    assert delta >= 1 and len({sum(b) for b in B}) == 1
    assert all(sum(b) == 2 * delta for b in B)
    discrete(D, n)
    assert all(tuple(int(i == k) for k in range(n)) in D for i in range(n))
    subsets = [tuple(i for i in range(n) if mask >> i & 1)
               for mask in range(1 << n)]
    rho = [max(sum(x[i] for i in A) for x in D) for A in subsets]
    for a in range(1 << n):
        for b in range(1 << n):
            assert rho[a] + rho[b] >= rho[a | b] + rho[a & b]
            if a & b == a:
                assert rho[a] <= rho[b]
    bounds = tuple(rho[1 << i] for i in range(n))
    assert D == {x for x in box(bounds)
                 if all(sum(x[i] for i in A) <= r
                        for A, r in zip(subsets, rho))}

    @lru_cache(None)
    def lattice(j):
        if j < 0:
            return frozenset()
        return frozenset(x for x in box(tuple(j * b for b in bounds))
                         if all(sum(x[i] for i in A) <= j * r
                                for A, r in zip(subsets, rho)))

    @lru_cache(None)
    def interior(j):
        if j <= 0:
            return frozenset()
        return frozenset(x for x in product(*(range(1, j*b) for b in bounds))
                         if all(sum(x[i] for i in subsets[m]) < j * rho[m]
                                for m in range(1, 1 << n)))

    def in_semigroup(j, x):
        return j >= 0 and all(t >= 0 for t in x) and all(
            sum(x[i] for i in A) <= j * r for A, r in zip(subsets, rho))

    s = next(j for j in range(1, n + 2) if interior(j))
    # Independently enumerate all minimal canonical monomials; degree-one
    # subtraction tests ideal generation, not the proposed conclusion.
    generators = []
    for j in range(1, n + 2):
        for u in sorted(interior(j)):
            if not any(sub(u, a) in interior(j - 1) for a in D):
                generators.append((j, u))
    assert generators
    level = all(j == s for j, _ in generators)
    gorenstein = len(generators) == 1

    def in_inverse(d, v):
        # f * omega <= R, checked on every canonical ideal generator.
        return all(in_semigroup(d + j, add(v, u)) for j, u in generators)

    witnesses = {}
    for a in sorted(D):
        for j, u in generators:
            v = sub(a, u)
            if in_inverse(1 - j, v):
                witnesses[a] = (j, u, v)
                break
    nearly = len(witnesses) == len(D)
    assert not nearly or level, ('CONJECTURE_FAILED', label, c)
    if nearly:
        # Check the claimed exact lattice sum, not only convex hull equality.
        for N in range(s, n + 2):
            rhs = {add(u, w) for u in interior(s) for w in lattice(N - s)}
            assert interior(N) == rhs, ('lattice_sum', label, c, N)
        if not gorenstein:
            # Lemma 2: every tested inverse monomial producing degree one
            # has degree >= 1-s and the canonical factor has degree exactly s.
            for a in D:
                for j, u in generators:
                    if in_inverse(1-j, sub(a, u)):
                        assert j == s
    # Complete graphs: independently reproduce the source's classification,
    # which is stronger than its Corollary 1.7 implication.
    if len(edges) == n * (n - 1) // 2:
        cube = len(B) == 1
        cube_lengths = next(iter(B)) if cube else None
        expected_nearly = (cube and (max(cube_lengths) <= 2 or
                                    min(cube_lengths) >= 2)) or (
                                        n == 3 and c == (1, 1, 1))
        expected_level = cube or (n == 3 and min(c) == 1)
        assert nearly == expected_nearly, ('complete_nearly', c, nearly)
        assert level == expected_level, ('complete_level', c, level)
        assert not nearly or level  # Corollary 1.7
    return dict(label=label, n=n, edges=[list(e) for e in edges], c=list(c),
                delta=delta, bases=len(B), D=len(D), codegree=s,
                canonical_generators=len(generators), gorenstein=gorenstein,
                nearly_gorenstein=nearly, level=level,
                last_canonical_degree_checked=n+1)


def main():
    cases = []
    # Every labelled graph without isolated vertices, n=3, capacities in {1,2}.
    # Every such graph, n=4, for both equal and unequal capacities.
    for n in (3, 4):
        E = list(combinations(range(n), 2))
        caps = list(product((1, 2), repeat=n)) if n == 3 else [
            (1, 1, 1, 1), (1, 2, 1, 2)]
        for mask in range(1 << len(E)):
            edges = [e for k, e in enumerate(E) if mask >> k & 1]
            if any(not any(i in e for e in edges) for i in range(n)):
                continue
            for c in caps:
                cases.append(check_case(n, edges, c, f'all-labelled-{n}-{mask}'))
    named = [
        (5, [(0,1),(1,2),(2,3),(3,4),(0,4)], (1,1,1,1,1), 'C5'),
        (5, [(0,1),(1,2),(2,3),(3,4),(0,4)], (1,2,1,2,1), 'C5-unequal'),
        (4, [(0,1),(1,2),(0,2),(2,3)], (2,1,2,1), 'triangle-pendant-unequal'),
        (4, list(combinations(range(4),2)), (2,2,2,2), 'K4-s1'),
        (3, list(combinations(range(3),2)), (2,2,2), 'K3-s1'),
        (4, [(0,1),(2,3)], (1,2,1,1), 'disconnected-unequal'),
        (4, [(0,2),(0,3),(1,2),(1,3)], (3,2,2,2), 'source-Example2.10-3'),
    ]
    for args in named:
        cases.append(check_case(*args))
    anchor = cases[-1]
    assert anchor['codegree'] == 1 and anchor['nearly_gorenstein']
    for r in cases:
        print(json.dumps(r, sort_keys=True))
    print(json.dumps({'cases':len(cases),
                      'nearly_gorenstein':sum(r['nearly_gorenstein'] for r in cases),
                      'non_gorenstein_nearly':sum(r['nearly_gorenstein'] and not r['gorenstein'] for r in cases),
                      'nonlevel':sum(not r['level'] for r in cases)}, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
