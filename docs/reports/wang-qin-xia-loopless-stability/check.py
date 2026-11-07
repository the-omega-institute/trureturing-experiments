#!/usr/bin/env python3
"""Exact finite checks for Wang--Qin--Xia, arXiv:2509.26170v2, Conjecture 1.3.

Python standard library only. These finite checks are not a proof for all m.
Relations include diagonal entries whenever the source permits loops.
"""
from collections import deque
from functools import lru_cache
from itertools import permutations, product
from math import factorial


def matrix(vertices, adjacent):
    return tuple(tuple(bool(adjacent(u, v)) for v in vertices) for u in vertices)


def neighbors(a):
    return tuple(frozenset(j for j, edge in enumerate(row) if edge) for row in a)


def connected(a):
    if not a:
        return False
    seen = {0}
    queue = deque([0])
    while queue:
        for v, edge in enumerate(a[queue.popleft()]):
            if edge and v not in seen:
                seen.add(v)
                queue.append(v)
    return len(seen) == len(a)


def bipartite(a):
    color = {}
    for start in range(len(a)):
        if start in color:
            continue
        color[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v, edge in enumerate(a[u]):
                if not edge:
                    continue
                if v in color:
                    if color[u] == color[v]:
                        return False
                else:
                    color[v] = 1 - color[u]
                    queue.append(v)
    return True


def twin_free(a):
    return len(set(neighbors(a))) == len(a)


def direct(a, b):
    vs = tuple(product(range(len(a)), range(len(b))))
    return matrix(vs, lambda u, v: a[u[0]][v[0]] and b[u[1]][v[1]])


def maximum_independent_set(a):
    assert all(not a[i][i] for i in range(len(a)))
    masks = tuple(sum(1 << j for j, edge in enumerate(row) if edge) for row in a)

    @lru_cache(None)
    def solve(mask):
        if not mask:
            return 0
        remaining = tuple(i for i in range(len(a)) if mask & (1 << i))
        v = max(remaining, key=lambda i: (mask & masks[i]).bit_count())
        without = mask & ~(1 << v)
        return max(solve(without), 1 + solve(without & ~masks[v]))

    return solve((1 << len(a)) - 1)


def two_fold_automorphisms(a):
    """Exhaust ALL degree-preserving row permutations, derive each column map.

    For a twin-free symmetric matrix each column neighborhood is unique. If
    alpha is a row permutation, beta(v) must have neighborhood alpha(N(v)).
    Therefore there is at most one beta for each tested alpha. Row degrees are
    necessary invariants of every two-fold automorphism. No heuristic pruning.
    """
    ns = neighbors(a)
    assert twin_free(a)
    lookup = {n: i for i, n in enumerate(ns)}
    degree_groups = {}
    for i, n in enumerate(ns):
        degree_groups.setdefault(len(n), []).append(i)
    groups = tuple(tuple(v) for v in degree_groups.values())
    tested = 0
    result = []
    for choices in product(*(permutations(g) for g in groups)):
        alpha = [None] * len(a)
        for group, image in zip(groups, choices):
            for u, v in zip(group, image):
                alpha[u] = v
        tested += 1
        image_columns = tuple(frozenset(alpha[u] for u in n) for n in ns)
        if any(n not in lookup for n in image_columns):
            continue
        beta = tuple(lookup[n] for n in image_columns)
        assert len(set(beta)) == len(a)
        assert all(a[u][v] == a[alpha[u]][beta[v]]
                   for u in range(len(a)) for v in range(len(a)))
        result.append((tuple(alpha), beta))
    assert tested == factorial_product(len(g) for g in groups)
    return tested, tuple(result)


def factorial_product(sizes):
    answer = 1
    for n in sizes:
        answer *= factorial(n)
    return answer


def isomorphic(a, b):
    if len(a) != len(b):
        return False
    if sorted(map(sum, a)) != sorted(map(sum, b)):
        return False
    for p in permutations(range(len(a))):
        if all(a[i][j] == b[p[i]][p[j]]
               for i in range(len(a)) for j in range(len(a))):
            return True
    return False


def symmetric_relations(n):
    positions = tuple((i, j) for i in range(n) for j in range(i, n))
    for bits in product((False, True), repeat=len(positions)):
        a = [[False] * n for _ in range(n)]
        for (i, j), bit in zip(positions, bits):
            a[i][j] = a[j][i] = bit
        yield tuple(map(tuple, a))


def check_source_cycle_skeletons():
    # WQX v2, Section 2, equation S(C_n)=B(C_n): adjacent indices differ by 2.
    for n in (3, 5, 6, 7, 8, 9, 10):
        a = matrix(tuple(range(n)), lambda i, j: (i - j) % n in (1, n - 1))
        ns = neighbors(a)
        b = matrix(tuple(range(n)), lambda u, v: u != v and bool(ns[u] & ns[v]))
        expected = matrix(tuple(range(n)), lambda u, v: (u - v) % n in (2, n - 2))
        assert b == expected
        dispensable = []
        for u in range(n):
            for v in range(u + 1, n):
                if not b[u][v]:
                    continue
                for w in range(n):
                    first = ((ns[u] & ns[v]) < (ns[u] & ns[w])
                             or ns[u] < ns[w] < ns[v])
                    second = ((ns[v] & ns[u]) < (ns[v] & ns[w])
                              or ns[v] < ns[w] < ns[u])
                    if first and second:
                        dispensable.append((u, v))
        assert not dispensable
        assert all(sum(row) == 2 for row in b)
        assert connected(b) == (n % 2 == 1)
        print(f'FIDELITY WQX cycle n={n}: B(C_n)=S(C_n), no dispensable edges')


def relation_r(a, b):
    return (a, b) in ((0, 1), (1, 0), (1, 1))


def mixer(v):
    (t, a, b), (side, i) = v
    if side == 0:
        return ((t, a, i), (0, b))
    return ((t, i, b), (1, a))


def main():
    check_source_cycle_skeletons()
    bits = (0, 1)
    delta_vertices = tuple(product(bits, bits))
    delta = matrix(delta_vertices,
                   lambda u, v: relation_r(u[0], v[1]) and relation_r(v[0], u[1]))
    assert neighbors(delta) == (frozenset((3,)), frozenset((1, 3)),
                                frozenset((2, 3)), frozenset((0, 1, 2, 3)))
    assert sum(delta[i][i] for i in range(4)) == 3
    tested, tf = two_fold_automorphisms(delta)
    assert tested == 2 and len(tf) == 2 and all(a == b for a, b in tf)
    assert connected(delta) and not bipartite(delta) and connected(direct(delta, ((0, 1), (1, 0))))
    assert all(x * y != 3 for x in range(3) for y in range(3))
    print('FIDELITY Wang--Gao Section 5: Delta has 3 loops, 2 diagonal TF automorphisms; loop-count excludes 2-by-2 factors')

    omega_vertices = tuple(product(bits, bits))  # (side, label): x_i or y_i
    omega = matrix(omega_vertices,
                   lambda u, v: u[0] != v[0] and relation_r(u[1], v[1]))
    path4 = matrix(tuple(range(4)), lambda i, j: abs(i - j) == 1)
    assert isomorphic(omega, path4) and connected(omega) and bipartite(omega) and twin_free(omega)

    # Enumerate every possible finite loop-allowed factor of P4: its order is 2 or 4.
    # In each decomposition both orders are fixed by the four-element vertex set.
    p4_factors = []
    for d in (2, 4):
        for f in symmetric_relations(d):
            if any(isomorphic(direct(j, f), path4) for j in symmetric_relations(4 // d)):
                p4_factors.append(f)
                assert connected(f) and twin_free(f)
                if not bipartite(f):
                    assert d == 2 and f[0][1] and sum(f[i][i] for i in range(2)) == 1
    eligible = [f for f in p4_factors if not bipartite(f)]
    assert len(eligible) == 2  # the two labelled versions of the one-loop edge
    print(f'P4 full loop-allowed factor screen: {len(p4_factors)} labelled factors; 2 nonbipartite factors, both one-loop edges')

    for m in range(3, 9):
        vs = tuple(product(range(m), bits, bits))

        def adj(u, v):
            return (u[0] != v[0] and relation_r(u[1], v[2])
                    and relation_r(v[1], u[2]))

        a = matrix(vs, adj)
        assert a == direct(matrix(tuple(range(m)), lambda i, j: i != j), delta)
        assert all(not a[i][i] for i in range(len(a)))
        assert connected(a) and twin_free(a) and not bipartite(a)
        assert all(adj((i, 1, 1), (j, 1, 1)) for i in range(3) for j in range(3) if i != j)
        for u, row in zip(vs, a):
            i, x, y = u
            expected = {v for v in vs if v[0] != i and (
                (x, y) == (1, 1) or
                ((x, y) == (0, 0) and v[1:] == (1, 1)) or
                ((x, y) == (0, 1) and v[2] == 1) or
                ((x, y) == (1, 0) and v[1] == 1))}
            assert {v for v, edge in zip(vs, row) if edge} == expected
        alpha = maximum_independent_set(a)
        assert alpha == m + 2 < 2 * m
        # Every remaining eligible common factor would create a half-size independent slice.
        for f in eligible:
            unlooped = next(i for i in range(2) if not f[i][i])
            assert all(not f[unlooped][unlooped] for _ in range(2 * m))
        assert connected(direct(a, ((0, 1), (1, 0))))
        pv = tuple(product(vs, omega_vertices))
        images = tuple(map(mixer, pv))
        assert len(set(images)) == len(pv)
        assert all(mixer(mixer(v)) == v for v in pv)
        assert all((adj(u[0], v[0]) and omega[omega_vertices.index(u[1])][omega_vertices.index(v[1])])
                   == (adj(mixer(u)[0], mixer(v)[0])
                       and omega[omega_vertices.index(mixer(u)[1])][omega_vertices.index(mixer(v)[1])])
                   for u in pv for v in pv)
        fixed = (0, 0, 0)
        assert mixer((fixed, (0, 0)))[0] != mixer((fixed, (0, 1)))[0]
        print(f'FAMILY m={m}: vertices={4*m}, independent_number={alpha}, connected/twin-free/triangle/cover OK, mixer ordered pairs={len(pv)**2}, noncomponentwise OK')
        if m == 3:
            tested, tf = two_fold_automorphisms(a)
            diagonal = sum(x == y for x, y in tf)
            assert tested == 25920 and len(tf) == diagonal == 12
            assert diagonal == 2 * factorial(m)
            print(f'STABILITY Gamma_3: degree-preserving row permutations={tested}, TF={len(tf)}, Aut={diagonal}, nondiagonal=0')
    print('SCOPE: finite family checks m=3..8; exhaustive stability check m=3; universal family remains a proposed argument, with no Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
