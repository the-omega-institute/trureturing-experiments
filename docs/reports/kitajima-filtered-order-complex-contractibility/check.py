#!/usr/bin/env python3
"""Exact finite checks for Kitajima, arXiv:2606.15241v2, Remark 2.13.

No numerical topology package or floating-point arithmetic is used. Homology
ranks are over F_2; nonzero H_2(F_2) obstructs contractibility. Order-map and
simplicial-collapse certificates supply independent contractibility checks.
The finite family checks do not constitute a proof for every t.
"""
from itertools import combinations
from collections import Counter


def subfaces(face):
    face = tuple(sorted(face))
    return {frozenset(c) for k in range(len(face) + 1) for c in combinations(face, k)}


def complex_from_facets(facets):
    return set().union(*(subfaces(f) for f in facets))


def order_complex(vertices, lt):
    result = {frozenset()}
    def extend(chain):
        for v in sorted(vertices):
            if not chain or lt[chain[-1]][v]:
                c = chain + (v,)
                result.add(frozenset(c))
                extend(c)
    extend(())
    return result


def filtered(simplices, rank, k):
    return {s for s in simplices if not s or max(rank[v] for v in s) - min(rank[v] for v in s) <= k}


def cone(apex, simplices):
    return simplices | {s | {apex} for s in simplices}


def gf2_rank(columns):
    pivots = {}
    for col in columns:
        while col:
            p = col.bit_length() - 1
            if p in pivots:
                col ^= pivots[p]
            else:
                pivots[p] = col
                break
    return len(pivots)


def homology(simplices):
    dim = max(map(len, simplices)) - 1
    faces = [sorted((s for s in simplices if len(s) == d + 1), key=lambda s: tuple(sorted(s))) for d in range(dim + 1)]
    boundary = [[]]
    for d in range(1, dim + 1):
        index = {s: i for i, s in enumerate(faces[d - 1])}
        boundary.append([sum(1 << index[s - {v}] for v in s) for s in faces[d]])
    for d in range(2, dim + 1):
        for col in boundary[d]:
            composite = 0
            while col:
                bit = col & -col
                composite ^= boundary[d - 1][bit.bit_length() - 1]
                col ^= bit
            assert composite == 0, 'boundary squared is not zero'
    ranks = [0] + [gf2_rank(boundary[d]) for d in range(1, dim + 1)] + [0]
    f = [len(fs) for fs in faces]
    betti = [f[d] - ranks[d] - ranks[d + 1] for d in range(dim + 1)]
    return f, betti


def collapse_to_point(simplices):
    """Return an explicitly checked sequence of elementary free-face collapses."""
    remaining = set(simplices)
    steps = 0
    while len(remaining) > 2:
        chosen = None
        for s in sorted((s for s in remaining if s), key=lambda s: (-len(s), tuple(sorted(s)))):
            strict_cofaces = [t for t in remaining if s < t]
            if len(strict_cofaces) == 1 and len(strict_cofaces[0]) == len(s) + 1:
                chosen = (s, strict_cofaces[0])
                break
        if chosen is None:
            return None
        remaining.difference_update(chosen)
        steps += 1
    assert frozenset() in remaining
    assert len(remaining) == 2 and max(map(len, remaining)) == 1
    return steps


def check_order_map(domain, target, function, below=False, above=False):
    for s in domain:
        assert function(s) in target
        if below:
            assert function(s) <= s
        if above:
            assert s <= function(s)
        for t in domain:
            if s <= t:
                assert function(s) <= function(t)


def witness(t):
    ds = set(range(3, 3 + t))
    K = complex_from_facets([{0, 1, d} for d in ds] + [{1, 2, d} for d in ds])
    L = subfaces({0, 1})
    B = complex_from_facets([{0, 1, d} for d in ds])
    C = K - {frozenset()}
    U, V = C - L, B - {frozenset()}
    J = U & V
    faces = sorted(C, key=lambda s: (len(s), tuple(sorted(s))))
    labels = faces + ['a', 'x', 'b', 'y']
    index = {s: i for i, s in enumerate(labels)}
    a, x, b, y = (index[s] for s in ['a', 'x', 'b', 'y'])
    n = len(labels)
    rank = [len(s) for s in faces] + [0, 0, 4, 4]
    lt = [[False] * n for _ in labels]
    for s in C:
        i = index[s]
        lt[a][i] = lt[i][b] = True
        if s in U:
            lt[x][i] = True
        if s in V:
            lt[i][y] = True
        for q in C:
            if s < q:
                lt[i][index[q]] = True
    for z in range(n):
        for i in range(n):
            if lt[i][z]:
                for j in range(n):
                    lt[i][j] = lt[i][j] or lt[z][j]
    assert not any(lt[i][i] for i in range(n))
    for i in range(n):
        for j in range(n):
            if lt[i][j]:
                assert rank[i] < rank[j]
                for z in range(n):
                    assert not lt[j][z] or lt[i][z]
                if not any(lt[i][z] and lt[z][j] for z in range(n)):
                    assert rank[j] == rank[i] + 1
    assert {i for i in range(n) if not any(lt[j][i] for j in range(n))} == {a, x}
    assert {i for i in range(n) if not any(lt[i][j] for j in range(n))} == {b, y}
    assert all(lt[i][j] for i in [a, x] for j in [b, y])
    assert n == 9 + 6 * t
    assert [Counter(rank)[i] for i in range(5)] == [2, 3 + t, 2 + 3 * t, 2 * t, 2]
    # Explicit comparable order maps: C,V -> faces containing the cone apex 1.
    for domain in [C, V]:
        check_order_map(domain, domain, lambda s: s | {1}, above=True)
        assert all({1} <= (s | {1}) for s in domain)
    W = {s & ({2} | ds) for s in U}
    check_order_map(U, W, lambda s: s & ({2} | ds), below=True)
    check_order_map(W, W, lambda s: s | {2}, above=True)
    assert all({2} <= (s | {2}) for s in W)
    # J is exactly t components, each with its distinct minimum {d}.
    assert J == {frozenset({d}) | s for d in ds for s in L}
    components = [{s for s in J if d in s} for d in sorted(ds)]
    for comp, d in zip(components, sorted(ds)):
        assert frozenset({d}) in comp and all({d} <= s for s in comp)
    assert all(not (s <= q or q <= s) for i, comp in enumerate(components) for other in components[i + 1:] for s in comp for q in other)
    def oc(domain):
        return order_complex({index[s] for s in domain}, lt)
    cC, cU, cV, cJ = (oc(domain) for domain in [C, U, V, J])
    A = cone(a, cC) | cone(x, cU)
    E = cone(b, cC) | cone(y, cV)
    T = cone(a, cV) | cone(x, cJ)
    full = order_complex(set(range(n)), lt)
    stage3 = filtered(full, rank, 3)
    assert stage3 == A | E and A & E == cC
    assert cone(a, cC) & cone(x, cU) == cU
    assert cone(b, cC) & cone(y, cV) == cV
    assert cone(a, cV) & cone(x, cJ) == cJ
    assert full == cone(b, A) | cone(y, T)
    assert cone(b, A) & cone(y, T) == T
    assert {s for s in full if all(lt[v][y] for v in s)} == T
    # Independent exact verification of Proposition 2.14, including k=0.
    for k in range(5):
        windows = set().union(*({s for s in full if all(i <= rank[v] <= i + k for v in s)} for i in range(5 - k)))
        assert filtered(full, rank, k) == windows
    f3, h3 = homology(stage3)
    f4, h4 = homology(full)
    assert h3 == [1, 0, 0, 0]
    assert h4 == [1, 0, t - 1, 0, 0]
    collapse = collapse_to_point(stage3)
    assert collapse is not None, 'no elementary collapse certificate was found'
    assert sum((-1) ** i * f3[i] for i in range(len(f3))) == 1
    assert sum((-1) ** i * f4[i] for i in range(len(f4))) == t
    if t == 2:
        assert f3 == [21, 99, 151, 72]
        assert f4 == [21, 103, 201, 169, 52]
    print(f'FAMILY t={t} elements={n} rank_sizes={[Counter(rank)[i] for i in range(5)]}')
    print(f'  Delta^3 f={f3} Betti_F2={h3} elementary_collapses={collapse}')
    print(f'  Delta^4 f={f4} Betti_F2={h4}')
    print('  order_maps, components, covers, extrema, rank_windows, cone_intersections: OK')


def main():
    # Example 2.12 (and the explicit range in Remark 3.5): chains, 1<=k<=n.
    # The isolated k=0 stage is discrete, not contractible for n>0.
    for n in range(1, 7):
        lt = [[i < j for j in range(n + 1)] for i in range(n + 1)]
        full = order_complex(set(range(n + 1)), lt)
        for k in range(1, n + 1):
            stage = filtered(full, list(range(n + 1)), k)
            assert homology(stage)[1] == [1] + [0] * k
            assert collapse_to_point(stage) is not None
    print('FIDELITY Example 2.12: chains n=1..6, all 1<=k<=n collapse to a point: OK')
    for t in range(2, 7):
        witness(t)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
