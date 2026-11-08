#!/usr/bin/env python3
"""Exact checks for Spiridon, arXiv:2606.19006v1, Question 3.29.

These finite and symbolic checks support a proposed paper argument. They do not
classify all subspaces over an algebraically closed field or prove the settlement.
Only the Python standard library and SymPy exact polynomial algebra are used.
"""
from fractions import Fraction
from itertools import combinations, product
from random import Random
import sympy as sp


def rref(rows, width):
    a = [[Fraction(x) for x in row] for row in rows if any(row)]
    pivot = 0
    for col in range(width):
        j = next((j for j in range(pivot, len(a)) if a[j][col]), None)
        if j is None:
            continue
        a[pivot], a[j] = a[j], a[pivot]
        x = a[pivot][col]
        a[pivot] = [v / x for v in a[pivot]]
        for j in range(len(a)):
            if j != pivot:
                x = a[j][col]
                a[j] = [u - x * v for u, v in zip(a[j], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return tuple(tuple(row) for row in a[:pivot])


def pairs(n):
    return list(combinations(range(n), 2))


def cycle(n):
    return {tuple(sorted((i, (i + 1) % n))) for i in range(n)}


def disconnected(n, edges, support):
    if len(support) < 2:
        return False
    reached = {min(support)}
    while True:
        new = reached | {v for u in reached for v in support
                         if tuple(sorted((u, v))) in edges and u != v}
        if new == reached:
            return len(reached) != len(support)
        reached = new


def components(n, edges):
    subsets = [frozenset(c) for m in range(2, n + 1)
               for c in combinations(range(n), m)]
    ds = [s for s in subsets if disconnected(n, edges, s)]
    return sorted([s for s in ds if not any(s < t for t in ds)],
                  key=lambda s: tuple(sorted(s)))


def coordinate_k(n, edges):
    return rref([[int(p == e) for p in pairs(n)] for e in sorted(edges)], len(pairs(n)))


def contracts(n, rows, a):
    result = []
    for row in rows:
        b = [Fraction(0)] * n
        for c, (i, j) in zip(row, pairs(n)):
            b[j] += c * a[i]
            b[i] -= c * a[j]
        result.append(b)
    return result


def resonant(n, rows, a):
    # a is always in the contraction kernel. A nonproportional b exists iff
    # nullity is at least two (with the source's explicit convention at a=0).
    return not any(a) or len(rref(contracts(n, rows, a), n)) <= n - 2


def locus_signature(n, edges):
    cs = components(n, edges)
    return tuple(any(frozenset(s) <= t for t in cs)
                 for m in range(1, n + 1) for s in combinations(range(n), m))


def partitions(items):
    if not items:
        yield []
        return
    first, *rest = items
    for p in partitions(rest):
        yield [[first]] + p
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]


def main():
    c4 = cycle(4)
    c5 = cycle(5)
    u4 = components(4, c4)
    u5 = components(5, c5)
    assert u4 == [frozenset((0, 2)), frozenset((1, 3))]
    expected = [(0, 1, 3), (0, 2, 3), (0, 2, 4), (1, 2, 4), (1, 3, 4)]
    assert [tuple(sorted(s)) for s in u5] == expected
    print('C4 maximal disconnected supports: 13, 24')
    print('C5 maximal disconnected supports: 124, 134, 135, 235, 245')
    counts = {p: sum(set(p) <= u for u in u5) for p in pairs(5)}
    assert all(counts[p] == (1 if p in c5 else 2) for p in counts)
    print('C5 pair coverage: every edge once; every nonedge twice; all 10 covered')

    # Reproduce the source graph-resonance formula from contraction ranks, on
    # every support, for every labelled graph on 4 or 5 vertices. The rank is
    # support-invariant by diagonal rescaling of the nonzero coordinates.
    total_graphs = total_supports = 0
    equal_counts = {}
    for n in (4, 5):
        ps = pairs(n)
        target = locus_signature(n, cycle(n))
        equal = []
        for mask in range(1 << len(ps)):
            edges = {p for i, p in enumerate(ps) if mask >> i & 1}
            sig = locus_signature(n, edges)
            rows = coordinate_k(n, edges)
            k = 0
            for m in range(1, n + 1):
                for support in combinations(range(n), m):
                    a = [int(i in support) for i in range(n)]
                    assert resonant(n, rows, a) == sig[k]
                    k += 1
                    total_supports += 1
            if sig == target:
                equal.append(edges)
            total_graphs += 1
        assert equal == [cycle(n)]
        equal_counts[n] = len(equal)
    print(f'Graph-resonance formula: {total_graphs} graphs, {total_supports} nonempty supports checked')
    print('Coordinate-subspace equality: C4 unique among 64; C5 unique among 1024')

    # Source fidelity: Proposition 3.28 on deterministic random rational
    # subspaces, plus exact alternative bases of K_C4. For each unequal row
    # space a rational vector witnesses unequal resonance; no sampled equality
    # is inferred merely from agreement on a grid.
    rng = Random(260619006)
    target_rows = coordinate_k(4, c4)
    grid = [a for a in product((-1, 0, 1), repeat=4)
            if any(a) and next(x for x in a if x) == 1]
    sampled = equal_samples = different_samples = 0
    for dim in range(7):
        for _ in range(20):
            rows = rref([[rng.randrange(-3, 4) for _ in range(6)] for _ in range(dim)], 6)
            same = rows == target_rows
            witness = next((a for a in grid
                            if resonant(4, rows, a) != resonant(4, target_rows, a)), None)
            assert same or witness is not None
            equal_samples += same
            different_samples += not same
            sampled += 1
    for _ in range(20):
        rows = [list(r) for r in target_rows]
        for _ in range(12):
            i, j = rng.sample(range(4), 2)
            c = rng.choice((-2, -1, 1, 2))
            rows[i] = [x + c * y for x, y in zip(rows[i], rows[j])]
        assert rref(rows, 6) == target_rows
        assert all(resonant(4, rows, a) == resonant(4, target_rows, a) for a in grid)
        equal_samples += 1
        sampled += 1
    print(f'Proposition 3.28 sample: {sampled} rational subspace presentations; '
          f'{equal_samples} exact equalities; {different_samples} witnessed inequalities')

    # The identity used in the source proof of Proposition 3.28.
    a12, a14, a23, a34 = sp.symbols('a12 a14 a23 a34')
    a = [a14, 1, a34, 0]
    b = [-a12, 0, a23, 1]
    wedge = [sp.expand(a[i] * b[j] - a[j] * b[i]) for i, j in pairs(4)]
    assert wedge == [a12, a12*a34 + a14*a23, a14, a23, 1, a34]
    print('Proposition 3.28 decomposable-bivector identity: exact polynomial identity checked')

    # The two omitted coordinate equations of each U_i have nonzero determinant
    # on A=span(p,q). All entries are integers, so this holds in characteristic 0.
    determinants = []
    for u in u5:
        r, t = sorted(set(range(5)) - u)
        d = (t + 1) - (r + 1)
        assert d != 0
        determinants.append(d)
    print(f'C5 pencil A intersection determinants: {determinants}; all nonzero')

    # Exact universal algebra for the dimension-four obstruction. The original
    # contraction matrix M has rank four iff B does at each nonzero (s,t).
    s, t = sp.symbols('s t')
    uv = sp.Matrix(sp.symbols('u0:4'))
    vv = [sp.Matrix([sp.Symbol(f'v{j}_{i}') for i in range(4)]) for j in range(3)]
    ww = [sp.Matrix([sp.Symbol(f'w{j}_{i}') for i in range(4)]) for j in range(3)]
    B = sp.Matrix.hstack(uv, *[s*v + t*w for v, w in zip(vv, ww)])
    C = sp.Matrix([[-t, s, 0, 0, 0], [0, 0, 1, 0, 0],
                   [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    M = sp.Matrix.hstack(-t*uv, s*uv, *[s*v + t*w for v, w in zip(vv, ww)])
    assert M == B*C
    f = sp.expand(B.det(method='domain-ge'))
    assert all(sum(m) == 3 for m, _ in sp.Poly(f, s, t).terms())
    assert sp.expand(f.subs(t, 0) - s**3 * sp.Matrix.hstack(uv, *vv).det()) == 0
    assert sp.expand(f.subs(s, 0) - t**3 * sp.Matrix.hstack(uv, *ww).det()) == 0
    # The maximal minors of C are s,-t,0,0,0. Multiplicativity
    # gives the maximal minors of M as s*f,-t*f,0,0,0.
    assert [C[:, [j for j in range(5) if j != i]].det() for i in range(5)] == [s, -t, 0, 0, 0]
    print('Dimension-four obstruction: contraction factorization and homogeneous cubic checked')

    # Every arrangement of zero/nonzero W_i and coincidences of nonzero lines.
    # Each nonzero pair coordinate must lie in every W_i containing that pair.
    patterns = 0
    for mask in range(32):
        active = [i for i in range(5) if mask >> i & 1]
        for blocks in partitions(active):
            classes = {i: j for j, block in enumerate(blocks) for i in block}
            allowed = []
            for p in pairs(5):
                owners = [i for i, u in enumerate(u5) if set(p) <= u]
                if all(i in classes for i in owners) and len({classes[i] for i in owners}) == 1:
                    allowed.append(p)
            # The span of all W_i has dimension at most the number of line classes.
            if len(blocks) == 5:
                assert set(allowed) == c5
            if set(allowed) - c5:
                assert len(blocks) <= 4
            patterns += 1
    assert patterns == 203
    print(f'C5 gluing: {patterns} zero/line-coincidence patterns; any surviving nonedge forces rank <=4')
    print('LIMIT: finite/symbolic checks only; arbitrary-K rigidity and algebraic-closure root step remain paper mathematics')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
