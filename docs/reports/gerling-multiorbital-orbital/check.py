#!/usr/bin/env python3
"""Exact checks for Gerling arXiv:2609.01866v4's converse question.

Standard library only. Polynomials are sparse dictionaries (x exponent, y
exponent) -> Fraction. No floating point or external symbolic engine is used.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb


def clean(p):
    return {k: Q(v) for k, v in p.items() if v}


def add(*ps):
    out = defaultdict(Q)
    for p in ps:
        for k, v in p.items():
            out[k] += v
    return clean(out)


def scale(c, p):
    return clean({k: Q(c) * v for k, v in p.items()})


def evaluate(p, x, y):
    return sum((c * x**i * y**j for (i, j), c in p.items()), Q(0))


def star_poly(m):
    """(x-y)x^m + y(x-1)^m, expanded over Q exactly."""
    return add({(m + 1, 0): 1, (m, 1): -1},
               {(k, 1): comb(m, k) * (-1)**(m-k) for k in range(m+1)})


def graph_poly(n, edges):
    """Inclusion-exclusion of equal restricted endpoint-color events.

    A selected connected component incident to an edge has y choices, and
    an isolated vertex has x choices. Loops are retained; parallel edges
    are coalesced. This routine is used only for the small fidelity cases.
    """
    edges = sorted({tuple(sorted(e)) for e in edges})
    terms = defaultdict(Q)
    for mask in range(1 << len(edges)):
        parent = list(range(n))
        active = set()

        def root(u):
            while u != parent[u]:
                u = parent[u]
            return u

        for i, (u, v) in enumerate(edges):
            if mask >> i & 1:
                active.update((u, v))
                parent[root(u)] = root(v)
        components = len({root(u) for u in active})
        terms[n-len(active), components] += (-1)**mask.bit_count()
    return clean(terms)


def direct_count(n, edges, x, y):
    assert 0 <= y <= x
    return sum(all(not (f[u] == f[v] and f[u] < y) for u, v in edges)
               for f in product(range(x), repeat=n))


def cycles(p):
    unseen = set(range(len(p)))
    out = []
    while unseen:
        start = min(unseen)
        orbit = []
        u = start
        while u not in orbit:
            assert u in unseen
            unseen.remove(u)
            orbit.append(u)
            u = p[u]
        assert u == start
        out.append(tuple(orbit))
    return out


def quotient(n, edges, p):
    cs = cycles(p)
    labels = {v: i for i, c in enumerate(cs) for v in c}
    qe = {tuple(sorted((labels[u], labels[v]))) for u, v in edges}
    return len(cs), qe


def source_fidelity():
    p3 = [(0, 1), (1, 2)]
    k2i = [(0, 1)]
    e = (0, 1, 2)
    quotients = []
    for edges, g in [(p3, (2, 1, 0)), (k2i, (1, 0, 2))]:
        polys = [graph_poly(*quotient(3, edges, p)) for p in (e, g)]
        for p in (e, g):
            n, qe = quotient(3, edges, p)
            for x in range(5):
                for y in range(x+1):
                    assert evaluate(graph_poly(n, qe), x, y) == direct_count(n, qe, x, y)
        quotients.append(polys)
    op = [scale(Q(1, 2), add(*ps)) for ps in quotients]
    assert op[0] == op[1] == clean({(3, 0): Q(1, 2), (2, 0): Q(1, 2), (1, 1): -1})
    f = [add(scale(Q(3, 2), ps[0]), scale(Q(1, 2), ps[1])) for ps in quotients]
    assert add(f[0], scale(-1, f[1])) == clean({(0, 1): 1, (1, 1): -1})
    for m in range(6):
        edges = [(0, v) for v in range(1, m+1)]
        assert graph_poly(m+1, edges) == star_poly(m)
        for x in range(5):
            for y in range(x+1):
                assert evaluate(star_poly(m), x, y) == direct_count(m+1, edges, x, y)
    print('SOURCE_FIDELITY: Section 3.3, OP(P3,C2)=OP(K2+K1,C2)=(x^3+x^2)/2-xy')
    print('SOURCE_FIDELITY: F(P3,C2)-F(K2+K1,C2)=y(1-x), including the loop quotient')
    print('STAR_FORMULA: exact expansion and literal coloring counts for m=0..5, 0<=y<=x<=4')


G = tuple(product(range(2), range(2), range(3)))
ONE = (0, 0, 0)


def plus(g, h):
    return ((g[0]+h[0]) % 2, (g[1]+h[1]) % 2, (g[2]+h[2]) % 3)


def group_and_subgroups():
    assert all(plus(ONE, g) == g == plus(g, ONE) for g in G)
    assert all(plus(plus(g, h), k) == plus(g, plus(h, k)) for g, h, k in product(G, repeat=3))
    assert all(any(plus(g, h) == ONE for h in G) for g in G)
    subs = []
    for mask in range(1 << len(G)):
        H = frozenset(g for i, g in enumerate(G) if mask >> i & 1)
        if ONE in H and all(plus(g, h) in H for g in H for h in H):
            subs.append(H)
    v4subs = [frozenset([ONE[:2]]),
              frozenset([(0, 0), (1, 0)]), frozenset([(0, 0), (0, 1)]),
              frozenset([(0, 0), (1, 1)]), frozenset(product(range(2), repeat=2))]
    split = {frozenset((u, v, t) for u, v in H for t in K)
             for H in v4subs for K in [(0,), (0, 1, 2)]}
    assert len(subs) == 10 and set(subs) == split
    weights = {g: sum((Q(1, len(H)) for H in subs if g in H), Q(0)) for g in G}
    assert sum(weights.values()) == len(subs)
    print('GROUP: all 1728 associativity triples, identities and inverses checked')
    print('SUBGROUPS: all 4096 subsets tested; exactly 10, each a V4/C3 product')
    return subs, weights


def leaves(spec, r=0):
    out = []
    copies_by_kind = dict(spec)
    copies_by_kind['I'] = copies_by_kind.get('I', 0) + r
    for kind, copies in copies_by_kind.items():
        qs = G if kind == 'R' else range({'T': 3, 'A': 2, 'B': 2, 'C': 2, 'I': 1}[kind])
        out.extend((kind, copy, q) for copy in range(copies) for q in qs)
    return out


SPEC1 = [('T', 4), ('A', 9), ('B', 9), ('C', 5)]
SPEC2 = [('R', 1), ('A', 8), ('B', 8), ('I', 14)]


def act(g, leaf):
    kind, copy, q = leaf
    u, v, t = g
    q1 = {'T': lambda: (q+t) % 3, 'A': lambda: (q+v) % 2,
          'B': lambda: (q+u) % 2, 'C': lambda: (q+u+v) % 2,
          'R': lambda: plus(g, q), 'I': lambda: q}[kind]()
    return kind, copy, q1


def star_action(spec, r):
    ll = leaves(spec, r)
    assert len(ll) == len(set(ll)) == 58+r
    idx = {v: i+1 for i, v in enumerate(ll)}
    perms = {g: (0,) + tuple(idx[act(g, v)] for v in ll) for g in G}
    n = len(ll)+1
    assert all(set(p) == set(range(n)) for p in perms.values())
    assert perms[ONE] == tuple(range(n))
    assert all(perms[plus(g, h)] == tuple(perms[g][perms[h][i]] for i in range(n))
               for g, h in product(G, repeat=2))
    assert len(set(perms.values())) == len(G)  # full faithfulness, not just nontriviality
    edges = {(0, i) for i in range(1, n)}
    assert all({tuple(sorted((p[u], p[v]))) for u, v in edges} == edges for p in perms.values())
    counts = {}
    for g, p in perms.items():
        nq, qe = quotient(n, edges, p)
        assert all(u != v for u, v in qe)  # no collapsed center/leaf edge, hence no loop
        assert qe == {(0, i) for i in range(1, nq)}
        counts[g] = nq-1
    return counts


def profile(counts, weights):
    out = defaultdict(Q)
    for g, m in counts.items():
        out[m] += weights[g]
    return clean(out)


def poly_from_profile(p):
    return add(*(scale(c, star_poly(m)) for m, c in p.items()))


def subgroup_average(counts, subs):
    return add(*(scale(Q(1, len(H)), add(*(star_poly(counts[h]) for h in H))) for H in subs))


def main():
    source_fidelity()
    subs, weights = group_and_subgroups()
    c1, c2 = star_action(SPEC1, 0), star_action(SPEC2, 0)
    classes = [('1', [ONE], 44, (58, 58)),
               ('z,z^2', [(0, 0, 1), (0, 0, 2)], 11, (50, 50)),
               ('a,b', [(1, 0, 0), (0, 1, 0)], 12, (44, 44)),
               ('c', [(1, 1, 0)], 12, (40, 36)),
               ('az,az^2,bz,bz^2', [(u, v, t) for u, v in [(1, 0), (0, 1)] for t in [1, 2]], 3, (36, 40)),
               ('cz,cz^2', [(1, 1, 1), (1, 1, 2)], 3, (32, 32))]
    assert {g for _, cl, _, _ in classes for g in cl} == set(G)
    print('ACTIONS: each has 58 leaves plus a fixed center, all 144 composition pairs checked')
    print('ACTIONS: faithful, connected K_(1,58), every quotient is a loop-free star')
    print('TABLE: class | multiplicity | individual 12*c_G | leaf cycles (L1,L2)')
    for name, cl, w, pair in classes:
        assert all(12*weights[g] == w and (c1[g], c2[g]) == pair for g in cl)
        print(name, len(cl), w, pair, sep=' | ')
    fp1, fp2 = profile(c1, weights), profile(c2, weights)
    expected12 = {58: 44, 50: 22, 44: 24, 40: 12, 36: 12, 32: 6}
    assert fp1 == fp2 and scale(12, fp1) == clean(expected12)
    f1, f2 = subgroup_average(c1, subs), subgroup_average(c2, subs)
    assert f1 == f2 == poly_from_profile(fp1)
    uniform = {g: Q(1, len(G)) for g in G}
    op1, op2 = poly_from_profile(profile(c1, uniform)), poly_from_profile(profile(c2, uniform))
    delta = add(op2, scale(-1, op1))
    assert delta == scale(Q(1, 4), add(star_poly(40), scale(-1, star_poly(36))))
    val1, val2 = evaluate(op1, 2, 1), evaluate(op2, 2, 1)
    assert val1 == 24209895276412929 and val2 == 24210152974450689
    assert val2-val1 == 257698037760 == 15*2**34
    print('POLYNOMIAL: 12F1=12F2=44S58+22S50+24S44+12S40+12S36+6S32')
    print('POLYNOMIAL: OP2-OP1=(S40-S36)/4, coefficient equality checked in Q[x,y]')
    print('PALETTE (2,1): OP1 =', val1, '; OP2 =', val2, '; difference =', val2-val1)
    # A fixed leaf adds exactly one singleton cycle for every g. The base
    # profile identity is therefore the symbolic identity after m -> m+r.
    # Check the action construction and expanded polynomials at additional r.
    for r in (1, 2, 10, 100):
        d1, d2 = star_action(SPEC1, r), star_action(SPEC2, r)
        assert all(d1[g] == c1[g]+r and d2[g] == c2[g]+r for g in G)
        assert profile(d1, weights) == profile(d2, weights) == {m+r: c for m, c in fp1.items()}
        dop = add(poly_from_profile(profile(d2, uniform)), scale(-1, poly_from_profile(profile(d1, uniform))))
        assert dop == scale(Q(1, 4), add(star_poly(40+r), scale(-1, star_poly(36+r))))
        assert evaluate(dop, 2, 1) == 15*2**(34+r)
    print('FAMILY: fixed leaves shift every exponent by r; symbolic cycle-profile equality established')
    print('FAMILY: direct expanded checks r=1,2,10,100; palette difference 15*2^(34+r)')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
