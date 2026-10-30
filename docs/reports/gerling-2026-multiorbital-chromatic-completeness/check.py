#!/usr/bin/env python3
"""Exact checks for Gerling multiorbital completeness; no numerical tolerances."""
from collections import Counter
from itertools import combinations, product
from math import comb
from functools import lru_cache
import sympy as sp

x, y = sp.symbols("x y")
z = x - y
H1 = ((0, 1), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
H2 = ((0, 1), (0, 4), (1, 4), (2, 3), (2, 4), (3, 4))
P = x**5 - 6*x**3*y + 8*x**2*y + 5*x*y**2 - 8*x*y - 4*y**2 + 4*y
T = x**3 - 3*x*y + 2*y


def equal(a, b):
    assert sp.Poly(sp.expand(a - b), x, y).is_zero, (a, b)


def falling(t, k):
    return sp.prod(t - j for j in range(k))


def join(edges, r):
    return tuple(edges) + tuple((u, v) for u in range(5) for v in range(5, 5+r)) + tuple(combinations(range(5, 5+r), 2))


def J(poly, r):
    return sp.expand(sum(comb(r, k)*falling(y, k)*z**(r-k)*poly.subs({x: x-k, y: y-k}, simultaneous=True) for k in range(r+1)))


def quotient(n, edges, perm):
    assert sorted(perm) == list(range(n))
    lab = [-1]*n
    m = 0
    for u in range(n):
        if lab[u] != -1:
            continue
        v = u
        while lab[v] == -1:
            lab[v] = m
            v = perm[v]
        m += 1
    return m, tuple(sorted({tuple(sorted((lab[u], lab[v]))) for u, v in edges}))


@lru_cache(None)
def partition_polynomial(n, edges):
    """Color-kernel partitions: a block containing an edge must use a Z color.

    With b total blocks and a edge-free blocks, choosing k of those a blocks
    to use Y gives binom(a,k) (y)_k (x-y)_(b-k). Distinct blocks have
    distinct colors. This independently enumerates vertex partitions.
    """
    hist = Counter()
    blocks = []
    def visit(v):
        if v == n:
            bad = sum(any((mask >> u)&1 and (mask >> w)&1 for u, w in edges) for mask in blocks)
            hist[len(blocks), len(blocks)-bad] += 1
            return
        for j in range(len(blocks)):
            blocks[j] |= 1 << v
            visit(v+1)
            blocks[j] ^= 1 << v
        blocks.append(1 << v)
        visit(v+1)
        blocks.pop()
    visit(0)
    return sp.expand(sum(count*sum(comb(a, k)*falling(y, k)*falling(z, b-k) for k in range(a+1)) for (b, a), count in hist.items()))


def inclusion_exclusion(n, edges):
    """Full edge-subset IE; each incident component is forced to one Y color."""
    out = Counter()
    m = len(edges)
    for mask in range(1 << m):
        parent = list(range(n))
        touched = 0
        def root(a):
            while parent[a] != a:
                a = parent[a]
            return a
        for j, (u, v) in enumerate(edges):
            if (mask >> j)&1:
                touched |= (1 << u) | (1 << v)
                parent[root(u)] = root(v)
        isolated = n - touched.bit_count()
        components = len({root(v) for v in range(n) if (touched >> v)&1})
        out[isolated, components] += -1 if mask.bit_count()%2 else 1
    return sp.expand(sum(c*x**a*y**b for (a, b), c in out.items()))


def colour_count(n, edges, xx, yy, fixed=None):
    return sum(1 for c in product(range(xx), repeat=n)
               if all(not(c[u] == c[v] and c[u] < yy) for u, v in edges)
               and (fixed is None or all(c[fixed[v]] == c[v] for v in range(n))))


def graph_checks(n, edges, r):
    es = set(edges)
    perm = list(range(n))
    perm[5], perm[6] = perm[6], perm[5]
    assert {tuple(sorted((perm[u], perm[v]))) for u, v in edges} == es
    assert any(perm[v] != v for v in range(n))
    assert all(perm[perm[v]] == v for v in range(n))
    deg = [sum(v in e for e in edges) for v in range(n)]
    seen = {0}
    while True:
        nxt = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if nxt == seen:
            break
        seen = nxt
    assert len(seen) == n and min(deg) >= 3
    return sorted(deg), sum(d == n-1 for d in deg), perm


def main():
    print("Exact arithmetic: integer enumeration and SymPy rational polynomial algebra")
    # Fidelity: source subsection 3.4, reflected P3 versus an edge plus an isolate.
    path = ((0, 1), (1, 2))
    edge = ((0, 1),)
    path_p = partition_polynomial(3, path)
    edge_p = partition_polynomial(3, edge)
    nq1, eq1 = quotient(3, path, (2, 1, 0))
    nq2, eq2 = quotient(3, edge, (1, 0, 2))
    q1 = partition_polynomial(nq1, eq1)
    q2 = partition_polynomial(nq2, eq2)
    equal(path_p, x**3 - 2*x*y + y)
    equal(edge_p, x**3 - x*y)
    equal(q1, x**2 - y)
    equal(q2, x**2 - x*y)
    equal((path_p+q1)/2, (edge_p+q2)/2)
    equal((3*path_p+q1-3*edge_p-q2)/2, y*(1-x))
    print("FIDELITY_OK: equal orbital averages; F(P3)-F(K2 disjoint K1)=y(1-x)")
    for edges in (H1, H2):
        equal(partition_polynomial(5, edges), P)
        equal(inclusion_exclusion(5, edges), P)
    equal(z*(x**2-y)**2+y*((x-1)**2-(y-1))**2, P)
    equal(z*x*T+y*(x-1)*((x-1)**2-(y-1)+T.subs({x:x-1,y:y-1}, simultaneous=True)), P)
    print("BASE_OK: both five-vertex polynomials and both conditioning identities")
    for r in range(6):
        for edges in (H1, H2):
            equal(partition_polynomial(5+r, join(edges, r)), J(P, r))
        if r >= 2:
            for i, edges in enumerate((H1, H2)):
                joined = join(edges, r)
                deg, universal, perm = graph_checks(5+r, joined, r)
                assert universal == r+i
                nq, eq = quotient(5+r, joined, perm)
                equal(partition_polynomial(nq, eq), z*J(P, r-2))
                equal((3*partition_polynomial(5+r, joined)+partition_polynomial(nq, eq))/2,
                      (3*J(P, r)+z*J(P, r-2))/2)
            print(f"FAMILY_OK: r={r}; connected, faithful C2; universal counts {r},{r+1}")
    for i, edges in enumerate((H1, H2)):
        joined = join(edges, 2)
        assert len(joined) == 17
        equal(inclusion_exclusion(7, joined), J(P, 2))
        print(f"IE_OK: seven-vertex graph {i+1}; all {1 << len(joined)} edge subsets")
    display = 3*x**7+x**6-52*x**5*y+144*x**4*y+207*x**3*y**2-388*x**3*y-699*x**2*y**2+820*x**2*y-158*x*y**3+1300*x*y**2-1184*x*y+262*y**3-1120*y**2+864*y
    equal(3*J(P, 2)+z*P, display)
    print("EXPANDED_2F =", sp.expand(display))
    palettes = 0
    for xx in range(5):
        for yy in range(xx+1):
            for r in (0, 2):
                rows = []
                for edges in (H1, H2):
                    n = 5+r
                    joined = join(edges, r)
                    perm = list(range(n))
                    if r:
                        perm[-2], perm[-1] = perm[-1], perm[-2]
                    c = colour_count(n, joined, xx, yy)
                    assert c == J(P, r).subs({x:xx, y:yy})
                    fixed = colour_count(n, joined, xx, yy, perm)
                    nq, eq = quotient(n, joined, perm)
                    assert fixed == colour_count(nq, eq, xx, yy)
                    if r:
                        assert fixed == (z*P).subs({x:xx, y:yy})
                    rows.append((c, fixed))
                assert rows[0] == rows[1]
            palettes += 1
    print(f"COLORINGS_OK: all {palettes} palettes 0<=y<=x<=4; base and r=2; fixed counts and quotient counts")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
