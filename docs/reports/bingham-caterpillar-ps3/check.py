#!/usr/bin/env python3
"""Exact preregistration checks for arXiv:2610.31785v1, Remark 3.11.

All coefficient arithmetic is over integers. SymPy is used only for a
polynomial identity in Z[q,x,y,z]. No floating point is used as evidence.
"""
from collections import Counter, deque
import hashlib
import platform
import sympy as sp


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, r):
    out = [0] * max(len(p), len(r))
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(r):
        out[i] += c
    return trim(out)


def mul(p, r):
    if not p or not r:
        return ()
    out = [0] * (len(p) + len(r) - 1)
    for i, c in enumerate(p):
        for j, d in enumerate(r):
            out[i + j] += c * d
    return trim(out)


def power(p, n):
    assert n >= 0
    r = (1,)
    for _ in range(n):
        r = mul(r, p)
    return r


def mono(n, c=1):
    return trim([0] * n + [c])


def source_graph(parts):
    """v_i indexed 0..l-1, then the source's alpha_i-1 distinct leaves."""
    assert parts and all(isinstance(a, int) and a >= 1 for a in parts)
    assert parts == (1,) or (parts[0] >= 2 and parts[-1] >= 2)
    l = len(parts)
    edges = [(i, i + 1) for i in range(l - 1)]
    n = l
    for i, a in enumerate(parts):
        for _ in range(a - 1):
            edges.append((i, n))
            n += 1
    assert n == sum(parts) and len(edges) == n - 1
    return n, edges


def direct_colorings(parts, k):
    """Enumerate actual proper colorings on the source vertex/edge list."""
    n, edges = source_graph(parts)
    previous = [[] for _ in range(n)]
    for u, v in edges:
        assert u < v
        previous[v].append(u)
    colors = [-1] * n
    coeff = Counter()

    def visit(v, weight):
        if v == n:
            coeff[weight] += 1
            return
        for c in range(k):
            if all(colors[u] != c for u in previous[v]):
                colors[v] = c
                visit(v + 1, weight + c)
        colors[v] = -1

    visit(0, 0)
    p = tuple(coeff[i] for i in range(max(coeff, default=-1) + 1))
    assert sum(p) == k * (k - 1) ** (n - 1)
    return trim(p)


def transfer(parts, k):
    """Sum leaves independently and propagate the endpoint spine color."""
    dp = None
    for a in parts:
        weights = []
        for c in range(k):
            leaf = tuple(int(d != c) for d in range(k))
            weights.append(mul(mono(c), power(leaf, a - 1)))
        if dp is None:
            dp = weights
        else:
            dp = [mul(weights[c], sum_polys(dp[d] for d in range(k) if d != c))
                  for c in range(k)]
    p = sum_polys(dp)
    assert sum(p) == k * (k - 1) ** (sum(parts) - 1)
    return p


def sum_polys(polys):
    out = ()
    for p in polys:
        out = add(out, p)
    return out


def bipartition(parts):
    n, edges = source_graph(parts)
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    sides = [-1] * n
    sides[0] = 0
    todo = deque([0])
    while todo:
        v = todo.popleft()
        for w in adj[v]:
            if sides[w] == -1:
                sides[w] = 1 - sides[v]
                todo.append(w)
            else:
                assert sides[w] != sides[v]
    assert -1 not in sides
    return sides.count(0), sides.count(1)


def family(j):
    assert j >= 1
    return (5, j, j, 4, 3), (4, j, j, 5, 3)


def main():
    print('Python:', platform.python_version(), 'SymPy:', sp.__version__)
    print('Source: arXiv:2610.31785v1, Remark 3.11, first question only')
    a, b = (5, 4, 3), (4, 5, 3)
    for k in (2, 3):
        for parts in (a, b):
            assert direct_colorings(parts, k) == transfer(parts, k)
    assert transfer(a, 3) == transfer(b, 3)
    assert transfer(a, 2) != transfer(b, 2)
    print('FIDELITY base pair (5,4,3)/(4,5,3): ps3 equal; ps2 unequal')
    print('Base common ps3 coefficients (low to high):', transfer(a, 3))
    print('Base ps2 coefficients:', transfer(a, 2), transfer(b, 2))

    for j in range(1, 4):
        for parts in family(j):
            for k in (2, 3):
                assert direct_colorings(parts, k) == transfer(parts, k)
    print('INDEPENDENT GRAPH ENUMERATION: both members, j=1..3, k=2,3 agree with transfer')

    for j in range(1, 31):
        a, b = family(j)
        sa, sb = bipartition(a), bipartition(b)
        assert sorted(sa) == [j + 5, j + 7]
        assert sorted(sb) == [j + 6, j + 6]
        pa, pb = transfer(a, 2), transfer(b, 2)
        assert pa == add(mono(j + 5), mono(j + 7))
        assert pb == mono(j + 6, 2)
        assert pa != pb
        p3 = transfer(a, 3)
        assert p3 == transfer(b, 3)
        digest = hashlib.sha256(','.join(map(str, p3)).encode()).hexdigest()
        print('j=%d vertices=%d bipartitions=%s/%s ps2=UNEQUAL ps3=EQUAL degree=%d sha256=%s'
              % (j, sum(a), sa, sb, len(p3) - 1, digest))
    print('FIDELITY reported family j=1..7 reproduced; extension j=8..30 checked')

    q, x, y, z = sp.symbols('q x y z')
    t = sp.ones(3) - sp.eye(3)
    one = sp.ones(3, 1)

    def diag(a):
        return sp.diag(*[q**c * sum(q**d for d in range(3) if d != c)**(a-1)
                         for c in range(3)])

    assert (diag(3) - sp.diag(q**2*(1+q)**2, q*(1+q**2)**2, q**2*(1+q)**2)).applyfunc(sp.expand) == sp.zeros(3)
    assert (diag(4) - sp.diag(q**3*(1+q)**3, q*(1+q**2)**3, q**2*(1+q)**3)).applyfunc(sp.expand) == sp.zeros(3)
    assert (diag(5) - sp.diag(q**4*(1+q)**4, q*(1+q**2)**4, q**2*(1+q)**4)).applyfunc(sp.expand) == sp.zeros(3)
    arbitrary = sp.diag(x, y, z)
    left = (one.T * diag(3) * t * diag(4) * t * arbitrary * t * arbitrary * t * diag(5) * one)[0]
    right = (one.T * diag(3) * t * diag(5) * t * arbitrary * t * arbitrary * t * diag(4) * one)[0]
    difference = sp.Poly(sp.expand(left - right), q, x, y, z, domain=sp.ZZ)
    assert difference.is_zero
    print('EXACT SYMBOLIC IDENTITY in Z[q,x,y,z]: arbitrary diagonal difference =', difference.as_expr())
    print('Scope: exact arithmetic evidence and a general algebraic route; no Lean settlement proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
