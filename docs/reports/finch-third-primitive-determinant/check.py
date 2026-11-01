#!/usr/bin/env python3
"""Exact checks for Finch arXiv:2607.07688v4, Section 4.

Python 3 plus SymPy. The finite screen visits every pair of rows, bounds the
remaining row by k*sum(max(cofactor,0)), and enumerates that row whenever the
bound permits a determinant at least the proposed third value. Thus omitted
matrices cannot affect the reported top three values. No floating point is used.
"""
from itertools import product
from math import gcd
import sympy as sp


def cross(p, q):
    return (p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2],
            p[0]*q[1]-p[1]*q[0])


def determinant(a):
    c = cross(a[1], a[2])
    return sum(x*y for x, y in zip(a[0], c))


def minors(a):
    return tuple(a[i][u]*a[j][v]-a[i][v]*a[j][u]
                 for i, j in ((0,1),(0,2),(1,2))
                 for u, v in ((0,1),(0,2),(1,2)))


def primitive(a):
    return gcd(*minors(a)) == 1


def symbolic_checks():
    k = sp.symbols('k', integer=True, positive=True)
    h = sp.symbols('h', integer=True, nonnegative=True)
    a, b = k*k, k*k-k
    B, S, T = 2*k**3-2*k**2+k, 2*k**3-2*k**2, k**3+(k-1)**3
    source = ((0,k,k-1),(k,1,k),(k,k,0))
    assert sp.expand(determinant(source)-B) == 0
    assert sp.expand(determinant(((k,0,k-1),(k-1,k,0),(0,k-1,k)))-T) == 0
    witnesses = [((k,k,0),(0,k,k-1),(k,1,k)),
                 ((k,k,0),(0,k-1,k),(k-1,0,k)),
                 ((k,0,k-1),(k-1,k,0),(0,k-1,k))]
    for A, value in zip(witnesses, (B,S,T)):
        assert sp.expand(determinant(A)-value) == 0
    assert sp.expand(B-S-k) == 0
    assert sp.expand(S-T-(k*(k-3)+1)) == 0
    # Six possible ordered factor pairs, written as deficits from k.
    pairs = [(0,0),(0,1),(1,0),(0,2),(2,0),(1,1)]
    def nonnegative_from_three(poly):
        return all(c >= 0 for c in sp.Poly(sp.expand(poly.subs(k,h+3)),h).all_coeffs())
    catalog = set()
    decoded = 0
    for (qa,xi),(pa,eta) in product(pairs, repeat=2):
        P = sp.expand((k-qa)*(k-xi)+(k-pa)*(k-eta))
        if not nonnegative_from_three(P-(2*k*k-3*k+3)):
            continue
        for p1,q2 in product(range(xi+1),range(eta+1)):
            p, q = (p1,p1+k-xi,k-pa),(q2+k-eta,q2,k-qa)
            u,v,nw = map(sp.expand,cross(p,q))
            # u,v can be interchanged in the loss identity.
            if sp.Poly(u,k).all_coeffs() < sp.Poly(v,k).all_coeffs():
                u,v = v,u
            catalog.add((u,v,-nw))
            decoded += 1
    expected = {
        (a,a,a),(a,b,b),(a,b,a),
        (a,a-2*k,a-2*k),(b,b,b),(a,a-2*k,a),
        (a,(k-1)**2,b),(b+1,b,a),
        (b,b,(k-1)**2),(a,a-2*k,b),(b,b,a-1),(b,b,a)}
    expected = {tuple(map(sp.expand,t)) for t in expected}
    assert catalog == expected
    assert decoded == 28 and len(catalog) == 12
    x,y,z = sp.symbols('x y z', integer=True, nonnegative=True)
    for u,v,w in catalog:
        assert nonnegative_from_three(u-1)
        assert nonnegative_from_three(v-1)
        assert nonnegative_from_three(w-1)
        H = sp.expand(k*(u+v))
        R = sp.expand(H-T-1)
        if H == sp.expand(2*k**3):
            assert (u,v,w) == (a,a,a)
            assert nonnegative_from_three(3*a-R-1)
            for j in range(3):
                assert sp.expand(H-j*a) in map(sp.expand,(2*k**3,2*k**3-k*k,S))
        elif H == sp.expand(2*k**3-k*k):
            # 2*b > R: at most one of x,y,z is nonzero, and it equals 1.
            assert all(nonnegative_from_three(t-b) for t in (u,v,w))
            assert nonnegative_from_three(2*b-R-1)
            for loss in (0,u,v,w):
                assert sp.expand(H-loss) in map(sp.expand,(2*k**3-k*k,S,B))
        else:
            assert H in map(sp.expand,(B,S))
            assert all(nonnegative_from_three(t-R-1) for t in (u,v,w))
    # This is algebraic checking of the proposed catalog and loss inequalities,
    # not a Lean proof that the catalog covers every k.
    print('SYMBOLIC source determinant B; three witness determinants B,S,T: OK')
    print('SYMBOLIC 28 decodings; 12 cofactor types; exact loss inequalities: OK')
    for i,(u,v,w) in enumerate(sorted(catalog,key=str),1):
        print(f'  type {i}: ({u}, {v}, -({w}))')


def finite_screen(k):
    rows = tuple(product(range(k+1),repeat=3))
    B, S, T = 2*k**3-2*k**2+k, 2*k**3-2*k**2, k**3+(k-1)**3
    high_all, high_primitive = set(), set()
    eligible = checked = maxP = 0
    for p in rows:
        for q in rows:
            c = cross(p,q)
            P = sum(max(t,0) for t in c)
            maxP = max(maxP,P)
            if k*P < T:
                continue
            eligible += 1
            for r in rows:
                checked += 1
                d = r[0]*c[0]+r[1]*c[1]+r[2]*c[2]
                if d < T:
                    continue
                A = (r,p,q)
                high_all.add(d)
                if primitive(A):
                    high_primitive.add(d)
                if d > B:
                    assert all(m % k == 0 for m in minors(A))
    # Fidelity: source's cube ceiling and primitive levels through k=7.
    assert k*maxP == 2*k**3
    assert sorted(high_primitive,reverse=True) == [B,S,T]
    assert high_all == {T,S,B,2*k**3-k*k,2*k**3}
    source = ((0,k,k-1),(k,1,k),(k,k,0))
    assert determinant(source) == B and primitive(source)
    for A,value in zip([((k,k,0),(0,k,k-1),(k,1,k)),
                        ((k,k,0),(0,k-1,k),(k-1,0,k)),
                        ((k,0,k-1),(k-1,k,0),(0,k-1,k))],(B,S,T)):
        assert determinant(A) == value and primitive(A)
    print(f'EXHAUSTIVE k={k}: row_pairs={len(rows)**2}; eligible={eligible}; '
          f'rows_checked={checked}; ceiling={k*maxP}; '
          f'primitive_top3={sorted(high_primitive,reverse=True)}; '
          f'all_values_ge_T={sorted(high_all,reverse=True)}')


def extended_witness_checks():
    for k in range(3,1001):
        values = (2*k**3-2*k**2+k,2*k**3-2*k**2,k**3+(k-1)**3)
        matrices = (((k,k,0),(0,k,k-1),(k,1,k)),
                    ((k,k,0),(0,k-1,k),(k-1,0,k)),
                    ((k,0,k-1),(k-1,k,0),(0,k-1,k)))
        for A,v in zip(matrices,values):
            assert all(0 <= t <= k for row in A for t in row)
            assert determinant(A) == v and primitive(A)
        assert values[0] > values[1] > values[2] > 0
    print('WITNESSES k=3..1000: bounds, determinant identities, gcd=1, strict ordering: OK')


if __name__ == '__main__':
    symbolic_checks()
    for k in range(3,10):
        finite_screen(k)
    extended_witness_checks()
    print('ALL_CHECKS_OK')
