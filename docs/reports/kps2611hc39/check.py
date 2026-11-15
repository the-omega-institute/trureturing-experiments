#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact independent checks for arXiv:2611.16894v1, Question 7(1).

Coordinates (a,b) mean a+b*exp(2*pi*i/3). The metric is the positive
quadratic form Q(a,b)=a*a-a*b+b*b. No floating-point arithmetic is used.
The finite support enumeration uses the minimum-enclosing-circle theorem:
its support is one point, a diameter pair, or three noncollinear points.
Translation by the honeycomb period lattice reduces a support anchor to +1
or -1. Every point of a radius <= sqrt(39)/2 disk is at squared distance
<=39 from that anchor, so Q>=3*a*a/4 and Q>=3*b*b/4 give |a|,|b|<=7.
This checks the counting and arithmetic inputs, not the infinite rigidity
lemma 3.14 or the final claim in Lean.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt

R2 = F(39, 4)

def Q(v):
    a, b = v
    return a*a-a*b+b*b

def add(u, v):
    return (u[0]+v[0], u[1]+v[1])

def sub(u, v):
    return (u[0]-v[0], u[1]-v[1])

def mul(u, v):
    a,b=u; c,d=v
    return (a*c-b*d, a*d+b*c-b*d)

def honey(v):
    return (v[0]+v[1]) % 3 != 0

def circumcenter(v, w):
    va,vb=v; wa,wb=w
    A,B=2*va-vb,2*vb-va
    C,D=2*wa-wb,2*wb-wa
    det=A*D-B*C
    if not det:
        return None
    assert det == 3*(va*wb-vb*wa)
    c=(F(Q(v)*D-Q(w)*B,det), F(A*Q(w)-C*Q(v),det))
    assert Q(c)==Q(sub(c,v))==Q(sub(c,w))
    return c

def integer_points_near_anchor(s):
    return [(a,b) for a,b in product(range(-7,8), repeat=2)
            if Q((a,b)) <= 39 and honey((s+a,b))]

def support_census():
    maxima={False:0, True:0}
    witnesses={}
    total_circles=0
    extremal_centers=set()
    triangle_at_endpoint=0
    for s in (1,-1):
        pool=integer_points_near_anchor(s)
        assert (0,0) in pool
        vs=[v for v in pool if v!=(0,0)]
        circles={((F(0),F(0)),F(0))}
        for v in vs:
            c=(F(v[0],2),F(v[1],2))
            circles.add((c,Q(c)))
        for v,w in combinations(vs,2):
            if Q(sub(v,w))>39:
                continue
            c=circumcenter(v,w)
            if c is None or Q(c)>R2:
                continue
            rr=Q(c)
            if rr==R2:
                triangle_at_endpoint+=1
                # Every endpoint support triangle must have a diameter.
                assert 39 in (Q(v),Q(w),Q(sub(v,w)))
            circles.add((c,rr))
        total_circles+=len(circles)
        for c,rr in sorted(circles):
            inside=[v for v in pool if Q(sub(v,c))<=rr]
            endpoint=rr==R2
            n=len(inside)
            if n>maxima[endpoint]:
                maxima[endpoint]=n
                witnesses[endpoint]=(add(c,(s,0)),rr)
            if n==26:
                assert endpoint
                diameter=[(u,v) for u,v in combinations(inside,2)
                          if Q(sub(u,v))==39]
                assert diameter
                for u,v in diameter:
                    assert add(u,v)==(2*c[0],2*c[1])
                extremal_centers.add(add(c,(s,0)))
        print(f"anchor {s}: relative honeycomb points={len(pool)}, support circles={len(circles)}")
    assert maxima=={False:25,True:26}, maxima
    print(f"fidelity: max count at radius^2=39/4 is {maxima[True]}")
    print(f"fidelity: max count over all strictly smaller radii is {maxima[False]}")
    for endpoint in (False,True):
        print(f"attainment endpoint={endpoint}: center={witnesses[endpoint][0]}, radius^2={witnesses[endpoint][1]}")
    print(f"support circles total={total_circles}; endpoint triangle supports={triangle_at_endpoint}")
    print(f"26-point disk centers witnessed={len(extremal_centers)}; all have antipodal honeycomb points")


def norm_and_triangle_checks():
    # These residue checks are exhaustive on the necessary modulus.
    for a,b in product(range(2),repeat=2):
        assert (Q((a,b))%2==0)==(a==b==0)
    for a,b in product(range(3),repeat=2):
        assert Q((a,b))%3==(a+b)**2%3
    vectors=[(a,b) for a,b in product(range(-7,8),repeat=2)
             if 0<Q((a,b))<=39]
    assert not any(Q(v)==26 for v in vectors)
    endpoints=0
    for v,w in combinations(vectors,2):
        k=v[0]*w[1]-v[1]*w[0]
        q1,q2,q3=Q(v),Q(w),Q(sub(v,w))
        if not k or q3>39:
            continue
        assert 4*q1*q2-(q1+q2-q3)**2==3*k*k
        if 4*q1*q2*q3==117*k*k:
            endpoints+=1
            assert any(q%13==0 for q in (q1,q2,q3))
            assert 26 not in (q1,q2,q3)
            assert 13 not in (q1,q2,q3)
            assert 39 in (q1,q2,q3)
    assert endpoints>0
    roots156=[(a,b) for a,b in product(range(-14,15),repeat=2) if Q((a,b))==156]
    assert roots156 and all(a%2==b%2==0 for a,b in roots156)
    print(f"arithmetic: endpoint triangle pairs={endpoints}; Q=26 impossible; every endpoint triangle has side^2=39")
    print(f"integrality: Q(2z)=156 has {len(roots156)} solutions, all coordinates even")


def quotient_checks():
    units=[(a,b) for a,b in product(range(-1,2),repeat=2) if Q((a,b))==1]
    assert len(units)==6
    nonzero_parities={(1,0),(0,1),(1,1)}
    assert {(a%2,b%2) for a,b in units}==nonzero_parities
    zs=[(a,b) for a,b in product(range(-7,8),repeat=2) if Q((a,b))==39]
    assert len(zs)==12
    # Classes of 2c modulo 2*Lambda0; CRT invariants are (a+b mod3,a mod2,b mod2).
    t_reps={}
    for a,b in product(range(6),repeat=2):
        t_reps.setdefault(((a+b)%3,a%2,b%2),(a,b))
    assert len(t_reps)==12
    ts=[t for key,t in sorted(t_reps.items()) if key[0]!=0 and key[1:]!=(0,0)]
    assert len(ts)==6
    hist=Counter()
    for z in zs:
        a,b=z
        assert (a+b)%3==0
        assert {(mul(e,z)[0]%2,mul(e,z)[1]%2) for e in units}==nonzero_parities
        def cls(v):
            x,y=v
            return (((a-b)*x+b*y)%39,(-b*x+a*y)%39)
        # 39*Lambda is contained in z*Lambda, so this is a complete quotient.
        classes={cls(v) for v in product(range(39),repeat=2)}
        hclasses={cls(v) for v in product(range(39),repeat=2) if honey(v)}
        assert len(classes)==39 and len(hclasses)==26
        assert cls(z)==cls(mul(z,(0,1)))==cls((0,0))
        for t in ts:
            c=(F(t[0],2),F(t[1],2))
            pts=[(x,y) for x in range(-8,12) for y in range(-8,12)
                 if honey((x,y)) and Q(sub((x,y),c))<=R2]
            assert len(pts)<=26
            choices=[mul(e,z) for e in units if all((mul(e,z)[i]-t[i])%2==0 for i in (0,1))]
            assert choices
            for v in choices:
                u=tuple((t[i]+v[i])//2 for i in (0,1))
                w=tuple((t[i]-v[i])//2 for i in (0,1))
                assert u!=w and honey(u) and honey(w)
                assert Q(sub(u,c))==Q(sub(w,c))==R2
                assert u in pts and w in pts and sub(u,w)==v
                assert cls(u)==cls(w)
            covered={cls(p) for p in pts}
            assert len(covered)<26
            assert hclasses-covered
            hist[(len(pts),len(covered))]+=1
    print(f"quotients: {len(zs)} norm-39 generators x {len(ts)} center classes; index=39, honeycomb classes=26")
    print("(disk points, represented classes) histogram:", sorted(hist.items()))
    print("residue collision: every tested disk omits a honeycomb class; translates by z*Lambda preserve that omission")


def conditional_prime_family_checks():
    # Only the arithmetic step is generalized: no counting hypothesis is asserted.
    primes=[p for p in range(7,44) if p%3==1 and all(p%d for d in range(2,isqrt(p)+1))]
    count=0
    for p in primes:
        m=3*p
        bound=isqrt(4*p)+1
        vs=[(a,b) for a,b in product(range(-bound,bound+1),repeat=2) if 0<Q((a,b))<=m]
        for v,w in combinations(vs,2):
            k=v[0]*w[1]-v[1]*w[0]
            qq=(Q(v),Q(w),Q(sub(v,w)))
            if not k or qq[2]>m:
                continue
            if 4*qq[0]*qq[1]*qq[2]==9*p*k*k:
                count+=1
                assert m in qq
    print(f"conditional family arithmetic: p={primes}, {count} endpoint triangle pairs; no counting hypotheses assumed")


def main():
    assert F(8, 1) / (3 * F(4, 39)) == 26
    print("density equality: d0^2=4/39, 8/(3*d0^2)=26")
    support_census()
    norm_and_triangle_checks()
    quotient_checks()
    conditional_prime_family_checks()
    print("NOT_CHECKED: infinite covering compactness/Voronoi equality reduction (Lemma 3.14); Lean settlement")
    print("ALL_CHECKS_OK")

if __name__=='__main__':
    main()
