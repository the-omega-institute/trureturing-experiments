#!/usr/bin/env python3
"""Exact checks for Hohn's OEIS A386671 bilateral-symmetry conjecture.

Original checker; Python standard library only. Numerical examples are attributed
in README.md. No floating-point operation is used.
"""
from fractions import Fraction
from functools import cmp_to_key
from math import gcd

K = 48612265
BASE = [(40723,3636),(39917,8844),(-39204,11603),(-14893,38076),
        (-16236,-37523),(-19667,35844),(-22116,-34387),(31476,-26093)]
MULTIPLIERS = [(-660,989),(1020,-611),(1020,611)]
EDGES = [(-30473184,37875287),(-35091936,33640873),(14399273,-46430736),
         (-27827784,-39859337),(47826007,8707776),(-22469496,-43107703),
         (48605303,822696),(5031817,48351144),(43759056,-21173033),
         (46119024,-15368407),(-32898647,35788704),(8073576,47937143),
         (-39487273,-28353264),(1840344,48577417),(-43568777,-21561864),
         (16162697,-45846696),(39315864,28590473),(35311656,33410167),
         (-47077513,-12118584),(-38455296,29737897),(6365833,-48193656),
         (-41961024,24544343),(-1547863,-48587616),(48048343,-7383024)]
ORDER = [7,5,17,18,12,8,14,1,11,2,20,22,19,15,13,4,6,23,21,3,16,9,10,24]
SOURCE_EXAMPLE = [(38076,14893),(35844,19667),(11603,39204),(-34387,22116),
                  (-37523,16236),(-26093,-31476),(3636,-40723),(8844,-39917)]


def add(z,w):
    return z[0]+w[0], z[1]+w[1]


def sub(z,w):
    return z[0]-w[0], z[1]-w[1]


def mul(z,w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conj(z):
    return z[0], -z[1]


def norm_sq(z):
    return z[0]*z[0]+z[1]*z[1]


def cross(z,w):
    return z[0]*w[1]-z[1]*w[0]


def angle_cmp(z,w):
    """Increasing angle in [0,2*pi), using integers and half-plane tests."""
    h = lambda v: 0 if v[1]>0 or (v[1]==0 and v[0]>=0) else 1
    if h(z)!=h(w):
        return -1 if h(z)<h(w) else 1
    c=cross(z,w)
    return -1 if c>0 else 1 if c<0 else 0


def vertices(edges):
    result=[(0,0)]
    for e in edges:
        result.append(add(result[-1],e))
    assert result[-1]==(0,0)
    return result[:-1]


def permitted_closed(edges,k):
    assert k>0 and len(edges)>0
    assert all(x and y and gcd(x,y)==1 and x*x+y*y==k*k for x,y in edges)
    assert len({tuple(sorted((abs(x),abs(y)))) for x,y in edges})==len(edges)
    assert (sum(x for x,y in edges),sum(y for x,y in edges))==(0,0)


def support_values(edges):
    vs=vertices(edges)
    return [cross(e,sub(vs[j],vs[i])) for i,e in enumerate(edges)
            for j in range(len(edges)) if j not in (i,(i+1)%len(edges))]


def reflection_failures(edges):
    """All possible linear reflections: image of edge 0 determines lambda.

    R(z) = lambda*conj(z). Equal length gives
    lambda = edge_j*edge_0 / k^2. Check numerators without division.
    Each returned pair (j,i) certifies that candidate j fails at edge i.
    """
    n2=norm_sq(edges[0])
    targets={(n2*x,n2*y) for x,y in edges}
    failures=[]
    for j,w in enumerate(edges):
        c=mul(edges[0],w)
        assert norm_sq(c)==n2*n2
        bad=[i for i,z in enumerate(edges) if mul(c,conj(z)) not in targets]
        if bad:
            failures.append((j+1,bad[0]+1))
    return failures


def vertex_reflections(edges):
    """Independent affine-axis test about the exact rational vertex centroid.

    A boundary reflection of a strictly convex polygon permutes its vertices,
    fixes their centroid, and maps the first centered vertex to some vertex.
    """
    vs=vertices(edges)
    c=(Fraction(sum(x for x,y in vs),len(vs)),Fraction(sum(y for x,y in vs),len(vs)))
    centered=[sub(v,c) for v in vs]
    a=centered[0]
    assert norm_sq(a)>0
    n2=norm_sq(a)
    target=set(centered)
    candidates=0
    symmetries=0
    for b in centered:
        if norm_sq(b)!=n2:
            continue
        candidates+=1
        lam=tuple(x/n2 for x in mul(a,b))
        if {mul(lam,conj(z)) for z in centered}==target:
            symmetries+=1
    return candidates,symmetries


def mod29_checks(edges):
    assert K%29==0 and 12*12%29==28
    plus=[z for z in edges if (z[1]-12*z[0])%29==0]
    minus=[z for z in edges if (z[1]+12*z[0])%29==0]
    assert len(plus)==8 and len(minus)==16
    assert all((x%29,y%29)!=(0,0) for x,y in edges)
    assert all((y+12*x)%29==0 for x,y in edges[:16])
    assert all((y-12*x)%29==0 for x,y in edges[16:])
    # Exhaust the nonzero isotropic vectors in F_29^2 and their product law.
    iso=[(x,y) for x in range(29) for y in range(29)
         if (x or y) and (x*x+y*y)%29==0]
    assert len(iso)==56
    for z in iso:
        ev=((z[0]+12*z[1])%29,(z[0]-12*z[1])%29)
        assert (ev[0]==0)!=(ev[1]==0)
        for w in iso:
            p=mul(z,w)
            product_zero=p[0]%29==0 and p[1]%29==0
            opposite=((z[1]-12*z[0])%29==0)!=((w[1]-12*w[0])%29==0)
            assert product_zero==opposite
    return len(plus),len(minus),len(iso)**2


def main():
    # Reproduce the displayed A386671 example (not its unproved minimality).
    permitted_closed(SOURCE_EXAMPLE,40885)
    se=sorted(SOURCE_EXAMPLE,key=cmp_to_key(angle_cmp))
    assert min(support_values(se))>0
    assert len(reflection_failures(se))<8
    assert vertex_reflections(se)[1]>0
    assert 8*40885==327080
    print('FIDELITY: source eight-edge example closes; all norms=40885^2; all gcds=1;')
    print('FIDELITY: eight distinct leg pairs; boundary has a reflection; total length=327080.')
    assert EDGES==[mul(b,u) for u in MULTIPLIERS for b in BASE]
    permitted_closed(BASE,40885)
    assert all(norm_sq(u)==1189**2 for u in MULTIPLIERS)
    assert K==40885*1189 and K*K==2363152308430225
    permitted_closed(EDGES,K)
    assert sorted(range(1,25),key=cmp_to_key(lambda i,j:angle_cmp(EDGES[i-1],EDGES[j-1])))==ORDER
    ordered=[EDGES[i-1] for i in ORDER]
    supports=support_values(ordered)
    assert len(supports)==528 and min(supports)==90286825602984
    print(f'WITNESS: k={K}; k^2={K*K}; edges=24; gcd/norm/nonzero/distinct/closure checks=OK.')
    print('ANGLE_ORDER [0,2*pi):',','.join(map(str,ORDER)))
    neg=[i for i in ORDER if EDGES[i-1][1]<0]
    pos=[i for i in ORDER if EDGES[i-1][1]>0]
    print('ANGLE_ORDER (-pi,pi]:',','.join(map(str,neg+pos)))
    print(f'STRICT_CONVEXITY: support_tests={len(supports)}; minimum={min(supports)}.')
    print('VERTICES [0,2*pi):',vertices(ordered))
    plus,minus,products=mod29_checks(EDGES)
    print(f'MOD29: y=12x count={plus}; y=-12x count={minus}; finite product tests={products}.')
    failures=reflection_failures(EDGES)
    assert len(failures)==24
    print('REFLECTION_EXHAUSTION: (candidate image, failing edge)=',failures)
    vc,vs=vertex_reflections(ordered)
    assert vs==0
    print(f'INDEPENDENT_AFFINE_AXES: norm-compatible images of vertex 0={vc}; symmetries={vs}.')
    # A symbolic infinite-family mechanism is documented in README. These are
    # exact sample instances, not computational evidence for an infinite quantifier.
    for m in [0,1,2,3,10,100]:
        t=2*K*(m+3)
        a,b=t*t-1,2*t
        c=t*t+1
        assert t>4*K+1 and gcd(a,b)==1 and gcd(c,K)==1
        assert a*a+b*b==c*c
        es=[mul(z,(a,b)) for z in EDGES]
        permitted_closed(es,K*c)
        assert len(reflection_failures(es))==24
        so=sorted(es,key=cmp_to_key(angle_cmp))
        assert min(support_values(so))==min(supports)*c*c
        assert vertex_reflections(so)[1]==0
        assert mod29_checks(es)[:2]==(8,16)
        print(f'FAMILY_SAMPLE: m={m}; length={K*c}; all exact checks=OK.')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
