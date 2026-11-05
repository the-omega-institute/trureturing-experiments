#!/usr/bin/env python3
"""Exact source-fidelity, admissibility, Bellman and exhaustive-tree checks.

Only the Python standard library is used. Affine vectors have coordinates
(constant, a, b, epsilon, sigma); every coordinate is a Fraction.
"""
from fractions import Fraction as Q
from functools import cache
from itertools import product
import hashlib
import json
from pathlib import Path

ZERO = (Q(0),) * 5

def affine(c=0, a=0, b=0, e=0, s=0):
    return tuple(map(Q, (c, a, b, e, s)))

def add(*vs):
    return tuple(sum(xs, Q(0)) for xs in zip(*vs))

def sub(v, u):
    return tuple(x-y for x,y in zip(v,u))

def at(v, p):
    return v[0] + sum((x*y for x,y in zip(v[1:],p)), Q(0))

BOX = ((Q(41,20),Q(43,20)), (Q(13,4),Q(67,20)),
       (Q(1,1000),Q(1,100)), (Q(1,1000),Q(1,100)))
CORNERS = tuple(product(*BOX))
BASE = (Q(21,10),Q(33,10),Q(1,100),Q(1,1000))
C = (61,81,95)

def inc(h,i):
    return tuple(x+(j==i) for j,x in enumerate(h))

def weight(h):
    n=sum(h)
    if n==0: return affine(1)
    if n==1:
        return affine(a=1,e=h[1]+2*h[2])
    if n==2:
        return {(2,0,0):affine(a=1), (1,1,0):affine(a=1,e=1),
                (1,0,1):affine(b=1), (0,2,0):affine(b=1),
                (0,1,1):affine(b=1,e=1),
                (0,0,2):affine(b=1,e=2)}[h]
    return affine(100+n-3,e=h[1]+2*h[2])


def histories(n):
    return [(i,j,n-i-j) for i in range(n+1) for j in range(n-i+1)]


def admissibility():
    comparisons=0
    for n in range(4):
        for h in histories(n):
            assert all(at(weight(h),p)>=1 for p in CORNERS)
            ws=[weight(inc(h,i)) for i in range(3)]
            for i in range(3):
                assert all(at(sub(ws[i],weight(h)),p)>=0 for p in CORNERS)
                comparisons+=len(CORNERS)
            for i in range(2):
                assert all(at(sub(ws[i+1],ws[i]),p)>=0 for p in CORNERS)
                comparisons+=len(CORNERS)
    assert weight((0,0,0))==affine(1)
    # These symbolic equalities cover all remaining infinitely many histories.
    # With n>=3, a coordinate increment changes w by 1+rank(i)*epsilon.
    # The difference between next-family j and i is (j-i)*epsilon.
    for i in range(3):
        assert all(at(affine(1,e=i),p)>0 for p in CORNERS)
    for i,j in product(range(3), repeat=2):
        if C[j]>=C[i]:
            assert all(at(affine(e=j-i),p)>=0 for p in CORNERS)
    assert 0<C[0]<C[1]<C[2]<=100<=110
    assert all(p[3]>=0 for p in CORNERS)
    print('ADMISSIBILITY: finite boundary comparisons',comparisons,
          '; all n>=3 covered by increments 1+i*epsilon and (j-i)*epsilon')


@cache
def solve(z,h):
    """Bellman implementation; choose policy by the rational base point.

    Stored affine inequalities are subsequently checked at EVERY box corner.
    A third cut costs at least 100 at its child and is strictly dominated by
    the pure branching continuation (at most 7*max(w(h))<100 here).
    """
    w=weight(h)
    if z<=0:
        return w,0,('L',),()
    lv,lc,lt,_=solve(z-100,h)
    rv,rc,rt,_=solve(z-110,h)
    options=[(add(w,lv,rv),lc+rc,('B',lt,rt),'branch')]
    if sum(h)<2:
        for i in range(3):
            v,c,t,_=solve(z-C[i],inc(h,i))
            options.append((add(w,affine(s=1),v),c+1,('C',i,t),'cut'+str(i+1)))
    best=min(options,key=lambda x:(at(x[0],BASE),x[1]))
    return best[0],best[1],best[2],tuple(options)


def bellman():
    v,n,t,actions=solve(250,(0,0,0))
    assert at(v,BASE)==Q(6307,500)
    assert [at(x[0],BASE) for x in actions]==[
        Q(6307,500),Q(12713,1000),Q(12743,1000),Q(6511,500)]
    assert v==affine(3,a=3,b=1,e=1,s=4)
    expected=('B',('C',0,('C',2,('L',))),('C',0,('C',1,('L',))))
    assert t==expected and n==4
    @cache
    def pure(z):
        return 1 if z<=0 else 1+pure(z-100)+pure(z-110)
    assert pure(250)==15
    print('PURE BRANCHING:',pure(250))
    # Populate all Bellman states before taking a snapshot of the cache keys.
    seen=set()
    def walk(z,h):
        if (z,h) in seen:return
        seen.add((z,h))
        if z<=0:return
        walk(z-100,h);walk(z-110,h)
        if sum(h)<2:
            for i in range(3):walk(z-C[i],inc(h,i))
    walk(250,(0,0,0))
    certificates=[]; margins=[]; ties=0;checks=0
    for z,h in sorted(seen):
        v,n,t,opts=solve(z,h)
        for u,k,tree,action in opts:
            delta=sub(u,v)
            vals=[at(delta,p) for p in CORNERS]
            assert min(vals)>=0,(z,h,action,vals)
            if tree!=t:
                if delta==ZERO:
                    assert k>=n;ties+=1
                else:
                    assert min(vals)>0; margins.extend(vals)
                checks+=len(vals)
            certificates.append({'z':z,'h':h,'action':action,
                                 'difference':list(map(str,delta))})
        if z>0 and sum(h)==2:
            # The pure branch subtree has depth <=2 and at most seven nodes.
            assert z<=128
            assert all(7*at(weight(h),p)<100 for p in CORNERS)
    print('BELLMAN: states',len(seen),'nonleaf',sum(z>0 for z,h in seen),
          'alternative corner comparisons',checks,'exact action ties',ties,
          'minimum positive margin',min(margins))
    print('BASE root actions:',[str(at(x[0],BASE)) for x in actions])
    print('UNIQUE root optimum:',at(v,BASE),'cuts',n,'tree',t)
    return expected,certificates


@cache
def enumerate_trees(z,h,n):
    """All trimmed, proving ordered trees of EXACTLY n nodes, not a minimum.

    The dummy children of unary cuts are not counted. No Bellman values or
    Bellman policy decisions are used here. Each tree carries its affine time.
    """
    if n<1:return ()
    if z<=0:
        return ((weight(h),0,('L',)),) if n==1 else ()
    w=weight(h);out=[]
    if sum(h)<2:
        for i in range(3):
            for v,k,t in enumerate_trees(z-C[i],inc(h,i),n-1):
                out.append((add(w,affine(s=1),v),k+1,('C',i,t)))
    for nl in range(1,n-1):
        nr=n-1-nl
        for lv,lc,lt in enumerate_trees(z-100,h,nl):
            for rv,rc,rt in enumerate_trees(z-110,h,nr):
                out.append((add(w,lv,rv),lc+rc,('B',lt,rt)))
    return tuple(out)


def exhaustive(expected):
    # Proposed tree costs <=257/20<13 on the entire box. Since every actual
    # node costs >=1, EVERY competitor no worse has at most twelve nodes.
    # A third cut has a child of weight >=100. Continuing after a proved leaf
    # strictly raises total time. These exclusions cover all omitted trees.
    v=affine(3,a=3,b=1,e=1,s=4)
    assert max(at(v,p) for p in CORNERS)==Q(257,20)
    counts=[];matches=0;comparisons=0;minmargin=None
    for n in range(1,13):
        ts=enumerate_trees(250,(0,0,0),n);counts.append(len(ts))
        for u,k,t in ts:
            if t==expected:
                matches+=1;assert u==v and k==4
            else:
                vals=[at(sub(u,v),p) for p in CORNERS]
                assert min(vals)>0,(n,t,vals)
                minmargin=min(vals) if minmargin is None else min(minmargin,min(vals))
                comparisons+=len(vals)
    assert matches==1
    print('EXHAUSTION: trees by size 1..12',counts,'total',sum(counts))
    print('EXHAUSTION: unique optimum across all 16 corners; comparisons',
          comparisons,'minimum competitor gap',minmargin)
    return {'counts_by_size':counts,'corner_comparisons':comparisons,
            'minimum_gap':str(minmargin)}


def fidelity():
    # Source Example, label ex:separation-placement, literal rationalized data.
    for sigma in [Q(0),Q(1,10),Q(1,5),Q(1,2),Q(1),Q(2)]:
        @cache
        def f(z,h):
            w=Q(1) if h==0 else Q(7,5) if h==1 else Q(100)
            if z<=0:return w,0
            l,lc=f(z-2,h);r,rc=f(z-4,h)
            opts=[(w+l+r,lc+rc)]
            if h<1:
                v,c=f(z-1,h+1);opts.append((w+sigma+v,c+1))
            return min(opts)
        expected=min((Q(39,5)+2*sigma,2),(8+sigma,1),(Q(9),0))
        assert f(5,0)==expected
        print('FIDELITY ex:separation-placement sigma',sigma,'(time,cuts)',
              tuple(map(str,expected)))


def gain_robustness():
    # Compare all possible path prefixes, up to the first threshold crossing,
    # under independent +-1/10 perturbations of Z,ell,r,c1,c2,c3.
    prefixes=set()
    def walk(counts):
        gain=sum(a*b for a,b in zip(counts,(100,110)+C))
        prefixes.add(counts)
        if gain>=250:return
        for i in range(5):
            if i>=2 and sum(counts[2:])>=2:continue
            nxt=list(counts);nxt[i]+=1;walk(tuple(nxt))
    walk((0,)*5)
    for counts in prefixes:
        residual=250-sum(a*b for a,b in zip(counts,(100,110)+C))
        radius=Q(1+sum(counts),10)
        assert abs(residual)>radius,(counts,residual,radius)
    print('GAIN ROBUSTNESS: prefix multisets',len(prefixes),
          '; all threshold signs persist under independent +/-1/10 perturbations')


def main():
    fidelity();admissibility();expected,cert=bellman()
    ex=exhaustive(expected);gain_robustness()
    data={'parameters':{'base':list(map(str,BASE)),
                        'box':[[str(x) for x in pair] for pair in BOX]},
          'affine_coordinates':['constant','a','b','epsilon','sigma'],
          'policy':expected,'bellman_action_differences':cert,'exhaustion':ex}
    Path(__file__).with_name('han_family_certificate.json').write_text(
        json.dumps(data,indent=2)+'\n')
    print('check.py sha256',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')

if __name__=='__main__':main()
