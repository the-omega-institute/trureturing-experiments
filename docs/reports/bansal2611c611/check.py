#!/usr/bin/env python3
"""Exact forest-labeling, Monk-basis and shelling checks; Python standard library only.

All integers are exact. Indices in exponent tuples are zero based; source labels
are one based. Trees are tuples (left,right); an integer is a labeled leaf.
No third-party program source is used. This is an experiment, not a Lean proof.
"""
from functools import lru_cache
from collections import defaultdict
from itertools import product
import json

N = 11
ZERO = (0,) * N

def clean(p):
    return {a: b for a, b in p.items() if b}

def add(p, q, scale=1):
    out = defaultdict(int, p)
    for a, b in q.items(): out[a] += scale * b
    return clean(out)

def mul(p, q):
    out = defaultdict(int)
    for a, b in p.items():
        for c, d in q.items(): out[tuple(x+y for x,y in zip(a,c))] += b*d
    return clean(out)

def monomial(c):
    return {tuple(c)+(0,)*(N-len(c)): 1}

def canon(c):
    c = list(c)
    while c and not c[-1]: c.pop()
    return tuple(c)

def forest(c):
    roots = list(range(1, len(c)+sum(c)+2))
    for i in reversed(range(len(c))):
        for _ in range(c[i]): roots[i:i+2] = [(roots[i], roots[i+1])]
    return tuple(roots)

@lru_cache(None)
def rho(t):
    return t if isinstance(t,int) else rho(t[0])

@lru_cache(None)
def tree_poly(t, lower=1):
    if isinstance(t,int): return {ZERO:1}
    out = {}
    for label in range(lower, rho(t)+1):
        a = [0]*N; a[label-1] = 1
        out = add(out, mul({tuple(a):1}, mul(tree_poly(t[0],label),tree_poly(t[1],label+1))))
    return out

@lru_cache(None)
def poly(c):
    out = {ZERO:1}
    for t in forest(c): out = mul(out, tree_poly(t))
    assert min(out) == tuple(c)+(0,)*(N-len(c)) and out[min(out)] == 1
    return out

@lru_cache(None)
def support(c):
    def leaves(t):
        return [t] if isinstance(t,int) else leaves(t[0])+leaves(t[1])
    s = set()
    for t in forest(c):
        if not isinstance(t,int):
            a = leaves(t); s.update(range(a[0],a[-1]))
    return frozenset(s)

def h(i):
    out={}
    for j in range(i):
        a=[0]*N;a[j]=1;out[tuple(a)]=1
    return out

def expand(p):
    """Independent triangular basis expansion, using the source's lex minimum."""
    p=dict(p);out={}
    while p:
        a=min(p);b=p[a];c=canon(a)
        out[c]=b;p=add(p,poly(c),-b)
    return out

@lru_cache(None)
def monk(c,i):
    out=expand(mul(h(i),poly(c)))
    assert all(v==1 for v in out.values()), (c,i,out)
    return out

def R(p,i):
    out=defaultdict(int)
    for a,b in p.items():
        if a[i-1]==0:out[a[:i-1]+a[i:]+(0,)] += b
    return clean(out)

def T(p,i):
    d=add(R(p,i+1),R(p,i),-1);out={}
    for a,b in d.items():
        assert a[i-1]>0
        x=list(a);x[i-1]-=1;out[tuple(x)]=b
    return out

def coefficient(p,c):
    for i in reversed(range(len(c))):
        for _ in range(c[i]):p=T(p,i+1)
    assert all(a==ZERO for a in p)
    return p.get(ZERO,0)

def weight(u,v,n):
    return sum(monk(u,i).get(v,0) for i in range(1,n+1))

E=(1,0,2,0,0,1)
Z=(2,2,1,1,1)
A=((1,0,2,0,2),(1,0,2,1,0,1),(1,0,3,0,1),(1,0,3,1),(1,2,1,0,0,1),(2,1,1,0,0,1))
B=((1,0,3,1,1),(2,0,2,1,1),(2,1,1,0,2),(2,1,1,1,0,1),(2,2,1,0,1),(2,2,1,1))
EXPECTED_U=[1,1,1,1,2,2]
EXPECTED_M=[[0,2,2,0,0,0],[0,1,0,3,0,0],[2,0,0,0,2,0],[1,0,0,0,0,3],[0,0,0,0,1,1],[0,0,1,1,0,0]]
EXPECTED_V=[4,4,2,1,2,1]

@lru_cache(None)
def window(n):
    if n==0:return ((),)
    out=[]
    def rec(i,tail,total):
        if i<0:out.append(canon(tail));return
        for k in range(n-i-total+1):rec(i-1,(k,)+tail,total+k)
    rec(n-1,(),0)
    return tuple(sorted(set(out),key=lambda c:(sum(c),c)))

@lru_cache(None)
def window_edges(n):
    codes=window(n);code_set=set(codes);down={c:set() for c in codes}
    for c in codes:
        if sum(c)==n:continue
        for i in range(1,n+1):
            for d in monk(c,i):
                if d in code_set:down[d].add(c)
    return down

def interval_down(top,down):
    found={top};todo=[top]
    while todo:
        t=todo.pop()
        for u in down[t]-found:found.add(u);todo.append(u)
    return found

def facets(top,down):
    elements=sorted(interval_down(top,down),key=lambda c:(sum(c),c))
    idx={c:i for i,c in enumerate(elements) if c and c!=top}
    @lru_cache(None)
    def go(v):
        if not v:return (0,)
        bit=0 if v==top else 1<<idx[v]
        return tuple(sorted(bit|f for u in down[v] for f in go(u)))
    fs=go(top);assert len(fs)==len(set(fs))
    assert len({x.bit_count() for x in fs})==1
    return fs

def subfaces(f):
    s=f
    while True:
        yield s
        if s==0:break
        s=(s-1)&f

def shelling(fs):
    """Find a shelling; independently verify all maximal old intersections."""
    remaining=set(fs);seen=set();ordered=[]
    while remaining:
        chosen=None
        for f in sorted(remaining):
            if not ordered:chosen=f;break
            bits=[1<<i for i in range(f.bit_length()) if f>>i&1]
            restriction=sum(b for b in bits if (f^b) in seen)
            if f.bit_count()<=1 or restriction and restriction not in seen:
                chosen=f;break
        if chosen is None:return None
        remaining.remove(chosen);ordered.append(chosen);seen.update(subfaces(chosen))
    assert set(ordered)==set(fs)
    for j,f in enumerate(ordered):
        if j==0:continue
        intersections={f&g for g in ordered[:j]}
        maximal={s for s in intersections if not any(s!=t and s&t==s for t in intersections)}
        assert maximal and all(s.bit_count()+1==f.bit_count() for s in maximal)
    return ordered

def main():
    # Source Figure 2 / first-covers caption: these exact polynomial identities.
    assert poly((0,1,1))==add(add(monomial((1,1)),monomial((1,0,1))),monomial((0,1,1)))
    assert poly((1,0,1))==add(add(monomial((2,)),monomial((1,1))),monomial((1,0,1)))
    G=(2,0,1)
    assert monk((0,1,1),1)=={G:1,(1,1,1):1}
    assert monk((1,0,1),1)=={G:1,(3,):1}
    assert monk((2,),3)=={G:1,(3,):1}
    print('FIDELITY source displayed polynomials and all three Monk identities: OK',flush=True)
    # All 188 source tree tops, with a separately checked shelling for each.
    for n,expected in [(4,14),(5,42),(6,132)]:
        down=window_edges(n);tops=[c for c in window(n) if sum(c)==n]
        assert len(tops)==expected
        largest=0
        for c in tops:
            fs=facets(c,down);s=shelling(fs)
            assert s is not None, ('shelling search failed',n,c,len(fs))
            largest=max(largest,len(fs))
        print(f'FIDELITY rank {n}: {len(tops)}/{expected} tree intervals shellable; largest facet count {largest}',flush=True)
    u=[weight(E,a,7) for a in A];m=[[weight(a,b,7) for b in B] for a in A];v=[weight(b,Z,7) for b in B]
    assert (u,m,v)==(EXPECTED_U,EXPECTED_M,EXPECTED_V)
    L={}
    for i in range(1,8):L=add(L,h(i))
    total=mul(mul(mul(poly(E),L),L),L)
    scalar=coefficient(total,Z);expanded=expand(total).get(Z,0)
    assert scalar==expanded==sum(u[i]*m[i][j]*v[j] for i in range(6) for j in range(6))==50
    assert poly(Z)==monomial(Z)
    for name, cc in [('E',E),('Z',Z)]+[(f'A{i+1}',cc) for i,cc in enumerate(A)]+[(f'B{i+1}',cc) for i,cc in enumerate(B)]:
        print('POLYNOMIAL '+name+' '+json.dumps([[canon(a),b] for a,b in sorted(poly(cc).items())]),flush=True)
    for names, lowers, uppers in [('E-A',(E,),A),('A-B',A,B),('B-Z',B,(Z,))]:
        labels=[[i+1,j+1,[k for k in range(1,8) if vv in monk(uu,k)]] for i,uu in enumerate(lowers) for j,vv in enumerate(uppers) if weight(uu,vv,7)]
        print('COVER_LABELS '+names+' '+json.dumps(labels),flush=True)
    a,b,c,d,e,f,g=[monomial((0,)*i+(1,)) for i in range(7)]
    def plus(*ps):
        q={}
        for p0 in ps:q=add(q,p0)
        return q
    def times(*ps):
        q={ZERO:1}
        for p0 in ps:q=mul(q,p0)
        return q
    def scale(p0,k):return {aa:bb*k for aa,bb in p0.items()}
    H=plus(times(a,a),times(a,b),times(a,c),times(b,b),times(b,c),times(c,c))
    assert poly(E)==times(a,H,plus(a,b,c,d,e,f))
    linA=plus(scale(a,7),scale(b,6),scale(c,5),scale(d,4),f)
    trims=[
        times(a,H,plus(a,b,c,d,e),plus(scale(times(linA,linA),3),scale(times(linA,e),15),scale(times(e,e),19))),
        times(a,H,plus(a,b,c,d),plus(scale(a,63),scale(b,54),scale(c,45),scale(d,29),scale(e,9))),
        times(a,plus(a,b,c),plus(scale(plus(times(a,a),times(a,b),times(b,b)),16),times(plus(a,b,c),plus(scale(a,63),scale(b,54),scale(c,45),scale(d,9))))),
        scale(times(a,plus(a,b),plus(a,b)),25),
        plus(scale(times(a,a),50),scale(times(a,b),25)),
        scale(a,50),{ZERO:50}]
    q=total
    for k, expected in zip([5,4,3,2,2,1,1],trims):
        q=T(q,k);assert q==expected
    print('TRIMMING all seven displayed intermediate polynomial identities: OK',flush=True)
    print('CERTIFICATE u='+json.dumps(u)+' M='+json.dumps(m)+' v='+json.dumps(v),flush=True)
    print('CERTIFICATE [P_Z](L^3 P_E)=50 by triangular expansion AND trimming; uMv=50',flush=True)
    # Directly enumerate the full lower interval in P_7, independent of exhaustion.
    down=window_edges(7);below=interval_down(Z,down)
    up={c:set() for c in below}
    for c in below:
        for d in down[c]:
            if d in below:up[d].add(c)
    above={E};todo=[E]
    while todo:
        for d in up[todo.pop()]-above:above.add(d);todo.append(d)
    between=(above&below)-{E,Z};assert between==set(A)|set(B)
    adj={c:(down[c]|up[c])&between for c in between}
    components=[];unseen=set(between)
    while unseen:
        component={min(unseen)};todo=list(component)
        while todo:
            for d in adj[todo.pop()]-component:component.add(d);todo.append(d)
        unseen-=component;components.append(component)
    assert sorted(map(len,components))==[6,6] and all(len(adj[c])==2 for c in between)
    print('CERTIFICATE full interval (E,Z): 12 vertices, 12 edges, two connected six-cycles',flush=True)
    print('COMPONENTS '+json.dumps([sorted(component) for component in components]),flush=True)
    # Actual link of a saturated chain from the bottom to E, in the proper part.
    chain=[];c=E
    while c:
        chain.append(c);c=min(down[c])
    face=set(chain)
    # Comparable to every face vertex, excluding face vertices, is exactly (E,Z).
    ancestors={c:interval_down(c,down) for c in below}
    link_vertices={c for c in below-{(),Z}-face if all(c in ancestors[d] or d in ancestors[c] for d in face)}
    assert link_vertices==between
    print('LINK chain='+json.dumps(list(reversed(chain)))+'; link is exactly the disconnected graph above',flush=True)
    for m0 in (0,1,2):
        z=(m0+2,2,1,1,1);zn=(m0+3,2,1,1,1)
        assert all(z) and len(support(z))==sum(z) and min(support(z))==1
        assert poly(z)==monomial(z) and mul(h(1),poly(z))==poly(zn)
        assert monk(z,1)=={zn:1}
        print(f'FAMILY m={m0}: rank={sum(z)}, zigzag tree, monomial polynomial, h1 cover to m={m0+1}: OK',flush=True)
    print('ALL_CHECKS_OK')

if __name__=='__main__':main()
