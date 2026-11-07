#!/usr/bin/env python3
"""Exact independent tableau and isobaric checks for arXiv:2608.11009v1.
All letters are bit masks; tensor factors follow the source's column reading.
No third-party program source is used. No finite check proves stabilization.
"""
from collections import Counter
from functools import lru_cache
from itertools import product


def cells(shape):
    return tuple((r,c) for c in range(max(shape,default=0)) for r in range(len(shape)-1,-1,-1) if c<shape[r])


def highest(shape):
    return tuple(1<<r for r,c in cells(shape))


def valid(shape,t):
    d=dict(zip(cells(shape),t))
    return all(s and ((r,c+1) not in d or s.bit_length()-1 <= (d[r,c+1]&-d[r,c+1]).bit_length()-1) and ((r+1,c) not in d or s.bit_length()-1 < (d[r+1,c]&-d[r+1,c]).bit_length()-1) for (r,c),s in d.items())


def f(t,i):
    """The square-root tensor rule, with epsilon/phi multiplied by two."""
    eps=delta=0; pos=None
    lo,hi=1<<i,1<<(i+1)
    for k,s in enumerate(t):
        a,b=bool(s&lo),bool(s&hi)
        e=1 if a and b else 2 if b else 0
        p=1 if a and b else 2 if a else 0
        if k==0 or eps<p:pos=k
        eps=max(e,eps-2*(int(a)-int(b)))
        delta+=int(a)-int(b)
    if pos is None or not t[pos]&lo:return None
    s=t[pos];u=s&~lo if s&hi else s|hi
    return t[:pos]+(u,)+t[pos+1:]


def ordinary_f(t,i):
    u=f(t,i)
    return None if u is None else f(u,i)


def closure(shape,word,ordinary=False):
    states={highest(shape)}; op=ordinary_f if ordinary else f
    for i in reversed(word):
        out=set(states)
        for t in states:
            seen={t}
            while (t:=op(t,i)) is not None:
                assert t not in seen,'nonterminating color string'
                seen.add(t);out.add(t)
        states=out
    assert all(valid(shape,t) for t in states)
    if ordinary:assert all(all(s.bit_count()==1 for s in t) for t in states)
    return states


def character(shape,states):
    n=len(shape);out=Counter()
    for t in states:
        wt=tuple(sum(bool(s&(1<<i)) for s in t)-shape[i] for i in range(n))
        excess=sum(s.bit_count() for s in t)-sum(shape)
        out[(excess,wt)]+=1
    return out


def pi(poly,i):
    out=Counter()
    for v,k in poly.items():
        d=v[i]-v[i+1]
        for r in range(0,d+1) if d>=0 else range(d+1,0):
            z=list(v);z[i]-=r;z[i+1]+=r
            out[tuple(z)]+=k*(1 if d>=0 else -1)
    return Counter({v:k for v,k in out.items() if k})


def ordinary_character(shape,word):
    p=Counter({tuple(shape):1})
    for i in reversed(word):p=pi(p,i)
    assert all(k>0 for k in p.values())
    return Counter({tuple(v[i]-shape[i] for i in range(len(shape))):k for v,k in p.items()})


def root_coeffs(n,rootset,cap=2):
    """Coefficients on the box [0,cap] in simple-root coordinates, beta=0."""
    p=Counter({(0,)*(n-1):1})
    for i,j in rootset:
        q=Counter()
        for a,v in p.items():
            for k in range(cap+1):
                b=tuple(a[t]+(k if i<=t<j else 0) for t in range(n-1))
                if max(b,default=0)<=cap:q[b]+=v
        p=q
    return p


def as_simple(wt):
    assert sum(wt)==0
    return tuple(-sum(wt[:i+1]) for i in range(len(wt)-1))


def main():
    # Source's explicit standard-crystal character following Definition 3.3.
    c=Counter((s.bit_count()-1,tuple(int(bool(s&(1<<i))) for i in range(3))) for s in range(1,8))
    assert c==Counter({(0,(1,0,0)):1,(0,(0,1,0)):1,(0,(0,0,1)):1,(1,(1,1,0)):1,(1,(1,0,1)):1,(1,(0,1,1)):1,(2,(1,1,1)):1})
    assert closure((1,0,0),(0,1,0))=={(s,) for s in range(1,8)}
    # Source Figure 3.1: the component of {1} tensor {1} has these 17 states.
    expected={(1,1),(1,3),(1,2),(1,7),(3,2),(1,6),(1,5),(2,2),(3,6),(1,4),(3,4),(2,6),(7,4),(2,4),(6,4),(5,4),(4,4)}
    actual=closure((2,0,0),(0,1,0))
    assert actual==expected
    print('FIDELITY: standard S_3 character and 17 tableaux in Figure 3.1 reproduced')
    word=(1,0,2,1);roots=[(0,1),(1,2),(2,3),(0,2),(1,3),(0,3)]
    dims=[30,184,638,1645,3536,6720];weights=[26,125,350,753,1386,2301]
    for m in range(1,7):
        shape=(4*m,3*m,2*m,m)
        ss=closure(shape,word,ordinary=True)
        c=character(shape,ss); assert all(e==0 for e,w in c)
        p=Counter({w:k for (e,w),k in c.items()});assert p==ordinary_character(shape,word)
        assert len(ss)==dims[m-1] and len(p)==weights[m-1]
        for A,B,C in product(range(3),repeat=3):
            if m>=A+C+2*B+1:assert p[(-A,A-B,B-C,C)]==min(B+1,A+C+1)
        print(f'ORDINARY m={m}: tableaux={len(ss)}, distinct_weights={len(p)}, signed-operator agreement')
        if m<=2:
            sv=closure(shape,word);svchar=character(shape,sv)
            assert Counter({w:k for (e,w),k in svchar.items() if e==0})==p
            print(f'SET_VALUED m={m}: tableaux={len(sv)}, max_excess={max(e for e,w in svchar)}, beta^0 agrees with singleton closure')
    p=ordinary_character((32,24,16,8),word)
    for A,B,C in product(range(3),repeat=3):assert p[(-A,A-B,B-C,C)]==min(B+1,A+C+1)
    target={a:min(a[1]+1,a[0]+a[2]+1) for a in product(range(3),repeat=3)}
    matches=[]
    for mask in range(64):
        q=root_coeffs(4,[r for k,r in enumerate(roots) if mask>>k&1])
        if all(q[a]==v for a,v in target.items()):matches.append(mask)
    assert matches==[]
    print('ROOT SUBSETS: all 64 rejected; actual [a,b,c,ab,bc,abc]=[1,1,1,2,2,2]')
    # Uniform-proof test: vary all four gaps, rather than only the cofinal ray.
    # The word acts on the three rightmost row boundaries; all these gaps are positive.
    for gaps in product(range(1,4), repeat=3):
        x,y,z=gaps;shape=(x+y+z+1,y+z+1,z+1,1)
        square=character(shape,closure(shape,word))
        for A,B,C in product(range(2), repeat=3):
            assert square[(0,(-A,A-B,B-C,C))]==min(B+1,A+C+1)
    print('GAPS: all 27 strictly dominant three-gap patterns in {1,2,3}^3 checked')
    # Mixed square-root strings: enumerate exponents up to two at each stage.
    # This is experimental evidence for the uniform local classification lemma.
    shape=(8,6,4,2);states={(highest(shape),())}
    for i in (1,2,0,1):
        out=set()
        for t,path in states:
            for k in range(3):
                if t is None:break
                out.add((t,path+(k,)));t=f(t,i)
        states=out
    paths=sorted(path for t,path in states if (0,(-1,0,0,1)) in character(shape,{t}))
    assert paths==[(0,2,2,2),(2,2,2,0)]
    print('ABC STRINGS: application-order exponents',paths,'; no odd mixed string survives')
    # Every four-coordinate placement in ranks 4 through 7, source operators.
    for n in range(4,8):
        for j in range(n-3):
            shape=tuple(n-r for r in range(n));w=tuple(j+i for i in word)
            ordinary=character(shape,closure(shape,w,ordinary=True))
            square=character(shape,closure(shape,w))
            assert Counter({v:k for (e,v),k in square.items() if e==0})==Counter({v:k for (e,v),k in ordinary.items()})
            for A,B,C in product(range(2),repeat=3):
                wt=[0]*n;wt[j:j+4]=[-A,A-B,B-C,C]
                assert ordinary[(0,tuple(wt))]==min(B+1,A+C+1)
            assert all(all(v[i]==0 for i in range(n) if not j<=i<j+4) for e,v in ordinary)
        print(f'FAMILY n={n}: all {n-3} consecutive blocks checked, including set-valued closure')
    print('ALL_CHECKS_OK')

if __name__=='__main__':main()
