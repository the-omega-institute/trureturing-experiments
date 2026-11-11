#!/usr/bin/env python3
"""Qi--Cui--Xu Problem 2: exact structural, polynomial and Gram checks."""
from collections import Counter
from itertools import combinations
import platform
import sympy as s


def opp(e):
    (i,j),(k,l)=e
    return {(i,l),(k,j)}


def nondeg(e):
    return e[0][0]!=e[1][0] and e[0][1]!=e[1][1]


def comp(e,f):
    return nondeg(e) and nondeg(f) and set(f)==opp(e)


def vertices(cells):
    return {i for i,j in cells},{j for i,j in cells}


def structure(one,two,m,n):
    counts=Counter(c for e in two for c in e)
    assert all(len(set(e))==2 for e in two)
    assert set(counts).isdisjoint(one) and all(v==1 for v in counts.values())
    H=one|set(counts)
    assert all(1<=i<=m and 1<=j<=n for i,j in H)
    for i,k in combinations(range(1,m+1),2):
        for j,l in combinations(range(1,n+1),2):
            assert not {(i,j),(i,l),(k,j),(k,l)}<=one
    D={(a,b) for a,e in enumerate(two) for b,f in enumerate(two)
       if a!=b and nondeg(e) and nondeg(f) and opp(e)<=H and opp(e)&set(f)}
    assert all(not opp(e)<=one for e in two if nondeg(e))
    pairs=0
    edges=[(False,{c}) for c in sorted(one)]+[(True,set(e)) for e in two]
    for (at,a),(bt,b) in combinations(edges,2):
        ar,ac=vertices(a);br,bc=vertices(b)
        if (at or bt) and ar.isdisjoint(br) and ac.isdisjoint(bc):
            pairs+=1
            O={(i,j) for i in ar for j in bc}|{(i,j) for i in br for j in ac}
            assert O-H
    return H,D,pairs


def fidelity():
    one={(1,3),(2,3),(3,1),(3,2),(4,2),(4,3),(5,1),(5,3),(6,1)}
    two=[((1,1),(2,2)),((1,2),(2,1)),((6,2),(3,3))]
    H,D,pairs=structure(one,two,6,3)
    assert len(one)+len(two)==12 and len(H)==15 and pairs==5
    assert D=={(0,1),(1,0)} and comp(two[0],two[1])
    assert {(i,j) for i in range(1,7) for j in range(1,4)}-H=={(4,1),(5,2),(6,3)}
    print('FIDELITY section 5: edges=12; occupied=15; complementary-pair arcs=2; W3_pairs=5')
    x=s.symbols('x1:5');y=s.symbols('y1:5')
    z=s.Matrix([x[0]*y[0],x[0]*y[1],x[0]*y[2],x[1]*y[0],x[1]*y[1]+x[2]*y[2],
                x[1]*y[3],x[2]*y[0]+x[3]*y[1],x[2]*y[1],x[2]*y[3],x[3]*y[2],x[3]*y[3]])
    N=s.zeros(11)
    for i,j,v in [(0,4,1),(1,3,-1),(1,9,1),(2,6,-1)]:N[i,j]=N[j,i]=v
    assert s.expand((z.T*N*z)[0])==0
    t=s.Symbol('t')
    assert s.expand(N.charpoly(t).as_expr()-t**5*(t*t-1)**2*(t*t-2))==0
    Q=s.eye(11)+N/s.sqrt(2)
    assert Q.rank()==10 and all(e.is_nonnegative is True for e in Q.eigenvals())
    assert s.expand((z.T*Q*z)[0]-(z.T*z)[0])==0
    print('FIDELITY section 8: z^T N z=0; charpoly=t^5(t^2-1)^2(t^2-2); exact PSD Gram rank=10<11')


def family(k):
    m=k+1;succ=lambda i:i%k+1;pred=lambda i:(i-2)%k+1
    one={(m,1)};two=[((i,1),(succ(i),2)) for i in range(1,k+1)]
    H,D,pairs=structure(one,two,m,2)
    assert one and len(H)==2*k+1 and pairs==0 and all(nondeg(e) for e in two)
    assert not any(comp(e,f) for e in two for f in two)
    assert D=={(i-1,succ(i)-1) for i in range(1,k+1)}|{(i-1,pred(i)-1) for i in range(1,k+1)}
    assert len(D)==2*k and all((i,(i+1)%k) in D for i in range(k))
    adj={}
    for i,j in H:
        a,b=('r',i),('c',j);adj.setdefault(a,set()).add(b);adj.setdefault(b,set()).add(a)
    seen=set();todo=[next(iter(adj))]
    while todo:
        v=todo.pop()
        if v not in seen:seen.add(v);todo.extend(adj[v]-seen)
    assert seen==set(adj)  # on its own occupied vertices
    x=s.symbols(f'x1:{m+1}');y1,y2=s.symbols('y1 y2')
    f=[x[-1]*y1]+[x[i-1]*y1+x[succ(i)-1]*y2 for i in range(1,k+1)]
    P=s.expand(sum(v*v for v in f))
    M=s.Matrix([[s.expand(v.subs({y1:1,y2:0})).coeff(w) for w in x] for v in f])
    assert M.T*M==s.eye(m) and M.rank()==m
    assert s.expand(P.subs({y1:1,y2:0})-sum(v*v for v in x))==0
    factor='not_requested'
    if k<=6:
        for opts in [{},{'extension':s.sqrt(2)}]:
            c,fs=s.factor_list(P,*x,y1,y2,**opts)
            assert len(fs)==1 and fs[0][1]==1
        # Exact certificate over C(y1,y2), hence R: an x-quadratic of rank >=3
        # cannot split into linear factors; content 1 excludes a pure y factor.
        A=s.hessian(P,x)/2
        assert A.subs({y1:1,y2:0})==s.eye(m) and s.expand(A.det())!=0
        assert A[-1,-1]==y1**2 and A[0,0]==y1**2+y2**2
        assert s.gcd(y1**2,y1**2+y2**2)==1
        factor='QQ:one;Q(sqrt(2)):one;R:absolute-rank-and-content-certificate'
    print(f'FAMILY k={k}: shape={m}x2; length={m}; arcs={len(D)}; complementary=0; '
          f'W3_pairs={pairs}; specialized_rank={m}; connected=True; factor={factor}')


if __name__=='__main__':
    print(f'Python {platform.python_version()}; SymPy {s.__version__}; exact integer/rational/algebraic arithmetic')
    fidelity()
    for k in range(3,25):family(k)
    print('ALL_CHECKS_OK')
