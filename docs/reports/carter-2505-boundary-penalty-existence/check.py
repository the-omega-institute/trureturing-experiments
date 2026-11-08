#!/usr/bin/env python3
"""Exact checks for Carter 2505.20005v2 section 7.2; no floating point.
Finite checks support an analytic proposal, not a proof of its universal claim.
"""
from fractions import Fraction as F
from itertools import combinations, product, permutations
from functools import lru_cache
from math import comb
from random import Random
rng = Random(250520005)
def mat(a): return tuple(tuple(F(x) for x in row) for row in a)
def trans(a): return tuple(zip(*a))
def mul(a,b): return tuple(tuple(sum((x*y for x,y in zip(row,col)),F(0)) for col in trans(b)) for row in a)
def add(a,b): return tuple(tuple(x+y for x,y in zip(ar,br)) for ar,br in zip(a,b))
def scale(c,a): return tuple(tuple(c*x for x in row) for row in a)
def eye(n): return mat([[int(i==j) for j in range(n)] for i in range(n)])
def block(a,inds): return tuple(tuple(a[i][j] for j in inds) for i in inds)
def det(a):
    if not a: return F(1)
    return sum(((-1)**j*a[0][j]*det(tuple(tuple(row[k] for k in range(len(a)) if k!=j) for row in a[1:])) for j in range(len(a))),F(0))
def det2(a):
    out=F(0)
    for p in permutations(range(len(a))):
        t=F((-1)**sum(p[i]>p[j] for i in range(len(a)) for j in range(i+1,len(a))))
        for i in range(len(a)): t*=a[i][p[i]]
        out+=t
    return out
def prod(xs):
    t=F(1)
    for x in xs: t*=x
    return t
def trace(a,b): return sum((a[i][j]*b[j][i] for i in range(len(a)) for j in range(len(a))),F(0))
def norm2(a): return sum((x*x for row in a for x in row),F(0))
def spd(n):
    l=mat([[rng.randint(1,3) if i==j else rng.randint(-2,2) if j<i else 0 for j in range(n)] for i in range(n)])
    # Scale rows, preserving triangular invertibility.
    l=tuple(tuple(c*x for x in row) for row in l for c in [rng.choice([F(1,4),F(1),F(3)])])
    h=mul(l,trans(l))
    assert all(det(block(h,range(k)))>0 for k in range(1,n+1))
    return h
@lru_cache(None)
def log_interval(q):
    assert q>0
    k=0; r=q
    while r>=2: r/=2; k+=1
    while r<1: r*=2; k-=1
    def series(x):
        u=(x-1)/(x+1)
        lo=2*sum((u**(2*j+1)/F(2*j+1) for j in range(40)),F(0))
        return lo,lo+2*u**81/(81*(1-u*u))
    lo,hi=series(r); l2,h2=series(F(2))
    return (lo+k*l2,hi+k*h2) if k>=0 else (lo+k*h2,hi+k*l2)
def objective(s,h,a,b):
    lo,hi=log_interval(det(h)); lo-=trace(s,h); hi-=trace(s,h)
    for i in range(len(h)):
        c=a if h[i][i]<=1 else b
        dl,dh=log_interval(h[i][i]); lo-=c*dh; hi-=c*dl
    return lo,hi

# Source fidelity: Proposition 4, det(I+t vv^T)=1+t for a unit null vector.
v=mat([[F(3,5)],[F(4,5)]]); w=mat([[F(4,5)],[F(-3,5)]])
s=mul(w,trans(w))
for t in [F(k,3) for k in range(1,65)]:
    h=add(eye(2),scale(t,mul(v,trans(v))))
    assert det(h)==1+t and trace(s,h)==trace(s,eye(2))==1
    assert all(h[i][i]<=1+t for i in range(2))
print('FIDELITY Proposition 4: 64 rational rank-one paths')
count=0
for n in range(1,5):
    for a in product([-1,0,1],repeat=n):
        for b in product([-1,0,1],repeat=n):
            assert 2*sum((x+y)**2 for x,y in zip(a,b))-sum(x*x for x in a)+2*sum(y*y for y in b)==sum((x+2*y)**2 for x,y in zip(a,b))
            count+=1
print(f'TRACE_IDENTITY all {count} vector pairs, lengths 1..4, entries -1,0,1')
count=0
for p in range(1,5):
    for _ in range(40):
        h=spd(p); u=mat([[rng.randint(-2,2) for _ in range(p)] for _ in range(p)]); s=mul(u,trans(u))
        assert det(h)==det2(h)
        b=tuple(i for i in range(p) if h[i][i]>=1); c=tuple(i for i in range(p) if i not in b)
        assert det(h)<=det(block(h,b))*det(block(h,c))
        assert det(block(h,c))<=prod(h[i][i] for i in c)
        assert trace(s,h)>=trace(block(s,b),block(h,b))/2-trace(block(s,c),block(h,c))
        assert trace(block(s,c),block(h,c))<=sum(abs(x) for row in s for x in row)
        count+=1
print(f'BLOCK_REDUCTION {count} rational pairs, dimensions 1..4')
count=0
for p in range(1,6):
    for r in range(1,p+1):
        a=mat([[F(j+1)**i for j in range(p)] for i in range(r)]); s=mul(trans(a),a); k=p-r
        for _ in range(4):
            h=spd(p); t=trace(s,h); lhs=rhs=F(1); sets=list(combinations(range(p),k))
            for inds in sets:
                l=a+tuple(eye(p)[i] for i in inds); dl=det(l)
                assert dl!=0
                left=dl**2*det(h); right=(t/r)**r*prod(h[i][i] for i in inds)
                assert left<=right
                lhs*=left; rhs*=right; count+=1
            exponent=comb(p-1,k-1) if k else 0
            assert rhs==(t/r)**(r*len(sets))*prod(h[i][i]**exponent for i in range(p))
            assert lhs<=rhs
print(f'GENERAL_POSITION {count} maximal-minor bounds and exact products')
count=0
for p in range(1,6):
    v=mat([[i+1] for i in range(p)]); q=add(eye(p),scale(-2/norm2(v),mul(v,trans(v))))
    assert mul(q,trans(q))==eye(p)
    for r in range(1,p+1):
        d=mat([[i+1 if i==j and i<r else 0 for j in range(p)] for i in range(p)]); s=mul(mul(q,d),trans(q))
        assert all(s[i][i]>0 for i in range(p))
        for _ in range(5):
            h=spd(p); t=max(h[i][i] for i in range(p)); m=p-r; T=trace(s,h)
            u=max([abs(h[i][j]) for i in range(p) for j in range(i+1,p)] or [F(0)])
            amin=min(s[i][i] for i in range(p)); b=2*sum(abs(s[i][j]) for i in range(p) for j in range(i+1,p)); c0=max(F(1),b)/amin
            assert T>=amin*t-b*u and t<=c0*(T+u)
            assert T>=t/(2*c0) or u>=t/(2*c0)
            assert det(h)<=(T/r)**r*((p*t/m)**m if m else 1)
            count+=1
print(f'JOINT_PENALTY {count} range/kernel and trace/off-diagonal bounds')

# Rational data in Gaussian support, not Gaussian-distributed draws.
# Source special case P=c log x: the analytic maximizer has x=y=d, z as below
# after rescaling by the nonzero data coordinates. Test stationarity and Hessian.
maxima=comparisons=0
for c in [F(3,5),F(2,3),F(3,4),F(4,5)]:
    for _ in range(5):
        a,b=(rng.choice([-3,-2,-1,1,2,3]) for _ in range(2)); s=mat([[a*a,a*b],[a*b,b*b]])
        d=c*(1-c)/(2*c-1); z=-(1-c)**2/(2*c-1); h=mat([[d/(a*a),z/(a*b)],[z/(a*b),d/(b*b)]])
        x,y,z0,delta=h[0][0],h[1][1],h[0][1],det(h)
        assert delta>0 and x>0 and y>0
        assert y/delta-a*a-c/x==x/delta-b*b-c/y==-2*z0/delta-2*a*b==0
        hes=mat([[-y*y/delta**2+c/x**2,1/delta-x*y/delta**2,2*y*z0/delta**2],
                 [1/delta-x*y/delta**2,-x*x/delta**2+c/y**2,2*x*z0/delta**2],
                 [2*y*z0/delta**2,2*x*z0/delta**2,-2/delta-4*z0*z0/delta**2]])
        assert all((-1)**k*det(block(hes,range(k)))>0 for k in range(1,4))
        flo,_=objective(s,h,c,c)
        for candidate in [spd(2) for _ in range(6)]:
            _,upper=objective(s,candidate,c,c); assert upper<flo; comparisons+=1
        maxima+=1
print(f'STRICT_THRESHOLD {maxima} analytic maximizers, {comparisons} certified log comparisons')
previous=F(0)
for t in range(1,65):
    equality=F(t*t)/(F(t)+F(1,4))**2
    below=F((1+2*t)**3,(1+t)**2)
    assert previous<equality<1 and below>=2*t
    previous=equality
print('BOUNDARY_AND_BELOW 64 exact path values: equality limit; subthreshold divergence')
s=mat([[1,1,0],[1,1,0],[0,0,1]]); a,b=F(3,4),F(2,3); count=0
for rate in [1,2,5]:
    previous=None
    for k in range(1,21):
        t=F(2**k); h=mat([[t+1,-t,0],[-t,t+1,0],[0,0,t**(-rate)]])
        lo,hi=objective(s,h,a,b)
        if previous is not None: assert hi<previous
        previous=lo; coeff=1-2*b-rate*(1-a); assert coeff<0
        lt,_=log_interval(t); _,l3=log_interval(F(3)); assert hi<l3+coeff*lt-2
        count+=1
print(f'MIXED_BOUNDARY {count} certified decreases, shrinking rates 1,2,5')
print('ALL_CHECKS_OK')
