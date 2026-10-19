"""Checks for Svinin, arXiv:1505.06394v2, Conjecture 1, Eq. (48), case (h,n)=(1,1): symbolic check with independent variables u_k,
for r = 1..5 and s = r-1..r+3 (the source proves r = 1, 2, 3)."""
import sympy as sp, itertools
from functools import lru_cache
U=sp.symbols('u0:60'); off=20
def u(k): return U[k+off]
@lru_cache(None)
def T(r,s,i):
    if r==0: return sp.Integer(1)
    if s<=0: return sp.Integer(0)
    tot=sp.Integer(0)
    for lam in itertools.combinations(range(s),r):
        p=sp.Integer(1)
        for a,l in enumerate(lam,1): p*=u(i+l+a-1)
        tot+=p
    return sp.expand(tot)
@lru_cache(None)
def Sp(r,s,i):
    if r==0: return sp.Integer(1)
    if s<=0: return sp.Integer(0)
    tot=sp.Integer(0)
    for lam in itertools.combinations_with_replacement(range(s),r):
        p=sp.Integer(1)
        for a,l in enumerate(lam,1): p*=u(i+l+r-a)
        tot+=p
    return sp.expand(tot)
def F(r,s,i):  # first argument of <,>: zero for s <= 0 when r >= 1
    return Sp(r,s,i)
ok=True; cnt=0
for r in range(1,6):
    for s in range(r-1, r+4):
        i=0
        lhs=sp.expand(sum((-1)**j*F(r-j,s+j,i)*T(j,s-r+j,i+r-j) for j in range(r+1)))
        rhs=sp.expand(sum((-1)**j*u(i+r-j-2)*u(i+s+j)*Sp(r-j-2,s+j+2,i)*T(j,s-r+j,i+r-j) for j in range(r-1)))
        good=sp.expand(lhs-rhs)==0; cnt+=1; ok&=good
        print('r=%d s=%d eq(48) holds: %s'%(r,s,good),flush=True)
print('cases',cnt,'all hold:',ok)
print("RESULT", "PASS" if ok else "FAIL")
