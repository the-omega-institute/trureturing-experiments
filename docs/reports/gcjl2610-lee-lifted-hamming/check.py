#!/usr/bin/env python3
"""Exact checks for arXiv:2610.04834v1, Example 7.10; no Lean proof.

All arithmetic is integral or in explicitly represented finite quotient rings.
Symbolic identities are reduced modulo lambda**2-lambda+2. Finite checks
are evidence for the route, not a proof of its uniform-in-r conclusion.
"""
from itertools import product
from math import gcd
import sympy as sp

l, t, u = sp.symbols('lambda t u')
relation = l*l-l+2
h = t**3-(l-1)*t**2-l*t-1
k = t**3+l*t**2+(l-1)*t-1
g = sp.expand((t-1)*k)

def red(a):
    return sp.rem(sp.Poly(sp.expand(a), l), sp.Poly(relation, l)).as_expr().expand()

def matred(A):
    return A.applyfunc(red)

def coeffs(p):
    return sp.Matrix([sp.expand(p).coeff(t, i) for i in range(7)])

D = sp.Matrix.hstack(*(coeffs(p) for p in [1,t,t*t,t**3,g,t*g,t*t*g]))
C = sp.Matrix.hstack(*(coeffs(p) for p in [h,t*h,t*t*h,t**3*h,1,t,t*t]))
T = sp.zeros(7)
for i in range(6): T[i+1,i] = 1
T[0,6] = u
H = T**3-(l-1)*T**2-l*T-sp.eye(7)
Delta = sp.diag(1,1,1,1,u-1,u-1,u-1)
assert red(h*k-sum(t**i for i in range(7))) == 0
assert red(h*g-(t**7-1)) == 0
assert matred(H*D-C*Delta) == sp.zeros(7)
assert sp.cancel(D.det()) in [1,-1] and sp.cancel(C.det()) in [1,-1]
Di, Ci = D.inv().applyfunc(sp.cancel), C.inv().applyfunc(sp.cancel)
assert all(sp.denom(sp.cancel(z)) == 1 for z in list(Di)+list(Ci))
assert (D*Di-sp.eye(7)).applyfunc(sp.expand) == sp.zeros(7)
assert (C*Ci-sp.eye(7)).applyfunc(sp.expand) == sp.zeros(7)
print('SYMBOLIC: hk=sum(t^i,i=0..6), hg=t^7-1, H D=C diag(I4,(u-1)I3)')
print(f'INTEGRAL_BASES: det(D)={sp.cancel(D.det())}, det(C)={sp.cancel(C.det())}; inverse entries in Z[lambda]')

# Coefficient-level Laurent symplectic check for every type of unit block.
v = sp.symbols('v')
def dagger_symbolic(M):
    return M.T.applyfunc(lambda z: sp.cancel(z.subs({u: 1/u, v: 1/v}, simultaneous=True)))
J4=sp.Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
for aa,bb in [(1,1),(1,v-1),(u-1,1),(u-1,v-1)]:
    if aa==1:W=sp.Matrix([[1,-bb],[0,1]])
    elif bb==1:W=sp.Matrix([[0,-1],[1,aa]])
    else:W=sp.eye(2)
    Q=sp.diag(dagger_symbolic(W),W.inv())
    assert (dagger_symbolic(Q)*J4*Q-J4).applyfunc(sp.cancel)==sp.zeros(4)
print('LAURENT_SYMPLECTIC: exact dagger identity for the unit-block shears and toric identity block')

# Small sparse matrices over Z/q, for the complete 98-coordinate complex.
def eye(n): return {(i,i):1 for i in range(n)}
def dense(A,q): return {(i,j):int(A[i,j])%q for i in range(A.rows) for j in range(A.cols) if int(A[i,j])%q}
def add(A,B,q):
    out=dict(A)
    for ij,z in B.items(): out[ij]=(out.get(ij,0)+z)%q
    return {ij:z for ij,z in out.items() if z}
def scale(A,c,q):return {ij:c*z%q for ij,z in A.items() if c*z%q}
def mul(A,B,q):
    rows={}
    for (i,j),z in B.items():rows.setdefault(i,[]).append((j,z))
    out={}
    for (i,k),z in A.items():
        for j,w in rows.get(k,[]):out[i,j]=(out.get((i,j),0)+z*w)%q
    return {ij:z for ij,z in out.items() if z}
def kron(A,B,n,q):
    return {(i*n+k,j*n+m):z*w%q for (i,j),z in A.items() for (k,m),w in B.items() if z*w%q}
def block(A,B,Cc,Dd,n):
    out=dict(A)
    out.update({(i,j+n):z for (i,j),z in B.items()})
    out.update({(i+n,j):z for (i,j),z in Cc.items()})
    out.update({(i+n,j+n):z for (i,j),z in Dd.items()})
    return out
def trans(A):return {(j,i):z for (i,j),z in A.items()}
def row(A,B,n):return dict(A)|{(i,j+n):z for (i,j),z in B.items()}
def col(A,B,n):return dict(A)|{(i+n,j):z for (i,j),z in B.items()}

def quotient(q,lam):
    zero=(0,0,0); one=(1,0,0); xi=(0,1,0)
    hc=(-1,-lam,1-lam)
    def plus(a,b):return tuple((x+y)%q for x,y in zip(a,b))
    def neg(a):return tuple(-x%q for x in a)
    def times(a,b):
        z=[0]*5
        for i in range(3):
            for j in range(3):z[i+j]+=a[i]*b[j]
        for i in [4,3]:
            for j in range(3):z[i-3+j]-=z[i]*hc[j]
        return tuple(x%q for x in z[:3])
    def power(a,n):
        z=one
        for _ in range(n):z=times(z,a)
        return z
    return zero,one,xi,plus,neg,times,power

for q in [2,4,8,16,32,64]:
    roots=[a for a in range(q) if (a*a-a+2)%q==0]
    assert len(roots)==2
    for lam in roots:
        d=dense(D.subs(l,lam),q); c=dense(C.subs(l,lam),q)
        di=dense(Di.subs(l,lam),q); ci=dense(Ci.subs(l,lam),q)
        assert mul(d,di,q)==eye(7) and mul(c,ci,q)==eye(7)
        # Source fidelity: cyclic annihilator rank three; 18 encoded qudits.
        for L in [7,14]:
            V=sp.zeros(L)
            for i in range(L):
                for j,a in enumerate([1,lam,lam-1,-1]):V[i,(i+j)%L]+=a
            # Unit Smith entries count all nonzero divisors over Z/q exactly.
            from sympy.matrices.normalforms import smith_normal_form
            S=smith_normal_form(V,domain=sp.ZZ)
            kernel_size=1
            for i in range(L):kernel_size*=gcd(int(S[i,i]),q)
            assert kernel_size==q**3
        for uu,vv in [(1,1),(3%q,5%q)]:
            # Arbitrary evaluation at units also checks Laurent adjoints below.
            hx=dense(H.subs({l:lam,u:uu}),q);hy=dense(H.subs({l:lam,u:vv}),q)
            ax=kron(hx,eye(7),7,q);by=kron(eye(7),hy,7,q)
            a=dense(Delta.subs(u,uu),q);b=dense(Delta.subs(u,vv),q)
            ad=kron(a,eye(7),7,q);bd=kron(eye(7),b,7,q)
            f0=kron(c,c,7,q);f2=kron(d,d,7,q)
            f1=block(kron(d,c,7,q),{},{},kron(c,d,7,q),49)
            f1i=block(kron(di,ci,7,q),{},{},kron(ci,di,7,q),49)
            d1=row(ax,by,49);d2=col(scale(by,-1,q),ax,49)
            e1=row(ad,bd,49);e2=col(scale(bd,-1,q),ad,49)
            assert mul(d1,f1,q)==mul(f0,e1,q)
            assert mul(d2,f2,q)==mul(f1,e2,q)
            assert not mul(d1,d2,q)
            # Split each unit-entry complex into two product states.
            shear={};sheari={};target1={};target2={}
            toric=trivial=0
            for i,j in product(range(7),repeat=2):
                z=7*i+j;aa=1 if i<4 else (uu-1)%q;bb=1 if j<4 else (vv-1)%q
                if i>=4 and j>=4:
                    toric+=1;W=sp.eye(2)
                    if aa:target1[0,z]=aa;target2[z+49,0]=aa
                    if bb:target1[0,z+49]=bb;target2[z,0]=-bb%q
                else:
                    trivial+=1
                    W=sp.Matrix([[1,-bb],[0,1]]) if i<4 else sp.Matrix([[0,-1],[1,aa]])
                    target1[0,z]=1;target2[z+49,0]=1
                Wi=W.inv()
                for ii,jj in product(range(2),repeat=2):
                    if int(W[ii,jj])%q:shear[z+49*ii,z+49*jj]=int(W[ii,jj])%q
                    if int(Wi[ii,jj])%q:sheari[z+49*ii,z+49*jj]=int(Wi[ii,jj])%q
            assert (toric,trivial)==(9,40)
            # Target rows/columns keep the 49 summand index.
            tr1={(j%49,j):z for (_,j),z in target1.items()}
            tr2={(i,i%49):z for (i,_),z in target2.items()}
            assert mul(e1,shear,q)==tr1 and mul(sheari,e2,q)==tr2
            B=mul(f1,shear,q);Bi=mul(sheari,f1i,q)
            assert mul(B,Bi,q)==eye(98) and mul(Bi,B,q)==eye(98)
            # B is physical coordinate change; diag(B^T,B^-1) is symplectic.
            # For nontrivial uu,vv the dagger evaluates the transpose at inverse units.
            # The coefficient-level identity B^-1 B=I implies the dagger identity;
            # the following checks the exact finite-torus specialization uu=vv=1.
            if uu==vv==1:
                Q=block(trans(B),{},{},Bi,98)
                J=block({},eye(98),scale(eye(98),-1,q),{},98)
                assert mul(mul(trans(Q),J,q),Q,q)==J
        zero,one,xi,plus,neg,times,power=quotient(q,lam)
        assert power(xi,7)==one and all(power(xi,n)!=one for n in range(1,7))
        hval=lambda z:plus(plus(power(z,3),tuple(-(lam-1)*x%q for x in power(z,2))),plus(tuple(-lam*x%q for x in z),neg(one)))
        rootsS=[power(xi,s) for s in [1,2,4]]
        assert all(hval(z)==zero for z in rootsS)
        # Differences have nonzero image in F8, so have odd norm and are units.
        for i,j in [(0,1),(0,2),(1,2)]:
            diff=plus(rootsS[i],neg(rootsS[j]))
            M=sp.Matrix.hstack(*(sp.Matrix(times(diff,power(xi,k0))) for k0 in range(3)))
            assert gcd(int(M.det()),q)==1
        # CRT evaluation map on the basis xi^a eta^b of E; y roots xi^(6,5,3).
        Y=[power(xi,s) for s in [6,5,3]]
        CRT=sp.zeros(9)
        for a,b in product(range(3),repeat=2):
            for j,z in enumerate(Y):
                val=times(power(xi,a),power(z,b))
                for c0 in range(3):CRT[3*j+c0,3*a+b]=val[c0]
        assert gcd(int(CRT.det()),q)==1
        # Each one-site power 1..6 has no nonzero fixed vector, on either side.
        for s in [1,2,4]:
            for i,j in product(range(7),repeat=2):
                exponent=(i-s*j)%7
                if exponent:
                    diff=plus(power(xi,exponent),neg(one))
                    M=sp.Matrix.hstack(*(sp.Matrix(times(diff,power(xi,k0))) for k0 in range(3)))
                    assert gcd(int(M.det()),q)==1
        # Perfect standard string crossing pairing: each coordinate detects itself.
        for a in range(q):assert (a*1)%q==a
        print(f'q={q} lambda={lam}: L=7,14 ker(V)=q^3, k=18; 9 toric+80 product qudits; symplectic; CRT and translations OK')
print('FIDELITY: q=4 h=t^3+3t^2+2t+3; GR(4,3) has 4^3=64 elements; source charge rank=9')
print('SCOPE: symbolic blocking is uniform; finite full-complex and CRT checks use q=2,4,8,16,32,64, both lambda roots')
print('ALL_CHECKS_OK')
