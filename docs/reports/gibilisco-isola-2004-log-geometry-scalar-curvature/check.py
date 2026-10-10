#!/usr/bin/env python3
"""Gibilisco–Isola math-ph/0407007v2, Conjecture 4.1, p = +infinity (log embedding).
(1) Intrinsic scalar curvature of the pull-back metric of rho -> log(rho) on the simplex, computed from the
    metric in coordinates (x_1..x_{n-1}), x_n = 1 - sum, via Christoffel symbols (sympy, exact rational points),
    for n = 3, 4, 5, compared with R = 2 (e_3 - 4 e_4) / s_2^2.
(2) Sanity: the same code for p = 2 (map 2 sqrt(rho)) gives (n-1)(n-2)/4 (source, Section 4).
(3) The identity K = sum z_i^2 (1-2z_i)^2 + u(3-4u)B + u^2(1-u) + v(4B - 4u^2 + 6u - 2) with K = (A - 4E_2)T + 4N,
    symbolically for n - 2 = 1..4 remaining coordinates.
(4) Strict Schur monotonicity on random positive vectors under random T-transforms (n = 3..8)."""
import sys, random, itertools
import sympy as sp
from fractions import Fraction as Fr
def scal(n, pt, embed):
    xs=sp.symbols(f'x1:{n}'); xn=1-sum(xs); X=list(xs)+[xn]
    F=[embed(v) for v in X]
    J=sp.Matrix([[sp.diff(f,xi) for xi in xs] for f in F])
    g=sp.simplify(J.T*J); m=n-1
    gi=g.inv()
    Gam=[[[sum(gi[k,l]*(sp.diff(g[l,i],xs[j])+sp.diff(g[l,j],xs[i])-sp.diff(g[i,j],xs[l]))/2 for l in range(m)) for j in range(m)] for i in range(m)] for k in range(m)]
    sub={xs[i]:pt[i] for i in range(m)}
    Gs=[[[sp.nsimplify(Gam[k][i][j].subs(sub)) for j in range(m)] for i in range(m)] for k in range(m)]
    dG=[[[[sp.nsimplify(sp.diff(Gam[k][i][j],xs[l]).subs(sub)) for l in range(m)] for j in range(m)] for i in range(m)] for k in range(m)]
    # Ricci R_ij = d_k G^k_ij - d_j G^k_ik + G^k_kl G^l_ij - G^k_jl G^l_ik
    Ric=sp.zeros(m,m)
    for i in range(m):
        for j in range(m):
            Ric[i,j]=sum(dG[k][i][j][k]-dG[k][i][k][j]+sum(Gs[k][k][l]*Gs[l][i][j]-Gs[k][j][l]*Gs[l][i][k] for l in range(m)) for k in range(m))
    gis=gi.subs(sub)
    return sp.nsimplify(sum(gis[i,j]*Ric[i,j] for i in range(m) for j in range(m)))
bad=0
for n,pt in [(3,[Fr(1,5),Fr(3,10)]),(4,[Fr(1,10),Fr(1,5),Fr(3,10)]),(5,[Fr(1,10),Fr(3,20),Fr(1,5),Fr(1,4)])]:
    pts=[sp.Rational(p.numerator,p.denominator) for p in pt]; full=pts+[1-sum(pts)]
    R=scal(n,pts,lambda v: sp.log(v))
    e=lambda k: sum(sp.prod(c) for c in itertools.combinations(full,k))
    s2=sum(v**2 for v in full); pred=2*(e(3)-4*e(4))/s2**2
    print('n',n,'R',R,'pred',sp.nsimplify(pred),'eq',sp.simplify(R-pred)==0); bad+= sp.simplify(R-pred)!=0
    if n<=4:
        R2=scal(n,pts,lambda v: 2*sp.sqrt(v)); print('  p=2 check',sp.nsimplify(R2),'expected',sp.Rational((n-1)*(n-2),4)); bad+= sp.simplify(R2-sp.Rational((n-1)*(n-2),4))!=0
# (3) identity
for r in range(1,5):
    z=sp.symbols(f'z1:{r+1}'); u,v=sp.symbols('u v')
    a_b=[u,v]  # a+b=u, ab=v
    A=sum(z); B=sum(t**2 for t in z)
    E=lambda k: sum(sp.prod(c) for c in itertools.combinations(z,k)) if k<=r else 0
    E2,E3,E4=E(2),E(3),E(4)
    N=v*A+u*E2+E3-4*(v*E2+u*E3+E4); T=u**2-2*v+B
    K=(A-4*E2)*T+4*N
    rhs=sum(t**2*(1-2*t)**2 for t in z)+u*(3-4*u)*B+u**2*(1-u)+v*(4*B-4*u**2+6*u-2)
    diff=sp.expand((K-rhs).subs(u,1-A))
    print('identity r=',r, diff==0); bad+= diff!=0
# (4) monotonicity
random.seed(1)
def R(x):
    s2=sum(t*t for t in x); e3=sum(a*b*c for a,b,c in itertools.combinations(x,3)); e4=sum(a*b*c*d for a,b,c,d in itertools.combinations(x,4))
    return 2*(e3-4*e4)/s2**2
viol=0
for _ in range(20000):
    n=random.randint(3,8); x=[random.random()**3+1e-3 for _ in range(n)]; s=sum(x); x=[t/s for t in x]
    i,j=random.sample(range(n),2); lam=random.uniform(0.01,0.99)
    y=list(x); y[i]=lam*x[i]+(1-lam)*x[j]; y[j]=(1-lam)*x[i]+lam*x[j]
    if abs(x[i]-x[j])>1e-6 and not R(y)>R(x): viol+=1
print('T-transform violations',viol); bad+=viol
print('bad=',bad); sys.exit(1 if bad else 0)
