#!/usr/bin/env python3
"""Fiol, Martinez-Montoya, Rios Fukelman, "The planar limit of N=2 superconformal quiver theories",
arXiv:2006.06379v2; JHEP 08 (2020) 161. Section 2.2, Conjecture 1: "For any tree with e edges, any fixed positive
integers n_1,...,n_e and arbitrary rho(k_i,n_i)>0 the polynomials P(tau,x,n_1,...,n_e) have all their x roots on the
unit circle", with P = sum_{k_1..k_e} rho(k_i,n_i) sum_states tau^{e_pm} prod_{up vertices} x^{m_j} (eq. polbeforesum).

Counterexample: path L-C-R, n_1 = n_2 = 4, splits a, c in {1,2,3}: m_L = a, m_C = 8-a-c, m_R = c; tau = -1;
weights W_M = [[1,1,1],[1,10,1],[1,1,M]] (rows a, columns c), M >= 31.
Checks:
 (1) exact (sympy): P_M = (M+17) - 6x + (M-23)x^2 - 2Mx^3 + 24x^4 - 2Mx^5 + (M-23)x^6 - 6x^7 + (M+17)x^8;
 (2) exact: P_M(x) = (1-x)^2 x^3 Q_M(x + 1/x), Q_M(y) = (M+17)y^3 + (2M+28)y^2 + (M-35)y - 52;
 (3) exact: disc(Q_M) = -24M^3 - 5348M^2 + 85664M + 2936364 < 0 for M >= 31 (as -24s^3-7580s^2-315104s-262464, s = M-31);
 (4) exact: for M = 31, q = Q_31/2 = 24y^3+45y^2-2y-26 < 0 on y <= 0 (the certificate identity and interval bound);
 (5) exact: disc(q) = -16404 < 0, so q (and Q_31) has a nonreal root; numerically P_31 has a root with | |x|-1 | ~ 0.019;
 (6) interior tau: M = 60, tau = -19/20 gives 400P = the stated octic and a nonreal root of R(y);
 (7) numeric: product weights rho = u_a v_c on the same path keep every root on the unit circle (2000 random trials)."""
import sys, itertools, random
import sympy as sp, numpy as np
bad=0
def rep(l,ok):
    global bad; print(l,ok); bad+=0 if ok else 1
x,y,M,s,tau=sp.symbols('x y M s tau')
def P_of(W,t):
    P=0
    for a in (1,2,3):
        for c in (1,2,3):
            m=(a,8-a-c,c)
            for sp_ in itertools.product((0,1),repeat=3):
                e=(sp_[0]!=sp_[1])+(sp_[1]!=sp_[2])
                P+=W[a-1][c-1]*t**e*x**(sum(si*mi for si,mi in zip(sp_,m)))
    return sp.expand(P)
W=[[1,1,1],[1,10,1],[1,1,M]]
PM=P_of(W,-1)
target=(M+17)-6*x+(M-23)*x**2-2*M*x**3+24*x**4-2*M*x**5+(M-23)*x**6-6*x**7+(M+17)*x**8
rep('(1) P_M formula',sp.expand(PM-target)==0)
QM=(M+17)*y**3+(2*M+28)*y**2+(M-35)*y-52
rep('(2) reciprocal reduction',sp.expand(PM-((1-x)**2*x**3*QM.subs(y,x+1/x)))==0 or sp.simplify(PM-((1-x)**2*x**3*QM.subs(y,x+1/x)))==0)
D=sp.discriminant(QM,y)
rep('(3) discriminant formula',sp.expand(D-(-24*M**3-5348*M**2+85664*M+2936364))==0 and sp.expand(D.subs(M,s+31)-(-24*s**3-7580*s**2-315104*s-262464))==0)
q=24*y**3+45*y**2-2*y-26; z=sp.symbols('z')
ok=sp.expand(q.subs(y,z-sp.Rational(5,4))-(24*z**2*(z-sp.Rational(3,4))-27*(z+sp.Rational(1,27))**2-sp.Rational(11,432)))==0
ok&=sp.expand(QM.subs(M,31)-2*q)==0
rep('(4) certificate identity for q and Q_31 = 2q',ok)
r=np.roots([float(c) for c in sp.Poly(PM.subs(M,31),x).all_coeffs()])
dev=max(abs(abs(r)-1)); rep(f'(5) disc(q) = {sp.discriminant(q,y)} < 0 (q has a nonreal root); P_31 numeric max | |x|-1 | = {dev:.4f}',sp.discriminant(q,y)==-16404 and dev>1e-3)
P60=P_of([[1,1,1],[1,10,1],[1,1,60]],sp.Rational(-19,20))
oct_=sp.Rational(1,400)*(30800-2280*x+12901*x**2-45676*x**3+8664*x**4-45676*x**5+12901*x**6-2280*x**7+30800*x**8)
R=30800*y**4-2280*y**3-110299*y**2-38836*y+44462
ok=sp.expand(P60-oct_)==0 and sp.discriminant(R,y)<0
rr=np.roots([float(c) for c in sp.Poly(P60,x).all_coeffs()]); dev60=max(abs(abs(rr)-1)); print('   P60 numeric max | |x|-1 | =',dev60); ok&=dev60>1e-4
rep('(6) interior tau=-19/20, M=60',ok)
random.seed(3); np.random.seed(3); worst=0
for _ in range(2000):
    u=np.exp(np.random.randn(3)*2); v=np.exp(np.random.randn(3)*2); t=random.choice([-1.0,1.0,random.uniform(-1,1)])
    Pn=np.zeros(9)
    for a in (1,2,3):
        for c in (1,2,3):
            m=(a,8-a-c,c)
            for sp_ in itertools.product((0,1),repeat=3):
                e=(sp_[0]!=sp_[1])+(sp_[1]!=sp_[2]); Pn[sum(si*mi for si,mi in zip(sp_,m))]+=u[a-1]*v[c-1]*t**e
    worst=max(worst,max(abs(abs(np.roots(Pn[::-1]))-1)))
rep(f'(7) product weights stay on the circle, max dev {worst:.2e}',worst<1e-4)
print('bad=',bad); sys.exit(1 if bad else 0)
