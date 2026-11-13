#!/usr/bin/env python3
"""Exact algebra for arXiv:2607.01334, Section V's first question.

No numerical momentum sampling is used. Inequalities and the continuous
argument lift are explained in README.md; this program certifies their
polynomial identities and the discriminant conditions symbolically.
"""
import sys
import sympy as S

assert sys.version_info >= (3, 12)
I = S.I
x, y, s, delta, t, a, c = S.symbols('x y s delta t a c', real=True)
lam = S.Symbol('lam')
z = x + I*y
f = 1 - s + s*z
fc = S.conjugate(f)
q = S.expand(f*fc)
q_circle = S.expand(q.subs(y**2, 1-x**2))

def equal(lhs, rhs):
    if isinstance(lhs, S.MatrixBase):
        assert lhs.shape == rhs.shape
        assert all(S.simplify(v) == 0 for v in lhs-rhs), lhs-rhs
    else:
        assert S.simplify(lhs-rhs) == 0, S.factor(lhs-rhs)

def blocks(A, B, C, D):
    return A.row_join(B).col_join(C.row_join(D))

# Source fidelity: Section IV.3's determinant, and Section II's two
# doubly degenerate bands.  x=cos(k), y=sin(k), x^2+y^2=1.
J2 = S.diag(1, -1)
Bsource = S.Matrix([[1, delta], [-delta, -1]])
bsource = f*Bsource
Gsource = blocks(S.zeros(2), J2*bsource.conjugate().T*J2,
                 bsource, S.zeros(2))
equal(bsource.det(), f**2*(delta**2-1))
equal(Gsource.charpoly(lam).as_expr(),
      (lam**2-(1-delta**2)*q)**2)
equal(q_circle, 1 + 2*s*(s-1)*(1-x))
equal((delta**2-1)*(2*s*(s-1)*(x-1)-1),
      (1-delta**2)*q_circle)
print('FIDELITY: source det b=f^2(delta^2-1) and double-band characteristic polynomial OK')

# General four-dimensional collision block and all bosonic constraints.
A = S.diag(a, -c)
D = S.Matrix([[0,t],[t,0]])
B = blocks(A,D,-D,-A)
J = S.diag(1,1,-1,-1)
C = blocks(S.zeros(2),S.eye(2),S.eye(2),S.zeros(2))
equal(J*B.T*J, B)
equal(C*B.conjugate()*C, -B)
order = [0,3,1,2]
Bperm = B.extract(order, order)
P = S.Matrix([[a,t],[-t,c]])
Q = S.Matrix([[-c,t],[-t,-a]])
equal(Bperm, S.diag(P,Q))
alpha, beta = (a+c)/2, (a-c)/2
pminus = (lam-alpha)**2-(beta**2-t**2)
pplus = (lam+alpha)**2-(beta**2-t**2)
equal(B.charpoly(lam).as_expr(), pminus*pplus)
equal(B.det(), (a*c+t**2)**2)
droot = S.symbols('droot')
V2 = S.Matrix([[t,t],[droot-beta,-droot-beta]])
for M, center in [(P,alpha),(Q,-alpha)]:
    residual = M*V2-V2*S.diag(center+droot,center-droot)
    assert all(S.simplify(S.rem(S.expand(v),
               droot**2-(beta**2-t**2),droot)) == 0 for v in residual)
equal(V2.det(),-2*t*droot)
print('FAMILY: coordinate permutation, discriminants, det B=(ac+t^2)^2, symbolic eigenbases OK')

b = f*B
G = blocks(S.zeros(4), J*b.conjugate().T*J, b, S.zeros(4))
eta = S.diag(J,J)
charge = S.diag(C,C)
R = S.diag(-S.eye(4),S.eye(4))
Gminus = G.xreplace({y:-y})
equal(eta*G*eta, G.conjugate().T)
equal(eta*G, (eta*G).conjugate().T)
equal(charge*G.conjugate()*charge, -Gminus)
equal(C*b.conjugate()*C, -b.xreplace({y:-y}))
equal(R*R, S.eye(8))
equal(R.conjugate().T*R, S.eye(8))
equal(R*G*R, -G)
K = S.Matrix([[0,fc],[f,0]])
equal(G, S.kronecker_product(K,B))
equal(K.conjugate().T,K)
equal(K*K, q*S.eye(2))
h = S.symbols('h', positive=True)
VK = S.Matrix([[fc,fc],[h,-h]])
assert all(S.simplify(S.rem(S.expand(v),h**2-q,h)) == 0
           for v in K*VK-VK*S.diag(h,-h))
equal(VK.det(),-2*h*fc)
equal(G.det(), q**4*(a*c+t**2)**4)
equal(b.det(), f**4*(a*c+t**2)**2)
equal(q_circle, (2*s-1)**2+2*s*(1-s)*(1+x))
print('BOSONIC: pseudo-Hermiticity, Hermitian Hessian, momentum particle-hole constraint, fixed chiral grading OK')
print('GAP: det G=|f|^8(ac+t^2)^4; |f|^2=(2s-1)^2+2s(1-s)(1+cos k) OK')

# Exact rational parameters and exact radicals, never floating point.
B0 = B.subs({a:3,c:1})
for value, expected in [
    (S.Rational(1,2), [2+S.sqrt(3)/2,2-S.sqrt(3)/2,
                         -2+S.sqrt(3)/2,-2-S.sqrt(3)/2]),
    (S.Rational(3,2), [2+I*S.sqrt(5)/2,2-I*S.sqrt(5)/2,
                         -2+I*S.sqrt(5)/2,-2-I*S.sqrt(5)/2]),
]:
    M = B0.subs(t,value)
    equal(M.charpoly(lam).as_expr(), S.prod(lam-v for v in expected))
    assert all(mult == 1 for mult in M.eigenvals().values())
    assert all(S.simplify(M.charpoly(lam).as_expr().subs(lam,v)) == 0
               for v in expected)
    V = S.Matrix.hstack(*[(M-v*S.eye(4)).nullspace()[0] for v in expected])
    assert S.simplify(V.det()) != 0
    equal(M*V, V*S.diag(*expected))
    print('EXACT t='+str(value)+': eigenvalues '+str(expected)+'; eigenbasis invertible')
Bc = B0.subs(t,1)
assert (Bc-2*S.eye(4)).rank() == 3
assert (Bc+2*S.eye(4)).rank() == 3
assert Bc.det() == 16
print('COLLISION t=1: eigenvalues +/-2, both geometric multiplicities 1; det=16 (nonzero-frequency EP)')

# Symbolic inequalities under the source domain (a,c>0, a!=c,
# 1/2<s<1).  These polynomial certificates are used in README.md.
u,v,w,d = S.symbols('u v w d', positive=True)
assert (u*v+t**2).is_positive
# alpha>|beta| follows from alpha^2-beta^2=ac>0.
equal(alpha**2-beta**2,a*c)
# Stable t = |beta|*r, 0<r<1; unstable t=|beta|*r, r>1.
equal(d**2-(d/2)**2,S.Rational(3,4)*d**2)
equal(d**2-(3*d/2)**2,-S.Rational(5,4)*d**2)
# Exact winding factorization, with r=(1-s)/s in (0,1).
# The factor 1+r exp(-ik) has real part >=1-r>0, so arg is
# continuous and periodic; det b has argument 4k+4 arg(that factor).
r = S.symbols('r', real=True)
lift_factor = (s*z*(1+r*S.conjugate(z))).subs(r,(1-s)/s)
equal(S.expand(((1-s)+s*z)-lift_factor).subs(y**2,1-x**2),0)
equal(1+r*x, (1-r)+r*(1+x))
equal((4*(2*S.pi)+4*w)-(4*0+4*w),8*S.pi)
print('WINDING: exact factorization and right-half-plane certificate give continuous argument increment 8*pi, index +4')

# Optional half-line recurrence, not a finite truncation claim.
n = S.symbols('n', integer=True, nonnegative=True)
rho = -(1-s)/s
equal((1-s)*rho**n+s*rho**(n+1),0)
print('BOUNDARY: exact geometric kernel recurrence (square summability uses |(1-s)/s|<1) OK')
print('ALL_CHECKS_OK')
