#!/usr/bin/env python3
"""Exact fidelity and counterexample checks for arXiv:1904.05732v1.

Requires Python >= 3.12 and SymPy 1.14. Mathematical assertions use only
integers, rationals and exact radicals. Floating-point spectral radii are
reported separately as a requested diagnostic, never as proof evidence.
"""
from fractions import Fraction as F
import cmath
import hashlib
from pathlib import Path
import sys
import sympy as sp

assert sys.version_info >= (3, 12)
w = sp.symbols('omega', real=True)
z = sp.symbols('z')


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def same(A, B):
    assert A.shape == B.shape
    for entry in A - B:
        zero(entry)


def outer(a):
    return a * a.T / (a.T * a)[0]


def fidelity():
    # Section 3.1, Example 1: alpha = pi/3; two projections are averaged.
    a = sp.Matrix([-sp.sqrt(3) / 2, sp.Rational(1, 2)])
    b = sp.Matrix([0, 1])
    B = sp.eye(2) - w * (outer(a) + outer(b)) / 2
    printed = sp.Matrix([[1 - 3*w/8, sp.sqrt(3)*w/8],
                         [sp.sqrt(3)*w/8, 1 - 5*w/8]])
    same(B, printed)
    zero(B.charpoly(z).as_expr() - (z - (1-w/4))*(z - (1-3*w/4)))
    endpoint = sp.Rational(8, 3)
    assert (1 - 3*w/4).subs(w, endpoint) == -1
    assert set(B.subs(w, 2).eigenvals()) == {-sp.Rational(1, 2), sp.Rational(1, 2)}
    print('FIDELITY: source Example 1 averaged matrix reproduced exactly at alpha=pi/3')
    print('FIDELITY: eigenvalues=1-omega/4,1-3*omega/4; Omega=8/3; omega_opt=2; rho_opt=1/2')


def exact_family(m):
    n = 2*m + 1
    eye = sp.eye(n)
    e0 = eye[:, 0]
    us = [eye[:, i+1] for i in range(m)]
    vs = [eye[:, m+i+1] for i in range(m)]
    rows = [e0] + [4*e0 + sign*7*us[i] + 4*vs[i]
                    for i in range(m) for sign in (1, -1)]
    A = sp.Matrix.vstack(*(a.T for a in rows))
    assert abs(A.det()) == 56**m
    assert [(a.T*a)[0] for a in rows] == [1] + [81]*(2*m)
    R = eye - w*outer(e0)
    # Direct source update: update root, then leaf, then pool, same omega.
    B = sum(((eye - w*outer(a))*R for a in rows[1:]), sp.zeros(n)) / (2*m)
    M = sum((outer(a) for a in rows[1:]), sp.zeros(n)) / (2*m)
    same(B, (eye - w*M)*R)
    for u in us:
        same(M*u, sp.Rational(49, 81*m)*u)
    s = sum(vs, sp.zeros(n, 1))
    same(M*e0, sp.Rational(16, 81)*e0 + sp.Rational(16, 81*m)*s)
    for v in vs:
        same(M*v, sp.Rational(16, 81*m)*(e0+v))
    lu = 1 - sp.Rational(49, 81*m)*w
    lv = 1 - sp.Rational(16, 81*m)*w
    C = sp.Matrix([[(1-16*w/81)*(1-w), -16*w/81],
                   [-16*w*(1-w)/(81*m), 1-16*w/(81*m)]])
    D = sp.diag(*([lu]*m + [lv]*(m-1)), C)
    P = sp.Matrix.hstack(*(us + [vs[i]-vs[-1] for i in range(m-1)] + [e0, s]))
    assert P.det() != 0
    same(B*P, P*D)
    alpha = sp.Rational(16*(m+1), 81*m)
    t = 2-(1+alpha)*w+16*w*w/81
    d = (1-w)*(1-alpha*w)
    zero(C.trace()-t)
    zero(C.det()-d)
    zero(B.charpoly(z).as_expr() - (z-lu)**m*(z-lv)**(m-1)*(z*z-t*z+d))
    endpoint = 1 + 1/alpha
    assert endpoint > 4
    assert endpoint < sp.Rational(162*m, 49) < sp.Rational(162*m, 16)
    assert d.subs(w, endpoint) == 1
    assert d.subs(w, endpoint+sp.Rational(1, 10000)) > 1
    # The affine source iteration has B as its error map for an arbitrary x*.
    xstar = sp.Matrix(sp.symbols('x0:'+str(n)))
    x = sp.Matrix(sp.symbols('y0:'+str(n)))
    rhs = A*xstar
    rootx = x+w*(rhs[0]-(rows[0].T*x)[0])/(rows[0].T*rows[0])[0]*rows[0]
    pooled = sum((rootx+w*(rhs[j]-(rows[j].T*rootx)[0])/81*rows[j]
                  for j in range(1,n)), sp.zeros(n,1))/(2*m)
    same(pooled-xstar, B*(x-xstar))
    print(f'EXACT m={m}: dimension={n}; abs(det A)={abs(A.det())}; det(P)={P.det()}; Omega={endpoint}')
    print(f'EXACT m={m}: full symbolic characteristic polynomial and similarity verified')
    if m == 2:
        fixed = sp.Rational(17, 4)
        assert lu.subs(w,fixed) == -sp.Rational(185,648)
        assert lv.subs(w,fixed) == sp.Rational(47,81)
        assert t.subs(w,fixed) == sp.Rational(19,324)
        assert d.subs(w,fixed) == sp.Rational(91,108)
        assert (t*t-4*d).subs(w,fixed) < 0
        print('EXACT m=2 omega=17/4: scalar eigenvalues=-185/648 (twice),47/81; block trace=19/324,det=91/108')
    return F(int(sp.numer(endpoint)), int(sp.denom(endpoint)))


def symbolic_family():
    m = sp.symbols('m', integer=True, positive=True)
    alpha = 16*(m+1)/(81*m)
    C = sp.Matrix([[(1-16*w/81)*(1-w), -16*w/81],
                   [-16*w*(1-w)/(81*m), 1-16*w/(81*m)]])
    t = 2-(1+alpha)*w+16*w*w/81
    d = (1-w)*(1-alpha*w)
    zero(C.trace()-t)
    zero(C.det()-d)
    a, b, c = 1-t+d, 1+t+d, 1-d
    zero(a-16*w*w/(81*m))
    zero(b-(4-2*(1+alpha)*w+16*(2*m+1)*w*w/(81*m)))
    zero(c-w*(1+alpha-alpha*w))
    discriminant = sp.discriminant(b,w)
    zero(discriminant+4*(7*m+16)*(137*m-16)/(6561*m*m))
    endpoint = 1+1/alpha
    zero(endpoint-(sp.Rational(97,16)-sp.Rational(81,16)/(m+1)))
    zero(162*m/sp.Integer(49)-endpoint-(2592*m*m-2161*m-784)/(784*(m+1)))
    zero((2592*m*m-2161*m-784)-(2592*m*(m-2)+3023*m-784))
    zero((endpoint-4)-(33*m-48)/(16*(m+1)))
    assert sp.limit(endpoint,m,sp.oo) == sp.Rational(97,16)
    print('SYMBOLIC all-m: trace, determinant, three Jury expressions, negative discriminant and scalar-threshold gap identities verified')
    print('SYMBOLIC all-m: Omega-4=(33*m-48)/(16*(m+1)); limit Omega=97/16')


def grid(m, endpoint):
    # Exact Schur/Jury and scalar inequalities at 4095 rational grid points.
    alpha = F(16*(m+1),81*m)
    maximum_rho = 0.0
    for k in range(1,4096):
        q = endpoint*F(k,4096)
        t = 2-(1+alpha)*q+F(16,81)*q*q
        d = (1-q)*(1-alpha*q)
        assert 1-t+d > 0 and 1+t+d > 0 and 1-d > 0
        lu,lv = 1-F(49,81*m)*q,1-F(16,81*m)*q
        assert -1 < lu < 1 and -1 < lv < 1
        # Requested numerical spectral-radius diagnostic; exact checks above
        # carry the evidence, not these floating-point values.
        tf,df = float(t),float(d)
        delta = cmath.sqrt(tf*tf-4*df)
        rho = max(abs((tf+delta)/2),abs((tf-delta)/2),abs(float(lu)),abs(float(lv)))
        assert rho < 1.0
        maximum_rho = max(maximum_rho,rho)
    q = endpoint+F(1,10000)
    t,d = 2-(1+alpha)*q+F(16,81)*q*q,(1-q)*(1-alpha*q)
    assert d > 1
    tf,df=float(t),float(d)
    delta=cmath.sqrt(tf*tf-4*df)
    rho=max(abs((tf+delta)/2),abs((tf-delta)/2),abs(float(1-F(49,81*m)*q)),abs(float(1-F(16,81*m)*q)))
    assert rho >= 1.0
    print(f'EXACT GRID m={m}: 4095 rational interior points pass; at Omega+1/10000 determinant>1')
    print(f'NUMERIC DIAGNOSTIC ONLY m={m}: max interior rho={maximum_rho:.12g}; rho(Omega+1/10000)={rho:.12g}')


def main():
    print(f'Python={sys.version_info.major}.{sys.version_info.minor}; SymPy={sp.__version__}')
    print('All mathematical evidence is exact; decimal spectral radii are diagnostics only.')
    fidelity()
    symbolic_family()
    for m in (2,3,4):
        grid(m,exact_family(m))
    print('check.py sha256='+hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
