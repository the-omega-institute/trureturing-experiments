#!/usr/bin/env python3
"""Exact checks for arXiv:2605.31491v1, Conjecture 8.7.

Independent implementation; Python 3 and SymPy. No floating-point evidence.
Source conventions: Borovik, de Korte, Meurrens, Pavlov, Section 8,
and their CC-BY-4.0 supplement https://doi.org/10.5281/zenodo.22208208.
No supplementary code is required to run this program.
"""

from fractions import Fraction
from itertools import product
import sympy as S


def zero(expression):
    value = S.cancel(expression)
    assert value == 0, value


def monomials(spins):
    """Nonnegative integer solutions of the source's three spin equations."""
    out = []
    for h12, h23, h13 in product(range(max(spins) + 1), repeat=3):
        v = (spins[0] - h12 - h13,
             spins[1] - h12 - h23,
             spins[2] - h23 - h13)
        if min(v) >= 0:
            out.append((*v, h12, h23, h13))
    return out


def source_count(spins):
    """Proposition 3.9 / equation (24), ordered spins."""
    a, b, c = sorted(spins)
    p = max(0, a + b - c)
    return (Fraction((a+1)*(a+2)*(3*b-a+3), 6)
            - Fraction(p*(p+2)*(2*p+5), 24)
            - Fraction(1-(-1)**p, 16))


def fidelity():
    # The two rows quoted are the source's formal monomial counts, before
    # Gram relations. In particular 112 (equivalently 211) has five, not ten.
    spins = [(0,0,1),(0,1,1),(1,1,1),(0,0,2),(0,1,2),(1,1,2),
             (0,0,3),(0,1,3),(1,1,3),(0,2,3),(1,2,3),(2,2,2),
             (2,2,3),(0,3,3),(1,3,3),(2,3,3),(3,3,3)]
    counts = [1,2,4,1,2,5,1,2,5,3,8,11,13,4,10,17,23]
    for s, n in zip(spins, counts):
        assert len(monomials(s)) == source_count(s) == n
    for s in product(range(7), repeat=3):
        assert len(monomials(s)) == source_count(s)
    print("FIDELITY: all 17 monomial counts in source Table 3 reproduced exactly")
    print("FIDELITY: Proposition 3.9 agrees for all 343 spins in [0,6]^3")


# Particle 1 varies; the four fixed vectors are P2,P3,Z2,Z3. Their Gram
# matrix uses only constraints on these fixed vectors. No differentiated
# particle-1 cone constraint is substituted before differentiating.
x, y, alpha, beta = X = S.symbols("x y alpha beta")
u, v, w, h = U = S.symbols("u v w h")
c, A, B, L = S.symbols("c A B L")
r = S.symbols("r", integer=True, positive=True)
G = S.Matrix([[0,c,0,B], [c,0,A,0], [0,A,0,L], [B,0,L,0]])
V1 = (u*y-v*x)/c
V2 = (A*x-alpha*c)/y
V3 = (beta*c-B*y)/x
H12 = -2*(w*x-u*alpha)
H13 = -2*(h*y-v*beta)
J = H23 = -2*(L*c-A*B)
T = V2*V3
basis = [V1**2*V2*V3, V1*H12*V3, V1*H13*V2,
         H12*H13, V1**2*H23]


def ppartial(q, i):
    # K = x^(-r-1) c^(-1) y^(-r-1). Exact K^-1 ∂(Kq).
    logterm = -(r+1)/X[i] if i < 2 else 0
    return S.cancel(S.diff(q, X[i]) + logterm*q)


def mixed(q):
    return S.cancel(sum(G[i,j]*ppartial(S.diff(q,U[j]), i)
                        for i,j in product(range(4), repeat=2)))


def lap_z(q):
    return S.cancel(sum(G[i,j]*S.diff(q,U[i],U[j])
                        for i,j in product(range(4), repeat=2)))


def lap_p(q):
    return S.cancel(sum(G[i,j]*ppartial(ppartial(q,j),i)
                        for i,j in product(range(4), repeat=2)))


def twice_diverged(q):
    # On quadratic Z1 inputs, the d=3 Thomas definition implies
    # C1^2 = (3/4)(∂P1·∂Z1)^2 - (1/4)□P1□Z1.
    return S.cancel(S.Rational(3,4)*mixed(mixed(q))
                    - S.Rational(1,4)*lap_p(lap_z(q)))


a = (2*r-1)*(r+1)
e = r*(2*r-3)
b = -2*(r+1)**2
k = 2*(r*r-3*r+1)
expected = [a*T+J, -a*T-J, -a*T-J, b*T+k*J, e*J]


def family():
    for i, (q, target) in enumerate(zip(basis, expected), 1):
        zero(twice_diverged(q)-target)
        print(f"FAMILY basis {i}: symbolic-r rational-function residual = 0", flush=True)
    vv1,vv2,vv3,hh12,hh13,jj = S.symbols("V1 V2 V3 H12 H13 H23")
    tt = vv2*vv3
    formal_basis = [vv1**2*vv2*vv3, vv1*hh12*vv3, vv1*hh13*vv2,
                    hh12*hh13, vv1**2*jj]
    exponents = [S.Poly(q,vv1,vv2,vv3,hh12,jj,hh13).monoms()[0]
                 for q in formal_basis]
    assert set(exponents) == set(monomials((2,1,1)))
    targets = [a*tt+jj, -a*tt-jj, -a*tt-jj, b*tt+k*jj, e*jj]
    for q, target in zip(formal_basis, targets):
        value = (e*S.diff(q,vv1,2)/2
                 + ((a-e)*tt+jj)*S.diff(q,vv1,2,vv2,vv3)/2
                 - (a*tt+jj)*(S.diff(q,vv1,hh12,vv3)
                              + S.diff(q,vv1,hh13,vv2))
                 + (b*tt+k*jj)*S.diff(q,hh12,hh13))
        zero(value-target)
    print("OPERATOR: all five polynomial images match; order <= 4 <= 6; no ∂H23")
    # The authors' D2Solution.m has 461 coefficient rules in variables
    # (V1,V2,V3,H12,H31). At Delta2=Delta3=2 precisely these seven
    # nonzero rules can act on this basis; all other derivative monomials
    # kill it or have zero coefficient. This is their critical-weight
    # operator, independently of our off-critical construction above.
    for q, target in zip(formal_basis, targets):
        source_d2 = (-S.diff(q,vv1,2)/2
            + (-2*jj-8*tt)*S.diff(q,hh12,hh13)
            + vv2*S.diff(q,vv1,2,vv2)/2
            + vv3*S.diff(q,vv1,2,vv3)/2
            + (-jj-2*tt)*(S.diff(q,vv1,vv2,hh13)
                          + S.diff(q,vv1,vv3,hh12))
            + (jj+tt)*S.diff(q,vv1,2,vv2,vv3)/2)
        zero(source_d2-target.subs(r,1))
    print("FIDELITY: supplementary critical-weight D2 at Δ=(2,2,2) agrees on all five inputs")
    for n in (2,3):
        coeffs = [(S.expand(t).coeff(tt), S.expand(t).coeff(jj)) for t in targets]
        print(f"r={n}: normalized images in basis (T,J) = "
              f"{[tuple(x.subs(r,n) for x in cs) for cs in coeffs]}")


def dot(q, t):
    return S.sympify(-q[0]*t[0] + sum(q[i]*t[i] for i in range(1,5)))


def ambient():
    """Independent actual-coordinate differentiation, before cone evaluation."""
    p = S.symbols("p0:5")
    z = S.symbols("z0:5")
    eta = (-1,1,1,1,1)
    p1 = (S.Rational(1,2),0,0,0,S.Rational(1,2))
    p2,p3 = (1,1,0,0,0), (1,0,1,0,0)
    z1,z2,z3 = (0,1,0,S.I,0), (1,1,0,S.I,-1), (1,0,1,-S.I,-1)
    for pp,zz in zip((p1,p2,p3),(z1,z2,z3)):
        assert dot(pp,pp) == dot(pp,zz) == S.expand(dot(zz,zz)) == 0
    xx,yy,cc = dot(p,p2),dot(p,p3),dot(p2,p3)
    av1 = (dot(z,p2)*yy-dot(z,p3)*xx)/cc
    av2 = (dot(z2,p3)*xx-dot(z2,p)*cc)/yy
    av3 = (dot(z3,p)*cc-dot(z3,p2)*yy)/xx
    ah12 = -2*(dot(z,z2)*xx-dot(z,p2)*dot(z2,p))
    ah13 = -2*(dot(z,z3)*yy-dot(z,p3)*dot(z3,p))
    aj = -2*(dot(z2,z3)*cc-dot(z2,p3)*dot(z3,p2))
    aq = [av1**2*av2*av3, av1*ah12*av3, av1*ah13*av2,
          ah12*ah13, av1**2*aj]
    point = dict(zip(p,p1)) | dict(zip(z,z1))
    blocks = [S.cancel(q.subs(point)) for q in (av1,av2,av3,ah12,ah13,aj)]
    assert blocks == [S.Rational(1,2),1,-1,-2,1,4]
    print("AMBIENT: real P and complex Z satisfy all nine cone constraints")
    print(f"AMBIENT: blocks (V1,V2,V3,H12,H13,H23) = {blocks}")
    for n in (2,3):
        kk = xx**(-n-1)*cc**(-1)*yy**(-n-1)
        assert not kk.has(S.Float)
        inputs,outputs = [],[]
        for j,q in enumerate(aq):
            F = S.cancel(kk*q)
            # Z-Hessian coefficients; p remains fully ambient throughout.
            hz = [[S.cancel(S.diff(F,z[i],z[l])) for l in range(5)]
                  for i in range(5)]
            div2 = sum(eta[i]*eta[l]*S.diff(hz[i][l],p[i],p[l])
                       for i,l in product(range(5), repeat=2))
            trace = sum(eta[i]*hz[i][i] for i in range(5))
            trace_lap = sum(eta[i]*S.diff(trace,p[i],2) for i in range(5))
            value = S.cancel((S.Rational(3,4)*div2
                              - S.Rational(1,4)*trace_lap).subs(point))
            target = expected[j].subs({r:n,x:dot(p1,p2),y:dot(p1,p3),
                        alpha:dot(p1,z2),beta:dot(p1,z3),c:cc,
                        A:dot(z2,p3),B:dot(z3,p2),L:dot(z2,z3)})
            zero(value - kk.subs(point)*target)
            inputs.append(S.cancel(F.subs(point)));outputs.append(value)
        print(f"AMBIENT r={n}: K={kk.subs(point)}, inputs={inputs}, C1² outputs={outputs}",flush=True)
        assert any(value != 0 for value in outputs)
        if n == 2:
            assert (inputs[0],outputs[0],inputs[4],outputs[4]) == (16,320,-64,-512)


if __name__ == "__main__":
    fidelity()
    family()
    ambient()
    print("ALL_CHECKS_OK")
