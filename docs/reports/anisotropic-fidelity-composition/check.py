#!/usr/bin/env python3
"""Exact sphere-moment and SU(2)-composition checks; no Lean proof is claimed."""

from fractions import Fraction as F
from functools import lru_cache
from math import factorial, prod
import sympy as S

q = S.symbols("q0:4", real=True)
r = S.symbols("r0:4", real=True)
k = S.symbols("k", integer=True, positive=True)
a = S.symbols("a", real=True)
R = S.Rational


def rising(x, count):
    return prod((x + j for j in range(count)), start=F(1))


@lru_cache(None)
def haar_monomial(exponents):
    """Normalized surface measure on S^3: Dirichlet(1/2,1/2,1/2,1/2)."""
    if any(e % 2 for e in exponents):
        return F(0)
    halves = [e // 2 for e in exponents]
    return prod((rising(F(1, 2), h) for h in halves), start=F(1)) / rising(
        F(2), sum(halves)
    )


def haar_integral(poly):
    result = S.Integer(0)
    for exponents, coefficient in S.Poly(S.expand(poly), *q).terms():
        moment = haar_monomial(exponents)
        result += coefficient * R(moment.numerator, moment.denominator)
    return S.factor(result)


def exact_equal(left, right):
    assert S.simplify(left - right) == 0, (left, right)


def U(x):
    return S.Matrix([[x[0] + S.I*x[1], -x[2] + S.I*x[3]],
                     [x[2] + S.I*x[3], x[0] - S.I*x[1]]])


def physical_fidelity(M, N):
    """Expand the actual product matrix before applying independent moments."""
    z = (U(r)*U(q))[0, 0]
    square = S.Poly(S.expand(z*S.conjugate(z)), *(r+q))
    value = S.Integer(0)
    for exponents, coefficient in square.terms():
        er, eq = exponents[:4], exponents[4:]
        ir = [i for i, e in enumerate(er) for _ in range(e)]
        iq = [i for i, e in enumerate(eq) for _ in range(e)]
        assert len(ir) == len(iq) == 2
        value += coefficient*N[ir[0], ir[1]]*M[iq[0], iq[1]]
    return S.factor(value)


def beta_x_moment(power, alpha, beta):
    """E[(2u-1)^power] for u ~ Beta(alpha,beta), by finite exact expansion."""
    return S.factor(sum(S.binomial(power, j)*2**j*(-1)**(power-j)
                        *S.rf(alpha, j)/S.rf(alpha+beta, j)
                        for j in range(power+1)))


def main():
    exact_equal(U(q).det(), sum(x*x for x in q))
    assert S.simplify(U(q).conjugate().T*U(q)
                      - sum(x*x for x in q)*S.eye(2)) == S.zeros(2)
    print("SU2: determinant and adjoint-product identities checked symbolically")

    # Reproduce the journal's Example 1 without using its answer as an input.
    g1 = (q[1]**2 + 2*q[3]**2)/3
    g2 = 1-q[0]**4
    Z1, Z2 = haar_integral(g1), haar_integral(g2)
    p1 = haar_integral((q[0]**2+q[1]**2)*g1)/Z1
    p2 = haar_integral((q[0]**2+q[1]**2)*g2)/Z2
    M1 = S.Matrix(4, 4, lambda i,j: haar_integral(q[i]*q[j]*g1)/Z1)
    M2 = S.Matrix(4, 4, lambda i,j: haar_integral(q[i]*q[j]*g2)/Z2)
    actual = physical_fidelity(M1, M2)
    exact_equal(Z1, R(1,4)); exact_equal(Z2, R(7,8))
    exact_equal(p1, R(4,9)); exact_equal(p2, R(13,28))
    exact_equal(actual, R(127,252))
    exact_equal(actual, 2*(p1-R(1,2))*(p2-R(1,2))+R(1,2))
    print(f"SOURCE_EXAMPLE_1: normalized masses {Z1}, {Z2}; F1^2={p1}; F2^2={p2}; Fsum^2={actual}")
    print("SOURCE_EXAMPLE_1: surface volumes are 2*pi^2 times the normalized masses")

    D = (k+3)*(k+4)
    A = (k*k+k+3)/D
    B = (2*k+3)/D
    kap = k/(k+5)
    var = 5*(2*k+5)/((k+5)**2*(k+6))
    m = a*B*var
    p = A+B
    base_mean = beta_x_moment(1, k+R(3,2), R(3,2))
    base_second = beta_x_moment(2, k+R(3,2), R(3,2))
    tilted_mean = beta_x_moment(1, k+R(5,2), R(5,2))
    tilted_second = beta_x_moment(2, k+R(5,2), R(5,2))
    exact_equal(base_mean, k/(k+3)); exact_equal(base_second, A)
    exact_equal((1-A)/3, B); exact_equal(tilted_mean, kap)
    exact_equal(tilted_second-tilted_mean**2, var)
    exact_equal(a*B*(tilted_mean-kap), 0)
    exact_equal(a*B*(tilted_second-kap*tilted_mean), m)
    exact_equal(p, (k*k+3*k+6)/D)
    M = S.diag(A,B,B,B); M[0,2] = M[2,0] = m
    composed = physical_fidelity(M,M)
    predicted = p*p+(1-p)**2
    exact_equal(composed, predicted-2*m*m)
    assert S.factor(composed-predicted) != 0
    print("SYMBOLIC_FAMILY: Beta moments, zero transverse mean, mixed second moment and composition checked")

    # Positivity is an analytic bound: |a q2(q0-kappa)| < 2a < 1/10.
    # All displayed differences have positive denominators for k>=2.
    h = p-R(1,2)
    exact_equal(h-R(1,30), (k-2)*(7*k+3)/(15*D))
    exact_equal(R(1,3)-B, A/3)
    exact_equal(1-var, (k**3+16*k**2+75*k+125)/((k+5)**2*(k+6)))
    exact_equal(composed-R(1,2), 2*(h*h-m*m))
    exact_equal(S.limit(p,k,S.oo), 1)
    exact_equal(S.limit(m,k,S.oo), 0)
    exact_equal(S.limit(k**3*m/a,k,S.oo),20)
    exact_equal(S.limit(k**6*(composed-predicted)/a**2,k,S.oo),-800)
    print("BOUNDS: h>=1/30, 0<m<a/3<1/60 imply composed F^2>1/2")
    print("LIMITS: F^2->1; m->0; k^3*m/a->20; k^6*defect/a^2->-800")

    # Independent direct polynomial integration checks the claimed Beta reduction.
    count = 0
    for integer_k in range(2,15):
        z = R(2*prod(range(1,2*integer_k+2,2)), factorial(integer_k+2))
        weight = (1+q[0])**integer_k
        exact_equal(haar_integral(weight),z)
        for rational_a in [R(1,100),R(1,40),R(1,21)]:
            density = weight*(1+rational_a*q[2]*(q[0]-R(integer_k,integer_k+5)))/z
            exact_equal(haar_integral(density),1)
            first = [haar_integral(qi*density) for qi in q]
            assert first == [R(integer_k,integer_k+3),0,0,0]
            direct_M = S.Matrix(4,4,lambda i,j:haar_integral(q[i]*q[j]*density))
            subs = {k:integer_k,a:rational_a}
            assert S.simplify(direct_M-M.subs(subs)) == S.zeros(4)
            actual = physical_fidelity(direct_M,direct_M)
            expected = S.factor(composed.subs(subs))
            exact_equal(actual,expected)
            assert R(1,2) < actual < 1
            assert R(1,2) < p.subs(k,integer_k) < 1
            assert actual < predicted.subs(k,integer_k)
            count += 1
        print(f"DIRECT_HAAR k={integer_k}: Z={z}; p={p.subs(k,integer_k)}; all 3 rational perturbations checked")
    exact_equal(p.subs(k,2),R(8,15))
    exact_equal(m.subs({k:2,a:R(1,40)}),R(3,4480))
    exact_equal(composed.subs({k:2,a:R(1,40)}),R(113,225)-R(9,10035200))
    print(f"DIRECT_HAAR: {count} independent exact polynomial density checks")
    print("EXACT_MEMBER k=2,a=1/40: p=8/15; m=3/4480; Fsum^2=113/225-9/10035200")
    print("SCOPE: verifies the stated SU(2) model; identification with the source's general sum remains a separate obligation")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
