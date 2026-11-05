#!/usr/bin/env python3
"""Exact fidelity and obstruction checks for arXiv:2610.39674v1.

No floating point is used. The finite checks support the proposed argument;
they do not certify Haar/Gaussian probabilities or replace the eventual Lean proof.
"""
from itertools import product
from fractions import Fraction as F
import sympy as s


def gaussian_expectation(poly, variables):
    """Exact standard-independent-Gaussian polynomial moments."""
    total = 0
    for powers, coeff in s.Poly(s.expand(poly), *variables).terms():
        moment = 1
        for power in powers:
            if power % 2:
                moment = 0
                break
            moment *= s.factorial2(power - 1) if power else 1
        total += coeff * moment
    return s.expand(total)


def fidelity():
    # Section 7.2: population MSE = squared coefficient residual, at all tested depths.
    cases = 0
    for D, M, L in product(range(1, 4), range(1, 5), range(1, 4)):
        dims = [D] + [M] * L
        x = s.Matrix(s.symbols(f'x0:{D}'))
        u = s.Matrix([j + 1 for j in range(D)])
        a = s.Matrix([(-1) ** j * (j + 1) for j in range(M)])
        P = s.eye(D)
        for l in range(L):
            W = s.Matrix(dims[l+1], dims[l], lambda i,j: (3*i+2*j+l) % 5 - 2)
            P = W * P
        residual = P.T * a - u
        scalar = (residual.T * x)[0]
        assert gaussian_expectation(scalar**2, list(x)) == (residual.T * residual)[0]
        cases += 1
    # Symbolic independent check of the same identity in dimension 3.
    x = s.symbols('x0:3'); r = s.symbols('r0:3')
    assert gaussian_expectation(sum(r[i]*x[i] for i in range(3))**2, x) == sum(v*v for v in r)
    print(f'FIDELITY population-MSE identity: {cases} integer networks, D=1..3, M=1..4, L=1..3; symbolic D=3')


def obstruction():
    cases = 0
    for D, M in product(range(1, 5), range(1, 7)):
        p = M*(D+1)
        B = s.zeros(p)
        for i in range(M):
            B[i*D, M*D+i] = B[M*D+i, i*D] = s.Rational(1, 2)
        # q = a^T W e_1; exact eigenvalues determine both signs and the radical.
        eig = B.eigenvals()
        want = {s.Rational(-1,2): M, s.Rational(1,2): M}
        if D > 1:
            want[s.Integer(0)] = M*(D-1)
        assert eig == want
        for d in range(M+1):
            T = s.zeros(p, d)
            for i in range(d):
                T[i*D, i] = -1
                T[M*D+i, i] = 1
            assert T.T*B*T == -s.eye(d)
            assert T.T*T == 2*s.eye(d)
            # Each translated plane with a small rational offset misses q=c,
            # certified by its exact quadratic maximum; includes all d=0..M.
            for denom in (4*p, 8*p):
                g = s.Matrix([s.Rational((-1)**i, denom) for i in range(p)])
                R = T.T*B*T
                b = T.T*B*g
                maximum = (g.T*B*g)[0] - ((b.T*R.inv()*b)[0] if d else 0)
                assert (g.T*g)[0] < s.Rational(1,3)
                assert maximum < s.Rational(1,2)
                cases += 1
        # Explicit nonempty solution and rank-D Jacobian at that solution.
        W = s.zeros(M, D); W[0,0] = 1
        a = s.zeros(M, 1); a[0] = 1
        assert W.T*a == s.eye(D)[:,0]
        J = s.zeros(D,p)
        for j in range(D):
            for i in range(M):
                J[j,i*D+j] = a[i]
                J[j,M*D+i] = W[i,j]
        assert J.rank() == D
    # Universal scalar estimates reduce to polynomial identities and positivity.
    r,t = s.symbols('r t', nonnegative=True)
    assert s.expand(s.Rational(3,2)*r*r - (r*r/2+r*t-t*t/4) - (t-2*r)**2/4) == 0
    # A nonempty open raw-matrix neighborhood is quantified without irrational bases:
    # T has norm sqrt(2), q(Tx)=-||x||^2. If ||E||_op<1/8, then
    # q((T+E)x) <= (-1+sqrt(2)/8+1/128)||x||^2 < -1/2 ||x||^2.
    assert F(2,8) + F(1,128) < F(1,2)  # sqrt(2) < 2 gives a rational upper bound.
    # Strong Grassmann neighborhood: projection distance <1/8 gives
    # q(t) <= [-1/2*(1-delta^2)+delta+delta^2/2] ||t||^2 < -1/4 ||t||^2.
    delta=F(1,8)
    assert -F(1,2)*(1-delta*delta)+delta+delta*delta/F(2) < -F(1,4)
    print(f'OBSTRUCTION exact signatures: 24 (D,M) pairs; {cases} translated-plane maxima with d=0..M')
    print('OBSTRUCTION universal square-completion identity and strict rational neighborhood margins: OK')


def scalar_threshold():
    # D=1: any (M+1)-plane in dimension 2M meets each signed M-plane.
    # Dimension inequality gives a nonzero vector in each intersection.
    for M in range(1, 30):
        assert (M+1)+M-2*M == 1
    w,a,b,c,v,z = s.symbols('w a b c v z')
    # Exact line criterion for M=1 and the family centered at the negative line.
    expr = s.expand((w+z*v)*(a+z*b)-c)
    assert s.Poly(expr,z).all_coeffs() == [b*v, a*v+b*w, a*w-c]
    for W,A in product(range(-2,3),repeat=2):
        for C in (1,2,3):
            discr=(A-W)**2+4*(W*A-C)
            assert discr == (A+W)**2-4*C
    print('SCALAR threshold dimension bounds: M=1..29; line polynomial/discriminant identity: OK')


if __name__ == '__main__':
    fidelity()
    obstruction()
    scalar_threshold()
    print('ALL_CHECKS_OK')
