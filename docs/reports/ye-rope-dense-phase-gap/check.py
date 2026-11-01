#!/usr/bin/env python3
"""Exact checks for Ye's resonant-ring spectrum and proposed dense-phase route.

No floating-point evaluations or asymptotic proof are used. SymPy rational
functions encode x = exp(beta) symbolically in the fidelity example.
"""
from fractions import Fraction as Q
from math import gcd
import sympy as s


def fidelity():
    # Figure 1(a), n=16, m=4: four distinct phases, each repeated four times.
    # At the figure's beta=2, substitute x=exp(2); identities hold for all x>1.
    x = s.symbols('x', positive=True)
    z = s.symbols('z')
    row = [x, s.Integer(1), 1/x, s.Integer(1)]
    W = s.Matrix(16, 16, lambda i, j: row[(i-j) % 4])
    Z = sum(W[0, j] for j in range(16))
    assert all(sum(W[i, j] for j in range(16)) == Z for i in range(16))
    B = s.Matrix(4, 4, lambda i, j: row[(i-j) % 4]) / sum(row)
    u = (x-1)/(x+1)
    assert s.factor(B.charpoly(z).as_expr() - (z-1)*(z-u)**2*(z-u**2)) == 0
    # Rank and a complete eigenbasis certify all twelve zero eigenvalues.
    assert W.rank() == 4
    roots = [s.Integer(1), s.I, s.Integer(-1), -s.I]
    expected = [1, u, u**2, u]
    for q in range(4):
        v = s.Matrix([roots[q]**j for j in range(16)])
        assert all(s.cancel(a) == 0 for a in W*v/Z - expected[q]*v)
    t = s.symbols('t', positive=True)
    assert s.cancel((u-u**2).subs(x, 1+t) - 2*t/(t+2)**2) == 0
    print('FIDELITY n=16,m=4: rank=4; twelve zeros; spectrum 1,u,u,u^2; u=(x-1)/(x+1)')
    # Theorem 5.3(ii), L=2, with x=exp(beta).
    B2 = s.Matrix([[x, 1/x], [1/x, x]])/(x+1/x)
    v = s.Matrix([1, -1])
    lam = (x*x-1)/(x*x+1)
    assert all(s.cancel(a) == 0 for a in B2*v-lam*v)
    assert s.cancel((1-lam)-2/(x*x+1)) == 0
    print('FIDELITY L=2: gamma=2/(x^2+1), exactly the source formula')


def quotient_checks():
    count = 0
    for n in range(2, 33):
        for m in range(-n, n+1):
            g = gcd(m, n)
            L = n//g
            phase = [(m//g*i) % L for i in range(n)]
            assert all(phase.count(a) == g for a in range(L))
            # A pair is labelled by a phase and one of its g repeated copies.
            order = sorted(range(n), key=lambda i: (phase[i], i))
            for i in range(n):
                for j in range(n):
                    original = (m//g*(order[i]-order[j])) % L
                    tensor = (i//g-j//g) % L
                    assert original == tensor
            count += 1
    print(f'QUOTIENT {count} signed resonances, 2<=n<=32: exact phase multiplicities and B_L tensor J_g/g labels')


def ode_checks():
    # Differentiate S_q: S_q'=(S_{q-1}+S_{q+1})/2. Reduce S_q=S_{L-q}.
    count = 0
    for L in range(4, 65):
        r = L//2
        C = [[Q(0) for _ in range(r+1)] for _ in range(r+1)]
        for i in range(r+1):
            for j in [(i-1) % L, (i+1) % L]:
                C[i][min(j, L-j)] += Q(1, 2)
        M = [[Q(0) for _ in range(r)] for _ in range(r)]
        M[0][0], M[0][1] = -Q(1, 2), Q(1, 2)
        for i in range(1, r-1):
            M[i][i-1] = M[i][i+1] = Q(1, 2)
        M[r-1][r-2] = Q(1, 2)
        if L % 2 == 0:
            M[r-1][r-1] = -Q(1, 2)
        for i in range(r):
            # D*C=M*D for D*S=(S_0-S_1,...,S_{r-1}-S_r).
            left = [C[i][j]-C[i+1][j] for j in range(r+1)]
            right = [Q(0) for _ in range(r+1)]
            for k in range(r):
                right[k] += M[i][k]
                right[k+1] -= M[i][k]
            assert left == right
            assert all(M[i][j] >= 0 for j in range(r) if j != i)
            assert all(M[i][j] + (Q(1, 2) if i == j else 0) >= 0 for j in range(r))
        # Initial data S(0)=(L,0,...,0) give d(0)=(L,0,...,0).
        S0 = [Q(L)] + [Q(0)]*r
        assert [S0[i]-S0[i+1] for i in range(r)] == [Q(L)]+[Q(0)]*(r-1)
        count += 1
    print(f'MODE_ORDER {count} periods, 4<=L<=64: exact cooperative difference identities, both parities, nonnegative shifted generator')


def quadrature_checks():
    count = 0
    for cells in range(1, 33):
        h = Q(1, cells)
        for power in range(1, 17):
            discrete = h*sum((j*h)**power for j in range(cells))
            integral = Q(1, power+1)
            variation = Q(1)
            assert abs(discrete-integral) <= h*variation
            count += 1
    print(f'QUADRATURE {count} exact polynomial cases: left-rule error <= mesh times total variation')


def main():
    fidelity()
    quotient_checks()
    ode_checks()
    quadrature_checks()
    print('SCOPE finite exact identities only; the unbounded joint limit remains for the post-filing probe')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
