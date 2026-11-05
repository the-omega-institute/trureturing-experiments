#!/usr/bin/env python3
"""Exact checks for arXiv:2606.27740v1, Section 2.3(E).
Finite checks support the general argument, not an infinite-dimensional theorem.
Requires Python 3 and SymPy. No floating point is used.
"""
from fractions import Fraction as F
from itertools import product
import sympy as s

I2 = s.eye(2)
Z = s.diag(1, -1)

def clean(z):
    return s.expand(z)

def coin(a, b):
    return s.Matrix([[a, b], [-s.conjugate(b), s.conjugate(a)]])

def equal_matrix(a, b):
    assert all(clean(z) == 0 for z in a - b)

# Symbolic identities modulo the source's unit-sphere equation.
x, y, u, v = s.symbols('x y u v', real=True)
a, b = x + s.I*y, u + s.I*v
C = coin(a, b)
A = (C.conjugate().T * Z * C + Z) / 2
expected = s.Matrix([[x*x+y*y, s.conjugate(a)*b],
                     [a*s.conjugate(b), -x*x-y*y]])
sphere = s.groebner([x*x+y*y+u*u+v*v-1], x, y, u, v, domain=s.QQ)
for residual in [C.conjugate().T*C-I2, A-expected,
                 A.conjugate().T-A, A*A-(x*x+y*y)*I2]:
    for z in residual:
        re, im = clean(z).as_real_imag()
        assert sphere.reduce(clean(re))[1] == 0
        assert sphere.reduce(clean(im))[1] == 0
print('SYMBOLIC: C* C = I; A=(C* sigma_z C+sigma_z)/2; A*=A; A^2=|a|^2 I')

moduli = [(F(0), F(1)), (F(1), F(0)), (F(3,5), F(4,5)),
          (F(4,5), F(3,5)), (F(5,13), F(12,13)), (F(12,13), F(5,13))]
phases = [s.Integer(1), s.Integer(-1), s.I, -s.I]
caps = [F(1), F(3,2), F(2), F(4), F(13)]
local_count = 0
for (am, bm), pa, pb, R in product(moduli, phases, phases, caps):
    a, b = s.Rational(am)*pa, s.Rational(bm)*pb
    C = coin(a, b)
    h = R if am == 0 else min(R, 1/am)
    A = (C.conjugate().T * Z * C + Z)/2
    B = s.Rational(h)*A
    equal_matrix(C.conjugate().T*C, I2)
    equal_matrix(B.conjugate().T, B)
    equal_matrix(B*B, s.Rational(h*am)**2 * I2)
    assert 0 <= h*am <= 1
    local_count += 1
print(f'LOCAL: {local_count} exact coin/cap cases, including zero and unit transmission')

# Full-line finite-support actions, without cyclic boundary artifacts.
def trim(z):
    return {k: clean(w) for k, w in z.items() if clean(w) != 0}

def C_action(z, coefficients, adjoint=False):
    out = {}
    for j in {j for j, _ in z}:
        mat = coin(*coefficients(j))
        if adjoint:
            mat = mat.conjugate().T
        w = mat*s.Matrix([z.get((j,0),0), z.get((j,1),0)])
        out[j,0], out[j,1] = w[0], w[1]
    return trim(out)

def shift(z, inverse=False):
    sign = -1 if inverse else 1
    return {(j + sign*(1 if k == 0 else -1), k): w for (j,k),w in z.items()}

def walk(z, coefficients):
    return shift(C_action(z, coefficients))

def adjoint_walk(z, coefficients):
    return C_action(shift(z, inverse=True), coefficients, adjoint=True)

def minus(z, w):
    return trim({k: z.get(k,0)-w.get(k,0) for k in z.keys() | w.keys()})

# Source fidelity: Section 2.3(D), delta=3/5, b=4/5 on the right half-line.
def halfline(j):
    return (s.Integer(1),s.Integer(0)) if j <= 0 else (s.Rational(3,5),s.Rational(4,5))
z = {(0,1):s.Integer(1)}
for t in range(21):
    assert z == {(-t,1):s.Integer(1)}
    qnorm_sq = sum(j*j*clean(s.conjugate(w)*w) for (j,_),w in z.items())
    assert qnorm_sq == t*t
    z = walk(z, halfline)
for N in range(1,13):
    H = F(2*N)/(N + N/F(3,5))
    assert H == F(3,4) < 1
print('FIDELITY: source remark (D), U^t delta_0^- = delta_-t^- for 0<=t<=20; Q-norm^2=t^2')
print('FIDELITY: 2N/(N+N/delta)=3/4<1 for 1<=N<=12, delta=3/5')

profile_count = increment_count = 0
alphabet = [moduli[0],moduli[2],moduli[1]]
for p in range(1,5):
    for profile in product(alphabet, repeat=p):
        def coefficients(j):
            am,bm = profile[(j-1)%p]
            return s.Rational(am)*phases[j%4], s.Rational(bm)*phases[(2*j+1)%4]
        for R in [F(1), F(2), F(4)]:
            def h(j):
                am = profile[(j-1)%p][0]
                return R if am == 0 else min(R,1/am)
            def r(j):
                if j >= 0:
                    return sum((h(k) for k in range(1,j+1)),F(0))
                return -sum((h(k) for k in range(j+1,1)),F(0))
            def X(z):
                return trim({(j,k):s.Rational(r(j-1 if k == 0 else j))*w for (j,k),w in z.items()})
            M = sum((h(j) for j in range(1,p+1)),F(0))/p
            for j in range(-3*p,3*p+1):
                assert r(j)-r(j-1) == h(j)
                assert r(j+p)-r(j) == p*M
                assert r(j+p)-M*(j+p) == r(j)-M*j
            for j,k in product(range(-2,3),range(2)):
                z = {(j,k):s.Integer(1)}
                lhs = minus(adjoint_walk(X(walk(z,coefficients)),coefficients),X(z))
                C = coin(*coefficients(j))
                A = (C.conjugate().T*Z*C+Z)/2
                w = s.Rational(h(j))*A*s.eye(2)[:,k]
                rhs = trim({(j,0):w[0],(j,1):w[1]})
                assert lhs == rhs
                increment_count += 1
            profile_count += 1
print(f'GLOBAL: {increment_count} basis-vector checks of U* X U-X=h_j A_j on the full integer line')
print(f'MEANS: {profile_count} profile/cap cases; all words of lengths 1..4 over transmissions {{0,3/5,1}}')
print('MEANS: exact positive/negative period increments and periodic remainder r_j-M*j')
print('SCOPE: general limit-periodic conclusion uses the ordinary proof outline; no Lean proof is claimed')
print('ALL_CHECKS_OK')
