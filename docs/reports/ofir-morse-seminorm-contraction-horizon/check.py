#!/usr/bin/env python3
"""Exact certificates for a fixed seminorm and a family of stochastic matrices.

Only integers and Fraction are used. Finite checks are arithmetic evidence,
not a proof of the universal horizon statement. See README.md for the
operator-norm interpretation of the bilinear certificates.
"""
from fractions import Fraction as F
from itertools import product

N = 3

def mat(rows):
    return tuple(tuple(F(x) for x in row) for row in rows)

def add(A, B):
    return mat([[A[i][j] + B[i][j] for j in range(N)] for i in range(N)])

def scale(c, A):
    return mat([[c * A[i][j] for j in range(N)] for i in range(N)])

def sub(A, B):
    return add(A, scale(-1, B))

def mul(A, B):
    return mat([[sum(A[i][r] * B[r][j] for r in range(N))
                 for j in range(N)] for i in range(N)])

def mv(A, x):
    return tuple(sum(a*b for a, b in zip(row, x)) for row in A)

def vm(x, A):
    return tuple(sum(x[i]*A[i][j] for i in range(N)) for j in range(N))

def dot(x, y):
    return sum(a*b for a, b in zip(x, y))

def transpose(A):
    return mat(zip(*A))

def frobenius_squared(A):
    return sum(x*x for row in A for x in row)

def coe(A):
    return max(sum(abs(A[i][r]-A[j][r]) for r in range(N))
               for i in range(N) for j in range(N))/2

def quotient_linf_corner_value(A):
    # A preserves constant vectors. Centering a vector of quotient norm one
    # places it in [-1,1]^3; the convex piecewise-linear objective attains its
    # maximum at a cube corner. Nonconstant corners have quotient norm one.
    values = []
    for x in product((F(-1), F(1)), repeat=N):
        y = mv(A, x)
        values.append((max(y)-min(y))/2)
    return max(values)

def stochastic(A):
    return all(x >= 0 for row in A for x in row) and all(sum(row) == 1 for row in A)

I = mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
J = scale(F(1, 3), mat([[1]*N for _ in range(N)]))
H = sub(I, J)
P = mat([[1, 0, 0], [0, 1, 0], [1, 0, 0]])
T = mul(H, P)
S = add(scale(F(29, 30), P), scale(F(1, 30), J))
x = (F(2), F(-3), F(1))
pi = (F(59, 90), F(30, 90), F(1, 90))
L = mat([pi]*N)
assert mul(P, P) == P
assert mul(H, H) == H and transpose(H) == H
assert mul(mul(H, P), H) == T
assert mul(T, T) == T and mul(H, T) == T and mul(T, H) == T
assert vm(pi, S) == pi and sum(pi) == 1 and min(pi) > 0
assert mul(sub(I, L), H) == sub(I, L)
assert mv(H, x) == x and dot(x, x) == 14 and dot(x, mv(T, x)) == 15
print('STRUCTURE: P^2=P; H^2=H=H^transpose; HPH=T; T^2=HT=TH=T')
print('STRUCTURE: pi=(59,30,1)/90; pi*S=pi; (I-L)H=I-L; x*x=14; x*T*x=15')

# Source fidelity: Section II, Proposition 1, |A|_infinity = coe(A),
# with strict threshold equivalent to scrambling. Include a non-scrambling
# boundary matrix and positive matrices. This is an independent corner method.
for name, A in [('P', P), ('S', S), ('A_1', add(scale(F(1,2), I), scale(F(1,2), S)))]:
    assert stochastic(A)
    coefficient = coe(A)
    assert quotient_linf_corner_value(A) == coefficient
    scrambling = all(any(A[i][r] > 0 and A[j][r] > 0 for r in range(N))
                     for i in range(N) for j in range(N))
    assert (coefficient < 1) == scrambling
    print(f'FIDELITY Proposition 1: {name}: coe=induced_linf={coefficient}; scrambling={scrambling}')

# Exact finite grid for the two scalar inequalities used in the outline.
for denominator in (1, 2, 3, 7, 30, 40):
    for numerator in range(denominator+1):
        u = F(numerator, denominator)
        for j in range(101):
            assert 1-j*u <= (1-u)**j <= 1-j*u+F(j*(j-1),2)*u*u
print('SCALAR: Bernoulli and quadratic upper bounds: exact grid, j=0..100')

ks = list(range(1, 13)) + [20, 50, 100]
power_checks = 0
for K in ks:
    e = F(1, 2*K)
    a, b = 1-e, 1-e/30
    A = add(scale(a, I), scale(e, S))
    assert stochastic(A) and min(z for row in A for z in row) >= e/90 > 0
    assert 0 < a < b < 1 and vm(pi, A) == pi
    assert mul(H, A) == add(scale(a, H), scale(F(29,30)*e, T))
    Aj = I
    for j in range(1, K+1):
        Aj = mul(Aj, A)
        HA = mul(H, Aj)
        assert HA == add(scale(a**j, H), scale(b**j-a**j, T))
        q = dot(x, mv(HA, x))
        lower = 14 + F(j, 2)*e*(1-(j-1)*e)
        assert q == 15*b**j-a**j and q >= lower > 14
        assert sub(Aj, L) == mul(sub(I, L), HA)
        power_checks += 1
    horizon_lower = 1 + F(K+1, 112*K)
    assert q/14 >= horizon_lower > F(113,112)
    print(f'FAMILY: K={K}: all j=1..{K} exact; positive; stationary; '
          f'nu(A^K)>={horizon_lower}>113/112')

# Source fidelity: an explicit instance of Corollary 2 (eventual contraction).
# nu(A^512)^2 <= ||H A^512||_F^2 < 1. No singular-value floats.
K = 1
e = F(1, 2*K)
a, b = 1-e, 1-e/30
HA512 = add(scale(a**512, H), scale(b**512-a**512, T))
assert frobenius_squared(HA512) < 1
print('FIDELITY Corollary 2: K=1, ||H A_1^512||_F^2<1, hence nu(A_1^512)<1')
print(f'CHECKED: {len(ks)} horizons; {power_checks} positive-length matrix powers; exact arithmetic only')
print('ALL_CHECKS_OK')
