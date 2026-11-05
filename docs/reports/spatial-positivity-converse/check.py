#!/usr/bin/env python3
"""Exact algebra checks for the spatial-positivity coefficient converse.

No numerical evaluation of analytic characters or PSD eigenvalues is used.
The infinite-series convergence and modular transformations require proof.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

# Q(sqrt(2)) with exact ordering; pairs encode a + b sqrt(2).
def add(x, y):
    return x[0] + y[0], x[1] + y[1]
def mul(x, y):
    return x[0]*y[0] + 2*x[1]*y[1], x[0]*y[1] + x[1]*y[0]
def nonnegative(x):
    a, b = x
    if a >= 0 and b >= 0:
        return True
    if a < 0 and b <= 0:
        return False
    return a*a >= 2*b*b if a >= 0 else 2*b*b >= a*a
Z = (F(0), F(0))
O = (F(1), F(0))
R = (F(0), F(1))
def total(xs):
    z = Z
    for x in xs:
        z = add(z, x)
    return z

# Fidelity: source section 3.2.1, labels eq:TTDLIsing and
# eq:projectors_verlindes_ising, projector
# decompositions and the three nonnegative combinations of simple TDLs.
T = [[O, O, O], [O, O, (-F(1), F(0))], [R, (F(0), -F(1)), Z]]
B = [[(F(1,4),0), (F(1,4),0), (0,F(1,4))],
     [(F(1,4),0), (F(1,4),0), (0,-F(1,4))],
     [(F(1,2),0), (-F(1,2),0), Z]]
assert [[total(mul(B[i][k], T[k][j]) for k in range(3))
         for j in range(3)] for i in range(3)] == [
             [O if i == j else Z for j in range(3)] for i in range(3)]
positive = []
for size in range(1,4):
    for support in combinations(range(3), size):
        coefficients = [total(B[i][j] for i in support) for j in range(3)]
        if all(map(nonnegative, coefficients)):
            positive.append(support)
assert set(positive) == {(0,), (0,1), (0,1,2)}
print('FIDELITY: Ising projector inverse and precisely the 3 source-positive supports:', positive)

# The normalized Walsh matrix is S[a,b] = (-1)^(a dot b)/sqrt(N).
# Use unnormalized integer rows and divide products by N, never approximate sqrt(N).
def sign(a, b):
    return -1 if (a & b).bit_count() % 2 else 1
A = {0, 4, 2, 5}
H = {0, 4, 2, 6}
Q = [[2,1,1,0,0,1,-1,0], [1,2,0,1,1,0,0,-1],
     [1,0,2,1,-1,0,0,1], [0,1,1,2,0,-1,1,0],
     [0,1,-1,0,2,1,1,0], [1,0,0,-1,1,2,0,1],
     [-1,0,0,1,1,0,2,1], [0,-1,1,0,0,1,1,2]]
C = [[F(sum(sign(a,x)*sign(b,x) for x in A),8)
      for b in range(8)] for a in range(8)]
assert C == [[F(z,4) for z in row] for row in Q]
assert C[0][6] == -F(1,4)
assert all(sum(sign(a,x)*sign(b,x) for x in range(8)) ==
           (8 if a == b else 0) for a in range(8) for b in range(8))
assert 5 != 6 and (5).bit_count() == (6).bit_count()
assert A ^ H == {5,6}
print('BASE: full 8-by-8 transformed coefficient matrix equals Q/4; C[000,110] =', C[0][6])
for row in C:
    print(' '.join(map(str, row)))

# Symbolic full-character check: x=|f0|^2, y=|f1|^2.
# The exact weight enumerator is x*(x+y)^(r+2).
# Fourier transforms are checked by independent direct summation, including
# every coefficient-matrix entry, for r=0,...,6.
for r in range(7):
    n = 3+r
    N = 1 << n
    support = {a*(1 << r)+y for a in A for y in range(1 << r)}
    assert 0 in support and len(support) == 4*(1 << r)
    enumerator = [sum(x.bit_count() == w for x in support) for w in range(n+1)]
    assert enumerator == [comb(r+2,w) if w <= r+2 else 0 for w in range(n+1)]
    ft = [F(sum(sign(k,x) for x in support),N) for k in range(N)]
    for a in range(N):
        for b in range(N):
            expected = C[a >> r][b >> r] if (a % (1 << r)) == (b % (1 << r)) else F(0)
            assert ft[a ^ b] == expected
    assert ft[6 << r] == -F(1,4)
    print(f'FAMILY r={r}: sectors={N}, selected={len(support)}, unique vacuum, '
          'weight polynomial x*(x+y)^(r+2), transform C tensor I, negative entry -1/4')

# All-r algebraic recurrence underlying the finite family checks:
# the appended-bit Fourier sum is 1+(-1)^d = 2*delta(d,0).
assert [sum(sign(d,y) for y in range(2)) for d in range(2)] == [2,0]
# S2 preserves E and sends |f0|^2 to |f0+f1|^2/2:
# exact coefficient matrices in independent holomorphic/antiholomorphic symbols.
assert [[F(sum(sign(i,k)*sign(j,k) for k in range(2)),2)
         for j in range(2)] for i in range(2)] == [[F(1),F(0)],[F(0),F(1)]]
assert [[F(sign(i,0)*sign(j,0),2) for j in range(2)] for i in range(2)] == [[F(1,2)]*2]*2

# Formal q-series evidence, not an analytic-convergence claim.  After factoring
# q^(-1/24) from f0 and q^(5/24) from f1, coefficients through degree 24
# are exact convolutions of partition numbers and the shifted A1 theta series.
L = 24
part = [1]+[0]*L
for d in range(1,L+1):
    for m in range(d,L+1):
        part[m] += part[m-d]
theta0 = [sum(z*z == m for z in range(-L-1,L+2)) for m in range(L+1)]
theta1 = [sum(z*(z+1) == m for z in range(-L-1,L+2)) for m in range(L+1)]
chars = [[sum(t[k]*part[m-k] for k in range(m+1)) for m in range(L+1)]
         for t in [theta0,theta1]]
assert chars[0][:6] == [1,3,4,7,13,19]
assert chars[1][:6] == [2,2,6,8,14,20]
assert all(c >= 0 for row in chars for c in row)
print('Q-SERIES: exact coefficients through degree 24, leading coefficients 1 and 2')
print('ANALYTIC BOUNDARY: convergence, theta/eta modular laws, and the full upper-half-plane Gram limit require proof.')
print('ALL_CHECKS_OK')
