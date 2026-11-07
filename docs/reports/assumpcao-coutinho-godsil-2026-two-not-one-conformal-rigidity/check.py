#!/usr/bin/env python3
"""Exact checks for the Section 9 two-not-one conformal-rigidity question.

Source: Assumpcao, Coutinho, Godsil, arXiv:2605.08508v2.
Only integer/rational/symbolic arithmetic is used. The script checks certificates;
the universal spectral variational argument is stated in README.md.
"""
from itertools import combinations
import sympy as s

Q = s.Rational
n = 5
edges = [e for e in combinations(range(n), 2) if e != (0, 1)]
assert len(edges) == 9
I = s.eye(n)
J = s.ones(n)
one = s.ones(n, 1)

def z(e):
    v = s.zeros(n, 1)
    v[e[0]] = 1
    v[e[1]] = -1
    return v

vectors = [z(e) for e in edges]

def laplacian(weights):
    return sum((w * v * v.T for w, v in zip(weights, vectors)), s.zeros(n))

def sorted_spectrum(A):
    values = A.eigenvals()
    assert all(x.is_Rational for x in values)
    return sorted(x for x, mult in values.items() for _ in range(mult))

# Fidelity: source equations (1) and (2), symbolically for arbitrary weights
# and every symmetric test matrix, not just the proposed certificate.
w = s.symbols('w0:9', nonnegative=True)
t = s.symbols('t0:15', real=True)
T = s.zeros(n)
for value, (a, b) in zip(t, combinations(range(n + 1), 2)):
    # Pairs in range(n+1) parametrise the upper triangle via (a,b-1).
    T[a, b - 1] = value
    T[b - 1, a] = value
L = laplacian(w)
for (a, b), v in zip(edges, vectors):
    assert s.expand((v.T * T * v)[0] - (T[a, a] + T[b, b] - 2*T[a, b])) == 0
assert s.expand(s.trace(L*T) - sum(weight*(v.T*T*v)[0] for weight, v in zip(w, vectors))) == 0
assert L*one == s.zeros(n, 1)
assert s.expand(s.trace(L) - 2*sum(w)) == 0
print('FIDELITY_OK: source equations (1),(2), symbolic adjoint identity, row sums, trace')

U = laplacian([s.Integer(1)]*len(edges))
assert sorted_spectrum(U) == [0, 3, 5, 5, 5]
print('uniform spectrum:', sorted_spectrum(U))

X = Q(1,45)*s.Matrix([
    [12,12,-8,-8,-8],
    [12,12,-8,-8,-8],
    [-8,-8,22,-3,-3],
    [-8,-8,-3,22,-3],
    [-8,-8,-3,-3,22],
])
assert X == X.T
assert sorted_spectrum(X) == [0, 0, Q(5,9), Q(5,9), Q(8,9)]
assert all(v >= 0 for v in sorted_spectrum(X))
assert all(v >= 0 for v in sorted_spectrum(I-X))
assert s.trace(X) == 2
assert X*one == s.zeros(n,1)
assert {(v.T*X*v)[0] for v in vectors} == {Q(10,9)}
assert s.expand(s.trace(L*X) - Q(10,9)*sum(w)) == 0
assert s.trace(U*X) == 10
print('X spectrum:', sorted_spectrum(X))
print('I-X spectrum:', sorted_spectrum(I-X))
print('X trace:', s.trace(X), '; every edge quadratic value:', Q(10,9))
print('SYMBOLIC_UPPER_CERTIFICATE_OK: tr(L(w) X) = (10/9) sum(w) = 10 on Delta_E')

P = I - Q(1,n)*J
Y = P-X
assert P*P == P and P*one == s.zeros(n,1)
assert Y == Y.T and Y*one == s.zeros(n,1)
assert s.trace(Y) == 2
assert all(v >= 0 for v in sorted_spectrum(Y))
assert all(v >= 0 for v in sorted_spectrum(I-Y))
assert {(v.T*Y*v)[0] for v in vectors} == {Q(8,9)}
assert s.expand(s.trace(L*Y) - Q(8,9)*sum(w)) == 0
assert s.trace(U*Y) == 8
print('Y spectrum:', sorted_spectrum(Y))
print('SYMBOLIC_LOWER_CERTIFICATE_OK: tr(L(w) Y) = (8/9) sum(w) = 8 on Delta_E')

perturbed_weights = [Q(6,5) if a < 2 else Q(3,5) for a, b in edges]
assert sum(perturbed_weights) == 9
assert min(perturbed_weights) > 0
W = laplacian(perturbed_weights)
assert sorted_spectrum(W) == [0, Q(18,5), Q(21,5), Q(21,5), 6]
assert sorted_spectrum(W)[1] > sorted_spectrum(U)[1]
print('perturbed weights:', perturbed_weights)
print('perturbed spectrum:', sorted_spectrum(W))
print('NOT_ONE_RIGID_OK: lambda_2 = 18/5 > 3')

# Exact boundary checks include disconnected admissible weighted graphs.
for j in range(len(edges)):
    weights = [s.Integer(9) if i == j else s.Integer(0) for i in range(len(edges))]
    A = laplacian(weights)
    ev = sorted_spectrum(A)
    assert ev == [0,0,0,0,18]
    assert sum(ev[-2:]) >= 10
    assert sum(ev[1:3]) <= 8
    assert s.trace(A*X) == 10 and s.trace(A*Y) == 8
print('BOUNDARY_OK: all nine simplex vertices; additional zero eigenvalues are allowed')
print('ALL_CHECKS_OK')
