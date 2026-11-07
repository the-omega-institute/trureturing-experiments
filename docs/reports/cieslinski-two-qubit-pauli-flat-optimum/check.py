#!/usr/bin/env python3
"""Exact Gaussian-integer certificates for the two-qubit Pauli-flat optimum.

Only Python's standard library is required. A Gaussian integer is an (int, int)
pair; rational matrices are represented by an integer numerator and denominator.
No numerical eigenvalue approximation or floating-point arithmetic is used.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from collections import Counter

ZERO = (0, 0)
ONE = (1, 0)

def add(a, b):
    return (a[0] + b[0], a[1] + b[1])

def neg(a):
    return (-a[0], -a[1])

def times(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def conj(a):
    return (a[0], -a[1])

def gsum(values):
    result = ZERO
    for a in values:
        result = add(result, a)
    return result

def matrix(rows):
    return tuple(tuple((x, 0) if isinstance(x, int) else x for x in row) for row in rows)

def madd(A, B):
    return tuple(tuple(add(a, b) for a, b in zip(ar, br)) for ar, br in zip(A, B))

def scale(k, A):
    return tuple(tuple((k * a[0], k * a[1]) for a in row) for row in A)

def mmul(A, B):
    return tuple(tuple(gsum(times(A[i][k], B[k][j]) for k in range(len(B)))
                       for j in range(len(B[0]))) for i in range(len(A)))

def dagger(A):
    return tuple(tuple(conj(A[j][i]) for j in range(len(A))) for i in range(len(A[0])))

def trace(A):
    return gsum(A[i][i] for i in range(len(A)))

def trace_product(A, B):
    return gsum(times(A[i][j], B[j][i]) for i in range(len(A)) for j in range(len(A)))

def kron(A, B):
    return tuple(tuple(times(A[i][j], B[k][l]) for j in range(2) for l in range(2))
                 for i in range(2) for k in range(2))

def outer(v):
    return tuple(tuple(times(a, conj(b)) for b in v) for a in v)

def determinant(A):
    result = ZERO
    for p in permutations(range(len(A))):
        parity = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
        term = ONE
        for i, j in enumerate(p):
            term = times(term, A[i][j])
        result = add(result, neg(term) if parity % 2 else term)
    return result

def principal_minors(A):
    result = {}
    for k in range(1, len(A) + 1):
        result[k] = []
        for indices in combinations(range(len(A)), k):
            minor = tuple(tuple(A[i][j] for j in indices) for i in indices)
            value = determinant(minor)
            assert value[1] == 0
            result[k].append(value[0])
    return result

I2 = matrix([[1, 0], [0, 1]])
X = matrix([[0, 1], [1, 0]])
Y = matrix([[0, (0, -1)], [(0, 1), 0]])
Z = matrix([[1, 0], [0, -1]])
letters = dict(zip('IXYZ', (I2, X, Y, Z)))
words = {a + b: kron(letters[a], letters[b]) for a in letters for b in letters}
I4 = words['II']
labels = [s for s in words if s != 'II']
for a, A in words.items():
    assert dagger(A) == A
    assert mmul(A, A) == I4
    for b, B in words.items():
        assert trace_product(A, B) == ((4 if a == b else 0), 0)
print('Pauli fidelity: all 16 words Hermitian, involutory and trace-orthogonal')

# The source's displayed rho^(2) at gamma = 1/3 has numerator S, denominator 6.
S = matrix([[0, 0, 0, 0],
            [0, 2, (1, 1), (-1, -1)],
            [0, (1, -1), 2, (-1, -1)],
            [0, (-1, 1), (-1, 1), 2]])
assert dagger(S) == S and trace(S) == (6, 0)
minors = principal_minors(S)
assert all(v >= 0 for row in minors.values() for v in row)
assert [sum(minors[k]) for k in range(1, 5)] == [6, 6, 0, 0]
source_expectations = {}
for label in labels:
    t = trace_product(S, words[label])
    assert t[1] == 0 and abs(t[0]) == 2
    source_expectations[label] = Fraction(t[0], 6)
assert trace_product(S, S) == (24, 0)
print('Source fidelity: rho^(2)(1/3), all 15 principal minors nonnegative, trace = 1')
print('Source fidelity: charpoly = t^2 * (t^2 - t + 1/6), purity = 2/3')
print('Source Pauli expectations:', ', '.join(f'{s}:{v}' for s, v in source_expectations.items()))

# The Mermin square is checked as operators, before enumerating sign choices.
contexts = [('XI', 'IX', 'XX', 1), ('IY', 'YI', 'YY', 1),
            ('XY', 'YX', 'ZZ', 1), ('XI', 'IY', 'XY', 1),
            ('IX', 'YI', 'YX', 1), ('XX', 'YY', 'ZZ', -1)]
assert Counter(s for c in contexts for s in c[:3]) == Counter({s: 2 for s in
    ('XI', 'IX', 'XX', 'IY', 'YI', 'YY', 'XY', 'YX', 'ZZ')})
certificates = {}
for ci, (p, q, r, eta) in enumerate(contexts):
    P, Q, R = (words[s] for s in (p, q, r))
    assert mmul(P, Q) == mmul(Q, P)
    assert mmul(P, R) == mmul(R, P)
    assert mmul(Q, R) == mmul(R, Q)
    assert mmul(mmul(P, Q), R) == scale(eta, I4)
    for a, b, c in product((-1, 1), repeat=3):
        if a * b * c != -eta:
            continue
        T = madd(madd(madd(I4, scale(-a, P)), scale(-b, Q)), scale(-c, R))
        assert dagger(T) == T and trace(T) == (4, 0)
        assert mmul(T, T) == scale(4, T)
        assert mmul(madd(I4, scale(-a, P)), madd(I4, scale(-b, Q))) == T
        assert all(v >= 0 for row in principal_minors(T).values() for v in row)
        certificates[ci, a, b, c] = T
print('Contexts: six commuting product identities and 24 PSD rank-one projectors checked exactly')

# Gray-code enumeration builds every signed Hamiltonian by one exact update.
H = scale(0, I4)
for s in labels:
    H = madd(H, words[s])
previous = 0
mismatch_histogram = Counter()
for index in range(1 << 15):
    mask = index ^ (index >> 1)
    if index:
        bit = (mask ^ previous).bit_length() - 1
        old_sign = -1 if previous & (1 << bit) else 1
        H = madd(H, scale(-2 * old_sign, words[labels[bit]]))
    previous = mask
    signs = {s: (-1 if mask & (1 << j) else 1) for j, s in enumerate(labels)}
    mismatches = [ci for ci, (p, q, r, eta) in enumerate(contexts)
                  if signs[p] * signs[q] * signs[r] == -eta]
    assert mismatches
    mismatch_histogram[len(mismatches)] += 1
    ci = mismatches[0]
    p, q, r, _ = contexts[ci]
    T = certificates[ci, signs[p], signs[q], signs[r]]
    assert trace(H) == ZERO and dagger(H) == H
    assert trace_product(T, H) == (-12, 0)
    assert trace_product(H, H) == (60, 0)
print('All 32768 sign assignments checked; mismatch-count histogram:', dict(sorted(mismatch_histogram.items())))
print('For every sign assignment: Tr(T)=4, Tr(T H)=-12; tau=T/4 is PSD')
print('Hence Tr(tau rho_gamma)=(1-3*gamma)/4; positivity forces gamma <= 1/3')
print('For every sign assignment: Tr(H)=0, Tr(H^2)=60; purity=(1+15*gamma^2)/4')

v = ((1, -1), ONE, ONE, ZERO)
w = (ONE, ZERO, ZERO, (0, 1))
N = madd(outer(v), outer(w))
Hplus = scale(0, I4)
for s in labels:
    Hplus = madd(Hplus, words[s])
assert scale(2, N) == madd(scale(3, I4), Hplus)
assert dagger(N) == N and trace(N) == (6, 0)
assert all(trace_product(N, words[s]) == (2, 0) for s in labels)
assert all(x >= 0 for row in principal_minors(N).values() for x in row)
assert Fraction(trace_product(N, N)[0], 36) == Fraction(2, 3)
assert Fraction(1 + 15 * Fraction(1, 3) ** 2, 4) == Fraction(2, 3)
print('Attainment: rho_*=(v v^*+w w^*)/6=(3 I+sum(P))/12; all 15 expectations = 1/3')
print('Exact optima: gamma=1/3 and purity=2/3')
print('ALL_CHECKS_OK')
