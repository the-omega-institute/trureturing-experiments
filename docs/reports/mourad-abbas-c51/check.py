#!/usr/bin/env python3
"""Exact checks for Mourad--Abbas, arXiv:1310.1273v1, Conjecture 5.1.

Finite checks and symbolic identities support a preregistered proof outline.
They do not quantify over all real matrices and are not a Lean proof of DS.
Run with Python 3 and sympy 1.14.0. No floating-point computation is used.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import sympy as sp


def matrix(n, f):
    return [[Q(f(i, j)) for j in range(n)] for i in range(n)]


def eye(n):
    return matrix(n, lambda i, j: int(i == j))


def uniform(n):
    return matrix(n, lambda i, j: Q(1, n))


def c_matrix(n):
    assert n >= 2
    return matrix(n, lambda i, j: 0 if i == j else Q(1, n - 1))


def affine(a, b, s):
    n = len(a)
    return matrix(n, lambda i, j: (1 - s) * a[i][j] + s * b[i][j])


def add(a, b, scale=Q(1)):
    return matrix(len(a), lambda i, j: a[i][j] + scale * b[i][j])


def mul(a, b):
    n = len(a)
    return matrix(n, lambda i, j: sum(a[i][h] * b[h][j] for h in range(n)))


def blocks(parts):
    sizes = [len(a) for a in parts]
    n = sum(sizes)
    ans = matrix(n, lambda i, j: 0)
    start = 0
    for a, size in zip(parts, sizes):
        for i in range(size):
            for j in range(size):
                ans[start + i][start + j] = a[i][j]
        start += size
    return ans


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def symmetric_stochastic(a):
    n = len(a)
    assert all(a[i][j] >= 0 and a[i][j] == a[j][i] for i in range(n) for j in range(n))
    assert all(sum(a[i]) == 1 for i in range(n))
    assert all(sum(a[i][j] for i in range(n)) == 1 for j in range(n))


def polynomial_check(a, lam):
    n = len(a)
    assert mul(add(a, eye(n), -1), add(a, eye(n), -lam)) == matrix(n, lambda i, j: 0)


def characteristic(a):
    return sp.Matrix(a).charpoly(x).as_expr().expand()


def four_cycle(n):
    h = matrix(n, lambda i, j: 0)
    for i, j, sign in [(0, 1, 1), (2, 3, 1), (1, 2, -1), (3, 0, -1)]:
        h[i][j] = h[j][i] = Q(sign)
    assert trace(h) == 0 and all(sum(row) == 0 for row in h)
    return h


def nonextreme_check(a, h, eps):
    assert eps > 0
    plus, minus = add(a, h, eps), add(a, h, -eps)
    symmetric_stochastic(plus)
    symmetric_stochastic(minus)
    assert plus != minus and affine(plus, minus, Q(1, 2)) == a


x = sp.Symbol('x')
# Fidelity: Section 2 explicitly gives this full spectrum of C_n.
for n in range(2, 9):
    c = c_matrix(n)
    symmetric_stochastic(c)
    assert characteristic(c) == sp.expand((x - 1) * (x + sp.Rational(1, n - 1)) ** (n - 1))
print('FIDELITY: C_n spectrum (1, -1/(n-1) repeated n-1), n=2..8: OK')
# Fidelity: the source warns that direct sums of DS matrices need not be DS.
for k in range(3, 7):
    a = blocks([uniform(2), uniform(2 * k - 2)])
    b = blocks([uniform(k), uniform(k)])
    assert characteristic(a) == characteristic(b)
    assert Q(1, k) in {v for row in b for v in row}
    assert Q(1, k) not in {v for row in a for v in row}
print('FIDELITY: J_2 + J_(2k-2) and J_k + J_k cospectral, different entries, k=3..6: OK')

# Symbolic scalar identities at arbitrary dimension and parameter.
k, t, s = sp.symbols('k t s', positive=True)
lam = (k * t - 1) / (k - 1)
w = (1 - t) / (k - 1)
assert sp.factor(lam + (1 - lam) / k - t) == 0
assert sp.factor((1 - lam) / k - w) == 0
assert sp.factor(t + (k - 1) * w - 1) == 0
assert sp.factor(lam + 1 / (k - 1) - k * t / (k - 1)) == 0
assert sp.factor(1 + k * lam - (k ** 2 * t - 1) / (k - 1)) == 0
# P_s = (A - (1-s)I)/s: diagonal, off-diagonal, and cycle witnesses.
pdiag, poff = (t - 1 + s) / s, w / s
assert sp.factor(pdiag + (k - 1) * poff - 1) == 0
assert sp.factor(poff - poff / 2 - poff / 2) == 0
print('SYMBOLIC: affine block, row sum, eigenvalue bounds, size obstruction, I-ray perturbation identities: OK')

cases = [(k0, r, Q(j, 4 * k0 * k0)) for k0 in range(4, 8) for r in (2, 3) for j in (1, 2, 3)]
cases.append((4, 2, Q(1, 40)))
for k0, r, t0 in cases:
    n = k0 * r
    lam0 = (k0 * t0 - 1) / (k0 - 1)
    assert 0 < t0 < Q(1, k0 * k0)
    assert -Q(1, k0 - 1) < lam0 < -Q(1, k0) < 0
    d = affine(c_matrix(k0), eye(k0), t0)
    a = blocks([d] * r)
    symmetric_stochastic(a)
    assert trace(a) == n * t0 > 0
    polynomial_check(a, lam0)
    e = matrix(n, lambda i, j: (a[i][j] - lam0 * int(i == j)) / (1 - lam0))
    assert e == blocks([uniform(k0)] * r)
    assert mul(e, e) == e and trace(e) == r
    # Enumerate ALL sorted component-size partitions with r parts summing to n.
    # Flor's theorem / the planned rigidity lemma reduces arbitrary cospectral B
    # to this form; that general reduction is not supplied by this computation.
    possible = [parts for parts in combinations_with_replacement(range(1, n + 1), r)
                if sum(parts) == n and all(1 + (m - 1) * lam0 >= 0 for m in parts)]
    assert possible == [(k0,) * r]
    # Cross-block zero and within-block positivity exclude [I_n,C_n].
    assert a[0][k0] == 0 and a[0][1] > 0
    h = four_cycle(n)
    nonextreme_check(a, h, a[0][1] / 2)
    # Full algebraic I-ray argument is above; also check its end and interior points.
    for s0 in (1 - t0, 1 - t0 / 2, Q(1)):
        p = matrix(n, lambda i, j: (a[i][j] - (1 - s0) * int(i == j)) / s0)
        symmetric_stochastic(p)
        assert affine(eye(n), p, s0) == a
        nonextreme_check(p, h, p[0][1] / 2)
print('FAMILY: 25 exact instances, k=4..7, r=2,3, three interior parameters plus (4,2,1/40): OK')
print('RIGIDITY REDUCTION: unique admissible component-size partition in every tested instance: OK')
print('SEGMENTS: support obstruction; exact nonextremality witnesses for A and three I-ray points each: OK')

k0, r, t0 = 4, 2, Q(1, 40)
a = blocks([affine(c_matrix(k0), eye(k0), t0)] * r)
assert characteristic(a) == sp.expand((x - 1) ** 2 * (x + sp.Rational(3, 10)) ** 6)
nonextreme_check(a, four_cycle(8), Q(1, 100))
print('NOMINATED: N=8, t=1/40, diagonal=1/40, within-block off-diagonal=13/40, trace=1/5, lambda=-3/10: OK')

# The upper bound is strict: at t=1/k^2, unequal block sizes k-1,k+1
# give an actual non-permutation-similar cospectral stochastic matrix.
for k0 in range(4, 7):
    lam0 = -Q(1, k0)
    a = blocks([affine(c_matrix(k0), eye(k0), Q(1, k0 * k0))] * 2)
    b = blocks([matrix(m, lambda i, j: lam0 * int(i == j) + (1 - lam0) / m)
                for m in (k0 - 1, k0 + 1)])
    symmetric_stochastic(a)
    symmetric_stochastic(b)
    assert characteristic(a) == characteristic(b)
    assert {a[i][i] for i in range(2 * k0)} != {b[i][i] for i in range(2 * k0)}
print('BOUNDARY CONTROL: at t=1/k^2, unequal blocks are cospectral mates, k=4..6: OK')
print('SCOPE: finite exact checks plus symbolic identities; arbitrary-real-matrix DS and full family remain prose proof obligations.')
print('ALL_CHECKS_OK')
