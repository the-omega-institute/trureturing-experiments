#!/usr/bin/env python3
"""Exact rational checks for Evert--Jury--McCullough, arXiv:2607.13251, Q1.

No floating point, external packages, random sampling, or Lean proof is used.
The finite tests supplement the general construction; they do not certify all n.
"""
from fractions import Fraction as Q
from itertools import product


def mat(n, f):
    return [[Q(f(i, j)) for j in range(n)] for i in range(n)]


def eye(n):
    return mat(n, lambda i, j: i == j)


def add(*args):
    return [[sum(a[i][j] for a in args) for j in range(len(args[0]))]
            for i in range(len(args[0]))]


def scale(c, a):
    return [[c * x for x in row] for row in a]


def mul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def kron(a, b):
    n, m = len(a), len(b)
    return mat(n * m, lambda i, j: a[i // m][j // m] * b[i % m][j % m])


def outer(v):
    return [[x * y for y in v] for x in v]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def rank(a):
    b = [row[:] for row in a]
    r = 0
    for j in range(len(b[0])):
        pivot = next((i for i in range(r, len(b)) if b[i][j]), None)
        if pivot is None:
            continue
        b[r], b[pivot] = b[pivot], b[r]
        b[r] = [x / b[r][j] for x in b[r]]
        for i in range(r + 1, len(b)):
            c = b[i][j]
            b[i] = [x - c * y for x, y in zip(b[i], b[r])]
        r += 1
    return r


def psd_pivots(a):
    """Exact symmetric elimination: congruence to a nonnegative diagonal."""
    n = len(a)
    assert all(a[i][j] == a[j][i] for i, j in product(range(n), repeat=2))
    b = [row[:] for row in a]
    pivots = []
    for k in range(n):
        d = b[k][k]
        assert d >= 0, (k, d)
        pivots.append(d)
        if d == 0:
            assert all(b[k][j] == 0 for j in range(k, n))
        else:
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    b[i][j] -= b[i][k] * b[k][j] / d
    return pivots


def determinant(a):
    """Independent row elimination, with swaps, for a determinant cross-check."""
    b = [row[:] for row in a]
    d = Q(1)
    for k in range(len(b)):
        p = next((i for i in range(k, len(b)) if b[i][k]), None)
        if p is None:
            return Q(0)
        if p != k:
            b[k], b[p] = b[p], b[k]
            d = -d
        v = b[k][k]
        d *= v
        for i in range(k + 1, len(b)):
            c = b[i][k] / v
            for j in range(k + 1, len(b)):
                b[i][j] -= c * b[k][j]
    return d


def basis(n):
    return [mat(n, lambda a, b: a == i and b == j)
            for i, j in product(range(n), repeat=2)]


def choi(n, phi):
    """Source order (i,s,a),(j,t,b), computed directly from all matrix units."""
    units = basis(n)
    c = mat(n ** 3, lambda i, j: 0)
    for i, j, s, t in product(range(n), repeat=4):
        out = phi(units[n * i + j], units[n * s + t])
        for a, b in product(range(n), repeat=2):
            c[(i * n + s) * n + a][(j * n + t) * n + b] = out[a][b]
    return c


def schur_weighted(p, a, b):
    return mat(len(p), lambda i, j: p[i][j] * a[i][j] * b[i][j])


def diagonal_support(n):
    vs = [[Q(int(k == (i * n + i) * n + i)) for k in range(n ** 3)]
          for i in range(n)]
    w = [sum(v[k] for v in vs) for k in range(n ** 3)]
    return add(*(outer(v) for v in vs), outer(w))


def source_fidelity():
    # Example 5.10: the two diagonal entries are doubled; the others are not.
    p = [[Q(2), Q(1)], [Q(1), Q(2)]]
    t = [[Q(1, 2), Q(1)], [Q(1), Q(1, 2)]]
    phi = lambda a, b: schur_weighted(p, a, b)
    units = basis(2)
    for a in units:
        assert phi(t, a) == a == phi(a, t)
    for a, b in product(units, repeat=2):
        assert phi(a, b) == [[2 * a[0][0] * b[0][0], a[0][1] * b[0][1]],
                             [a[1][0] * b[1][0], 2 * a[1][1] * b[1][1]]]
    c = choi(2, phi)
    assert c == diagonal_support(2)
    assert rank(c) == 2
    assert sum(x > 0 for x in psd_pivots(c)) == 2
    assert determinant(t) == Q(-3, 4)
    print('FIDELITY Example 5.10: all 64 Choi entries, both unit identities, '
          'PSD rank 2, det(unit)=-3/4: OK')


def expected_settlement(n):
    ident = eye(n)
    if n == 1:
        phi = lambda a, b: [[a[0][0] * b[0][0]]]
        c = choi(1, phi)
        assert c == [[Q(1)]]
        for a in basis(1):
            assert phi(ident, a) == a == phi(a, ident)
        assert rank(c) == 1 and determinant(c) == 1
        print('SETTLEMENT n=1: scalar product, unit 1, Choi rank 1: OK')
        return
    ones = mat(n, lambda i, j: 1)
    p = add(ident, ones)
    t = add(ones, scale(Q(-1, 2), ident))
    h = add(scale(2 * n - 1, ident), scale(-1, ones))
    k = add(scale(n, ident), scale(-1, ones))
    f = lambda a: trace(mul(h, a))
    phi = lambda a, b: add(schur_weighted(p, a, b), scale(f(a) * f(b), ident))
    assert trace(mul(h, t)) == 0
    assert mul(h, t) == scale(Q(2 * n - 1, 2), add(ones, scale(-1, ident)))
    assert mat(n, lambda i, j: p[i][j] * t[i][j]) == ones
    for a in basis(n):
        assert phi(t, a) == a == phi(a, t)
    # K = sum_{i<j} (e_i-e_j)(e_i-e_j)^T is a rational Gram certificate.
    edges = [[Q(int(a == i) - int(a == j)) for a in range(n)]
             for i in range(n) for j in range(i + 1, n)]
    assert k == add(*(outer(v) for v in edges))
    assert h == add(scale(n - 1, ident), k)
    # Independent definition of C by matrix-unit evaluation vs. closed formula.
    c = choi(n, phi)
    gram = diagonal_support(n)
    htensor = kron(kron(h, h), ident)
    assert c == add(gram, htensor)
    bound = (n - 1) ** 2
    gap = add(c, scale(-bound, eye(n ** 3)))
    expanded_gap = add(gram,
                       scale(n - 1, kron(kron(k, ident), ident)),
                       scale(n - 1, kron(kron(ident, k), ident)),
                       kron(kron(k, k), ident))
    assert gap == expanded_gap
    psd_pivots(gap)
    pivots = psd_pivots(c)
    assert all(x > 0 for x in pivots)
    assert rank(c) == n ** 3
    d = determinant(c)
    d_from_pivots = Q(1)
    for x in pivots:
        d_from_pivots *= x
    assert d == d_from_pivots
    if n == 2:
        assert d == 22113
    print(f'SETTLEMENT n={n}: {n ** 6} Choi entries; {2 * n ** 2} unit-basis identities; '
          f'Tr(HT)=0; K Gram identity; C>={bound}I; rank={n ** 3}; det(C)={d}: OK')


def main():
    source_fidelity()
    for n in (1, 2, 3):
        expected_settlement(n)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
