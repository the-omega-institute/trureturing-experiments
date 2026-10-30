#!/usr/bin/env python3
"""Exact algebra checks for KMV arXiv:2610.22248v1, Question 31.

Checks support a proposed argument; they do not classify arbitrary evolution
algebras or prove non-realizability of an algebraic group. No floating-point
arithmetic is used. Requires Python 3 and SymPy.
"""
from itertools import product
import sympy as s


def unit(d, i, j):
    a = s.zeros(d)
    a[i, j] = 1
    return a


def derivation_equations(c):
    """Leibniz constraints, columns D[a,b], on a multiplication c[i][j][k]."""
    d = len(c)
    rows = []
    for i, j, k in product(range(d), repeat=3):
        row = [0] * (d * d)
        for a in range(d):
            row[k * d + a] += c[i][j][a]
            row[a * d + i] -= c[a][j][k]
            row[a * d + j] -= c[i][a][k]
        rows.append(row)
    return s.Matrix(rows)


def fidelity():
    # Proposition 14: e_i^2=u, u^2=u, all distinct products zero.
    readings = []
    for m in range(1, 8):
        d = m + 1
        c = [[[int(i == j and k == m) for k in range(d)]
              for j in range(d)] for i in range(d)]
        equations = derivation_equations(c)
        skew = [unit(d, i, j) - unit(d, j, i)
                for i in range(m) for j in range(i + 1, m)]
        expected = m * (m - 1) // 2
        assert d * d - equations.rank() == expected
        if skew:
            columns = s.Matrix.hstack(*(a.reshape(d * d, 1) for a in skew))
            assert columns.rank() == expected
            assert equations * columns == s.zeros(equations.rows, expected)
        else:
            assert expected == 0
        readings.append((m, expected))
    # Nontrivial exact orthogonal automorphisms for the m=2 example.
    q = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5)],
                  [s.Rational(4, 5), s.Rational(3, 5)]])
    for q2 in (q, s.diag(-1, 1)):
        assert q2.T * q2 == s.eye(2)
        g = s.diag(q2, 1)
        for i, j in product(range(3), repeat=2):
            before = s.Matrix([0, 0, int(i == j)])
            after = s.Matrix([0, 0, sum(g[a, i] * g[a, j] for a in range(3))])
            assert g * before == after
    print('FIDELITY Proposition 14: (m, dim Der(A_m)) =', readings)
    print('FIDELITY full Leibniz kernels equal the embedded skew spaces; exact orthogonal examples verified')


def tangent_identities():
    eps = s.Symbol('eps')
    for d in range(1, 4):
        D = s.Matrix(d, d, s.symbols('d:' + str(d * d)))
        F = s.eye(d) + eps * D
        inverse_product = F * (s.eye(d) - eps * D)
        for i, j in product(range(d), repeat=2):
            p = s.Poly(inverse_product[i, j], eps)
            assert p.nth(0) == int(i == j) and p.nth(1) == 0
        c = [[[s.Symbol(f'c_{i}_{j}_{k}') for k in range(d)]
              for j in range(d)] for i in range(d)]
        for i, j, k in product(range(d), repeat=3):
            residual = (sum(F[k, a] * c[i][j][a] for a in range(d))
                        - sum(F[a, i] * F[b, j] * c[a][b][k]
                              for a, b in product(range(d), repeat=2)))
            coefficient = s.Poly(residual, eps).nth(1)
            leibniz = (sum(D[k, a] * c[i][j][a] for a in range(d))
                       - sum(D[a, i] * c[a][j][k] + D[a, j] * c[i][a][k]
                             for a in range(d)))
            assert s.expand(coefficient - leibniz) == 0
            assert s.Poly(residual, eps).nth(0) == 0
        assert s.Poly(F.det(), eps).nth(1) == s.trace(D)
    print('TANGENT generic bilinear products d=1..3: dual-number equations = Leibniz; det(I+eps D) = 1+eps tr(D) mod eps^2')


def sl3_killing():
    n = 3
    basis = [unit(n, i, j) for i in range(n) for j in range(n) if i != j]
    basis += [unit(n, i, i) - unit(n, n - 1, n - 1) for i in range(n - 1)]
    cols = s.Matrix.hstack(*(a.reshape(n * n, 1) for a in basis))
    ads = []
    for a in basis:
        brackets = [(a * b - b * a).reshape(n * n, 1) for b in basis]
        ads.append(s.Matrix.hstack(*(cols.gauss_jordan_solve(b)[0] for b in brackets)))
    K = s.Matrix(8, 8, lambda i, j: s.trace(ads[i] * ads[j]))
    G = s.Matrix(8, 8, lambda i, j: 6 * s.trace(basis[i] * basis[j]))
    assert K == G and K.det() != 0
    print('SL3 exact adjoint Killing matrix: dimension=8, determinant=' + str(K.det()))
    allowed = [(m, m * (m - 1) // 2) for m in range(3, 10) if m * (m - 1) // 2 <= 8]
    assert allowed == [(3, 3), (4, 6)]
    assert not [(a, b) for a in range(3) for b in range(2) if 3 * a + 6 * b == 8]
    print('SL3 dimension obstruction: so factors of dimension <=8 have dimensions 3 or 6; no sum equals 8')


def family():
    n = s.Symbol('n', integer=True, positive=True)
    delta_even = s.factor(n**2 - 1 - (2*n - 2) * (2*n - 3) / 2)
    delta_odd = s.factor(n**2 - 1 - (2*n - 1) * (2*n - 2) / 2)
    assert s.expand(delta_even + (n - 1) * (n - 4)) == 0
    assert s.expand(delta_odd + (n - 1) * (n - 2)) == 0
    assert s.solve(delta_even, n) == [1, 4]
    assert s.solve(delta_odd, n) == [1, 2]
    # floor(m/2)=n-1 permits exactly m=2n-2 or 2n-1.
    simultaneous = []
    dimension_only = []
    for a in range(2, 257):
        for m in range(3, 2*a + 2):
            equal_rank = m // 2 == a - 1
            assert equal_rank == (m in (2*a - 2, 2*a - 1))
            equal_dim = a*a - 1 == m*(m - 1) // 2
            if equal_dim:
                dimension_only.append((a, m))
            if equal_rank and equal_dim:
                simultaneous.append((a, m))
    assert simultaneous == [(2, 3), (4, 6)]
    assert (11, 16) in dimension_only
    print('FAMILY exact all-n polynomial reductions:', delta_even, ';', delta_odd)
    print('FAMILY n=2..256 simultaneous (rank,dimension) matches:', simultaneous)
    print('CONTROL dimension alone has extra matches:', dimension_only)
    print('BOUNDARY n=2 and n=4 are not excluded; these checks do not assert their group realization')


def main():
    print('Arithmetic: integers and exact SymPy polynomial/rational/matrix algebra')
    fidelity()
    tangent_identities()
    sl3_killing()
    family()
    print('SCOPE classification, differentiation of natural isomorphisms, and maximal-Cartan rank formulas require the mathematical argument')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
