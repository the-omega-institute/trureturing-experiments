#!/usr/bin/env python3
"""Exact checks for Lin--Miao--Yu--Sheng Conjecture 5.11 (F2 uniqueness).

Requires Python 3 and SymPy. Matrices have integer/rational entries. No numeric
approximation is used. Finite checks support, but do not prove, the infinite
family argument in the preregistration.
"""
from fractions import Fraction
from itertools import combinations
import sympy as sp

X = sp.Symbol('x')


def graph(n, edges):
    e = {tuple(sorted((u, v))) for u, v in edges}
    assert all(0 <= u < v < n for u, v in e)
    a = sp.zeros(n)
    for u, v in e:
        a[u, v] = a[v, u] = 1
    d = [int(sum(a.row(i))) for i in range(n)]
    reached = {0}
    while True:
        enlarged = reached | {v for u in reached for v in range(n) if a[u, v] == 1}
        if enlarged == reached:
            break
        reached = enlarged
    return a, sp.diag(*d) - a, e, sorted(d, reverse=True), len(reached) == n


def gmn(n, m):
    r = m - n + 1
    return graph(n, [(0, i) for i in range(1, n)] + [(1, i) for i in range(2, r + 2)])


def h(n):
    return graph(n, [(0, i) for i in range(1, n)] + [(1, 2), (3, 4)])


def characteristic_check(matrix, eigenvalues):
    expected = sp.Poly(sp.prod(X - sp.sympify(v) for v in eigenvalues), X)
    actual = matrix.charpoly(X).as_poly()
    assert actual == expected, (actual, expected)
    return sp.factor(actual.as_expr())


def exact_centered_norms(lap, edge_count, spectrum):
    n = lap.rows
    shift = Fraction(2 * edge_count, n)
    centered_eigs = [Fraction(int(v)) - shift for v in spectrum]
    li = lap - sp.Rational(shift.numerator, shift.denominator) * sp.eye(n)
    characteristic_check(li, [sp.Rational(v.numerator, v.denominator) for v in centered_eigs])
    # Independent matrix check of all singular-value squares, not just trace.
    singular = sorted([abs(v) for v in centered_eigs], reverse=True)
    characteristic_check(li.T * li, [sp.Rational(v.numerator, v.denominator) ** 2 for v in singular])
    assert all(sp.sqrt(sp.Rational(v.numerator, v.denominator) ** 2) ==
               sp.Rational(v.numerator, v.denominator) for v in singular)
    return singular, [sum(singular[:k], Fraction(0)) for k in range(1, n + 1)]


def main():
    print('Arithmetic: Python integers, fractions.Fraction, SymPy exact matrices/polynomials')
    print('SymPy:', sp.__version__)
    # Source fidelity: introduction, regular graphs satisfy ||LI||_Fk = ||A||_Fk.
    for n in range(3, 11):
        adj, lap, edges, deg, connected = graph(n, combinations(range(n), 2))
        assert connected and deg == [n - 1] * n
        spectrum = [n] * (n - 1) + [0]
        characteristic_check(lap, spectrum)
        characteristic_check(adj, [n - 1] + [-1] * (n - 1))
        singular, norms = exact_centered_norms(lap, len(edges), spectrum)
        assert lap - (n - 1) * sp.eye(n) == -adj
        assert norms == [Fraction(n + k - 2) for k in range(1, n + 1)]
        print(f'FIDELITY K_{n}: regular identity for every k=1,...,{n}; F2={norms[1]}')

    # Figure-1 reconstruction throughout every small admissible (n,m).
    cases = 0
    for n in range(3, 10):
        for m in range(n, 2 * n - 2):
            adj, lap, edges, deg, connected = gmn(n, m)
            r = m - n + 1
            spectrum = [n, r + 2] + [2] * (r - 1) + [1] * (n - r - 2) + [0]
            assert connected and len(edges) == m and n <= m <= 2 * n - 3
            assert deg.count(1) == 2 * n - m - 3
            characteristic_check(lap, spectrum)
            cases += 1
    print(f'FIGURE_FIDELITY: {cases} admissible (n,m) pairs, 3<=n<=9; correct pendant count 2n-m-3')

    # Two separately constructed cones, including independent Gram checks.
    for n in range(6, 19):
        b, hh = gmn(n, n + 1), h(n)
        spectra = ([n, 4, 2] + [1] * (n - 4) + [0],
                   [n, 3, 3] + [1] * (n - 4) + [0])
        results = []
        for data, spec in zip((b, hh), spectra):
            adj, lap, edges, deg, connected = data
            assert connected and len(edges) == n + 1 and n <= len(edges) <= 2 * n - 3
            characteristic_check(lap, spec)
            singular, norms = exact_centered_norms(lap, len(edges), spec)
            assert singular[:2] == [Fraction(n) - Fraction(2 * (n + 1), n), Fraction(2 * (n + 1), n)]
            assert norms[1] == n
            results.append(norms)
        assert b[3] == [n - 1, 3, 2, 2] + [1] * (n - 4)
        assert hh[3] == [n - 1, 2, 2, 2, 2] + [1] * (n - 5)
        assert b[3] != hh[3]  # isomorphism preserves the degree multiset
        assert results[0][2] == Fraction(n + 2) - Fraction(2, n)
        assert results[1][2] == Fraction(n + 1) + Fraction(2, n)
        assert results[1][2] < results[0][2]
        print(f'FAMILY n={n}, m={n+1}: connected, admissible, nonisomorphic; '
              f'F2(B)=F2(H)={n}; F3(B)={results[0][2]}, F3(H)={results[1][2]}')
        if n == 6:
            print('n=6 chi_B =', characteristic_check(b[1], spectra[0]))
            print('n=6 chi_H =', characteristic_check(hh[1], spectra[1]))
    print('Scope: finite checks only; no exhaustive extremal-graph classification is asserted by this program.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
