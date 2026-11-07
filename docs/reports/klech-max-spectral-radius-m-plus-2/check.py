#!/usr/bin/env python3
"""Exact diagnostics for Klech, arXiv:2609.18367v1, Conjecture 5.15.

Python 3 and SymPy 1.14.0. No approximate eigenvalues are used.
These finite diagnostics are not a proof for arbitrary m.
"""
from collections import defaultdict
from itertools import combinations, permutations
import sympy as sp

X, Z = sp.symbols('X Z')


def connected(n, edges):
    rows = [0] * n
    for u, v in edges:
        rows[u] |= 1 << v
    for start in range(n):
        seen = 1 << start
        while True:
            new = seen
            for v in range(n):
                if seen & (1 << v):
                    new |= rows[v]
            if new == seen:
                break
            seen = new
        if seen != (1 << n) - 1:
            return False
    return True


def matrix(n, edges):
    es = set(edges)
    return sp.Matrix(n, n, lambda i, j: int((i, j) in es))


def bouquet(lengths):
    """Three directed cycles, common vertex 0, disjoint interiors."""
    edges, nxt = set(), 1
    for length in lengths:
        assert length >= 1
        vertices = [0] + list(range(nxt, nxt + length - 1))
        nxt += length - 1
        for u, v in zip(vertices, vertices[1:] + [0]):
            assert (u, v) not in edges
            edges.add((u, v))
    return nxt, edges


def relabelings(n, edges):
    return {frozenset((p[u], p[v]) for u, v in edges)
            for p in permutations(range(n))}


def perron(poly):
    """PF: for a nonnegative irreducible matrix, largest real root is rho.

    real_roots(radicals=False) isolates algebraic roots with exact rationals.
    SymPy exact algebraic comparisons refine isolating intervals as necessary.
    """
    roots = sp.Poly(poly, X).real_roots(radicals=False)
    assert roots and roots[-1] > 0
    return roots[-1]


def exhaustive(n):
    groups = defaultdict(list)
    for edges in combinations([(i, j) for i in range(n) for j in range(n)], n + 2):
        if connected(n, edges):
            cp = tuple(matrix(n, edges).charpoly(X).all_coeffs())
            groups[cp].append(frozenset(edges))
    expected_n, general_edges = bouquet([1, 2, n - 1])
    assert expected_n == n
    general_isos = relabelings(n, general_edges)
    general_max = perron(matrix(n, general_edges).charpoly(X).as_expr())
    # Source fidelity: Theorem 5.14's minimum polynomial 1-2z^(n-1)-z^n.
    source_min = perron(X**n - 2*X - 1)
    general_equal, loopless_equal, minimum_equal = [], [], []
    count, loopless_count = 0, 0
    if n >= 4:
        ln, loopless_edges = bouquet([2, 2, n - 2])
        assert ln == n
        loopless_isos = relabelings(n, loopless_edges)
        loopless_max = perron(matrix(n, loopless_edges).charpoly(X).as_expr())
    for coeffs, graphs in groups.items():
        rho = perron(sp.Poly.from_list(list(coeffs), X).as_expr())
        assert rho <= general_max
        assert source_min <= rho
        if rho == source_min:
            minimum_equal.extend(graphs)
        if rho == general_max:
            general_equal.extend(graphs)
        count += len(graphs)
        if n >= 4:
            loopless = [g for g in graphs if all(u != v for u, v in g)]
            if loopless:
                loopless_count += len(loopless)
                assert rho <= loopless_max
                if rho == loopless_max:
                    loopless_equal.extend(loopless)
    assert set(general_equal) == general_isos
    assert minimum_equal
    assert count == {3: 39, 4: 688}[n]
    print(f'm={n}: strongly connected labelled graphs={count}; characteristic polynomials={len(groups)}')
    print(f'  source fidelity: minimum rho is largest real root of X^{n}-2X-1; attaining graphs={len(minimum_equal)}')
    print(f'  clause 1: exact maximum and complete bouquet equality class; labelled maximizers={len(general_equal)}')
    if n >= 4:
        assert loopless_count == 316
        assert set(loopless_equal) == loopless_isos
        print(f'  clause 2: loopless graphs={loopless_count}; exact maximum and complete bouquet equality class; labelled maximizers={len(loopless_equal)}')


def polynomial_checks():
    # Section 4, HH type: T=1-z^k1-z^k2-z^(l+1).
    # det(I-zA) is computed independently, with exact integer coefficients.
    checked = 0
    for m in range(3, 13):
        tuples = [(2, 1, 1, m-2), (m-1, 1, 1, 1), (m-1, 2, 1, 0)]
        if m >= 4:
            tuples += [(2, 2, 1, m-3), (m-2, 2, 1, 1)]
        for k1, k2, t, ell in tuples:
            assert t == 1
            n, edges = bouquet([k1, k2, ell+1])
            assert n == m and len(edges) == m+2 and connected(n, edges)
            cp = matrix(n, edges).charpoly(X).as_expr()
            det = sp.Poly(sp.expand(Z**n * cp.subs(X, 1/Z)), Z).as_expr()
            hh = 1-Z**k1-Z**k2-Z**(ell+1)
            assert sp.expand(det-hh) == 0
            expected = (1-Z-Z**2-Z**(m-1) if min(k1, k2, ell+1) == 1
                        else 1-2*Z**2-Z**(m-2))
            assert sp.expand(det-expected) == 0
            checked += 1
    print(f'HH source formula and all conjectured tuple polynomials: {checked} exact cases, 3<=m<=12')
    r = sp.symbols('r')
    # Algebraic deficit used in the general two-self-return endpoint comparison.
    assert sp.expand((1-3*r+2*r**2+r**3)-((1-r)**3+r**2*(2*r-1))) == 0
    print('Endpoint deficit identity: exact polynomial expansion agrees')


def main():
    print('Exact arithmetic: Python integers and SymPy rational/algebraic root isolation; no floating-point evidence')
    polynomial_checks()
    for n in (3, 4):
        exhaustive(n)
    print('Scope: finite diagnostics only; the unbounded positive settlement requires the preregistered general argument')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
