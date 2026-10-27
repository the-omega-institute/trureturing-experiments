#!/usr/bin/env python3
"""Exact checks for Elphick--Linz, journal Question 4.3 (arXiv v3 Question 2).

Requires Python 3 and SymPy. Every matrix entry, polynomial coefficient,
eigenvalue, energy and ratio used as evidence is an integer or rational.
The finite checks do not establish the unbounded-family assertion.
"""

from collections import deque
from itertools import combinations, product
import sys

import sympy as sp


X = sp.Symbol("x")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def exact_spectrum(matrix):
    """Extract every root with multiplicity from the independently computed charpoly."""
    polynomial = matrix.charpoly(X).as_poly()
    require(polynomial.domain == sp.ZZ, "nonintegral characteristic polynomial")
    spectrum = polynomial.ground_roots()
    require(sum(spectrum.values()) == matrix.rows, "omitted characteristic roots")
    require(all(value.is_Rational for value in spectrum), "nonrational eigenvalue")
    require(not polynomial.as_expr().has(sp.Float), "floating-point evidence")
    return polynomial, spectrum


def energies(spectrum):
    positive = sum((multiplicity * value**2 for value, multiplicity in spectrum.items()
                    if value > 0), sp.Integer(0))
    negative = sum((multiplicity * value**2 for value, multiplicity in spectrum.items()
                    if value < 0), sp.Integer(0))
    require(positive.is_Rational and negative.is_Rational, "inexact square energy")
    return positive, negative


def connected(matrix):
    if matrix.rows == 0:
        return False
    reached = {0}
    todo = deque([0])
    while todo:
        u = todo.popleft()
        for v in range(matrix.rows):
            if matrix[u, v] == 1 and v not in reached:
                reached.add(v)
                todo.append(v)
    return len(reached) == matrix.rows


def check_source_fidelity():
    # Section 4, immediately after the definition of squared spread:
    # K_n has s^+ - s^- = (n - 1)^2 - (n - 1) = (n - 1)(n - 2).
    for n in range(1, 13):
        matrix = sp.ones(n) - sp.eye(n)
        polynomial, spectrum = exact_spectrum(matrix)
        expected = sp.Poly((X - (n - 1)) * (X + 1)**(n - 1), X)
        require(polynomial == expected, f"K_{n} characteristic polynomial")
        positive, negative = energies(spectrum)
        require(positive == (n - 1)**2, f"K_{n} positive square energy")
        require(negative == n - 1, f"K_{n} negative square energy")
        require(positive - negative == (n - 1) * (n - 2), f"K_{n} squared spread")
    print("FIDELITY: complete graphs K_n, 1 <= n <= 12: exact source identity OK")


def check_petersen_family():
    pairs = [frozenset(s) for s in combinations(range(1, 6), 2)]
    base = sp.Matrix([[int(s.isdisjoint(t)) for t in pairs] for s in pairs])
    require(len(pairs) == 10, "wrong base order")
    require(base == base.T and all(base[i, i] == 0 for i in range(10)), "not simple")
    require(all(sum(base.row(i)) == 3 for i in range(10)), "wrong base degree")
    require(sum(base) == 30 and connected(base), "wrong edges or connectivity")
    require(base**2 == 2 * sp.eye(10) - base + sp.ones(10), "base matrix identity")
    polynomial, spectrum = exact_spectrum(base)
    require(spectrum == {sp.Integer(3): 1, sp.Integer(1): 5, sp.Integer(-2): 4},
            "base multiplicities")
    require(polynomial == sp.Poly((X - 3) * (X - 1)**5 * (X + 2)**4, X),
            "base characteristic polynomial")
    print("BASE: n=10, m=15, degree=3, connected; A^2 = 2I - A + J OK")
    print("BASE: charpoly=(x-3)(x-1)^5(x+2)^4; spectrum=3^1,1^5,(-2)^4")

    for t in range(1, 9):
        vertices = list(product(pairs, range(t)))
        # Construct the adjacency relation directly from disjoint two-subsets.
        matrix = sp.Matrix([[int(u[0].isdisjoint(v[0])) for v in vertices]
                            for u in vertices])
        n = len(vertices)
        require(n == 10 * t and matrix == matrix.T, f"t={t} order or symmetry")
        require(all(matrix[i, i] == 0 for i in range(n)), f"t={t} loop")
        require(all(sum(matrix.row(i)) == 3 * t for i in range(n)), f"t={t} degree")
        require(sum(matrix) == 30 * t**2 and connected(matrix), f"t={t} edges/connectivity")
        require(matrix == sp.kronecker_product(base, sp.ones(t)), f"t={t} Kronecker model")
        polynomial, spectrum = exact_spectrum(matrix)
        expected = {sp.Integer(3*t): 1, sp.Integer(t): 5, sp.Integer(-2*t): 4}
        if t > 1:
            expected[sp.Integer(0)] = 10 * (t - 1)
        require(spectrum == expected, f"t={t} full multiplicities")
        require(polynomial == sp.Poly((X-3*t)*(X-t)**5*(X+2*t)**4*X**(10*(t-1)), X),
                f"t={t} characteristic polynomial")
        positive, negative = energies(spectrum)
        require(positive == 14 * t**2 and negative == 16 * t**2, f"t={t} energies")
        spread = positive - negative
        require(positive + negative == sum(matrix), f"t={t} second moment")
        require(sp.Rational(-spread, n**2) == sp.Rational(1, 50), f"t={t} coefficient")
        print(f"FAMILY: t={t}, n={n}, m={15*t*t}, s+={positive}, s-={negative}, "
              f"spread={spread}, (-spread)/n^2=1/50; exact charpoly and connectivity OK")
    print("SCOPE: finite exact checks t=1,...,8; unbounded-family claim requires the paper argument")


def main():
    print(f"ENVIRONMENT: Python {sys.version.split()[0]}, SymPy {sp.__version__}")
    check_source_fidelity()
    check_petersen_family()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
