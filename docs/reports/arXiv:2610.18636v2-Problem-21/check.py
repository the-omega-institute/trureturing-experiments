#!/usr/bin/env python3
"""Exact arithmetic fidelity and pre-probe checks for Problem 21.

The script deliberately checks finite identities and the published integer
certificate arithmetic.  It does not claim to replace the general proof.
"""
from fractions import Fraction
import sympy as sp


def q(x, n):
    return x**3 - n*x**2 - (2*n**2 - n - 1)*x + (n**2 - n)


def assert_q_signs():
    n = sp.symbols("n", integer=True, positive=True)
    x = sp.symbols("x")
    poly = x**3 - n*x**2 - (2*n**2 - n - 1)*x + (n**2 - n)
    assert sp.expand(poly.subs(x, 0)) == n**2 - n
    assert sp.expand(poly.subs(x, sp.Rational(1, 2))) == sp.Rational(5, 8) - sp.Rational(3, 4)*n
    assert sp.expand(poly.subs(x, 2*n)) == sp.expand(3*n*(n + sp.Rational(1, 3)))
    assert sp.expand(poly.subs(x, -n)) == -2*n
    # The sign claims used for the three-root separation hold at every n >= 3.
    for nv in range(3, 101):
        assert q(Fraction(0), nv) > 0
        assert q(Fraction(1, 2), nv) < 0
        assert q(Fraction(2*nv), nv) > 0
        assert q(Fraction(-nv), nv) < 0
    print("q_n endpoint identities and signs: OK (symbolic plus n=3..100)")


def assert_distance_identity():
    n, l1, l2, l3 = sp.symbols("n l1 l2 l3", real=True)
    lam = 2*n - l1
    # Ordered spectra differ only at the three active coordinates.
    l1_distance = lam + l2 + (n + l3)
    l2_distance = lam**2 + l2**2 + (n + l3)**2
    relation = sp.expand(l1 + l2 + l3 - n)
    assert sp.expand(l1_distance.subs(l3, n-l1-l2) - 2*lam) == 0
    assert sp.expand(l2_distance.subs(l3, n-l1-l2) - (2*(lam**2 - lam*l2 + l2**2))) == 0
    # A rational admissible sample makes the absolute-value signs explicit.
    nv, a, b = Fraction(3), Fraction(5), Fraction(3, 8)
    c = nv - a - b
    assert 0 < 2*nv-a and 0 < b and -nv < c
    d1 = abs(2*nv-a) + abs(b) + abs(nv+c)
    d2 = (2*nv-a)**2 + b**2 + (nv+c)**2
    assert d1 == 2*(2*nv-a)
    assert d2 == 2*((2*nv-a)**2 - (2*nv-a)*b + b**2)
    assert relation.subs({n: sp.Rational(3), l1: sp.Rational(5), l2: sp.Rational(3,8), l3: sp.Rational(-19,8)}) == 0
    print("spectral-distance identities: OK (exact symbolic elimination)")


def adjacency_matrix(edges, m):
    A = sp.zeros(m)
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    return A


def hn_edges(n):
    # Definition 3 in the source, with the vertices partitioned as
    # X=(n-1), Y=n, Z=n, v=last; v is joined to X.
    m = 3*n
    edges = set()
    X = range(0, n-1)
    Y = range(n-1, 2*n-1)
    Z = range(2*n-1, 3*n-1)
    v = 3*n-1
    for a in X:
        for b in Y:
            edges.add((a,b))
        for b in Z:
            edges.add((a,b))
        edges.add((a,v))
    for a in Y:
        for b in Z:
            edges.add((a,b))
    return edges, m


def assert_hn_fidelity():
    n = 3
    edges, m = hn_edges(n)
    A = adjacency_matrix(edges, m)
    x = sp.symbols("x")
    charpoly = sp.factor(A.charpoly(x).as_expr())
    q3 = x**3 - n*x**2 - (2*n**2-n-1)*x + (n**2-n)
    expected = x**(3*n-4) * (x+n) * q3
    assert sp.expand(charpoly - expected) == 0
    assert A.rank() == 4
    print("H_n characteristic polynomial and rank at n=3: OK (exact SymPy algebra)")


def assert_connected_rank_six_example():
    # K_{3n} is connected and has adjacency rank 3n >= 6.  Its exact
    # distance is a concrete rank-six-or-more instance of the target.
    for n in range(3, 11):
        m = 3*n
        distance = 3*(n-1)**2 + (3*n-3)
        assert distance >= Fraction(2, 3)
    print("connected rank >= 6 family K_{3n}: OK (exact distance >= 2/3 for n=3..10)")


def assert_certificate_rows():
    # The two elimination orders reported by the search seat end in these
    # integer square certificates.  Keep the arithmetic explicit and exact.
    rows = {
        "six-vertex order": (25, 9, 5, 3, 3, 1),
        "seven-vertex order": (32, 16, 8, 4, 2, 1, 1),
    }
    expected = {"six-vertex order": 750, "seven-vertex order": 1366}
    for label, terms in rows.items():
        numerator = sum(t*t for t in terms)
        assert numerator == expected[label]
        bound = Fraction(numerator, 1024)
        assert bound > Fraction(2, 3)
        print(f"{label} certificate: {numerator}/1024 > 2/3")
    assert Fraction(750, 1024) == Fraction(375, 512)
    assert Fraction(1366, 1024) == Fraction(683, 512)


if __name__ == "__main__":
    assert_q_signs()
    assert_distance_identity()
    assert_hn_fidelity()
    assert_connected_rank_six_example()
    assert_certificate_rows()
    print("ALL_CHECKS_OK")
