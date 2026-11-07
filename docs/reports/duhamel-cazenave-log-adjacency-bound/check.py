#!/usr/bin/env python3
"""Exact spectral checks for the Duhamel--Cazenave Table 3 logarithmic row.

Only integers, fractions and symbolic polynomial algebra are used. Exponential
series and log(1+x) <= x are the analytic implications of the exact certificates;
this program is independent evidence, not a Lean certificate.
"""

from fractions import Fraction
from math import factorial

import sympy as sp


t = sp.Symbol("t")


def connected(a):
    seen = {0}
    frontier = [0]
    while frontier:
        u = frontier.pop()
        for v in range(a.rows):
            if a[u, v] and v not in seen:
                seen.add(v)
                frontier.append(v)
    return len(seen) == a.rows


def matrices(m):
    vertices = [(s, i) for s in range(2) for i in range(m)]
    a = sp.Matrix(2 * m, 2 * m, lambda u, v:
                  int(u != v and (vertices[u][0] == vertices[v][0]
                                  or vertices[u][1] == vertices[v][1])))
    degrees = [sum(a.row(u)) for u in range(2 * m)]
    l = sp.diag(*degrees) - a
    return a, l, degrees


def prism_check(m):
    a, l, degrees = matrices(m)
    assert a == a.T and all(a[i, i] == 0 for i in range(2 * m))
    assert connected(a) and degrees == [m] * (2 * m)
    assert sum(degrees) // 2 == m * m
    one = sp.ones(m, 1)
    columns = [one.col_join(one), one.col_join(-one)]
    eigs = [m, m - 2]
    differences = [sp.eye(m)[:, i] - sp.eye(m)[:, m - 1]
                   for i in range(m - 1)]
    columns += [x.col_join(x) for x in differences]
    eigs += [0] * (m - 1)
    columns += [x.col_join(-x) for x in differences]
    eigs += [-2] * (m - 1)
    basis = sp.Matrix.hstack(*columns)
    assert basis.det() != 0
    assert a * basis == basis * sp.diag(*eigs)
    lap_eigs = [m - e for e in eigs]
    assert l * basis == basis * sp.diag(*lap_eigs)
    assert sp.expand(a.charpoly(t).as_expr()
                     - t ** (m - 1) * (t + 2) ** (m - 1)
                     * (t - m) * (t - (m - 2))) == 0
    assert sp.expand(l.charpoly(t).as_expr()
                     - t * (t - 2) * (t - m) ** (m - 1)
                     * (t - (m + 2)) ** (m - 1)) == 0
    assert sorted(eigs, reverse=True)[1] == m - 2
    assert max(lap_eigs) == m + 2
    if m == 7:
        print("SEED m=7: vertices=14 edges=49 degrees=7 basis_det=" + str(basis.det()))
        print("SEED adjacency spectrum:", sorted(eigs, reverse=True))
        print("SEED Laplacian spectrum:", sorted(lap_eigs))
        print("SEED characteristic polynomials: t^6*(t+2)^6*(t-7)*(t-5); "
              "t*(t-2)*(t-7)^6*(t-9)^6")


def fidelity_check():
    # Source equation (1): lambda_max(A) <= lambda_max(L) - 1.
    # Reproduce equality on K_n; do not claim to reproduce the paper's timings.
    for n in range(2, 11):
        a = sp.ones(n) - sp.eye(n)
        l = (n - 1) * sp.eye(n) - a
        assert sp.expand(a.charpoly(t).as_expr()
                         - (t - (n - 1)) * (t + 1) ** (n - 1)) == 0
        assert sp.expand(l.charpoly(t).as_expr()
                         - t * (t - n) ** (n - 1)) == 0
        adjacency_roots = sp.Poly(a.charpoly(t).as_expr(), t).ground_roots()
        laplacian_roots = sp.Poly(l.charpoly(t).as_expr(), t).ground_roots()
        assert max(adjacency_roots) == max(laplacian_roots) - 1
    print("FIDELITY: source equation (1) is an equality on K_n, 2 <= n <= 10; "
          "exact adjacency and D-A characteristic polynomials verified")


def family_certificates():
    # Generic action on (a*1,b*1), (x,x), (x,-x), with sum(x)=0.
    m, a, b, x = sp.symbols("m a b x", real=True)
    assert sp.expand(((m - 1) * a + b) - m * a).subs(b, a) == 0
    assert sp.expand(((m - 1) * a + b) - (m - 2) * a).subs(b, -a) == 0
    assert -x + x == 0 and -x - x == -2 * x
    assert 1 + 1 + (m - 1) + (m - 1) == 2 * m
    # For every integer m>=7, 2 ln(m+2) <= 2 ln9 + 2(m-7)/9.
    # exp(5/2) strictly exceeds its positive cubic Taylor partial sum.
    lower = sum(Fraction(5, 2) ** k / factorial(k) for k in range(4))
    assert lower == Fraction(443, 48) and lower > 9
    slack = sp.expand((m - 2) - (5 + sp.Rational(2, 9) * (m - 7)))
    assert sp.expand(slack - sp.Rational(7, 9) * (m - 7)) == 0
    print("FAMILY: invariant-subspace formulas checked symbolically; "
          "exact matrices/eigenbases/characteristic polynomials checked for 2 <= m <= 14")
    print("EXP CERTIFICATE: sum(k=0..3) (5/2)^k/k! = 443/48 > 9")
    print("FAMILY LOG CERTIFICATE: (m-2)-(5+2*(m-7)/9)=7*(m-7)/9 >= 0 for m>=7")
    print("CONCLUSION (using exp Taylor lower bound and log(1+x)<=x): "
          "m-2 > 2*ln(m+2) for every integer m>=7")


def tree_check():
    # The optional tree witness has two joined centers, each with d leaves.
    d = sp.Symbol("d", integer=True, positive=True)
    q = sp.Matrix([[0, 1, d, 0], [1, 0, 0, d],
                   [1, 0, 0, 0], [0, 1, 0, 0]])
    r = sp.Matrix([[d + 1, -1, -d, 0], [-1, d + 1, 0, -d],
                   [-1, 0, 1, 0], [0, -1, 0, 1]])
    assert sp.expand(q.charpoly(t).as_expr()
                     - (t*t - t - d) * (t*t + t - d)) == 0
    assert sp.expand(r.charpoly(t).as_expr()
                     - t * (t - d - 1) * (t*t - (d + 3)*t + 2)) == 0
    leaves = 256
    n = 2 + 2 * leaves
    a = sp.zeros(n)
    a[0, 1] = a[1, 0] = 1
    for side in range(2):
        for i in range(leaves):
            v = 2 + side * leaves + i
            a[side, v] = a[v, side] = 1
    degrees = [sum(a.row(u)) for u in range(n)]
    assert connected(a) and sum(degrees) // 2 == n - 1 == 513
    l = sp.diag(*degrees) - a
    columns = [sp.eye(n)[:, 0], sp.eye(n)[:, 1]]
    for side in range(2):
        columns.append(sp.Matrix([int(2 + side*leaves <= i < 2 + (side+1)*leaves)
                                  for i in range(n)]))
    u = sp.Matrix.hstack(*columns)
    assert a * u == u * q.subs(d, leaves)
    assert l * u == u * r.subs(d, leaves)
    # The remaining independent leaf differences have eigenvalues A=0, L=1.
    # Verify every column equation without a dense 514x514 characteristic polynomial.
    for side in range(2):
        anchor = 2 + side*leaves + leaves - 1
        for i in range(leaves - 1):
            v = 2 + side*leaves + i
            assert a[:, v] - a[:, anchor] == sp.zeros(n, 1)
            expected = sp.zeros(n, 1)
            expected[v] = 1
            expected[anchor] = -1
            assert l[:, v] - l[:, anchor] == expected
    assert 1025 > 32**2 and 67073 == 259**2 - 8 < 259**2
    lower = sum(Fraction(7)**k / factorial(k) for k in range(6))
    assert lower == Fraction(9893, 30) and lower > 259
    print("OPTIONAL TREE: n=514 edges=513; quotient and all 510 leaf-difference equations verified")
    print("TREE spectra: adjacency roots of (t^2-t-256)*(t^2+t-256), plus 0^510; "
          "Laplacian roots of t*(t-257)*(t^2-259*t+2), plus 1^510")
    print("TREE bounds: (sqrt(1025)-1)/2 > 31/2 > 14; "
          "(259+sqrt(67073))/2 < 259 < exp(7), since 9893/30 > 259")


if __name__ == "__main__":
    fidelity_check()
    for size in range(2, 15):
        prism_check(size)
    family_certificates()
    tree_check()
    print("ALL_CHECKS_OK")
