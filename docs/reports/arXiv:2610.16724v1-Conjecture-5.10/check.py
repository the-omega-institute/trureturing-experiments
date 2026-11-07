#!/usr/bin/env python3
"""Exact rational checks for Améndola–Riffo–Schmitz Conjecture 5.10.

Dependency: SymPy 1.14.0. No floating point or modular-rank evidence is used.
Pointwise Jacobian ranks are lower bounds on generic rank, not upper bounds.
"""

from itertools import product
from math import comb
from random import Random

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def s(d, r):
    return comb(r, 2) if r <= d else d * r - comb(d + 1, 2)


def expected(d, m, n):
    assert m > 0 and n > 0
    return d * m * n - s(d, m) - s(d, n) - s(d, (m - 1) * (n - 1))


def b(r):
    """e e^T + K: diagonal 1, strict upper triangle 2, lower triangle 0."""
    return sp.Matrix(r, r, lambda i, j: 1 if i == j else 2 if i < j else 0)


def moment(r):
    return sp.Matrix(r, r, lambda i, j: sp.Rational(j + 1, i + j + 2))


def pair_cores(m, n):
    return sp.kronecker_product(b(m), b(n)), sp.kronecker_product(b(m), b(n).T)


def exact_rank(matrix):
    return DomainMatrix.from_Matrix(matrix).convert_to(sp.QQ).rank()


def jacobian(a, cores):
    """Derivative of the pair A -> (A M_12 A^T, A M_21 A^T).

    These cores multiply the actual integral coordinates by 4. This nonzero
    scalar preserves rank. All 2*d*d output coordinates are retained.
    """
    d, q = a.shape
    rows = []
    for core in cores:
        left, right = core * a.T, a * core
        for i, j in product(range(d), repeat=2):
            rows.append([
                (left[t, j] if k == i else 0) + (right[i, t] if k == j else 0)
                for k, t in product(range(d), range(q))
            ])
    return sp.Matrix(2 * d * d, d * q, lambda i, j: rows[i][j])


def check_fidelity():
    # Corollary labelled cor:productMembranes: direct polynomial integration.
    x, y = sp.symbols("x y")
    for m in range(1, 4):
        integrated = sp.Matrix(m, m, lambda i, j:
            sp.integrate(sp.integrate(
                (i + 1) * x**i * (j + 1) * y**j, (x, 0, y)), (y, 0, 1)))
        assert integrated == moment(m)
    print("FIDELITY moment-curve integral S[p,q]=q/(p+q), 1<=p,q<=3: OK")

    # Example labelled ex:momentMembrane: the displayed sigma^(21)_11.
    a, c, e, f = sp.symbols("a c e f")
    row = sp.Matrix([[a, c, e, f]])
    value = (row * sp.kronecker_product(moment(2), moment(2).T) * row.T)[0]
    source = (a*a + c*c + e*e + f*f) / 4 + (
        a*c + a*e + c*f + e*f) / 2 + sp.Rational(4, 9)*a*f + sp.Rational(5, 9)*c*e
    assert sp.expand(value - source) == 0
    print("FIDELITY ex:momentMembrane displayed sigma^(21)_11 polynomial: OK")

    # Differentiate symbolically once to check the Jacobian implementation.
    aa = sp.Matrix(2, 4, sp.symbols("x0:8"))
    cores = pair_cores(2, 2)
    outputs = sp.Matrix([z for core in cores for z in aa * core * aa.T])
    assert outputs.jacobian(list(aa)) == jacobian(aa, cores)
    print("FIDELITY analytic quadratic-map Jacobian vs symbolic derivative: OK")


def check_ranks():
    # Table labelled fig:dimsMdmn gives the lower triangles for d=3 and d=4.
    source_tables = {
        3: [[3], [5, 10], [6, 13, 15], [6, 14, 15, 15]],
        4: [[4], [7, 14], [9, 19, 24], [10, 22, 25, 26], [10, 23, 25, 26, 26]],
    }
    checked = 0
    fidelity_entries = 0
    for d, m, n in product(range(1, 5), range(1, 6), range(1, 6)):
        cores = pair_cores(m, n)
        ranks = []
        for sample in range(3):
            rng = Random(261016724 + 100000*d + 1000*m + 10*n + sample)
            point = sp.Matrix(d, m*n, lambda i, j: rng.randrange(-7, 8))
            ranks.append(exact_rank(jacobian(point, cores)))
        target = expected(d, m, n)
        assert max(ranks) == target, (d, m, n, ranks, target)
        assert all(r <= target for r in ranks)
        print(f"RANK_Q d={d} m={m} n={n} samples={','.join(map(str, ranks))} expected={target}")
        checked += 1
        if d in source_tables and max(m, n) <= len(source_tables[d]):
            source = source_tables[d][max(m, n)-1][min(m, n)-1]
            assert max(ranks) == source
            fidelity_entries += 1
    # The empty target space d=0 is a genuine defined boundary.
    for m, n in product(range(1, 6), repeat=2):
        assert expected(0, m, n) == 0
        assert exact_rank(jacobian(sp.zeros(0, m*n), pair_cores(m, n))) == 0
    assert all(s(d, 0) == 0 for d in range(5))
    print(f"SUMMARY positive_d_cases={checked} samples_per_case=3 rational_ranks=300")
    print(f"FIDELITY source dimension table entries (including reflected entries)={fidelity_entries}: OK")
    print("BOUNDARIES d=0 with 1<=m,n<=5: 25 zero-dimensional cases; s_d(0)=0: OK")
    print("SCOPE m=0 or n=0 excluded: source uniform-grid denominators are undefined")
    print("LIMITATION sampled ranks do not certify a generic upper bound or the universal conjecture")


def main():
    print(f"EXACT_ARITHMETIC SymPy={sp.__version__} field=QQ seed=261016724")
    check_fidelity()
    check_ranks()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
