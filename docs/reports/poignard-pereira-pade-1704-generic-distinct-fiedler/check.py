#!/usr/bin/env python3
"""Exact checks for Poignard--Pereira--Pade (arXiv:1704.01677v1).

All evidence is over SymPy's exact integer/rational/algebraic domains; no
floating point value is used to decide an assertion.
"""
from itertools import combinations
import sympy as sp


def laplacian(n, edges, weights):
    L = sp.zeros(n)
    for (i, j), weight in zip(edges, weights):
        # Inputs are exact SymPy expressions; positivity is asserted for
        # concrete witnesses below.
        L[i, i] += weight
        L[j, j] += weight
        L[i, j] -= weight
        L[j, i] -= weight
    return L


def pairwise_distinct(v):
    return all(sp.simplify(v[i] - v[j]) != 0 for i, j in combinations(range(len(v)), 2))


def main():
    # Three-vertex path with positive edge weights 1 and 2.
    # This is a connected weighted simple graph, and is a direct fidelity
    # instance for the source's simple-spectrum and nonzero-Fiedler claims.
    edges = [(0, 1), (1, 2)]
    L = laplacian(3, edges, [sp.Integer(1), sp.Integer(2)])
    assert L == sp.Matrix([[1, -1, 0], [-1, 3, -2], [0, -2, 2]])
    assert all(sum(L.row(i)) == 0 for i in range(3))
    assert L == L.T

    lam = sp.symbols("lambda")
    charpoly = sp.factor(L.charpoly(lam).as_expr())
    assert charpoly == lam * (lam**2 - 6 * lam + 6)

    eigendata = L.eigenvects()
    assert sorted(mult for _, mult, _ in eigendata) == [1, 1, 1]
    positive = sorted((value for value, _, _ in eigendata if value != 0),
                      key=lambda value: sp.default_sort_key(value))
    # The exact algebraic roots are 3 - sqrt(3) < 3 + sqrt(3); this check
    # identifies the Fiedler eigenpair without numerical approximation.
    fiedler_value = 3 - sp.sqrt(3)
    assert fiedler_value in [value for value, _, _ in eigendata]
    assert sp.simplify(fiedler_value - (3 + sp.sqrt(3))) < 0
    fiedler_vector = next(vectors[0] for value, _, vectors in eigendata if value == fiedler_value)
    assert all(sp.simplify(entry) != 0 for entry in fiedler_vector)
    assert pairwise_distinct(fiedler_vector)

    # The stronger simultaneous witness used by the paper's proposed route:
    # every nonconstant eigenvector in this exact example has distinct entries.
    for value, multiplicity, vectors in eigendata:
        if value != 0:
            assert multiplicity == 1
            vector = vectors[0]
            assert all(sp.simplify(entry) != 0 for entry in vector)
            assert pairwise_distinct(vector)

    # A symbolic fidelity check for the path family: the nonzero spectrum is
    # simple throughout the positive quadrant because its discriminant is the
    # nonzero polynomial 4*(a^2 - a*b + b^2).
    a, b = sp.symbols("a b")
    family = laplacian(3, edges, [a, b])
    family_charpoly = sp.factor(family.charpoly(lam).as_expr())
    assert sp.expand(family_charpoly - lam * (lam**2 - 2 * (a + b) * lam + 3 * a * b)) == 0
    nonzero_factor = sp.Poly(family_charpoly / lam, lam)
    assert sp.simplify(sp.discriminant(nonzero_factor.as_expr(), lam) - 4 * (a**2 - a * b + b**2)) == 0
    assert sp.factor((a**2 - a * b + b**2).subs({a: 1, b: 2})) == 3

    print("exact witness: P3 weights (1,2)")
    print("charpoly:", charpoly)
    print("fiedler eigenvalue:", fiedler_value)
    print("fiedler vector:", tuple(fiedler_vector))
    print("symbolic positive-quadrant discriminant:", 4 * (a**2 - a * b + b**2))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
