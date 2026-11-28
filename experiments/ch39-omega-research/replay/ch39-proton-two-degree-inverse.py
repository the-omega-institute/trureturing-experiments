"""Exact identifiability tests for two shell sites of unequal degree."""

from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp


def positive_on_unit_interval(function, variable):
    """Map the open unit interval to the positive ray and certify coefficient signs."""
    x = sp.Symbol("x", positive=True)
    numerator, denominator = sp.cancel(function.subs(variable, x/(1+x))).as_numer_denom()
    for poly in (numerator, denominator):
        coefficients = sp.Poly(poly, x).all_coeffs()
        assert all(c >= 0 for c in coefficients) and any(c > 0 for c in coefficients)


def main():
    replay = Path(__file__).resolve().parent
    fold = runpy.run_path(str(replay / "ch39-proton-shell-overlap-fold.py"))
    w = fold["w"]
    ell, h = fold["bond_coefficients"]()
    # At rho=1/2 the two normalized coefficients are ell/u² and h/u⁴.
    # Their quotient removes the unknown scale u. All inequalities are strict.
    positive_on_unit_interval(-ell, w)
    positive_on_unit_interval(-h, w)
    positive_on_unit_interval(-sp.diff(sp.cancel(h/ell**2), w), w)

    # Unnormalized degree differences must survive at order t⁴ at every pair.
    # The 3x3 derivative uses (rho,u,w), and constants S'_AB and z are positive.
    z1, z2, u = sp.symbols("z1 z2 u", positive=True)
    normalized = sp.Matrix([ell/u**2, h/u**4])
    determinant = sp.factor(normalized.jacobian([u, w]).det())
    target = -2*(ell*sp.diff(h, w)-2*h*sp.diff(ell, w))/u**7
    assert sp.cancel(determinant-target) == 0
    shell = runpy.run_path(str(replay / "ch39-proton-overlap-identifiability-graphs.py"))
    vertices, _ = shell["shell_graph"]()
    edges = [(a, b) for a, b in combinations(vertices, 2)
             if sum((x-y)**2 for x, y in zip(a, b)) <= 804]
    assert len(edges) == 135
    degrees = [sum(vertex in edge for edge in edges) for vertex in vertices]
    assert set(degrees) == {10,14,16,18}
    raw_factor = (z2-z1)*z1
    for a, b in combinations(sorted(set(degrees)), 2):
        assert raw_factor.subs({z1:a, z2:b}) > 0
        print(f"Certified degree pair {a}, {b}")


if __name__ == "__main__":
    main()
