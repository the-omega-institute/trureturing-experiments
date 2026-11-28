"""Exact triangle response and noise orders for the finite contact probe."""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import perm
from pathlib import Path
import runpy

import flint
import sympy as sp


HERE = Path(__file__).resolve().parent
RING = flint.fmpq_mpoly_ctx.get(("w", "rho"))


def amplitude(terms, rho=None, der=(0, 0), uv=1):
    """Insert first. Set u=1, v=w and evaluate the full polynomial in rho."""
    i, j = der
    coefficients = Counter()
    for (n11, n22, n21, _), c in terms.items():
        if i <= n11 and j <= n22:
            coefficients[n22-j+n21, n21] += c*perm(n11, i)*perm(n22, j)
    return RING.from_dict({k: flint.fmpq(c.numerator, c.denominator)
                           for k, c in coefficients.items() if c})


def triangle_coefficient(graphs, amplitudes):
    tri, count = graphs["triangle_t6"](amplitudes, None, 1, amplitude=amplitude)
    assert count == 3216
    fixed, fixed_count = graphs["triangle_t6"](amplitudes, Fraction(1, 2), 1)
    assert fixed_count == count
    for d, poly in tri.items():
        sliced = poly.subs({"rho": flint.fmpq(1, 2)})
        projected = flint.fmpq_poly([0])
        for (j, k), c in sliced.to_dict().items():
            assert k == 0
            projected[int(j)] = c
        assert projected == fixed[d]
    static = {d: amplitude(amplitudes[0, 0], der=d) for d in graphs["DERS"]}
    w, r = sp.symbols("w rho")

    def expression(p):
        return sum(sp.Rational(int(c.p), int(c.q))*w**int(i)*r**int(j)
                   for (i, j), c in p.to_dict().items())

    # At order six only the connected triangle contributes. Zeroth order is
    # the product of three single sites; the two unobserved sites carry no insertions.
    s = {d: expression(v) for d, v in static.items()}
    correction = {d: expression(v) for d, v in tri.items()}
    n0 = s[0, 0]*s[1, 1]
    d0 = s[1, 0]*s[0, 1]
    ratio = sp.cancel(n0/d0)
    n6 = correction[0, 0]*s[1, 1]+s[0, 0]*correction[1, 1]
    d6 = correction[1, 0]*s[0, 1]+s[1, 0]*correction[0, 1]
    return sp.factor((n6-ratio*d6)/(s[0, 0]**2*d0*sp.diff(ratio, r)))


def main():
    graphs = runpy.run_path(str(HERE / "ch39-proton-overlap-identifiability-graphs.py"))
    points, _ = graphs["shell_graph"]()
    adjacent = {p: {q for q in points if p != q and sum((a-b)**2 for a,b in zip(p,q)) <= 804}
                for p in points}
    types = Counter((len(neighbors), sum(q in adjacent[p] for p, q in combinations(neighbors, 2)))
                    for neighbors in adjacent.values())
    assert types == {(10, 45): 6, (14, 83): 3, (16, 100): 6, (18, 117): 4}
    determinants = [sp.Matrix([[1, z, n] for z, n in triple]).det()
                    for triple in combinations(sorted(types), 3)]
    assert determinants == [-8, -16, -8, 0]
    print("Site (degree, triangle count) types and multiplicities:", sorted(types.items()))
    print("Three-site determinants:", determinants)
    amps = graphs["site_amplitudes"]()
    triangle = triangle_coefficient(graphs, amps)
    print("Triangle displacement", triangle, flush=True)
    low = runpy.run_path(str(HERE / "ch39-proton-two-degree-signs.py"))
    _, ell, _, _ = low["bond_responses"]()
    w, r = sp.symbols("w rho")
    cross = sp.factor(3*sp.diff(ell,w)*triangle-2*ell*sp.diff(triangle,w))
    half = sp.Rational(1, 2)
    tiny = sp.Rational(1, 100)
    assert cross.subs({w: half, r: half}) == -sp.Rational(101786027, 52502068096)
    assert cross.subs({w: tiny, r: tiny}) < 0 < cross.subs({w: half, r: tiny})
    print("3 ell_w m - 2 ell m_w =", cross)
    print("At w=rho=1/2 the coefficient is -101786027/52502068096.")
    print("At rho=1/100 a zero lies between w=1/100 and w=1/2.")


if __name__ == "__main__":
    main()
