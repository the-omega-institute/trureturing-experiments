"""Exact checks for the convex Gibbs family in arXiv:2410.02470v2, Conjecture 1."""

from fractions import Fraction
from math import comb

import sympy as sp


def check_semicircle_moments():
    # Substitute x = 2t in the density of Definition 2. The remaining
    # integral is B(m + 1/2, 3/2), multiplied by 2 * 4**m / pi.
    # The base moment is 1 by the unit disk area formula.
    moment = Fraction(1)
    for m in range(9):
        assert moment == Fraction(comb(2 * m, m), m + 1)
        moment *= 4 * Fraction(2 * m + 1, 2 * m + 4)
    print("Fidelity: Definition 2 yields Catalan moments m = 0,...,8; odd moments vanish by symmetry.")


def check_boundary_and_curvature():
    u, v, e = sp.symbols("u v e", real=True)
    w = u - sp.I * v
    k = 1 / w + w + e * w**2
    p = 2 * e * u**3 + u**2 - 1
    height_squared = 1 / (1 + 2 * e * u) - u**2
    assert sp.factor(height_squared + p / (1 + 2 * e * u)) == 0
    assert sp.factor(sp.im(k) - v * (1 / (u**2 + v**2) - 1 - 2 * e * u)) == 0
    assert sp.factor(k.diff(u) - p.subs(u, w) / w**2) == 0
    assert tuple(p.subs(u, x) for x in (-2, -1, 0, 1)) == (3 - 16 * e, -2 * e, -1, 2 * e)
    assert 3 - 16 * Fraction(1, 16) == 2
    x = 2 * u + 4 * e * u**2 - e / (1 + 2 * e * u)
    re_boundary = u * (1 + 2 * e * u) + u + e * (u**2 - height_squared)
    assert sp.factor(re_boundary - x) == 0
    slope = 2 + 8 * e * u + 2 * e**2 / (1 + 2 * e * u)**2
    assert sp.factor(slope - x.diff(u)) == 0
    bound = Fraction(1, 16)
    denominator_lower = 1 - 4 * bound
    assert denominator_lower == Fraction(3, 4)
    assert 2 - 16 * bound == 1
    slope_upper = 2 + 8 * bound + 2 * bound**2 / denominator_lower**2
    assert slope_upper < 4
    potential = 2 * u**2 + sp.Rational(16, 3) * e * u**3
    potential += sp.log(1 + 2 * e * u) + 1 / (1 + 2 * e * u)
    assert sp.factor(potential.diff(u) - 2 * u * slope) == 0
    distance, endpoint, endpoint_value, endpoint_slope = sp.symbols(
        "distance endpoint endpoint_value endpoint_slope", real=True
    )
    extension = endpoint_value + 2 * endpoint * distance + distance**2 / endpoint_slope
    assert extension.subs(distance, 0) == endpoint_value
    assert extension.diff(distance).subs(distance, 0) == 2 * endpoint
    assert extension.diff(distance, 2) == 2 / endpoint_slope
    print(f"Family: exact boundary and endpoint identities, 1 < x' <= {slope_upper} < 4.")


def check_convolution():
    w, e, t = sp.symbols("w e t")
    plus = w + e * w**2
    minus = -plus.subs(w, -w)
    assert sp.expand(plus + minus) == 2 * w
    # Inversion of zG = 1 + G² + eG³ at infinity, using the
    # Lagrange coefficient formula for G(w) = w(1 + G(w)² + eG(w)³).
    phi = sp.Poly(1 + t**2 + e * t**3, t)
    moments = [(phi**(n + 1)).nth(n) / (n + 1) for n in range(4)]
    assert moments == [1, 0, 1, e]
    assert [m.subs(e, -e) for m in moments] == [1, 0, 1, -e]
    print("Family: exact moments (1, 0, 1, +/- e) and R-transform sum 2w.")


if __name__ == "__main__":
    check_semicircle_moments()
    check_boundary_and_curvature()
    check_convolution()
    print("ALL_CHECKS_OK")
