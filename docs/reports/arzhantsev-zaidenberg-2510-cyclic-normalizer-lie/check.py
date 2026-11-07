#!/usr/bin/env python3
"""Exact arithmetic checks for Question 2 of arXiv:2510.17223.

The checks use integers, Fraction, and SymPy's exact polynomial ring.  They
verify the source's cyclic-weight calculation and the proposed jet witness
family; the script is an independent consistency check, not a proof of the
ind-group tangent statement.
"""

from fractions import Fraction
from math import gcd

import sympy as sp


x, y = sp.symbols("x y")


def inverse_mod(e: int, d: int) -> int:
    return pow(e, -1, d)


def witness_data(d: int, e: int):
    assert d > 1 and 1 <= e < d and gcd(d, e) == 1
    ep = inverse_mod(e, d)
    k = (e * ep - 1) // d
    a, b = divmod(d, e)
    m = e + ep - 1
    return ep, k, a, b, m


def check_weight_fidelity():
    # Proposition 5.2 / Lemma 5.10: the two triangular families are invariant.
    for d in range(2, 31):
        for e in range(1, d):
            if gcd(d, e) != 1:
                continue
            ep = inverse_mod(e, d)
            for j in range(4):
                # y^(e+jd) d/dx has coefficient weight e+jd-e = jd.
                assert (e + j * d - e) % d == 0
                # x^(e'+jd) d/dy has coefficient weight (e'+jd)e-1.
                assert ((ep + j * d) * e - 1) % d == 0


def truncate(poly, degree):
    poly = sp.Poly(sp.expand(poly), x, y, domain=sp.ZZ)
    kept = 0
    for (i, j), coeff in poly.terms():
        if i + j < degree:
            kept += int(coeff) * x**i * y**j
    return sp.expand(kept)


def compose(f, g, degree):
    # f after g, with exact integer coefficients.
    return (truncate(f[0].subs({x: g[0], y: g[1]}), degree),
            truncate(f[1].subs({x: g[0], y: g[1]}), degree))


def check_jet_constraint(d: int, e: int, ep: int):
    m = e + ep - 1
    # Several nontrivial shears, including higher d-periodic terms.
    maps = [
        (x + y**e * (1 + 2 * y**d), y),
        (x, y + x**ep * (2 - x**d)),
        (x + 3 * y ** (e + d), y),
        (x, y - x ** (ep + d)),
    ]
    current = (x, y)
    for g in maps:
        current = compose(g, current, m)
        for coordinate in current:
            p = sp.Poly(coordinate, x, y, domain=sp.ZZ)
            for (i, j), coeff in p.terms():
                if i > 0 and j > 0:
                    assert i + j >= m, (d, e, ep, coordinate, i, j, m)


def check_counterexample_family():
    bad = []
    for d in range(2, 101):
        for e in range(1, d):
            if gcd(d, e) != 1:
                continue
            ep, k, a, b, m = witness_data(d, e)
            if k >= 2:
                # Lemma 2: e' >= 2a+1, and Lemma 3's mixed degree is below m.
                assert ep >= 2 * a + 1, (d, e, ep, a, b)
                assert a + b + 1 < m, (d, e, ep, a, b, m)
                assert (a * e + b) % d == 0
                D = sp.expand(x**a * y**b * ((b + 1) * x - (a + 1) * y))
                assert D != 0
                bad.append((d, e, ep, k, a, b, m))
    assert bad, "the k >= 2 family must contain examples"
    # The paper's explicit cubic missing direction.
    ep, k, a, b, m = witness_data(4, 3)
    assert (ep, k, a, b, m) == (3, 2, 1, 1, 5)
    assert (sp.expand((b + 1) * x ** (a + 1) * y**b),
            sp.expand(-(a + 1) * x**a * y ** (b + 1))) == (2 * x**2 * y, -2 * x * y**2)
    return bad


def check_arithmetic_classification():
    # Source's positive cases are exactly k = 0 (ee'=1) or k = 1 (ee'=d+1).
    for d in range(2, 101):
        for e in range(1, d):
            if gcd(d, e) != 1:
                continue
            ep = inverse_mod(e, d)
            k = (e * ep - 1) // d
            assert (k == 0) == (e * ep == 1)
            assert (k == 1) == (e * ep == d + 1)


def main():
    check_weight_fidelity()
    check_arithmetic_classification()
    bad = check_counterexample_family()
    # The first five pairs keep the exact truncated compositions small while
    # covering the requested (4,3) case and several k >= 2 instances.
    for d, e, ep, *_ in bad[:5]:
        check_jet_constraint(d, e, ep)
    print(f"checked {len(bad)} counterexample-family pairs for d <= 100")
    print("fidelity: cyclic weights and triangular normalizer generators")
    print("witness: (d,e,e')=(4,3,3), D=2*x*y*(x*d/dx-y*d/dy)")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
