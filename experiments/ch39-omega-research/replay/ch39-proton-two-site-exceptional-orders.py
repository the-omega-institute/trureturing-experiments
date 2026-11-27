"""Exceptional hopping orders of the three-ratio Jacobian on a single bond.

Run with Python 3.12, python-flint and sympy. The insertion moments come from
the four-flavour partition polynomial in ch39-proton-current-overlap-hopping.py.
"""

import json
import runpy
from collections import defaultdict
from math import perm, prod
from pathlib import Path

import flint
import sympy as sp


HERE = Path(__file__).resolve().parent
RING = flint.fmpz_mpoly_ctx.get(("u", "v", "rho", "T"))
U, V, RHO, T = range(4)


def insertion_moments():
    source = runpy.run_path(str(HERE / "ch39-proton-current-overlap-hopping.py"))
    variables = [source[key] for key in ("u0", "u1", "v0", "v1", "x", "y", "t")]
    partition = source["Z_of"](source["generating_series"](), 4)
    terms = sp.Poly(partition, *variables, domain=sp.ZZ).terms()
    assert len(terms) == 225

    moments = {}
    for i, j in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)):
        coefficients = defaultdict(int)
        for (a, b, c, d, h, g, k), value in terms:
            assert g == h and k % 4 == 0
            assert a + b + c + d + 2*h + k//2 == 16
            if a < i or c < j:
                continue
            # Differentiate at site 0 holding x,y fixed, then set the sites
            # equal and impose x=rho*u*v, y=1, T=t^4.
            key = (a + b + h - i, c + d + h - j, h, k//4)
            coefficients[key] += int(value) * perm(a, i) * perm(c, j)
        moments[i, j] = RING.from_dict(coefficients)
    return moments


def jacobian_fraction(moments, point):
    a, b = moments[1, 0], moments[0, 1]
    rows, denominators = [], []
    for derivative, den in (((1, 1), a*b), ((2, 0), a*a), ((0, 2), b*b)):
        num = moments[0, 0] * moments[derivative]
        n, d = num.subs(point), den.subs(point)
        rows.append([d*num.derivative(j).subs(point) - n*den.derivative(j).subs(point)
                     for j in (U, V, RHO)])
        denominators.append(d**2)
    a, b, c = rows
    determinant = (a[0]*(b[1]*c[2] - b[2]*c[1]) + a[1]*(b[2]*c[0] - b[0]*c[2])
                   + a[2]*(b[0]*c[1] - b[1]*c[0]))
    return determinant, prod(denominators)


def coefficient(poly, order, variable, symbol):
    """Extract [T^order] from a slice with one remaining parameter."""
    coefficients = {}
    for powers, value in poly.to_dict().items():
        assert all(power == 0 for j, power in enumerate(powers) if j not in (T, variable))
        if powers[T] == order:
            coefficients[(int(powers[variable]),)] = int(value)
    return sp.Poly.from_dict(coefficients, symbol, domain=sp.ZZ)


def main():
    moments = insertion_moments()
    w, r = sp.symbols("w r")

    # On rho=0 leave w=v/u symbolic at u=1. The homogeneity checked above
    # restores a factor u^(-2-2k) to the coefficient of T^k.
    num, den = jacobian_fraction(moments, {"u": 1, "rho": 0})
    assert all(coefficient(num, k, V, w).is_zero for k in range(6))
    n6 = coefficient(num, 6, V, w)
    d0 = coefficient(den, 0, V, w)
    assert d0.terms()[0][1] > 0 and len(d0.terms()) == 1
    positive = (w**2 - 14)**2 + 60
    assert (n6.as_expr() * 30198988800*w**9 + d0.as_expr()*positive).expand() == 0
    print("[t^24] J(1,w,0) = -((w^2-14)^2+60)/(30198988800*w^9)")

    num, den = jacobian_fraction(moments, {"u": 1, "v": 1})
    assert all(coefficient(num, k, RHO, r).is_zero for k in range(4))
    n4 = coefficient(num, 4, RHO, r)
    n5 = coefficient(num, 5, RHO, r)
    d0 = coefficient(den, 0, RHO, r)
    data = json.loads((HERE / "ch39-proton-overlap-jacobian-coordination-law.json").read_text())
    diagonal = sp.Poly(sum(int(c)*r**int(key.split(",")[1]) for key, c in data["P"].items()), r)
    p = 1 + 24*r + 90*r**2 + 80*r**3 + 15*r**4
    q = 1 + 18*r + 45*r**2 + 20*r**3
    # Independent comparison with the published t^16 coefficient on w=1.
    scale = sp.Rational(data["C"])
    assert n4*sp.Poly(scale.q*q**8*p**2, r) == scale.p*d0*diagonal*sp.Poly(r**2, r)
    lo, hi = sp.Rational("0.929583"), sp.Rational("0.929584")
    assert diagonal.degree() == 26 and diagonal.count_roots(0, 1) == 1
    assert diagonal.eval(lo)*diagonal.eval(hi) < 0
    assert diagonal.gcd(n5).degree() == diagonal.gcd(d0).degree() == 0
    print("On w=1, P(1,r) has one root in (0,1), isolated by", lo, hi)
    print("At this root [t^16] J = 0 and [t^20] J != 0 (both gcd degrees zero)")


if __name__ == "__main__":
    main()
