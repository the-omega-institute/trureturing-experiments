"""Certify the two-degree inverse on 0<w,rho<1 using exact Bernstein coefficients."""

from math import comb
from pathlib import Path
import runpy

import sympy as sp
from sympy.polys.fields import field


def bond_responses():
    source = Path(__file__).with_name("ch39-proton-two-site-exceptional-orders.py")
    moments = runpy.run_path(str(source))["insertion_moments"]()
    fractions, w, rho = field("w,rho", sp.QQ)
    series = {}
    for key, poly in moments.items():
        coefficients = [fractions.zero] * 3
        # u=1, v=w; the insertion moments were differentiated before restricting.
        for (_, v_power, rho_power, order), value in poly.to_dict().items():
            if order <= 2:
                coefficients[order] += int(value)*w**int(v_power)*rho**int(rho_power)
        series[key] = coefficients

    def product(left, right):
        return [sum(left[j]*right[k-j] for j in range(k+1)) for k in range(3)]

    def ratio(insertion, first, second):
        n = product(series[0, 0], series[insertion])
        d = product(series[first], series[second])
        static = n[0]/d[0]
        fourth = (n[1]-d[1]*static)/d[0]
        eighth = (n[2]-d[2]*static-d[1]*fourth)/d[0]
        return static, fourth, eighth

    mixed = ratio((1, 1), (1, 0), (0, 1))
    aa = ratio((2, 0), (1, 0), (1, 0))
    bb = ratio((0, 2), (0, 1), (0, 1))
    assert aa[:2] == bb[:2]
    slope = mixed[0].diff(rho)
    shift = mixed[1]/slope
    h = bb[2]-aa[2]
    drift = 2*h*shift.diff(w)-shift*h.diff(w)
    return [f.as_expr() for f in (slope, shift, h, drift)]


def bernstein_bounds(poly):
    """Bounds on the tensor Bernstein coefficients over the closed unit square."""
    m, n = poly.degree_list()
    terms = poly.terms()
    coefficients = [
        sum(c*sp.Rational(comb(k, i), comb(m, i))*sp.Rational(comb(l, j), comb(n, j))
            for (i, j), c in terms if i <= k and j <= l)
        for k in range(m+1) for l in range(n+1)
    ]
    return min(coefficients), max(coefficients)


def main():
    w, rho = sp.symbols("w rho")
    slope, shift, h, drift = bond_responses()
    for factor in sp.fraction(slope):
        assert all(c > 0 for c in sp.Poly(factor, rho).coeffs())
        assert factor.subs(rho, 0) > 0

    q = 20*rho**3+45*rho**2+18*rho+1
    p = 15*rho**4+80*rho**3+90*rho**2+24*rho+1
    h_reduced = sp.Poly(sp.cancel(h*49152*w**4*q**3*p/(rho*(16-w**2))), w, rho)
    drift_reduced = sp.Poly(sp.cancel(drift*294912*w**6*q**3*p**2/rho**2), w, rho)
    assert h_reduced.degree_list() == (2, 10)
    assert drift_reduced.degree_list() == (4, 13)
    assert bernstein_bounds(h_reduced) == (-22127616, 0)
    assert bernstein_bounds(drift_reduced) == (-17774007552, -768)

    # A different normalization by the square of the shift would be singular.
    assert shift.subs({w: sp.Rational(9, 10), rho: sp.Rational(1, 100)}) > 0
    assert shift.subs({w: sp.Rational(9, 10), rho: sp.Rational(1, 2)}) < 0
    print("0<w,rho<1: S'_AB>0, h<0, 2h ell_w - ell h_w < 0.")
    print("Bernstein bounds after clearing positive factors:", bernstein_bounds(h_reduced),
          bernstein_bounds(drift_reduced))


if __name__ == "__main__":
    main()
