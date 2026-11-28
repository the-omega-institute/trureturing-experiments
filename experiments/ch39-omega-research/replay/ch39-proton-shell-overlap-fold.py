"""Exact tests for an overlap-changing fold of three contact ratios at small hopping."""

import json
import runpy
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
w, rho = sp.symbols("w rho")
half = sp.Rational(1, 2)


def bond_coefficients():
    source = runpy.run_path(str(HERE / "ch39-proton-two-site-exceptional-orders.py"))
    moments = source["insertion_moments"]()
    series = {}
    for key, poly in moments.items():
        # T=t^4. Set u=1 and v=w after taking the local insertion derivatives.
        coefficients = [sp.S.Zero] * 3
        for (_, j, r, k), c in poly.to_dict().items():
            if k < 3:
                coefficients[k] += int(c)*w**j*rho**r
        series[key] = coefficients

    def product(a, b, k):
        return sum(a[i]*b[k-i] for i in range(k+1))

    def quotient(i, j, left, right):
        numer = [product(series[0, 0], series[i, j], k) for k in range(3)]
        denom = [product(series[left], series[right], k) for k in range(3)]
        static = sp.cancel(numer[0]/denom[0])
        n = [c.subs(rho, half) for c in numer]
        d = [c.subs(rho, half) for c in denom]
        x0 = static.subs(rho, half)
        x4 = sp.cancel((n[1]-x0*d[1])/d[0])
        x8 = sp.cancel((n[2]-x0*d[2]-x4*d[1])/d[0])
        return static, x4, x8

    ab = quotient(1, 1, (1, 0), (0, 1))
    aa = quotient(2, 0, (1, 0), (1, 0))
    bb = quotient(0, 2, (0, 1), (0, 1))
    assert sp.cancel(aa[0]-bb[0]) == 0 and sp.cancel(aa[1]-bb[1]) == 0
    slope = sp.diff(ab[0], rho).subs(rho, half)
    assert slope > 0
    return sp.factor(ab[1]/slope), sp.factor(bb[2]-aa[2])


def main():
    data = json.loads((HERE / "ch39-proton-overlap-jacobian-coordination-law.json").read_text())
    terms = [(int(k.split(",")[0]), int(k.split(",")[1]), int(v)) for k, v in data["P"].items()]
    determinant = sp.Poly(sum(c*w**i*half**j for i, j, c in terms), w)
    lo, hi = sp.Rational("0.4000162"), sp.Rational("0.4000163")
    assert determinant.count_roots(0, 1) == 1
    assert determinant.eval(lo)*determinant.eval(hi) < 0
    assert determinant.gcd(determinant.diff()).degree() == 0

    shift, diff8 = bond_coefficients()
    # Holding the leading transverse readout u^-4 h(w) fixed gives
    # du/dw = u h'/(4h); the derivative of u^-2 l(w) then has this numerator.
    drift = sp.factor(2*diff8*sp.diff(shift, w)-shift*sp.diff(diff8, w))
    for function in (shift, diff8, drift):
        numerator, denominator = sp.fraction(function)
        assert determinant.gcd(sp.Poly(numerator, w)).degree() == 0
        assert determinant.gcd(sp.Poly(denominator, w)).degree() == 0
    print("fold at rho=1/2, w in", lo, hi)
    print("ell =", shift)
    print("h =", diff8)
    print("2 h ell' - ell h' =", drift)


if __name__ == "__main__":
    main()
