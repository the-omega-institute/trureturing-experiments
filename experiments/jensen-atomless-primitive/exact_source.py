"""Exact compact positive-law source and actual-node first-defect certificate."""
import argparse
import json
from pathlib import Path

import sympy as sp


def run():
    x, t, c, y = sp.symbols("x t c y")
    k0 = (1 - x**2)**5
    raw = sp.cancel(-(sp.diff(k0, x) + 2*sp.diff(k0, x, 3)/c
                      + sp.diff(k0, x, 5)/c**2)/x)
    mass = sp.factor(sp.integrate(raw, (x, -1, 1)))
    assert sp.cancel(mass - 512*(c**2 - 36*c + 504)/(63*c**2)) == 0
    bernstein = [
        10*(c*c-48*c+720)/c**2, 5*(c*c-36*c+600)/c**2,
        sp.Rational(15, 7)*(c*c-56)/c**2,
        sp.Rational(5, 7)*(c*c+120*c-3024)/c**2,
        (c*c+720*c-21840)/(7*c**2),
        60*(9*c-350)/(7*c**2), 120*(2*c-105)/(7*c**2),
        480/c**2, 3840/c**2,
    ]
    py = sp.Poly(raw, x).as_expr().subs(x*x, y)
    assert sp.expand(py - sum(bernstein[j]*sp.binomial(8,j)*y**j*(1-y)**(8-j)
                              for j in range(9))) == 0
    v = sp.symbols("v")
    for b in bernstein:
        assert all(a >= 0 for a in sp.Poly(sp.expand((c*c*b).subs(c,v+64)),v).all_coeffs())
        assert b.subs(c,64) > 0
    density = sp.Poly(raw.subs(c,64)/mass.subs(c,64),x)
    moment = lambda j: sum(2*a/sp.Integer(e[0]+2*j+1) for e,a in density.terms())
    coefficient = lambda j: (-1)**j*moment(j)/sp.factorial(2*j)
    falling = lambda n,j: sp.factorial(n)/sp.factorial(n-j)
    def row(n, shift=0):
        return sp.Poly(sum(falling(n,j)*sp.factorial(j+shift)/sp.factorial(j)
                           *coefficient(j+shift)*t**j/sp.Integer(n)**j
                           for j in range(n+1)),t)
    rows = []
    for n in range(1,14):
        poly = row(n)
        positive = int(poly.count_roots(0,sp.oo))
        negative = int(poly.count_roots(-sp.oo,0))
        gcd = sp.gcd(poly,poly.diff()).degree()
        assert positive == (n if n <= 12 else 11)
        assert negative == 0 and gcd == 0
        rows.append({"degree":n,"positive_roots":positive,"negative_roots":negative,
                     "gcd_with_derivative_degree":gcd})
    j, g = row(12), row(12,1)
    lo, hi = sp.Rational(131,2), sp.Rational(6551,100)
    assert j.count_roots(lo,hi) == 1 and j.eval(lo)*j.eval(hi) < 0
    signs = {}
    for name,poly,sign in [
        ("G12",g,-1), ("J12prime",j.diff(),1),
        ("G12_plus_7over5_J12prime",g+sp.Rational(7,5)*j.diff(),-1),
        ("G12_plus_3over2_J12prime",g+sp.Rational(3,2)*j.diff(),1),
    ]:
        count = int(poly.count_roots(lo,hi))
        endpoint_sign = int(sp.sign(poly.eval(lo)))
        assert count == 0 and endpoint_sign == sign
        signs[name] = {"root_count":count,"sign":endpoint_sign}
    # Check the source derivative normalization by exact coefficient extraction.
    kappa = sp.Rational(256,3157)
    bcoef = lambda j: (-1)**j/(4**j*sp.factorial(j)*sp.rf(sp.Rational(13,2),j))
    for j in range(15):
        rhs = bcoef(j)
        if j >= 1:
            rhs -= sp.Rational(1,32)*bcoef(j-1)
        if j >= 2:
            rhs += sp.Rational(1,4096)*bcoef(j-2)
        assert (j+1)*coefficient(j+1) == -kappa*rhs
    return {
        "source":{"c":64,"density":str(sp.factor(density.as_expr())),
                  "raw_mass":"41/9","derivative_kappa":"256/3157",
                  "moments_0_through_3":[str(moment(j)) for j in range(4)],
                  "positive_bernstein_at_64":[str(b.subs(c,64)) for b in bernstein]},
        "rows":rows,
        "first_defect":13,
        "actual_node":{"row":12,"interval":["131/2","6551/100"],"root_count":1,
                       "ratio_strict_interval":["-3/2","-7/5"],"sign_certificates":signs},
        "derivative_identity_orders_checked":list(range(15)),
        "evidence_limit":"Exact finite algebra and Sturm data, not a Lean proof of the "
                         "infinite shifted-row or smooth-source existence statements, "
                         "the actual Xi cofactor signs, RH or Robin."
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print("Exact source certificate: first defect 13; actual J12 node ratio in (-3/2,-7/5).")
