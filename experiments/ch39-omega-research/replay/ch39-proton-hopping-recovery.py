"""Exact low- and high-hopping coefficients for global two-site contact recovery."""

from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    source = Path(__file__).with_name("ch39-proton-two-site-exceptional-orders.py")
    u, v, rho, T = sp.symbols("u v rho T")
    moments = {}
    for insertion, raw in run_path(str(source))["insertion_moments"]().items():
        terms = {powers: int(c) for powers, c in raw.to_dict().items()}
        assert all(c > 0 for c in terms.values())
        assert any(r == 0 and k == 0 for _, _, r, k in terms)
        assert all(i + j + 2*k == 16 - sum(insertion) for i, j, _, k in terms)
        poly = sp.Poly.from_dict(terms, (u, v, rho, T)).as_expr()
        moments[insertion] = sp.Poly(poly, T)
    z, a, b = [moments[key] for key in ((0, 0), (1, 0), (0, 1))]
    assert [poly.degree() for poly in (z, a, b)] == [8, 7, 7]
    assert a.eval(0).subs(rho, 0) != 0 and b.eval(0).subs(rho, 0) != 0
    aa_slope = -3 * (24*rho*u*v + 16*u**2 - v**2)
    bb_slope = -3 * (24*rho*u*v - 16*u**2 + v**2)
    p = 20*rho**3 + 45*rho**2 + 18*rho + 1
    h = 15*rho**4 + 80*rho**3 + 90*rho**2 + 24*rho + 1
    static_aa = 3*(15*rho**2 + 12*rho + 1)*h / (4*p**2)
    sign = sp.Poly(sp.cancel(-2*p**3 * sp.diff(static_aa, rho) / (63*rho)), rho)
    assert static_aa.subs(rho, 0) == sp.Rational(3, 4)
    assert static_aa.subs(rho, 1) == sp.Rational(5, 8)
    assert sign.TC() == 1 and all(c > 0 for c in sign.all_coeffs())
    aa, bb, ab = [moments[key] for key in ((2, 0), (0, 2), (1, 1))]

    def at_infinity(numerator, denominator):
        assert numerator.degree() == denominator.degree()
        d = numerator.degree()
        limit = sp.cancel(numerator.LC() / denominator.LC())
        slope = sp.cancel((numerator.nth(d-1) - limit*denominator.nth(d-1))
                          / denominator.LC())
        return limit, slope

    assert [poly.degree() for poly in (aa, bb, ab)] == [6, 6, 6]
    q_limit, q_slope = at_infinity(ab**2, aa*bb)
    assert sp.expand(q_limit - (1+rho)**2) == 0
    total = 16*u**2 + v**2
    expected = (1+3*rho-9*rho**2)*total + 40*rho*u*v
    assert sp.cancel(q_slope/(8*(1+rho)) - expected) == 0
    # Each quotient compares the mixed readout to a diagonal readout.
    quot_a = (ab*a, aa*b)
    quot_b = (ab*b, bb*a)
    slope_a, slope_b = [at_infinity(*pair)[1] for pair in (quot_a, quot_b)]
    assert sp.cancel(slope_a - slope_b - 8*(1+rho)*(16*u**2-v**2)) == 0
    for num, den in (quot_a, quot_b):
        static = sp.cancel(num.eval(0)/den.eval(0))
        assert sp.cancel(num.LC()/den.LC() - (1+rho)) == 0
        assert static.free_symbols <= {rho} and den.degree() == 13

    # At rho=1, v=4u the AA curve as a function of T/u² has a single minimum.
    unit_ray = {u: 1, v: 4, rho: 1}
    num = sp.Poly((z*aa).as_expr().subs(unit_ray), T)
    den = sp.Poly((a*a).as_expr().subs(unit_ray), T)
    slope_num = num.diff()*den - num*den.diff()
    signs = [sp.sign(slope_num.nth(i)) for i in range(slope_num.degree()+1)]
    assert signs == [0] + [-1]*9 + [1]*17

    reference = (sp.Rational(1, 4), sp.Rational(1, 8), sp.Rational(1, 2))
    second = [sp.Rational(s) for s in ("0.1476388875849972", "0.0560800415954564",
                                       "0.7716008713962137")]
    margin = sp.Rational(1, 10**11)
    for key, den in (((1, 1), a*b), ((2, 0), a*a), ((0, 2), b*b)):
        num = z * moments[key]
        assert num.degree() == den.degree() == 14
        lead = sp.cancel(num.LC() / den.LC())
        slope = sp.cancel((num.nth(13) - lead*den.nth(13)) / den.LC())
        if key == (1, 1):
            assert sp.cancel(lead - 3*(rho + 1)/4) == 0
            assert slope.subs(rho, 0).expand() == 3*(16*u**2 + v**2)
            assert sp.cancel((num.eval(0)/den.eval(0)).subs(rho, 0)) == 1
            # All coefficients are nonnegative. Opposite corners bound the ratio on the box.
            n, d = num.eval(16), den.eval(16)
            variables = (u, v, rho)
            def evaluate(poly, point):
                return poly.subs(dict(zip(variables, point, strict=True)))
            first_readout = evaluate(n, reference) / evaluate(d, reference)
            lo = [value - margin for value in second]
            hi = [value + margin for value in second]
            assert 1095 < 1000 * first_readout < 1096
            assert 1222 < 1000 * evaluate(n, lo) / evaluate(d, hi)
            assert 1000 * evaluate(n, hi) / evaluate(d, lo) < 1223
        else:
            assert lead == sp.Rational(3, 4)
            assert sp.cancel(num.eval(0)/den.eval(0) - static_aa) == 0
            assert sp.expand(slope - (aa_slope if key == (2, 0) else bb_slope)) == 0
    print("Moment degrees, normalized quotients, exceptional-ray derivative, static limits and collision bounds verified.")


if __name__ == "__main__":
    verify()
