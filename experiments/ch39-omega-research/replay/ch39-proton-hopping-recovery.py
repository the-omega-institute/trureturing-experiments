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
    assert sign.TC() == 1 and all(c > 0 for c in sign.all_coeffs())
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
    print("Degrees (14,14), static monotonicity, three high-hopping limits, slopes and collision bounds verified.")


if __name__ == "__main__":
    verify()
