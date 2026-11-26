"""Relate the shell mirror probability to fermion thresholds and test exact cancellation."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Rational, cancel, symbols

COMMON = run_path(str(Path(__file__).with_name("ch39-common-dirac-threshold.py")))
MATCH, a = COMMON["MATCH"], COMMON["a"]


def algebra():
    p, k, u = symbols("p k u", positive=True)
    ell = 2*p/(1-p)**2
    height2 = ell**2-Rational(1,4)
    inverse = (1-p)**2/(2*p)
    assert cancel(16*height2/(4*height2+1)**2-inverse**2*(1-inverse**2/4)) == 0
    positive, negative = [p*u, k*u/p], [u, k*u]
    assert cancel(positive[0]*positive[1]-negative[0]*negative[1]) == 0
    assert cancel(sum(positive)-sum(negative)-u*(1-p)*(k/p-1)) == 0
    b, z = symbols("b z", positive=True)
    # Pointwise maximum of the Euclidean integrand derivative for p=b^2.
    assert cancel((1-b)/(1+b)-z*(1-b*b)/((z+b*b)*(z+1))
                  -(1-b)*(z-b)**2/((1+b)*(z+b*b)*(z+1))) == 0
    for U, V in ((9,10), (10,9)):
        example = Rational(7,U)
        assert example*U == 7 and V/example == Rational(90,7)
        assert 2*example/(1-example)**2 > 1


def verify():
    ctx.prec = 180
    algebra()
    ratio = MATCH["RATIO"]
    polys = ratio["product_polynomials"]()
    for n, kind, r, v, state in MATCH["matched_points"]():
        if kind != "center":
            continue
        gamma = acb.zeta_zero(n).imag
        ell = ((gamma/624)**2+a("1/4")).sqrt()
        root = (2*ell+1).sqrt()
        p = (root-1)/(root+1)
        z = 1/state[0]
        d = r*(a("13/10")*z).sqrt()
        h = d*p.sqrt()
        D = d*acb(((ell+a("1/2"))/2).sqrt(), ((ell-a("1/2"))/2).sqrt())
        E, W, P = D-d, -h*(D+d)/d, -h*D/d
        assert abs(E).overlaps(abs(P)) and abs(W).overlaps(abs(D))
        ratio_ED = abs(E)**2/abs(D)**2
        assert (1-p.sqrt())/(6*(1+p.sqrt())) < a(".014")
        inverse = (1-ratio_ED)**2/(2*ratio_ED)
        predicted = inverse**2*(1-inverse**2/4)
        rho = acb(arb(1)/2, gamma)
        w = (rho-1)/rho
        half_angle = 624*(w-1)/(w+1)
        reflection = (1+half_angle)/(1-half_angle)
        assert predicted.overlaps(reflection.imag**2) and predicted > a(".0031")
        offsets, determinant = COMMON["masses"](d, h, gamma, 624)
        epsilon = (state[5]/state[2]).sqrt()
        scalar = MATCH["SCATTER"]["scalar_thresholds"](polys, r, a("1/12"))
        correction = MATCH["threshold_correction"](r*r, a("1/256"))
        low, new_v, lx = ratio["MODEL"]["matched"](
            [s+c+ds for s,c,ds in zip(scalar, correction, offsets)], epsilon, a("8/5"), r)
        assert new_v.overlaps(v) and low.overlaps(1/(state[0]+3*state[1]).sqrt())
        bridge = lx.exp()
        assert bridge.overlaps(state[2].sqrt()*(-determinant.log()/9).exp())
        assert bridge**2 < 4*arb.pi() and epsilon*bridge*v > 1000
        assert abs(E)*v > 1000 and abs(D)*v > 1000
        print(n, "cancelling xi5", (a("13/10")*p).str(11), "bridge", bridge.str(7),
              "mirror probability", predicted.str(12))


if __name__ == "__main__":
    verify()
