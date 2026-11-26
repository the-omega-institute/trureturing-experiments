"""Recover the integer port gain from independent geometric and Dirac threshold ratios."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

BOSON = run_path(str(Path(__file__).with_name("ch39-boson-shell-invariant.py")))
GRID, a = BOSON["GRID"], BOSON["a"]
SOURCE = GRID["MATCH"]["SOURCE"]
TARGET = SOURCE["PLANCK"]["TARGET"]
OFFSET = F(1,4)+30+F(2,15)


def calibration(products, epsilon):
    scalar = SOURCE["MATCH"]["SCATTER"]["scalar_thresholds"](products, arb(1), a("1/12"))
    g, v, _ = SOURCE["RATIO"]["MODEL"]["matched"](scalar, epsilon, arb(1))
    return g, v/g*(2*arb.pi()/(3*arb(5).sqrt()/2)).sqrt()


def effective(p):
    ell = 2*p/(1-p)**2
    return (ell**2-a("1/4")).sqrt()


def height(g, anchor, eps, scale, invariant, product):
    power = (scale/anchor)**88*(invariant/g**4)*(eps**6*product)**532
    return (power.root(44)-a(OFFSET)).sqrt()


def verify():
    ctx.prec = 220
    assert 134*F(63,67*7) == 18
    assert 18+F(134,7) == F(260,7)
    assert -4*F(134,7)+88*F(260,7) == 6*532
    products = SOURCE["RATIO"]["product_polynomials"]()
    references = []
    for eps in (a(F(6,19**4)), a(F(1,21700))):
        references.append((*calibration(products, eps), eps))
    # Enclose the logarithmic sensitivity of the boson readout on the whole branch.
    q = a(".16").union(a(".17"))**2
    t = q*(1+a(-GRID["MATCH"]["WINDOW"]).union(a(GRID["MATCH"]["WINDOW"])))/256
    x = 6*q/(q+6+((q+6)**2-12*q).sqrt())
    e = a("1e-6")
    stretch = ((1+e)/(1-e))
    x_box, t_box = [v/stretch.union(1/stretch) for v in (x,t)]
    ax, at = BOSON["derivative_bounds"](x_box,t_box)
    assert abs(ax)+abs(at)+2+abs(ax+at+2) < 1760
    boson_error = 1760*(-(1-e).log())
    product_error = 5*stretch.log()
    rounding = a(SOURCE["PLANCK"]["HALF_UNIT"])/TARGET
    scale_error = -(1-rounding).log()
    log_geometry_error = 2*scale_error+(boson_error+532*product_error)/44
    assert log_geometry_error < a(".000244")
    estimates = []
    zero_data = SOURCE["RATIO"]["high_zero_readings"]()
    for (n, center), (label, gamma, eps_n) in zip(GRID["models"](),zero_data,strict=True):
        assert n == label
        _, vectors, scalars = GRID["mass_ranges"](center,(F(0),F(0)),arb(1))
        invariant = BOSON["reconstruct"](a("2/3"),scalars[5][0],scalars[-3][0],vectors[-1][0])
        product = eps_n**(-6)
        for g, anchor, eps in references:
            value = height(g, anchor, eps, TARGET, invariant, product)
            assert value.overlaps(gamma) and value.rad() < a("1e-6")
        ell = ((gamma/624)**2+a("1/4")).sqrt()
        p = ((2*ell+1).sqrt()-1)/((2*ell+1).sqrt()+1)
        # Work in ratios to avoid any assumption on the sign of individual errors.
        p_noise = p/stretch.union(1/stretch)
        x_noise = effective(p_noise.lower()).union(effective(p_noise.upper()))
        geom_noise = (gamma**2+a(OFFSET))*(log_geometry_error.exp()).union((-log_geometry_error).exp())
        gamma_noise = (geom_noise-a(OFFSET)).sqrt()
        gain = gamma_noise/x_noise
        assert gain > a("623.91") and gain < a("624.09")
        assert abs(624*x_noise-gamma) < a(".14")
        estimates.append(624*x_noise)
        print(n, "gain at 1 ppm and rounded Planck scale", gain.str(7))
    assert all(u < v for u,v in zip(estimates,estimates[1:]))


if __name__ == "__main__":
    verify()
