"""Finite Bloch spectra and mean scalar effective resistance on the five-periodic box."""
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))


def scalar_resistance():
    # Pair conjugate Fourier values. The square-root terms cancel under i <-> j.
    total = Fraction(0)
    for i in range(5):
        for j in range(5-i):
            if i+j:
                multiplicity = comb(4, i)*comb(4-i, j)*2**(i+j)
                total += multiplicity*Fraction(10*(i+j), 25*(i+j)**2-5*(j-i)**2)
    assert 2*total/625 == Fraction(586384, 1959375)


def laplacian_spectrum():
    ctx.prec = 160
    links = base["field"]()[:4]
    c = (arb(5).sqrt()-1)/4
    root = acb(c, (1-c*c).sqrt())
    phases = [[u*root**j for j in range(5)] for u in links]
    # Four rational cuts, all measured in units of 10^-8.
    endpoints = (182679944, 182679945, 208329614, 208329615)
    cuts = [arb(k)/100000000 for k in endpoints]
    counts = [[] for _ in cuts]
    origin = 8*base["identity"](5)
    for p in product(range(5), repeat=4):
        forward = sum((phases[a][p[a]] for a in range(1, 4)), phases[0][p[0]])
        block = origin-forward-base["adjoint"](forward)
        coeffs = block.charpoly().coeffs()
        assert all(z.imag.contains(0) for z in coeffs)
        poly = arb_poly([z.real for z in coeffs])
        for points, cut in zip(counts, cuts):
            below = 5-base["positive_roots"](poly, cut)
            points.extend([p]*below)
    first, second = (4, 1, 2, 1), (4, 0, 3, 2)
    assert counts == [[], [first], [first], sorted([first, second])]
    lower = Fraction(endpoints[0], 100000000)
    assert lower > 4*Fraction(67579, 100000)**2
    assert lower-Fraction(8, 100000) > 4*Fraction(67578, 100000)**2
    print("Lowest two Bloch momenta and bounds in 10^-8 units:", first, second, endpoints)


if __name__ == "__main__":
    scalar_resistance()
    laplacian_spectrum()
