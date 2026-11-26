"""Intertwiner singular gap and rational margins for a rough open Higgs domain."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

import sympy as sp

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))


def intertwiner_gap():
    identity = sp.eye(5)
    rows = []
    for permutation, powers in zip(base["ROWS"], base["POWERS"]):
        link = sp.zeros(5)
        for i, (j, exponent) in enumerate(zip(permutation, powers)):
            link[i, j] = sp.I**exponent
        loop = link**5
        rows.append(sp.kronecker_product(loop, identity)-sp.kronecker_product(identity, loop.T))
    K = sp.Matrix.vstack(*rows)
    polynomial = (K.H*K).charpoly().as_poly().set_domain(sp.ZZ)
    x = polynomial.gens[0]
    positive = polynomial.exquo(sp.Poly(x, x, domain=sp.ZZ))
    assert positive.eval(0) != 0 and positive.count_roots(0, 1) == 0


def margins():
    delta = F(2, 10**12)
    target = F(508983788, 10**9)
    assert (2*40)**2*5 < 180**2
    assert 2*(180+8)+1 == 377
    assert target + 800*delta < F(50898379, 10**8)
    assert F(508983786062, 10**12) + 2*delta < target
    assert 400*delta*32 < F(1, 2000)
    assert 10*delta < 1-10*delta and 40**2*5*delta**2 < 1


if __name__ == "__main__":
    intertwiner_gap()
    margins()
    print("The rooted loop intertwiner has second singular value > 1; open-domain margins passed.")
