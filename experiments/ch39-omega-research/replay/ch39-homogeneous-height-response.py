"""Uniform-height response from the face incidence counts and conditional mixing."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def pair_error():
    radius = F(1, 2**49)
    kappa = 2**115
    scaled2 = kappa*radius**2
    sixth = 10395*108**6
    bad = sixth/scaled2**6
    density = 4*1440*radius
    off = F(1, 2**30)
    assert scaled2/20 > 23*14+182+14*115+30
    def truncation(mult, degree, moment):
        return mult*(sixth/scaled2**(6-degree)+moment*bad)/(1-bad)
    dx, dy = 27*kappa*radius**3, 2500*kappa*radius**3
    ex = off + dx + 34*density + truncation(F(3,10),1,108)
    ey = off + dy + 1621*density + truncation(15,1,108)
    exy = (off+3242*dx+68*dy+dx*dy+160000*density
           +truncation(F(9,2),2,3*108**2))
    mean_x, mean_y = F(54,10), F(1092,5)
    largest_cov = F(3294,25)
    u = F(32,11134**2)
    denom_loss = (1/(1-u)**2-1)*(mean_x*mean_y+largest_cov+exy+(mean_x+ex)*(mean_y+ey))
    assert exy+mean_x*ey+mean_y*ex+ex*ey+denom_loss < F(1,100)
    assert 2*(F(1,2**60)+6*4*672*radius) < F(1,2**33)
    assert 2*(F(1,2**60)+256*4*672*radius) < F(1,2**28)
    assert 2*(F(1,2**60)+4*224*radius) < F(1,2**36)
    # In the first oscillation estimate the actual off-chart cost is much smaller.
    assert scaled2/20 > 23*14+182+14*115+60


def geometry():
    faces = [(a,b) for a,b in product(range(1,16),repeat=2) if not a&b]
    edge_counts = Counter(edge for a,b in faces for edge in (a,b,a|b))
    types = Counter(6+sum(edge_counts[e] for e in (a,b,a|b)) for a,b in faces)
    assert types == {36:30,42:20}
    assert 2*28 < 2**36
    disjoint_error = F(84*14*256,2**33) + F(2*84*392,2**(33+28))
    assert disjoint_error < F(1,1000)
    assert F(34,100)+disjoint_error < F(7,20)
    return sorted(types)


def gaps(types):
    ctx.prec = 160
    setup = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
    q = setup["q"]
    power = (setup["high"]["geo"]["BOOK_PLANCK"]/(1000*q(setup["L0"]).exp()))**2/(8*arb.pi())
    strengths = [1/(acb.zeta_zero(i).imag**2+q(1)/4) for i in range(11485,11491)]
    for n in types:
        coeff = q(366)*n/25
        intervals = [(power*(a-b)*(coeff-q("0.35")),power*(a-b)*(coeff+q("0.35")))
                     for a,b in zip(strengths,strengths[1:])]
        floor = min(x[0] for x in intervals)
        threshold = {36: "0.01075", 42: "0.01254"}[n]
        assert floor > q(threshold)
        print(n, "five lower bounds:", " ".join(x[0].str(9) for x in intervals))


if __name__ == "__main__":
    pair_error()
    gaps(geometry())
