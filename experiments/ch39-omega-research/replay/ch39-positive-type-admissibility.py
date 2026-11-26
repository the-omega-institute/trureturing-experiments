"""Positive-type curvature and a compact semialgebraic support for the chiral smoothness bound."""
from fractions import Fraction as F
from itertools import combinations, permutations
from json import loads
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EPS = F(1, 2**26)
RADIUS = EPS/4
F0 = 30+F(2, 15)


def cayley(x):
    return (1-x*x)/(1+x*x), 2*x/(1+x*x)


def verify():
    x = EPS/2+EPS**3/8
    re, im = cayley(x)
    assert 2*(1-re) > EPS**2
    assert 1/(8*x*x) > 2**51-F(1, 4)
    prefix = loads((ROOT/"ch39-tick-separators.json").read_text())["prefix"]
    for n, height, ceiling, rotation, lower in [
        (30484477, F(11134), 1, F(1, 3), F(1, 9)),
        (144196081, F(prefix[3199][0], 2**400), 42, F(1, 20), F(4, 25)),
    ]:
        stiffness = 4*n/(height**2+F(1, 4)+F0)
        assert stiffness < ceiling
        assert 1-8*ceiling*rotation**2 >= lower
    re, im = cayley(RADIUS/8)
    third = re**3-3*re*im**2
    assert 4*(1-re) < RADIUS**2/4 and 4*(1-third) > RADIUS**2
    assert 1-4*42*(1-re) > F(8, 9)
    assert 2*(2*RADIUS) == EPS and 2*(3*RADIUS) < F(1, 1024)
    size = 250000
    assert 2*(6*8)**2*4 < RADIUS**2*size**4
    # At most six adjacent swaps put the two ordered blocks in their canonical order.
    assert max(sum(a>b for a,b in combinations(order, 2)) for order in permutations(range(4))) == 6
    print("Hard-support curvature exceeds 2^51 - 1/4; the semialgebraic cutoff violates a four-point Gram test.")


if __name__ == "__main__":
    verify()
