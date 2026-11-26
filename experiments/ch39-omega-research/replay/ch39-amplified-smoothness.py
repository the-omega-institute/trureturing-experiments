"""Conditional SU(5) Haar small balls for the Newton-calibrated and repeated RH plaquette weights."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from json import loads
from flint import arb

ROOT = Path(__file__).resolve().parent
geo = run_path(str(ROOT/"ch39-regge-gauge-scale.py"))
q = geo["q"]


def curvature(height, power):
    s = 1/(height**2+q("1/4"))
    denominator = 1-q("452/15")*s
    return 14*power*(96*s/denominator+(q("160/3")*s/denominator)**2)+q("32/49")


def log_bound_factor(k, radius, epsilon):
    jac = 1-radius**2/6
    return (k*radius**2/2).exp()/jac**20*(q("257/256")*epsilon/radius)**24


def certify(k, radius, median, probability, bits, tail):
    radius, median, probability = q(radius), q(median), q(probability)
    eps = q(F(1, 2**26))
    assert log_bound_factor(k, radius, eps) < q(F(1, 2**bits))
    assert log_bound_factor(k, radius, median) < probability < q("1/2")
    assert 4*probability*(1-probability) < q(tail)**2
    # On these radii 2 sin(257 eps/512) >= eps by the cubic sine bound.
    for e in [eps, median]:
        a = q("257/512")*e
        assert 2*(a-a**3/6) > e
    print("Hessian bound, median threshold, probability bound, small-ball bits:",
          k.str(10), median.str(5), probability.str(5), bits)


def verify():
    mu = 1000*q(F(429, 16)).exp()
    prefactor = (geo["BOOK_PLANCK"]/mu)**2/(8*arb.pi())
    assert prefactor < 30484477
    first = q(11134)
    k_planck = curvature(first, q(30484477))
    certify(k_planck, "27/100", "3/20", ".184", 560, ".78")
    prefix = loads((ROOT/"ch39-tick-separators.json").read_text())["prefix"]
    k_repeat = curvature(q(F(prefix[3199][0], 1 << 400)), q(144196081))
    certify(k_repeat, "1/24", "3/125", ".349", 496, ".954")


if __name__ == "__main__":
    verify()
