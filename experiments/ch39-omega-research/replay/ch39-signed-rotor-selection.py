"""A frustrated four-step return potential, its exact tuning interval, and grid localization."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def q(x):
    return arb(x.numerator)/x.denominator


def verify():
    ctx.prec = 192
    source = run_path(str(Path(__file__).with_name("ch39-rotor-shell-returns.py")))
    _, adj, root = source["charged_shell"]()
    _, _, three, four = source["returns"](adj,root)
    b, a, c = three[1], four[1], four[2]
    heights = [z.imag for z in acb.zeta_zeros(10874,3)]
    x = [1-1/(2*(g*g+arb(1)/4)) for g in heights]
    low,high = [(q(a)+2*q(c)*(s+t))/q(b) for s,t in zip(x,x[1:])]
    parameter = Fraction(7460124,259581545)
    parents = (Fraction(583423,20300714), Fraction(6876701,239280831))
    left,right = parents
    assert right.numerator*left.denominator-left.numerator*right.denominator == 1
    assert parameter == Fraction(left.numerator+right.numerator,left.denominator+right.denominator)
    assert q(left) < low < q(parameter) < high < q(right)
    assert arb("7.60031e-18") < high-low < arb("7.60032e-18")

    center = q((parameter*b-a)/(4*c))
    total = 2**32
    center_index = center.acos()*total/(2*arb.pi())
    assert 64361 < center_index < 64362
    def energy(j):
        return 4*q(c)*((2*arb.pi()*j/total).cos()-center)**2
    ground = energy(64362)
    assert ground < energy(64361)
    outside_excess = min(energy(64357),energy(64364))-ground
    assert outside_excess > arb("5.08e-38")
    assert 2*arb("2.54e-40") < arb(".01")*outside_excess
    print("Unique selection interval width", (high-low).str(15), "least-denominator choice", parameter)
    print("Lowest phase-grid pair +/-64362, outside cost above minimum", outside_excess.str(15))


if __name__ == "__main__":
    verify()
