"""本原极点的最小分母正方向、单位面积对角线与留出强耦合。"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

model = run_path(str(Path(__file__).with_name("ch39-native-center-parity-crossing.py")))


def first_arc():
    for n in (1, 2, 3):
        radius2 = 3*n*n
        points = [point for point in combinations_with_replacement(range(1, isqrt(radius2)+1), 3)
                  if sum(x*x for x in point) == radius2 and point != (n, n, n)]
        assert points == ([(1, 1, 5)] if n == 3 else [])
    a, b = (3, 3, 3), (1, 1, 5)
    plus = [x+y for x, y in zip(a, b)]
    minus = [x-y for x, y in zip(a, b)]
    assert sum(x*y for x, y in zip(plus, minus)) == 0
    ratio = Fraction(sum(x*x for x in plus), sum(x*x for x in minus))
    assert ratio == 8 and Fraction(sum(x*y for x, y in zip(a, b)), 27) == Fraction(7, 9)
    return arb(ratio.numerator).sqrt()


def readout(levels, radius, u, fermion_weight=3, higgs_levels=None, min_mu="1e13"):
    # For every L >= Lmin, the intercept intervals enclose the finite M_Z subtraction.
    mz, beta = model["mz"], model["beta"]
    lower_log = (arb(min_mu)/mz).log()
    y, w, c = model["intercept"](levels, radius, u, lower_log, color_inverse=arb(0),
                                  fermion_weight=fermion_weight, higgs_levels=higgs_levels)
    slope = (beta[0]-beta[1])/(2*arb.pi())
    log = (y-w)/slope
    inverse_strong = -(115*y-333*w+218*c)/218
    # Bound the derivative in L of the finite reference terms of Y-W.
    scales = [energy/u if charge in (1,3) else energy*u
              for charge, row in levels.items() for energy in row]
    if higgs_levels is not None:
        scales.extend(energy/u for energy in higgs_levels)
    coeff = sum((7*2/(4*arb.pi()),
                 sum(abs(v) for v in model["higgs"][:2])/(2*arb.pi()),
                 fermion_weight*sum(sum(row[:2]) for row in model["weyl"])/(2*arb.pi())))
    assert all(len(row) <= 18 for row in levels.values())
    assert higgs_levels is None or len(higgs_levels) <= 18
    assert 18*coeff*radius**2*(-2*lower_log).exp()/(4*min(x.lower() for x in scales)) < slope
    return inverse_strong, mz*log.exp(), y-beta[0]*log/(2*arb.pi())


def verify():
    ctx.prec = 224
    u = first_arc()
    middle = arb(".1180")
    lower, upper = arb(".1171"), arb(".1189")
    actual, radius = model["setup"]()
    flat = {q: actual[0] for q in actual}
    sites = model["shell"]["soft"]["geo"]["sites"]
    distances = Counter(sum((x-y)**2 for x,y in zip(a,b)) for a,b in combinations(sites, 2))
    assert len(distances) == 30 and sum(distances.values()) == 171
    for name, levels, target in (("实际", actual, 392), ("平坦", flat, 372)):
        selected = []
        for d in sorted(distances):
            inv, scale, common = readout(levels, radius, (arb(3468)/d - 1).sqrt())
            if inv < 1/upper or inv > 1/lower:
                continue
            assert 1/upper < inv < 1/lower and arb("1e13") < scale < arb("1e17") and common > 0
            selected.append(d)
            print(name, "原壳入窗平方弦长、弦数、强耦合", d, distances[d], (1/inv).str(13))
        assert selected == [target]
    flat_inverse, flat_scale, _ = readout(flat, radius, u)
    assert arb(".117285206") < 1/flat_inverse < arb(".117285208")
    assert arb("5.7144e13") < flat_scale < arb("5.7145e13")
    errors = {}
    for index in range(11427, 11436):
        levels, r = (actual, radius) if index == 11431 else model["setup"](index)
        inv, scale, common = readout(levels, r, u)
        strong = 1/inv
        assert lower < strong < upper and arb("1e13") < scale < arb("1e17") and common > 0
        errors[index] = abs(strong-middle)
        if index == 11431:
            assert arb(".11793471703") < strong < arb(".11793471706")
            assert arb("5.95966e13") < scale < arb("5.95968e13")
        print(index, "固定几何刚度的强耦合、交汇能标", strong.str(13), scale.str(12))
    assert all(errors[11432] < err for index, err in errors.items() if index != 11432)


if __name__ == "__main__":
    verify()
