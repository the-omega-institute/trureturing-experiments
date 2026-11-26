"""Tagged zero-height response of a rank-four pinned fundamental triangle."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def moment_budget():
    kappa = 2**115
    radius = F(1, 2**49)
    R2 = kappa*radius**2
    sixth = 10395 * 54**6
    bad = sixth/R2**6
    assert R2/20 > 23*14+180+14*115+30
    off = F(1, 2**30)
    density = 4*720*radius
    dx = 27*kappa*radius**3
    dy = 2500*kappa*radius**3
    def truncation(multiplier, power, moment):
        return multiplier*(sixth/R2**(6-power)+moment*bad)/(1-bad)
    ex = off + dx + 17*density + truncation(F(3, 10), 1, 54)
    ey = off + dy + 811*density + truncation(15, 1, 54)
    exy = (off + 1622*dx + 34*dy + dx*dy + 40000*density
           + truncation(F(9, 2), 2, 3*54**2))
    mean_x = F(54, 10)
    mean_y = 3*(72+12*F(1,15))
    moment_error = exy + mean_x*ey + mean_y*ex + ex*ey
    # The tagged-height derivative uses Y/[(1+s*f0)(1+s*f)].
    u = F(32, 11134**2)
    distortion = 1/(1-u)**2-1
    target = 9*(4*36+8*F(1,15)*F(9,2))/10
    assert target == F(3294,25)
    absolute_xy = mean_x*mean_y + target + exy
    # Both signed errors also bound the corresponding absolute moments,
    # since the reference quadratic forms are nonnegative.
    ratio_error = distortion*(absolute_xy+(mean_x+ex)*(mean_y+ey))
    assert moment_error+ratio_error < F(1,100)
    print("Scaled response coefficient", target, "error budget", float(moment_error+ratio_error))
    return target


def height_gaps(target):
    ctx.prec = 160
    settings = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
    q = settings["q"]
    power = (settings["high"]["geo"]["BOOK_PLANCK"]/(1000*q(settings["L0"]).exp()))**2/(8*arb.pi())
    strengths = [1/(acb.zeta_zero(i).imag**2+q(1)/4) for i in range(11485,11491)]
    factor = arb(target.numerator)/target.denominator
    spans = []
    for first, second in zip(strengths, strengths[1:]):
        coefficient = power*(first-second)
        lower, upper = coefficient*(factor-q("0.01")), coefficient*(factor+q("0.01"))
        assert lower > 0
        spans.append((lower, upper))
        print(lower.str(12), upper.str(12))
    intervals = [(410653,410716), (626777,626873), (269641,269683),
                 (569463,569551), (268976,269018)]
    for (lower, upper), (lo, hi) in zip(spans, intervals, strict=True):
        assert lower > q(lo)/10**8 and upper < q(hi)/10**8
    return spans


if __name__ == "__main__":
    height_gaps(moment_budget())
