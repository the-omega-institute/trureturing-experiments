"""对实际 724 弦钟的全部有限二值配置作条件计数，复核平直胞的力与复 Hessian。"""
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def face_statistics(n, tests, x, prior):
    root = x.sqrt()
    sum_w, sum_d, sum_d2 = arb(0), arb(0), arb(0)
    sum_score, sum_cross = arb(0), arb(0)
    prior_mass = (1-prior)**n
    odds = prior/(1-prior)
    for k in range(n+1):
        p = arb(k)/n
        mismatch = p**2-x
        likelihood = 1-mismatch**2/2
        weight = prior_mass*(likelihood**tests)
        delta = p-root
        score = tests*mismatch/likelihood
        sum_w += weight
        sum_d += weight*delta
        sum_d2 += weight*delta**2
        sum_score += weight*score
        sum_cross += weight*score*delta
        prior_mass *= odds*(n-k)/(k+1)
    shift = sum_d/sum_w
    variance = sum_d2/sum_w-shift**2
    slope = (sum_cross-shift*sum_score)/sum_w
    return root+shift, slope, n*variance, sum_w


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    ctx.prec = 160
    u = (1+model["phases"](acb.zeta_zero(11366).imag)[724].real)/2
    q = 2*arb(3).sqrt()-3
    b = 4*q/(arb(5)/3).sqrt()
    expected = -54*b/arb(5).sqrt()
    for n in (390,5432,32768):
        tests = n*isqrt(n)
        rows = [(30, 1-q, u, b*b/6, arb(18)/5, -arb(108)/25),
                (20, q, 1-u, b*b/4, -18*arb(3).sqrt()/5, 162*arb(3).sqrt()/25)]
        force, hessian, imaginary = arb(0), arb(0), arb(0)
        for count, area, prior, dx, angle1, angle2 in rows:
            mean, slope, noise, accepted = face_statistics(n,tests,area**2,prior)
            assert 0 < accepted <= 1 and 0 < noise and 0 < slope
            force += count*mean*angle1
            hessian += count*(mean*angle2+2*slope*dx*angle1)
            imaginary += count*noise*angle1**2/(8*arb.pi())
        print("槽位、每面检验数、形状力乘 8π、实 Hessian 与 Regge 之比、虚 Hessian 乘 8π",
              n, tests, force.str(12), (hessian/expected).str(12), imaginary.str(12))
        if n == 32768:
            assert arb(".0004864") < force < arb(".0004865")
            assert arb(".93748") < hessian/expected < arb(".93749")
            assert 0 < imaginary < arb(".3")


if __name__ == "__main__":
    verify()
