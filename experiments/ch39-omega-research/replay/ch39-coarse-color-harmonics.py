"""First two color-circle Fourier coefficients and a one-bit small-step construction."""
from fractions import Fraction as F
from json import loads
from math import factorial
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

ROOT = Path(__file__).resolve().parent
geo = run_path(str(ROOT/"ch39-regge-gauge-scale.py"))
q = geo["q"]


def harmonics(height, power, degree, count=2):
    amplitude = 16/(height**2+q("1/4")+q("452/15"))
    b = amplitude/(1-2*amplitude)
    # Positive coefficients of the powers of 1+cos(phi).
    term = q(1)
    sums = [q(1)]+[q(0) for _ in range(count)]
    for k in range(1, degree+1):
        term *= (power-k+1)*b*(2*k-1)/(k*k)
        sums[0] += term
        harmonic = term
        for j in range(1, min(k, count)+1):
            harmonic = harmonic*(k-j+1)/(k+j)
            sums[j] += harmonic
    alpha = 2*power*b
    err = alpha.exp()*alpha**(degree+1)/factorial(degree+1)
    assert err < q("1e-15") and power > degree+1
    integrals = [v+arb(0, err) for v in sums]
    return tuple(v/integrals[0] for v in integrals[1:])


def verify():
    ctx.prec = 200
    mu = 1000*q(F(429,16)).exp()
    power = (geo["BOOK_PLANCK"]/mu)**2/(8*arb.pi())
    eps = F(1, 2**26)
    eta = q(2*eps**2/(1-2*eps**2))
    for index in [3200, *range(11485, 11491)]:
        if index == 3200:
            cell = loads((ROOT/"ch39-tick-separators.json").read_text())["prefix"][3199][0]
            height = q(F(cell,1 << 400)).union(q(F(cell+1,1 << 400)))
            n, degree = q(144196081), 1500
            gap, radius, spread = q("1.09e-7"), q(".0042"), q(".000134")
        else:
            height = acb.zeta_zero(index).imag
            n, degree = power, 128
            gap, radius, spread = q(".0128"), q(".258"), q(".046")
        b1, b2 = harmonics(height, n, degree)
        diff = b2-b1**4
        log_ratio = (b2.log()-4*b1.log())/(-b1.log())
        assert diff > gap and (log_ratio/(2*(1+log_ratio))).sqrt() > radius
        assert (diff-eta)/6 > spread*spread
        assert q("1/2") < b1 < 1 and b1**4+eta < b2 < 2*b1-1
        print("Zero, harmonics and fourth-power discrepancy:", index,
              b1.str(11), b2.str(11), diff.str(9))
    x = eps/3
    d = 2*x*x/(1+x*x)
    assert 2*d < eps*eps and 2**56*d/2 > 1


if __name__ == "__main__":
    verify()
