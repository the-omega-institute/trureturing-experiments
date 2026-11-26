"""Finite soft-pinning comparison for parallel fundamental Higgs links."""
from fractions import Fraction as F

from flint import arb


def verify():
    t = F(1,15)
    strength = F(1,11134**2)
    power = 30484477
    d = (30+2*t)*strength
    derivative = power*strength*50*(1+t)/(1-d)
    assert derivative**2*27 < 70**2
    lower = (1-arb(d.numerator)/d.denominator)/(1+arb(d.numerator)/d.denominator)
    assert lower**power > arb(2)**-23
    a = F(1,2**52)
    edges = faces = length = 2**20
    degree = 14
    kappa = 2**115
    exponent = F(3,4)*kappa*a*a
    # (5/a)^10 < 2^(10*(52+3)), and exp(-x) <= 2^(-x).
    assert exponent > 20+23*degree+10*(52+3)+30
    eta = 70*faces*a
    assert 2*eta < F(1,2)
    error = 2*length*a + 4*eta + F(2,2**30)
    assert error < F(1,10**7)
    print("Soft-pinning comparison below 1e-7 for the stated finite volume, perimeter and coefficient.")


if __name__ == "__main__":
    verify()
