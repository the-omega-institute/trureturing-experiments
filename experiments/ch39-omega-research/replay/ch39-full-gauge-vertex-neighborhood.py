"""Uniform determinant continuity budget transverse to pure weak gauge links."""
from fractions import Fraction as F


def verify():
    # Relative support radius in units of the smaller of 1 and the core gap.
    width = F(1, 2**9600)
    assert 96*width < F(1,2)
    sign_cost = 192*width
    rotation_cost = sign_cost
    row_cost = sign_cost+2*rotation_cost
    assert row_cost < 2**10*width
    n = 1440+12
    terms = 16*(6*2**4)**3
    absolute_error = terms*n*2**(n-1)*row_cost
    assert absolute_error < 2**1498*width
    r = F(1,2**20)
    numerator = F(216,16**5)*(r/4)**147 * F(1,10)**1440
    denominator = (r/4)**147
    assert absolute_error < min(numerator,denominator)/4
    assert numerator/(4*2**1440) > F(1,2**9474)
    assert 32*r/(1-32*r) < F(1,100)
    print("Full H-link integral on positive transverse fibers: |A| > 2^-9474.")
    print("Common upper bound on the change of the numerator integrand < 2^-8102.")


if __name__ == "__main__":
    verify()
