"""Coarse-edge influence and center-charged area bounds after independent face refinement."""
from fractions import Fraction as F
from math import factorial


def verify():
    color_lo,weak_lo = F(7180,100000), F(7777,100000)
    color_hi,weak_hi = F(7188,100000), F(7786,100000)
    other = F(448,10000)
    rest = F(9,4)-1-10*(3*color_lo**2+2*weak_lo**2)
    def sup(n):
        return 10*(3*color_hi**n+2*weak_hi**n)+other**(n-2)*rest

    width = sup(3)
    var = 10*(3*color_hi**6+2*weak_hi**6)+other**4*rest
    density_ratio = ((1+width)/(1-width))**13
    assert width < F(6424,100000) and var < F(354,100000)**2
    assert 28**2*density_ratio*var/(1-width)**2 < F(245,1000)**2
    # Sufficient bound on the number of ordered child slots in a rooted face animal.
    e_upper = sum(F(1,factorial(k)) for k in range(9)) + F(10,9*factorial(9))
    for n, target in ((4,F(74,100)),(5,F(43,1000))):
        p = 2*sup(n)/(1+sup(n))
        assert 39*e_upper*p < target
    print("Three-face coarse-edge contraction < 0.245; four/five-face area factors < 0.74 and 0.043.")


if __name__ == "__main__":
    verify()
