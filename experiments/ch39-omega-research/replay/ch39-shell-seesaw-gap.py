"""Sharp gap bound for arbitrary clock angles, and minimum common intensity to reach oscillation data."""
from fractions import Fraction as F
from flint import arb

# NuFIT 6.0 normal and inverted, both variants, three-sigma projections (Table 1).
floor=arb("0.0692")/arb("2.606")
assert floor>arb(".02655")
for c2,h,l in ((F(81,98),"0.9689519862424","0.01833231954330"),
               (F(8741547,12563572),"0.8816252200039","0.06272919306043")):
    c=(arb(c2.numerator)/c2.denominator).sqrt()
    assert c*(2+c)>1
    bound=(1-c)**2*(2+c)/4
    assert bound<floor
    # At the first crossing, the adjacent low-mass gap is the smaller gap.
    def strength(high,low):
        return ((1-floor)/high**2+floor/low**2).sqrt()
    opt=strength(1+c,1-c)
    fixed=strength(arb(h)+arb(0,"1e-12"),arb(l)+arb(0,"1e-13"))
    print("axis cosine",c.str(12),"uniform gap supremum",bound.str(12),
          "minimum common intensity for three-sigma lower bound, all angles / at 719",
          opt.str(12),fixed.str(10))
