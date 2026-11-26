"""A CKM obstruction for round-sphere quadratic zero modes sampled at three cyclic points."""
from fractions import Fraction as F
from flint import arb

def rational(x):
    value=F(x)
    return arb(value.numerator)/value.denominator
# G=6I-V^*V after rephasing the degree-two Lagrange sections and a common scaling.
# The conditioned two-by-two minor of G has diagonal >=5-1/5, off-diagonal <=1+1/5.
mu = (1+F(1,5))/(5-F(1,5))
assert mu==F(1,4)
low = 1-2*mu
assert low>0
# Enlarge the quoted central values by the listed input errors at 1 TeV.
up = rational('.00316')/rational('.8573')
down = rational('.000271')/rational('.01392')
def direction_error(r):
    x=rational(mu)
    return (2*(1+x)).sqrt()*x*r/(rational(low)-(1+x)*r)
upper=direction_error(up)+direction_error(down)
different=rational(mu)+arb(2).sqrt()*upper
print('aligned-site off-diagonal upper, distinct-site diagonal upper',upper.upper().str(12),different.upper().str(12))
assert upper<rational('.020') and different<rational('.28')
