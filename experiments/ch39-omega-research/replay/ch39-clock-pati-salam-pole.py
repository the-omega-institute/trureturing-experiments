"""Lower bounds from a color- and internal-clock-singlet elementary condensate."""
from fractions import Fraction as F
from math import comb
from flint import arb

# h=6n, D=-3k, R=k-8n for n>=1 (conjugate all weights if necessary).
# A complex X representation is forced by h != 0. For k!=0 bound spectator dimension.
small = [comb(k+3,3)*(9-k) for k in range(1,4)]
assert min(small) == 32 and comb(7,3) > 32
tx, dx = F(15,2), 10  # SU(3) (3,0) or (0,3); dimension and index monotone in labels
def spin_index(label):
    return F(label*(label+1)*(label+2),12)
assert spin_index(8) == 60
weight_x = F(16,3)
bx0 = -11 + tx*9/3
br0 = -F(22,3) + dx*spin_index(8)/3 + 4  # three ordinary right families
bx1 = -11 + tx*32/3
assert (bx0,br0,bx1) == (F(23,2),F(590,3),69)
def pole_upper(budget):
    # Tree-level matching: (16/3)/gX^2+1/gR^2=1/gU^2.
    return (2*arb.pi()*arb("40.9906") * budget.denominator / budget.numerator).exp()
bounds = [pole_upper(weight_x*bx0+br0),pole_upper(weight_x*bx1)]
assert bounds[0] < arb("2.72") and bounds[1] < arb("2.02")
print("upper bounds for the one-loop pole ratio, k=0 and k!=0",
      [value.upper().str(13) for value in bounds])
