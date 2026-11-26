"""Exact light branching, residual baryon symmetry and charge-lattice checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd

# sign distinguishes (4,2,1,3_F) and its right-handed partner (bar4,1,2,bar3_F).
spectrum = Counter()
for sign, color_index, family in product((1, -1), range(4), range(3)):
    colored = color_index < 3
    d = sign * (1 if colored else -3)
    phase = (7 * d + 3 * sign) % 9
    assert phase == (sign if colored else 0) % 9
    for component in (-1, 1):
        r = component if sign == -1 else 0
        spectrum[sign, colored, d + 3 * r] += 1
assert spectrum == {(1, True, 1): 18, (1, False, -3): 6,
                    (-1, True, -4): 9, (-1, True, 2): 9,
                    (-1, False, 0): 3, (-1, False, 6): 3}
assert 4 * (1 - (-2)) == 12
# Cubic indices of fundamental minus antifundamental (multiplicities at left and right).
assert 2 * 3 - 2 * 3 == 0 and 4 * 2 - 4 * 2 == 0
assert 4 * 3 % 2 == 0  # both SU(2) groups
assert (3 * 7) % 9 == 3 and (3 * 7 * -3) % 9 == 0  # g^3 is the color center

# Columns H,D,R,S,F. The last entry is the integer triality of the family representation.
breakers = ((0,-3,1,0,1), (2,0,0,-4,0), (0,0,-2,3,0))
for h,d,r,s,f in breakers:
    assert 4*h + d + 3*r + 2*s == 0
    assert (7*d + 3*f) % 9 == 0
# After modding out connected color, beta=3 theta_D is solved by the first row.
# Independent of family component, the remaining two rows act on H,R,S.
a, b = (2,0,-4), (0,-2,3)
minors = [a[i]*b[j]-a[j]*b[i] for i,j in combinations(range(3),2)]
assert gcd(*minors) == 2 and gcd(*a,*b) == 1
assert [abs(n)//2 for n in reversed(minors)] == [4,3,2]
assert F(7,9)*(-3) + F(1,3) == -2
print("SM Weyl components after branching", sorted(spectrum.items()))
print("residual components modulo connected color and hypercharge: Z/2 x Z/3")
