"""Inverse-coupling bound orthogonal to the charged X vector tower in the bridge model."""
from pathlib import Path
from fractions import Fraction as F
import runpy
from flint import arb

arc = runpy.run_path(str(Path(__file__).with_name("ch39-clock-abelian-pole.py")))
left, right = arc["levels"](False), arc["levels"](True)
def ball(q):
    q = F(q)
    return arb(q.numerator)/q.denominator
# On a multiplet of fixed Y, K=Y-2h/3; two extrema suffice for the lattice distance.
def distance(y):
    n = 3*y//2
    return min(abs(y-F(2*h,3)) for h in (n,n+1))
assert [distance(y) for y in (1,F(9,2),5,4)] == [F(1,3),F(1,6),F(1,3),0]
pure = -F(11,3)*(F(2,3)*4+2+F(4,9)*2)
sector = -F(29,3)*F(2,3)+F(5,3)+F(4,9)*F(23,3)
assert (pure,sector) == (-F(550,27),-F(37,27))
budget = 2*arb.pi()*arb("40.9906")
def crossing(base,coeffs,matching):
    def deficit(x):
        return (ball(base)*x.log()
                + ball(coeffs[0])*arc["exposure"](x,left)
                + ball(coeffs[1])*arc["exposure"](x,right) - ball(matching)*budget)
    lo,hi = arb(1),arb(64)
    assert deficit(hi)>0
    for _ in range(40):
        mid = (lo+hi)/2
        if deficit(mid)>0:
            hi=mid
        else:
            lo=mid
    assert deficit(hi)>0
    return hi.upper()
for count,y in ((4,5),(6,4)):
    coeffs = (F(2,3)*count*distance(1)**2,
              (2*distance(F(9,2))**2+distance(y)**2)/3)
    values = [crossing(base,coeffs,match) for base,match in
              ((pure,F(5,3)),(sector,F(5,3)),(sector,1))]
    print("Y-scalar",y,"coeffs",coeffs,"upper bounds",*[v.str(13) for v in values])
    assert max(values)<32
