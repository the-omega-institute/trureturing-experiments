# Check of the cubic C([0,6],C) counterexample to Krishna's C*-algebraic Sendov conjectures
# (arXiv:2203.06916v1, Conjectures 2.4 and 2.5): branch exchange of the critical points along a loop.
from fractions import Fraction as F
import cmath, math
V = [complex(-0.9,0), complex(-0.1,0.8), complex(0.1,0.8), complex(0.1,0.9), complex(-0.1,0.9), complex(-0.1,0.8), complex(-0.9,0)]
def c(t):
    k = min(int(math.floor(t)), 5); u = t - k
    return (1-u)*V[k] + u*V[k+1]
a = [0.5, -0.5]
ok = True
# exact vertex moduli
Vq = [(F(-9,10),F(0)),(F(-1,10),F(4,5)),(F(1,10),F(4,5)),(F(1,10),F(9,10)),(F(-1,10),F(9,10)),(F(-1,10),F(4,5)),(F(-9,10),F(0))]
assert max(x*x+y*y for x,y in Vq) == F(82,100)
# discriminant zero would need c = +-i sqrt(3)/2; imaginary-axis crossings at 4/5 i and 9/10 i straddle sqrt(3)/2
assert F(4,5)**2 < F(3,4) < F(9,10)**2
N = 60000
prev = None; worst_dist = {}
maxmod = 0; min_disc = 1e9; min_gap = 1e9; min_aj_gap = 1e9
for i in range(N+1):
    t = 6*i/N
    cc = c(t); D = cc*cc + 0.75
    min_disc = min(min_disc, abs(D))
    r = cmath.sqrt(D)
    roots = [(cc + r)/3, (cc - r)/3]
    if prev is None:
        bp, bm = (roots[0], roots[1]) if roots[0].real > roots[1].real else (roots[1], roots[0])
    else:
        # continue branches by nearest root
        if abs(roots[0]-prev[0]) + abs(roots[1]-prev[1]) <= abs(roots[1]-prev[0]) + abs(roots[0]-prev[1]):
            bp, bm = roots[0], roots[1]
        else:
            bp, bm = roots[1], roots[0]
    prev = (bp, bm)
    if i == 0: start = (bp, bm)
    for b in (bp, bm):
        assert abs(3*b*b - 2*cc*b - 0.25) < 1e-9
        maxmod = max(maxmod, abs(b))
        # barycentric positive weights
        q = [1/abs(b-aj)**2 for aj in (0.5, -0.5, cc)]
        z = sum(qj*aj for qj, aj in zip(q, (0.5, -0.5, cc)))/sum(q)
        assert abs(z - b) < 1e-9
    min_gap = min(min_gap, abs(bp-bm))
    min_aj_gap = min(min_aj_gap, abs(cc-0.5), abs(cc+0.5))
    maxmod = max(maxmod, abs(cc))
end = prev
print('start b+ b-:', start)
print('end   b+ b-:', end)
assert abs(end[0]-start[1]) < 1e-6 and abs(end[1]-start[0]) < 1e-6   # monodromy: the branches exchange
d0 = abs(start[1]-0.5); d6 = abs(end[0]-0.5)
print('|b-(0)-1/2| =', d0, ' |b+(6)-1/2| =', d6, ' (12+sqrt39)/15 =', (12+math.sqrt(39))/15)
assert d0 > 1 and d6 > 1
print('max |a_j|,|b| =', maxmod, '< 1;  min |D| =', min_disc, ' min |b+ - b-| =', min_gap, ' min |c -+ 1/2| =', min_aj_gap)
assert maxmod < 1 and min_disc > 0.05 and min_gap > 0.05 and min_aj_gap > 0.05
print('ALL_OK')
