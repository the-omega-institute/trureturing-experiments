"""Isospectral vectorlike doublets and complex singlets with an all-scale internal trace identity."""
from fractions import Fraction
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-slope-matched-one-loop.py'))
q = m['q']
g = m['g']['g']
tau, net = m['g']['base']
xs = m['xs']
L = list(map(Fraction,('2/5',2,1)))
R = list(map(Fraction,('8/5',0,1)))
states = [q(a*5/3+b) for a,b in zip(L,R)]
assert all(3*a+3*b == 6 for a,b in zip(L,R))
assert all(a*5/3+b == 2+2*a/3 for a,b in zip(L,R))
target = net/2
r_step = ((target+sum(x.log() for x in xs)/2)/18).exp()
assert all(x < r_step**2 for x in xs) and m['edge19']>r_step**2
assert r_step/min(x.sqrt().lower() for x in xs)<tau.exp()
# Smooth exposure J is strictly increasing in r for each positive mass.
def smooth_interval(r):
    ys = [x/q(r)**2 for x in xs]
    tail = (-q(100)).exp()/200
    assert min(y.lower() for y in ys)*(2*tau).exp()>100
    result = sum(y.upper().expint(1).union(y.lower().expint(1))/2 for y in ys)
    # Use 38 boundary conditions and the lower Dirichlet floor c*j^2/3.
    # All omitted y>1, so E1(y)/2 < exp(-y).
    first = m['edge19']/q(r)**2
    c = m['d']/(3*q(r)**2)
    assert first>1 and c>1
    tail_high = 20*(-first).exp()+270*(-c).exp()/(1-(-c).exp())
    return result+(-18*tail).union(tail_high)
def bracket():
    # Independently bisect a lower and upper endpoint for all allowed geometries.
    ends=[]
    for direction in ('lower','upper'):
        lo,hi=Fraction(6),Fraction(22)
        for _ in range(18):
            mid=(lo+hi)/2
            val=smooth_interval(mid)
            if (val.upper() if direction=='lower' else val.lower()) < target:
                lo=mid
            else:
                hi=mid
        ends.append(lo if direction=='lower' else hi)
    assert smooth_interval(ends[0]) < target and smooth_interval(ends[1])>target
    return q(ends[0]).union(q(ends[1]))
r_smooth = bracket()
# For either choice each full family contributes a common vector.
# Sharp scheme: three lines, first mass >=Q/e, each reaches at most seven modes.
# Smooth scheme: on each line y_j>=j/e^2 and sum_i J <= e^2/2.
core = 4*arb.pi()*g['G0'][0]-(q('41/10')*tau+states[0]*target)/(2*arb.pi())
for label,r,family_budget in [('step',r_step,q(21)),('smooth',r_smooth,3*q(2).exp()/2)]:
    alpha = core-q('8/3')*q(0).union(family_budget)/(2*arb.pi())
    assert alpha>1
    print(label,'scale ratio',r.lower().str(13),r.upper().str(13),'inverse coupling',alpha.lower().str(12),alpha.upper().str(12))
