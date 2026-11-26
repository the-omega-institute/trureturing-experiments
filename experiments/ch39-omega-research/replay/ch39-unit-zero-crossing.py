"""An actual zero with a prespecified pentagon weight and nearby shell ratio."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import acb, arb

root = Path(__file__).resolve().parent
strip = runpy.run_path(str(root / 'ch39-unit-clock-strip.py'))
unit = strip['unit']
m = unit['model']
flow = runpy.run_path(str(root / 'ch39-spin-phase-closure.py'))
q = m['real']
g = acb.zeta_zero(288).imag
u, v = (m['shell'][s]*g for s in ('a','b'))
width = q('.00001')
terms = list(strip['weighted_terms']())
slope = max(strip['gradient'](terms, u+arb(0,width), v+arb(0,width),
                              q(Fraction(i-4,4)).union(q(Fraction(i-3,4)))).upper()
            for i in range(8))
mesh = runpy.run_path(str(root/'ch39-gaussian-mesh-momentum.py'))
transfer = q(0)
for a,b in [(0,1)] + [(2**j,2**(j+1)) for j in range(10)]:
    for j in range(256):
        left = q(Fraction(a) + Fraction((b-a)*j,256))
        right = q(Fraction(a) + Fraction((b-a)*(j+1),256))
        transfer += (right*right-left*left)*mesh['mesh_slope_error'](left,right)/2
transfer += 2*m['tail_bound'](1024)
print('zero 288 continuum slope upper and mesh correction upper', slope.str(12), transfer.upper().str(12), flush=True)
assert transfer < q('.2192') and slope + transfer < 0
assert q('-.3389') + transfer < q('-.1197')
xs = flow['spectrum'](u,width,True,v,q(0))
r = flow['r']
eps = q('.001')
lo,hi = r*(1-eps),r*(1+eps)
t0,t1 = q('26.994'),q('27.008')
T = t0.union(t1)
for time,scale,k,sign in ((t0,lo.union(hi),0,1),(t1,lo.union(hi),0,-1),(T,lo,1,1),(T,hi,1,-1)):
    margin = sign*flow['run'](time,xs,8192,scale)[k]
    print('zero 288 gauge face margin',margin.lower().str(12),flush=True)
    assert margin>0
