"""Pure Gaussian ground-state selection at an actual zero and on a nearby phase strip."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import acb, arb

root = Path(__file__).resolve().parent
strip = runpy.run_path(str(root/'ch39-unit-clock-strip.py'))
mesh = runpy.run_path(str(root/'ch39-gaussian-mesh-momentum.py'))
m, q = strip['m'], strip['q']
terms = list(strip['weighted_terms']())
cutoff = 8192
error = q(0)
for a,b in [(0,1)] + [(2**j,2**(j+1)) for j in range(13)]:
    for j in range(256):
        left = q(Fraction(a) + Fraction((b-a)*j,256))
        right = q(Fraction(a) + Fraction((b-a)*(j+1),256))
        error += (right*right-left*left)*mesh['mesh_slope_error'](left,right)/2
tail = 2*m['tail_bound'](cutoff)
assert error < q('.05938') and tail < q('.01316')


def upper_slopes(u, v):
    continuous = max(strip['gradient'](terms, u, v,
                     q(Fraction(i-4,4)).union(q(Fraction(i-3,4))), 0).upper()
                     for i in range(8))
    # Each moving 804 edge differentiates its axis through sin(v).
    discrete = continuous + abs(v.sin())*error + tail
    return continuous, discrete.upper()


g = acb.zeta_zero(288).imag
width = arb(0,q('.00001'))
u0, v0 = (m['shell'][s]*g + width for s in ('a','b'))
values = upper_slopes(u0, v0)
print('zero 288 pure Gaussian slope upper, continuum and mesh', *(v.str(12) for v in values), flush=True)
assert values[0] < q('-.06037') and values[1] < q('-.01474')
u = q('2.8') + arb(0,q('.0001'))
worst = q(-10)
for j in range(32):
    v = q('.40') + q(Fraction(j,100)).union(q(Fraction(j+1,100)))
    c, d = upper_slopes(u, v)
    worst = max(worst,d)
print('strip pure finite-mesh slope upper',worst.str(12),flush=True)
assert worst < q('-.00113')

flow = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
for phase, sign in ((q('.40'),1),(q('.72'),-1)):
    xs = flow['spectrum'](q('2.8'), q('.0001'), True, phase, q(0))
    time = q('26.97').union(q('27.03'))
    for j in range(3):
        target, res, _, _ = flow['run'](time,xs,8192 if j<2 else 16384,return_target=True)
        assert time.overlaps(target)
        time = time.intersection(target)
    print('pure strip gauge phase margin', (sign*res).lower().str(12), flush=True)
    assert sign*res > 0
