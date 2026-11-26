"""One-loop crossing with proper-time heat thresholds and arbitrary complementary scalar towers."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-slope-matched-one-loop.py'))
q = m['q']
tau, required = m['g']['base']
rlo,rhi = q(127),q(450)
def Jinf(y):
    return y.expint(1)/2
u = q(36)
def multiplier(u):
    # For y>=1 integrate from v=1 to v=2 in J = (1/2) int exp(-vy) dv/v.
    high = 2*(2-u).exp()/q(2).log()
    # On y >= exp(-8) the omitted upper-limit exponential integral is tiny.
    assert (2*tau-8).exp()>100
    tail = (-q(100)).exp()/200
    tiny = 1/(Jinf((-q(8)).exp())-tail)
    candidates = [high.upper(),tiny.upper()]
    for j in range(4000):
        t = (q(j)/1000).union(q(j+1)/1000)
        y = (-2*t).exp()
        lower = Jinf(y.upper())-tail
        candidates.append(((-u*y.lower()).exp()/lower).upper())
    return max(candidates)
def exposure_upper(r):
    low = sum(Jinf((x/r**2).lower()) for x in m['xs'])
    high = 20*Jinf(q('1.2710')*(400/r)**2)
    for n,factor in m['factors']:
        c = factor/r**2
        for j in range(1,6):
            high += n*Jinf(c*j*j/(1+arb.pi()**2*j*j*m['h']**2/6))
        # j>=6: y>= c*j^2/3 >= 2*c*j >1/2, hence E1(y)/2 <= exp(-y).
        a = (-2*c).exp()
        assert 12*c > q('1/2')
        high += n*a**6/(1-a)
    return low+high
price = multiplier(u)
ratio = m['theta'](804)/m['theta'](12)
arc_u = ratio**2
arc_price = multiplier(arc_u)
arc_gaps = []
rows=[]
for lo in range(127,450,10):
    hi = min(450,lo+10)
    budget = exposure_upper(q(hi))-required
    left = sum((-u*x/q(lo)**2).exp() for x in m['xs'])
    gap = left-price*budget
    arc_left = sum((-arc_u*x/q(lo)**2).exp() for x in m['xs'])
    arc_gaps.append(arc_left-arc_price*budget)
    rows.append(gap)
    print('r',lo,hi,'J_R upper',budget.upper().str(9),'heat difference lower',gap.lower().str(9))
print('universal price and gap lower',price.str(12),min(v.lower() for v in rows).str(12))
assert min(v.lower() for v in rows) > 0

# Construct a common crossing and two matching heat moments at Q with 30 charged left levels.
band = q('1.2710').union(q('1.273'))
low = [x/400**2 for x in m['xs']]
heat = sum((-y).exp() for y in low)+12*(-band).exp()
moment = sum(y*(-y).exp() for y in low)+12*band*(-band).exp()
assert min(y.lower() for y in low)*(2*tau).exp()>100
tail_error = (-q(100)).exp()/200
left_J = sum(Jinf(y.lower()).union(Jinf(y.upper())) for y in low)+12*Jinf(band.lower()).union(Jinf(band.upper()))
left_J += q(0).union(-30*tail_error)
def residual(s):
    y = (-2*s).exp()
    weight = 18*(-y).exp()
    w, mu = heat-weight, moment-y*weight
    assert w>0 and w<12*(-q(1)).exp()
    assert mu-w > (-q(2)).exp() and w*(12/w).log()>mu
    low_R_J = 18*(Jinf(y.lower()).union(Jinf(y.upper()))+q(0).union(-tail_error))
    assert y*(2*tau).exp()>100
    return left_J-low_R_J-required-q(0).union(w.upper()/2)
for j in range(1980,2150):
    residual((q(j)/1000).union(q(j+1)/1000))
left_residual, right_residual = residual(q('1.98')),residual(q('2.15'))
assert left_residual>0 and right_residual<0
# All three complete-family towers have y_j >= j/e^2. Integral test for a decreasing E1:
# sum_{j>=1} (1/2)E1(j/e^2) <= e^2/2 * integral_0^infty E1(x) dx = e^2/2.
fermion_J = 3*q(2).exp()/2
inverse_lower = 4*arb.pi()*m['g']['g']['G0'][0]-(q('41/10')*tau+
    q('2/15')*left_J.upper()+q('8/15')*(left_J.upper()-required)+q('8/3')*fermion_J)/(2*arb.pi())
assert inverse_lower>1
print('proper-time crossing residual at endpoints',left_residual.lower().str(12),right_residual.upper().str(12))
print('inverse unified coupling lower',inverse_lower.lower().str(12))

print('proper time at the extreme arc ratio, gap lower',min(v.lower() for v in arc_gaps).str(12))
assert min(v.lower() for v in arc_gaps)>3
