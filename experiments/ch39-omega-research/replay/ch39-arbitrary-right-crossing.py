"""A priori enclosures at r=400 allowing arbitrary positive complementary scalar spectra."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root/'ch39-right-scalar-rge.py'))
g = model['g']
m = model['data']['g']
q = g['q']
D = m['D']
cap = 4*arb.pi()
dot = m['dot']
# Retain the previously certified energy and common inverse coupling window.
tau = q('26.92').union(q('26.97'))
alpha_inverse = q('33.1367').union(q('33.2625'))
# Gauge-square lower bounds follow from the three inverse-coupling equations,
# discarding positive heavy contributions. The SM color self coefficient is -26.
def gauge_lower(t):
    lower = []
    for i in range(3):
        adverse_gauge = sum(min(q(0), b) for b in g['B'][i])*cap
        beta_floor = g['beta'][i]+(adverse_gauge-3*cap*sum(g['D'][i]))/D
        upper_G = g['G0'][i]-2*t*beta_floor/D
        assert upper_G > 0
        lower.append(1/upper_G.upper())
    return lower
# On each step a log-squared-Yukawa supersolution ignores the negative
# opposite-sector term and keeps the sample's actual gauge lower bounds.
y_upper = [v.upper() for v in g['y0']]
integral = arb(0)
step = tau.upper()/256
for j in range(256):
    t = (step*j).union(step*(j+1))
    gs = gauge_lower(t)
    ceiling = [v+q('.01') for v in y_upper]
    yukawa = [v.exp() for v in ceiling]
    trace = [sum(yukawa[i:i+3]) for i in (0,3,6)]
    total = 3*(trace[0]+trace[1])+trace[2]
    for i in range(9):
        family = i//3
        rate = (3*yukawa[i]+2*(total-dot(g['C'][family],gs)))/D
        # Uniform upper derivative up to the first possible exit of this tube.
        assert step*rate < q('.01')
        y_upper[i] = (y_upper[i]+step*rate).upper()
    integral += step*dot(g['D'][0],trace)
# Hypercharge has nonnegative gauge coefficients. Lightweight Yukawas can
# only increase this upper bound by the amount just enclosed.
WR = ((g['G0'][0]-alpha_inverse/(4*arb.pi()))*D/2-g['beta'][0]*tau
      -g['d'][0]*model['WL']-q('8/3')+integral/D)/q('8/15')
print('arbitrary-right integrated hypercharge Yukawa upper',integral.upper().str(12))
print('arbitrary-right clipped total exposure upper',WR.upper().str(12))
assert WR < q(53)

read = runpy.run_path(str(root/'ch39-second-cutoff.py'))
def trace_price(u, tau_floor):
    u = q(u)
    # t = -log y / 2; for y>=1 the numerator is nonpositive.
    # On [0,1/4] bound the difference quotient by a derivative bound.
    ymin = (-q('1/2')).exp()
    assert u*ymin > 1
    first = 2*u*ymin*(-u*ymin).exp()
    a = (1-u).exp()
    bounds = [first.upper(), q('1/4').upper(), (1/tau_floor).upper()]
    for j in range(250,4000):
        t = (q(j)/1000).union(q(j+1)/1000)
        y = (-2*t).exp()
        bounds.append((((-u*y).exp()-a*(-y).exp())/t).upper())
    return max(bounds)
def required_exposure(u, tau_floor):
    price = trace_price(u,tau_floor)
    total = read['low_left']+read['high']
    low = sum((-u*y).exp() for y in read['low'])
    numerator = low-(1-q(u)).exp()*total
    return price,numerator/price
# Even states below the input obey 1 <= price*tau; arbitrary multiplicities are allowed.
one_time, net = model['data']['g']['base']
one_W = model['WL']-net
for u, tau_floor, budget in ((8,one_time.lower(),one_W),(24,tau.lower(),WR)):
    price, needed = required_exposure(u,tau_floor)
    gap = (needed-budget)*price
    print('inverse cutoff',u,'right logs needed above',needed.lower().str(12),
          'available below',budget.upper().str(12),'trace deficit above',gap.lower().str(12))
    assert gap > 0
