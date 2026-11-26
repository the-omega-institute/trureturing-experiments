"""Joint logarithmic exposure and second-cutoff obstruction across a continuous scalar scale range."""
from pathlib import Path
import runpy

root = Path(__file__).resolve().parent
flow = runpy.run_path(str(root/'ch39-all-right-crossings.py'))
q = flow['q']
data = flow['model']['data']
xs = data['xs']
rlo, rhi = q(127),q(450)
assert q('1.2710')*400**2 > rhi**2 and max(x.upper() for x in xs) < rlo**2
smin = q('3.5')
assert smin*(1+smin**2)*q(1).exp() < rlo
assert q(2).sqrt()*(1+2*smin**2)/(1+smin**2) > q(1).exp()
for a,b in ((724,216),(804,292)):
    k = data['theta'](a)/data['theta'](b)
    assert k*(1+k**2*smin**2)/(1+smin**2) > q(1).exp()
# First-exit bound used in crossings must cover the largest allowed WL.
g = flow['g']
loss = ((g['d'][1]+flow['value'](g['BD'][1],flow['trial']))*
        sum(rhi.log()-x.log()/2 for x in xs)+(8+flow['value'](g['BF'][1],flow['trial']))/3)
assert g['G0'][1]-2*loss/flow['D'] > q('1/2')
# Upper bound for the trace at u=1 with all 270 Dirichlet edge towers.
heat = q(38)
for n,factor in data['factors']:
    a = (-factor/(3*rhi**2)).exp()
    heat += n*a/(1-a)
u = 20
price = flow['read']['trace_price'](u,q(25))
ratio = data['theta'](804)/data['theta'](12)
arc_u = ratio**2
# Use a rational lower endpoint for u in the multiplier's comparison and
# low-left numerator at the upper endpoint, with the same exponential tail decomposition.
arc_floor = 72
assert arc_u > arc_floor and arc_u < 73
arc_price = flow['read']['trace_price'](arc_floor,q(25))
arc_gaps = []
rows = []
for lo in range(127,450,10):
    hi = min(450,lo+10)
    r = q(lo).union(q(hi))
    logs = [r.log()-x.log()/2 for x in xs]
    tau,w = flow['crossings'](logs)
    assert tau > 25
    numerator = sum((-u*x/q(lo)**2).exp() for x in xs)-(1-q(u)).exp()*heat
    gap = numerator-price*w
    # The high-right contribution is bounded at u=72; positivity gives S_R(arc_u)<=S_R(72).
    arc_low = sum((-arc_u*x/q(lo)**2).exp() for x in xs)
    arc_gaps.append(arc_low-(1-q(arc_floor)).exp()*heat-arc_price*w)
    rows.append((tau,w,gap))
    print('r interval',lo,hi,'right log upper',w.upper().str(9),'trace difference lower',gap.lower().str(9))
assert all(row[2]>0 for row in rows)
print('uniform trace deficit',min(row[2].lower() for row in rows).str(12))

print('extreme arcs Q ratio',ratio.str(14),'heat difference lower',min(v.lower() for v in arc_gaps).str(12))
assert min(v.lower() for v in arc_gaps)>3
