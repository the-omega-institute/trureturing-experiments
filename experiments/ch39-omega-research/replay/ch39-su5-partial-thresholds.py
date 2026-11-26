"""Exclude even partially active low levels below the first possible additional internal mode."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-su5-all-clocks.py'))
v = m['local']['vectors']
q,f = m['q'],m['f']
# Eliminate the three degree-four outer pairs on the static low block.
# Their eigenvalues are 3 and 5, for both flavors. Each inner site has six outer neighbors.
shape_cap = 3*(q(15)*q('6.8')).log()
pivot = (shape_cap/18).exp()
# Two static moments on the eighteen spin components. The true eigenvalues are bounded
# above, in order, by the static low-block eigenvalues by the endpoint extension principle.
degrees = [int(m['lap'][i,i]) for i in range(9)]
first = 2*sum(degrees)
second = 2*(sum(d*d for d in degrees)+sum(degrees))
assert (first,second)==(96,672)
mean = q(first)/18
sigma = (q(second)/18-mean*mean).sqrt()
excess = max(q(0),*((q(j)/2*((mean+sigma*(q(18-j)/j).sqrt())/pivot**2).log()).upper() for j in range(1,19)))
# For every cutoff R, W_q - max(0,18 log R - shape_cap) is between 0 and this bound
# for the positive-projection branches. The maximum of the difference is at the pivot.
spread = shape_cap-m['shape_lower'][2]+excess
required = m['base']-q('67/3')*spread
print('all cutoffs required vector contrast lower',required.lower().str(12))
assert required>q('4.47')
neutral = []
for z,mult in v['neutral'].charpoly().complex_roots():
    if z==0:continue
    assert z.real>0
    neutral.extend([z.real.upper().sqrt().upper()]*mult)
charged = []
for z in v['g']['matrix'].eig(multiple=True):
    s=z.real.lower()
    charged.append((s/(1+36*v['rho']+2*v['rho']*s)).sqrt().lower())
assert len(neutral)==18 and len(charged)==19
def exposure(r,masses):
    return sum(max(q(0),(r/s).log().upper()) for s in masses)
def contrast(lo,hi):
    low = sum(max(q(0),(lo/s).log().lower()) for s in neutral)
    return exposure(hi,charged)-low
# All these piecewise-log-linear functions can be enclosed by their boundary knots.
knots = [q(6),*charged,*neutral]
maximum = max(contrast(x.lower(),x.upper()).upper() for x in knots)
print('partial vector contrast upper to 6',maximum.str(12))
assert maximum<q(2) and min(charged)>1
dirichlet = 1/(v['rho']*(1/arb.pi()**2+q('.01')**2/6))
ceiling = min(dirichlet.sqrt().lower(),f['data']['edge19'].sqrt().lower())
largest = ceiling.log()-v['log_product'].lower()/2
assert largest<m['flux']
mass_min = min(m['floors'].values())
time_floor = (ceiling/mass_min.sqrt()).log()
print('low-band scale ceiling, saturation input time, full-activation vector upper',
      ceiling.str(12),time_floor.upper().str(12),largest.upper().str(12))
