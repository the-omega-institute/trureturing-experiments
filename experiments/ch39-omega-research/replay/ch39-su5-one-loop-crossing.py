"""One-loop running with the SU(5) vector thresholds and an independently specified matter spectrum."""
from pathlib import Path
import runpy
from flint import arb, acb

root = Path(__file__).resolve().parent
vectors = runpy.run_path(str(root/'ch39-su5-vector-threshold.py'))
f = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
q,r,D,flo = f['q'],f['r'],f['D'],f['flow']
u,v = (f['shell'][key]*acb.zeta_zero(8).imag for key in ('a','b'))
g = flo['g']
dot, project = flo['m']['dot'],flo['m']['project']
pt,pw = project(flo['d'],flo['beta']),project(flo['beta'],flo['d'])
target = [D*x/2-q(8)/3 for x in g['G0']]
indices = [5,3,2]
tp,wp = (q('7/2')*dot(vector,indices) for vector in (pt,pw))
base_t,base_w = dot(pt,target),dot(pw,target)
shape = None
for j in range(48):
    c0,c1 = (q('-.1855') + q('.0585')*i/48 for i in (j,j+1))
    left,right = c0.acos(),c1.acos()
    xs = f['spectrum'](u,(left-right)/2+q('.00001'),True,v,(left+right)/2)
    cell = sum(x.log()/2 for x in xs)
    shape = cell if shape is None else shape.union(cell)
# S_m=18 log(r)-shape, H=log(r)-log_product/2.
ell = vectors['log_product']
logr = (base_w+shape-wp*ell/2)/(18-wp)
ratio = logr.exp()
exposure = logr-ell/2
tau = base_t+tp*exposure
assert ratio**2>36 and ratio**2<1/(2*vectors['rho'])
assert q(26) < ratio**2 and ratio**2 < f['data']['edge19']
levels = f['data']['m']['model']['geometry'](ratio)
assert min(levels)>0 and max(levels)<1 and sum(levels).contains(1)
assert q(26)<tau and logr < tau
w = 18*logr-shape
s0 = 18*logr-(q(int(19*vectors['trees'])).log()+(-18*vectors['error']).union(q(0)))/2
inv = 4*arb.pi()*(g['G0'][0]-2*(flo['beta'][0]*tau+flo['d'][0]*w+q(8)/3
      -q('7/2')*(5*s0+5*exposure))/D)
mass = 1000*tau.exp()/ratio
for name,val in [('scale ratio',ratio),('relative arc ratio',ratio/r),('log time',tau),
                 ('unified inverse',inv),('mass scale in GeV',mass)]:
    print(name,val.lower().str(12),val.upper().str(12))
assert ratio/r>q('1.001') and inv>1
# Dropping vector antiscreening leaves a convex cumulative loss, maximal at an endpoint.
for i in range(3):
    upper_loss = max(q(0),(flo['beta'][i]*tau+flo['d'][i]*w+q(8)/3).upper())
    assert 4*arb.pi()*(g['G0'][i]-2*upper_loss/D)>1

