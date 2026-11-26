"""Minimal flat neutral portal size for the two-loop spin-balanced crossing."""
from pathlib import Path
import math
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-spin-mixed-singlets.py'))
h,q = m['h'],m['q']
D,lh = h['D'],h['lh']
light = h['light']
def evolve(state,delta,whole=False):
    n=math.ceil(abs(float(delta))/.00025)
    dt=delta/n
    pads=[q('.0003')/8]*3+[q('.0006')/8]*9+[q('.00009')/8,q('.00015')/8]
    hull=None
    for _ in range(n):
        deriv=light['slopes']([v+arb(0,p) for v,p in zip(state,pads)])
        assert all(abs(dt*d)<p for d,p in zip(deriv,pads))
        if whole:
            tube=[v+arb(0).union(dt)*d for v,d in zip(state,deriv)]
            hull=tube if hull is None else [a.union(b) for a,b in zip(hull,tube)]
        state=[v+dt*d for v,d in zip(state,deriv)]
    return hull if whole else state
mz=evolve(h['start'],-h['dt'])
init=evolve(mz+[h['higgs'],(h['higgs']*q('248.404')**2).log()],h['dt'])
low=evolve(h['start']+init[12:],q('26.99'))
low=evolve(low,q('.023'),True)
delta=lh-low[12]
p=(D*delta/(q('5/4')*m['overlap'])).sqrt()
neg=p/4
hard=m['pure']+q('5/4')*p**2*m['sun']
msq=low[13].exp()
b=(hard+msq*(D/h['Q'])**2)/(29*D*h['response'])
lim=8*arb.pi()
c=neg**2/(4*lh)+q('.0001')
diag=10*c
for sign in (-1,1):
    gap=(lim-sign*diag)*(lim-sign*6*lh-(2*p**2+1044*b**2)/lim)-8*neg**2
    print('five directions Schur determinant sign',sign,gap.lower().str(12),gap.upper().str(12))
    assert lim-sign*diag>0 and gap>0
for name,x in [('positive portal',p),('charged portal',b),('light lambda difference',delta),('radial diagonal',diag)]:
    print(name,x.lower().str(12),x.upper().str(12))
