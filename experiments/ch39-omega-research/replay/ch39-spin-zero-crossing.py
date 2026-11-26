"""A certified crossing near the arc ratio at an actual arithmetic clock."""
from pathlib import Path
import runpy
from flint import acb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
q=m['q']
r=m['r']
shell=m['shell']
gamma=acb.zeta_zero(78).imag
width=q('0.00005')
xs=m['spectrum'](shell['a']*gamma,width,True,shell['b']*gamma)
t0,t1=q('26.994'),q('27.008')
T=t0.union(t1)
eps=q('.001')
lo,hi=r*(1-eps),r*(1+eps)
for time,ratio,k,sign in ((t0,lo.union(hi),0,1),(t1,lo.union(hi),0,-1),(T,lo,1,1),(T,hi,1,-1)):
    row=m['run'](time,xs,8192,ratio)
    print('clock crossing face margin',k,sign,(row[k]*sign).lower().str(12),flush=True)
    assert sign*row[k]>0
print('zero, gamma, scale bounds',78,gamma.str(18),lo.lower().str(13),hi.upper().str(13))
