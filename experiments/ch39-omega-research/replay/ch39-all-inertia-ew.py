"""Uniform low SM matching target and pole floor for the broad all-inertia crossing."""
from pathlib import Path
import runpy
root=Path(__file__).parent
s=(root/'ch39-ew.py').read_text().split('interval=(q(1000)')[0]
g={}
exec(compile(s,'SM_low','exec'),g)
q=g['q'];ev=g['evolve'];base=g['low_input']
import sys
global_crossings='--all-crossings' in sys.argv
lo,hi=(q('26.85181'),q('26.88167')) if global_crossings else (q('26.837'),q('26.899'))
interval=(q(1000)/q('91.1876')).log()
mz=ev(base,-interval)
lam=q('.5611')/4
v=q('248.404')
te=ev(mz+[lam,(lam*v*v).log()],interval)
start=ev(base+te[12:],lo)
state=ev(start,hi-lo,True)
ll=state[12];m2=state[13].exp()
print('SM quartic and m2',ll.lower().str(13),ll.upper().str(13),m2.lower().str(13),m2.upper().str(13))
lambda_upper='-.02980' if global_crossings else '-.0295'
assert ll<q(lambda_upper) and m2>0
inertia=runpy.run_path(str(root/'ch39-quaternion-all-inertias.py'))
width='3.272' if global_crossings else '3.29'
_,_,L=inertia['interval'](q('49.01286') if global_crossings else q('49.222'))
assert L<q(width)
p=(root/'ch39-spectrum-width-obstruction.py').read_text().split("if __name__")[0]
alpha_ceiling='38.3451' if global_crossings else '38.364'
p=p.replace('38.3316',alpha_ceiling).replace('-.029827',lambda_upper)
bound={}
exec(compile(p,'width_bound','exec'),bound)
_,_,floor=bound['bound'](width)
assert floor>0
Q=1000*lo.exp()
result=floor*(Q/(16*g['arb'].pi()**2))**2
print('all inertias mass floor',result.lower().str(15),'dimensionless',floor.lower().str(12))
