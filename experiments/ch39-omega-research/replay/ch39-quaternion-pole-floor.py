"""Mass floor for the full-flavor clock example, allowing all Hermitian trace-free portals."""
from pathlib import Path
from flint import arb
import runpy
m=runpy.run_path(str(Path(__file__).with_name('ch39-quaternion-gauge.py')))
q=m['q'];r=q('31.0740').union(q('31.6548'))
alpha=q('38.3183').union(q('38.3273'))
g2=4*arb.pi()/alpha
Q=1000*q('26.865').union(q('26.873')).exp()
M=Q/r
xs=m['xs']
tau=m['geometry'](r)
D=16*arb.pi()**2
scalar=arb(0)
for x in xs:
    l=x.log()-2*r.log()
    R=(l-2)**2+1
    bracket=q('78/25')*(q('3/2')*l*l+2*l+q('7/2'))-q('216/25')*(3-l)**2/R
    assert l<-4 and bracket>0
    scalar+=x*bracket
fermion=q('36/5')*r*r*sum((1+12*t)*(-2*t).exp() for t in tau)
mass=g2**2*(scalar+fermion)*(M/D)**2
print('spectral scale M',M.lower().str(13),M.upper().str(13))
print('positive pole floor GeV^2',mass.lower().str(14))
assert mass>q('1e25')
