"""Bounds on the scalar logarithmic width for every positive common inertia given total exposure."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-scalar-levels.py'))
q=m['m']['m']['src']['theta']
ctxsrc=m['m']['m']['real_energy_matrix'].__globals__['src']
dstar=arb(256)
ctxsrc['rho']=arb.pi()**2/(dstar*ctxsrc['theta'](804))
xmin=arb('1.5')
shift=arb('.041')
assert m['count'](xmin+shift)==0
k=((xmin+shift)*arb.pi()**2/dstar).sqrt()
assert 4*arb('.01')*k/k.sin()<shift
# All inertias have eigenvalues bounded above by the fixed-endpoint static trial levels.
raw=(root/'ch39-quaternion-band-boxes.txt').read_text().splitlines()
vmax=sum((arb(int(line.split()[1]))/10**6).log()/2 for line in raw)/18
def interval(w):
    r_global=(w/18+vmax).exp()
    lowst=(r_global/xmin.sqrt()).log()
    s=arb('1.62')
    def term(d):
        x=s*(q(d)/q(804)).sqrt()
        return x.log() if x>1 else arb(0)
    WD=2*sum(term(d) for _,_,d,_ in ctxsrc['edges'])
    assert WD>w
    dense=(s*(dstar/xmin).sqrt()).log()
    upper=max(lowst.upper(),dense.upper())
    # With a single active fermion its lowest mass has to be Q/e.
    # s<3 when r/e<30, hence m2/m1 <= sqrt(2)*19/10 < e.
    assert r_global/arb(1).exp()<30 and arb(2).sqrt()*19/10<arb(1).exp()
    return r_global,WD,upper
if __name__=='__main__':
    for W in ['49.5','50','51','52']:
        r,wd,L=interval(arb(W))
        print(W,r.upper().str(11),wd.lower().str(11),L.upper().str(11))
