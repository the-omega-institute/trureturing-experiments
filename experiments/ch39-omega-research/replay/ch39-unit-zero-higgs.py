"""Joint gauge and Higgs solution in the zero-288 window with the fixed pentagon coefficient."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
z=runpy.run_path(str(root/'ch39-unit-zero-crossing.py'))
h=runpy.run_path(str(root/'ch39-five-singlet-floor.py'))
q=h['q'];D=h['D'];lam=h['lh'];delta=h['delta']
xs=z['xs'];T=z['T'];r=z['lo'].union(z['hi']);y=[x/r**2 for x in xs]
assert T>q('26.99') and T<q('27.013') and y[1].upper()<y[2].lower()
_,_,W,alpha=z['flow']['run'](T,xs,8192,r)
g2=4*arb.pi()/alpha
tF=z['flow']['data']['m']['model']['geometry'](r)
hard=q('36/5')*g2**2*sum((1+12*t)*(-2*t).exp() for t in tF)
resp=q(0)
for v in y:
    l=v.log()
    hard+=v*g2**2*((q('78/25')+3*q('12/25'))*(3*l*l/2+2*l+q('7/2'))+q('156/25')*(1-6*l))
    resp+=v*(1-l)
pair=min(v.lower() for v in y[:2]).union(max(v.upper() for v in y[:2]))
zeta=(pair.lower()/pair.upper()).lower().union(q(1))
l=pair.log()
sun=pair*((1+zeta)*((l-2)**2+1)/2-zeta*(2-l)*zeta.log()+
          (1-zeta)*q(0).union(arb.pi()**2/6))
d=6
p=(D*delta/(q(d)/(d-1)*(-l)/2)).sqrt()
neg=p/(d-1)
hard+=q(d)/(d-1)*p**2*sun
Q=1000*T.exp();ms=h['msq']
b=(hard+ms*(D/Q)**2)/(29*D*resp)
lim=8*arb.pi()
diag=2*d*(neg**2/(4*lam)+q('.0001'))
for sign in (-1,1):
    schur=(lim-sign*diag)*(lim-sign*6*lam-(2*p**2+1044*b**2)/lim)-2*(d-1)*neg**2
    assert lim-sign*diag>0 and schur>0
small=D*delta/(8*lam*(-pair.lower().log()/2))
assert small > lim and p<lim and b>0 and b<lim and diag/d<lim
sens=1+hard*(Q/D)**2/ms
for name,v in [('positive portal',p),('charged portal',b),('unified inverse',alpha),('sensitivity',sens),('rank five lower',small)]:
    print(name,v.lower().str(12),v.upper().str(12))
