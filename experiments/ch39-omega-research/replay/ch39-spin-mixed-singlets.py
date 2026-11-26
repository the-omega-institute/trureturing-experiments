"""Match with only the existing three neutral copies and two low spectral levels."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
h=runpy.run_path(str(root/'ch39-spin-balanced-higgs.py'))
q,D=h['q'],h['D']
y=h['y'];lh=h['lh']
assert max(v.upper() for v in y[:2]) < min(v.lower() for v in y[2:])
pair=min(v.lower() for v in y[:2]).union(max(v.upper() for v in y[:2]))
# Equal-modulus rank-one vector on six neutral directions; leading eigenvalues +p,-p/5.
overlap=-pair.log()/2
norm=q('6/5')
p=(D*(lh-h['lam'])/(norm*overlap)).sqrt()
pure=h['ferm']*h['g2']**2
for v in y:
    l=v.log()
    pure+=v*((q('78/25')+3*q('12/25'))*(3*l*l/2+2*l+q('7/2'))+
             q('156/25')*(1-6*l))*h['g2']**2
# Bound the divided mass sunset on both orders of the pair, including equal masses.
z=(pair.lower()/pair.upper()).lower().union(q(1))
l=pair.log()
sun=pair*((1+z)*((l-2)**2+1)/2-z*(2-l)*z.log()+
          (1-z)*q(0).union(arb.pi()**2/6))
hard=pure+norm*p**2*sun
b=(hard+h['msq']*(D/h['Q'])**2)/(29*D*h['response'])
sens=1+hard*(h['Q']/D)**2/h['msq']
neg=p/5;c=neg**2/(4*lh)+q('.0005')
diag=12*c
lim=q('24')
for sign in (-1,1):
    assert lim-sign*diag>0
    schur=lim-sign*6*lh-10*neg**2/(lim-sign*diag)-(2*p**2+18*58*b**2)/lim
    assert schur>0
frac=29*b**2*h['W']/(D*(lh-h['lam']))
for name,x in [('positive eigenvalue',p),('charged subleading portal',b),('radial negative block',diag),
               ('mass sensitivity',sens),('omitted portal-square / leading quartic',frac)]:
    print(name,x.lower().str(12),x.upper().str(12))
assert 0<b and 2*c<lim and p<lim and b<lim
# A flat traceless block with d<=4 cannot realize this quartic and the 8 pi contact bound.
# Its two-norm is d*p^2/(d-1), while the negative radial state tests 2*d*(p/(d-1))^2/(4*lambda_H).
logmax = -min(v.lower() for v in y).log()/2
smallblock = D*(lh-h['lam'])/(2*lh*3*logmax)
assert smallblock>8*arb.pi()
print('at most four neutral directions, necessary radial eigenvalue above',smallblock.lower().str(12))
