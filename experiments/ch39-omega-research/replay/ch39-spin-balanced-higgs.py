"""Gauge-singlet portals in the spin-balanced shell at its certified two-loop crossing."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-spin-balanced-rge.py'))
q,D = m['q'],m['D']
T=q('26.99').union(q('27.013'))
W=q('24.925').union(q('25.07'))
r=q('8.37369').union(q('8.53494'))
g2=4*arb.pi()/q('29.99894').union(q('30.07107'))
Q=1000*T.exp()
source=(root/'ch39-ew.py').read_text().split('interval=(q(1000)')[0]
light={}
exec(compile(source,'SM_light','exec'),light)
evolve=light['evolve']
start=light['low_input']
dt=(q(1000)/q('91.1876')).log()
mz=evolve(start,-dt)
higgs=q('.5611')/4
init=evolve(mz+[higgs,(higgs*q('248.404')**2).log()],dt)
endpoint=evolve(start+init[12:],q('26.99'))
low=evolve(endpoint,q('.023'),True)
lam,msq=low[12],low[13].exp()
lh=q('.02335')
# A_L=2 sum C_H T(L)=78/25, A_R=12/25, Dirac L carries 4 sum C_H T(L).
y=[x/r**2 for x in m['low']]
tF=m['m']['data']['m']['model']['geometry'](r)
ferm=q('36/5')*sum((1+12*t)*(-2*t).exp() for t in tF)
hard=ferm*g2**2
response=q(0)
for v in y:
    l=v.log()
    pure=(q('78/25')+3*q('12/25'))*(3*l*l/2+2*l+q('7/2'))
    pure+=q('156/25')*(1-6*l)
    # Per-layer sum of real portal squares divided by two is fixed by the quartic equation.
    portal=D*(lh-lam)/W
    hard+=v*(pure*g2**2+portal*((l-2)**2+1))
    response+=v*(1-l)
b=(hard+msq*(D/Q)**2)/(29*D*response)
assert b>0
sens=1+hard*(Q/D)**2/msq
for n in range(2,8):
    p=(D*(lh-lam)*(n-1)/(W*n)).sqrt()
    neg=p/(n-1)
    count=36*(n-1)
    minimal=(count+2)*neg**2/(4*lh)
    if n<=5:
        assert minimal>8*arb.pi()
    else:
        # c is chosen from the computed parameters; both normalized radial blocks
        # and all their orthogonal channels are bounded in the stated construction.
        c=neg**2/(4*lh)+q('.0005')
        diag=(count+2)*c
        lim=8*arb.pi()
        assert diag<lim
        for sign in (-1,1):
            margin=lim-sign*6*lh-count*neg**2/(lim-sign*diag)-(36*p**2+18*58*b**2)/lim
            assert margin>0
        assert 2*c<lim and p<lim and b<lim and neg<lim
    print('singlets',n,'positive leading portal',p.lower().str(10),p.upper().str(10),
          'negative-channel stability floor',minimal.lower().str(10),minimal.upper().str(10))
print('subleading charged scalar portal',b.lower().str(12),b.upper().str(12))
print('pole mass sensitivity lower',sens.lower().str(12))

ratio=29*b**2/(D*(lh-lam)/W)
print('relative portal-square correction to leading quartic upper',ratio.upper().str(12))
assert ratio<q('.0006')
