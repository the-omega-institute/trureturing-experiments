"""Uniform six-direction neutral Higgs matching at the eighth zero's exact-arc-ratio crossings."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
z=runpy.run_path(str(root/'ch39-eight-zero-axis.py'))
m=z['m'];q,D,r=m['q'],m['D'],m['r'];T=z['T']
h=runpy.run_path(str(root/'ch39-spin-balanced-higgs.py'))
low=h['evolve'](h['start']+h['init'][12:],q('26.98'))
low=h['evolve'](low,q('.045'),True)
lh=h['lh'];delta=lh-low[12];ms=low[13].exp();Q=1000*T.exp()
tF=m['data']['m']['model']['geometry'](r)
result=[]
for j in range(z['ncells']):
    half=(z['angle1']-z['angle0'])/(2*z['ncells'])
    xs=m['spectrum'](z['u'],z['clock_width']+half,True,z['v'],z['angle0']+(2*j+1)*half)
    y=[x/r**2 for x in xs]
    a=z['vs'][j]
    g2=4*arb.pi()/a
    rho=0
    hard=q('36/5')*g2**2*sum((1+12*t)*(-2*t).exp() for t in tF)
    for v in y:
        l=v.log();rho+=v*(1-l)
        hard+=v*g2**2*((q('78/25')+3*q('12/25'))*(3*l*l/2+2*l+q('7/2'))+q('156/25')*(1-6*l))
    pair=min(v.lower() for v in y[:2]).union(max(v.upper() for v in y[:2]))
    zz=(pair.lower()/pair.upper()).lower().union(q(1));l=pair.log()
    sun=pair*((1+zz)*((l-2)**2+1)/2-zz*(2-l)*zz.log()+(1-zz)*q(0).union(arb.pi()**2/6))
    p=(D*delta/(q('6/5')*(-l)/2)).sqrt()
    hard+=q('6/5')*p**2*sun
    b=(hard+ms*(D/Q)**2)/(29*D*rho)
    neg=p/5
    diag=12*(neg**2/(4*lh)+q('.0001'))
    lim=8*arb.pi()
    margins=[]
    for sign in (-1,1):
        schur=(lim-sign*diag)*(lim-sign*6*lh-(2*p**2+1044*b**2)/lim)-10*neg**2
        assert lim-sign*diag>0 and schur>0
        margins.append(schur.lower())
    assert b>0 and p<lim and b<lim and diag/6<lim
    result.append([p,b,1+hard*(Q/D)**2/ms,min(margins)])
for i,name in enumerate(('positive portal','charged portal','mass sensitivity','Schur margin')):
    print(name,min(row[i].lower() for row in result).str(12),max(row[i].upper() for row in result).str(12))
