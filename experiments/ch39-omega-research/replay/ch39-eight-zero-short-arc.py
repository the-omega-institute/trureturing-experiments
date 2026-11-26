"""Test the separately prescribed shortest-arc rotation at the eighth zero, over every candidate time."""
from pathlib import Path
import runpy
from flint import acb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
a=runpy.run_path(str(root/'ch39-spin-arc-scale.py'))
q,r=m['q'],m['r']
gamma=acb.zeta_zero(8).imag
u,v=m['shell']['a']*gamma,m['shell']['b']*gamma
axis=m['data']['theta'](12)*gamma
def project_fixed_axis(axis,width):
    xs=m['spectrum'](u,width,True,v,axis)
    shapes=[x.log()/2 for x in xs]
    w=18*r.log()-sum(shapes)
    T=a['narrow'](a['T'],a['gs'],shapes,w)
    assert T>5 and T<46
    survivors=[]
    for j in range(5,46):
        trial=a['narrow'](q(j).union(q(j+1)),a['gs'],shapes,w)
        if trial is None:continue
        t,*_=m['run'](trial,xs,4096,return_target=True)
        if trial.overlaps(t):survivors.append(a['base']['cut'](trial,t))
    assert survivors
    T=min(v.lower() for v in survivors).union(max(v.upper() for v in survivors))
    for _ in range(3):
        t,res,_,_=m['run'](T,xs,16384,return_target=True)
        T=a['base']['cut'](T,t)
    return T,res
if __name__=='__main__':
    T,res=project_fixed_axis(axis,q('.00001'))
    print('cos of shortest-arc axis angle',axis.cos().str(13),'candidate times',T.lower().str(12),T.upper().str(12))
    print('shortest-arc exposure residual',res.lower().str(12),res.upper().str(12))
    assert res>0
