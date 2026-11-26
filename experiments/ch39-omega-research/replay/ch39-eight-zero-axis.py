"""An actual zero with exact arc ratio, with the relative clock axis determined by matching."""
from pathlib import Path
import runpy
from flint import acb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
q,pi=m['q'],m['pi']
gamma=acb.zeta_zero(8).imag
u,v=m['shell']['a']*gamma,m['shell']['b']*gamma
t0,t1=q('26.98'),q('27.025')
angle0,angle1=q('.5').acos(),q('-.5').acos()
T=t0.union(t1)
clock_width=q('.0001')
for angle,sign in ((angle0,1),(angle1,-1)):
    xs=m['spectrum'](u,clock_width,True,v,angle)
    _,res,w,_=m['run'](T,xs,8192)
    print('axis endpoint, required residual margin',angle.str(10),(sign*res).lower().str(12))
    assert sign*res>0
ncells=16
margins=[]
vs=[]
for j in range(ncells):
    half=(angle1-angle0)/(2*ncells)
    # A clock-angle change has block norm <=delta; the axis change <=half on 804 only.
    # At most four clocked low-low edges in each row.
    xs=m['spectrum'](u,clock_width+half,False,v,angle0+(2*j+1)*half)
    ml=m['run'](t0,xs)[0]
    mr=-m['run'](t1,xs)[0]
    margins.append(min(ml.lower(),mr.lower()))
    vs.append(m['run'](T,xs)[3])
print('opposite time face margin',min(margins).str(12))
assert min(margins)>0
print('eighth zero',gamma.str(18),'unified inverse',min(x.lower() for x in vs).str(12),max(x.upper() for x in vs).str(12))
