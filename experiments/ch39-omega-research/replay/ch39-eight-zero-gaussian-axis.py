"""Test whether the zero-energy shell Gaussian determinant can dynamically select the crossing axis."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
fixed=runpy.run_path(str(root/'ch39-eight-zero-short-arc.py'))
m=runpy.run_path(str(root/'ch39-quaternion-spectral-polynomial.py'))
q=fixed['q']; pi=arb.pi()
clock_width=q('.00001')
u=fixed['u']+arb(0,clock_width)
v=fixed['v']+arb(0,clock_width)
su,x=u.sin_cos();sv,y=v.sin_cos()
for shift in (0,1,10):
    poly=m['polynomial'](shift)
    derivatives=[]
    for j in range(128):
        c=(q(j-64)/64).union(q(j+1-64)/64)
        z=x*y+c*su*sv
        ratio=1+m['evaluate'](poly,(x,y,z))
        grad=sum(a*k*x**i*y**n*z**(k-1) for (i,n,k),a in poly.items() if k)*su*sv
        assert ratio>0
        derivatives.append(grad/ratio)
    lo=min(w.lower() for w in derivatives)
    hi=max(w.upper() for w in derivatives)
    print('shift',shift,'log determinant slope in the axis dot product',lo.str(12),hi.str(12),flush=True)
    assert lo>0
for beta,sign in ((q(0),1),(pi,-1)):
    T,res=fixed['project_fixed_axis'](beta,clock_width)
    print('collinear axis',beta.str(9),'all-time oriented crossing gap',(sign*res).lower().str(12),flush=True)
    assert sign*res>0
