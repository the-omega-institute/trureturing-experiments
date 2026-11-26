"""A zero-shift finite-graph decision at the known ordinate square, and direct block replay."""
from pathlib import Path
import runpy
from fractions import Fraction
from flint import arb,acb,acb_mat
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-spectral-polynomial.py'))
src=m['src'];geom=runpy.run_path(str(root/'ch39-arc-quaternion-clocks.py'))
c=geom['c']
x0,y0=geom['p']
def delta_reading(poly,x,y,z,z0):
    # Factor z-z0; use the same x,y for both physical hypotheses.
    return (z-z0)*sum(v*x**i*y**j*sum(z**k*z0**(r-1-k) for k in range(r))
                          for (i,j,r),v in poly.items() if r)
def window(width,poly,count):
    vals=[]
    for i in range(count):
        for j in range(count):
            offsets=[arb(width)*(2*arb(k+0.5)/count-1)+arb(0,(arb(width)/count).upper()) for k in (i,j)]
            (s,x),(t,y)=[(a+e).sin_cos() for a,e in zip((x0,y0),offsets)]
            vals.append(delta_reading(poly,x,y,x*y+c*s*t,x*y+s*t))
    return min(v.lower() for v in vals).union(max(v.upper() for v in vals))
poly=m['polynomial'](0)
for w,n in [('0',1),('.01',12)]:
    val=window(w,poly,n)
    print('change at zero shift in normalized determinant',w,val.lower().str(13),val.upper().str(13))
    assert not val.contains(0)
# Replay on the actual nonperpendicular quaternion clock matrices.
s,x=x0.sin_cos();t,y=y0.sin_cos();perp=(1-c*c).sqrt()
U=acb_mat([[acb(x,s),0],[0,acb(x,-s)]])
V=acb_mat([[acb(y,c*t),t*perp],[-t*perp,acb(y,-c*t)]])
for u in [0,15000,-1000,-204000,-215000]:
    predicted=m['evaluate'](m['polynomial'](u),[x,y,x*y+c*s*t])
    full=m['direct'](u,U,V)
    assert predicted.overlaps(full) and full.rad()<arb('1e-30')
