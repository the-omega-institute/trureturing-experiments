"""Four-momentum concavity and a finite-mesh global-minimum obstruction to triangle actions."""
from pathlib import Path
from fractions import Fraction
import runpy

root=Path(__file__).resolve().parent
mesh=runpy.run_path(str(root/'ch39-gaussian-mesh-momentum.py'))
m=mesh['base']
real=m['real']
def curvature_tail(K):
    K=real(K)
    adj=[{b if i==a else a for a,b,d,z in m['shell']['edges'] if d in (724,804) and i in (a,b)}
         for i in range(19)]
    four=sum(len(a&b)**2 for a in adj for b in adj)
    # n>=5, 2 d^2(d^2-1)/(3n(n-2)) <= n^2/12, n^2 <= 25*(6/5)^(2(n-5)).
    large=47970/K**2*(real(25)/12)*(18/K)/(1-real('1.44')*18/K)
    return four/K**2+large

if __name__=='__main__':
    lo,hi=m['quadrature'](order=2,tail_control=curvature_tail)
    # Retain the continuum finite-range concavity, then use the mesh's own tail bound.
    error=real(0)
    for a,b in [(0,1)]+[(2**j,2**(j+1)) for j in range(10)]:
        for i in range(256):
            l=real(Fraction(a)+Fraction((b-a)*i,256))
            r=real(Fraction(a)+Fraction((b-a)*(i+1),256))
            error+=(r*r-l*l)*mesh['mesh_slope_error'](l,r)/2
    # The truncation used in deriving lo already allows one worst-case curvature tail.
    gap=lo/2-error
    print('finite mesh gap above endpoint chord divided by 1-c^2',gap.lower().str(12))
    assert gap>0
