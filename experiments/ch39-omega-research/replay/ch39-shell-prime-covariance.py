"""Exact windowed covariance of the six readouts, including finite prime terms."""
from itertools import product
from pathlib import Path
import runpy
from fractions import Fraction as F
from flint import arb,acb

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-local-prime.py")))
T=acb.zeta_zero(719).imag
D=m["denominator"](T)
def observables(r,ds,c2):
    angles=[(1-arb(d)/(2*r)).acos() for d in ds]
    p,q=[2*(u-angles[0]) for u in angles[1:]]
    c=(arb(c2.numerator)/c2.denominator).sqrt()
    harmonics=(p,q,q-p,p+q)
    rows=((-arb(1)/2,0,0,0),(0,-arb(1)/2,0,0),
          (-(1-c**2)/4,-(1-c**2)/4,-(1+c)**2/8,-(1-c)**2/8))
    return harmonics,rows
l,a = observables(17,(6,12,14),F(81,98))
r,b = observables(867,(244,724,804),F(8741547,12563572))
def mean(f):
    arch,prime,hits=m["cosine"](T,f)
    return (arch+prime)/D,arch/D,hits
L,R = [[mean(z) for z in side] for side in (l,r)]
cov, smooth, hits = {},{},{}
for i,j in product(range(4),repeat=2):
    plus,minus=mean(l[i]+r[j]),mean(abs(l[i]-r[j]))
    for k,e in ((1,plus),(-1,minus)):
        if e[2]:
            hits[i,j,k]=e[2][0]
    cov[i,j]=(plus[0]+minus[0])/2-L[i][0]*R[j][0]
    smooth[i,j]=(plus[1]+minus[1])/2-L[i][1]*R[j][1]
assert hits=={(1,1,1):5,(1,3,-1):3,(3,3,1):19,(3,2,-1):3,
              (2,0,-1):2,(2,1,1):3,(2,3,-1):5}
out,rem={},{}
for i,j in product(range(3),repeat=2):
    out[i,j],rem[i,j]=[sum(a[i][s]*b[j][u]*matrix[s,u]
                            for s,u in product(range(4),repeat=2)) for matrix in (cov,smooth)]
for i in range(3):
    print("cross covariance row",i,*[out[i,j].str(13) for j in range(3)])
assert max(abs(v).upper() for v in rem.values())<arb("0.000001")
assert all(v.rad()<arb("1e-13") for v in out.values())
print("maximum archimedean-only covariance magnitude",max(abs(v).upper() for v in rem.values()).str(12))
print("sum-frequency prime 5 contribution to covariance of second readouts",
      (out[1,1]-rem[1,1]+(L[1][0]-L[1][1])*R[1][0]/4).str(13))
