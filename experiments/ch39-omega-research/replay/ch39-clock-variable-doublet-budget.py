"""Dual certificate even when each matter layer chooses a different group lift."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from flint import arb
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-large-doublet-budget.py")))
weights=[F(n,10837) for n in (7917,280,0,2184,456)]
opt=F(16380,10837)
assert sum(weights)==1
def coeff(p,q,u,v,r,s,relaxed=False):
    # Cast all sums in rational arithmetic.
    _,b=m["examine"](p,q,u,v,r,s)
    fac=F(1,2) if relaxed and (p!=q or u!=v or (1+r+s)%2!=0) else F(1)
    return [fac*6*x/(11*n) for x,n in zip(b,m["ca"])]
def weighted(vec):
    return sum(w*x for w,x in zip(weights,vec))
limits=[]
for axis in range(6):
    z=[0]*6
    while weighted(coeff(*z,relaxed=True))<=opt:
        z[axis]+=1
    limits.append(z[axis])
# Each Casimir and dimension is nondecreasing in each Dynkin label. Use a relaxed real field
# for every partial product, then restore the actual type at the leaves.
visited=0
def verify(prefix):
    global visited
    padded=prefix+(0,)*(6-len(prefix))
    if weighted(coeff(*padded,relaxed=True))>opt:
        return
    if len(prefix)<6:
        for n in range(limits[len(prefix)]):
            verify(prefix+(n,))
    elif m["target"](*prefix):
        visited+=1
        assert weighted(coeff(*prefix))>=opt
verify(())
mix=((495,(0,0,0,0,1,12)),(1053,(0,0,0,0,9,0)),(280,(0,0,0,9,0,0)),(9009,(3,0,0,0,1,0)))
assert sum(n for n,_ in mix)==10837
assert all(m["target"](*rep) for _,rep in mix)
mean=[sum(F(n,10837)*coeff(*rep)[i] for n,rep in mix) for i in range(5)]
assert max(mean)==opt and weighted(mean)==opt
print("dual bound and saturating layer fractions",opt,[n for n,_ in mix],"out of 10837")
print("exclusive upper label bounds in relaxed enumeration",limits,"target leaves",visited)

# One rank-one adjoint Yang-Mills tower per factor, optimistically pure YM at every
# vector threshold. Upper bound includes one additional copy for ordinary gauge fields.
def ball(q):
    q=F(q)
    return arb(q.numerator)/q.denominator
c=[6*w/(11*n) for w,n in zip(weights,m["ca"])]
allowance=c[1]+max(c[0]*F(3,16),c[3],c[4]*F(9,4))
budget=2*arb.pi()*arb("40.9906")*ball(allowance)
# Both scalar graphs have equal total length. Use a certified lower Weyl slope.
spectrum=runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
kappa=2*spectrum["r"]*sum(spectrum["total"])/arb.pi()
def integrated(x):
    return ball(opt-1)*kappa*(x-1)-ball(270*opt+40)*x.log()
lo,hi=arb(1),arb(25)
assert integrated(hi)>budget
for _ in range(40):
    mid=(lo+hi)/2
    if integrated(mid)>budget:
        hi=mid
    else:
        lo=mid
assert integrated(hi)>budget and hi<15
print("weighted inverse-coupling budget and pole upper E/Q",
      budget.upper().str(13),hi.upper().str(13))
