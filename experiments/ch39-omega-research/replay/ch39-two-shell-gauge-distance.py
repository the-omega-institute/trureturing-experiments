"""Gauge-conjugacy mismatch for the embedded R=17 and R=867 neutral clock pairs."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-book-shell-clock-embedding.py")))
q,pi=m["q"],arb.pi()
c1,c2=m["axis"],m["m"]["m"]["axis"]
s1,s2=(1-c1*c1).sqrt(),(1-c2*c2).sqrt()
th1,th2=c1.acos(),c2.acos()
delta=abs(th1-th2)
total=th1+th2
tau=2*(delta/4).sin()
area=(pi-delta)/delta.sin()+total/total.sin()
coef=8*area/pi**3
p,q1,P,Q=[m["h"]*v for v in m["f"]]
def angle_data(x,y):
    a=x.cos()*y.cos()
    b=x.sin()*y.sin()
    return a,b,abs(a)+abs(b)
ap,bp,mp=angle_data(p,P)
aq,bq,mq=angle_data(q1,Q)
def positive(x):
    assert x>0 or x<0
    return 1 if x>0 else -1
def feasible(value):
    for ep in (-1,1):
        for eq in (-1,1):
            cosines=[(value-e*a)/abs(b) for e,a,b in ((ep,ap,bp),(eq,aq,bq))]
            if any(x>1 for x in cosines):
                continue
            assert all(x<1 for x in cosines)
            caps=[pi if x < -1 else x.acos() for x in cosines]
            angle=th2 if ep*eq*positive(bp*bq)>0 else pi-th2
            gap=sum(caps)-abs(th1-angle)
            assert not gap.contains(0)
            if gap>0:
                return True
    return False
lo,hi=q(0),q(1)
for _ in range(50):
    mid=(lo+hi)/2
    if feasible(mid):
        lo=mid
    else:
        hi=mid
assert feasible(lo) and not feasible(hi)
best=lo.union(hi)
distance=(2-2*best).sqrt()
print("unoriented angle discrepancy and minimax axis chord",delta.str(13),tau.str(13))
print("small-tolerance coefficient in the two-shell conjugacy density",coef.str(13))
print("exact minimax projective-quaternion distance at 719",distance.str(13))
