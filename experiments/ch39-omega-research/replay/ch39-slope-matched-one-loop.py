"""One-loop necessary net log exposure, and a uniform moment obstruction below a scale ratio."""
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-kinetic-trace-crossing.py"))
g=runpy.run_path(str(root/"ch39-global-gauge-crossings.py"))
shell=runpy.run_path(str(root/"ch39-arc-relative-determinant.py"))
q=m["q"]
d=m["d"]
required=g["base"][1]
h=q(".01")
alpha=arb.pi()**2/d
beta=alpha*h*h/6
ret=(1-16*alpha*h*h/arb.pi()**2)**2-9*alpha*h
xs=[(x.lower()/(1+beta*(x.lower()+10))).lower().union((x.upper()/ret).upper())
    for x in m["model"]["xs"]]
E=q(446)**2
k=arb.pi()*(E/d).sqrt()
assert k<arb.pi() and 959990-E-36*k/k.sin()>0
edge19=(E/(1+beta*(E+10))).lower()
edges={}
for _,_,ds,_ in shell["edges"]:
    edges[ds]=edges.get(ds,0)+2
theta=shell["theta"]
factors=[(n,d*theta(804)/theta(ds)) for ds,n in edges.items()]

def maximum(rlo,rhi):
    r=q(rlo).union(q(rhi))
    yy=[x/r**2 for x in xs]
    heat=sum((-y).exp() for y in yy).lower()
    moment=sum(y*(-y).exp() for y in yy).upper()
    base=edge19/r.upper()**2
    assert base>1
    moment+=20*base*(-base).exp()
    for n,factor in factors:
        k0=factor/r.upper()**2
        for j in range(1,5):
            low=k0*j*j/(1+arb.pi()**2*j*j*h*h/6)
            assert low>1
            moment+=n*low*(-low).exp()
        # y exp(-y) <= exp(-y/2), y >= k0*j^2/3; j^2 >= 5*j for j>=5.
        decay=(-5*k0/6).exp()
        assert 25*k0/3>2
        moment+=n*decay**5/(1-decay)
    ratio=moment/heat+(-2*g["base"][0]).exp()
    right=max(q(0),(-heat*ratio.log()/2).lower())
    left=sum(max(q(0),(-y.log()/2).upper()) for y in yy)
    return (left-right).upper()

if __name__=="__main__":
    bounds=[]
    assert 18*q(5).log()<required
    for j in range(5,189):
        bounds.append(maximum(j,j+1))
    assert max(bounds)<required
    print("r <= 189 net upper",max(bounds).str(13),"one-loop requirement",required.str(13))
