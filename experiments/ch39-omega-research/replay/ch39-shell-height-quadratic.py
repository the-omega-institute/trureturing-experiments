"""Exact interval of height penalties for a coarse rational Wilson target to pick zero 719."""
from pathlib import Path
from flint import arb,acb
import runpy

data=runpy.run_path(str(Path(__file__).with_name("ch39-shell-earlier-zero-dominance.py")))
observe=data["obs"]
h=acb.zeta_zero(719).imag
read=observe(h)[1]
def rounded(d):
    return [int((d*x+arb(1)/2).floor().unique_fmpz()) for x in read]
def loss(values,nums,d):
    return sum((x-arb(k)/d)**2 for x,k in zip(values,nums))
for d,n in ((1,583),(2,657)):
    v=rounded(d)
    assert loss(observe(acb.zeta_zero(n).imag)[1],v,d)<loss(read,v,d)
v=rounded(3)
u=loss(read,v,3)

earlier=[]
for start in range(1,719,80):
    for n,z in enumerate(acb.zeta_zeros(start,min(80,719-start)),start=start):
        zz=z.imag
        val=(loss(observe(zz)[1],v,3)-u)/(h*h-zz*zz)
        earlier.append((val,n))
earlier.sort(key=lambda t:float(t[0].mid()))
hi,nhi=earlier[0]
assert nhi==207 and all(hi<x for x,_ in earlier[1:])
print("third-grid numerators",v,"target loss",u.str(15),"upper penalty",nhi,hi.str(16))
assert arb(1)/51**5<hi<arb(1)/51**4
