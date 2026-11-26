"""Count all coarse Wilson grid targets choosing zero 719 at the geometric height penalty."""
from itertools import product
from pathlib import Path
import runpy
import numpy as np
from flint import arb,acb

source=runpy.run_path(str(Path(__file__).with_name('ch39-shell-earlier-zero-dominance.py')))
observe=source["obs"]
lam=arb(1)/51**5
hstar=acb.zeta_zero(719).imag
xstar=observe(hstar)[1]
def difference(x,h):
    return (sum(y*y for y in x)-sum(v*v for v in xstar)+lam*(h-hstar)*(h+hstar),
            [2*(t-v) for t,v in zip(x,xstar)])
first=[]
for k in range(1,1701,96):
    for z in acb.zeta_zeros(k,min(96,1701-k)):
        first.append(difference(observe(z.imag)[1],z.imag))
const=np.array([float(a.mid()) for a,_ in first])
coef=np.array([[float(v.mid()) for v in b] for _,b in first])
def check_grid(d):
    nums=list(product(range(d+1),repeat=6))
    dirs=np.array(nums,dtype=float).T/d
    cheapest=(const[:,None]-coef@dirs).argmin(axis=0)
    survivors=[]
    for pt,i in zip(nums,cheapest):
        if i==718:
            survivors.append(pt)
            continue
        a,b=first[int(i)]
        assert d*a < sum(v*k for v,k in zip(b,pt))
    return survivors
grids=[check_grid(d) for d in (1,2,3)]
print("targets not eliminated by first 1700 zeros",list(map(len,grids)))
print("provisional candidates",grids)
def remaining(pt):
    assert all(3*a>sum(v*k for v,k in zip(b,pt)) for i,(a,b) in enumerate(first) if i!=718)
    r=[arb(z)/3 for z in pt]
    def energy(x,t):
        return sum((u-v)**2 for u,v in zip(x,r))+lam*t**2
    winner=energy(xstar,hstar)
    start=int(acb.zeta_zero(1700).imag.floor().unique_fmpz())
    cap=int((winner/lam).sqrt().upper().ceil().unique_fmpz())
    def visit(a,b):
        # Use the lower endpoint for the increasing term.
        if energy(observe(a.union(b))[1],a)>winner:
            return None
        if b-a>arb(1)/8:
            mid=(a+b)/2
            return visit(a,mid) or visit(mid,b)
        left=a.zeta_nzeros().unique_fmpz()
        right=b.zeta_nzeros().unique_fmpz()
        assert left is not None and right is not None
        for idx in range(max(1701,int(left)+1),int(right)+1):
            h=acb.zeta_zero(idx).imag
            candidate=energy(observe(h)[1],h)
            if candidate<winner:
                return idx
            assert candidate>winner
        return None
    for i in range(start,cap):
        witness=visit(arb(i),arb(i+1))
        if witness:
            return witness
    assert lam*cap**2>winner
    return None
final={}
for pt in grids[2]:
    final[pt]=remaining(pt)
print("global result (None denotes a strict 719 win)",final)

assert len(grids[0])==len(grids[1])==0 and len(final)==3 and all(z is None for z in final.values())
