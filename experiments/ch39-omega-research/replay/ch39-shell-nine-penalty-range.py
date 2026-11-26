"""The order-nine ground index over an interval of height penalties, with its sharp upper endpoint."""
from functools import lru_cache
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb,arb

home=Path(__file__).resolve().parent
m=run_path(str(home/"ch39-shell-four-free-references.py"))
ref=m["m"]["reference"]
distance=m["m"]["distance"]
layout=m["layout"]
small=[tuple(ref(r,t[i],t[j]) for r,(i,j) in layout.items())
       for t in product(range(9),repeat=3) if t<=tuple(-v%9 for v in t)]
big=[ref(867,p,q) for p,q in product(range(9),repeat=2) if (p,q)<=(-p%9,-q%9)]


def minimum(vals):
    lo=min(x.lower() for x in vals)
    hi=min(x.upper() for x in vals)
    return lo.union(hi)


@lru_cache(None)
def point(n):
    h=acb.zeta_zero(n).imag
    obs=m["state"](h)
    a=minimum([sum(distance(obs[r],ref) for r,ref in zip(layout,target)) for target in small])
    b=minimum([distance(obs[867],target) for target in big])
    return a+b,h


e,t=point(25)
lines={}
for n in range(1,25):
    en,tn=point(n)
    lines[n]=(en-e)/(t*t-tn*tn)
idx=min(lines,key=lambda n:float(lines[n].mid()))
cap=lines[idx]
assert cap>0 and all(v>cap for n,v in lines.items() if n!=idx)
print("sharp upper endpoint",cap.str(16),"tied earlier zero",idx)
lam=arb(1)/51**7
assert arb(1)/51**4<cap and idx==8
E=e+lam*t*t
hi=int((E/lam).sqrt().upper().ceil().unique_fmpz())
assert lam*hi**2>E
found=[]
def visit(lo,hi):
    h=lo.union(hi)
    floor=m["paired_lower"](m["state"](h))+lam*h*h
    if floor>E:
        return
    if hi-lo>arb(1)/4:
        mid=(lo+hi)/2
        visit(lo,mid)
        visit(mid,hi)
        return
    nl,nh=(x.zeta_nzeros().unique_fmpz() for x in (lo,hi))
    assert nl is not None and nh is not None
    for j in range(max(25,int(nl))+1,int(nh)+1):
        t2=acb.zeta_zero(j).imag
        if m["paired_lower"](m["state"](t2))+lam*t2*t2>E:
            continue
        score,t2=point(j)
        if score+lam*t2*t2<E:
            found.append(j)
        else:
            assert score+lam*t2*t2>E
for lower in range(int(t.floor().unique_fmpz()),hi):
    visit(arb(lower),arb(lower+1))
assert not found
print("lower test",lam.str(13),"height cut",hi,"later strictly lower competitors",found)
