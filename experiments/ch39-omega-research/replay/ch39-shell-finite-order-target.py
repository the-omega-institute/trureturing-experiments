"""First common order of the reference rotations yielding the 719 ground state at lambda=51^-5."""
from itertools import product
from pathlib import Path
import numpy as np
from flint import arb,acb
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-selector-grid.py")))
first,hstar,xstar,lam,observe=[m[i] for i in ("first","hstar","xstar","lam","observe")]
def triples(n,c):
    out=[]
    for i,j in product(range(n),repeat=2):
        if (i,j)>((-i)%n,(-j)%n):
            continue
        a,b=(arb.pi()*i/n).sin_cos()
        t,v=(arb.pi()*j/n).sin_cos()
        out.append(((i,j),(a*a,t*t,1-(b*v+c*a*t)**2)))
    return out
def later(ref):
    def value(x,h):
        return sum((a-b)**2 for a,b in zip(x,ref))+lam*h*h
    win=value(xstar,hstar)
    start=int(acb.zeta_zero(1700).imag.floor().unique_fmpz())
    cap=int((win/lam).sqrt().upper().ceil().unique_fmpz())
    def scan(a,b):
        if value(observe(a.union(b))[1],a)>win:
            return None
        if b-a>arb(1)/8:
            mid=(a+b)/2
            return scan(a,mid) or scan(mid,b)
        no=a.zeta_nzeros().unique_fmpz()
        nh=b.zeta_nzeros().unique_fmpz()
        assert no is not None and nh is not None
        for idx in range(max(int(no)+1,1701),int(nh)+1):
            h=acb.zeta_zero(idx).imag
            vv=value(observe(h)[1],h)
            if vv<win:
                return idx
            assert vv>win
        return None
    for i in range(start,cap):
        if lam*i*i>win:
            break
        result=scan(arb(i),arb(i+1))
        if result:
            return result
    return None
def check(n):
    al,be=[triples(n,c) for c in m["source"]["cs"]]
    winners=[]
    pts=[(i+j,x+y) for i,x in al for j,y in be]
    for s in range(0,len(pts),256):
        batch=pts[s:s+256]
        arr=np.array([[float(a.mid()) for a in v] for _,v in batch])
        cheapest=(m["const"][:,None]-m["coef"]@arr.T).argmin(axis=0)
        for (angles,r),k in zip(batch,cheapest):
            if k!=718:
                a,b=first[int(k)]
                assert a<sum(x*y for x,y in zip(b,r))
            else:
                assert all(a>sum(x*y for x,y in zip(b,r)) for j,(a,b) in enumerate(first) if j!=718)
                if later(r) is None:
                    winners.append(angles)
    print("common order",n,"inequivalent targets",len(pts),"719 winners",winners)
    return winners
if __name__=="__main__":
    for n in range(1,16):
        result=check(n)
        if result:
            assert n==9 and result==[(2,1,4,6)]
            break
    else:
        raise AssertionError("no ninth-order reference found")
