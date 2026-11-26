"""The exposed edge of the readout body gives physically realizable reference triples."""
from fractions import Fraction as F
from pathlib import Path
from flint import arb,acb
import runpy

base=runpy.run_path(str(Path(__file__).with_name("ch39-shell-earlier-zero-dominance.py")))
see=base["obs"]
h=acb.zeta_zero(719).imag
c2=F(8741547,12563572)
s0=see(h)[1]
penalty=arb(1)/51**5
def score(x,t,ref):
    return sum((a-b)**2 for a,b in zip(x,ref))+penalty*t*t
def verify(tval):
    r=[F(1,3),0,F(1,3),1,tval,1-c2*tval]
    ref=[arb(x.numerator)/x.denominator for x in r]
    e=score(s0,h,ref)
    stop=int((e/penalty).sqrt().upper().ceil().unique_fmpz())
    best,idx=e,719
    def scan(a,b):
        nonlocal best,idx
        if score(see(a.union(b))[1],a,ref)>best:
            return
        if b-a>arb(1)/8:
            mid=(a+b)/2
            scan(a,mid)
            scan(mid,b)
            return
        lo=a.zeta_nzeros().unique_fmpz()
        hi=b.zeta_nzeros().unique_fmpz()
        assert lo is not None and hi is not None
        for n in range(int(lo)+1,int(hi)+1):
            if n==idx:
                continue
            t=acb.zeta_zero(n).imag
            v=score(see(t)[1],t,ref)
            if v<best:
                best,idx=v,n
            else:
                assert v>best
    for i in range(stop):
        if penalty*i**2>best:
            break
        scan(arb(i),arb(i+1))
    return idx,best,e
for t in (F(2,3),F(1)):
    n,value,e=verify(t)
    print("compatible target t19",t,"global winning index",n,value.str(14),"719 score",e.str(14))
