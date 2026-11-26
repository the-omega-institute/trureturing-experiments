"""Independent truncation on the zero side at height 3200 and a xi-resolvent tail bound."""
from pathlib import Path
import runpy
from flint import arb, acb, ctx

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-local-prime.py")))
ctx.prec=72
pi, a = arb.pi(), m["a"]
A=arb(3200)
T=acb.zeta_zero(719).imag
n=A.zeta_nzeros().unique_fmpz()
assert n is not None and 1000<int(n)<3000
# xi'/xi at real s>1: zeta'/zeta is negative, so dropping it gives an upper bound.
s=A+arb(1)/2
logder=1/s+1/(s-1)-pi.log()/2+(s/2).digamma()/2
tail=(2*A/(a**4*(A-T)**4)*logder).upper()
def clock(r,d):
    active=(1-arb(d)/(2*r)).acos()
    ref=(1-arb(6 if r==17 else 244)/(2*r)).acos()
    return active-ref
p,q=clock(17,14),clock(867,804)
sums=[arb(0) for _ in range(4)]
for lo in range(1,int(n)+1,80):
    for z in acb.zeta_zeros(lo,min(80,int(n)-lo+1)):
        t=z.imag
        assert t<A
        wt=(a*(t-T)).sinc()**4+(a*(t+T)).sinc()**4
        x,y=(p*t).sin()**2,(q*t).sin()**2
        for i,value in enumerate((1,x,y,x*y)):
            sums[i]+=wt*value
weight,*values=[v+arb(0).union(tail) for v in sums]
x,y,xy=[v/weight for v in values]
direct=xy-x*y
explicit=m["denominator"](T)
assert weight.contains(explicit)
assert direct>arb(".01") and direct<arb(".014")
print("positive zeros, rigorous tail bound on total weight",n,tail.str(12))
print("direct zero-side covariance bounds",direct.lower().str(12),direct.upper().str(12))
