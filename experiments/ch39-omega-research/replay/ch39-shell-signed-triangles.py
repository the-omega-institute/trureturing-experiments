"""Exposed points of the three Wilson-deficit surface and the signs that support zero 719."""
from fractions import Fraction as F
from flint import arb,acb,ctx

ctx.prec=120
T=acb.zeta_zero(719).imag
for r,ds,c2 in ((17,(6,12,14),F(81,98)),(867,(244,724,804),F(8741547,12563572))):
    z=[(1-arb(d)/(2*r)).acos() for d in ds]
    p,q=[2*T*(x-z[0]) for x in z[1:]]
    s,x=p.sin_cos()
    t,y=q.sin_cos()
    c2=arb(c2.numerator)/c2.denominator
    k=2*c2.sqrt()/(1+c2)
    sigma=1 if s*t>0 else -1
    v=sigma*y
    b=k*v*abs(s)/abs(t)-x
    a=k*x*abs(t)/abs(s)-v
    margin=b*b-(1-k*k)
    assert margin>0
    delta=(1-c2)/(1+c2)
    weights=[2*(a-sigma*delta),2*sigma*(b-delta),4*sigma/(1+c2)]
    assert all(not w.contains(0) for w in weights)
    print("radius squared",r,"normalized deficit weights",*[w.str(15) for w in weights],
          "concavity margin",margin.str(15))
