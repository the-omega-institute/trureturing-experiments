"""Gauge regeneration of a zero portal mass trace before the second heavy threshold."""
from pathlib import Path
from runpy import run_path
from flint import arb

base=run_path(str(Path(__file__).with_name("ch39-shell-gw-inert-quartic-band.py")))
model,flow,q,loop=(base[k] for k in ("model","flow","q","loop"))


def check():
    _,m=flow["masses"](43)
    eps=q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b=((1+eps)/(1-eps)).sqrt()
    spread=(1/b).union(b)
    s0,s1=q(".36"),q(".391")
    c=[v*spread for v in m[:2]]
    first=model["shifted"](c[0],s0.union(s1))
    def ratio(s):
        return (model["shifted"](c[1],s)/model["shifted"](c[0],s)).log()
    x=ratio(s0).union(ratio(s1))/loop
    assert x>0
    y,_=base["matched"]["get_boundary"](first.lower(),first.upper())
    g=[q(5,3)*y[2]**2,y[1]**2,y[0]**2]
    trace=sum(v*v for i,v in enumerate(y) if i not in (0,1,2,6))
    cap=trace/(1-10*x*trace)
    z=[v+arb(0,".001") for v in g]
    for i in range(3):
        b=flow["beta"][i]+model["coeff"]["db"][i]
        B=[v+d for v,d in zip(model["coeff"]["base"][i],model["coeff"]["extra"][i])]
        f=2*z[i]**2*(b+(sum(a*g2 for a,g2 in zip(B,z))-q(0).union(max(model["coeff"]["yukawa"][i])*cap))/loop)
        assert abs(x*f)<q(".00099")
    W,Y=z[1],q(3,5)*z[0]
    assert 3*(3*W+Y)<4
    forcing=q(3,2)*(3*W**2+Y**2)
    assert forcing>q(".45")
    # If V remains bounded below, its two single-field quartics are nonnegative.
    # The positive source and a lower diagonal bound -4 imply p(x)>=c*x*exp(-4*x).
    pmin=forcing*x*(-4*x).exp()
    # Off-diagonal mass beta is 2*p. The normalized mass response matrix is
    # cooperative, with each diagonal >=-2. Hence r_01 >= c*x²*exp(-6*x).
    response=forcing*x**2*(-6*x).exp()
    induced=response*first**2
    print("log mass ratio",float((loop*x).lower()),float((loop*x).upper()),
          "portal trace lower",pmin.lower().str(8),"mass matrix response lower",response.lower().str(8),
          "normalized heavy mass contribution GeV² above",induced.lower().str(8))
    assert induced>q("1.5e10")


if __name__=="__main__":
    check()
