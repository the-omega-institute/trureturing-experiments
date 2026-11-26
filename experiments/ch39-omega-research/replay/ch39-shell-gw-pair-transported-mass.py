"""Transport the two signed finite mass steps to one scale using the two-doublet mass RGE."""
from pathlib import Path
from runpy import run_path
from flint import arb

folder=Path(__file__).resolve().parent
pair=run_path(str(folder/"ch39-shell-gw-paired-mass-obstruction.py"))
model,flow,q,loop=(pair[k] for k in ("model","flow","q","loop"))


def check():
    _,m=flow["masses"](43)
    eps=q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b=((1+eps)/(1-eps)).sqrt()
    masses=[v*(1/b).union(b) for v in m[:2]]
    s0,s1=q(".36"),q(".391")
    first=model["shifted"](masses[0],s0.union(s1))
    def interval(s):
        return (model["shifted"](masses[1],s)/model["shifted"](masses[0],s)).log()/loop
    dt=interval(s0).union(interval(s1))
    y,_=pair["base"]["matched"]["get_boundary"](first.lower(),first.upper())
    light=q(".01")
    scalar=4*arb.pi()/3
    trace=8*arb.pi()
    gap=2*(light*scalar).sqrt()
    total=sum(v*v for i,v in enumerate(y) if i not in (0,1,2,6))
    total/=1-10*total*dt
    # The preceding two-doublet gauge-band certificate bounds 3*(3W+Y)<4.
    # In particular the positive constant gauge quartic is below 2/3.
    gauge4=q(2,3)
    rate=24*light**2+12*total*light+gauge4+trace**2/2+2*(trace/2+gap)**2
    assert y[6]+dt*rate<light
    pmin=-4*(light*scalar).sqrt()
    max_a00=12*light+6*total+2
    max_a11=12*scalar+2
    norm=arb.max(max_a00,max_a11)+2*trace
    difference=(norm*dt).exp()-1
    assert difference<q(".001")
    # Duhamel formula for U01, with U11 positive on the entire short interval.
    u01=2*pmin*dt*(max_a00*dt).exp()*(1+difference)
    u00_error=dt*(max_a00*(1+difference)+2*trace*difference)
    result=q(".020")-trace*u00_error+loop*u01
    needed=result*first**2/loop
    assert result>q(".018") and needed>q("6.3e17")
    print("U00 absolute displacement below",u00_error.upper().str(8),"U01 above",u01.lower().str(8),
          "transported dimensionless signed lower bound",result.lower().str(8),
          "light mass-squared correction at the second threshold exceeds",needed.lower().str(8),"GeV²")


if __name__=="__main__":
    check()
