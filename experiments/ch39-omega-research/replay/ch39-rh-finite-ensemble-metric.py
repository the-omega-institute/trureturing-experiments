"""All 8^24 fields contribute to a finite two-arc effective curvature."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

directory=Path(__file__).parent
arc=run_path(str(directory/"ch39-rh-internal-arc-calibration.py"))
fixed=arc["fixed"]


def main():
    ctx.prec,ctx.cap=1536,3
    z,weights,_=fixed["shell"]["quad"](21)
    def response(t):
        return sum(w/(1+t*n) for w,n in zip(weights,z))
    def root(target):
        lo,hi=arb(0),arb(256)
        for _ in range(164):
            mid=(lo+hi)/2
            gap=response(mid)-target
            if gap>0:lo=mid
            elif gap<0:hi=mid
            else:break
        return lo.union(hi)
    a,b=arb(33)/34,arb(14)/17
    t1,t2=root(a),root(b)
    goal=t2/(3*t1)
    share=(33*goal-89)/64
    h=(144*share/(1-share)-3).root(3)
    rate=(2/(1+h+((1+h)**2-4).sqrt())).root(32)
    den=2**115
    xnum=int((rate*den).mid().floor().unique_fmpz())
    tnum=int((t1*den).mid().floor().unique_fmpz())
    x,t=arb(xnum)/den,arb(tnum)/den
    assert 0<x<1
    u=x**32
    c=8*((1+u**3)**3+3*u**3*(1+u)**3)
    triad=1152*u**3*(1+u)**3
    r=Q(1,256)
    tail=442368*r**16+884736*r**20+(8**24-9472-442368-884736)*r**24
    p_out=arb(tail.numerator)/tail.denominator/(c+triad)
    assert p_out<arb(2)**(-112)
    # The quadratic coefficient for any unit vector is bounded in absolute value by 16.
    err=32*p_out
    sign=(16*c+arb(80)/9*triad)/(c+triad)
    color=arb(8)/3
    uniform=arb(-1).union(arb(1))
    d=(1+sign+uniform*err)/(1+color+uniform*err)
    assert 1+color>err and abs(response(t)-a)<arb(2)**(-96)
    assert abs(response(3*t*d)-b)<arb(2)**(-96)
    print("Numerators of (x,t) over 2^115:",xnum,tnum)
    print("Uniform mass outside the exact zero sectors below",p_out.str(10))
    print("Energy-curvature ratio",d.str(16),"two records within 2^-96.")
    # Neither the 7 shell readings nor any one of the 42 arithmetic records
    # is the required individual nearest-step weight.
    readings=[response(arb(j*j)) for j in range(42)]
    closest=min(range(42),key=lambda j:float(abs(readings[j]-rate).mid()))
    assert all(abs(v-rate)>arb("0.002") for v in readings)
    print("Closest integer offset",closest,"readout",readings[closest].str(12),"target",rate.str(12))


if __name__=="__main__":
    main()
