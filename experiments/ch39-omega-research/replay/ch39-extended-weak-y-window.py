"""Uniform scalar obstruction including the extra five and ten Dirac fermions."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
BASE = run_path(str(ROOT/"ch39-matrix-all-gauge-ratios.py"))
a = BASE["a"]


def exclusion(low, ratio, initial_radial, coeff):
    pi = arb.pi()
    step = a("1/256")
    u = 1/(8*pi**2)
    floor = (1/pi).lower()
    g,h = BASE["gauge_intervals"](low,ratio)
    r = initial_radial.lower()
    small = a(".01")
    trace_max = 12*pi+3*small
    fmax = 4*pi*(a("12/5")+a("54/5"))
    aa,bb,ab = coeff
    assert 4*aa*bb>ab**2

    def dg(g,h,tr,extra):
        return u*(a("2/3")-u/2*(a("8192/15")/g+594/h-15*tr-extra))

    def dh(g,h,tr):
        return u*(TEN_ONE+u/2*(TEN_SELF/h-TEN_CROSS/g+5*tr))

    def constrained_interval(lo,hi):
        return max(floor,lo.lower()).union(hi.upper())

    for i in range(1280):
        if g.upper()<floor or h.upper()<floor:
            return i*step
        gl,gu = max(floor,g.lower()),g.upper()
        hl,hu = max(floor,h.lower()),h.upper()
        ga,gb = max(floor,gl-step),gu+step
        ha,hb = max(floor,hl-step),hu+step
        assert -1 < dg(ga,ha,0,0) and dg(gb,hb,trace_max,fmax) < 1
        assert -1 < dh(ga,hb,0) and dh(gb,ha,trace_max) < 1
        z,w = (1/gb).lower(),(1/hb).lower()
        box = r+(-step).union(step)
        c5,c10 = a("18/5"),a("99/20")
        base = aa*z*z+bb*w*w+ab*z*w-12*(c5*z+c10*w)*box
        tangent_z = 2*aa*z+ab*w-12*c5*box
        tangent_w = 2*bb*w+ab*z-12*c10*box
        quad = (base+min(arb(0),tangent_z.lower())*(1/ga-z)
                    +min(arb(0),tangent_w.lower())*(1/ha-w))
        ferm = 12*small*min(arb(0),box.lower())-a("120/101")*small**2
        speed = (u/2*(416*box*box+quad+ferm)).lower()
        assert 0 < speed < 1, (i,ratio,speed,g,h,box)
        r = (r+step*speed).lower()
        if r>8*pi/202:
            return (i+1)*step
        # Endpoint comparisons use the indicated monotonicity of both gauge equations.
        gl2 = gl+step*dg(ga,ha,0,0)
        gu2 = gu+step*dg(gb,hb,trace_max,fmax)
        hl2 = hl+step*dh(ga,hl+step,0)
        hu2 = hu+step*dh(gb,max(floor,hu-step),trace_max)
        if gu2.upper()<floor or hu2.upper()<floor:
            return (i+1)*step
        g,h = constrained_interval(gl2,gu2),constrained_interval(hl2,hu2)
    raise AssertionError("No exclusion within the certified range")


def verify():
    ctx.prec = 120
    low = a(".739").union(a(".741"))
    pi = arb.pi()
    assert a(".001")**2*((13*a(".01")+40*pi)*5/(8*pi**2)).exp()<a(".01")
    assert (3/(pi/low**2-1)).sqrt()>a(".78") and (pi/low**2-3).sqrt()<a("1.67")
    coeff = BASE["MATRIX"]["gauge_source"]()
    worst = max(exclusion(low,a(F(i,100)).union(a(F(i+1,100))),arb(0),coeff)
                for i in range(78,167))
    ext = run_path(str(ROOT/"ch39-selected-extended-flow.py"))
    v,(g,h,x,y,_,_) = ext["initial"]()
    assert low.contains(1/(g+3*h).sqrt())
    assert (4*pi*y/x).sqrt()<a(".001")
    fixed = exclusion(1/(g+3*h).sqrt(),a("8/5"),a("40/2424"),coeff)
    assert worst <= a("925/256") and worst.exp() < a("37.088")
    assert fixed <= a("774/256") and v*fixed.exp() < a("1.219e16")
    print("All 89 gauge-ratio boxes by log scale",worst.str(10),"factor",worst.exp().str(10))
    print("Fixed 624 scalar boundary: factor",fixed.exp().str(10),"scale / GeV",(v*fixed.exp()).str(11))


if __name__ == "__main__":
    verify()
