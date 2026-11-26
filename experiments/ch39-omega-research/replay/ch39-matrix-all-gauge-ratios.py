"""Uniform comparison for all gauge ratios and bounded bridge couplings."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
MATRIX = run_path(str(ROOT/"ch39-quartic-matrix-crossing.py"))
a, MATCH = MATRIX["a"], MATRIX["MATCH"]


def gauge_intervals(low, k):
    return (1/(low**2*(1+3/(k*k))), 1/(low**2*(3+k*k)))


def exclude(low, k, gauge_coefficients):
    pi = arb.pi()
    step = a("1/256")
    unit = 1/(8*pi**2)
    floor = (1/pi).lower()
    g,h = gauge_intervals(low,k)
    radial = arb(0)
    bmax = a(".01")  # upper bound for each squared Y Yukawa
    tmax = 12*pi+3*bmax
    cap = 8*pi/202
    aa,bb,ab = gauge_coefficients
    c5,c10 = a("18/5"),a("99/20")
    assert 4*aa*bb>ab**2
    diag5, diag10 = a("6796/15"),TEN_SELF
    def dg(g,h,tr):
        return unit*(a("10/3")-(diag5/g+594/h-15*tr)*unit/2)
    def dh(g,h,tr):
        return unit*(TEN_ONE+(diag10/h-TEN_CROSS/g+5*tr)*unit/2)
    for j in range(1280):
        if g.upper()<floor or h.upper()<floor:
            return j*step
        gl,gu = max(floor,g.lower()),g.upper()
        hl,hu = max(floor,h.lower()),h.upper()
        # Bootstrap a unit speed tube on the conditional gauge box.
        g0,g1 = max(floor,gl-step),gu+step
        h0,h1 = max(floor,hl-step),hu+step
        assert dg(g0,h0,0)>-1 and dg(g1,h1,tmax)<1
        assert dh(g0,h1,0)>-1 and dh(g1,h0,tmax)<1
        z,w = (1/g1).lower(),(1/h1).lower()
        zspan,wspan = 1/g0-z,1/h0-w
        r = radial+(-step).union(step)
        # A tangent plane minorizes the positive definite gauge quadratic.
        base = aa*z*z+bb*w*w+ab*z*w-12*(c5*z+c10*w)*r
        sz = 2*aa*z+ab*w-12*c5*r
        sw = 2*bb*w+ab*z-12*c10*r
        lower_gauge = base+min(arb(0),sz.lower())*zspan+min(arb(0),sw.lower())*wspan
        # The external-leg term is nonnegative at R>=0; allow negative tube values.
        fermion = min(arb(0),r.lower())*12*bmax - a("120/101")*bmax*bmax
        speed = (unit/2*(416*r*r+lower_gauge+fermion)).lower()
        assert 0 < speed < 1
        radial = (radial+step*speed).lower()
        if radial>cap:
            return (j+1)*step
        # G' is increasing in G,H. H' increases in G and decreases in H.
        # At each moving endpoint only a single step of self-displacement is needed.
        gl2 = gl+step*dg(max(floor,gl-step),h0,0)
        gu2 = gu+step*dg(gu+step,h1,tmax)
        hl2 = hl+step*dh(g0,hl+step,0)
        hu2 = hu+step*dh(g1,max(floor,hu-step),tmax)
        if gu2.upper()<floor or hu2.upper()<floor:
            return (j+1)*step
        g = max(floor,gl2.lower()).union(gu2.upper())
        h = max(floor,hl2.lower()).union(hu2.upper())
    raise AssertionError("No crossing by the end of the comparison")


def verify():
    ctx.prec = 120
    low = a(".739").union(a(".741"))
    for _,_,_,_,state in MATCH["matched_points"]():
        assert low.contains(1/(state[0]+3*state[1]).sqrt())
    pi = arb.pi()
    assert (3/(pi/low**2-1)).sqrt()>a(".78") and (pi/low**2-3).sqrt()<a("1.67")
    # Bootstrap Y_a <= .01 for t <= 5 from |y_a(v)| <= .001.
    assert a(".001")**2*((13*a(".01")+40*pi)*5/(8*pi**2)).exp() < a(".01")
    coefficients = MATRIX["gauge_source"]()
    times = [exclude(low,a(F(i,100)).union(a(F(i+1,100))),coefficients)
             for i in range(78,167)]
    worst = max(times)
    print("All 89 ratio intervals excluded by t",worst.str(9),"scale ratio",worst.exp().str(9))
    assert worst < 5


if __name__ == "__main__":
    verify()
