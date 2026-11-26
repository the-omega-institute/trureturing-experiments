"""Certified breakdown scales for the rounded 11490,624 matching point."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb,arb_series,ctx

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
SELECT = run_path(str(ROOT/"ch39-rounded-phase-selection.py"))
MATCH,a = SELECT["MATCH"],SELECT["a"]


def initial():
    n,gain,sign,r,mass = next(SELECT["matches"]())
    assert (n,gain,sign)==(11490,624,1)
    ratio = MATCH["RATIO"]
    eps = list(ratio["high_zero_readings"]())[-1][2]
    s = MATCH["SCATTER"]["scalar_thresholds"](ratio["product_polynomials"](),r,a("1/12"))
    ds = MATCH["threshold_correction"](r*r,a("1/256"))
    k = a("8/5")
    low,v,lx = ratio["MODEL"]["matched"]([c+d for c,d in zip(s,ds)],eps,k,r)
    g,h = 1/(low*low*(1+3/k**2)),1/(low*low*(3+k*k))
    x = (2*lx).exp()
    return v,[g,h,x,x*eps**2,a("7/2")/g,a("13/10")/g]


def verify():
    ctx.prec = 150
    b = -F(10,3)+F(4,3)*(F(1,2)+F(3,2))
    diagonal = F(6796,15)+(2*F(12,5)+F(50,3))+3*(2*F(18,5)+F(50,3))
    assert (b,diagonal)==(F(-2,3),F(8192,15))
    v,state = initial()
    pi = arb.pi()
    u = 1/(8*pi**2)
    off = a("54/5"),a("297/20")
    def rhs(state):
        g,h,x,y,f,j = state
        z,w = 1/g,1/h
        guard = off[0]*z+off[1]*w
        return [u*(a("2/3")-(a(diagonal)*z+594*w-45*(x+y)-a("12/5")*f-a("54/5")*j)*u/2),
                u*(TEN_ONE+(TEN_SELF*w-TEN_CROSS*z+15*(x+y))*u/2),
                u*x*(13*x+10*y-guard),u*y*(10*x+13*y-guard),
                u*f*(a("16/5")*f+3*j-a("72/5")*z),
                u*j*(f+a("44/5")*j-a("108/5")*z)]
    def jet(state,order):
        series = [[s] for s in state]
        for n in range(order):
            slope = rhs([arb_series(seq,prec=n+1) for seq in series])
            for seq,item in zip(series,slope):
                seq.append(item[n]/(n+1))
        return series
    matrix = run_path(str(ROOT/"ch39-quartic-matrix-crossing.py"))
    aa,bb,ab = matrix["gauge_source"]()
    radial = a("242/2424").lower()
    step = a("1/1024")
    width = 3*step
    margin = (-width).union(width)
    gauge_crossings = {}
    scalar_end = None
    for i in range(6200):
        tube = [s+margin for s in state]
        slopes = rhs(tube)
        assert all(abs(d)<3 for d in slopes) and slopes[0]<0<slopes[1]
        if scalar_end is None:
            rbox = radial+(-step).union(step)
            z,w,x = 1/tube[0],1/tube[1],tube[2]
            bound = u/2*(416*rbox*rbox+(-12*(a("18/5")*z+a("99/20")*w)+12*x)*rbox
                         +aa*z*z+bb*w*w+ab*z*w-a("120/101")*x*x)
            speed = bound.lower()
            assert abs(speed)<1
            radial = (radial+step*speed).lower()
            if radial > 8*pi/202:
                scalar_end = v*((i+1)*step).exp()
        coeff,rem = jet(state,4),jet(tube,5)
        new = []
        for seq,tail in zip(coeff,rem):
            value = tail[-1]
            for c in reversed(seq):
                value = c+step*value
            new.append(value)
        for label,target in ((4,1/pi),(1,1/(4*pi))):
            if state[0].upper()>=target.lower() and new[0].lower()<=target.upper():
                window = (v*(i*step).exp()).union(v*((i+1)*step).exp())
                gauge_crossings[label] = window if label not in gauge_crossings else gauge_crossings[label].union(window)
        state = new
        if state[0] < 1/(4*pi):
            break
    assert scalar_end is not None and scalar_end < a("1.127e15") and scalar_end < 2*v
    assert a("6.12e16") < gauge_crossings[4] < a("6.15e16")
    assert a("1.72e17") < gauge_crossings[1] < a("1.74e17")
    print("X contact bound by",scalar_end.str(9),"GeV")
    print("Gauge alpha5 crossings",[(k,win.str(8)) for k,win in gauge_crossings.items()])


if __name__=="__main__":
    verify()
