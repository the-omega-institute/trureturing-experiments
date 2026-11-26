"""No perturbative gauge/Yukawa trajectory reaches the given finite horizon."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import Matrix, eye, ones

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
BASE = run_path(str(ROOT/"ch39-matrix-all-gauge-ratios.py"))
a = BASE["a"]


def check_trace_form():
    diag = 10*eye(3)+ones(3)
    cross = Matrix(3,3,lambda i,j: 5*int(j==i or j==(i+1)%3))
    form = diag.row_join(cross).col_join(cross.T.row_join(diag))
    assert form.eigenvals()=={23:1,3:1,5:2,15:2} and form*ones(6,1)==23*ones(6,1)


def trace_caps(horizon,step,pi,u):
    count = int((horizon/step).unique_fmpz())
    inv = (1/(24*pi)).lower()
    caps = []
    for i in reversed(range(count)):
        C = 54*pi/5+a("297/20")/(1/pi+a(".36")*i*step)
        C = C.upper()
        decay = (-u*C*step).exp()
        prev = (inv*decay+a("23/6")/C*(1-decay)).lower()
        caps.append(min(24*pi,1/min(prev,inv)).upper())
        inv = prev
    return caps[::-1]


def check_horizon(low,ratio,step,caps,pi,u):
    g,h = (x.upper() for x in BASE["gauge_intervals"](low,ratio))
    floor = (1/pi).lower()
    extra = a("66/5")*4*pi
    for i,tr in enumerate(caps):
        if min(g,h)<floor:
            return i
        # Each upper comparison endpoint can move at most one step in the tube.
        dg = u*(a("2/3")-u/2*(a("8192/15")/(g+step)+594/(h+step)-15*tr-extra))
        dh = u*(TEN_ONE+u/2*(TEN_SELF/max(floor,h-step)-TEN_CROSS/(g+step)+5*tr))
        g,h = (g+step*dg).upper(),(h+step*dh).upper()
    assert min(g,h)<floor, (low,ratio,g,h)
    return len(caps)


def verify():
    ctx.prec = 128
    check_trace_form()
    pi = arb.pi()
    u = 1/(8*pi*pi)
    step = a("1/256")
    assert u*(TEN_ONE-u*TEN_CROSS*pi/2)>a(".36")
    assert u*(a("2/3")+u/2*((a("8192/15")+594)*pi+15*24*pi+a("66/5")*4*pi))<1
    assert u*(TEN_ONE+u/2*((TEN_SELF+TEN_CROSS)*pi+5*24*pi))<1
    low = a(".739").union(a(".741"))
    assert (3/(pi/low**2-1)).sqrt()>a(".78") and (pi/low**2-3).sqrt()<a("1.67")
    horizon = a("17/2")
    caps = trace_caps(horizon,step,pi,u)
    worst = max(check_horizon(low,a(F(i,100)).union(a(F(i+1,100))),step,caps,pi,u)
                for i in range(78,167))
    ext = run_path(str(ROOT/"ch39-selected-extended-flow.py"))
    v,(g,h,*_) = ext["initial"]()
    selected = a("15/2")
    fixed = check_horizon(1/(g+3*h).sqrt(),a("8/5"),step,
                          trace_caps(selected,step,pi,u),pi,u)
    print("All gauge ratios excluded for horizon",horizon.str(6),"conditional step",worst)
    print("Selected horizon",selected.str(6),"conditional step",fixed,
          "scale /GeV",(v*selected.exp()).str(11))


if __name__ == "__main__":
    verify()
