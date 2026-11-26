"""Quantitative inertia and energy requirements from two thermal moments and the gauge exposure."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/"ch39-scalar-overlap-bound.py"))
q = m["q"]
r = q(400)
e = q(1).exp()

def heat_coeff(a):
    # f(y)=(1-a*y)e^-y.  On 0<=y<=1 bound its slope as a function of -log(y)/2.
    if a>0:
        minimum = -a*(-1-1/a).exp()
        y = 2*(1+a)/(1+3*a+((1+3*a)**2-4*a*(1+a)).sqrt())
    else:
        minimum,y=q(0),q(1)
    c = (2*y*(1+a-a*y)*(-y).exp()).upper()
    f1 = (1-a)/e
    return 1-f1,c,f1-minimum,-minimum

def inertia_table(low):
    T,W = m["ranges"](low)
    log_bounds = [(r.log()-v.log()/2).lower() for v in low]
    assert sum(log_bounds)>W and T>6
    test = []
    for i in range(101):
        a = q(i)/20
        h,c,bminus,bplus = heat_coeff(a)
        assert c*6>h and bminus>0
        f1 = 1-h
        hleft = sum((1-a*x/r**2)*(-(x/r**2)).exp()-f1 for x in low)
        test.append((h,c,bminus,bplus,hleft))
    rows = []
    for p in range(9,19):
        n = 18-p
        required = sum(log_bounds[p:])
        spare = W.upper()-required
        assert spare>0
        qmin = 1
        for h,c,bminus,bplus,hleft in test:
            need = (hleft-n*h-c*spare-p*bplus)/bminus
            for candidate in range(qmin,39):
                if need>candidate:
                    qmin = candidate+1
                else:
                    break
        # Of the top p of the first 18 right eigenvalues, small exposures force a large
        # sum of masses. For eigenvalues below the input floor, clip y at exp(-2T).
        upper_left_sum = sum(x.upper() for x in low)
        trace_min = p*r*r*((-2*spare/p).exp()-(-2*T.lower()).exp())-upper_left_sum
        trace_min = max(trace_min.lower(),q(0))
        # The global Jensen bound on all 18 is independent of p.
        global_trace = (18*r*r*((-2*W.upper()/18).exp()-(-2*T.lower()).exp())-upper_left_sum).lower()
        rows.append((p,qmin,max(global_trace,trace_min)))
    print("positive index; necessary negative index; positive trace >")
    for p,k,x in rows:
        print(p,k,int(float(x))-1)
    assert rows[0][0]+rows[0][1]>38
    return rows

if __name__=="__main__":
    bare = inertia_table(m["profile"]["low"])
    shared = inertia_table(m["source"]["mass_bounds"]())
