"""An entire positive-density region below the fitted clock energy at every nonnegative four-loop coupling."""
from pathlib import Path
import runpy
from flint import arb,acb,arb_mat

root = Path(__file__).resolve().parent
loop = runpy.run_path(str(root/"ch39-clock-wilson-length.py"))
m = loop["m"]
q,sh,short = m["q"],m["sh"],m["short"]
radius = q(".0025")
pi = arb.pi()
U0,a0 = m["rotations"](719)
Uc = {d:-loop["I"] for d in m["axes"]}

def end(side):
    raw,V = m["rows"](side)
    cycles = m["short_cycles"](raw,U0,Uc)
    ds = [4*abs(((pi-x)/2).sin()) for x in a0]
    # The nominal two-pi representatives may differ; the sine absolute is periodic.
    lo = q(".125")
    ratio = q(2)**(q(1)/8)
    out = q(0)
    extra = q(0)
    start = q(0)
    for j in range(8*33+1):
        hi = lo if j==0 else lo*ratio
        out += m["interval"](start,hi,raw,V,U0,Uc,ds,cycles)
        gd,weight,_,_ = m["bounds"](start,True,raw,V)
        K = arb_mat(19,19)
        for i,v in enumerate(gd):
            K[i,i]=v
        for (i,k),w in weight.items():
            K[i,k]=K[k,i]=-w
        G=K.inv()
        marks=2*sum(weight[i,k]*(G[i,k]-weight[i,k]/(gd[i]*gd[k])) for i,k,d,_,_,_ in raw if d in m["axes"])
        extra += (hi*hi-start*start)*radius*marks
        start=lo=hi
    degree = short["base"]["arc_budget"]()[1].upper()/2
    tail = 152*(2*degree)**3/lo
    return (out+arb(0,(extra+tail).upper()))/(32*pi**2*short["r"]**4)

def check():
    L,R = end("left"),end("right")
    neutral_upper = q("9.9583e-8")*2*radius**2
    total = 16*(R-L)+3*(neutral_upper-q("8.37729e-12"))
    assert total<0
    # All selected quadrilaterals have one edge of each clock type.
    # At (pi,pi) both matrices are -I, so t=2; the trace changes by <=4*radius.
    zmin = min(z.lower() for _,ds,z,word in loop["cycles"](4) if ds==(244,244,724,804))
    tmin = 2-4*radius
    deficit = 8-tmin**2-2*zmin*tmin
    assert deficit<q(".062") and deficit<q("5.41178")
    h = acb.zeta_zero(54542).imag
    for d in m["axes"]:
        a = (sh["theta"](d)-sh["theta"](244))*h
        assert a.cos()<0 and abs(a.sin())<radius.sin()
    print("uniform spectral energy difference upper",total.upper().str(12))
    print("four-edge action upper on square",deficit.upper().str(12))
    print("zero lower density of competitors",(radius**2/pi**2).lower().str(12))
if __name__=="__main__":
    check()
