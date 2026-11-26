"""Enclose a nonlinear two-loop crossing for every phase in two arithmetic windows."""
from pathlib import Path
from runpy import run_path
from flint import arb

folder = Path(__file__).resolve().parent
model = run_path(str(folder/"ch39-shell-gw-nonlinear-crossing.py"))
density = run_path(str(folder/"ch39-shell-gw-positive-density-cp.py"))
flow,q = model["flow"], model["q"]


def check(n,width):
    Q,m = flow["masses"](n)
    raw = flow["full"]["model"]
    if n==43:
        s,f,_,_ = raw["m"]["transported_spectra"]()
    else:
        s,f,_ = raw["scalar_models"]()
    gaps = flow["full"]["all_levels"](s,f)
    delta = 4*width
    assert model["count"](model["rows"](f),q(".04"))==model["count"](model["rows"](f),q("1.05"))==2
    assert q("1.05")-q(".04")-delta>1
    assert gaps[1]+delta < 1 < gaps[2]-delta and gaps[-1]+delta < 7
    eps = delta/gaps[0]
    factor = ((1+eps)/(1-eps)).sqrt()
    stretch = (1/factor).union(factor)
    mass = [v*stretch for v in m]
    s0,s1,t0,t1 = map(q,(".311",".338","-.114","-.094"))
    s,t = s0.union(s1),t0.union(t1)
    def light(mass,ss):
        x = (model["shifted"](mass,ss)/1000).log()/model["loop"]
        return model["low"]["light"](x)
    varied = light(mass[0],s)
    center = light(m[0],(s0+s1)/2)
    e = density["window"](n,width)
    trace = model["continuous_light"](n,varied[0],center[0],e,3*delta)
    def face(ss,tt):
        return model["run_face"](Q,mass,ss,tt,(s0,s1),varied,trace)
    assert face(s,t0)[0]>0 and face(s,t1)[0]<0
    assert face(s0,t)[1]<0 and face(s1,t)[1]>0
    _,_,inverse,cap = face(s,t)
    assert q("38.86") < inverse < q("38.94") and cap < q("2.43")
    print("zero",n,"square half-width",width.str(2),"inverse coupling",inverse.str(6),
          "Yukawa squared norm upper",cap.upper().str(7))


if __name__=="__main__":
    for n,w in ((43,"3e-8"),(195,"1e-8")):
        check(n,arb(w))
