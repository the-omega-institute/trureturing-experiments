"""Higgs matching using the same 719th-zero quark and two graph-generated scalar towers."""
from pathlib import Path
from types import FunctionType
import runpy
from flint import arb

root = Path(__file__).resolve().parent
fixed = runpy.run_path(str(root/"ch39-fixed-inertia-scale.py"))
base = fixed["base"]
profile = fixed["profile"]
q = fixed["q"]
T0,T1 = fixed["T0"],fixed["T1"]
T = T0.union(T1)
R = fixed["r0"].union(fixed["r1"])
V = fixed["v0"].union(fixed["v1"])
D = profile["old"]["D"]

def light():
    # Fit uncertainty across the physical clock and metric square, including the
    # down-type product of a within-sector ratio and the bottom/top normalization.
    env = {}
    source = (root/"ch39-ew.py").read_text().split("interval=(q(1000)")[0]
    exec(compile(source,"light_SM","exec"),env)
    gs = fixed["inputs"]["G0"]
    ys = []
    for value,rel in zip(profile["gauge"]["yc"][:6],
                                   (".002",".002","0",".002132",".002132",".000131")):
        ys.append(2*(q(value)*(1+arb(0,rel))).log())
    ys.extend(fixed["inputs"]["y0"][6:])
    initial = gs+ys
    dt = (q(1000)/q("91.1876")).log()
    evolve = env["evolve"]
    at_z = evolve(initial,-dt)
    lam0 = q(".5573").union(q(".5611"))/4
    vev0 = q("248.372").union(q("248.436"))
    at_tev = evolve(at_z+[lam0,(lam0*vev0**2).log()],dt)
    evolved = evolve(initial+at_tev[12:],T0)
    final = evolve(evolved,T1-T0,True)
    return final[12],final[13].exp()

def spectrum():
    mu = 400**2*(-2*V).exp()
    top = (mu+fixed["norm"]).upper()
    d = base["data"]["d"]
    eta = fixed["eta"]
    k = arb.pi()*(eta*top/d).sqrt()
    slope = 18*k*(k/2).tan()/top
    beta = base["data"]["beta"]*eta
    xright = (mu.lower()/((1+slope)*(1+beta*(top+10)))).lower().union(top)
    yl = [x/R**2 for x in base["left"]]
    yr = xright/R**2
    wl = sum(-y.log()/2 for y in yl)
    s = -yr.log()/2
    calc = profile["integrate"]
    logrs = [(R/x.sqrt()).log() for x in base["left"]]
    calc = FunctionType(calc.__code__,dict(calc.__globals__,g=fixed["inputs"],Ls=logrs,WL=wl),
                        argdefs=calc.__defaults__)
    alpha = calc(T,s)[2]
    return yl,yr,wl,s,alpha

def match():
    yl,yr,wl,s,alpha = spectrum()
    ll,msq = light()
    Q = 1000*T.exp()
    gsq = 4*arb.pi()/alpha
    lh = q(".02335")
    f = (2*arb.pi()**2*(lh-ll)/(3*wl)).sqrt()
    hard = q("36/5")*gsq**2*13*(-q(2)).exp()
    def poly(l):
        return 3*l*l/2+2*l+q("7/2")
    for y in yl:
        l=y.log()
        hard+=y*(q("78/25")*poly(l)*gsq**2+q("144/5")*(3-l)*gsq*f+24*((l-2)**2+1)*f**2)
    lr = yr.log()
    hard+=18*yr*q("12/25")*poly(lr)*gsq**2
    response = 8*18*yr*(1-lr)
    b = (hard+msq*(D/Q)**2)/(D*response)
    margin = 4*lh*q(".06")-f**2
    contact = 218*q(".06")+(216*f**2+72*(3*f)**2+288*b**2).sqrt()
    sensitivity = 1+hard*(Q/D)**2/msq
    ratio = b/f
    omitted = ratio.lower()**2*18*s.lower()/(3*wl.upper())
    omitted = omitted.union(ratio.upper()**2*18*s.upper()/(3*wl.lower()))
    assert 0<f and 0<b and margin>0 and contact<8*arb.pi()
    for name,x in (("quartic low",ll),("mass low GeV^2",msq),("inverse alpha",alpha),
                   ("left portal",f),("right portal",b),("stability",margin),("contact",contact),
                   ("relative mass sensitivity",sensitivity),("right quartic squared ratio",omitted)):
        print(name,x.lower().str(10),x.upper().str(10))

if __name__=="__main__":
    match()
