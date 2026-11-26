"""Integrated continuum relative curvature term for a common minimal coupling."""
from pathlib import Path
import runpy
from flint import arb

vac = runpy.run_path(str(Path(__file__).with_name("ch39-relative-vacuum-value.py")))
m,w,q,r = (vac[k] for k in ("m","w","q","r"))
delta = q(".0004")
T = delta.sqrt()
ref = q(6)
linear = -(ref+delta).log()-1
a = vac["a"]
length = vac["left_d"]

def primitive(a):
    return -(a+delta)*(a+delta).log()

def rem(a,d):
    _,hlog = vac["root_log"](a,d,T,delta)
    return d*hlog-(a+delta)*(T*d/(a+delta)).log1p()

left_rem = [sum(rem(a[i],d) for i in g) for g,d in zip(m["groups"],length)]

def estimate(u,v,s):
    shifts = [x/r**2 for x in m["branch"](u,v,s)]
    fs = w["factors"](u,v)
    value = -q(3).zeta()*sum(k*(1/f**2-1) for k,f in zip(vac["inv_length2"],fs))
    for j,(g,d,p) in enumerate(zip(m["groups"],length,shifts)):
        e = (m["du"][j]*(u-1)+m["dv"][j]*(v-1))/r
        f = d+e
        for i in g:
            x = a[i]
            if j<2:
                part = primitive(x+p)-primitive(x)-linear*p
            else:
                axis = x.union(x+p)
                part = -p*((axis+delta)/(ref+delta)).log()
            value += 2*(part + linear*e*(d+f)/2 + rem(x+p,f))
        value -= 2*left_rem[j]
    return value

def errors():
    c = 1-10/(delta*r*r)
    lens = [(r**2*m["shell"]["rho"]*m["shell"]["theta"](d)).sqrt()
            for _,_,d,_ in m["shell"]["edges"]]
    arc = -8*sum((T*l*f/2).tanh().log() for l in lens for f in (1,q(".979")))/(r*r*c)
    hs = 4*sum(m["base"]["weights"][idx]**2 for idx in m["source_ids"].values())/r**4
    p = m["branch"](*w["box"], w["s0"].union(w["s1"]))
    pa = max(abs(x).upper() for x in p)/r**2
    pd = max((abs(u*(w["box"][0]-1)+v*(w["box"][1]-1))/r).upper()
             for u,v in zip(m["du"],m["dv"]))
    off = hs*(pa/(2*delta**2)+2*pd/(3*T**3))/c**3
    # The bound K(x)<10 holds throughout [0,.01] by the vacuum-value estimate.
    low = 38*delta*((10/delta).log()+1)
    return (arc+off+low).upper()

def check():
    err = errors()
    values = [estimate(*p)+arb(0,err) for p in vac["matched_cells"]()]
    lo,hi = min(x.lower() for x in values),max(x.upper() for x in values)
    print("uniform error",err.str(12),"integrated curvature moment",lo.str(12),hi.str(12))
    assert lo>0

if __name__ == "__main__":
    check()
