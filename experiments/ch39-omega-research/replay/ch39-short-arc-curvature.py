"""Sign of the induced Einstein coefficient on the shortest-edge zero-volume solutions."""
from pathlib import Path
import runpy
from flint import arb

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
grav = runpy.run_path(str(Path(__file__).with_name("ch39-relative-einstein-value.py")))
m,w,q,r = (short[k] for k in ("m","w","q","r"))
delta,T = grav["delta"],grav["T"]
ref,linear = grav["ref"],grav["linear"]
a,left_rem,rem,primitive = (grav[k] for k in ("a","left_rem","rem","primitive"))

def vertex(inc,shifts):
    value = q(0)
    for j,(g,d,p) in enumerate(zip(m["groups"],grav["length"],shifts)):
        p = p/r**2
        e = inc[j]/r
        f = d+e
        for i in g:
            x = a[i]
            if j<2:
                diff = primitive(x+p)-primitive(x)-linear*p
            else:
                diff = -p*((x.union(x+p)+delta)/(ref+delta)).log()
            value += 2*(diff+linear*e*(f+d)/2+rem(x+p,f))
        value -= 2*left_rem[j]
    return value

def error():
    c = 1-10/(delta*r*r)
    shell = m["shell"]
    raw = [(r*(shell["rho"]*shell["theta"](d)).sqrt(),short["fsbox"][short["band"](d)].lower())
           for _,_,d,_ in shell["edges"]]
    arc = -8*sum((T*l*f/2).tanh().log() for l,g in raw for f in (1,g))/(r*r*c)
    off = short["hs"]*(short["pa"]/(2*delta**2)+2*short["pd"]/(3*T**3))/c**3
    # On [0,delta], the right low static eigenvalues exceed mu-10. The left trial
    # has 18 eigenvalues below 200 and all its ordered eigenvalues grow by at most
    # (1+beta)x, using linear interpolation of the endpoint data on each arc.
    beta = 9*shell["rho"]*shell["theta"](804)
    left_top = max(m["V"])+short["norm"]+(1+beta)*delta*r*r
    lows = short["shiftbox"][0].lower()-10
    highs = min((a+b).lower() for a,b in zip(m["potential"][2:],short["shiftbox"][2:]))
    left_low = 200+(1+beta)*delta*r*r
    base = m["base"]
    assert max(z.upper() for z in w["source"]["enclose_static"](
        w["source"]["static"]([v.upper() for v in base["weights"]]))) < 200
    assert 18*(lows/left_low).log()+20*(highs/left_top).log()>0
    # 0 < log(det K_R/K_L) < 38 log(10/x) throughout this interval.
    low = 38*delta*((10/delta).log()+1)
    rad = (arc+off).upper()
    return (-rad).union(rad+low)

def assess(box,s,err):
    x = short["coords"](box,s)
    f,inc = short["lengths"](x),short["increments"](x)
    p = short["branch"](inc,s)
    arcs = -q(3).zeta()*sum(k*(1/b**2-1) for k,b in zip(short["inv2"],f))
    return arcs+vertex(inc,p)+err

def check():
    err = error()
    s0,s1 = w["s0"],w["s1"]
    initial = [[short["full"]]*3,s0.union(s1),0]
    pending = [initial]
    bounds = []
    while pending:
        box,s,depth = pending.pop()
        val = assess(box,s,err)
        if val>q(".05") and val<q("2"):
            bounds.append(val)
            continue
        if any(not short["margin"](i,1,box,s).contains(0) for i in range(3)):
            continue
        assert depth<24
        pos = depth%4
        ax = s if pos==3 else box[pos]
        mid = ax.mid()
        for seg in (ax.lower().union(mid),mid.union(ax.upper())):
            child = box[:]
            if pos<3:
                child[pos] = seg
            pending.append((child,seg if pos==3 else s,depth+1))
    assert bounds
    lo = min(v.lower() for v in bounds)
    hi = max(v.upper() for v in bounds)
    print("relative Einstein integral",lo.str(10),hi.str(10),"accepted boxes",len(bounds))
    return lo,hi

if __name__=="__main__":
    check()
