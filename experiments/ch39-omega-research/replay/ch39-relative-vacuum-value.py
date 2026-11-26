"""Certified continuum relative determinant integral on the thermally matched branch."""
from itertools import product
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-relative-vacuum.py")))
w, q, r = m["w"], m["q"], m["r"]
shell = m["shell"]
delta, T = q(".01"), q(".1")
a0 = q(6)

def prime(a):
    return a*(a+delta).log()+(a-delta)/2

slope = (a0+delta).log()+a0/(a0+delta)+q(".5")
linear = prime(a0)-a0*slope

def root_log(a,d,T,delta):
    # Primitive correction when the diagonal contains d sqrt(x).
    log = (delta+T*d+a).log()
    discr = d*d/4-a
    if discr > 0:
        h = discr.sqrt()
        hlog = sum(k*(T+k).log() for k in (d/2+h,d/2-h))
    else:
        assert discr < 0
        h = (-discr).sqrt()
        hlog = d*log/2-2*h*(h/(T+d/2)).atan()
    return log,hlog

def remainder(a,d):
    log,hlog = root_log(a,d,T,delta)
    return ((a*a-delta**2)*(T*d/(a+delta)).log1p()-d*d*a*log+(d**3-2*d*a)*hlog)/2

# Arc correction, common-endpoint off-diagonal correction, and low-x enclosure.
left_d = [d/r for d in m["length"]]
a = [(v+b)/r**2 for v,b in zip(m["V"],m["added"])]
left_rem = [sum(remainder(a[i],d) for i in group) for group,d in zip(m["groups"],left_d)]
inv_length2 = [sum(1/(r*r*shell["rho"]*shell["theta"](d))
                   for _,_,d,_ in shell["edges"] if w["band"](d)==j) for j in range(3)]
inv_length4 = [sum(1/(r*r*shell["rho"]*shell["theta"](d))**2
                   for _,_,d,_ in shell["edges"] if w["band"](d)==j) for j in range(3)]

def uniform_error():
    # Both the exact and asymptotic endpoint matrices dominate x-10/r^2:
    # all edge energies (including the full hyperbolic response) are nonnegative.
    c = 1-10/(delta*r*r)
    arc = 28*q(3).zeta()*sum(inv_length2)*(1+1/q(".979")**2)/(r*r*c)
    # The common off-diagonal B has Hilbert-Schmidt norm squared 4 sum w_e^2/r^4.
    hs = 4*sum(m["base"]["weights"][idx]**2 for idx in m["source_ids"].values())/r**4
    shifts = m["branch"](*w["box"], w["s0"].union(w["s1"]))
    pa = max(abs(x).upper() for x in shifts)/r**2
    pd = max((abs(u*(w["box"][0]-1)+v*(w["box"][1]-1))/r).upper()
             for u,v in zip(m["du"],m["dv"]))
    off = hs*(pa/delta+2*pd/T)/c**3
    # K(0) is positive from the previous static mass certificate and the right onsite gap.
    # y*coth(y/2) <= y+2 bounds the dynamic norm on 0<=x<=delta.
    eta = q("1.018")
    xmax = max((v+p).upper() for v,p in zip(m["potential"],shifts))/r**2
    emax = max((r*r*shell["rho"]*shell["theta"](d)*delta).sqrt().upper() for _,_,d,_ in shell["edges"])
    _,norm = m["base"]["arc_budget"]()
    assert xmax+delta+(norm+18*(eta*emax+2))/r**2 < 10
    low = 38*delta**2*((10/delta).log()+q(".5"))/2
    return (arc+off+low).upper()

def integral(u,v,s):
    p = [x/r**2 for x in m["branch"](u,v,s)]
    fs = w["factors"](u,v)
    result = -q("1.5")*q(5).zeta()*sum(k*(1/f**4-1) for k,f in zip(inv_length4,fs))
    for j,(g,d,pj) in enumerate(zip(m["groups"],left_d,p)):
        e = (m["du"][j]*(u-1)+m["dv"][j]*(v-1))/r
        f = d+e
        diff2 = e*(f+d)
        # All vertices of an orbit share d,p, but can have distinct fitted stiffnesses.
        for i in g:
            x = a[i]
            derivative = q(0)
            for k in range(4):
                b = x+pj*q(k).union(q(k+1))/4
                if j < 2:
                    derivative += (prime(b)-linear-slope*b)/4
                else:
                    axis = b.union(a0)+delta
                    diff = b-a0
                    derivative += diff*diff*(1/axis+delta/(axis*axis))/8
            # Restore the two subtracted polynomial moments using C2=C4=0.
            result += 2*(pj*derivative + linear*diff2/2
                        +slope*(diff2*(x-(f*f+d*d)/4)+f*f*pj)
                        +remainder(x+pj,f))
        result -= 2*left_rem[j]
    return result

def segment(lo,hi,k,n):
    return lo+(hi-lo)*q(k).union(q(k+1))/n

def matched_cells():
    s0,s1 = w["s0"],w["s1"]
    trace,_ = m["base"]["arc_budget"]()
    # Here only the continuum is needed. The inherited Dirichlet enclosures also cover it.
    es = (1+q("1.018")**2)*m["base"]["errS"]+w["tail"]+q("1.001")*trace
    em = (1+q("1.018")**2)*m["base"]["errM"]+w["tail"]+q("1.137")*trace
    n = 8
    params = [[segment(c-d,c+d,k,n) for k in range(n)] for c,d in zip(w["center"],w["radii"])]
    for u,v in product(*params):
        inner = w["dirichlet"](u,v)
        uc,vc = u.mid(),v.mid()
        mid = w["dirichlet"](uc,vc)
        jac = w["jac"](u,v)
        inner = [z.intersection(h+jac[i,0]*(u-uc)+jac[i,1]*(v-vc)) for i,(z,h) in enumerate(zip(inner,mid))]
        for k in range(n):
            s = segment(s0,s1,k,n)
            onsite = m["endpoints"](m["branch"](u,v,s))
            if any(abs(i+h)>bound for i,h,bound in zip(inner,onsite,(es,em))):
                continue
            yield u,v,s

def check():
    err = uniform_error()
    estimates = [integral(u,v,s)+arb(0,err) for u,v,s in matched_cells()]
    assert estimates
    lo,hi = min(x.lower() for x in estimates),max(x.upper() for x in estimates)
    print("uniform determinant error",err.str(12),"cells",len(estimates))
    print("integrated S_L-S_R",lo.str(12),hi.str(12))
    assert lo>q(".09") and hi<q("1.88")

if __name__ == "__main__":
    check()
