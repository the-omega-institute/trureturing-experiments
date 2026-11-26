"""Dimensionless log moment and gravitational coefficients of the same relative spectrum."""
from pathlib import Path
import runpy
from flint import arb

grav = runpy.run_path(str(Path(__file__).with_name("ch39-relative-einstein-value.py")))
vac,m,w,q,r = (grav[k] for k in ("vac","m","w","q","r"))
source,base = w["source"],m["base"]
weights = base["weights"]
lo = source["enclose_static"](source["static"]([x.lower() for x in weights]))
hi = source["enclose_static"](source["static"]([x.upper() for x in weights]))
assert min(v.lower() for v in lo)>0
log_left = sum(a.lower().union(b.upper()).log() for a,b in zip(lo,hi))
_,norm = base["arc_budget"]()
edges = [sum(w["band"](d)==j for _,_,d,_ in m["shell"]["edges"]) for j in range(3)]

def zero_shift(u,v,s):
    fs = w["factors"](u,v)
    p = m["branch"](u,v,s)
    mu = p[0]
    value = 2*sum(n*f.log() for n,f in zip(edges,fs))
    value += 18*(mu.lower()-10).union((mu+norm).upper()).log()-log_left
    floors = []
    for n,onsite,h in zip(m["count"][2:],m["potential"][2:],p[2:]):
        value += 2*n*(h/onsite).log1p()
        floors.extend((onsite.lower(),(onsite+h).lower()))
    # Positive stiffness moves every ordered static eigenvalue by between 0 and norm.
    return value+arb(0,(20*(norm/min(floors)).log1p()).upper())

def check():
    spans = [[],[],[],[],[]]
    ev,eg = vac["uniform_error"](),grav["errors"]()
    for u,v,s in vac["matched_cells"]():
        volume = vac["integral"](u,v,s)+arb(0,ev)
        einstein = grav["estimate"](u,v,s)+arb(0,eg)
        log = zero_shift(u,v,s)
        assert einstein>0 and volume>0 and log>0
        # 16 real components per level, minimal coupling, Euclidean convention.
        rho = volume/(2*arb.pi()**2)
        invG = 4*einstein/(3*arb.pi())
        R2 = log/(144*arb.pi()**2)
        ratio = (rho.lower()/invG.upper()**2).union(rho.upper()/invG.lower()**2)
        for row,val in zip(spans,(log,rho,invG,R2,ratio)):
            row.append(val)
    for name,row in zip(("log moment","abs(relative volume)/Q^4","relative inverse G/Q^2",
                          "relative square coefficient","volume/(relative inverse G)^2"),spans):
        assert all(v.is_finite() and v>0 for v in row)
        print(name,min(v.lower() for v in row).str(12),max(v.upper() for v in row).str(12))

if __name__=="__main__":
    check()
