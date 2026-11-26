"""A kinetic readout outside the two fitting conditions on the nine-to-ten shell branch."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/"ch39-fixed-inertia-scale.py"))
b,q,eta = m["base"],m["q"],m["eta"]
v = m["v0"].union(m["v1"])
u = m["u0"].union(m["u1"])
r0,r1 = m["r0"],m["r1"]

def compare(r,t,rows,z,y):
    scale = (400/r)**2
    gs = 18*(-t*y).exp()+20*(-t*z).exp()
    left = sum((-t*x/r**2).exp() for x in m["left_static"])
    # For t>=1 all omitted higher Dirichlet levels share the previously bounded tail.
    ymin = b["ymin"]*(400/r1)**2
    ratio = (-8*ymin).exp()
    tail = 2*270*ratio**8/(1-ratio)
    diff = sum(n*(-t*a).exp()*(-t*(1/eta-1)*a).expm1() for n,a in rows)
    err = 2*m["errS"]+tail+t*m["trace"]*scale
    return gs-left+diff+arb(0,err.upper())

def check():
    readings = {2:[],4:[]}
    for i in range(256):
        r = r0+(r1-r0)*q(i).union(q(i+1))/256
        rows = [(n,a*(400/r)**2) for n,a in b["drows"]]
        for j in range(4):
            vv = m["v0"]+(m["v1"]-m["v0"])*q(j).union(q(j+1))/4
            z,_ = m["chart"](r,vv,u)
            y = (-2*vv).exp()*(400/r)**2
            for t in readings:
                readings[t].append(compare(r,t,rows,z,y))
    for t,intervals in readings.items():
        lo = min(z.lower() for z in intervals)
        hi = max(z.upper() for z in intervals)
        print("u, left minus right bounds",t,(-hi).str(12),(-lo).str(12))
        assert hi < (-q(".02") if t==2 else -q(".35"))

if __name__=="__main__":
    check()
