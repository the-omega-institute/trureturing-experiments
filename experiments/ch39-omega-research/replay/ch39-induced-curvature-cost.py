"""Induced gravity, scalar vacuum term and the quadratic-curvature cost of matching Newton."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/"ch39-induced-newton.py"))
q,pi,base,chart = m["q"],m["pi"],m["base"],m["m"]
Q = 1000*chart["T0"].union(chart["T1"]).exp()
MP = q("1.220890e19")
v,u = chart["v0"].union(chart["v1"]),chart["u0"].union(chart["u1"])

def integrated():
    values = {2:[],3:[]}
    for i in range(12):
        r = m["r0"]+(m["r1"]-m["r0"])*q(i).union(q(i+1))/12
        for n in values:
            left,right = m["traces"](r,v,u,n)
            values[n].append(16*(left+right))
    return [min(x.lower() for x in vals).union(max(x.upper() for x in vals)) for vals in values.values()]

def before_breakdown():
    # For a Dirichlet subspace of codimension 38, the full count is at most 38 higher.
    # Bound E2 by integrating the Gaussian integral bound for every Dirichlet ladder.
    s=q("3.85")
    total=q(0)
    for eta in (q(1),m["eta"]):
        total+=38+sum(n*(3*pi*eta/d).sqrt()*m["r1"]*s/3 for n,d in base["data"]["factors"])
    bound=16*total*(s*Q/MP)**2/(12*pi)
    print("fraction through perturbative range",bound.upper().str(12))
    assert bound<q(".00003")

def verify():
    w2,w3 = integrated()
    minimal = w2/(12*pi)*(Q/MP)**2
    xi = (1-1/minimal)/6
    density = w3/(32*pi**2)
    curvature = (MP/Q)**4/(16*w3)
    equality_scale = MP/(16*pi*curvature).sqrt()
    print("universal xi for induced Newton",xi.lower().str(12),xi.upper().str(12))
    print("vacuum magnitude / Q^4",density.lower().str(12),density.upper().str(12))
    print("square curvature coefficient lower",curvature.lower().str(12))
    print("equal-term sqrt(|R|) upper GeV",equality_scale.upper().str(12))
    assert xi<q("-3e6") and curvature>q("1e13") and density>1
    before_breakdown()

if __name__=="__main__":
    verify()
