"""Balance the first endpoint-dependent heat coefficient while preserving the fitted crossing."""
from pathlib import Path
import runpy
from types import FunctionType
from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-weyl-balanced-arcs.py"))
q,r,base,shell=m["q"],m["r"],m["b"],m["shell"]
source=base["source"]
V=shell["V"]
dL=[q(0) for _ in V]
dR=dL[:]
for i,j,d,_ in shell["edges"]:
    ell=(shell["rho"]*shell["theta"](d)).sqrt()
    f=m["fbox"][m["band"](d)]
    for k in (i,j):
        dL[k]+=ell
        dR[k]+=f*ell
# At each vertex D_i is a scalar times the internal 2-by-2 identity.
dhalf=sum(a*a-b*b for a,b in zip(dR,dL))
y=(-2*m["s0"].union(m["s1"])).exp()
mu=r**2*y
nu=-q("9/10")*mu+dhalf/20
assert min(V[9:])+nu>q(900000)
assert abs(dhalf/(20*mu))<q(".000001")

def gram():
    ids=source["edge_indices"]()
    added=[q(0) for _ in V]
    for i,j,d,_ in shell["edges"]:
        if (i,j) in ids:
            w=base["weights"][ids[i,j]]
            added[i]+=w
            added[j]+=w
    # Large-z coefficient of z^(-3/2) in the endpoint log-determinant difference.
    return 2*sum((a**3-b**3)/3-(a-b)*(pot+w)-a*(mu if i<9 else nu)
                 for i,(a,b,pot,w) in enumerate(zip(dR,dL,V,added)))

def refit():
    # Additional high onsite shift in static comparisons. Bound its effect with the derivative
    # on y >= 5.9, since the comparison before adding positive stiffness is diagonal.
    yy=min(V[9:])/r**2+nu.lower()/r**2
    assert yy>q("5.9")
    epsS=20*abs(nu/r**2).upper()*(-yy).exp()
    epsM=epsS*(yy-1)
    # The original proof already charges a full shear error to each graph.
    # Inject the small extra perturbation only in the two preconditioned face inequalities.
    fn=m["heat"]
    cfg=dict(base)
    cfg["errS"]+=epsS.upper()/2
    cfg["errM"]+=epsM.upper()/2
    FunctionType(fn.__code__,dict(fn.__globals__,b=cfg))()
    # For the gap estimate the minimum high potential has moved by nu, while the low
    # lower bound mu-10 holds after restricting to energies below 2000 M^2 pointwise.
    maxeta=q("1.018")**2
    beta=base["data"]["beta"]*maxeta
    E=r*r*(1+10*beta)/(1-r*r*beta)
    k=arb.pi()*(E*maxeta/base["data"]["d"]).sqrt()
    assert k<arb.pi() and min(V[9:])+nu-E-36*k/k.sin()>0
    assert min(V[9:])+nu>mu.upper()
    m["gauge"]()
    print("positive unbalanced heat coefficient",((18*mu-dhalf)/r**2).lower().str(12),
          ((18*mu-dhalf)/r**2).upper().str(12))
    print("required high/low shift",(-q("9/10")+dhalf/(20*mu)).str(8),"heat perturbation bounds",epsS.str(10),epsM.str(10))
    c3=gram()
    print("coefficient after compensation, dimensionless z^-3/2 interval",(c3/r**3).lower().str(10),
          (c3/r**3).upper().str(10))

if __name__=="__main__":
    refit()
