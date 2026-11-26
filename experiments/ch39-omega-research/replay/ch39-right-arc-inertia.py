"""Coupled gauge and two-moment solution from two onsite masses and one arc inertia on one shell."""
from pathlib import Path
from types import FunctionType
import runpy
from flint import arb

root = Path(__file__).resolve().parent
comparison = runpy.run_path(str(root/"ch39-common-arc-heat-obstruction.py"))
source = comparison["source"]
model = source["gauge"]
q,r = model["q"],model["r"]
data = model["data"]
errS,errM = comparison["heat_constants"]()
weights = source["source_bounds"]()
left = source["mass_bounds"]()
logs = [r.log()-x.log()/2 for x in left]
wl = sum(logs)
T0,T1 = q("26.86"),q("27.01")
s0,s1 = q("2.15"),q("2.38")
z0,z1 = q("1.06"),q("1.36")
u0,u1 = q("-.035"),q(".035")
eta0,eta1 = q(".84"),q(".95")
Sref = q("18.044")
h = q(".01")

# Finite-grid Dirichlet eigenvalue bounds also enclose the continuous graph.
drows = []
for mult,d in data["factors"]:
    for k in range(1,8):
        cont = k*k*d/r**2
        lo = cont/(1+arb.pi()**2*k*k*h**2/6)
        hi = cont/(1-k*k*h**2)**2
        assert lo>1
        drows.append((mult,lo.lower().union(hi.upper())))
# y >= k^2*(min d_e/r^2)/3; for y e^-y use e^(-.8y) when y>=26.
ymin = min(d.lower() for _,d in data["factors"])/r**2/3
assert 64*ymin>26 and (q(26)/5).exp()>26
tail_ratio = (-q(".8")*8*ymin).exp()
tail = 2*270*tail_ratio**8/(1-tail_ratio)

def dchange(eta,first):
    out = q(0)
    for n,y in drows:
        if first:
            out += n*y*(-y).exp()*((-(1/eta-1)*y).exp()/eta-1)
        else:
            out += n*(-y).exp()*((-(1/eta-1)*y).expm1())
    return out+arb(0,tail.upper())

def reference(s,z):
    y = (-2*s).exp()
    return Sref-18*(-y).exp()-20*(-z).exp()

def eta_range(s,z,u):
    target = reference(s,z)+u
    assert dchange(eta0,False)<target<dchange(eta1,False)
    def root_bound(lower):
        a,b = eta0,eta1
        for _ in range(26):
            mid = (a+b)/2
            difference = dchange(mid,False)-target
            if (difference<0 if lower else not difference>0):
                a=mid
            else:
                b=mid
        return a.lower() if lower else b.upper()
    return root_bound(True).union(root_bound(False))

def static_bounds():
    low = comparison["all_bounds"](source["static"]([x.lower() for x in weights]))
    high = comparison["all_bounds"](source["static"]([x.upper() for x in weights]))
    ys = [a.lower().union(b.upper())/r**2 for a,b in zip(low,high)]
    return sum((-v).exp() for v in ys),sum(v*(-v).exp() for v in ys)

def arc_budget():
    shell = source["shell"]
    ids = source["edge_indices"]()
    incident = [q(0)]*19
    for i,j,_,z in shell["edges"]:
        w = 1+weights[ids[i,j]].upper() if (i,j) in ids else q(1)
        incident[i]+=w
        incident[j]+=w
    # The static right graph Laplacian is PSD, with or without each added energy.
    trace = 2*sum(incident)/r**2
    norm = 2*max(incident)
    return trace,norm

def gauge_faces():
    span = s0.union(s1)
    mu = r**2*(-2*span).exp()
    _,norm = arc_budget()
    ceiling = (mu+norm).upper()
    k = arb.pi()*(ceiling*eta1/data["d"]).sqrt()
    slope = 18*k*(k/2).tan()/ceiling
    # Energy matrix dominates the onsite-minus-energy matrix up to this relative correction.
    beta = data["beta"]*eta1
    ratio_lo = 1/((1+slope)*(1+beta*(ceiling+10)))
    ratio_hi = 1+norm/mu.lower()
    dlog = (-ratio_hi.log()/2).lower().union((-ratio_lo.log()/2).upper())
    # No nineteenth level up to the cutoff, including the finite-grid mass correction.
    energy = r**2*(1+10*beta)/(1-r**2*beta)
    k = arb.pi()*(energy*eta1/data["d"]).sqrt()
    assert k<arb.pi() and z0*r**2-energy-36*k/k.sin()>0
    Y = model["integrate"]
    integrate = FunctionType(Y.__code__,dict(Y.__globals__,Ls=logs,WL=wl),argdefs=Y.__defaults__)
    T = T0.union(T1)
    faces = [integrate(t,s+dlog)[i] for t,s,i in ((T0,span,0),(T1,span,0),(T,s0,1),(T,s1,1))]
    alpha = integrate(T,span+dlog)[2]
    print("right low log shift",dlog.lower().str(9),dlog.upper().str(9))
    print("gauge faces",[(float(v.lower()),float(v.upper())) for v in faces])
    print("inverse alpha",alpha.lower().str(10),alpha.upper().str(10))
    assert faces[0]>0 and faces[1]<0 and faces[2]>0 and faces[3]<0

def heat_faces():
    s,z,u = s0.union(s1),z0.union(z1),u0.union(u1)
    eta = eta_range(s,z,u)
    SL,ML = static_bounds()
    trace,_ = arc_budget()
    # Monotonicity and |f'|<=1 for f=e^-x,x e^-x bound the onsite/Laplacian difference.
    s_res = Sref-SL+arb(0,(trace+2*errS).upper())
    assert s_res+u0<0 and s_res+u1>0
    y = (-2*s).exp()
    for z,sign in ((z0,-1),(z1,1)):
        band = eta_range(s,z,u)
        # The first Dirichlet moment increases with eta when all y/eta>1.
        diff = dchange(band.lower(),True).lower().union(dchange(band.upper(),True).upper())
        value = 18*y*(-y).exp()+20*z*(-z).exp()+diff-ML+arb(0,(trace+2*errM).upper())
        print("first moment face",z.str(6),value.lower().str(10),value.upper().str(10))
        assert value*sign>0
    print("inertia ratio",eta.lower().str(12),eta.upper().str(12))
    print("value face margin",min(-(s_res+u0).upper(),(s_res+u1).lower()).str(10))

if __name__=="__main__":
    gauge_faces()
    heat_faces()
