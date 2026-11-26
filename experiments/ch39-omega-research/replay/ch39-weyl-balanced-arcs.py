"""Two low-scale heat moments and the same ultraviolet Weyl coefficient on the fitted shell."""
from pathlib import Path
from types import FunctionType
import runpy

from flint import arb,arb_mat

root=Path(__file__).resolve().parent
b=runpy.run_path(str(root/"ch39-right-arc-inertia.py"))
q,source=b["q"],b["source"]
shell=source["shell"]
s0,s1=b["s0"],b["s1"]
r=q(400)
classes=sorted(set(d for _,_,d,_ in shell["edges"]))
def band(d):
    return 0 if d<=72 else 1 if d<=244 else 2
length=[q(0) for _ in range(3)]
for i,j,d,_ in shell["edges"]:
    length[band(d)]+=shell["theta"](d).sqrt()
a0,a1=[l/length[2] for l in length[:2]]
def factors(u,v):
    return u,v,1-a0*(u-1)-a1*(v-1)
center=[q(".9965"),q(".9888")]
radii=[q(".009"),q(".008")]
box=[x+arb(0,dx.upper()) for x,dx in zip(center,radii)]
fbox=factors(*box)
assert all(q(".979")<f<q("1.018") for f in fbox)
# Reference and finite mesh intervals for the first seven harmonics, classified by shell edge.
rows=[]
for d in classes:
    mult=2*sum(ds==d for _,_,ds,_ in shell["edges"])
    E=b["data"]["d"]*shell["theta"](804)/shell["theta"](d)/r**2
    for k in range(1,8):
        y=E*k*k
        lo=y/(1+arb.pi()**2*k*k*b["h"]**2/6)
        hi=y/(1-k*k*b["h"]**2)**2
        rows.append((band(d),mult,lo.lower().union(hi.upper())))
# Uniform absolute tail in the two heat functions; evaluated directly, not differentiated.
coef=b["ymin"]/q("1.018")**2
tailratio=(-q(".75")*8*coef).exp()
assert 64*coef>20
tail=2*270*tailratio**8/(1-tailratio)

def dirichlet(u,v):
    fs=factors(u,v)
    S=M=q(0)
    for j,n,k in rows:
        inv=1/fs[j]**2
        g=(-(inv-1)*k)
        S+=n*(-k).exp()*g.expm1()
        M+=n*k*(-k).exp()*(inv*g.exp()-1)
    return [S,M]

def jac(u,v):
    fs=factors(u,v)
    ds=[q(0) for _ in range(3)]
    dm=ds[:]
    for j,n,k in rows:
        y=k/fs[j]**2
        val=2*n*y*(-y).exp()/fs[j]
        ds[j]+=val
        dm[j]+=val*(y-1)
    return arb_mat([[z[0]-a0*z[2],z[1]-a1*z[2]] for z in (ds,dm)])

def heat():
    # At fixed mass measure the equilibrium static matrix is unaffected by the length rescaling.
    # Left/right use identical static stiffness and potential, plus a common positive low-site shift.
    # Relative to the vertex potential alone each graph Laplacian is positive with this trace cap.
    trace,_=b["arc_budget"]()
    errors=[(q("1.018")**2+1)*v+2*trace+tail for v in (b["errS"],b["errM"])]
    assert min(v/r**2 for v in shell["V"])>q("-.00007")
    errors=[w+q(".001")*2*trace for w in errors]
    J=jac(*center)
    P=arb_mat([[J[i,j].mid() for j in range(2)] for i in range(2)]).inv()
    for i in range(2):
        margins=[]
        j=1-i
        for sign in (-1,1):
            for k in range(32):
                f=center[:]
                f[i]+=sign*radii[i]
                f[j]+=radii[j]*(q(k).union(q(k+1))/16-1)
                dc=dirichlet(*f)
                for part in range(12):
                    s=s0+(s1-s0)*q(part).union(q(part+1))/12
                    y=(-2*s).exp()
                    sr=mr=q(0)
                    for pot in shell["V"][:9]:
                        x=pot/r**2
                        sr+=2*(-x).exp()*(-y).expm1()
                        mr+=2*(-x).exp()*((x+y)*(-y).exp()-x)
                    R=arb_mat([[w+z+arb(0,e.upper())] for w,z,e in zip(dc,(sr,mr),errors)])
                    margins.append((sign*(P*R)[i,0]).lower())
        print("heat residual face margin",i,min(margins).str(12))
        assert min(margins)>0

def gauge():
    mu=r**2*(-2*s0.union(s1)).exp()
    _,norm=b["arc_budget"]()
    # Keep the original potentials. Nine low vertices have V>=-10.
    ceiling=(mu+norm).upper()
    d=b["data"]["d"]
    max_eta=q("1.018")**2
    beta=b["data"]["beta"]*max_eta
    k=arb.pi()*(ceiling*max_eta/d).sqrt()
    slope=18*k*(k/2).tan()/ceiling
    lower=(1-10/mu.lower())/((1+slope)*(1+beta*(ceiling+10)))
    upper=1+norm/mu.lower()
    shift=(-upper.log()/2).lower().union((-lower.log()/2).upper())
    E=r**2*(1+10*beta)/(1-r*r*beta)
    k=arb.pi()*(E*max_eta/d).sqrt()
    assert k<arb.pi() and min(shell["V"][9:])-E-36*k/k.sin()>0
    t=b["T0"].union(b["T1"])
    s=s0.union(s1)
    fn=b["model"]["integrate"]
    run=FunctionType(fn.__code__,dict(fn.__globals__,Ls=b["logs"],WL=b["wl"]),argdefs=fn.__defaults__)
    vals=[run(z,w+shift)[j] for z,w,j in
          ((b["T0"],s,0),(b["T1"],s,0),(t,s0,1),(t,s1,1))]
    assert vals[0]>0 and vals[1]<0 and vals[2]>0 and vals[3]<0
    print("gauge face margins",*(z.str(10) for z in (vals[0].lower(),-vals[1].upper(),vals[2].lower(),-vals[3].upper())))
    print("optical length factors",[(float(f.lower()),float(f.upper())) for f in fbox])

if __name__=="__main__":
    heat()
    gauge()
