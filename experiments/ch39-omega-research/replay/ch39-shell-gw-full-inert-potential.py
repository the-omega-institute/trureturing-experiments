"""Positive complete scalar potential to the gauge crossing on the 43rd-zero window."""
from collections import Counter
from itertools import combinations_with_replacement
from math import ceil
from pathlib import Path
from runpy import run_path
from flint import arb

folder = Path(__file__).resolve().parent
inert = run_path(str(folder/"ch39-shell-gw-inert-quartic-band.py"))
matched = inert["matched"]
model,flow,q,loop = (matched[k] for k in ("model","flow","q","loop"))


def pair(i,j):
    return (min(i,j),max(i,j))


def rates(groups,state,g,ytot):
    W,Y=g[1],q(3,5)*g[0]
    damping=3*(3*W+Y)
    trace=q(0).union(3*ytot)
    trace4=q(0).union(3*ytot**2)
    # Every member of a class has its own complex rephasing. Permutations
    # inside a class relate coefficients, also for inequivalent heavy masses.
    rhs={}
    for i,num in groups.items():
        l=state["l",i]
        value=24*l*l+q(3,8)*(2*W**2+(W+Y)**2)-damping*l
        if i==0:
            value+=4*trace*l-2*trace4
        for j,count in groups.items():
            multiplicity=count-int(j==i)
            if multiplicity:
                ij=pair(i,j)
                a,b=state["a",ij],state["b",ij]
                value+=multiplicity*(2*a*a+2*a*b+b*b)
        rhs["l",i]=value
    for i,j in combinations_with_replacement(groups,2):
        if i==j and groups[i]==1:
            continue
        ij=pair(i,j)
        a,b=state["a",ij],state["b",ij]
        l=state["l",i]+state["l",j]
        f=(2*trace if i==0 else 0)-damping
        da=4*l*(3*a+b)+4*a*a+2*b*b+q(3,4)*(3*W**2+Y**2-2*W*Y)+f*a
        db=4*l*b+8*a*b+4*b*b+3*W*Y+f*b
        for h,count in groups.items():
            multiplicity=count-int(h==i)-int(h==j)
            if multiplicity:
                ih,jh=pair(i,h),pair(j,h)
                ai,aj=state["a",ih],state["a",jh]
                bi,bj=state["b",ih],state["b",jh]
                da+=multiplicity*(4*ai*aj+2*ai*bj+2*bi*aj)
                db+=2*multiplicity*bi*bj
        rhs["a",ij]=da
        rhs["b",ij]=db
    return rhs


def activate(k,groups,state):
    g=1 if k<=2 else 2
    groups[g]+=1
    if ("l",g) not in state:
        state["l",g]=q(".10" if g==1 else ".02")
    for i in groups:
        if i==g and groups[g]==1:
            continue
        ij=pair(i,g)
        if ("a",ij) not in state:
            key=(i,g)
            state["a",ij]=q(".60" if key==(0,1) else ".10" if key==(1,1) else ".01")
            state["b",ij]=q(".01" if g==1 else ".002")


def check(heavy_gain=arb(0), threshold_steps=None, mass_path=None, bulk_portal=None):
    Q,m=flow["masses"](43)
    s0,s1,t0,t1=map(q,(".36",".391","-.139","-.112"))
    eps=q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b=((1+eps)/(1-eps)).sqrt()
    m=[v*(1/b).union(b) for v in m]
    first=model["shifted"](m[0],s0.union(s1))
    ys,_=matched["get_boundary"](first.lower(),first.upper())
    state={("l",0):ys[6]}
    gauge=[q(5,3)*ys[2]**2,ys[1]**2,ys[0]**2]
    first_total=sum(v*v for i,v in enumerate(ys) if i not in (0,1,2,6))
    F,D=flow["fit"]["readouts"](43)
    raw=sum((flow["fit"]["hs"](A,A).real.sqrt()+3*q(".000062"))**2 for A in (F,D))*q("7.9")**2
    # Combined central-vertex perturbation over the 43rd-zero square is below 0.000062.
    if heavy_gain>0:
        assert matched["density"]["window"](43,q("3e-8"))<q(".000062")
    ytot=first_total
    Vmargin=q(1)
    Hnorm=q(0)
    classes=Counter({0:1})
    margin=q(".001")
    largest=q(0)
    scalar_norm=q(0)
    for k in range(1,54):
        activate(k,classes,state)
        if k==3 and bulk_portal is not None:
            state["a",(0,2)]=bulk_portal
        if mass_path is not None:
            mass=model["shifted"](m[k-1],s0.union(s1))
            mass_path.append(("threshold",k,dict(classes),state.copy(),mass))
        if threshold_steps is not None:
            kind=1 if k<=2 else 2
            portal=2*state["a",(0,kind)]+state["b",(0,kind)]
            mass=model["shifted"](m[k-1],s0.union(s1))
            threshold_steps.append((mass,portal))
        if k==3:
            ytot+=heavy_gain**2*raw
        if k==53:
            duration=(Q*t0.union(t1).exp()/model["shifted"](m[-1],s0.union(s1))).log()/loop
        else:
            def log_ratio(s):
                return (model["shifted"](m[k],s)/model["shifted"](m[k-1],s)).log()/loop
            duration=arb.max(q(0),log_ratio(s0).union(log_ratio(s1)))
        steps=max(1,ceil(float(duration.upper())*8000))
        dt=duration/steps
        Bs=[[a+k*b for a,b in zip(row,ext)] for row,ext in zip(model["coeff"]["base"],model["coeff"]["extra"])]
        bs=[a+k*b for a,b in zip(flow["beta"],model["coeff"]["db"])]
        for _ in range(steps):
            assert 10*dt*ytot<1
            ym=ytot/(1-10*dt*ytot)
            gs=[v+arb(0,".0016") for v in gauge]
            assert min(gs)>0
            dg=[]
            for i in range(3):
                damp=q(0).union(max(model["coeff"]["yukawa"][i])*ym)
                dg.append(2*gs[i]**2*(bs[i]+(sum(x*y for x,y in zip(Bs[i],gs))-damp)/loop))
            assert all(abs(dt*v)<q(".00159") for v in dg)
            tube={j:v+arb(0,".002") for j,v in state.items()}
            d=rates(classes,tube,gs,ym)
            assert all(abs(dt*v)<q(".00199") for v in d.values())
            full={j:v+q(0).union(dt)*d[j] for j,v in state.items()}
            if mass_path is not None:
                mass_path.append(("step",dict(classes),full,gs,ym,dt))
            lower=min(v.lower() for (tag,_),v in full.items() if tag=="l")
            assert lower>0
            for (tag,ij),v in full.items():
                if tag=="a" and (bulk_portal is None or ij!=(0,2)):
                    assert v>0 and v+full["b",ij]>0
            if k>=3:
                lower_v=lower
                frob2=288*sum(classes[i]*abs(full["l",i]).upper()**2 for i in classes)
                for (tag,ij),v in full.items():
                    if tag=="a":
                        i,j=ij
                        num=classes[i]*classes[j] if i!=j else classes[i]*(classes[i]-1)//2
                        b=full["b",ij]
                        lower_v=min(lower_v,v.lower()/2,(v+b).lower()/2)
                        frob2+=96*num*(abs(v).upper()+2*abs(b).upper())**2
                if bulk_portal is not None:
                    lower_v=q(".00165")
                    bulk=full["l",2].lower()-lower_v
                    if classes[2]>1:
                        au,bu=full["a",(2,2)],full["b",(2,2)]
                        bulk=min(bulk,min(au.lower(),(au+bu).lower())/2-lower_v)
                    assert bulk>0 and lower>lower_v
                    au,bu=full["a",(0,2)],full["b",(0,2)]
                    mixed=min(au.lower(),(au+bu).lower())-2*lower_v
                    assert mixed+2*((full["l",0].lower()-lower_v)*bulk).sqrt()>0
                    for (tag,ij),v in full.items():
                        if tag=="a" and ij!=(0,2):
                            assert min(v.lower(),(v+full["b",ij]).lower())>2*lower_v
                Vmargin=arb.min(Vmargin,lower_v)
                Hnorm=arb.max(Hnorm,frob2.sqrt().upper())
            for i in classes:
                row=6*abs(full["l",i])
                for j,num in classes.items():
                    if num>int(i==j):
                        ij=pair(i,j)
                        a,b=abs(full["a",ij]),abs(full["b",ij])
                        row+=(num-int(i==j))*(2*a+b)
                        scalar_norm=arb.max(scalar_norm,(a+3*b).upper())
                scalar_norm=arb.max(scalar_norm,row.upper())
            margin=arb.min(margin,lower)
            largest=arb.max(largest,max(abs(v).upper() for v in full.values()))
            state={j:v+dt*d[j] for j,v in state.items()}
            gauge=[v+dt*f for v,f in zip(gauge,dg)]
            ytot=ym
    assert scalar_norm<q("5.35") and ytot<q(".295")
    assert 72*ytot<8*arb.pi()
    print("scalar contact norm",scalar_norm.str(8),"total Yukawa",ytot.str(8))
    assert margin>q(".0008") and largest<q("1.5")
    end=(Q*t1.exp()/model["shifted"](m[2],s0)).log()/loop
    if heavy_gain>0:
        return ytot,Hnorm,scalar_norm,Vmargin,raw,end
    print("all 54 doublets to the upper crossing face: minimum diagonal quartic above",
          margin.lower().str(8),"largest component below",largest.upper().str(8),
          "light quartic at crossing",float(state["l",0].lower()),float(state["l",0].upper()))


if __name__=="__main__":
    check()
