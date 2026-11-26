"""The nine-to-ten shell-count inertia ratio, with the common scale determined by the four equations."""
from pathlib import Path
from types import FunctionType
import runpy

from flint import arb

root = Path(__file__).resolve().parent
base = runpy.run_path(str(root/"ch39-right-arc-inertia.py"))
src = base["source"]
profile = base["model"]
q = base["q"]
r0,r1 = q(356),q(430)
v0,v1 = q("2.205"),q("2.326")
T0,T1 = q("26.88"),q("27.025")
eta = q("9/10")
u0,u1 = q("-.030"),q(".030")
trace,norm = base["arc_budget"]()
errS,errM = base["errS"],base["errM"]
comp = base["comparison"]
lower = comp["all_bounds"](src["static"]([w.lower() for w in base["weights"]]))
upper = comp["all_bounds"](src["static"]([w.upper() for w in base["weights"]]))
left_static = [a.lower().union(b.upper()) for a,b in zip(lower,upper)]
inputs = dict(profile["inputs"])
# The flavor fit uses the measured top normalization. Keep all nine Yukawa boxes;
# the three gauge measurements are taken in a small box around the Table 2 centers.
inputs["G0"] = [1/arb(v,rad)**2 for v,rad in
                ((".467453",".000003"),(".63811",".000003"),("1.0583",".00025"))]

def chart(r,v,u):
    y = (-2*v).exp()*(400/r)**2
    ys = [x/r**2 for x in left_static]
    SL = sum((-x).exp() for x in ys)
    ML = sum(x*(-x).exp() for x in ys)
    rows = [(n,b*(400/r)**2) for n,b in base["drows"]]
    ymin = base["ymin"]*(400/r1)**2
    ratio = (-q(".75")*8*ymin).exp()
    assert 64*ymin>20 and (q(20)/4).exp()>20
    tail = 2*270*ratio**8/(1-ratio)
    fn = base["dchange"]
    change = FunctionType(fn.__code__,dict(fn.__globals__,drows=rows,tail=tail))
    diffS,diffM = change(eta,False),change(eta,True)
    w = SL-diffS-18*(-y).exp()+u
    assert w>0
    z = (20/w).log()
    moment = 18*y*(-y).exp()+w*z+diffM-ML+arb(0,(trace*(400/r)**2+2*errM).upper())
    return z,moment

def heat():
    v,u = v0.union(v1),u0.union(u1)
    zs=[]
    for j in range(128):
        r=r0+(r1-r0)*q(j).union(q(j+1))/128
        z,_=chart(r,v,u)
        zs.append(z)
    zmin=min(z.lower() for z in zs)
    zmax=max(z.upper() for z in zs)
    print("high right potential",zmin.str(10),zmax.str(10))
    assert zmin>q("1.004")
    margin=u1-trace*(400/r0)**2-2*errS
    assert margin>0
    vals=[[],[]]
    for i,rr in enumerate((r0,r1)):
        for j in range(12):
            vpart=v0+(v1-v0)*q(j).union(q(j+1))/12
            for k in range(8):
                upart=u0+(u1-u0)*q(k).union(q(k+1))/8
                vals[i].append(chart(rr,vpart,upart)[1])
    m0=min(w.lower() for w in vals[0])
    m1=max(w.upper() for w in vals[1])
    print("moment faces",m0.lower().str(10),m1.upper().str(10),"heat margin",margin.lower().str(10))
    assert m0>0 and m1<0
    return zmin

def gauge(zmin):
    v=v0.union(v1)
    T=T0.union(T1)
    r=r0.union(r1)
    # The low onsite mass in units of M is independent of r in these coordinates.
    mu=400**2*(-2*v).exp()
    ceiling=(mu+norm).upper()
    d,beta=base["data"]["d"],base["data"]["beta"]*eta
    k=arb.pi()*(ceiling*eta/d).sqrt()
    slope=18*k*(k/2).tan()/ceiling
    dlo=(-(1+norm/mu.lower()).log()/2).lower()
    dhi=(((1+slope)*(1+beta*(ceiling+10))).log()/2).upper()
    assert r1**2<base["data"]["edge19"]
    assert r0>5*26*q(1).exp()
    # Check the high onsite sector after removing the 18 low components.
    factor=(1+10*beta)/(1-r1**2*beta)
    k=arb.pi()*(r1**2*factor*eta/d).sqrt()
    assert k<arb.pi() and (zmin-factor)*r0**2>36*k/k.sin()
    # The same log(r/400) occurs in the two matter exposures. Remove it analytically
    # from the endpoint projections while keeping the true threshold positions in the flow.
    fn=profile["integrate"]
    old=profile["old"]
    coeffT=-18*old["dot"](old["PT"],old["effL"])
    coeffV=-old["dot"](old["PS"],old["effL"])-1
    # Evaluate the integral with true threshold positions, but carry only the reference
    # left exposure in the projected boundary residual. The shift is restored below.
    def run(t,w,logr):
        rz=w+logr+dlo.union(dhi)
        ls=[a+logr for a in base["logs"]]
        calc=FunctionType(fn.__code__,dict(fn.__globals__,Ls=ls,WL=base["wl"],g=inputs),
                          argdefs=fn.__defaults__)
        A,B,_=calc(t,rz)
        return A+coeffT*logr,B+logr+coeffV*logr
    # Avoid subtracting the uncertain shift twice in the second residual.
    # Split the range before computing endpoints; unions cover every physical scale.
    results=[[] for _ in range(4)]
    for j in range(10):
        r=r0+(r1-r0)*q(j).union(q(j+1))/10
        logr=(r/400).log()
        for i,(t,w) in enumerate(((T0,v),(T1,v),(T,v0),(T,v1))):
            results[i].append(run(t,w,logr)[i//2])
    # Arb's comparison is partial. Compare individual endpoints for the extrema.
    bounds=[min(z.lower() for z in results[0]),max(z.upper() for z in results[1]),
            min(z.lower() for z in results[2]),max(z.upper() for z in results[3])]
    print("gauge face margins",*(z.str(10) for z in bounds))
    assert bounds[0]>0 and bounds[1]<0 and bounds[2]>0 and bounds[3]<0

if __name__=="__main__":
    gauge(heat())
