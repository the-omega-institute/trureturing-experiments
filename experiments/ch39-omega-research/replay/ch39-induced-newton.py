"""Proper-time induced Einstein coefficient of the two matched heavy scalar towers."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/"ch39-fixed-inertia-scale.py"))
base,q = m["base"],m["q"]
r0,r1 = q(390),q(429)
eta = m["eta"]
pi = arb.pi()

def en(x,n):
    # E_n is decreasing on the positive real axis.
    return x.upper().expint(n).lower().union(x.lower().expint(n).upper())

def traces(r,v,u,n):
    z,_ = m["chart"](r,v,u)
    scale = (400/r)**2
    y = (-2*v).exp()*scale
    left = sum(en(x/r**2,n) for x in m["left_static"])
    right = 18*en(y,n)+20*en(z,n)
    # Monotonicity also gives the error sign of the positive right graph Laplacian.
    right -= arb(0).union((en(y,n-1)*m["trace"]*scale).upper())
    for count,bb in base["drows"]:
        a = bb*scale
        left += count*en(a,n)
        right += count*en(a/eta,n)
    ymin = base["ymin"]*(400/r1)**2
    ratio = (-8*ymin).exp()
    tail = 270*ratio**8/(1-ratio)/(n-1)
    err = arb(0,(m["errS"]/(n-1)).upper())+arb(0).union(tail.upper())
    return left+err,right+err

def split_difference(r,v,u,n):
    # Integrate the Dirichlet difference at the same squared mass to retain its correlation.
    z,_ = m["chart"](r,v,u)
    scale = (400/r)**2
    y = (-2*v).exp()*scale
    value = sum(en(x/r**2,n) for x in m["left_static"])-18*en(y,n)-20*en(z,n)
    for count,bb in base["drows"]:
        a = bb*scale
        # h(a)=E_n(a)-E_n(a/eta). Bound the variation around the midpoint.
        mid = a.mid()
        slope = -en(a,n-1)+en(a/eta,n-1)/eta
        value += count*(en(mid,n)-en(mid/eta,n)+(a-mid)*slope)
    ymin=base["ymin"]*(400/r1)**2
    ratio=(-8*ymin).exp()
    value += arb(0,(2*m["errS"]/(n-1)+540*ratio**8/((1-ratio)*(n-1))).upper())
    return value+arb(0).union((en(y,n-1)*m["trace"]*scale).upper())

def verify():
    v = m["v0"].union(m["v1"])
    u = m["u0"].union(m["u1"])
    lefts,rights,diffs = [],[],[]
    for i in range(32):
        r = r0+(r1-r0)*q(i).union(q(i+1))/32
        a,b = traces(r,v,u,2)
        lefts.append(a); rights.append(b)
        diffs.append(split_difference(r,v,u,2))
    total = min(a.lower()+b.lower() for a,b in zip(lefts,rights)).union(
            max(a.upper()+b.upper() for a,b in zip(lefts,rights)))
    ratio = 16*total/(12*pi)
    Q=1000*m["T0"].union(m["T1"]).exp()
    planck=q("1.220890e19")
    fraction=ratio*(Q/planck)**2
    d=min(x.lower() for x in diffs).union(max(x.upper() for x in diffs))
    print("one-component proper-time difference",d.lower().str(12),d.upper().str(12))
    print("induced inverse Newton / Q^2",ratio.lower().str(12),ratio.upper().str(12))
    print("fraction of inverse measured Newton",fraction.lower().str(12),fraction.upper().str(12))
    assert ratio>0 and fraction<q("1e-6")
    return d,ratio

if __name__=="__main__":
    verify()
