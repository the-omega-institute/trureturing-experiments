"""One-loop gauge intersections with a zero flat-space volume and positive relative Einstein term."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

charge = runpy.run_path(str(Path(__file__).with_name("ch39-spin-gauge-charged.py")))
m,b = charge["m"],charge["b"]
q,wl,wr,beta = (m[k] for k in ("q","wl","wr","beta"))
pi,scale = arb.pi(),8*arb.pi()**2
s0,s1 = b["base"]["s0"],b["base"]["s1"]
def fit(lh,rh,hmax,gs):
    G = [1/g**2 for g in gs]
    l1,r1 = q(F(lh,90)),q(F(rh-6048,360))
    a,c = beta[1]-beta[2], beta[0]-beta[2]
    det = c-2*a*r1
    t = (scale*(G[0]-2*r1*G[1]+(2*r1-1)*G[2])-l1*wl)/det
    z = (2*scale*(c*G[1]-a*G[0]+(a-c)*G[2])+2*a*l1*wl)/det
    # All left thresholds precede the right cluster, which also precedes Q/e.
    assert min(b["base"]["logs"])>s1+b["shift"].upper()
    assert s0+b["shift"].lower()>1
    assert 18*(s0+b["shift"].upper())<z<18*(s1+b["shift"].lower())
    # The negative contribution of the eighteen right vectors and scalars starts
    # no later in aggregate than if all eighteen masses equalled their geometric mean.
    assert beta[0]+18*(l1+r1)+q(8)/3<0
    tpre = t-z/18
    floor = G[0]-(beta[0]*tpre+l1*(wl-z))/scale
    # SU(2) inverse coupling first increases, then decreases (the slope only rises).
    # SU(3) remains asymptotically free; the extra complete family turns on at Q/e.
    invU = G[2]-(beta[2]*t+q(8)/3)/scale
    assert beta[2]+q(8)/3<0 and min(G[1].lower(),invU.lower())>2/(4*pi)
    hyper = q(F(3*hmax*hmax,5*36))/floor
    return t,z,4*pi*invU,hyper

def check():
    s = b["model"]["gauge_inputs"]
    # The full one-sigma input box, followed by a center-neighborhood that lies
    # inside the original scalar-theory scale interval.
    core = [g.mid()+arb(0,(g.rad()/4).upper()) for g in s]
    for lh,rh,hmax in ((144,4716,30),(216,4068,27),(360,2916,24)):
        for label,gs in (("full",s),("central quarter",core)):
            t,z,alpha,y2 = fit(lh,rh,hmax,gs)
            print(lh,label,"T,WR,inv alpha, effective hypercharge squared cap",
                  [(float(v.lower()),float(v.upper())) for v in (t,z,alpha,y2)])
            if lh in (144,216):
                assert y2<4*pi
                if label=="central quarter":
                    assert t>b["base"]["T0"] and t<b["base"]["T1"]

if __name__=="__main__":
    check()
