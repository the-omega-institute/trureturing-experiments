"""Saturate the nine-overlap bound with explicit copied masses and two balanced heat moments."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root/"ch39-common-source-tower.py"))
data = source["gauge"]
old,g,q,r = (data[k] for k in ("old","inputs","q","r"))
D = old["D"]
dot,active,positive = (old[k] for k in ("dot","active","positive"))

def residual(T,s,logs,inputs,copied,steps=2048):
    changed = 18-copied
    X = inputs["G0"][:]
    y = inputs["y0"][:]
    dt = T/steps
    for j in range(steps):
        frac = q(Fraction(j,steps)).union(q(Fraction(j+1,steps)))
        t,rel = T*frac,T*(frac-1)
        thr = [rel+v for v in logs]
        ls = [positive(v) for v in thr]
        ns = [active(v) for v in thr]
        el,nl = sum(ls),sum(ns)
        er,nr = sum(ls[-copied:])+changed*positive(rel+s),sum(ns[-copied:])+changed*active(rel+s)
        ef,nf = positive(rel+1)/3,active(rel+1)/3
        ys = [v+arb(0,".006") for v in y]
        G = [x+arb(0,".001")-2*(b*t+l*el+rgt*er+f*ef)/D for x,b,l,rgt,f in
             zip(X,g["beta"],old["effL"],old["effR"],old["hf"])]
        assert min(G)>1/(4*arb.pi())
        gs = [1/v for v in G]
        xp,yp,fp = inputs["rates"](gs,ys,nl)
        shift = [v-ref for v,ref in zip(gs,old["ref"])]
        xp = [v-2*nr*dot(row,shift)/D**2+nf*f for v,row,f in zip(xp,old["BR"],fp)]
        assert all(abs(dt*v)<q(".001") for v in xp) and all(abs(dt*v)<q(".006") for v in yp)
        X = [v+dt*z for v,z in zip(X,xp)]
        y = [v+dt*z for v,z in zip(y,yp)]
    wl,wr0 = sum(logs),sum(logs[-copied:])
    target = [x*D/2-l*wl-f/3 for x,l,f in zip(X,old["effL"],old["hf"])]
    return dot(old["PT"],target)-T,(dot(old["PS"],target)-wr0)/changed-s

def heat(low,s,copied):
    changed = 18-copied
    d = data["data"]
    y19 = d["edge19"]/r**2
    upper = sorted(bound/(1-d["h"]**2)**2/r**2 for count,bound in d["factors"] for _ in range(count))
    def cost(x):
        return (x-1)*(-x).exp()
    # Min-max: the first Dirichlet trial on each edge has zero endpoint term.
    assert len(upper)==270 and min(upper)>y19>1
    reserve = sum(min(cost(y19).lower(),cost(v).lower()) for v in upper[18:])
    y = (-2*s).exp()
    low_difference = sum((-(x/r**2)).exp() for x in low[:changed])-changed*(-y).exp()
    mass = low_difference + sum((-v).exp() for v in upper[18:])
    gap = reserve-changed+changed*(1-y)*(-y).exp()
    assert mass>0 and gap>2/q(1).exp()
    return mass.lower(),gap.lower()

def verify(low,inputs,t0,t1,s0,s1,copied=9):
    logs = [r.log()-v.log()/2 for v in low]
    t0,t1,s0,s1 = map(q,(t0,t1,s0,s1))
    T,S = t0.union(t1),s0.union(s1)
    assert t0>max(logs) and min(logs)>s1 and s0>0
    faces = [residual(t,s,logs,inputs,copied)[i] for t,s,i in ((t0,S,0),(t1,S,0),(T,s0,1),(T,s1,1))]
    print("faces",[(float(v.lower()),float(v.upper())) for v in faces])
    print("remaining moments",*(v.str(12) for v in heat(low,S,copied)))
    assert faces[0]>0 and faces[1]<0 and faces[2]>0 and faces[3]<0

if __name__=="__main__":
    verify(data["low"],g,"26.85","27.03",".011",".30")
    verify(data["low"],g,"26.85","27.03",".42",".80",copied=8)
    inputs = dict(g)
    inputs["G0"] = g["G0"][:2]+[1/arb("1.0555",".00015")**2]
    verify(source["mass_bounds"](),inputs,"26.86","27.04",".001",".21")
