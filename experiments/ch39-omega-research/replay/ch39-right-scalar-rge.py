"""Two-loop gauge and one-loop Yukawa crossing with the slope-matched right scalar representations."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
data=runpy.run_path(str(root/"ch39-slope-matched-one-loop.py"))
g=data["g"]["g"]
q=g["q"]
D=g["loop"]
beta=g["beta"]
BL=g["BD"]
BF=g["BF"]
ref=g["gref"]
rf=g["gf"]
effL=g["eff"]
hf=g["hf"]
states=[(["4/5",0,"1/2"],["4/15",0,"4/3"]),
        (["1/5",0,"1/2"],["1/15",0,"4/3"]),
        (["3/5",0,0],["3/5",0,0])]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def project(null,target):
    w=[null[(i+2)%3]-null[(i+1)%3] for i in range(3)]
    return [v/dot(w,target) for v in w]
BR=[[sum(q(T[i])*(4*q(c[j])+q("2/3")*(0,2,3)[i]*int(i==j))
         for T,c in states) for j in range(3)] for i in range(3)]
effR=[q(w)+dot(row,ref)/D for w,row in zip(("8/15",0,"1/3"),BR)]
PT,PS=project(effR,beta),project(beta,effR)
r=q(400)
Ls=[r.log()-x.log()/2 for x in data["xs"]]
WL=sum(Ls)
# At r=400 exactly the first mode of the longest line contributes.
s0=q(5)
assert s0*(1+s0**2) < r/q(1).exp()
assert q(2).sqrt()*(1+2*s0**2)/(1+s0**2)>q(1).exp()
theta=data["theta"]
for a,b in ((724,216),(804,292)):
    k=theta(a)/theta(b)
    assert k*(1+k**2*s0**2)/(1+s0**2)>q(1).exp()


def positive(z):
    return g["positive"](z)


def active(z):
    if z>0:return q(1)
    if z<0:return q(0)
    return q(0).union(q(1))


def run(tau,s,steps=2048):
    X=g["G0"][:]
    y=g["y0"][:]
    dt=tau/steps
    for j in range(steps):
        rel=tau*q(Fraction(j-steps,steps)).union(q(Fraction(j+1-steps,steps)))
        t=tau*q(Fraction(j,steps)).union(q(Fraction(j+1,steps)))
        threshold=[rel+li for li in Ls]
        EL=sum(positive(v) for v in threshold)
        nL=sum(active(v) for v in threshold)
        ER=18*positive(rel+s)
        nR=18*active(rel+s)
        Z=positive(rel+1)/3
        nF=active(rel+1)/3
        tube=[v+arb(0,".001") for v in X]
        ys=[v+arb(0,".006") for v in y]
        G=[v-2*(bt*t+le*EL+re*ER+f*Z)/D for v,bt,le,re,f in zip(tube,beta,effL,effR,hf)]
        assert all(v>0 for v in G)
        gs=[1/v for v in G]
        xp,yp,fp=g["rates"](gs,ys,nL)
        shifted=[v-k for v,k in zip(gs,ref)]
        xp=[v-2*nR*dot(row,shifted)/D**2+nF*f for v,row,f in zip(xp,BR,fp)]
        assert all(abs(dt*v)<q(".001") for v in xp) and all(abs(dt*v)<q(".006") for v in yp)
        X=[v+dt*w for v,w in zip(X,xp)]
        y=[v+dt*w for v,w in zip(y,yp)]
    v=[D*x/2-e*WL-f/3 for x,e,f in zip(X,effL,hf)]
    R=dot(PT,v)-tau
    S=dot(PS,v)/18-s
    GU=X[0]-2*(beta[0]*tau+effL[0]*WL+effR[0]*18*s+hf[0]/3)/D
    return R,S,4*arb.pi()*GU


if __name__=="__main__":
    t0,t1=q("26.92"),q("26.97")
    s0,s1=q("2.55"),q("2.605")
    T,S=t0.union(t1),s0.union(s1)
    faces=[run(t0,S),run(t1,S),run(T,s0),run(T,s1),run(T,S)]
    for i,row in enumerate(faces):
        print("face",i,[(z.lower().str(10),z.upper().str(10)) for z in row])
    assert faces[0][0]>0 and faces[1][0]<0 and faces[2][1]>0 and faces[3][1]<0
    right_y=(-2*S).exp()
    y19=data["edge19"]/r**2
    count=sum(n for n,factor in data["factors"] if factor/(1-data["h"]**2)**2/r**2<2)
    low_weight=18*(-right_y).exp()
    Srest=sum((-(x/r**2)).exp() for x in data["xs"])+(count-18)*(-q(2)).exp()-low_weight
    excess=-18+(count-18)*(y19-1)*(-y19).exp()+(1-right_y)*low_weight
    assert Srest>0 and excess>1/q(1).exp()
    print("Q GeV, complementary weak-singlet m/Q", (1000*T.exp()).str(12),(-S).exp().str(12))
    print("high S, M-S lower",Srest.lower().str(12),excess.lower().str(12))
