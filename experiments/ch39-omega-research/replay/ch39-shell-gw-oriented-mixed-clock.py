"""Recover the missing mixing coordinate and test a coupling to shortest shell loops."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

folder=Path(__file__).resolve().parent
fiber=run_path(str(folder/"ch39-shell-gw-pencil-fiber.py"))
m=fiber["m"]
adj=m["adj"]


def mass_clock():
    omega=acb(0,2*arb.pi()/3).exp()
    S=acb_mat([[omega**i if i==j else 0 for j in range(3)] for i in range(3)])
    V=m["ckm"](arb("1.139"))
    T=V*S*adj(V)
    W=S*T*adj(S)*adj(T)
    trace=sum(W[i,i] for i in range(3))
    J=(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag
    B=fiber["fiber"]()
    even=sum(z*z for row in B for z in row)*(arb(9)/4)-arb(15)/4
    assert (trace.real-even(arb(0))).contains(0) and (trace.imag+27*J).contains(0)
    slope=even.derivative()
    # The derivative is affine everywhere.
    assert min(slope(arb("-.020")),slope(arb(".045")))>7
    print("mass commutator trace",trace.str(12),"orientation-even slope above seven")


def shell_sign(index,width):
    w=run_path(str(folder/"ch39-shell-gw-wilson-coordinates.py"))
    height=acb.zeta_zero(index).imag
    p,q=[height*(a-w["angles"][0]) for a in w["angles"][1:]]
    pts=w["pts"]
    short={}
    for i in (0,2,4):
        j=next(j for j in range(6) if sum(a*b for a,b in zip(pts[i],pts[j]))==25)
        short[pts[i].index(1)]=acb_mat(w["model"]["link"](pts[i],pts[j],p,q))
    for k in range(3):
        v=adj(short[k])*short[(k+1)%3]
        odd=sum(v[j,j] for j in range(3)).imag
        assert (odd+2*(1-p.cos()*q.cos())*p.sin()*q.cos()).contains(0)
    p,q=p+arb(0,width),q+arb(0,width)
    odd=-2*(1-p.cos()*q.cos())*p.sin()*q.cos()
    if index==43:
        assert odd<arb("-.01396")
    else:
        assert odd>arb(".02160")
    print("shell orientation",index,odd.str(9))


def selected():
    B=fiber["fiber"]()
    a,b,c=[B[0][i]*B[1][i] for i in range(3)]
    square=4*a*b-(c-a-b)**2
    deriv=square.derivative()
    interval=arb("-.020").union(arb(".045"))
    assert deriv.derivative()(interval)<0
    t=fiber["isolate"](deriv,".001",".03")
    J=square(t).sqrt()/2
    a,b,c=fiber["angles"](B,t)
    cosine=(B[1][0](t)-(a.sin()*c.cos())**2-(a.cos()*c.sin()*b.sin())**2)
    cosine/=2*a.sin()*a.cos()*c.sin()*c.cos()*b.sin()
    assert arb("1.42148")<cosine.acos()<arb("1.42149")
    print("phase-magnitude maximizer",t.str(10),J.str(12),"phase magnitude",cosine.acos().str(11))


if __name__=="__main__":
    mass_clock()
    shell_sign(43,arb("3e-8"))
    shell_sign(195,arb("1e-8"))
    selected()
