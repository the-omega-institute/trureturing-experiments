"""Test the convex combination of complete one-loop spectra on both short-loop comparison squares."""
from pathlib import Path
import runpy
from flint import arb,arb_mat

root=Path(__file__).resolve().parent
cone=runpy.run_path(str(root/"ch39-clock-short-wilson-cone.py"))
neutral=runpy.run_path(str(root/"ch39-neutral-clock-consecutive.py"))
m,q,pi=cone["m"],cone["q"],cone["pi"]
pad=cone["half"]
r=m["short"]["r"]

def neutral_upper(p,t,steps):
    # Evaluate the actual twist on the second square by positive walk comparison.
    axes=(724,804)
    c=m["m"]["c"]
    v=m["m"]["v"]
    I=arb_mat([[int(i==j) for i in range(3)] for j in range(3)])
    raw=neutral["raw"]
    weights=[sum(1+w.upper() for _,_,d,_,w in raw if d==a) for a in axes]
    sigma=[abs((h+arb(0,pad.upper())).sin()).upper()**2 for h in (p,t)]
    lo=q(".25")
    high=8*sum(w*h for w,h in zip(weights,sigma))*(lo/19+q(18)*lo**2/(19*16))
    ratio=q(2)**(q(1)/steps)
    # Use nominal rotations in the determinant and the absolute marked-walk
    # bound for each finite angular displacement (the walk series respects unitarity).
    exact={}
    for d,y,z,theta in ((724,q(0),q(1),p),(804,v,c,t)):
        s,b=(2*theta).sin_cos()
        J=arb_mat([[0,-z,y],[z,0,0],[-y,0,0]])
        exact[d]=I-s*J+(1-b)*(J*J)
    raw2=[(i,j,d,0,ell,w) for i,j,d,ell,w in raw]
    for _ in range(steps*32):
        hi=lo*ratio
        val=neutral["point"](lo,True,exact)
        diag,bs,_,_=m["bounds"](lo,True,raw2,[q(0)]*19)
        marks=total_marks(diag,bs,raw2)
        high+=(hi*hi-lo*lo)/2*(val+4*pad*marks)
        lo=hi
    degree=m["short"]["base"]["arc_budget"]()[1].upper()/2
    return (high+608*degree**3/lo*sum(sigma))/(32*pi**2*r**4)

def total_marks(diag,bs,raw):
    K=arb_mat(19,19)
    for i,v in enumerate(diag):
        K[i,i]=v
    for (i,j),w in bs.items():
        K[i,j]=K[j,i]=-w
    G=K.inv()
    return 2*sum(bs[i,j]*(G[i,j]-bs[i,j]/(diag[i]*diag[j])) for i,j,d,_,_,_ in raw if d in m["axes"])

def charged(side,steps):
    raw,V=m["rows"](side)
    U0,a0=m["rotations"](719)
    rings=[]
    for pair in cone["angles"]:
        U=cone["rotations"](*pair)
        ds=[4*abs(((x-y)/2).sin()).upper() for x,y in zip(pair,a0)]
        rings.append((U,ds,m["short_cycles"](raw,U0,U)))
    lo=q(0)
    hi=q(".125")
    ratio=q(2)**(q(1)/steps)
    out=q(0)
    for _ in range(steps*33+1):
        gd,bs,low,sl=m["bounds"](hi,False,raw,V)
        D,B,high,su=m["bounds"](lo,True,raw,V)
        span=2*pad*total_marks(D,B,raw)
        vals=[]
        for U,ds,cycles in rings:
            val=m["log_ratio"](gd,bs,raw,U0,U)
            extra=sum(d*(h-l) for d,h,l in zip(ds,high,low)).upper()
            pos=neg=q(0)
            for cycle,c in cycles:
                a=b=q(1)
                for i,j in cycle:
                    a*=sl[i,j]
                    b*=su[i,j]
                diff=(b-a)*c
                pos+=max(q(0),diff.upper())
                neg+=min(q(0),diff.lower())
            vals.append(val+neg.lower().union(pos.upper())+arb(0,extra))
        out+=(hi*hi-lo*lo)/2*((vals[0]+2*vals[1])/3+arb(0,span.upper()))
        lo,hi=hi,hi*ratio
    degree=m["short"]["base"]["arc_budget"]()[1].upper()/2
    tail=152*(2*degree)**3/lo
    return (out+arb(0,tail.upper()))/(32*pi**2*r**4)

def check():
    steps=16
    N0=q("9.9583e-8")*2*pad**2
    N1=neutral_upper(*cone["angles"][1],steps)
    left,right=charged("left",steps),charged("right",steps)
    neutral_change=N0+2*N1-3*q("8.37729e-12")
    ferm_bose=16*(right-left)
    upper=(neutral_change+ferm_bose).upper()
    print("neutral mixture subtotal upper",neutral_change.upper().str(12),"spin-pair upper",ferm_bose.upper().str(12))
    print("convex one-loop energy difference upper",upper.str(12))
    assert upper<0
if __name__=="__main__":
    check()
