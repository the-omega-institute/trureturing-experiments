"""Explicit strength sufficient to expose the calibrated height in the presence of all fitted spin spectra."""
from pathlib import Path
import runpy
from flint import arb,arb_mat,acb_mat

root=Path(__file__).resolve().parent
pent=runpy.run_path(str(root/"ch39-five-clock-selection.py"))
m,q= pent["m"],pent["q"]
short=m["short"]
pi,r=arb.pi(),short["r"]
axes=(724,804)
threshold=q(".25")
degree=short["base"]["arc_budget"]()[1].upper()/2

def mark_sums(z,raw,V):
    # The positive upper envelope has g from lower input bounds, b from upper ones.
    diag=[z+v.lower() for v in V]
    K=arb_mat(19,19)
    C=arb_mat(19,19)
    pairs=[]
    for i,j,d,_,ell,w in raw:
        if z.is_zero():
            a=b=q(1)
        else:
            y=(z.sqrt()*ell.lower()).lower()
            a,b=y/y.tanh(),y/y.sinh()
        v=(a+w.lower()).lower()
        b=(b+w.upper()).upper()
        diag[i]+=v
        diag[j]+=v
        K[i,j]=K[j,i]=-b
        if d in axes:
            C[i,j]=C[j,i]=b
            pairs.append((i,j,b))
    for i,v in enumerate(diag):
        K[i,i]=v
    G=K.inv()
    assert all(sum(G[i,j] for j in range(19))>0 for i in range(19))
    GC=G*C
    bounce=2*sum(b*b/(diag[i]*diag[j]) for i,j,b in pairs)
    first=GC.trace()-bounce
    second=GC.trace()+(GC*GC).trace()-2*bounce
    return first.upper(),second.upper()

def integrate(raw,V,start,rank,speed):
    lo=start
    hi=max(q(".25"),start*q(2)**(q(1)/8))
    result=[q(0),q(0)]
    ratio=q(2)**(q(1)/8)
    while lo<q(2)**30:
        result=[u+(hi*hi-lo*lo)*v/2 for u,v in zip(result,mark_sums(lo,raw,V))]
        lo,hi=hi,hi*ratio
    # n>=3 with rho <= 2*degree/z < 1/2. Use k<=n for every marked walk.
    result[0]+=38*(2*degree)**3/lo
    result[1]+=152*(2*degree)**3/lo
    norm=32*pi**2*r**4
    a,b=[u/norm for u in result]
    return 2*rank*a,rank*speed*a,rank*speed**2*b

def neutral_gap(raw):
    # A constant basis realizes both adjoint rotations as real 3x3 matrices.
    axis=m["m"]
    Js=[]
    for y,z in ((q(0),q(1)),(axis["v"],axis["c"])):
        Js.append(arb_mat([[0,-z,y],[z,0,0],[-y,0,0]]))
    I=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    target=pent["target"]
    maxdeg=max(sum(1+w.upper() for i,j,d,_,_,w in raw if d in axes and k in (i,j)) for k in range(19))
    half=pent["radius"]/8
    # Use a certified center eigenvalue and the angular Lipschitz bound on the block row sums.
    floor=q(".04")
    for x in range(8):
        for y in range(8):
            ps=[h+half*(2*j-7) for h,j in zip(target,(x,y))]
            R={}
            for d,u,J in zip(axes,ps,Js):
                ss,cc=(2*u).sin_cos()
                R[d]=I-ss*J+(1-cc)*(J*J)
            A=acb_mat(57,57)
            for i,j,d,_,_,w in raw:
                b=1+w.lower()
                for k in range(3):
                    A[3*i+k,3*i+k]+=b
                    A[3*j+k,3*j+k]+=b
                    for l in range(3):
                        A[3*i+k,3*j+l]=A[3*j+l,3*i+k]=-b*R.get(d,I)[k,l]
            cut=floor+2*half*maxdeg
            assert axis["source"]["inertia"](A,cut)==0
    return floor,maxdeg

def check():
    left,Vl=m["rows"]("left")
    right,Vr=m["rows"]("right")
    neutral=right
    gap,deg=neutral_gap(neutral)
    N=list(integrate(neutral,[q(0)]*19,threshold,3,2))
    norm=32*pi**2*r**4
    W=sum(1+w.upper() for i,j,d,_,_,w in neutral if d in axes)
    # For the global oscillation on the omitted interval, 0 <= F <= 8 W*(1/19z+18/152).
    N[0]+=8*W*(threshold/19+18*threshold**2/304)/norm
    # Near the fitted angles every eigenvalue of K(z) exceeds gap+z.
    # The sharper trace-norm constants are unnecessary here.
    K1,K2=2*deg,4*deg
    N[1]+=57*K1*threshold**2/(2*gap*norm)
    N[2]+=57*(K2/gap+(K1/gap)**2)*threshold**2/(2*norm)
    L=integrate(left,Vl,q(0),2,1)
    R=integrate(right,Vr,q(0),2,1)
    osc,grad,hess=[16*(l+r)+3*n for l,r,n in zip(L,R,N)]
    inv=pent["C"].inv()
    shift=[sum(inv[i,j]**2 for j in range(2)).sqrt()*grad for i in range(2)]
    shift=[v.upper() for v in shift]
    positivity=max(shift[i]/pent["weights"][i+1].lower() for i in range(2))
    outside=(osc+9*sum(shift))/q(".015")
    inside=(hess+80*sum(shift))/5
    needed=max(x.upper() for x in (positivity,outside,inside))
    print("dimensionless spectrum oscillation, gradient and Hessian caps",*[v.upper().str(12) for v in (osc,grad,hess)])
    print("individual coefficient corrections",*[v.str(12) for v in shift])
    print("sufficient common weight above",needed.str(12))
    assert needed<q("9e-8")
if __name__=="__main__":
    check()
