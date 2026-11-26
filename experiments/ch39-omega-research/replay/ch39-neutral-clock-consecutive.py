"""All-walk integral enclosure at consecutive zero ordinates, with the same frozen other inputs."""
from pathlib import Path
import runpy
from flint import acb,arb,arb_mat

m = runpy.run_path(str(Path(__file__).with_name("ch39-neutral-clock-susceptibility.py")))
pot,sh,raw,q,r = (m[k] for k in ("pot","sh","raw","q","r"))
qz = pot["m"]["source"]
axis = qz["gauge"]["spectrum"]
axes = (724,804)

def rotations(n):
    height = acb.zeta_zero(n).imag
    R = {}
    angles = []
    for k,y,z in ((724,q(0),q(1)),(804,axis["sin_axis"],axis["cos_axis"])):
        theta = height*(sh["theta"](k)-sh["theta"](244))
        angles.append(theta)
        s,c = (2*theta).sin_cos()
        # Ad(exp(i theta sigma)) rotates by -2 theta.
        J = arb_mat([[0,-z,y],[z,0,0],[-y,0,0]])
        I = arb_mat([[int(i==j) for i in range(3)] for j in range(3)])
        R[k] = I-s*J+(1-c)*(J*J)
    return R,angles

def point(z,upper,R):
    K = arb_mat(19,19)
    A = arb_mat(57,57)
    diag = [z for _ in range(19)]
    for i,j,d,ell,w in raw:
        if upper:
            y = (z.sqrt()*ell.lower()).lower()
            h = (y/y.tanh()+w.lower()).lower()
            b = (y/y.sinh()+w.upper()).upper()
        else:
            y = (z.sqrt()*ell.upper()).upper()
            h = (y/y.tanh()+w.upper()).upper()
            b = (y/y.sinh()+w.lower()).lower()
        diag[i]+=h
        diag[j]+=h
        K[i,j]=K[j,i]=-b
        if d in R:
            for a in range(3):
                for k in range(3):
                    A[3*i+a,3*j+k]=A[3*j+k,3*i+a]=-b*R[d][a,k]
        else:
            for a in range(3):
                A[3*i+a,3*j+a]=A[3*j+a,3*i+a]=-b
    for i,v in enumerate(diag):
        K[i,i]=v
        for a in range(3):
            A[3*i+a,3*i+a]=v
    # Positive witness for convergence of the path expansion.
    inv = K.solve(arb_mat([[q(1)] for _ in range(19)]))
    assert all(inv[i,0]>0 for i in range(19))
    return (A.det()/K.det()**3).log()

def cert(n):
    R,angle = rotations(n)
    weights = [sum(1+w.upper() for _,_,d,_,w in raw if d==k) for k in axes]
    s = [t.sin()**2 for t in angle]
    lo = q(".25")
    high = 8*sum(w*e for w,e in zip(weights,s))*(lo/19+q(18)*lo**2/(19*16))
    low = q(0)
    ratio = q(2)**(q(1)/32)
    for _ in range(32*32):
        hi = lo*ratio
        width = (hi*hi-lo*lo)/2
        low += width*point(hi,False,R).lower()
        high += width*point(lo,True,R).upper()
        lo=hi
    degree = pot["short"]["base"]["arc_budget"]()[1].upper()/2
    high += 608*degree**3/lo*sum(s)
    return low.lower().union(high.upper())/(32*arb.pi()**2*r**4)

def check():
    vals = [cert(n) for n in (718,719,720)]
    for n,val in zip((718,719,720),vals):
        print(n,val.lower().str(12),val.upper().str(12),flush=True)
    assert vals[2]<vals[0]<vals[1]
    return vals

if __name__=="__main__":
    check()
