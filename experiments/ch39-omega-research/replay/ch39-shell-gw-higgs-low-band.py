"""Verified commutator bound for singlet and adjoint Higgs modes below a common spectral cut."""
from itertools import product
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb,acb_mat,arb,ctx


folder=Path(__file__).resolve().parent
yuk=run_path(str(folder/"ch39-shell-gw-local-yukawa.py"))
model=yuk["model"]["model"]
pts=model["points"]
ctx.prec=256
adj=yuk["adjoint"]
trace=lambda a: sum(a[i,i] for i in range(a.nrows()))


def hs(a,b):
    return trace(adj(a)*b)


def basis():
    third=arb(3).sqrt()
    out=[acb_mat([[1,0,0],[0,-1,0],[0,0,0]])/arb(2).sqrt(),
         acb_mat([[1,0,0],[0,1,0],[0,0,-2]])/(arb(2).sqrt()*third)]
    for i in range(3):
        for j in range(i):
            for phase in (acb(1),acb(0,1)):
                v=acb_mat(3,3)
                v[i,j]=phase/arb(2).sqrt()
                v[j,i]=v[i,j].conjugate()
                out.append(v)
    return out


def transported_spectra():
    h=acb.zeta_zero(43).imag
    angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
    p,q=h*(angles[1]-angles[0]),h*(angles[2]-angles[0])
    spin=[];phases=[];G=basis()
    for i in range(6):
        x,y,z=pts[i]
        a,b=arb(26).sqrt()+z,acb(x,y)
        chi=(arb(3).sqrt()+1)*a+acb(1,-1)*b
        t=chi/abs(chi)
        spin.append((a/t,b/t));phases.append(t.conjugate()**2)
    sing=yuk["identity"](6)*4
    family=yuk["identity"](48)*4
    for i,j in product(range(6),repeat=2):
        if i==j or j==(i^1):continue
        v=sum(spin[j][r].conjugate()*spin[i][r] for r in range(2))
        edge=-(v.conjugate()/abs(v))**2
        sing[i,j]=edge
        U=acb_mat(model["link"](pts[j],pts[i],p,q))
        for a,src in enumerate(G):
            dest=U*src*adj(U)
            for b,tgt in enumerate(G):
                family[8*i+b,8*j+a]=edge*hs(tgt,dest)
    return sing,family,phases,G


def low_frame(k,d):
    n=k.nrows();rows=[[k[i,j] for j in range(n)] for i in range(n)]
    count=model["model"]["count"]
    assert count(rows,arb(".02"))==count(rows,arb(1))==d
    vals,vectors=np.linalg.eigh(np.array([[complex(z) for z in row] for row in rows]))
    frame=[]
    for i in range(d):
        v=acb_mat([[acb(float(z.real),float(z.imag))] for z in vectors[:,i]])
        trial=arb(float(vals[i]))
        assert trial<arb(".02")
        residual=k*v-v*trial
        error=hs(residual,residual).real.sqrt()/(1-trial)
        assert error<arb("1e-11")
        e=arb(0,error.upper())
        frame.append(acb_mat([[v[j,0]+acb(e,e)] for j in range(n)]))
    gram=acb_mat([[hs(a,b) for b in frame] for a in frame])
    assert gram.det().real>0
    return frame


def low_yukawas():
    sing,high,phases,g=transported_spectra()
    left,right=yuk["zero_basis"](1),yuk["zero_basis"](-1)
    blocks=[]
    for v,basis0 in [(v,[yuk["identity"](3)]) for v in low_frame(sing,1)]+[(v,g) for v in low_frame(high,2)]:
        out=acb_mat(3,3)
        for i in range(6):
            vertex=acb_mat(3,3)
            for j,h in enumerate(basis0):
                vertex+=h*v[len(basis0)*i+j,0]
            a=acb_mat([[left[3*i+j,b] for b in range(3)] for j in range(3)])
            b=acb_mat([[right[3*i+j,t] for t in range(3)] for j in range(3)])
            out+=phases[i]*a.transpose()*vertex.transpose()*b
        blocks.append(out)
    return blocks


if __name__=="__main__":
    mats=[]
    for x in low_yukawas():
        for y in mats:x-=hs(y,x)*y
        mats.append(x/hs(x,x).real.sqrt())
    pairs=[a*adj(b) for a in mats for b in mats]
    commutators=[a*b-b*a for a in pairs for b in pairs]
    bound=sum(hs(c,c).real for c in commutators).sqrt()
    assert bound<arb(".0031")
    # Common-scale reference matrices normalized by their own Hilbert-Schmidt norm.
    def normed_squared(weights):
        t=[arb(w)**2 for w in weights];return acb_mat([[z/sum(t) if i==j else 0 for j in range(3)] for i,z in enumerate(t)])
    u=normed_squared((".00000615",".00311",".8616"))
    d=normed_squared((".0000135",".000268",".01401"))
    v=yuk["identity"](3)
    for i,j,angle,phase in ((1,2,".04275","0"),(0,2,".00377","-1.139"),(0,1,".22704","0")):
        a=arb(angle); mat=yuk["identity"](3);mat[i,i]=mat[j,j]=a.cos()
        mat[i,j]=acb(0,arb(phase)).exp()*a.sin();mat[j,i]=-mat[i,j].conjugate()
        v=v*mat
    d=v*d*adj(v);c=u*d-d*u
    obs=hs(c,c).real.sqrt()
    assert obs>arb(".0602")
    print("low-band commutator <",bound.upper().str(10),"input",obs.str(12),
          "necessary fractional component outside the low span >",((obs-bound)/4).lower().str(9))
