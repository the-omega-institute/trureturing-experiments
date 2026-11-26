"""原壳 Weyl 点延续至零软化，并由四枚实际零点产生整数一。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def flat_certificate(shell, previous):
    soft=shell["soft"]
    spoke=acb_mat(19,19)
    for j in soft["star"]["orbit"]:
        i=soft["pole"]
        spoke[i,i]+=1
        spoke[j,j]+=1
        spoke[i,j]-=1
        spoke[j,i]-=1
    flat_soft=dict(soft,star=dict(soft["star"],orbit=()))
    def data(h,e,tau):
        v=h*h+soft["point"]["study"]["offset"]
        delta=arb(23)/25*arb(2)**-20/v
        ddelta=-2*h*delta/v
        def lap(q,height):
            return shell["laplacian"](q,height)+(1-tau)*delta*spoke
        proxy=dict(shell,soft=flat_soft,laplacian=lap)
        G,J,(H,D,B,dB,V)=previous["calculate"](proxy,h,e)
        correction=-tau*ddelta*spoke
        k=-(adj(V)*correction*V)[0,1]
        J[0,0]+=k.real
        J[1,0]+=k.imag
        return G,J,(H,D+correction,B,V)
    h,e=arb("80.9284256910337"),arb("8.92931408205483")
    eye=arb_mat([[1,0],[0,1]])
    for _ in range(2):
        g,j,_=data(h,e,arb(0))
        correction=j.solve(arb_mat([[g[0,1].real],[g[0,1].imag]]))
        h,e=(h-correction[0,0]).mid(),(e-correction[1,0]).mid()
    for tau,r in ((arb(0),arb("1e-23")),(arb(0).union(arb(1)),arb("2e-9"))):
        G0,J0,_=data(h,e,arb(0))
        pre=J0.inv()
        G,J,(H,D,B,V)=data(h+arb(0,r),e+arb(0,r),tau)
        defect=eye-pre*J
        rho=max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
        res=max(abs(z).upper() for z in (pre*arb_mat([[G0[0,1].real],[G0[0,1].imag]])).entries())
        if not tau.is_zero():
            Vt=data(h,e,tau)[2][3]
            delta=arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
            step=delta*(adj(Vt)*spoke*Vt)[0,1]
            res+=max(abs(z).upper() for z in (pre*arb_mat([[step.real],[step.imag]])).entries())
        assert rho<1 and res+rho*r<r, (res,rho)
        t=[-1/G[i,i].real for i in range(2)]
        assert -1<t[0]<0<t[1]
        center=data(h,e,arb(0))[2][0]
        vals=center.eig(multiple=True)
        change=sum(abs(z).upper()**2 for z in (H-center).entries()).sqrt()
        assert len(vals)==19 and all(abs(z.real-e)>change+r for z in vals) and sum(z.real<e for z in vals)==3
        v=acb_mat([[V[i,0]] for i in range(19)])
        v=v/(adj(v)*v)[0,0].real.sqrt()
        w=acb_mat([[V[i,1]] for i in range(19)])
        w=w-v*(adj(v)*w)[0,0]
        w=w/(adj(w)*w)[0,0].real.sqrt()
        basis=acb_mat([[v[i,0],w[i,0]] for i in range(19)])
        ders=[D]
        for i in range(2):
            col=acb_mat([[B[j,i]] for j in range(19)])
            ders.append(col*adj(col))
        pauli=[]
        for m in ders:
            A=adj(basis)*m*basis
            pauli.append([A[0,1].real,-A[0,1].imag,(A[0,0].real-A[1,1].real)/2])
        det=arb_mat(pauli).det()
        assert det>0
        print("软化乘子、求根半径、压缩率、速度行列式严格下界",
              "0" if tau.is_zero() else "[0,1]",r.str(7),rho.str(8),det.lower().str(10))
        if tau.is_zero():
            flat_root = (h+arb(0,r), e+arb(0,r), [1+a for a in t])
            print("零软化的高度、交叉本征值与两条正电导",
                  (h+arb(0,r)).str(18),(e+arb(0,r)).str(18),*[(1+a).str(16) for a in t])
    return flat_root


def tetrahedron(shell):
    table=((9,".0005","27"),(4,".0005","27"),(16,".016","7"),(10,".0005","27"))
    frames=[]
    for n,x,y in table:
        H=shell["laplacian"](-3,acb.zeta_zero(n).imag)
        for (i,j),w in zip(((0,2),(5,14)),map(arb,(x,y))):
            H[i,i]+=w-1
            H[j,j]+=w-1
            H[i,j]*=w
            H[j,i]=H[i,j].conjugate()
        vals,v=H.eig(right=True,algorithm="rump")
        ix=sorted(range(19),key=lambda i:float(vals[i].real))
        assert all(vals[i].real<vals[j].real for i,j in zip(ix[:3],ix[1:4]))
        cols=[]
        for k in ix[:3]:
            col=[v[i,k] for i in range(19)]
            norm=sum(abs(z)**2 for z in col).sqrt()
            cols.append([z/norm for z in col])
        frames.append(acb_mat(cols).transpose())
    links={}
    for i in range(4):
        for j in range(i+1,4):
            z=(adj(frames[j])*frames[i]).det()
            assert abs(z) > 0
            links[i,j]=z/abs(z)
            links[j,i]=links[i,j].conjugate()
    angles=[]
    for i,j,k in ((1,2,3),(0,3,2),(0,1,3),(0,2,1)):
        angle=(links[i,j]*links[j,k]*links[k,i]).arg()
        assert 1<angle<2
        angles.append(angle)
    value=sum(angles)/(2*arb.pi())
    assert value.unique_fmpz()==1
    print("四枚实际零点及两权",table)
    print("四张有向面的主值与总整数",*[z.str(15) for z in angles],value.unique_fmpz())


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    previous=run_path(str(root/"ch39-first-chern-edge.py"))
    ctx.prec=256
    flat_certificate(shell,previous)
    tetrahedron(shell)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
