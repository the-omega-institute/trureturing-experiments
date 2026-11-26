"""十五个循环块的前二谱投影在全有效双钟环面与物理软亏额上平凡。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, acb_mat, arb, ctx


def adj(A):
    return A.conjugate().transpose()


def sub(A,rows,cols):
    return acb_mat([[A[i,j] for j in cols] for i in rows])


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    light=run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec=160
    bases=list(light["rotation_bases"](shell))
    # 内站是 (13,13,23) 的三排列。软辐边两端都在外站。
    sites=shell["geo"]["sites"]
    inner={j for j,p in enumerate(sites) if sorted(p)==[13,13,23]}
    softened={shell["soft"]["pole"],*shell["soft"]["star"]["orbit"]}
    assert len(inner)==3 and inner.isdisjoint(softened)
    eps=6*arb("1e-8")
    for q in (1,-4,2,-3,6):
        base,links=light["clock_matrices"](shell,q)
        corners=[[],[]]
        for k,Q in enumerate(bases):
            A,B,C=[adj(Q)*X*Q for X in (base,*links)]
            for point,v in enumerate((-1,1)):
                H=A-B-adj(B)+v*(C+adj(C))
                sample=H.eig(multiple=True)
                assert all(e.imag.contains(0) for e in sample)
                values=sorted((e.real for e in sample),key=float)
                assert values[0]<values[1]<values[2]
                corners[point].append(values[1]+arb(0,eps.upper()))
            n=A.nrows()
            z,=[j for j in range(n) if any(abs(Q[i,j])>0 for i in inner)]
            outside=[i for i in range(n) if i!=z]
            assert all(Q[i,z]==0 for i in range(19) if i not in inner)
            assert all(B[i,j]==C[i,j]==0 for i in range(n) for j in range(n) if j!=z or i==z)
            D=sub(A,outside,outside)
            c=A[z,z].real
            eig,vec=D.eig(right=True,algorithm="rump")
            order=sorted(range(len(eig)),key=lambda j:float(eig[j].real))
            levels=[eig[i].real for i in order]
            assert all(x.imag.contains(0) for x in eig)
            assert all(a+2*eps<b for a,b in zip(levels,levels[1:]))
            cols=[sub(X,outside,[z]) for X in (A,B,C)]
            norms=[sum(abs(x)**2 for x in b.entries()).sqrt().upper() for b in cols]
            coeff=[]
            for i,j in enumerate(order):
                v=sub(vec,range(len(levels)),[j])
                v=v/(adj(v)*v)[0,0].real.sqrt()
                sep=min(abs(levels[i]-b).lower() for t,b in enumerate(levels) if t!=i)
                # 对齐两单位本征向量的相位，再用 sin 主角估计对所有软亏额的分量误差。
                err=arb(2).sqrt()*eps/(sep-eps)
                abc=[]
                for col,norm in zip(cols,norms):
                    rad=arb(0,(err*norm).upper())
                    abc.append((adj(v)*col)[0,0]+acb(rad,rad))
                coeff.append(abc)
            mu=[v+arb(0,eps.upper()) for v in levels]
            assert c-mu[1]>0
            a,b,d=coeff[1]
            ra,rb,rd=abs(a),abs(b),abs(d)
            sides=[ra+rb-rd,ra+rd-rb,rb+rd-ra]
            # 若三段不能闭合，则对所有钟角第2个外态的耦合非零。
            if k:
                excess=max((-v).lower() for v in sides)
                assert excess>arb(".0017")
                print("荷、特征标、不能闭合三角的长度差下界",q,k,excess.str(12))
                continue
            assert all(v>0 for v in sides)
            cs=(rd**2-ra**2-rb**2)/(2*ra*rb)
            assert -1<cs<1
            sn=(1-cs**2).sqrt()
            bounds=[]
            for s in (-1,1):
                u=(a/ra)*(rb/b)*acb(cs,s*sn)
                v=-(a+b*u)/d
                # 这两个解覆盖 b₂=0 的所有钟角。其余外本征方向给出标量 Schur 补。
                schur=c-mu[1]-sum(abs(x+y*u+h*v)**2/(mu[i]-mu[1])
                                   for i,(x,y,h) in enumerate(coeff) if i!=1)
                assert schur>arb("1.99")
                bounds.append(schur.lower())
            print("荷、特征标、两闭合相位的 Schur 余量共同下界",q,k,min(bounds).str(12))
        winners=[]
        for values in corners:
            order=sorted(range(3),key=lambda i:float(values[i]))
            assert all(values[j]-values[order[0]]>arb(".10") for j in order[1:])
            winners.append(order[0])
        assert winners==[0,1 if q in (1,2,6) else 2]
        print("荷、有效钟点 (-1,-1) 与 (-1,+1) 上第四态的特征标",q,*winners)


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root)
