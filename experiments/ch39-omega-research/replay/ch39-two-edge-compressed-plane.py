"""两条原边正交补上十七态压缩的全钟面谱隙与低二平面。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys
from flint import acb,acb_mat,arb,ctx


def adj(A):
    return A.conjugate().transpose()


def sub(A,rows,cols):
    return acb_mat([[A[i,j] for j in cols] for i in rows])


def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


@cache
def clock(i,n):
    # 整段角而非格心，两个钟的差角在同宽的二倍区间内。
    mid=arb.pi()*(2*i+1)/n
    return acb(0,mid+arb(0,(arb.pi()/n).upper())).exp()


def positive(S):
    # 对 Hermitian 矩阵，正顺序主子式等价于正定。
    return (S[0,0].real>0 and (S[0,0].real*S[1,1].real-abs(S[0,1])**2)>0
            and S.det().real>0)


def verify(root):
    light=run_path(str(root/"ch39-native-light-projector-pump.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=144
    edges=((0,2),(5,14))
    other=sorted(set(range(19))-{v for e in edges for v in e})
    inner={i for i,p in enumerate(shell["geo"]["sites"]) if sorted(p)==[13,13,23]}
    assert len(inner)==3 and inner.isdisjoint(v for e in edges for v in e)
    assert inner.isdisjoint({shell["soft"]["pole"],*shell["soft"]["star"]["orbit"]})
    inner_cols=[j for j,i in enumerate(other) if i in inner]
    outer_cols=sorted(set(range(17))-set(inner_cols))
    eps=6*arb("1e-8")
    results=[]
    for q in (1,-4,2,-3,6):
        base,links=light["clock_matrices"](shell,q)
        basis=acb_mat(19,17)
        B=acb_mat(19,2)
        for j,i in enumerate(other):
            basis[i,j]=1
        for k,(i,j) in enumerate(edges):
            assert all(m[i,j]==0 and m[j,i]==0 for m in links)
            B[i,k]=1
            B[j,k]=base[j,i]
            basis[i,15+k]=1/arb(2).sqrt()
            basis[j,15+k]=-B[j,k]/arb(2).sqrt()
        assert all(z.contains(0) for z in (adj(B)*basis).entries())
        assert all(z.contains(0) for z in (adj(basis)*basis-eye(17)).entries())
        A,X,Y=[adj(basis)*m*basis for m in (base,*links)]
        assert all(M[i,j]==0 for M in (X,Y) for i in range(17) for j in range(17)
                   if i not in outer_cols or j not in inner_cols)
        D=sub(A,outer_cols,outer_cols)
        C=sub(A,inner_cols,inner_cols)
        triples=[sub(m,outer_cols,inner_cols) for m in (A,X,Y)]
        # 单位辐边处的阈值加上整段软亏额的范数，以便保留 8.22 的严格界。
        t=arb("8.22")+eps
        vals=D.eig(multiple=True)
        assert all(z.imag.contains(0) for z in vals)
        levels=sorted([z.real for z in vals],key=float)
        assert levels[0]<levels[1]<arb("8.178")<t<levels[2]
        R=(D-t*eye(14)).inv()
        S0=C-t*eye(3)
        assert positive(S0)
        # Schur 补的非常数 Fourier 项只有 u、v 和 v*conj(u) 及其共轭。
        schur=S0
        for H in triples:
            schur-=adj(H)*R*H
        H0,H1,H2=triples
        M1,M2,M3=[-adj(a)*R*b for a,b in ((H0,H1),(H0,H2),(H1,H2))]
        stack=[(i,j,8) for i in range(8) for j in range(8)]
        leaves=0
        depth=0
        while stack:
            i,j,n=stack.pop()
            u,v=clock(i,n),clock(j,n)
            T=u*M1+v*M2+u.conjugate()*v*M3
            if positive(schur+T+adj(T)):
                leaves+=1
                depth=max(depth,n)
            else:
                assert n<2**12
                stack.extend((2*i+x,2*j+y,n*2) for x in (0,1) for y in (0,1))
        print("荷、外块第二值上界、全两钟 Schur 证书叶格数与最细分格",q,
              levels[1].upper().str(15),leaves,depth,flush=True)
        results.append(levels[1].upper())
    print("全物理亏额与全耦合缩放下压缩第二值上界、第三值严格下界、统一隙下界",
          max(results).str(15),"8.22",(arb("8.22")-max(results)).lower().str(15))


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root)
